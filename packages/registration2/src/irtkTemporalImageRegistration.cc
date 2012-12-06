#include <irtkRegistration2.h>

#include <irtkResamplingWithPadding.h>

#include <irtkGaussianBlurring2D.h>
#include <irtkGaussianBlurring.h>

#include <irtkGradientImageFilter.h>

#ifdef WIN32
#include <time.h>
#else
#include <sys/resource.h>
#endif

#define MAX_NO_LINE_ITERATIONS 20

irtkGreyImage **tmp_ttarget;
extern irtkGreyImage *tmp_source;

inline double GetBasisSplineValue(double x)
{
    x = fabs(x);
    double value = 0.0;
    if (x < 2.0) {
        if (x < 1.0) {
            value = (double)(2.0f/3.0f + (0.5f*x-1.0)*x*x);
        } else {
            x -= 2.0f;
            value = -x*x*x/6.0f;
        }
    }
    return value;
}

inline double GetBasisSplineDerivativeValue(double x)
{
    double y = fabs(x);
    double value = 0.0;
    if(y < 2.0) {
        if(y < 1.0) {
            value = (double)((1.5*y-2.0)*x);
        }  else {
            y -= 2.0;
            value = -0.5 * y * y;
            if(x < 0.0) value =-value;
        }
    }
    return value;
}

irtkTemporalImageRegistration::irtkTemporalImageRegistration()
{
    int i;

    for (i = 0; i < MAX_NO_RESOLUTIONS; i++) {
        // Default parameters for target image
        _TargetBlurring[i]      = 0;
        _TargetResolution[i][0] = 0;
        _TargetResolution[i][1] = 0;
        _TargetResolution[i][2] = 0;

        // Default parameters for source image
        _SourceBlurring[i]      = 0;
        _SourceResolution[i][0] = 0;
        _SourceResolution[i][1] = 0;
        _SourceResolution[i][2] = 0;

        // Default parameters for optimization
        _NumberOfIterations[i]  = 20;
        _MinStep[i]             = 0.01;
        _MaxStep[i]             = 1;

        // Default dimension of subdivision (= 3D)
        _SubdivisionDim[i]      = 3;
    }

    // Default parameters for registration
    _NumberOfLevels     = 1;
    _NumberOfBins       = 64;

    // Default parameters for optimization
    _SimilarityMeasure  = SSD;
    _InterpolationMode  = Interpolation_Linear;
    _Epsilon            = 0;

    // Default parameters for debugging
    _DebugFlag = false;

    // Set parameters
    _TargetPadding   = MIN_GREY;
    _SourcePadding   = MIN_GREY;

    // Set inputs
    _target = NULL;
    _source = NULL;

    // Set output
    _transformation = NULL;

    // Set metric
    _histogram = NULL;

    // Allocate interpolation object
    _interpolator = NULL;

    _N_target = 1;
}

irtkTemporalImageRegistration::~irtkTemporalImageRegistration()
{}

void irtkTemporalImageRegistration::Initialize()
{
    // Check that the t-dimensions are the same for both images.
    // Do not currently support different t-dimensions.
    for (int n=0; n<_N_target; n++)
        if ((_target[n]->GetT() > 1) || (_source->GetT() > 1)) {
            cerr << this->NameOfClass() << "::Initialize() : Not implemented for images with t > 1" << endl;
            exit(1);
        }
}

void irtkTemporalImageRegistration::Initialize(int level)
{
    double dxT, dyT, dzT, dxS, dyS, dzS, tempT, tempS;
    int i, j, k, t, n, target_nbins, source_nbins;

    int blurrZ = false;
    for (n = 0; n < _N_target; n++) {
        if (_target[n]->GetZ() == 1) {
            blurrZ = true;
            break;
        }
    }

    // Swap source with temp space copies
    tmp_source = new irtkGreyImage(*_source);

    // Swap source with temp space copies
    swap(tmp_source, _source);

    // Blur images if necessary
    if (_SourceBlurring[level] > 0) {
        cout << "Blurring source ... "; cout.flush();
        irtkGaussianBlurringWithPadding<irtkGreyPixel> blurring(_SourceBlurring[level], _SourcePadding);
        blurring.SetInput (_source);
        blurring.SetOutput(_source);
        blurring.Run();
        cout << "done" << endl;
    }

    _source->GetPixelSize(&dxS, &dyS, &dzS);
    tempS = fabs(_SourceResolution[level][0]-dxS) 
        + fabs(_SourceResolution[level][1]-dyS) 
        + fabs(_SourceResolution[level][2]-dzS);

    if (level > 0 || tempS > 0.000001) {
        cout << "Resampling source ... "; cout.flush();
        // Create resampling filter
        irtkResamplingWithPadding<irtkGreyPixel> resample(_SourceResolution[level][0],
            _SourceResolution[level][1],
            _SourceResolution[level][2], _SourcePadding);

        resample.SetInput (_source);
        resample.SetOutput(_source);
        resample.Run();
        cout << "done" << endl;
    }

    tmp_ttarget = new irtkGreyImage *[_N_target];
    for (n = 0; n < _N_target; n++) {
        // Copy target to temp space
        tmp_ttarget[n] = new irtkGreyImage(*_target[n]);
        // Swap target with temp space copies
        swap(tmp_ttarget[n], _target[n]);
    }

    _target_max = MIN_GREY;
    _target_min = MAX_GREY;

    irtkResamplingWithPadding<irtkGreyPixel> resample(_TargetResolution[level][0],
        _TargetResolution[level][1],
        _TargetResolution[level][2],
        _TargetPadding);

    // loop over all target images
    cout << "Blurring (and resampling) target ... "; cout.flush();
    for (n=0; n<_N_target; n++) {

        // Blur images if necessary
        if (_TargetBlurring[level] > 0) {
            if (blurrZ) {
                irtkGaussianBlurringWithPadding2D<irtkGreyPixel> blurring(_TargetBlurring[level], _TargetPadding);
                blurring.SetInput (_target[n]);
                blurring.SetOutput(_target[n]);
                blurring.Run();
            } else {
                irtkGaussianBlurringWithPadding<irtkGreyPixel> blurring(_TargetBlurring[level], _TargetPadding);
                blurring.SetInput (_target[n]);
                blurring.SetOutput(_target[n]);
                blurring.Run();
            }
        }

        _target[n]->GetPixelSize(&dxT, &dyT, &dzT);
        tempT = fabs(_TargetResolution[level][0]-dxT) 
            + fabs(_TargetResolution[level][1]-dyT) 
            + fabs(_TargetResolution[level][2]-dzT);

        if (level > 0 || tempT > 0.000001) {
            // Create resampling filter
            resample.SetInput (_target[n]);
            resample.SetOutput(_target[n]);
            resample.Run();
        }

        // Find out the min and max values in target image, ignoring padding
        for (t = 0; t < _target[n]->GetT(); t++) {
            for (k = 0; k < _target[n]->GetZ(); k++) {
                for (j = 0; j < _target[n]->GetY(); j++) {
                    for (i = 0; i < _target[n]->GetX(); i++) {
                        if (_target[n]->Get(i, j, k, t) > _TargetPadding) {
                            if (_target[n]->Get(i, j, k, t) > _target_max)
                                _target_max = _target[n]->Get(i, j, k, t);
                            if (_target[n]->Get(i, j, k, t) < _target_min)
                                _target_min = _target[n]->Get(i, j, k, t);
                        } else {
                            _target[n]->Put(i, j, k, t, _TargetPadding);
                        }
                    }
                }
            }
        }

    }
    cout << "done" << endl;

    // Find out the min and max values in source image
    _source_max = MIN_GREY;
    _source_min = MAX_GREY;
    for (t = 0; t < _source->GetT(); t++) {
        for (k = 0; k < _source->GetZ(); k++) {
            for (j = 0; j < _source->GetY(); j++) {
                for (i = 0; i < _source->GetX(); i++) {
                    if (_source->Get(i, j, k, t) > _SourcePadding) {
                        if (_source->Get(i, j, k, t) > _source_max)
                            _source_max = _source->Get(i, j, k, t);
                        if (_source->Get(i, j, k, t) < _source_min)
                            _source_min = _source->Get(i, j, k, t);
                    } else {
                        _source->Put(i, j, k, t, _SourcePadding);
                    }
                }
            }
        }
    }
    if(_source_max == MIN_GREY && _source_min == MAX_GREY){
        _source_max = _SourcePadding;
        _source_min = _SourcePadding;
    }

    // Check whether dynamic range of data is not to large
    if (_target_max - _target_min > MAX_GREY) {
        cerr << this->NameOfClass()
            << "::Initialize: Dynamic range of target is too large" << endl;
        exit(1);
    } else {
        for (n = 0; n < _N_target; n++) {
            for (t = 0; t < _target[n]->GetT(); t++) {
                for (k = 0; k < _target[n]->GetZ(); k++) {
                    for (j = 0; j < _target[n]->GetY(); j++) {
                        for (i = 0; i < _target[n]->GetX(); i++) {
                            if (_target[n]->Get(i, j, k, t) > _TargetPadding) {
                                _target[n]->Put(i, j, k, t, _target[n]->Get(i, j, k, t) - _target_min);
                            } else {
                                _target[n]->Put(i, j, k, t, -1);
                            }
                        }
                    }
                }
            }
        }
    }

    // Compute the maximum possible difference across target and source
    _maxDiff = (_target_min - _source_max) * (_target_min - _source_max) > (_target_max - _source_min) * (_target_max - _source_min) ?
        (_target_min - _source_max) * (_target_min - _source_max) : (_target_max - _source_min) * (_target_max - _source_min);

    if (_SimilarityMeasure == SSD) {
        if (_source_max - _target_min > MAX_GREY) {
            cerr << this->NameOfClass()
                << "::Initialize: Dynamic range of source is too large" << endl;
            exit(1);
        } else {
            for (t = 0; t < _source->GetT(); t++) {
                for (k = 0; k < _source->GetZ(); k++) {
                    for (j = 0; j < _source->GetY(); j++) {
                        for (i = 0; i < _source->GetX(); i++) {
                            if ( _source->Get(i, j, k, t) > _SourcePadding) {
                                _source->Put(i, j, k, t, _source->Get(i, j, k, t) - _target_min);
                            }else {
                                _source->Put(i, j, k, t, -1);
                            }
                        }
                    }
                }
            }
        }
    } else {
        if (_source_max - _source_min > MAX_GREY) {
            cerr << this->NameOfClass()
                << "::Initialize: Dynamic range of source is too large" << endl;
            exit(1);
        } else {
            for (t = 0; t < _source->GetT(); t++) {
                for (k = 0; k < _source->GetZ(); k++) {
                    for (j = 0; j < _source->GetY(); j++) {
                        for (i = 0; i < _source->GetX(); i++) {
                            if ( _source->Get(i, j, k, t) > _SourcePadding) {
                                _source->Put(i, j, k, t, _source->Get(i, j, k, t) - _source_min);
                            }else {
                                _source->Put(i, j, k, t, -1);
                            }
                        }
                    }
                }
            }
        }
    }

    // Pad target image if necessary
    for (n = 0; n < _N_target; n++)
        irtkPadding(*_target[n], _TargetPadding);

    irtkPadding(*_source, _SourcePadding);

    // Allocate memory for histogram if necessary
    switch (_SimilarityMeasure) {
    case SSD:
        break;
    case JE:
    case MI:
    case NMI:
        //      // Rescale images by an integer factor if necessary
        //      target_nbins = irtkCalculateNumberOfBins(_target, _NumberOfBins, _target_min, _target_max);
        //      source_nbins = irtkCalculateNumberOfBins(_source, _NumberOfBins, _source_min, _source_max);
        //      _histogram = new irtkHistogram_2D<double>(target_nbins, source_nbins);
        //      break;
    default:
        cerr << this->NameOfClass() << "::Initialize(int): No such metric implemented" << endl;
        exit(1);
    }


    // Compute spatial gradient of source image
    irtkGradientImageFilter<double> gradient(irtkGradientImageFilter<double>::GRADIENT_VECTOR);
    irtkGenericImage<double> tmp = *_source;
    gradient.SetInput (&tmp);
    gradient.SetOutput(&_sourceGradient);
    gradient.SetPadding(-1);
    gradient.Run();

    // Determine attributes of target image
    irtkImageAttributes attr_t;

    // Set up transformed source image
    _transformedSource = new irtkGenericImage<double> [_N_target];
    for (n = 0; n < _N_target; n++) {
        attr_t = _target[n]->GetImageAttributes();
        _transformedSource[n].Initialize(attr_t);
    }

    // Set up gradient of source image (transformed)
    _transformedSourceGradient = new irtkGenericImage<double> [_N_target];
    for (n = 0; n < _N_target; n++) {
        attr_t = _target[n]->GetImageAttributes();
        attr_t._t = 3;
        _transformedSourceGradient[n].Initialize(attr_t);
    }

    // Set up gradient of similarity metric
    _similarityGradient = new irtkGenericImage<double> [_N_target];
    for (n = 0; n < _N_target; n++) {
        attr_t = _target[n]->GetImageAttributes();
        attr_t._t = 3;
        _similarityGradient[n].Initialize(attr_t);
    }

    // Setup the interpolator
    _interpolator = irtkInterpolateImageFunction::New(_InterpolationMode, _source);

    // Setup interpolation for the source image
    _interpolator->SetInput(_source);
    _interpolator->Initialize();

    // Calculate the source image domain in which we can interpolate
    _interpolator->Inside(_source_x1, _source_y1, _source_z1,
        _source_x2, _source_y2, _source_z2);

    // Setup the interpolator
    _interpolatorGradient = irtkInterpolateImageFunction::New(_InterpolationMode, &_sourceGradient);

    // Setup interpolation for the source image
    _interpolatorGradient->SetInput(&_sourceGradient);
    _interpolatorGradient->Initialize();

    // Calculate the source image domain in which we can interpolate
    _interpolatorGradient->Inside(_source_x1, _source_y1, _source_z1,
        _source_x2, _source_y2, _source_z2);

    //Print some debugging information
    for (int n=0; n<_N_target; n++){
        cout << "Target images"<<n<<" (reference)" << endl;
        _target[n]->Print();
    }
    cout << "Targets range is from " << _target_min << " to " << _target_max << endl;

    cout << "Source image (transform)" << endl;
    _source->Print();
    cout << "Range is from " << _source_min << " to " << _source_max << endl;

    // Print initial transformation
    cout << "Initial transformation for level = " << level+1 << endl;;
    _transformation->Print();

}

void irtkTemporalImageRegistration::Finalize()
{}

void irtkTemporalImageRegistration::Finalize(int level)
{
    int i;
    // Print final transformation
    cout << "Final transformations for level = " << level+1 << endl;
    _transformation->Print();

    // Swap source, target and time array back with temp space copies (see Initialize)
    for (i = 0; i < _N_target; i++) {
        swap(tmp_ttarget[i], _target[i]);
    }
    swap(tmp_source, _source);

    for (i = 0; i < _N_target; i++) {
        delete tmp_ttarget[i];
        tmp_ttarget[i] = NULL;
    }
    delete [] tmp_ttarget;
    tmp_ttarget = NULL;
    delete tmp_source;
    tmp_source = NULL;
    delete _interpolator;
    _interpolator = NULL;
    delete _interpolatorGradient;
    _interpolatorGradient = NULL;
    if (_histogram != NULL) {
        delete _histogram;
        _histogram = NULL;
    }
    delete []_similarityGradient;
    _similarityGradient = NULL;
    delete []_transformedSourceGradient;
    _transformedSourceGradient = NULL;
    delete []_transformedSource;
    _transformedSource = NULL;
}

void irtkTemporalImageRegistration::Update(bool updateGradient)
{
    // Update
    if (updateGradient == true) {
        this->UpdateSourceAndGradient();
    } else {
        this->UpdateSource();
    }
}

void irtkTemporalImageRegistration::Run()
{
    int i, k;
    char buffer[256];
    double *gradient, delta, step, min_step, max_step, max_length, best_similarity, new_similarity, old_similarity;

    // Print debugging information
    this->Debug("irtkTemporalImageRegistration::Run");

    if (_source == NULL) {
        cerr << "irtkTemporalImageRegistration::Run: Filter has no source input" << endl;
        exit(1);
    }

    if (_target == NULL) {
        cerr << "irtkTemporalImageRegistration::Run: Filter has no target input" << endl;
        exit(1);
    }

    if (_transformation == NULL) {
        cerr << "irtkTemporalImageRegistration::Run: Filter has no transformation output" << endl;
        exit(1);
    }

    // Do the initial set up for all levels
    this->Initialize();

    // Loop over levels
    for (_CurrentLevel = _NumberOfLevels-1; _CurrentLevel >= 0; _CurrentLevel--) {

        // Initial step size
        min_step = _MinStep[_CurrentLevel];
        max_step = _MaxStep[_CurrentLevel];

        // Print resolution level
        cout << "Resolution level no. " << _CurrentLevel+1 << " (step sizes " << min_step << " to " << max_step  << ")\n";

        // Initialize for this level
        this->Initialize(_CurrentLevel);

        // Save pre-processed images if we are debugging
        sprintf(buffer, "source_%d.nii.gz", _CurrentLevel);
        if (_DebugFlag == true)
            _source->Write(buffer);
        if (_DebugFlag == true)
            for (int n = 0; n < _N_target; n++) {
                sprintf(buffer, "target_%d_t%d.nii.gz", _CurrentLevel, n);
                _target[n]->Write(buffer);
            }

            // Allocate memory for gradient vector
            gradient = new double[_transformation->NumberOfDOFs()];

            // Run the registration filter at this resolution
            _CurrentIteration = 0;
            while (_CurrentIteration < _NumberOfIterations[_CurrentLevel]) {
                cout << "Iteration = " << _CurrentIteration + 1 << " (out of " << _NumberOfIterations[_CurrentLevel] << ")"<< endl;

                // Update source image
                this->Update(true);

                // Compute current metric value
                best_similarity = old_similarity = this->Evaluate();
                cout << "Current best metric value is " << best_similarity << endl;

                // Compute gradient of similarity metric
                max_length = this->EvaluateGradient(gradient);

                // Step along gradient direction until no further improvement is necessary
                i = 0;
                delta = 0;
                step = max_step;
                do {
                    double current = step / max_length;

                    // Move along gradient direction
                    for (k = 0; k < _transformation->NumberOfDOFs(); k++) {
                        _transformation->Put(k, _transformation->Get(k) + current * gradient[k]);
                    }

                    // We have just changed the transformation parameters, so we need to update
                    this->Update(false);

                    // Compute new similarity
                    new_similarity = this->Evaluate();

                    if (new_similarity > best_similarity + _Epsilon) {
                        cout << "New metric value is " << new_similarity << "; step = " << step << endl;
                        best_similarity = new_similarity;
                        delta += step;
                        break; // recalculate gradient if current one was succesful
                        step = step * 1.1;
                        if (step > max_step) step = max_step;

                    } else {
                        // Last step was no improvement, so back track
                        cout << "Rejected metric value is " << new_similarity << "; step = " << step << endl;
                        for (k = 0; k < _transformation->NumberOfDOFs(); k++) {
                            _transformation->Put(k, _transformation->Get(k) - current * gradient[k]);
                        }
                        step = step * 0.5;
                    }
                    i++;
                } while ((i < MAX_NO_LINE_ITERATIONS) && (step > min_step));
                _CurrentIteration++;

                // Check for convergence
                if (delta == 0) break;
            }

            // Delete gradient
            delete gradient;

            // Do the final cleaning up for this level
            this->Finalize(_CurrentLevel);
    }

    // Do the final cleaning up for all levels
    this->Finalize();
}

double irtkTemporalImageRegistration::EvaluateSSD()
{
    int i, n, nV;
    double norm, ssd;

    // Print debugging information
    this->Debug("irtkTemporalImageRegistration::EvaluateSSD");

    // Initialize metric
    nV = 0;
    ssd = 0;
    //cout<<"irtkTemporalImageRegistration::EvaluateSSD start"<<endl;

    for (n = 0; n < _N_target; n++) {
        // Pointer to voxels in images
        short  *ptr2target = _target[n]->GetPointerToVoxels();
        double *ptr2source = _transformedSource[n].GetPointerToVoxels();
        for (i = 0; i < _target[n]->GetNumberOfVoxels(); i++) {
            if ((*ptr2target >= 0) && (*ptr2source >= 0)) {
                ssd += (*ptr2target - *ptr2source) * (*ptr2target - *ptr2source);
                nV++;
            }
            ptr2target++;
            ptr2source++;
        }
    }
    //cout<<"irtkTemporalImageRegistration::EvaluateSSD end"<<endl;

    // Normalize similarity measure by number of voxels and maximum SSD
    norm = 1.0 / ((double)nV * (double)_maxDiff);

    // Return similarity measure
    if (norm > 0) {
        return -ssd * norm;
    } else {
        cerr << "irtkTemporalImageRegistration::EvaluateSSD: No samples available" << endl;
        return 0;
    }
}

double irtkTemporalImageRegistration::EvaluateNMI()
{
    cerr << "irtkTemporalImageRegistration::EvaluateNMI: NMI not implemented yet!" << endl;
    exit(1);
    //  int i;
    //
    //  // Print debugging information
    //  this->Debug("irtkTemporalImageRegistration::EvaluateNMI");
    //
    //  // Pointer to voxels in images
    //  short  *ptr2target = _target->GetPointerToVoxels();
    //  double *ptr2source = _transformedSource.GetPointerToVoxels();
    //
    //  // Initialize metric
    //  _histogram->Reset();
    //
    //  // Compute metric
    //  for (i = 0; i < _target->GetNumberOfVoxels(); i++) {
    //    if ((*ptr2target >= 0) && (*ptr2source >= 0)) {
    //      _histogram->Add(*ptr2target, round(*ptr2source));
    //    }
    //    ptr2target++;
    //    ptr2source++;
    //  }
    //
    //  // Smooth histogram if appropriate
    //  _histogram->Smooth();
    //
    //  // Evaluate similarity measure
    return _histogram->NormalizedMutualInformation();
}

double irtkTemporalImageRegistration::Evaluate()
{
    double metric;

#ifdef USE_TIMING
    // Start timing
    clock_t start, end;
    double cpu_time_used;
    start = clock();
#endif

    // Print debugging information
    this->Debug("irtkTemporalImageRegistration::Evaluate");

    // Allocate memory for metric
    switch (_SimilarityMeasure) {
    case SSD:
        metric = this->EvaluateSSD();
        break;
    case NMI:
        metric = this->EvaluateNMI();
        break;
    default:
        metric = 0;
        cerr << this->NameOfClass() << "::Evaluate: No such metric implemented" << endl;
        exit(1);
    }

#ifdef USE_TIMING
    // Stop timing
    end = clock();
    cpu_time_used = ((double) (end - start)) / CLOCKS_PER_SEC;
    cout << "CPU time for irtkTemporalImageRegistration::Evaluate() = " << cpu_time_used << endl;
#endif

    // Evaluate similarity measure
    return metric;
}

void irtkTemporalImageRegistration::EvaluateGradientSSD()
{
    double ssd;
    int i, j, k, n;

    // Print debugging information
    this->Debug("irtkTemporalImageRegistration::EvaluateGradient");
    //cout<<"irtkTemporalImageRegistration::EvaluateGradientSSD start"<<endl;

    for (n = 0; n < _N_target; n++) {
        // Pointer to voxels in images
        short  *ptr2target = _target[n]->GetPointerToVoxels();
        double *ptr2source = _transformedSource[n].GetPointerToVoxels();
        // Compute gradient
        for (k = 0; k < _target[n]->GetZ(); k++) {
            for (j = 0; j < _target[n]->GetY(); j++) {
                for (i = 0; i < _target[n]->GetX(); i++) {
                    if ((*ptr2target >= 0) && (*ptr2source >= 0)) {
                        ssd = 2 * (*ptr2target - *ptr2source) / _maxDiff;
                        _similarityGradient[n](i, j, k, 0) = ssd * _transformedSourceGradient[n](i, j, k, 0);
                        _similarityGradient[n](i, j, k, 1) = ssd * _transformedSourceGradient[n](i, j, k, 1);
                        _similarityGradient[n](i, j, k, 2) = ssd * _transformedSourceGradient[n](i, j, k, 2);
                    } else {
                        _similarityGradient[n](i, j, k, 0) = 0;
                        _similarityGradient[n](i, j, k, 1) = 0;
                        _similarityGradient[n](i, j, k, 2) = 0;
                    }
                    ptr2target++;
                    ptr2source++;
                }
            }
        }
    }
    //cout<<"irtkTemporalImageRegistration::EvaluateGradientSSD end"<<endl;
}

void irtkTemporalImageRegistration::EvaluateGradientNMI()
{
    cerr << "irtkTemporalImageRegistration::EvaluateGradientNMI: NMI not implemented yet!" << endl;
    exit(1);
    //  int i, j, k, l, t, r;
    //  double w, je, nmi, tmp, targetEntropyGrad[3], sourceEntropyGrad[3], jointEntropyGrad[3];
    //
    //  // Allocate new histograms
    //  irtkHistogram_1D<double> logMarginalXHistogram(_histogram->NumberOfBinsX());
    //  irtkHistogram_1D<double> logMarginalYHistogram(_histogram->NumberOfBinsY());
    //  irtkHistogram_2D<double> logJointHistogram(_histogram->NumberOfBinsX(), _histogram->NumberOfBinsY());
    //
    //  // Recompute joint histogram
    //  for (j = 0; j < _histogram->NumberOfBinsY(); j++) {
    //    for (i = 0; i < _histogram->NumberOfBinsX(); i++) {
    //
    //      logJointHistogram.Add(i, j, _histogram->irtkHistogram_2D<double>::operator()(i, j));
    //    }
    //  }
    //
    //  // Smooth joint histogram
    //  //logJointHistogram.Smooth();
    //  je  = logJointHistogram.JointEntropy();
    //  nmi = logJointHistogram.NormalizedMutualInformation();
    //
    //  // Recompute marginal histogram for X
    //  for (i = 0; i < _histogram->NumberOfBinsX(); i++) {
    //    for (j = 0; j < _histogram->NumberOfBinsY(); j++) {
    //      logMarginalXHistogram.Add(i, _histogram->irtkHistogram_2D<double>::operator()(i, j));
    //    }
    //  }
    //  // Recompute marginal histogram for Y
    //  for (j = 0; j < _histogram->NumberOfBinsY(); j++) {
    //    for (i = 0; i < _histogram->NumberOfBinsX(); i++) {
    //      logMarginalYHistogram.Add(j, _histogram->irtkHistogram_2D<double>::operator()(i, j));
    //    }
    //  }
    //
    //  // Log transform histograms
    //  logJointHistogram.Log();
    //  logMarginalXHistogram.Log();
    //  logMarginalYHistogram.Log();
    //
    //  // Pointer to voxels in images
    //  short *ptr2target  = _target->GetPointerToVoxels();
    //  double *ptr2source = _transformedSource.GetPointerToVoxels();
    //
    //  // Loop over images
    //  for (k = 0; k < _target->GetZ(); k++) {
    //    for (j = 0; j < _target->GetY(); j++) {
    //      for (i = 0; i < _target->GetX(); i++) {
    //
    //        // This code is based on an idea from Marc Modat for computing the NMI derivative as suggested in his niftyreg package
    //        if ((*ptr2target >= 0) && (*ptr2source >= 0)) {
    //          irtkGreyPixel targetValue = *ptr2target;
    //          irtkRealPixel sourceValue = *ptr2source;
    //
    //          if ((targetValue >= 0) && (sourceValue >= 0.0f) && (targetValue < _histogram->NumberOfBinsX()) && (sourceValue < _histogram->NumberOfBinsY())) {
    //            // The two is added because the image is resample between 2 and bin +2
    //            // if 64 bins are used the histogram will have 68 bins et the image will be between 2 and 65
    //
    //            sourceValue = (float)floor((double)sourceValue);
    //
    //            for (l = 0; l < 3; l++) {
    //              jointEntropyGrad [l] = 0;
    //              targetEntropyGrad[l] = 0;
    //              sourceEntropyGrad[l] = 0;
    //            }
    //
    //            for (t = targetValue-1; t < targetValue+2; t++) {
    //              if ((t >= 0) && (t < _histogram->NumberOfBinsX())) {
    //                for (r = (int)(sourceValue-1.0); r < (int)(sourceValue+2.0); r++) {
    //                  if ((r >= 0) && (r < _histogram->NumberOfBinsY())) {
    //                    w = GetBasisSplineValue((double)t-(double)targetValue) * GetBasisSplineDerivativeValue((double)r-(double)sourceValue);
    //
    //                    double jointLog  = logJointHistogram(t, r);
    //                    double targetLog = logMarginalXHistogram(t);
    //                    double resultLog = logMarginalYHistogram(r);
    //
    //                    tmp = -w * _transformedSourceGradient(i, j, k, 0);
    //                    jointEntropyGrad[0]  -= tmp * jointLog;
    //                    targetEntropyGrad[0] -= tmp * targetLog;
    //                    sourceEntropyGrad[0] -= tmp * resultLog;
    //
    //                    tmp = -w * _transformedSourceGradient(i, j, k, 1);
    //                    jointEntropyGrad[1]  -= tmp * jointLog;
    //                    targetEntropyGrad[1] -= tmp * targetLog;
    //                    sourceEntropyGrad[1] -= tmp * resultLog;
    //
    //                    tmp = -w * _transformedSourceGradient(i, j, k, 2);
    //                    jointEntropyGrad[2]  -= tmp * jointLog;
    //                    targetEntropyGrad[2] -= tmp * targetLog;
    //                    sourceEntropyGrad[2] -= tmp * resultLog;
    //
    //                  }
    //                }
    //              }
    //            }
    //
    //            _similarityGradient(i, j, k, 0) = ((targetEntropyGrad[0] + sourceEntropyGrad[0] - nmi * jointEntropyGrad[0]) / je);
    //            _similarityGradient(i, j, k, 1) = ((targetEntropyGrad[1] + sourceEntropyGrad[1] - nmi * jointEntropyGrad[1]) / je);
    //            _similarityGradient(i, j, k, 2) = ((targetEntropyGrad[2] + sourceEntropyGrad[2] - nmi * jointEntropyGrad[2]) / je);
    //          } else {
    //            _similarityGradient(i, j, k, 0) = 0;
    //            _similarityGradient(i, j, k, 1) = 0;
    //            _similarityGradient(i, j, k, 2) = 0;
    //          }
    //        }else {
    //            _similarityGradient(i, j, k, 0) = 0;
    //            _similarityGradient(i, j, k, 1) = 0;
    //            _similarityGradient(i, j, k, 2) = 0;
    //        }
    //        ptr2target++;
    //        ptr2source++;
    //      }
    //    }
    //  }
}

double irtkTemporalImageRegistration::EvaluateGradient(double *)
{
    int i, j, k, n;
    double x, y, z;

    // Start timing
    clock_t start, end;
    double cpu_time_used;
    start = clock();

    // Allocate memory for metric
    switch (_SimilarityMeasure) {
    case SSD:
        this->EvaluateGradientSSD();
        break;
    case NMI:
        this->EvaluateGradientNMI();
        break;
    default:
        cerr << this->NameOfClass() << "::Evaluate: No such metric implemented" << endl;
        exit(1);
    }

    // Extract matrix for reorientation of gradient
    irtkMatrix m = _source->GetImageToWorldMatrix();
    //cout<<"irtkTemporalImageRegistration::EvaluateGradient start"<<endl;
    for (n = 0; n < _N_target; n++) {
        // Pointer to voxels in images
        short  *ptr2target = _target[n]->GetPointerToVoxels();
        double *ptr2source = _transformedSource[n].GetPointerToVoxels();
        // Reorient gradient
        for (k = 0; k < _target[n]->GetZ(); k++) {
            for (j = 0; j < _target[n]->GetY(); j++) {
                for (i = 0; i < _target[n]->GetX(); i++) {
                    if ((*ptr2target >= 0) && (*ptr2source >= 0)) {
                        x = m(0, 0) * _similarityGradient[n](i, j, k, 0) + m(0, 1) * _similarityGradient[n](i, j, k, 1) + m(0, 2) * _similarityGradient[n](i, j, k, 2);
                        y = m(1, 0) * _similarityGradient[n](i, j, k, 0) + m(1, 1) * _similarityGradient[n](i, j, k, 1) + m(1, 2) * _similarityGradient[n](i, j, k, 2);
                        z = m(2, 0) * _similarityGradient[n](i, j, k, 0) + m(2, 1) * _similarityGradient[n](i, j, k, 1) + m(2, 2) * _similarityGradient[n](i, j, k, 2);
                        _similarityGradient[n](i, j, k, 0) = x;
                        _similarityGradient[n](i, j, k, 1) = y;
                        _similarityGradient[n](i, j, k, 2) = z;
                    }
                    ptr2target++;
                    ptr2source++;
                }
            }
        }
    }
    //cout<<"irtkTemporalImageRegistration::EvaluateGradient end"<<endl;
    // Stop timing
    end = clock();
    cpu_time_used = ((double) (end - start)) / CLOCKS_PER_SEC;
    //cout << "CPU time for irtkTemporalImageRegistration::EvaluateGradient() = " << cpu_time_used << endl;

    // This function always returns 0
    return 0;
}

bool irtkTemporalImageRegistration::Read(char *buffer1, char *buffer2, int &level)
{
    int i, n, ok = false;
    double dx, dy, dz;

    // Resolution level
    if (strstr(buffer1, "Resolution level") != NULL) {
        level = atoi(buffer2)-1;
        ok = true;
    }
    // Target blurring
    if (strstr(buffer1, "Target blurring (in mm)") != NULL) {
        if (level == -1) {
            for (i = 0; i < MAX_NO_RESOLUTIONS; i++) {
                this->_TargetBlurring[i] = pow(2.0, double(i)) * atof(buffer2);
            }
        } else {
            this->_TargetBlurring[level] = atof(buffer2);
        }
        ok = true;
    }

    // Target resolution
    if (strstr(buffer1, "Target resolution (in mm)") != NULL) {
        _target[0]->GetPixelSize(&dx, &dy, &dz);
        if (level == -1) {
            n = sscanf(buffer2, "%lf %lf %lf", &(this->_TargetResolution[0][0]),  &(this->_TargetResolution[0][1]),  &(this->_TargetResolution[0][2]));
            if (n == 1) {
                this->_TargetResolution[0][1] = this->_TargetResolution[0][0];
                this->_TargetResolution[0][2] = this->_TargetResolution[0][0];
            }
            if (this->_TargetResolution[0][0] == 0) this->_TargetResolution[0][0] = dx;
            if (this->_TargetResolution[0][1] == 0) this->_TargetResolution[0][1] = dy;
            if (this->_TargetResolution[0][2] == 0) this->_TargetResolution[0][2] = dz;
            for (i = 0; i < MAX_NO_RESOLUTIONS; i++) {
                this->_TargetResolution[i][0] = pow(2.0, double(i)) * this->_TargetResolution[0][0];
                this->_TargetResolution[i][1] = pow(2.0, double(i)) * this->_TargetResolution[0][1];
                this->_TargetResolution[i][2] = pow(2.0, double(i)) * this->_TargetResolution[0][2];
            }
        } else {
            n = sscanf(buffer2, "%lf %lf %lf", &(this->_TargetResolution[level][0]),  &(this->_TargetResolution[level][1]),  &(this->_TargetResolution[level][2]));
            if (n == 1) {
                this->_TargetResolution[level][1] = this->_TargetResolution[level][0];
                this->_TargetResolution[level][2] = this->_TargetResolution[level][0];
            }
            if (this->_TargetResolution[level][0] == 0) this->_TargetResolution[level][0] = dx;
            if (this->_TargetResolution[level][1] == 0) this->_TargetResolution[level][1] = dy;
            if (this->_TargetResolution[level][2] == 0) this->_TargetResolution[level][2] = dz;
        }
        ok = true;
    }

    // Source blurring
    if (strstr(buffer1, "Source blurring (in mm)") != NULL) {
        if (level == -1) {
            for (i = 0; i < MAX_NO_RESOLUTIONS; i++) {
                this->_SourceBlurring[i] = pow(2.0, double(i)) * atof(buffer2);
            }
        } else {
            this->_SourceBlurring[level] = atof(buffer2);
        }
        ok = true;
    }
    // Source resolution
    if (strstr(buffer1, "Source resolution (in mm)") != NULL) {
        _source->GetPixelSize(&dx, &dy, &dz);
        if (level == -1) {
            n = sscanf(buffer2, "%lf %lf %lf", &(this->_SourceResolution[0][0]),  &(this->_SourceResolution[0][1]),  &(this->_SourceResolution[0][2]));
            if (n == 1) {
                this->_SourceResolution[0][1] = this->_SourceResolution[0][0];
                this->_SourceResolution[0][2] = this->_SourceResolution[0][0];
            }
            if (this->_SourceResolution[0][0] == 0) this->_SourceResolution[0][0] = dx;
            if (this->_SourceResolution[0][1] == 0) this->_SourceResolution[0][1] = dy;
            if (this->_SourceResolution[0][2] == 0) this->_SourceResolution[0][2] = dz;
            for (i = 0; i < MAX_NO_RESOLUTIONS; i++) {
                this->_SourceResolution[i][0] = pow(2.0, double(i)) * this->_SourceResolution[0][0];
                this->_SourceResolution[i][1] = pow(2.0, double(i)) * this->_SourceResolution[0][1];
                this->_SourceResolution[i][2] = pow(2.0, double(i)) * this->_SourceResolution[0][2];
            }
        } else {
            n = sscanf(buffer2, "%lf %lf %lf", &(this->_SourceResolution[level][0]),  &(this->_SourceResolution[level][1]),  &(this->_SourceResolution[level][2]));
            if (n == 1) {
                this->_SourceResolution[level][1] = this->_SourceResolution[level][0];
                this->_SourceResolution[level][2] = this->_SourceResolution[level][0];
            }
            if (this->_SourceResolution[level][0] == 0) this->_SourceResolution[level][0] = dx;
            if (this->_SourceResolution[level][1] == 0) this->_SourceResolution[level][1] = dy;
            if (this->_SourceResolution[level][2] == 0) this->_SourceResolution[level][2] = dz;
        }
        ok = true;
    }
    if (strstr(buffer1, "No. of resolution levels") != NULL) {
        this->_NumberOfLevels = atoi(buffer2);
        ok = true;
    }
    if (strstr(buffer1, "No. of bins") != NULL) {
        this->_NumberOfBins = atoi(buffer2);
        ok = true;
    }
    if (strstr(buffer1, "No. of iterations") != NULL) {
        if (level == -1) {
            for (i = 0; i < MAX_NO_RESOLUTIONS; i++) {
                this->_NumberOfIterations[i] = atoi(buffer2);
            }
        } else {
            this->_NumberOfIterations[level] = atoi(buffer2);
        }
        ok = true;
    }
    if (strstr(buffer1, "Maximum length of steps") != NULL) {
        if (level == -1) {
            for (i = 0; i < MAX_NO_RESOLUTIONS; i++) {
                this->_MaxStep[i] = atof(buffer2);
            }
        } else {
            this->_MaxStep[level] = atof(buffer2);
        }
        ok = true;
    }
    if (strstr(buffer1, "Minimum length of steps") != NULL) {
        if (level == -1) {
            for (i = 0; i < MAX_NO_RESOLUTIONS; i++) {
                this->_MinStep[i] = atof(buffer2);
            }
        } else {
            this->_MinStep[level] = atof(buffer2);
        }
        ok = true;
    }
    if (strstr(buffer1, "Dimension subdivision") != NULL) {
        if (level == -1) {
            for (i = 0; i < MAX_NO_RESOLUTIONS; i++) {
                this->_SubdivisionDim[i] = atoi(buffer2);
            }
        } else {
            this->_SubdivisionDim[level] = atoi(buffer2);
        }
        ok = true;
    }
    if (strstr(buffer1, "Epsilon") != NULL) {
        this->_Epsilon = atof(buffer2);
        ok = true;
    }
    if (strstr(buffer1, "Padding value") != NULL) {
        this->_TargetPadding = atoi(buffer2);
        ok = true;
    }
    if (strstr(buffer1, "Padding value") != NULL) {
        this->_TargetPadding = atoi(buffer2);
        ok = true;
    }
    if (strstr(buffer1, "Source padding value") != NULL) {
        this->_SourcePadding = atoi(buffer2);
        ok = true;
    }
    if (strstr(buffer1, "Similarity measure") != NULL) {
        if (strstr(buffer2, "CC") != NULL) {
            this->_SimilarityMeasure = CC;
            ok = true;
        } else {
            if (strstr(buffer2, "JE") != NULL) {
                this->_SimilarityMeasure = JE;
                ok = true;
            } else {
                if (strstr(buffer2, "NMI") != NULL) {
                    this->_SimilarityMeasure = NMI;
                    ok = true;
                } else {
                    if (strstr(buffer2, "MI") != NULL) {
                        this->_SimilarityMeasure = MI;
                        ok = true;
                    } else {
                        if (strstr(buffer2, "SSD") != NULL) {
                            this->_SimilarityMeasure = SSD;
                            ok = true;
                        } else {
                            if (strstr(buffer2, "CR_XY") != NULL) {
                                this->_SimilarityMeasure = CR_XY;
                                ok = true;
                            } else {
                                if (strstr(buffer2, "CR_YX") != NULL) {
                                    this->_SimilarityMeasure = CR_YX;
                                    ok = true;
                                } else {
                                    if (strstr(buffer2, "LC") != NULL) {
                                        this->_SimilarityMeasure = LC;
                                        ok = true;
                                    } else {
                                        if (strstr(buffer2, "K") != NULL) {
                                            this->_SimilarityMeasure = K;
                                            ok = true;
                                        } else {
                                            if (strstr(buffer2, "ML") != NULL) {
                                                this->_SimilarityMeasure = ML;
                                                ok = true;
                                            }
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
    }
    if (strstr(buffer1, "Interpolation mode") != NULL) {
        if (strstr(buffer2, "NN") != NULL) {
            this->_InterpolationMode = Interpolation_NN;
            cout << "Interpolation Mode is ... NN" << endl;
            ok = true;
        } else {
            if (strstr(buffer2, "Linear") != NULL) {
                this->_InterpolationMode = Interpolation_Linear;
                cout << "Interpolation Mode is ... Linear" << endl;
                ok = true;
            } else {
                if (strstr(buffer2, "CSpline") != NULL) {
                    this->_InterpolationMode = Interpolation_CSpline;
                    cout << "Interpolation Mode is ... CSpline" << endl;
                    ok = true;
                } else {
                    if (strstr(buffer2, "BSpline") != NULL) {
                        this->_InterpolationMode = Interpolation_BSpline;
                        cout << "Interpolation Mode is ... BSpline" << endl;
                        ok = true;
                    } else {
                        if (strstr(buffer2, "Sinc") != NULL) {
                            this->_InterpolationMode = Interpolation_Sinc;
                            cout << "Interpolation Mode is ... Sinc" << endl;
                            ok = true;
                        }
                    }
                }
            }
        }
    }

    if (ok == false) {
        cerr << "irtkTemporalImageRegistration::Read: Can't parse line " << buffer1 << endl;
        exit(1);
    }

    return ok;
}

void irtkTemporalImageRegistration::Write(ostream &to)
{
    int i;

    to << "\n#\n# Registration parameters\n#\n\n";
    to << "No. of resolution levels          = " << this->_NumberOfLevels << endl;
    to << "No. of bins                       = " << this->_NumberOfBins << endl;
    to << "Epsilon                           = " << this->_Epsilon << endl;
    to << "Padding value                     = " << this->_TargetPadding << endl;
    to << "Source padding value              = " << this->_SourcePadding << endl;

    switch (this->_SimilarityMeasure) {
    case K:
        to << "Similarity measure                = K" << endl;
        break;
    case LC:
        to << "Similarity measure                = LC" << endl;
        break;
    case CC:
        to << "Similarity measure                = CC" << endl;
        break;
    case JE:
        to << "Similarity measure                = JE" << endl;
        break;
    case MI:
        to << "Similarity measure                = MI" << endl;
        break;
    case NMI:
        to << "Similarity measure                = NMI" << endl;
        break;
    case SSD:
        to << "Similarity measure                = SSD" << endl;
        break;
    case CR_XY:
        to << "Similarity measure                = CR_XY" << endl;
        break;
    case CR_YX:
        to << "Similarity measure                = CR_YX" << endl;
        break;
    case ML:
        to << "Similarity measure                = ML" << endl;
        break;
    }

    switch (this->_InterpolationMode) {
    case Interpolation_NN:
        to << "Interpolation mode                = NN" << endl;
        break;
    case Interpolation_Linear:
        to << "Interpolation mode                = Linear" << endl;
        break;
    case Interpolation_CSpline:
        to << "Interpolation mode                = CSpline" << endl;
        break;
    case Interpolation_BSpline:
        to << "Interpolation mode                = BSpline" << endl;
        break;
    case Interpolation_Sinc:
        to << "Interpolation mode                = Sinc" << endl;
        break;
    case Interpolation_Gaussian:
        to << "Interpolation mode                = Gaussian" << endl;
        break;
    default:
        cerr << "irtkTemporalImageRegistration::Write: Interpolation mode not supported" << endl;
        exit(1);
    }

    for (i = 0; i < this->_NumberOfLevels; i++) {
        to << "\n#\n# Registration parameters for resolution level " << i+1 << "\n#\n\n";
        to << "Resolution level                  = " << i+1 << endl;
        to << "Target blurring (in mm)           = " << this->_TargetBlurring[i] << endl;
        to << "Target resolution (in mm)         = " << this->_TargetResolution[i][0] << " " << this->_TargetResolution[i][1] << " " << this->_TargetResolution[i][2] << endl;
        to << "Source blurring (in mm)           = " << this->_SourceBlurring[i] << endl;
        to << "Source resolution (in mm)         = " << this->_SourceResolution[i][0] << " " << this->_SourceResolution[i][1] << " " << this->_SourceResolution[i][2] << endl;
        to << "No. of iterations                 = " << this->_NumberOfIterations[i] << endl;
        to << "Minimum length of steps           = " << this->_MinStep[i] << endl;
        to << "Maximum length of steps           = " << this->_MaxStep[i] << endl;
        to << "Dimension subdivision             = " << this->_SubdivisionDim[i] << endl;
    }
}

void irtkTemporalImageRegistration::Read(char *filename)
{
    int level;
    char buffer1[255], *buffer2;

    ifstream from(filename);

    if (!from) {
        cerr << "irtkTemporalImageRegistration::Read: Can't open file " << filename
            << endl;
        exit(1);
    }

    level = -1;
    while (from.eof() != true) {
        if (read_line(from, buffer1, buffer2) != 0) {
            if (this->Read(buffer1, buffer2, level) == false) {
                cerr << "Couldn't parse line: " << buffer1 << endl;
            }
        }
    }
}

void irtkTemporalImageRegistration::Write(char *filename)
{
    ofstream to(filename);

    if (!to) {
        cerr << "irtkTemporalImageRegistration::Write: Can't open file " << filename;
        exit(1);
    }

    this->Write(to);
}
