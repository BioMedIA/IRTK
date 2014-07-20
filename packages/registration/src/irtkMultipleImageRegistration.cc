/*=========================================================================

Library   : Image Registration Toolkit (IRTK)
Module    : $Id$
Copyright : Imperial College, Department of Computing
Visual Information Processing (VIP), 2008 onwards
Date      : $Date$
Version   : $Revision$
Changes   : $Author$

Copyright (c) IXICO LIMITED
All rights reserved.
See COPYRIGHT for details

=========================================================================*/

#ifdef HAS_VTK

#include <irtkRegistration.h>

#include <irtkGradientDescentConstrainedOptimizer.h>

#include <irtkResamplingWithPadding.h>

#include <irtkGaussianBlurring.h>

#define HISTORY

#ifdef HAS_TBB

extern tbb::deprecated::concurrent_queue<irtkSimilarityMetric *> sim_queue;

#endif

irtkGreyImage **tmp_mtarget, **tmp_msource;

irtkMultipleImageRegistration::irtkMultipleImageRegistration()
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
        _NumberOfIterations[i] = 20;
        _NumberOfSteps[i]      = 5;
        _LengthOfSteps[i]      = 2;
        _Delta[i]              = 0;
        _Beta[i]               = 1;
    }

    // Default parameters for registration
    _NumberOfLevels     = 1;
    _NumberOfBins       = 64;

    // Default parameters for optimization
    _SimilarityMeasure  = NMI;
    _OptimizationMethod = DownhillDescent;
    _InterpolationMode  = Interpolation_Linear;
    _Epsilon            = 0;
    _Lregu              = 0;

    // Default parameters for debugging
    _DebugFlag = false;

    // Set parameters
    _TargetPadding   = MIN_GREY;

    // Set inputs
    _target = NULL;
    _source = NULL;
    _locator = NULL;

#ifdef HAS_VTK

    // Set landmarks
    _ptarget = NULL;
    _psource = NULL;

#endif

    // Set output
    _transformation = NULL;

    // Set metric
    _metric = NULL;

    // Allocate interpolation object
    _interpolator = NULL;

    // Allocate optimizer object
    _optimizer = NULL;

#ifdef HISTORY
    history = new irtkHistory;
#endif
}

irtkMultipleImageRegistration::~irtkMultipleImageRegistration()
{
#ifdef HISTORY
    delete history;
#endif
}

void irtkMultipleImageRegistration::Initialize()
{
    int i;

    // Check that the t-dimensions are the same for both images.
    // Do not currently support different t-dimensions.
    for (i = 0; i < _numberOfImages; i++) {
        if (_target[i]->GetT() != _source[i]->GetT()) {
            cerr << this->NameOfClass() << "::Initialize() : Images no. " << i << " have different t-dimensions." << endl;
            exit(1);
        }
    }

    _source_x1 = new double[_numberOfImages];
    _source_y1 = new double[_numberOfImages];
    _source_z1 = new double[_numberOfImages];
    _source_x2 = new double[_numberOfImages];
    _source_y2 = new double[_numberOfImages];
    _source_z2 = new double[_numberOfImages];
    _interpolator = new irtkInterpolateImageFunction *[_numberOfImages];

#ifdef HAS_VTK

    if(_ptarget != NULL && _psource != NULL){
        //calculate surface from polydata
        vtkDelaunay2D *delny = vtkDelaunay2D::New();
#if VTK_MAJOR_VERSION >= 6
     delny->SetInputData(_psource);
#else
     delny->SetInput(_psource);
#endif
   
        delny->SetTolerance(0.00001);
        delny->Update();
        _psource->DeepCopy(delny->GetOutput());
        delny->Delete();

        _locator = new irtkLocator;
        _locator->SelectLocatorType(0);
        _locator->SetDataSet(_psource); // data represents the surface  
    }

#endif

}

void irtkMultipleImageRegistration::Initialize(int level)
{
    int i, j, k, t, n;
    double dx, dy, dz, temp;
    irtkGreyPixel target_min, target_max, target_nbins;
    irtkGreyPixel source_min, source_max, source_nbins;

    target_max = MIN_GREY;
    target_min = MAX_GREY;
    source_max = MIN_GREY;
    source_min = MAX_GREY;

    // Allocate memory for temporary images
    tmp_mtarget = new irtkGreyImage*[_numberOfImages];
    tmp_msource = new irtkGreyImage*[_numberOfImages];

    // Allocate memory for metric
    _metric = new irtkSimilarityMetric*[_numberOfImages];

    for (n = 0; n < _numberOfImages; n++) {

        target_max = MIN_GREY;
        target_min = MAX_GREY;
        source_max = MIN_GREY;
        source_min = MAX_GREY;

        // Copy source and target to temp space
        tmp_mtarget[n] = new irtkGreyImage(*_target[n]);
        tmp_msource[n] = new irtkGreyImage(*_source[n]);

        // Swap source and target with temp space copies
        swap(tmp_mtarget[n], _target[n]);
        swap(tmp_msource[n], _source[n]);

        // Blur images if necessary
        if (_TargetBlurring[level] > 0) {
            cout << "Blurring target ... "; cout.flush();
            irtkGaussianBlurringWithPadding<irtkGreyPixel> blurring(_TargetBlurring[level], _TargetPadding);
            blurring.SetInput (_target[n]);
            blurring.SetOutput(_target[n]);
            blurring.Run();
            cout << "done" << endl;
        }

        if (_SourceBlurring[level] > 0) {
            cout << "Blurring source ... "; cout.flush();
            irtkGaussianBlurringWithPadding<irtkGreyPixel> blurring(_SourceBlurring[level], _TargetPadding);
            blurring.SetInput (_source[n]);
            blurring.SetOutput(_source[n]);
            blurring.Run();
            cout << "done" << endl;
        }

        _target[n]->GetPixelSize(&dx, &dy, &dz);
        temp = fabs(_TargetResolution[0][0]-dx) + fabs(_TargetResolution[0][1]-dy) + fabs(_TargetResolution[0][2]-dz);

        if (level > 0 || temp > 0.000001) {
            cout << "Resampling target ... "; cout.flush();
            // Create resampling filter
            irtkResamplingWithPadding<irtkGreyPixel> resample(_TargetResolution[level][0],
                _TargetResolution[level][1],
                _TargetResolution[level][2],
                _TargetPadding);
            resample.SetInput (_target[n]);
            resample.SetOutput(_target[n]);
            resample.Run();
            cout << "done" << endl;
        }

        _source[n]->GetPixelSize(&dx, &dy, &dz);
        temp = fabs(_SourceResolution[0][0]-dx) + fabs(_SourceResolution[0][1]-dy) + fabs(_SourceResolution[0][2]-dz);

        if (level > 0 || temp > 0.000001) {
            cout << "Resampling source ... "; cout.flush();
            // Create resampling filter
            irtkResamplingWithPadding<irtkGreyPixel> resample(_SourceResolution[level][0],
                _SourceResolution[level][1],
                _SourceResolution[level][2], 
                _TargetPadding);

            resample.SetInput (_source[n]);
            resample.SetOutput(_source[n]);
            resample.Run();
            cout << "done" << endl;
        }

        // Find out the min and max values in target image, ignoring padding
        for (t = 0; t < _target[n]->GetT(); t++) {
            for (k = 0; k < _target[n]->GetZ(); k++) {
                for (j = 0; j < _target[n]->GetY(); j++) {
                    for (i = 0; i < _target[n]->GetX(); i++) {
                        if (_target[n]->Get(i, j, k, t) > _TargetPadding) {
                            if (_target[n]->Get(i, j, k, t) > target_max)
                                target_max = _target[n]->Get(i, j, k, t);
                            if (_target[n]->Get(i, j, k, t) < target_min)
                                target_min = _target[n]->Get(i, j, k, t);
                        } else {
                            _target[n]->Put(i, j, k, t, _TargetPadding);
                        }
                    }
                }
            }
        }

        // Find out the min and max values in source image, ignoring padding
        for (t = 0; t < _source[n]->GetT(); t++) {
            for (k = 0; k < _source[n]->GetZ(); k++) {
                for (j = 0; j < _source[n]->GetY(); j++) {
                    for (i = 0; i < _source[n]->GetX(); i++) {
                        if (_source[n]->Get(i, j, k, t) > _TargetPadding) {
                            if (_source[n]->Get(i, j, k, t) > source_max)
                                source_max = _source[n]->Get(i, j, k, t);
                            if (_source[n]->Get(i, j, k, t) < source_min)
                                source_min = _source[n]->Get(i, j, k, t);
                        }else{
                            _source[n]->Put(i, j, k, t, _TargetPadding);
                        }
                    }
                }
            }
        }
        // Check whether dynamic range of data is not to large
        if (target_max - target_min > MAX_GREY) {
            cerr << this->NameOfClass()
                << "::Initialize: Dynamic range of target is too large" << endl;
            exit(1);
        } else {
            for (t = 0; t < _target[n]->GetT(); t++) {
                for (k = 0; k < _target[n]->GetZ(); k++) {
                    for (j = 0; j < _target[n]->GetY(); j++) {
                        for (i = 0; i < _target[n]->GetX(); i++) {
                            if (_target[n]->Get(i, j, k, t) > _TargetPadding) {
                                _target[n]->Put(i, j, k, t, _target[n]->Get(i, j, k, t) - target_min);
                            } else {
                                _target[n]->Put(i, j, k, t, _TargetPadding - 10);
                            }
                        }
                    }
                }
            }
        }

        if ((_SimilarityMeasure == SSD) || (_SimilarityMeasure == CC) ||
            (_SimilarityMeasure == LC)  || (_SimilarityMeasure == K) || (_SimilarityMeasure == ML)) {
                if (source_max - target_min > MAX_GREY) {
                    cerr << this->NameOfClass()
                        << "::Initialize: Dynamic range of source is too large" << endl;
                    exit(1);
                } else {
                    for (t = 0; t < _source[n]->GetT(); t++) {
                        for (k = 0; k < _source[n]->GetZ(); k++) {
                            for (j = 0; j < _source[n]->GetY(); j++) {
                                for (i = 0; i < _source[n]->GetX(); i++) {
                                    if (_source[n]->Get(i, j, k, t) > _TargetPadding) {
                                        _source[n]->Put(i, j, k, t, _source[n]->Get(i, j, k, t) - target_min);
                                    } else {
                                        _source[n]->Put(i, j, k, t, _TargetPadding - 10);
                                    }
                                }
                            }
                        }
                    }
                }
        } else {
            if (source_max - source_min > MAX_GREY) {
                cerr << this->NameOfClass()
                    << "::Initialize: Dynamic range of source is too large" << endl;
                exit(1);
            } else {
                for (t = 0; t < _source[n]->GetT(); t++) {
                    for (k = 0; k < _source[n]->GetZ(); k++) {
                        for (j = 0; j < _source[n]->GetY(); j++) {
                            for (i = 0; i < _source[n]->GetX(); i++) {
                                if (_source[n]->Get(i, j, k, t) > _TargetPadding) {
                                    _source[n]->Put(i, j, k, t, _source[n]->Get(i, j, k, t) - source_min);
                                } else {
                                    _source[n]->Put(i, j, k, t, _TargetPadding - 10);
                                }
                            }
                        }
                    }
                }
            }
        }

        switch (_SimilarityMeasure) {
        case SSD:
            _metric[n] = new irtkSSDSimilarityMetric;
            break;
        case CC:
            // Rescale images by an integer factor if necessary
            _metric[n] = new irtkCrossCorrelationSimilarityMetric;
            break;
        case JE:
            // Rescale images by an integer factor if necessary
            target_nbins = irtkCalculateNumberOfBins(_target[n], _NumberOfBins,
                target_min, target_max);
            source_nbins = irtkCalculateNumberOfBins(_source[n], _NumberOfBins,
                source_min, source_max);
            _metric[n] = new irtkJointEntropySimilarityMetric(target_nbins, source_nbins);
            break;
        case MI:
            // Rescale images by an integer factor if necessary
            target_nbins = irtkCalculateNumberOfBins(_target[n], _NumberOfBins,
                target_min, target_max);
            source_nbins = irtkCalculateNumberOfBins(_source[n], _NumberOfBins,
                source_min, source_max);
            _metric[n] = new irtkMutualInformationSimilarityMetric(target_nbins, source_nbins);
            break;
        case NMI:
            // Rescale images by an integer factor if necessary
            target_nbins = irtkCalculateNumberOfBins(_target[n], _NumberOfBins,
                target_min, target_max);
            source_nbins = irtkCalculateNumberOfBins(_source[n], _NumberOfBins,
                source_min, source_max);
            _metric[n] = new irtkNormalisedMutualInformationSimilarityMetric(target_nbins, source_nbins);
            break;
        case CR_XY:
            // Rescale images by an integer factor if necessary
            target_nbins = irtkCalculateNumberOfBins(_target[n], _NumberOfBins,
                target_min, target_max);
            source_nbins = irtkCalculateNumberOfBins(_source[n], _NumberOfBins,
                source_min, source_max);
            _metric[n] = new irtkCorrelationRatioXYSimilarityMetric(target_nbins, source_nbins);
            break;
        case CR_YX:
            // Rescale images by an integer factor if necessary
            target_nbins = irtkCalculateNumberOfBins(_target[n], _NumberOfBins,
                target_min, target_max);
            source_nbins = irtkCalculateNumberOfBins(_source[n], _NumberOfBins,
                source_min, source_max);
            _metric[n] = new irtkCorrelationRatioYXSimilarityMetric(target_nbins, source_nbins);
            break;
        case LC:
            _metric[n] = new irtkLabelConsistencySimilarityMetric;
            break;
        case K:
            // Rescale images by an integer factor if necessary
            target_nbins = irtkCalculateNumberOfBins(_target[n], _NumberOfBins,
                target_min, target_max);
            source_nbins = irtkCalculateNumberOfBins(_source[n], _NumberOfBins,
                source_min, source_max);
            _metric[n] = new irtkKappaSimilarityMetric(target_nbins, source_nbins);
            break;
        case ML:
            // Rescale images by an integer factor if necessary
            _metric[n] = new irtkMLSimilarityMetric(classification);
            if (_metric[n]==NULL) {
                cerr<<"Please, do not forget to set the ML metric!!!"<<endl;
            }
            break;
        }
        // Print some debugging information
        cout << "Target image (reference) no. " << n << endl;
        _target[n]->Print();
        cout << "Range is from " << target_min << " to " << target_max << endl;

        cout << "Source image (transform) no. " << n << endl;
        _source[n]->Print();
        cout << "Range is from " << source_min << " to " << source_max << endl;

        // Setup the interpolator
        _interpolator[n] = irtkInterpolateImageFunction::New(_InterpolationMode, _source[n]);

        // Setup interpolation for the source image
        _interpolator[n]->SetInput(_source[n]);
        _interpolator[n]->Initialize();

        // Calculate the source image domain in which we can interpolate
        _interpolator[n]->Inside(_source_x1[n], _source_y1[n], _source_z1[n],
            _source_x2[n], _source_y2[n], _source_z2[n]);

    }

    // Setup the optimizer
    switch (_OptimizationMethod) {
    case DownhillDescent:
        _optimizer = new irtkDownhillDescentOptimizer;
        break;
    case GradientDescent:
        _optimizer = new irtkGradientDescentOptimizer;
        break;
    case GradientDescentConstrained:
        _optimizer = new irtkGradientDescentConstrainedOptimizer;
        break;
    case SteepestGradientDescent:
        _optimizer = new irtkSteepestGradientDescentOptimizer;
        break;
    case ConjugateGradientDescent:
        _optimizer = new irtkConjugateGradientDescentOptimizer;
        break;
    default:
        cerr << "Unkown optimizer" << endl;
        exit(1);
    }
    _optimizer->SetTransformation(_transformation);
    _optimizer->SetRegistration(this);

    if(_Lregu > 0){
        _Lregu = _Lregu*_Beta[level];
    }

    // Print initial transformation
    cout << "Initial transformation for level = " << level+1 << endl;;
    _transformation->Print();

}

void irtkMultipleImageRegistration::Finalize()
{
    delete []_source_x1;
    delete []_source_y1;
    delete []_source_z1;
    delete []_source_x2;
    delete []_source_y2;
    delete []_source_z2;
    delete []_interpolator;

#ifdef HAS_VTK
    if(_ptarget != NULL && _psource != NULL){
        delete _locator;
    }

#endif
}

void irtkMultipleImageRegistration::Finalize(int level)
{
    int n;
    // Print final transformation
    cout << "Final transformation for level = " << level+1 << endl;;
    _transformation->Print();
    // Swap source and target back with temp space copies (see Initialize)
    swap(tmp_mtarget, _target);
    swap(tmp_msource, _source);

#ifdef HAS_TBB
    irtkSimilarityMetric *metric;
    while (sim_queue.size() > 0) {
        sim_queue.pop(metric);
        delete metric;
    }
#endif
    for (n = 0; n < _numberOfImages; n++) {
        delete tmp_mtarget[n];
        delete tmp_msource[n];
        delete _interpolator[n];
        delete _metric[n];
    }
    delete []tmp_mtarget;
    delete []tmp_msource;
    delete []_metric;
    delete _optimizer;
}

void irtkMultipleImageRegistration::Run()
{
    int i, j, level;
    char buffer[256];
    double step, epsilon, delta, maxChange;

    // Print debugging information
    this->Debug("irtkMultipleImageRegistration::Run");

    if (_source == NULL) {
        cerr << "Registration::Run: Filter has no source input" << endl;
        exit(1);
    }

    if (_target == NULL) {
        cerr << "Registration::Run: Filter has no target input" << endl;
        exit(1);
    }

    if (_transformation == NULL) {
        cerr << "irtkMultipleImageRegistration::Run: Filter has no transformation output" << endl;
        exit(1);
    }

    // Do the initial set up for all levels
    this->Initialize();

    // Loop over levels
    for (level = _NumberOfLevels-1; level >= 0; level--) {


        // Initial step size
        step = _LengthOfSteps[level];

        // Print resolution level
        cout << "Resolution level no. " << level+1 << " (step sizes ";
        cout << step << " to " << step / pow(2.0, static_cast<double>(_NumberOfSteps[level]-1)) << ")\n";

        // Initial Delta
        delta = _Delta[level];
        cout << "Delta values : " << delta << " to ";
        cout << delta / pow(2.0, static_cast<double>(_NumberOfSteps[level]-1)) << "\n";

#ifdef HISTORY
        history->Clear();
#endif

        // Initialize for this level
        this->Initialize(level);

#ifdef HAS_TBB
        task_scheduler_init init(tbb_no_threads);

        tick_count t_start = tick_count::now();
#endif

        // Run the registration filter at this resolution
        for (i = 0; i < _NumberOfSteps[level]; i++) {
            for (j = 0; j < _NumberOfIterations[level]; j++) {

                cout << "Iteration = " << j + 1 << " (out of " << _NumberOfIterations[level];
                cout << "), step size = " << step << endl;

                // Optimize at lowest level of resolution
                _optimizer->SetStepSize(step);
                _optimizer->SetEpsilon(_Epsilon);
                _optimizer->Run(epsilon, maxChange);

                // Check whether we made any improvement or not
                if (epsilon > _Epsilon && maxChange > delta) {
                    sprintf(buffer, "log_%.3d_%.3d_%.3d.dof", level, i+1, j+1);
                    if (_DebugFlag == true) _transformation->Write(buffer);
                    this->Print();
                } else {
                    sprintf(buffer, "log_%.3d_%.3d_%.3d.dof", level, i+1, j+1);
                    if (_DebugFlag == true) _transformation->Write(buffer);
                    this->Print();
                    break;
                }
            }
            step = step / 2;
            delta = delta / 2.0;
        }

#ifdef HAS_TBB

        tick_count t_end = tick_count::now();
        if (tbb_debug) cout << this->NameOfClass() << " = " << (t_end - t_start).seconds() << " secs." << endl;
        init.terminate();

#endif

        // Do the final cleaning up for this level
        this->Finalize(level);

#ifdef HISTORY
        history->Print();
#endif

    }

    // Do the final cleaning up for all levels
    this->Finalize();
}

double irtkMultipleImageRegistration::LandmarkPenalty()
{

#ifdef HAS_VTK

    int i;
    double d = 0,distance = 0, p[3],q[3];

    if (_ptarget == NULL || _psource == NULL){
        return 0;
    }

    for (i = 0; i < _ptarget->GetNumberOfPoints(); i++) {
        _ptarget->GetPoints()->GetPoint(i,p);
        _transformation->Transform(p[0],p[1],p[2]);
        q[0] = p[0]; q[1] = p[1]; q[2] = p[2]; 
        _locator->FindClosestPoint (q);
        d = sqrt((p[0] - q[0]) * (p[0] - q[0]) +
            (p[1] - q[1]) * (p[1] - q[1]) +
            (p[2] - q[2]) * (p[2] - q[2]));
        distance += d;
    }
    return -(distance/double(_ptarget->GetNumberOfPoints()));

#else

    return 0;

#endif

}

double irtkMultipleImageRegistration::EvaluateGradient(float step, float *dx)
{
    int i;
    double s1, s2, norm, parameterValue;

    for (i = 0; i < _transformation->NumberOfDOFs(); i++) {
        if (_transformation->irtkTransformation::GetStatus(i) == _Active) {
            parameterValue = _transformation->Get(i);
            _transformation->Put(i, parameterValue + step);
            s1 = this->Evaluate();
            _transformation->Put(i, parameterValue - step);
            s2 = this->Evaluate();
            _transformation->Put(i, parameterValue);
            dx[i] = s1 - s2;
        } else {
            dx[i] = 0;
        }
    }

    // Calculate norm of vector
    norm = 0;
    for (i = 0; i < _transformation->NumberOfDOFs(); i++) {
        norm += dx[i] * dx[i];
    }

    // Normalize vector
    norm = sqrt(norm);
    if (norm > 0) {
        for (i = 0; i < _transformation->NumberOfDOFs(); i++) {
            dx[i] /= norm;
        }
    } else {
        for (i = 0; i < _transformation->NumberOfDOFs(); i++) {
            dx[i] = 0;
        }
    }

    return norm;
}

bool irtkMultipleImageRegistration::Read(char *buffer1, char *buffer2, int &level)
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
        _source[0]->GetPixelSize(&dx, &dy, &dz);
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
    if (strstr(buffer1, "No. of steps") != NULL) {
        if (level == -1) {
            for (i = 0; i < MAX_NO_RESOLUTIONS; i++) {
                this->_NumberOfSteps[i] = atoi(buffer2);
            }
        } else {
            this->_NumberOfSteps[level] = atoi(buffer2);
        }
        ok = true;
    }
    if (strstr(buffer1, "Length of steps") != NULL) {
        if (level == -1) {
            for (i = 0; i < MAX_NO_RESOLUTIONS; i++) {
                this->_LengthOfSteps[i] = pow(2.0, double(i)) * atof(buffer2);
            }
        } else {
            this->_LengthOfSteps[level] = atof(buffer2);
        }
        ok = true;
    }
    if (strstr(buffer1, "Epsilon") != NULL) {
        this->_Epsilon = atof(buffer2);
        ok = true;
    }
    if (strstr(buffer1, "Lregulation") != NULL) {
        this->_Lregu = atof(buffer2);
        ok = true;
    }
    if (strstr(buffer1, "Delta") != NULL) {
        if (level == -1) {
            for (i = 0; i < MAX_NO_RESOLUTIONS; i++) {
                this->_Delta[i] = pow(2.0, double(i)) * atof(buffer2);
            }
        } else {
            this->_Delta[level] = atof(buffer2);
        }
        ok = true;
    }
    if (strstr(buffer1, "Beta") != NULL) {
        if (level == -1) {
            for (i = 0; i < MAX_NO_RESOLUTIONS; i++) {
                this->_Beta[i] = pow(2.0, double(i)) * atof(buffer2);
            }
        } else {
            this->_Beta[level] = atof(buffer2);
        }
        ok = true;
    }
    if (strstr(buffer1, "Padding value") != NULL) {
        this->_TargetPadding = atoi(buffer2);
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

    if (strstr(buffer1, "Optimization method") != NULL) {
        if (strstr(buffer2, "DownhillDescent") != NULL) {
            this->_OptimizationMethod = DownhillDescent;
            ok = true;
        } else {
            if (strstr(buffer2, "SteepestGradientDescent") != NULL) {
                this->_OptimizationMethod = SteepestGradientDescent;
                ok = true;
            } else {
                if (strstr(buffer2, "ConjugateGradientDescent") != NULL) {
                    this->_OptimizationMethod = ConjugateGradientDescent;
                    ok = true;
                } else {
                    if (strstr(buffer2, "GradientDescentConstrained") != NULL) {
                        this->_OptimizationMethod = GradientDescentConstrained;
                        ok = true;
                    } else {
                        if (strstr(buffer2, "GradientDescent") != NULL) {
                            this->_OptimizationMethod = GradientDescent;
                            ok = true;
                        }
                    }
                }
            }
        }
    }

    if (ok == false) {
        cerr << "irtkMultipleImageRegistration::Read: Can't parse line " << buffer1 << endl;
        exit(1);
    }

    return ok;
}

void irtkMultipleImageRegistration::Write(ostream &to)
{
    int i;

    to << "\n#\n# Registration parameters\n#\n\n";
    to << "No. of resolution levels          = " << this->_NumberOfLevels << endl;
    to << "No. of bins                       = " << this->_NumberOfBins << endl;
    to << "Epsilon                           = " << this->_Epsilon << endl;
    to << "Lregulation                       = " << this->_Lregu << endl;
    to << "Padding value                     = " << this->_TargetPadding << endl;

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
        cerr << "irtkMultipleImageRegistration::Write: Interpolation mode not supported" << endl;
        exit(1);
    }

    switch (this->_OptimizationMethod) {
    case DownhillDescent:
        to << "Optimization method               = DownhillDescent" << endl;
        break;
    case GradientDescent:
        to << "Optimization method               = GradientDescent" << endl;
        break;
    case SteepestGradientDescent:
        to << "Optimization method               = SteepestGradientDescent" << endl;
        break;
    case ConjugateGradientDescent:
        to << "Optimization method               = ConjugateGradientDescent" << endl;
        break;
    case GradientDescentConstrained:
        to << "Optimization method               = GradientDescentConstrained" << endl;
        break;
    case ClosedForm:
        to << "Optimization method               = ClosedForm" << endl;
        break;
    }

    for (i = 0; i < this->_NumberOfLevels; i++) {
        to << "\n#\n# Registration parameters for resolution level " << i+1 << "\n#\n\n";
        to << "Resolution level                  = " << i+1 << endl;
        to << "Target blurring (in mm)           = " << this->_TargetBlurring[i] << endl;
        to << "Target resolution (in mm)         = " << this->_TargetResolution[i][0] << " " << this->_TargetResolution[i][1] << " " << this->_TargetResolution[i][2] << endl;
        to << "Source blurring (in mm)           = " << this->_SourceBlurring[i] << endl;
        to << "Source resolution (in mm)         = " << this->_SourceResolution[i][0] << " " << this->_SourceResolution[i][1] << " " << this->_SourceResolution[i][2] << endl;
        to << "No. of iterations                 = " << this->_NumberOfIterations[i] << endl;
        to << "No. of steps                      = " << this->_NumberOfSteps[i] << endl;
        to << "Length of steps                   = " << this->_LengthOfSteps[i] << endl;
        to << "Delta                             = " << this->_Delta[i] << endl;
        to << "Beta                              = " << this->_Beta[i] << endl;
    }
}

void irtkMultipleImageRegistration::Read(char *filename)
{
    int level;
    char buffer1[255], *buffer2;

    ifstream from(filename);

    if (!from) {
        cerr << "irtkMultipleImageRegistration::Read: Can't open file " << filename
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

void irtkMultipleImageRegistration::Write(char *filename)
{
    ofstream to(filename);

    if (!to) {
        cerr << "irtkMultipleImageRegistration::Write: Can't open file " << filename;
        exit(1);
    }

    this->Write(to);
}

#endif // HAS_VTK
