/*=========================================================================

  Library   : Image Registration Toolkit (IRTK)
  Module    : $Id$
  Copyright : Imperial College, Department of Computing
              Visual Information Processing (VIP), 2009 onwards
  Date      : $Date$
  Version   : $Revision$
  Changes   : $Author$

=========================================================================*/

#include <irtkRegistration2.h>

#include <irtkResamplingWithPadding.h>

#include <irtkGaussianBlurring.h>

#include <irtkGradientImageFilter.h>

#ifdef WIN32
#include <time.h>
#else
#include <sys/resource.h>
#endif

#define MAX_NO_LINE_ITERATIONS 20

#ifdef HAS_VTK

extern irtkGreyImage **tmp_mtarget, **tmp_msource;

#else

irtkGreyImage **tmp_mtarget, **tmp_msource;

#endif

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

irtkMultipleImageRegistration2::irtkMultipleImageRegistration2()
{
  int i;

  for (i = 0; i < MAX_NO_RESOLUTIONS; i++) {
    // Default parameters for target image
    _TargetBlurring[i]      = 0;

    // Default parameters for source image
    _SourceBlurring[i]      = 0;

    // Default parameters for optimization
    _NumberOfIterations[i] = 40;
    _MinStep[i]            = 0.01;
    _MaxStep[i]            = 1;
  }

  // Default parameters for registration
  _NumberOfLevels     = 1;
  _NumberOfBins       = 64;

  // Default parameters for optimization
  _SimilarityMeasure  = NMI;
  _InterpolationMode  = Interpolation_BSpline;
  _Epsilon            = 0;
  _Lregu              = 0;

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

  // weight
  _weight = NULL;

#ifdef HAS_VTK

  // Set landmarks
  _ptarget = NULL;
  _psource = NULL;
  _pgradient = NULL;

#endif
}

irtkMultipleImageRegistration2::~irtkMultipleImageRegistration2()
{

#ifdef HAS_VTK
    if(_pgradient != NULL){
        _pgradient->Delete();
    }
#endif
}

void irtkMultipleImageRegistration2::Initialize()
{
  // Check that the t-dimensions are the same for both images.
  // Do not currently support different t-dimensions.
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
    _maxDiff = new int[_numberOfImages];
    _interpolator = new irtkInterpolateImageFunction *[_numberOfImages];
    _interpolatorGradient = new irtkInterpolateImageFunction *[_numberOfImages];
}

void irtkMultipleImageRegistration2::Initialize(int level)
{
  int i, j, k, t, n, target_nbins, source_nbins;
  double dx, dy, dz;

  // Allocate memory for temporary images
  tmp_mtarget = new irtkGreyImage*[_numberOfImages];
  tmp_msource = new irtkGreyImage*[_numberOfImages];
  _transformedSource = new irtkGenericImage<double>[_numberOfImages];

  /// Gradient of the original source
  _sourceGradient = new irtkGenericImage<double>[_numberOfImages];

  /// Gradient of the transformed source
  _transformedSourceGradient = new irtkGenericImage<double>[_numberOfImages];

  /// Gradient of the similarity metric
  _similarityGradient = new irtkGenericImage<double>[_numberOfImages];

  // Allocate memory for metric
  _histogram = new irtkHistogram_2D<double>*[_numberOfImages];

  _totalVoxels = 0;

  _totalWeight = 0;

  _level = level;

  for (n = 0; n < _numberOfImages; n++) {

      _histogram[n] = NULL;

      // Copy source and target to temp space
      tmp_mtarget[n] = new irtkGreyImage(*_target[n]);
      tmp_msource[n] = new irtkGreyImage(*_source[n]);

      // Swap source and target with temp space copies
      swap(tmp_mtarget[n], _target[n]);
      swap(tmp_msource[n], _source[n]);

      // Resample image according to level
      if (level > 0 ) {
    	  cout << "Resampling target "<< n <<"... "; cout.flush();
    	  // Guess resolution do not resample z if z is way larger then x y
    	  // For example cardiac SA image
    	  _target[n]->GetPixelSize(&dx, &dy, &dz);
    	  dx = GuessResolution(dx,dy);
    	  dx = pow(2.0,level)*dx;
    	  if(dx > dz) dz = dx;
    	  // Create resampling filter
    	  irtkResamplingWithPadding<irtkGreyPixel> resample(dx, dx, dz, _TargetPadding);
    	  resample.SetInput (_target[n]);
    	  resample.SetOutput(_target[n]);
    	  resample.Run();
    	  cout << "done" << endl;
      }

      // Blur images if necessary
      if (_TargetBlurring[level] > 0) {
    	  cout << "Blurring target "<< n <<"... "; cout.flush();
    	  irtkGaussianBlurringWithPadding<irtkGreyPixel> blurring(_TargetBlurring[level], _TargetPadding);
    	  blurring.SetInput (_target[n]);
    	  blurring.SetOutput(_target[n]);
    	  blurring.Run();
    	  cout << "done" << endl;
      }

      // Resample image according to level do not resample for level 0
      if (level > 0 ) {
    	  cout << "Resampling source "<< n <<"... "; cout.flush();
    	  // Guess resolution do not resample z if z is way larger then x y
    	  // For example cardiac SA image
    	  _source[n]->GetPixelSize(&dx, &dy, &dz);
    	  dx = GuessResolution(dx,dy);
    	  dx = pow(2.0,level)*dx;
    	  if(dx > dz) dz = dx;
    	  // Create resampling filter
    	  irtkResamplingWithPadding<irtkGreyPixel> resample(dx, dx, dz, _SourcePadding);
    	  resample.SetInput (_source[n]);
    	  resample.SetOutput(_source[n]);
    	  resample.Run();
    	  cout << "done" << endl;
      }

      if (_SourceBlurring[level] > 0) {
    	  cout << "Blurring source "<< n <<"... "; cout.flush();
    	  irtkGaussianBlurringWithPadding<irtkGreyPixel> blurring(_SourceBlurring[level], _SourcePadding);
    	  blurring.SetInput (_source[n]);
    	  blurring.SetOutput(_source[n]);
    	  blurring.Run();
    	  cout << "done" << endl;
      }

      //calculate total weight with or without external weights
      if(_weight != NULL) {
          cout << "image " << n << "th's external weight is " << _weight[_level][n] << endl;
          _totalWeight += _target[n]->GetNumberOfVoxels()*_weight[_level][n];
      }

      // without external weights images are weighted according to number of voxels
      _totalVoxels += _target[n]->GetNumberOfVoxels();

      // Find out the min and max values in target image, ignoring padding
      _target_max = MIN_GREY;
      _target_min = MAX_GREY;
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

      // Find out the min and max values in source image
      _source_max = MIN_GREY;
      _source_min = MAX_GREY;
      for (t = 0; t < _source[n]->GetT(); t++) {
          for (k = 0; k < _source[n]->GetZ(); k++) {
              for (j = 0; j < _source[n]->GetY(); j++) {
                  for (i = 0; i < _source[n]->GetX(); i++) {
                      if (_source[n]->Get(i, j, k, t) > _SourcePadding) {
                          if (_source[n]->Get(i, j, k, t) > _source_max)
                              _source_max = _source[n]->Get(i, j, k, t);
                          if (_source[n]->Get(i, j, k, t) < _source_min)
                              _source_min = _source[n]->Get(i, j, k, t);
                      } else {
                          _source[n]->Put(i, j, k, t, _SourcePadding);
                      }
                  }
              }
          }
      }

      // Check whether dynamic range of data is not to large
      if (_target_max - _target_min > MAX_GREY) {
          cerr << this->NameOfClass()
              << "::Initialize: Dynamic range of target is too large" << endl;
          exit(1);
      } else {
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

      // Compute the maximum possible difference across target and source
      _maxDiff[n] = (_target_min - _source_max) * (_target_min - _source_max) > (_target_max - _source_min) * (_target_max - _source_min) ?
          (_target_min - _source_max) * (_target_min - _source_max) : (_target_max - _source_min) * (_target_max - _source_min);

      if (_SimilarityMeasure == SSD) {
          if (_source_max - _target_min > MAX_GREY) {
              cerr << this->NameOfClass()
                  << "::Initialize: Dynamic range of source is too large" << endl;
              exit(1);
          } else {
              for (t = 0; t < _source[n]->GetT(); t++) {
                  for (k = 0; k < _source[n]->GetZ(); k++) {
                      for (j = 0; j < _source[n]->GetY(); j++) {
                          for (i = 0; i < _source[n]->GetX(); i++) {
                              if (_source[n]->Get(i, j, k, t) > _SourcePadding) {
                                  _source[n]->Put(i, j, k, t, _source[n]->Get(i, j, k, t) - _target_min);
                              } else {
                                  _source[n]->Put(i, j, k, t, -1);
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
              for (t = 0; t < _source[n]->GetT(); t++) {
                  for (k = 0; k < _source[n]->GetZ(); k++) {
                      for (j = 0; j < _source[n]->GetY(); j++) {
                          for (i = 0; i < _source[n]->GetX(); i++) {
                              if (_source[n]->Get(i, j, k, t) > _SourcePadding) {
                                  _source[n]->Put(i, j, k, t, _source[n]->Get(i, j, k, t) - _source_min);
                              } else {
                                  _source[n]->Put(i, j, k, t, -1);
                              }
                          }
                      }
                  }
              }
          }
      }

      // Pad target image if necessary
      irtkPadding(*_target[n], _TargetPadding);
      irtkPadding(*_source[n], _SourcePadding);

      // Allocate memory for histogram if necessary
      switch (_SimilarityMeasure) {
      case SSD:
          break;
      case JE:
      case MI:
      case NMI:
          // Rescale images by an integer factor if necessary
          target_nbins = irtkCalculateNumberOfBins(_target[n], _NumberOfBins, _target_min, _target_max);
          source_nbins = irtkCalculateNumberOfBins(_source[n], _NumberOfBins, _source_min, _source_max);
          _histogram[n] = new irtkHistogram_2D<double>(target_nbins, source_nbins);
          break;
      default:
          cerr << this->NameOfClass() << "::Initialize(int): No such metric implemented" << endl;
          exit(1);
      }

      // Compute spatial gradient of source image
      irtkGradientImageFilter<double> gradient(irtkGradientImageFilter<double>::GRADIENT_VECTOR);
      irtkGenericImage<double> tmp = *_source[n];
      gradient.SetInput (&tmp);
      gradient.SetOutput(&_sourceGradient[n]);
      gradient.SetPadding(_SourcePadding);
      gradient.Run();

      // Determine attributes of source image
      irtkImageAttributes attr = _target[n]->GetImageAttributes();

      // Set up gradient of source image (transformed)
      attr._t = 3;
      _transformedSourceGradient[n].Initialize(attr);

      // Set up gradient of similarity metric
      attr._t = 3;
      _similarityGradient[n].Initialize(attr);

      // Setup the interpolator
      _interpolator[n] = irtkInterpolateImageFunction::New(_InterpolationMode, _source[n]);

      // Setup interpolation for the source image
      _interpolator[n]->SetInput(_source[n]);
      _interpolator[n]->Initialize();

      // Calculate the source image domain in which we can interpolate
      _interpolator[n]->Inside(_source_x1[n], _source_y1[n], _source_z1[n],
          _source_x2[n], _source_y2[n], _source_z2[n]);

      // Setup the interpolator
      _interpolatorGradient[n] = irtkInterpolateImageFunction::New(_InterpolationMode, &_sourceGradient[n]);

      // Setup interpolation for the source image
      _interpolatorGradient[n]->SetInput(&_sourceGradient[n]);
      _interpolatorGradient[n]->Initialize();

      // Calculate the source image domain in which we can interpolate
      _interpolatorGradient[n]->Inside(_source_x1[n], _source_y1[n], _source_z1[n],
          _source_x2[n], _source_y2[n], _source_z2[n]);

      // Print some debugging information
      cout << "Target image (reference) "<< n << endl;
      _target[n]->Print();
      cout << "Range is from " << _target_min << " to " << _target_max << endl;

      cout << "Source image (transform) "<< n << endl;
      _source[n]->Print();
      cout << "Range is from " << _source_min << " to " << _source_max << endl;
  }

  // Print initial transformation
  cout << "Initial transformation for level = " << level+1 << endl;;
  _transformation->Print();
}

void irtkMultipleImageRegistration2::Finalize()
{
    delete []_source_x1;
    delete []_source_y1;
    delete []_source_z1;
    delete []_source_x2;
    delete []_source_y2;
    delete []_source_z2;
    delete []_maxDiff;
    delete []_interpolator;
    delete []_interpolatorGradient;
}

void irtkMultipleImageRegistration2::Finalize(int level)
{
  // Print final transformation
  cout << "Final transformation for level = " << level+1 << endl;;
  _transformation->Print();

  for (int n = 0; n < _numberOfImages; n++) {
      // Swap source and target back with temp space copies (see Initialize)
      swap(tmp_mtarget[n], _target[n]);
      swap(tmp_msource[n], _source[n]);

      delete tmp_mtarget[n];
      delete tmp_msource[n];
      delete _interpolator[n];
      delete _interpolatorGradient[n];
      if (_histogram[n] != NULL) {
          delete _histogram[n];
          _histogram[n] = NULL;
      }
  }

  delete []_transformedSource;

  /// Gradient of the original source
  delete []_sourceGradient;

  /// Gradient of the transformed source
  delete []_transformedSourceGradient;

  /// Gradient of the similarity metric
  delete []_similarityGradient;

  delete []_histogram;
  delete []tmp_mtarget;
  delete []tmp_msource;
}

void irtkMultipleImageRegistration2::Update(bool updateGradient)
{
  // Update
  if (updateGradient == true) {
    this->UpdateSourceAndGradient();
  } else {
    this->UpdateSource();
  }
}

void irtkMultipleImageRegistration2::Run()
{
  int i, k;
  char buffer[256];
  double *gradient, delta, step, min_step, max_step, max_length, best_similarity, new_similarity, old_similarity;

  // Print debugging information
  this->Debug("irtkMultipleImageRegistration2::Run");

  if (_source == NULL) {
    cerr << "irtkMultipleImageRegistration2::Run: Filter has no source input" << endl;
    exit(1);
  }

  if (_target == NULL) {
    cerr << "irtkMultipleImageRegistration2::Run: Filter has no target input" << endl;
    exit(1);
  }

  if (_transformation == NULL) {
    cerr << "irtkMultipleImageRegistration2::Run: Filter has no transformation output" << endl;
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
    for(int n = 0; n < _numberOfImages; n++){
        sprintf(buffer, "source_%d_%d.nii.gz", _CurrentLevel, n);
        if (_DebugFlag == true) _source[n]->Write(buffer);
        sprintf(buffer, "target_%d.nii.gz", _CurrentLevel);
        if (_DebugFlag == true) _target[n]->Write(buffer);
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

double irtkMultipleImageRegistration2::EvaluateSSD()
{
  int i, n, k, m;
  double norm, ssd, tmpssd;

  // Print debugging information
  this->Debug("irtkMultipleImageRegistration2::EvaluateSSD");

  m = 0;
  ssd = 0;

  for (k = 0; k < _numberOfImages; k++){

  // Pointer to voxels in images
  short  *ptr2target = _target[k]->GetPointerToVoxels();
  double *ptr2source = _transformedSource[k].GetPointerToVoxels();

  // Initialize metric
  n = 0;
  tmpssd = 0;

  // Compute metric
  for (i = 0; i < _target[k]->GetNumberOfVoxels(); i++) {
    if ((*ptr2target >= 0) && (*ptr2source >= 0)) {
      tmpssd += (*ptr2target - *ptr2source) * (*ptr2target - *ptr2source);
      n++;
      m++;
    }
    ptr2target++;
    ptr2source++;
  }

  // Normalize similarity measure by number of voxels and maximum SSD
  if(_totalWeight == 0)
      norm = 1.0 / ((double)n * (double)_maxDiff[k]);
  else
      norm = _weight[_level][k] / ((double)n * (double)_maxDiff[k]);

  ssd += tmpssd*norm*_target[k]->GetNumberOfVoxels();
  }

  if(_totalWeight == 0)
      ssd = ssd/_totalVoxels;
  else
      ssd = ssd/_totalWeight;

  // Return similarity measure
  if (m > 0) {
    return -ssd;
  } else {
    cerr << "irtkMultipleImageRegistration2::EvaluateSSD: No samples available" << endl;
    return 0;
  }
}

double irtkMultipleImageRegistration2::EvaluateNMI()
{
    int i,n;

    double similarity;
    similarity = 0;

    // Print debugging information
    this->Debug("irtkMultipleImageRegistration2::EvaluateNMI");
    for(n = 0; n < _numberOfImages; n++){
        // Pointer to voxels in images
        short  *ptr2target = _target[n]->GetPointerToVoxels();
        double *ptr2source = _transformedSource[n].GetPointerToVoxels();

        // Initialize metric
        _histogram[n]->Reset();

        // Compute metric
        for (i = 0; i < _target[n]->GetNumberOfVoxels(); i++) {
            if ((*ptr2target >= 0) && (*ptr2source >= 0)) {
                _histogram[n]->Add(*ptr2target, round(*ptr2source));
            }
            ptr2target++;
            ptr2source++;
        }

        // Smooth histogram if appropriate
        _histogram[n]->Smooth();

        // Evaluate similarity measure
        if(_totalWeight == 0)
            similarity += _histogram[n]->NormalizedMutualInformation()
            *_target[n]->GetNumberOfVoxels();
        else
            similarity += _histogram[n]->NormalizedMutualInformation()
            *_target[n]->GetNumberOfVoxels()*_weight[_level][n];
    }
    if(_totalWeight == 0)
        return similarity/_totalVoxels;
    else
        return similarity/_totalWeight;
}

double irtkMultipleImageRegistration2::Evaluate()
{
  double metric;

  IRTK_START_TIMING();

  // Print debugging information
  this->Debug("irtkMultipleImageRegistration2::Evaluate");

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

#ifdef HAS_VTK
  if(_Lregu > 0){
      metric += _Lregu*this->LandMarkPenalty();
  }
#endif

  IRTK_END_TIMING("irtkMultipleImageRegistration2::Evaluate");

  // Evaluate similarity measure
  return metric;
}

void irtkMultipleImageRegistration2::EvaluateGradientSSD()
{
  double ssd;
  int i, j, k, n;

  // Print debugging information
  this->Debug("irtkMultipleImageRegistration2::EvaluateGradient");

  for(n = 0; n < _numberOfImages; n++){
      // Pointer to voxels in images
      short  *ptr2target = _target[n]->GetPointerToVoxels();
      double *ptr2source = _transformedSource[n].GetPointerToVoxels();

      // Compute gradient
      for (k = 0; k < _target[n]->GetZ(); k++) {
          for (j = 0; j < _target[n]->GetY(); j++) {
              for (i = 0; i < _target[n]->GetX(); i++) {
                  if ((*ptr2target >= 0) && (*ptr2source >= 0)) {
                      if(_totalWeight == 0)
                          ssd = 2 * (*ptr2target - *ptr2source) 
                                  / (_maxDiff[n]*_totalVoxels);
                      else
                          ssd = 2 * (*ptr2target - *ptr2source) * _weight[_level][n] 
                                  / (_maxDiff[n]*_totalWeight);
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
}

void irtkMultipleImageRegistration2::EvaluateGradientNMI()
{
  int i, j, k, l, t, r, n;
  double w, je, nmi, tmp, targetEntropyGrad[3], sourceEntropyGrad[3], jointEntropyGrad[3], total,current;

  if(_totalWeight == 0)
      total = _totalVoxels;
  else
      total = _totalWeight;

  for(n = 0; n < _numberOfImages; n++){
  // Allocate new histograms
  irtkHistogram_1D<double> logMarginalXHistogram(_histogram[n]->NumberOfBinsX());
  irtkHistogram_1D<double> logMarginalYHistogram(_histogram[n]->NumberOfBinsY());
  irtkHistogram_2D<double> logJointHistogram(_histogram[n]->NumberOfBinsX(), _histogram[n]->NumberOfBinsY());

  if(_totalWeight == 0)
      current = _target[n]->GetNumberOfVoxels();
  else
      current = _target[n]->GetNumberOfVoxels()* _weight[_level][n];

  // Recompute joint histogram
  for (j = 0; j < _histogram[n]->NumberOfBinsY(); j++) {
    for (i = 0; i < _histogram[n]->NumberOfBinsX(); i++) {

      logJointHistogram.Add(i, j, _histogram[n]->irtkHistogram_2D<double>::operator()(i, j));
    }
  }

  // Smooth joint histogram
  //logJointHistogram.Smooth();
  je  = logJointHistogram.JointEntropy();
  nmi = logJointHistogram.NormalizedMutualInformation();

  // Recompute marginal histogram for X
  for (i = 0; i < _histogram[n]->NumberOfBinsX(); i++) {
      for (j = 0; j < _histogram[n]->NumberOfBinsY(); j++) {
      logMarginalXHistogram.Add(i, _histogram[n]->irtkHistogram_2D<double>::operator()(i, j));
    }
  }
  // Recompute marginal histogram for Y
  for (j = 0; j < _histogram[n]->NumberOfBinsY(); j++) {
    for (i = 0; i < _histogram[n]->NumberOfBinsX(); i++) {
      logMarginalYHistogram.Add(j, _histogram[n]->irtkHistogram_2D<double>::operator()(i, j));
    }
  }

  // Log transform histograms
  logJointHistogram.Log();
  logMarginalXHistogram.Log();
  logMarginalYHistogram.Log();

  // Pointer to voxels in images
  short *ptr2target  = _target[n]->GetPointerToVoxels();
  double *ptr2source = _transformedSource[n].GetPointerToVoxels();

  // Loop over images
  for (k = 0; k < _target[n]->GetZ(); k++) {
    for (j = 0; j < _target[n]->GetY(); j++) {
      for (i = 0; i < _target[n]->GetX(); i++) {

        // This code is based on an idea from Marc Modat for computing the NMI derivative as suggested in his niftyreg package
        if ((*ptr2target >= 0) && (*ptr2source >= 0)) {
          irtkGreyPixel targetValue = *ptr2target;
          irtkRealPixel sourceValue = *ptr2source;

          if ((targetValue < _histogram[n]->NumberOfBinsX()) && (sourceValue < _histogram[n]->NumberOfBinsY())) {
            // The two is added because the image is resample between 2 and bin +2
            // if 64 bins are used the histogram will have 68 bins et the image will be between 2 and 65

            sourceValue = (float)floor((double)sourceValue);

            for (l = 0; l < 3; l++) {
              jointEntropyGrad [l] = 0;
              targetEntropyGrad[l] = 0;
              sourceEntropyGrad[l] = 0;
            }

            for (t = targetValue-1; t < targetValue+2; t++) {
              if ((t >= 0) && (t < _histogram[n]->NumberOfBinsX())) {
                for (r = (int)(sourceValue-1.0); r < (int)(sourceValue+2.0); r++) {
                  if ((r >= 0) && (r < _histogram[n]->NumberOfBinsY())) {
                    w = GetBasisSplineValue((double)t-(double)targetValue) * GetBasisSplineDerivativeValue((double)r-(double)sourceValue);

                    double jointLog  = logJointHistogram(t, r);
                    double targetLog = logMarginalXHistogram(t);
                    double resultLog = logMarginalYHistogram(r);

                    tmp = -w * _transformedSourceGradient[n](i, j, k, 0);
                    jointEntropyGrad[0]  -= tmp * jointLog;
                    targetEntropyGrad[0] -= tmp * targetLog;
                    sourceEntropyGrad[0] -= tmp * resultLog;

                    tmp = -w * _transformedSourceGradient[n](i, j, k, 1);
                    jointEntropyGrad[1]  -= tmp * jointLog;
                    targetEntropyGrad[1] -= tmp * targetLog;
                    sourceEntropyGrad[1] -= tmp * resultLog;

                    tmp = -w * _transformedSourceGradient[n](i, j, k, 2);
                    jointEntropyGrad[2]  -= tmp * jointLog;
                    targetEntropyGrad[2] -= tmp * targetLog;
                    sourceEntropyGrad[2] -= tmp * resultLog;

                  }
                }
              }
            }

            _similarityGradient[n](i, j, k, 0) = current/total*((targetEntropyGrad[0] + sourceEntropyGrad[0] - nmi * jointEntropyGrad[0]) / je);
            _similarityGradient[n](i, j, k, 1) = current/total*((targetEntropyGrad[1] + sourceEntropyGrad[1] - nmi * jointEntropyGrad[1]) / je);
            _similarityGradient[n](i, j, k, 2) = current/total*((targetEntropyGrad[2] + sourceEntropyGrad[2] - nmi * jointEntropyGrad[2]) / je);
          } else {
            _similarityGradient[n](i, j, k, 0) = 0;
            _similarityGradient[n](i, j, k, 1) = 0;
            _similarityGradient[n](i, j, k, 2) = 0;
          }
        }else {
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
}

double irtkMultipleImageRegistration2::LandMarkPenalty ()
{

#ifdef HAS_VTK

    int i;
    double d = 0,distance = 0, p[3],q[3];

    if (_ptarget == NULL || _psource == NULL){
        return 0;
    }

    int tnumber = _ptarget->GetNumberOfPoints();
    int snumber = _psource->GetNumberOfPoints();
    int number = min(tnumber,snumber);

    // pairwise distance
    for (i = 0; i < number; i++) {
        _ptarget->GetPoints()->GetPoint(i,p);
        _transformation->Transform(p[0],p[1],p[2]);
        _psource->GetPoints()->GetPoint(i,q);
        d = (p[0] - q[0]) * (p[0] - q[0]) +
            (p[1] - q[1]) * (p[1] - q[1]) +
            (p[2] - q[2]) * (p[2] - q[2]);
        distance += d;
    }
    return -(distance/double(number));

#else

    return 0;

#endif

}

#ifdef HAS_VTK

void irtkMultipleImageRegistration2::EvaluateGradientLandmark ()
{
    int i;
    double p[3], q[3];

    if (_ptarget == NULL || _psource == NULL){
        return;
    }

    for (i = 0; i < _ptarget->GetNumberOfPoints(); i++) {
        _ptarget->GetPoints()->GetPoint(i,p);
        _transformation->Transform(p[0],p[1],p[2]);
        _psource->GetPoints()->GetPoint(i,q);
        _pgradient->GetPoints()->SetPoint(i,
            -2.0*(p[0] - q[0])/_ptarget->GetNumberOfPoints(),
            -2.0*(p[1] - q[1])/_ptarget->GetNumberOfPoints(),
            -2.0*(p[2] - q[2])/_ptarget->GetNumberOfPoints());
    }

}

#endif

double irtkMultipleImageRegistration2::EvaluateGradient(double *)
{
  int i, j, k;
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

#ifdef HAS_VTK
  if(_Lregu > 0){
      this->EvaluateGradientLandmark();
  }
#endif
  for (int n = 0; n < _numberOfImages; n++){

      // Extract matrix for reorientation of gradient
      irtkMatrix m = _source[n]->GetImageToWorldMatrix();

      // Pointer to voxels in images
      short  *ptr2target = _target[n]->GetPointerToVoxels();
      double *ptr2source = _transformedSource[n].GetPointerToVoxels();

      // Reorient gradient
      for (k = 0; k < _target[n]->GetZ(); k++) {
          for (j = 0; j < _target[n]->GetY(); j++) {
              for (i = 0; i < _target[n]->GetX(); i++) {
                  if ((*ptr2target >= 0) && (*ptr2source >= 0)) {
                      x = m(0, 0) * _similarityGradient[n](i, j, k, 0) 
                          + m(0, 1) * _similarityGradient[n](i, j, k, 1)
                          + m(0, 2) * _similarityGradient[n](i, j, k, 2);
                      y = m(1, 0) * _similarityGradient[n](i, j, k, 0)
                          + m(1, 1) * _similarityGradient[n](i, j, k, 1)
                          + m(1, 2) * _similarityGradient[n](i, j, k, 2);
                      z = m(2, 0) * _similarityGradient[n](i, j, k, 0)
                          + m(2, 1) * _similarityGradient[n](i, j, k, 1)
                          + m(2, 2) * _similarityGradient[n](i, j, k, 2);
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

  // Stop timing
  end = clock();
  cpu_time_used = ((double) (end - start)) / CLOCKS_PER_SEC;
  //cout << "CPU time for irtkMultipleImageRegistration2::EvaluateGradient() = " << cpu_time_used << endl;

  // This function always returns 0
  return 0;
}

bool irtkMultipleImageRegistration2::Read(char *buffer1, char *buffer2, int &level)
{
  int i, ok = false;

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
  if (strstr(buffer1, "Epsilon") != NULL) {
    this->_Epsilon = atof(buffer2);
    ok = true;
  }
  if (strstr(buffer1, "Lregulation") != NULL) {
      this->_Lregu = atof(buffer2);
      ok = true;
  }
  if (strstr(buffer1, "Padding value") != NULL) {
    this->_TargetPadding = atoi(buffer2);
    if(this->_SourcePadding == MIN_GREY){
        this->_SourcePadding = this->_TargetPadding;
    }
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
    cerr << "irtkMultipleImageRegistration2::Read: Can't parse line " << buffer1 << endl;
    exit(1);
  }

  return ok;
}

void irtkMultipleImageRegistration2::Write(ostream &to)
{
  int i;

  to << "\n#\n# Registration parameters\n#\n\n";
  to << "No. of resolution levels          = " << this->_NumberOfLevels << endl;
  to << "No. of bins                       = " << this->_NumberOfBins << endl;
  to << "Epsilon                           = " << this->_Epsilon << endl;
  to << "Lregulation                       = " << this->_Lregu << endl;
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
      cerr << "irtkMultipleImageRegistration2::Write: Interpolation mode not supported" << endl;
      exit(1);
  }

  for (i = 0; i < this->_NumberOfLevels; i++) {
    to << "\n#\n# Registration parameters for resolution level " << i+1 << "\n#\n\n";
    to << "Resolution level                  = " << i+1 << endl;
    to << "Target blurring (in mm)           = " << this->_TargetBlurring[i] << endl;
    to << "Source blurring (in mm)           = " << this->_SourceBlurring[i] << endl;
    to << "No. of iterations                 = " << this->_NumberOfIterations[i] << endl;
    to << "Minimum length of steps           = " << this->_MinStep[i] << endl;
    to << "Maximum length of steps           = " << this->_MaxStep[i] << endl;
  }
}

void irtkMultipleImageRegistration2::Read(char *filename)
{
  int level;
  char buffer1[255], *buffer2;

  ifstream from(filename);

  if (!from) {
    cerr << "irtkMultipleImageRegistration2::Read: Can't open file " << filename
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

void irtkMultipleImageRegistration2::Write(char *filename)
{
  ofstream to(filename);

  if (!to) {
    cerr << "irtkMultipleImageRegistration2::Write: Can't open file " << filename;
    exit(1);
  }

  this->Write(to);
}
