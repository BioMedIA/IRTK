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

#include <irtkResampling.h>

#include <irtkGaussianBlurring.h>

#include <irtkGradientImageFilter.h>

#ifdef WIN32
#include <time.h>
#else
#include <sys/resource.h>
#endif

#define MAX_NO_LINE_ITERATIONS 20

extern irtkRealImage *tmp_target, *tmp_source;

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

irtkImageRegistration2::irtkImageRegistration2()
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
    _MinStep[i]            = 0.01;
    _MaxStep[i]            = 1;
  }

  // Default parameters for registration
  _NumberOfLevels     = 1;
  _NumberOfBins       = 64;

  // Default parameters for optimization
  _SimilarityMeasure  = NMI;
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
}

irtkImageRegistration2::~irtkImageRegistration2()
{}

void irtkImageRegistration2::Initialize()
{
  // Check that the t-dimensions are the same for both images.
  // Do not currently support different t-dimensions.
  if ((_target->GetT() > 1) || (_source->GetT() > 1)) {
    cerr << this->NameOfClass() << "::Initialize() : Not implemented for images with t > 1" << endl;
    exit(1);
  }
}

void irtkImageRegistration2::Initialize(int level)
{
  double dx, dy, dz, temp;
  int i, j, k, t;

  // Copy source and target to temp space
  tmp_target = new irtkRealImage(*_target);
  tmp_source = new irtkRealImage(*_source);

  // Swap source and target with temp space copies
  swap(tmp_target, _target);
  swap(tmp_source, _source);

  // Blur images if necessary
  if (_TargetBlurring[level] > 0) {
      cout << "Blurring target ... "; cout.flush();
      irtkGaussianBlurring<irtkRealPixel> blurring(_TargetBlurring[level]);
      blurring.SetInput (_target);
      blurring.SetOutput(_target);
      blurring.Run();
      cout << "done" << endl;
  }

  if (_SourceBlurring[level] > 0) {
      cout << "Blurring source ... "; cout.flush();
      irtkGaussianBlurring<irtkRealPixel> blurring(_SourceBlurring[level]);
      blurring.SetInput (_source);
      blurring.SetOutput(_source);
      blurring.Run();
      cout << "done" << endl;
  }

  _target->GetPixelSize(&dx, &dy, &dz);
  temp = fabs(_TargetResolution[level][0]-dx) 
      + fabs(_TargetResolution[level][1]-dy) 
      + fabs(_TargetResolution[level][2]-dz);

  if (level > 0 || temp > 0.000001) {
    cout << "Resampling target ... "; cout.flush();
    // Create resampling filter
    irtkResampling<irtkRealPixel> resample(_TargetResolution[level][0],
        _TargetResolution[level][1],
        _TargetResolution[level][2]);
    resample.SetInput (_target);
    resample.SetOutput(_target);
    if (_target->GetZ() == 1) {
      resample.SetInterpolator(new irtkLinearInterpolateImageFunction2D);
    } else {
      resample.SetInterpolator(new irtkLinearInterpolateImageFunction);
    }
    resample.Run();
    cout << "done" << endl;
  }

  _source->GetPixelSize(&dx, &dy, &dz);
  temp = fabs(_SourceResolution[level][0]-dx) 
      + fabs(_SourceResolution[level][1]-dy) 
      + fabs(_SourceResolution[level][2]-dz);

  if (level > 0 || temp > 0.000001) {
    cout << "Resampling source ... "; cout.flush();
    // Create resampling filter
    irtkResampling<irtkRealPixel> resample(_SourceResolution[level][0],
        _SourceResolution[level][1],
        _SourceResolution[level][2]);
    resample.SetInput (_source);
    resample.SetOutput(_source);
    if (_source->GetZ() == 1) {
	  resample.SetInterpolator(new irtkLinearInterpolateImageFunction2D);
	} else {
	  resample.SetInterpolator(new irtkLinearInterpolateImageFunction);
	}
    resample.Run();
    cout << "done" << endl;
  }

  _target->GetMinMaxAsDouble(&_target_min, &_target_max);
  _source->GetMinMaxAsDouble(&_source_min, &_source_max);

  // Compute the maximum possible difference across target and source
  _maxDiff = (_target_min - _source_max) * (_target_min - _source_max) > (_target_max - _source_min) * (_target_max - _source_min) ?
             (_target_min - _source_max) * (_target_min - _source_max) : (_target_max - _source_min) * (_target_max - _source_min);

  // Initialise distance mask
  _distanceMask.Initialize(_target->GetImageAttributes());

  //Compute distance mask image if necessary
  irtkPadding(*_target, _TargetPadding, &_distanceMask);

  // Allocate memory for histogram if necessary
  switch (_SimilarityMeasure) {
    case SSD:
    case NGD:
    case NGP:
    case NGS:
      break;
    case JE:
    case MI:
    case NMI:
      //Create histogram
      _histogram = new irtkHistogram_2D<double>(_NumberOfBins, _NumberOfBins);
      break;
    default:
      cerr << this->NameOfClass() << "::Initialize(int): No such metric implemented" << endl;
      exit(1);
  }

  // Compute spatial gradient of source image
  irtkGradientImageFilter<double> gradient(irtkGradientImageFilter<double>::GRADIENT_VECTOR);
  irtkGenericImage<double> tmp = *_source;
  gradient.SetInput (&tmp);
  gradient.SetOutput(&_sourceGradient);
  gradient.SetPadding(_SourcePadding);
  gradient.Run();

  // Set up gradient of source image (transformed)
  irtkImageAttributes attr = _target->GetImageAttributes();
  attr._t = 3;
  _transformedSourceGradient.Initialize(attr);

  // Set up gradient of similarity metric
  attr._t = 3;
  _similarityGradient.Initialize(attr);

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

  // Print some debugging information
  cout << "Target image (reference)" << endl;
  _target->Print();
  cout << "Range is from " << _target_min << " to " << _target_max << endl;

  cout << "Source image (transform)" << endl;
  _source->Print();
  cout << "Range is from " << _source_min << " to " << _source_max << endl;

  // Print initial transformation
  cout << "Initial transformation for level = " << level+1 << endl;;
  _transformation->Print();
}

void irtkImageRegistration2::Finalize()
{}

void irtkImageRegistration2::Finalize(int level)
{
  // Print final transformation
  cout << "Final transformation for level = " << level+1 << endl;;
  _transformation->Print();

  // Swap source and target back with temp space copies (see Initialize)
  swap(tmp_target, _target);
  swap(tmp_source, _source);

  delete tmp_target;
  delete tmp_source;
  delete _interpolator;
  delete _interpolatorGradient;
  if (_histogram != NULL) {
    delete _histogram;
    _histogram = NULL;
  }
}

void irtkImageRegistration2::Update(bool updateGradient)
{
  // Update
  if (updateGradient == true) {
    this->UpdateSourceAndGradient();
  } else {
    this->UpdateSource();
  }
}

void irtkImageRegistration2::Run()
{
  int i, k;
  char buffer[256];
  double *gradient, delta, step, min_step, max_step, max_length, best_similarity, new_similarity, old_similarity;

  // Print debugging information
  this->Debug("irtkImageRegistration2::Run");

  if (_source == NULL) {
    cerr << "irtkImageRegistration2::Run: Filter has no source input" << endl;
    exit(1);
  }

  if (_target == NULL) {
    cerr << "irtkImageRegistration2::Run: Filter has no target input" << endl;
    exit(1);
  }

  if (_transformation == NULL) {
    cerr << "irtkImageRegistration2::Run: Filter has no transformation output" << endl;
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
    if (_DebugFlag == true) _source->Write(buffer);
    sprintf(buffer, "target_%d.nii.gz", _CurrentLevel);
    if (_DebugFlag == true) _target->Write(buffer);

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

double irtkImageRegistration2::EvaluateSSD()
{
  int i, n;
  double norm, ssd;

  // Print debugging information
  this->Debug("irtkImageRegistration2::EvaluateSSD");

  // Pointer to voxels in images
  double *ptr2target  = _target->GetPointerToVoxels();
  double *ptr2source  = _transformedSource.GetPointerToVoxels();
  short  *ptr2mask    = _distanceMask.GetPointerToVoxels();

  // Initialize metric
  n = 0;
  ssd = 0;

  // Compute metric
  for (i = 0; i < _target->GetNumberOfVoxels(); i++) {
	if ((*ptr2mask == 0) && (*ptr2source > _SourcePadding)) {
      ssd += (*ptr2target - *ptr2source) * (*ptr2target - *ptr2source);
      n++;
    }
    ptr2target++;
    ptr2source++;
    ptr2mask++;
  }

  // Normalize similarity measure by number of voxels and maximum SSD
  norm = 1.0 / ((double)n * (double)_maxDiff);

  // Return similarity measure
  if (norm > 0) {
    return -ssd * norm;
  } else {
    cerr << "irtkImageRegistration2::EvaluateSSD: No samples available" << endl;
    return 0;
  }
}

double irtkImageRegistration2::EvaluateNMI()
{
  int i, targetBinValue, sourceBinValue;

  // Print debugging information
  this->Debug("irtkImageRegistration2::EvaluateNMI");

  // Pointer to voxels in images
  double *ptr2target = _target->GetPointerToVoxels();
  double *ptr2source = _transformedSource.GetPointerToVoxels();
  short  *ptr2mask   = _distanceMask.GetPointerToVoxels();

  // Initialize metric
  _histogram->Reset();

  // Compute metric
  for (i = 0; i < _target->GetNumberOfVoxels(); i++) {
	if ((*ptr2mask == 0) && (*ptr2source > _SourcePadding)) {
	  targetBinValue = irtkGetBinIndex(*ptr2target, _target_min, _target_max, _histogram->NumberOfBinsX());
	  sourceBinValue = irtkGetBinIndex(*ptr2source, _source_min, _source_max, _histogram->NumberOfBinsY());
      _histogram->Add(targetBinValue, sourceBinValue);
    }
    ptr2target++;
    ptr2source++;
    ptr2mask++;
  }

  // Smooth histogram if appropriate
  _histogram->Smooth();

  // Evaluate similarity measure
  return _histogram->NormalizedMutualInformation();
}

double irtkImageRegistration2::Evaluate()
{
  IRTK_START_TIMING();

  this->Debug("irtkImageRegistration2::Evaluate");

  double metric;
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

  IRTK_END_TIMING("irtkImageRegistration2::Evaluate");
  return metric;
}

void irtkImageRegistration2::EvaluateGradientSSD()
{
  double ssd;
  int i, j, k;

  // Print debugging information
  this->Debug("irtkImageRegistration2::EvaluateGradient");

  // Pointer to voxels in images
  double *ptr2target = _target->GetPointerToVoxels();
  double *ptr2source = _transformedSource.GetPointerToVoxels();
  short  *ptr2mask   = _distanceMask.GetPointerToVoxels();

  // Compute gradient
  for (k = 0; k < _target->GetZ(); k++) {
    for (j = 0; j < _target->GetY(); j++) {
      for (i = 0; i < _target->GetX(); i++) {
    	if ((*ptr2mask == 0) && (*ptr2source > _SourcePadding)) {
          ssd = 2 * (*ptr2target - *ptr2source) / _maxDiff;
          _similarityGradient(i, j, k, 0) = ssd * _transformedSourceGradient(i, j, k, 0);
          _similarityGradient(i, j, k, 1) = ssd * _transformedSourceGradient(i, j, k, 1);
          _similarityGradient(i, j, k, 2) = ssd * _transformedSourceGradient(i, j, k, 2);
        } else {
          _similarityGradient(i, j, k, 0) = 0;
          _similarityGradient(i, j, k, 1) = 0;
          _similarityGradient(i, j, k, 2) = 0;
        }
        ptr2target++;
        ptr2source++;
        ptr2mask++;
      }
    }
  }
}

void irtkImageRegistration2::EvaluateGradientNMI()
{
  int i, j, k, l, t, r;
  double w, je, nmi, tmp, targetEntropyGrad[3], sourceEntropyGrad[3], jointEntropyGrad[3];

  //Compute constant values
  je  = _histogram->JointEntropy();
  nmi = _histogram->NormalizedMutualInformation();

  // Allocate new histograms
  irtkHistogram_1D<double> logMarginalXHistogram(_histogram->NumberOfBinsX());
  irtkHistogram_1D<double> logMarginalYHistogram(_histogram->NumberOfBinsY());
  irtkHistogram_2D<double> logJointHistogram(_histogram->NumberOfBinsX(), _histogram->NumberOfBinsY());

  // Recompute joint histogram
  for (j = 0; j < _histogram->NumberOfBinsY(); j++) {
    for (i = 0; i < _histogram->NumberOfBinsX(); i++) {
      logJointHistogram.Add(i, j, (double)_histogram->irtkHistogram_2D<double>::operator()(i, j));
    }
  }

  // Recompute marginal histogram for X
  for (i = 0; i < _histogram->NumberOfBinsX(); i++) {
    for (j = 0; j < _histogram->NumberOfBinsY(); j++) {
      logMarginalXHistogram.Add(i, (double)_histogram->irtkHistogram_2D<double>::operator()(i, j));
    }
  }
  // Recompute marginal histogram for Y
  for (j = 0; j < _histogram->NumberOfBinsY(); j++) {
    for (i = 0; i < _histogram->NumberOfBinsX(); i++) {
      logMarginalYHistogram.Add(j, (double)_histogram->irtkHistogram_2D<double>::operator()(i, j));
    }
  }

  // Log transform histograms
  logJointHistogram.Log();
  logMarginalXHistogram.Log();
  logMarginalYHistogram.Log();

  // Pointer to voxels in images
  double *ptr2target = _target->GetPointerToVoxels();
  double *ptr2source = _transformedSource.GetPointerToVoxels();
  short  *ptr2mask   = _distanceMask.GetPointerToVoxels();

  // Loop over images
  for (k = 0; k < _target->GetZ(); k++) {
    for (j = 0; j < _target->GetY(); j++) {
      for (i = 0; i < _target->GetX(); i++) {

        // This code is based on an idea from Marc Modat for computing the NMI derivative as suggested in his niftyreg package
    	if ((*ptr2mask == 0) && (*ptr2source > _SourcePadding)) {
    	  int targetBinValue = irtkGetBinIndex(*ptr2target, _target_min, _target_max, _histogram->NumberOfBinsX());
    	  int sourceBinValue = irtkGetBinIndex(*ptr2source, _source_min, _source_max, _histogram->NumberOfBinsY());

          for (l = 0; l < 3; l++) {
            jointEntropyGrad [l] = 0;
            targetEntropyGrad[l] = 0;
            sourceEntropyGrad[l] = 0;
          }

          for (t = targetBinValue-1; t <= targetBinValue+1; t++) {
            if ((t >= 0) && (t < _histogram->NumberOfBinsX())) {
              for (r = sourceBinValue-1; r <= sourceBinValue+1; r++) {
                if ((r >= 0) && (r < _histogram->NumberOfBinsY())) {
                  w = GetBasisSplineValue((double)t-(double)targetBinValue) * GetBasisSplineDerivativeValue((double)r-(double)sourceBinValue);

                  double jointLog  = logJointHistogram(t, r);
                  double targetLog = logMarginalXHistogram(t);
                  double resultLog = logMarginalYHistogram(r);

                  tmp = -w * _transformedSourceGradient(i, j, k, 0);
                  jointEntropyGrad[0]  -= tmp * jointLog;
                  targetEntropyGrad[0] -= tmp * targetLog;
                  sourceEntropyGrad[0] -= tmp * resultLog;

                  tmp = -w * _transformedSourceGradient(i, j, k, 1);
                  jointEntropyGrad[1]  -= tmp * jointLog;
                  targetEntropyGrad[1] -= tmp * targetLog;
                  sourceEntropyGrad[1] -= tmp * resultLog;

                  tmp = -w * _transformedSourceGradient(i, j, k, 2);
                  jointEntropyGrad[2]  -= tmp * jointLog;
                  targetEntropyGrad[2] -= tmp * targetLog;
                  sourceEntropyGrad[2] -= tmp * resultLog;
                }
              }
            }
          }

          _similarityGradient(i, j, k, 0) = ((targetEntropyGrad[0] + sourceEntropyGrad[0] - nmi * jointEntropyGrad[0]) / je);
          _similarityGradient(i, j, k, 1) = ((targetEntropyGrad[1] + sourceEntropyGrad[1] - nmi * jointEntropyGrad[1]) / je);
          _similarityGradient(i, j, k, 2) = ((targetEntropyGrad[2] + sourceEntropyGrad[2] - nmi * jointEntropyGrad[2]) / je);
        } else {
          _similarityGradient(i, j, k, 0) = 0;
          _similarityGradient(i, j, k, 1) = 0;
          _similarityGradient(i, j, k, 2) = 0;
        }
        ptr2target++;
        ptr2source++;
        ptr2mask++;
      }
    }
  }
}

double irtkImageRegistration2::EvaluateGradient(double *)
{
  int i, j, k;
  double x, y, z;

  IRTK_START_TIMING();

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

  // Pointer to voxels in mask
  double *ptr2source = _transformedSource.GetPointerToVoxels();
  short  *ptr2mask   = _distanceMask.GetPointerToVoxels();

  // Reorient gradient
  for (k = 0; k < _target->GetZ(); k++) {
    for (j = 0; j < _target->GetY(); j++) {
      for (i = 0; i < _target->GetX(); i++) {
      	if ((*ptr2mask == 0) && (*ptr2source > _SourcePadding)) {
      	  x = m(0, 0) * _similarityGradient(i, j, k, 0) + m(0, 1) * _similarityGradient(i, j, k, 1) + m(0, 2) * _similarityGradient(i, j, k, 2);
      	  y = m(1, 0) * _similarityGradient(i, j, k, 0) + m(1, 1) * _similarityGradient(i, j, k, 1) + m(1, 2) * _similarityGradient(i, j, k, 2);
      	  z = m(2, 0) * _similarityGradient(i, j, k, 0) + m(2, 1) * _similarityGradient(i, j, k, 1) + m(2, 2) * _similarityGradient(i, j, k, 2);
      	  _similarityGradient(i, j, k, 0) = x;
      	  _similarityGradient(i, j, k, 1) = y;
      	  _similarityGradient(i, j, k, 2) = z;
        }
        ptr2mask++;
        ptr2source++;
      }
    }
  }

  IRTK_END_TIMING("irtkImageRegistration2::EvaluateGradient()");

  // This function always returns 0
  return 0;
}

bool irtkImageRegistration2::Read(char *buffer1, char *buffer2, int &level)
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
    _target->GetPixelSize(&dx, &dy, &dz);
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
  if (strstr(buffer1, "Epsilon") != NULL) {
    this->_Epsilon = atof(buffer2);
    ok = true;
  }
  if (strstr(buffer1, "Padding value") != NULL) {
    this->_TargetPadding = atof(buffer2);
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
                      } else {
                        if (strstr(buffer2, "NGD") != NULL) {
                          this->_SimilarityMeasure = NGD;
                          ok = true;
                        } else {
                          if (strstr(buffer2, "NGP") != NULL) {
                            this->_SimilarityMeasure = NGP;
                            ok = true;
                          } else {
                            if (strstr(buffer2, "NGS") != NULL) {
                              this->_SimilarityMeasure = NGS;
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
    cerr << "irtkImageRegistration2::Read: Can't parse line " << buffer1 << endl;
    exit(1);
  }

  return ok;
}

void irtkImageRegistration2::Write(ostream &to)
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
    case NGD:
      to << "Similarity measure                = NGD" << endl;
      break;
    case NGP:
      to << "Similarity measure                = NGP" << endl;
      break;
    case NGS:
      to << "Similarity measure                = NGS" << endl;
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
      cerr << "irtkImageRegistration2::Write: Interpolation mode not supported" << endl;
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
  }
}

void irtkImageRegistration2::Read(char *filename)
{
  int level;
  char buffer1[255], *buffer2;

  ifstream from(filename);

  if (!from) {
    cerr << "irtkImageRegistration2::Read: Can't open file " << filename
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

void irtkImageRegistration2::Write(char *filename)
{
  ofstream to(filename);

  if (!to) {
    cerr << "irtkImageRegistration2::Write: Can't open file " << filename;
    exit(1);
  }

  this->Write(to);
}
