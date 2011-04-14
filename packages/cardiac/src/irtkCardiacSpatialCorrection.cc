/*=========================================================================

  Library   : Image Registration Toolkit (IRTK)
  Module    : $Id: irtkCardiacSpatialCorrection.cc 247 2010-10-28 22:07:41Z dr $
  Copyright : Imperial College, Department of Computing
              Visual Information Processing (VIP), 2008 onwards
  Date      : $Date: 2010-10-28 23:07:41 +0100 (ËÄ, 28 Ê®ÔÂ 2010) $
  Version   : $Revision: 247 $
  Changes   : $Author: dr $

=========================================================================*/

#include <irtkCardiac.h>

#include <irtkGradientDescentConstrainedOptimizer.h>

#include <irtkResamplingWithPadding.h>

#include <irtkGaussianBlurring.h>

#include <irtkHomogeneousTransformationIterator.h>

#ifdef HAS_TBB

concurrent_queue<irtkSimilarityMetric *> queue;

#endif

irtkGreyImage *tmp_satarget, *tmp_3dsource;

irtkCardiacSpatialCorrection::irtkCardiacSpatialCorrection()
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
  }

  // Default parameters for registration
  _NumberOfLevels     = 1;
  _NumberOfBins       = 64;

  // Default parameters for optimization
  _SimilarityMeasure  = NMI;
  _InterpolationMode  = Interpolation_Linear;
  _Epsilon            = 0;
  _IntLambda          = 0;
  _ConstrainLambda    = 0;

  // Default parameters for debugging
  _DebugFlag = false;

  // Set parameters
  _TargetPadding   = MIN_GREY;

  // Set inputs
  _target = NULL;
  _source = NULL;

  // Set output
  _transformation = NULL;
  _origintransformation = NULL;

  // Set metric
  _metric = NULL;
  _intmetric = NULL;

  // Allocate interpolation object
  _interpolator = NULL;

#ifdef HISTORY
  history = new irtkHistory;
#endif
}

irtkCardiacSpatialCorrection::~irtkCardiacSpatialCorrection()
{
#ifdef HISTORY
  delete history;
#endif
}

void irtkCardiacSpatialCorrection::Initialize()
{
  // Check that the t-dimensions are the same for both images.
  // Do not currently support different t-dimensions.
  if (_target->GetT() != _source->GetT()) {
    cerr << this->NameOfClass() << "::Initialize() : Images have different t-dimensions." << endl;
    exit(1);
  }
}

void irtkCardiacSpatialCorrection::Initialize(int level)
{
  int i, j, k, t;
  double dx, dy, dz, temp;
  irtkGreyPixel target_min, target_max, target_nbins;
  irtkGreyPixel source_min, source_max, source_nbins;

  // Copy source and target to temp space
  tmp_satarget = new irtkGreyImage(*_target);
  tmp_3dsource = new irtkGreyImage(*_source);

  // Swap source and target with temp space copies
  swap(tmp_satarget, _target);
  swap(tmp_3dsource, _source);

  // Blur images if necessary
  if (_TargetBlurring[level] > 0) {
    cout << "Blurring target ... "; cout.flush();
    irtkGaussianBlurringWithPadding<irtkGreyPixel> blurring(_TargetBlurring[level], _TargetPadding);
    blurring.SetInput (_target);
    blurring.SetOutput(_target);
    blurring.Run();
    cout << "done" << endl;
  }

  if (_SourceBlurring[level] > 0) {
    cout << "Blurring source ... "; cout.flush();
    irtkGaussianBlurring<irtkGreyPixel> blurring(_SourceBlurring[level]);
    blurring.SetInput (_source);
    blurring.SetOutput(_source);
    blurring.Run();
    cout << "done" << endl;
  }

  _target->GetPixelSize(&dx, &dy, &dz);
  temp = fabs(_TargetResolution[0][0]-dx) + fabs(_TargetResolution[0][1]-dy) + fabs(_TargetResolution[0][2]-dz);

  if (level > 0 || temp > 0.000001) {
    cout << "Resampling target ... "; cout.flush();
    // Create resampling filter
    irtkResamplingWithPadding<irtkGreyPixel> resample(_TargetResolution[level][0],
        _TargetResolution[level][1],
        _TargetResolution[level][2],
        _TargetPadding);
    resample.SetInput (_target);
    resample.SetOutput(_target);
    resample.Run();
    cout << "done" << endl;
  }

  _source->GetPixelSize(&dx, &dy, &dz);
  temp = fabs(_SourceResolution[0][0]-dx) + fabs(_SourceResolution[0][1]-dy) + fabs(_SourceResolution[0][2]-dz);

  if (level > 0 || temp > 0.000001) {
    cout << "Resampling source ... "; cout.flush();
    // Create resampling filter
    irtkResamplingWithPadding<irtkGreyPixel> resample(_SourceResolution[level][0],
        _SourceResolution[level][1],
        _SourceResolution[level][2], MIN_GREY);

    resample.SetInput (_source);
    resample.SetOutput(_source);
    resample.Run();
    cout << "done" << endl;
  }

  // Find out the min and max values in target image, ignoring padding
  target_max = MIN_GREY;
  target_min = MAX_GREY;
  for (t = 0; t < _target->GetT(); t++) {
    for (k = 0; k < _target->GetZ(); k++) {
      for (j = 0; j < _target->GetY(); j++) {
        for (i = 0; i < _target->GetX(); i++) {
          if (_target->Get(i, j, k, t) > _TargetPadding) {
            if (_target->Get(i, j, k, t) > target_max)
              target_max = _target->Get(i, j, k, t);
            if (_target->Get(i, j, k, t) < target_min)
              target_min = _target->Get(i, j, k, t);
          } else {
            _target->Put(i, j, k, t, _TargetPadding);
          }
        }
      }
    }
  }

  // Find out the min and max values in source image, ignoring padding
  source_max = MIN_GREY;
  source_min = MAX_GREY;
  for (t = 0; t < _source->GetT(); t++) {
    for (k = 0; k < _source->GetZ(); k++) {
      for (j = 0; j < _source->GetY(); j++) {
        for (i = 0; i < _source->GetX(); i++) {
          if (_source->Get(i, j, k, t) > source_max)
            source_max = _source->Get(i, j, k, t);
          if (_source->Get(i, j, k, t) < source_min)
            source_min = _source->Get(i, j, k, t);
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
    for (t = 0; t < _target->GetT(); t++) {
      for (k = 0; k < _target->GetZ(); k++) {
        for (j = 0; j < _target->GetY(); j++) {
          for (i = 0; i < _target->GetX(); i++) {
            if (_target->Get(i, j, k, t) > _TargetPadding) {
              _target->Put(i, j, k, t, _target->Get(i, j, k, t) - target_min);
            } else {
              _target->Put(i, j, k, t, -1);
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
      for (t = 0; t < _source->GetT(); t++) {
        for (k = 0; k < _source->GetZ(); k++) {
          for (j = 0; j < _source->GetY(); j++) {
            for (i = 0; i < _source->GetX(); i++) {
              _source->Put(i, j, k, t, _source->Get(i, j, k, t) - target_min);
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
      for (t = 0; t < _source->GetT(); t++) {
        for (k = 0; k < _source->GetZ(); k++) {
          for (j = 0; j < _source->GetY(); j++) {
            for (i = 0; i < _source->GetX(); i++) {
              _source->Put(i, j, k, t, _source->Get(i, j, k, t) - source_min);
            }
          }
        }
      }
    }
  }

  // Pad target image if necessary
  irtkPadding(*_target, _TargetPadding);

  // Allocate memory for metric
  switch (_SimilarityMeasure) {
  case SSD:
    _metric = new irtkSSDSimilarityMetric;
	_intmetric = new irtkSSDSimilarityMetric;
    break;
  case CC:
    // Rescale images by an integer factor if necessary
    _metric = new irtkCrossCorrelationSimilarityMetric;
	_intmetric = new irtkCrossCorrelationSimilarityMetric;
    break;
  case JE:
    // Rescale images by an integer factor if necessary
    target_nbins = irtkCalculateNumberOfBins(_target, _NumberOfBins,
                   target_min, target_max);
    source_nbins = irtkCalculateNumberOfBins(_source, _NumberOfBins,
                   source_min, source_max);
    _metric = new irtkJointEntropySimilarityMetric(target_nbins, source_nbins);
	_intmetric = new irtkJointEntropySimilarityMetric(target_nbins, target_nbins);
    break;
  case MI:
    // Rescale images by an integer factor if necessary
    target_nbins = irtkCalculateNumberOfBins(_target, _NumberOfBins,
                   target_min, target_max);
    source_nbins = irtkCalculateNumberOfBins(_source, _NumberOfBins,
                   source_min, source_max);
    _metric = new irtkMutualInformationSimilarityMetric(target_nbins, source_nbins);
	_intmetric = new irtkMutualInformationSimilarityMetric(target_nbins, target_nbins);
    break;
  case NMI:
    // Rescale images by an integer factor if necessary
    target_nbins = irtkCalculateNumberOfBins(_target, _NumberOfBins,
                   target_min, target_max);
    source_nbins = irtkCalculateNumberOfBins(_source, _NumberOfBins,
                   source_min, source_max);
    _metric = new irtkNormalisedMutualInformationSimilarityMetric(target_nbins, source_nbins);
	_intmetric = new irtkNormalisedMutualInformationSimilarityMetric(target_nbins, target_nbins);
    break;
  default:
  	  cerr<<"Can not recognize similarity metric type!!!"<<endl;
  	  exit(1);
  	break;
  }

  // Setup the interpolator
  _interpolator = irtkInterpolateImageFunction::New(_InterpolationMode, _source);

  // Setup interpolation for the source image
  _interpolator->SetInput(_source);
  _interpolator->Initialize();

  // Calculate the source image domain in which we can interpolate
  _interpolator->Inside(_source_x1, _source_y1, _source_z1,
                        _source_x2, _source_y2, _source_z2);

  // Print some debugging information
  cout << "Target image (reference)" << endl;
  _target->Print();
  cout << "Range is from " << target_min << " to " << target_max << endl;

  cout << "Source image (transform)" << endl;
  _source->Print();
  cout << "Range is from " << source_min << " to " << source_max << endl;

  // Print initial transformation
  cout << "Initial transformation for level = " << level+1 << endl;
  for(k = 0; k < _target->GetZ(); k++){
   _transformation[k]->Print();
   cout << endl;
  }
}

void irtkCardiacSpatialCorrection::Finalize()
{}

void irtkCardiacSpatialCorrection::Finalize(int level)
{
	int k;
  // Print final transformation
    cout << "Final transformation for level = " << level+1 << endl;
	for(k = 0; k < _target->GetZ(); k++){
		_transformation[k]->Print();
		cout << endl;
	}
    // Swap source and target back with temp space copies (see Initialize)
    swap(tmp_satarget, _target);
    swap(tmp_3dsource, _source);

#ifdef HAS_TBB
  irtkSimilarityMetric *metric;
  while (queue.size() > 0) {
    queue.pop(metric);
    delete metric;
  }
#endif

  delete tmp_satarget;
  delete tmp_3dsource;
  delete _metric;
  delete _intmetric;
  delete _interpolator;
}

void irtkCardiacSpatialCorrection::Run()
{
  int i, j, level;
  char buffer[256];
  double step, epsilon = 0, delta, maxChange = 0;

  // Print debugging information
  this->Debug("irtkCardiacSpatialCorrection::Run");

  if (_source == NULL) {
    cerr << "Registration::Run: Filter has no source input" << endl;
    exit(1);
  }

  if (_target == NULL) {
    cerr << "Registration::Run: Filter has no target input" << endl;
    exit(1);
  }

  if (_transformation == NULL) {
    cerr << "irtkCardiacSpatialCorrection::Run: Filter has no transformation output" << endl;
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

    // Save pre-processed images if we are debugging
    sprintf(buffer, "source_%d.nii.gz", level);
    if (_DebugFlag == true) _source->Write(buffer);
    sprintf(buffer, "target_%d.nii.gz", level);
    if (_DebugFlag == true) _target->Write(buffer);

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
        this->OptimizerRun(epsilon, maxChange, step);

        // Check whether we made any improvement or not
        if (epsilon > _Epsilon && maxChange > delta) { 
        } else {
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

void irtkCardiacSpatialCorrection::OptimizerRun(double &epsilon,double &maxChange,double &step){
	int i, j, k, n;
	double diff;

	n = 0;
	for(k = 0; k < _target->GetZ(); k++){
		n += _transformation[k]->NumberOfDOFs();
	}

	double *transformationParams = new double[n];

	j=0;
	for(k = 0; k < _target->GetZ(); k++){
		n = _transformation[k]->NumberOfDOFs();

		// Store the current state of the transformation.
		for (i = 0; i < n; ++i) {
			transformationParams[j] = _transformation[k]->Get(i);
			j++;
		}
	}

	// Epsilon is the change in similarity over an iteration of the optimizer.
	epsilon = this->RunStep(step);

	// MaxChange is the maximum change over the transformation parameters.
	maxChange = 0; j=0;
	for(k = 0; k < _target->GetZ(); k++){
		n = _transformation[k]->NumberOfDOFs();

		// Store the current state of the transformation.
		for (i = 0; i < n; ++i) {
			diff = fabs(_transformation[k]->Get(i) - transformationParams[j]);

			if (maxChange < diff) {
				maxChange = diff;
			}
			j++;
		}
	}

	delete []transformationParams;
}

double irtkCardiacSpatialCorrection::RunStep(double &step)
{
	int i, j, k, n;
	double similarity, new_similarity, old_similarity;

	// Assume that the transformation is the optimal transformation
	old_similarity = new_similarity = similarity = this->Evaluate();

	// Number of variables we have to optimize
	n = 0;
	for(k = 0; k < _target->GetZ(); k++){
		n += _transformation[k]->NumberOfDOFs();
	}

	// Convert some stuff to NR
	float *dx = new float[n];

	// Step along gradient direction until no further improvement is necessary
	do {
		new_similarity = similarity;

		this->EvaluateGradient(step, dx);

		j = 0;
		for(k = 0; k < _target->GetZ(); k++){
			for (i = 0; i < _transformation[k]->NumberOfDOFs(); i++) {
				_transformation[k]->Put(i, _transformation[k]->Get(i) + step * dx[j]);
				j++;
			}
		}
		similarity = this->Evaluate();
		if (similarity > new_similarity + _Epsilon) cout << similarity << endl;
	} while (similarity > new_similarity + _Epsilon);

	j = 0;
	for(k = 0; k < _target->GetZ(); k++){
		// Last step was no improvement, so back track
		for (i = 0; i < _transformation[k]->NumberOfDOFs(); i++) {
			_transformation[k]->Put(i, _transformation[k]->Get(i) - step * dx[j]);
			j++;
		}
	}

	// Delete NR memory
	delete []dx;

	if (new_similarity > old_similarity) {
		return new_similarity - old_similarity;
	} else {
		return 0;
	}
}

double irtkCardiacSpatialCorrection::Evaluate()
{

#ifndef HAS_TBB
	int i, j, k, t, n;
	double similarity;
#endif

	// Pointer to reference data
	irtkGreyPixel *ptr2target;

	// Print debugging information
	this->Debug("irtkCardiacSpatialCorrection::Evaluate");

	// Invert transformation
	//((irtkRigidTransformation *)_transformation)->Invert();

	// Create iterator
	irtkHomogeneousTransformationIterator *iterator
		= new irtkHomogeneousTransformationIterator[_target->GetZ()];

	for (k = 0; k < _target->GetZ(); k++) {
		iterator[k].SetTransformation((irtkHomogeneousTransformation *)_transformation[k]);
	}

	// Initialize metric
	_metric->Reset();

	// Pointer to voxels in target image
	ptr2target = _target->GetPointerToVoxels();

#ifdef HAS_TBB
	irtkMultiThreadedImageRigidRegistrationEvaluate evaluate(this);
	parallel_reduce(blocked_range<int>(0, _target->GetZ(), 20), evaluate);
#else

	for (t = 0; t < _target->GetT(); t++) {

		// Loop over all voxels in the target (reference) volume
		for (k = 0; k < _target->GetZ(); k++) {
			// Initialize iterator
			iterator[k].Initialize(_target, _source);
			// move to correct z slice
			for (n = 0; n < k; n++){
				iterator[k].NextZ();
			}
			for (j = 0; j < _target->GetY(); j++) {
				for (i = 0; i < _target->GetX(); i++) {
					// Check whether reference point is valid
					if (*ptr2target >= 0) {
						// Check whether transformed point is inside source volume
						if ((iterator[k]._x > _source_x1) && (iterator[k]._x < _source_x2) &&
							(iterator[k]._y > _source_y1) && (iterator[k]._y < _source_y2) &&
							(iterator[k]._z > _source_z1) && (iterator[k]._z < _source_z2)) {
								// Add sample to metric
								_metric->Add(*ptr2target, round(_interpolator->EvaluateInside(iterator[k]._x, iterator[k]._y, iterator[k]._z, t)));
						}
						iterator[k].NextX();
					} else {
						// Advance iterator by offset
						iterator[k].NextX(*ptr2target * -1);
						i          -= (*ptr2target) + 1;
						ptr2target -= (*ptr2target) + 1;
					}
					ptr2target++;
				}
				iterator[k].NextY();
			}
		}
	}

#endif


	// Invert transformation
	//((irtkRigidTransformation *)_transformation)->Invert();
	// Evaluate similarity measure
	delete []iterator;
	similarity = _metric->Evaluate();
	if(_IntLambda > 0)
		similarity += _IntLambda*this->InternalSimilarity();
	if(_ConstrainLambda > 0)
		similarity += _ConstrainLambda*this->ConstrainTransformation();

	return similarity;
}
	
double irtkCardiacSpatialCorrection::EvaluateGradient(float step, float *dx)
{
  int i, j, k;
  double s1, s2, norm, parameterValue;

  j = 0;
  for (k = 0; k < _target->GetZ(); k++){
	  for (i = 0; i < _transformation[k]->NumberOfDOFs(); i++) {
		  if (_transformation[k]->irtkTransformation::GetStatus(i) == _Active) {
			  parameterValue = _transformation[k]->Get(i);
			  _transformation[k]->Put(i, parameterValue + step);
			  s1 = this->Evaluate();
			  _transformation[k]->Put(i, parameterValue - step);
			  s2 = this->Evaluate();
			  _transformation[k]->Put(i, parameterValue);
			  dx[j] = s1 - s2;
		  } else {
			  dx[j] = 0;
		  }
		  j++;
	  }
  }

  // Calculate norm of vector
  j = 0; norm = 0;
  for (k = 0; k < _target->GetZ(); k++){
	  for (i = 0; i < _transformation[k]->NumberOfDOFs(); i++) {
		  norm += dx[j] * dx[j];
		  j++;
	  }
  }
  // Normalize vector
  norm = sqrt(norm);
  j = 0;
  for (k = 0; k < _target->GetZ(); k++){
	  if (norm > 0) {
		  for (i = 0; i < _transformation[k]->NumberOfDOFs(); i++) {
			  dx[j] /= norm;
			  j++;
		  }
	  } else {
		  for (i = 0; i < _transformation[k]->NumberOfDOFs(); i++) {
			  dx[j] = 0;
			  j++;
		  }
	  }
  }
  return norm;
}

double irtkCardiacSpatialCorrection::InternalSimilarity(){
	int i, j, k, t;
	double x, y, z;

	// Print debugging information
	this->Debug("irtkCardiacSpatialCorrection::InternalSimilarity()");
	_intmetric->Reset();

	
	for (k = 0; k < _target->GetZ() - 1; k++){
		((irtkAffineTransformation*)_transformation[k+1])->Invert();
		for (j = 0; j < _target->GetY(); j++){
			for (i = 0; i < _target->GetX(); i++){
				x = i; y = j; z = k;
				_target->ImageToWorld(x,y,z);
				_transformation[k]->Transform(x,y,z);
				_transformation[k+1]->Transform(x,y,z);
				_target->WorldToImage(x,y,z);
				z = z + 1;
				if(round(x)>=0 && round(x)<_target->GetX()
					&& round(y)>=0 && round(y)<_target->GetY()
					&& round(z)>=0 && round(z)<_target->GetZ()){
						for (t = 0; t < _target->GetT(); t++){
							_intmetric->Add(_target->GetAsDouble(i,j,k,t),
							_target->GetAsDouble(round(x),round(y),round(z),t));
						}
				}
			}
		}
		((irtkAffineTransformation*)_transformation[k+1])->Invert();
	}

	return _intmetric->Evaluate();
}

double irtkCardiacSpatialCorrection::ConstrainTransformation(){
	int i,k;
	double result;

	// Print debugging information
	this->Debug("irtkCardiacSpatialCorrection::ConstrainTransformation()");

	result = 0;

	for (k = 0; k < _target->GetZ(); k++){
		for (i = 0; i < 6; i++) {
				result += fabs(_transformation[k]->Get(i) - _origintransformation->Get(i));
		}
	}

	return -result;
}

bool irtkCardiacSpatialCorrection::Read(char *buffer1, char *buffer2, int &level)
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
  if (strstr(buffer1, "IntLambda") != NULL) {
	  this->_IntLambda = atof(buffer2);
	  ok = true;
  }
  if (strstr(buffer1, "TCLambda") != NULL) {
	  this->_ConstrainLambda = atof(buffer2);
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

  if (ok == false) {
    cerr << "irtkCardiacSpatialCorrection::Read: Can't parse line " << buffer1 << endl;
    exit(1);
  }

  return ok;
}

void irtkCardiacSpatialCorrection::Write(ostream &to)
{
	int i;

  to << "\n#\n# Registration parameters\n#\n\n";
  to << "No. of resolution levels          = " << this->_NumberOfLevels << endl;
  to << "No. of bins                       = " << this->_NumberOfBins << endl;
  to << "Epsilon                           = " << this->_Epsilon << endl;
  to << "IntLambda                         = " << this->_IntLambda << endl;
  to << "TCLambda                          = " << this->_ConstrainLambda << endl;
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
  default:
  	  cerr<<"Can not recognize similarity metric type!!!"<<endl;
  	  exit(1);
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
  	cerr << "irtkCardiacSpatialCorrection::Write: Interpolation mode not supported" << endl;
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
    to << "No. of steps                      = " << this->_NumberOfSteps[i] << endl;
    to << "Length of steps                   = " << this->_LengthOfSteps[i] << endl;
    to << "Delta                             = " << this->_Delta[i] << endl;
  }
}

void irtkCardiacSpatialCorrection::Read(char *filename)
{
  int level;
  char buffer1[255], *buffer2;

  ifstream from(filename);

  if (!from) {
    cerr << "irtkCardiacSpatialCorrection::Read: Can't open file " << filename
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

void irtkCardiacSpatialCorrection::Write(char *filename)
{
  ofstream to(filename);

  if (!to) {
    cerr << "irtkCardiacSpatialCorrection::Write: Can't open file " << filename;
    exit(1);
  }

  this->Write(to);
}

void irtkCardiacSpatialCorrection::GuessParameter()
{
	int i;
	double xsize, ysize, zsize;

	if ((_target == NULL) || (_source == NULL)) {
		cerr << "irtkImageRigidRegistration::GuessParameter: Target and source image not found" << endl;
		exit(1);
	}

	// Default parameters for registration
	_NumberOfLevels     = 1;
	_NumberOfBins       = 64;

	// Default parameters for optimization
	_SimilarityMeasure  = NMI;
	_Epsilon            = 0.0001;
	_IntLambda			= 0.3;
	_ConstrainLambda    = 0;

	// Read target pixel size
	_target->GetPixelSize(&xsize, &ysize, &zsize);

	// Default target parameters
	_TargetBlurring[0]      = 0;
	_TargetResolution[0][0] = xsize;
	_TargetResolution[0][1] = ysize;
	_TargetResolution[0][2] = zsize;

	// Read source pixel size
	_source->GetPixelSize(&xsize, &ysize, &zsize);

	// Default source parameters
	_SourceBlurring[0]      = 0;
	_SourceResolution[0][0] = xsize;
	_SourceResolution[0][1] = ysize;
	_SourceResolution[0][2] = zsize;

	// Remaining parameters
	
	_NumberOfIterations[0] = 20;
	_NumberOfSteps[0]      = 4;
	_LengthOfSteps[0]      = 4;

	// Try to guess padding by looking at voxel values in all eight corners of the volume:
	// If all values are the same we assume that they correspond to the padding value
	_TargetPadding = MIN_GREY;
	if ((_target->Get(_target->GetX()-1, 0, 0)                                 == _target->Get(0, 0, 0)) &&
		(_target->Get(0, _target->GetY()-1, 0)                                 == _target->Get(0, 0, 0)) &&
		(_target->Get(0, 0, _target->GetZ()-1)                                 == _target->Get(0, 0, 0)) &&
		(_target->Get(_target->GetX()-1, _target->GetY()-1, 0)                 == _target->Get(0, 0, 0)) &&
		(_target->Get(0, _target->GetY()-1, _target->GetZ()-1)                 == _target->Get(0, 0, 0)) &&
		(_target->Get(_target->GetX()-1, 0, _target->GetZ()-1)                 == _target->Get(0, 0, 0)) &&
		(_target->Get(_target->GetX()-1, _target->GetY()-1, _target->GetZ()-1) == _target->Get(0, 0, 0))) {
			_TargetPadding = _target->Get(0, 0, 0);
	}
}
