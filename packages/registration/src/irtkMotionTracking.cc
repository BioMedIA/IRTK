/*=========================================================================

  Library   : Image Registration Toolkit (IRTK)
  Module    : $Id$
  Copyright : Imperial College, Department of Computing
              Visual Information Processing (VIP), 2008 onwards
  Date      : $Date$
  Version   : $Revision$
  Changes   : $Author$

=========================================================================*/

#include <irtkRegistration.h>

#include <irtkMotionTracking.h>

#include <irtkResamplingWithPadding.h>

#include <irtkGaussianBlurring.h>

// Used as temporary memory for transformed intensities
irtkGreyImage *_tmpImage;

irtkGreyImage *tmp_image;

irtkMotionTracking::irtkMotionTracking()
{
  int i;

  for (i = 0; i < MAX_NO_RESOLUTIONS; i++) {
    // Default parameters for image
    _Blurring[i]      = 0;
    _Resolution[i][0] = 0;
    _Resolution[i][1] = 0;
    _Resolution[i][2] = 0;

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

  // Default optimization
  _OptimizationMethod = GradientDescent;

  // Default speedup factor
  _SpeedupFactor = 1;

  // Default parameters for motion tracking
  _DX          = 20;
  _DY          = 20;
  _DZ          = 20;

  // Default parameters for debugging
  _DebugFlag = false;

  // Set parameters
  _TargetPadding   = MIN_GREY;

  // Set inputs
  _image = NULL;
  _mask  = NULL;

  // Set output
  _transformation = NULL;

  // Set metric
  _metric = NULL;

  // Allocate interpolation object
  _interpolator = NULL;

  // Allocate optimizer object
  _optimizer = NULL;
}

irtkMotionTracking::~irtkMotionTracking()
{}

void irtkMotionTracking::Initialize()
{}

void irtkMotionTracking::Initialize(int level)
{
  int i, j, k, t;
  double dx, dy, dz, temp;
  irtkGreyPixel min, max, nbins;

  // Copy image to temp space
  tmp_image = new irtkGreyImage(*_image);

  // Swap source and target with temp space copies
  swap(tmp_image, _image);

  // Blur images if necessary
  if (_Blurring[level] > 0) {
    cout << "Blurring image ... "; cout.flush();
    irtkGaussianBlurring<irtkGreyPixel> blurring(_Blurring[level]);
    blurring.SetInput (_image);
    blurring.SetOutput(_image);
    blurring.Run();
    cout << "done" << endl;
  }

  _image->GetPixelSize(&dx, &dy, &dz);
  temp = fabs(_Resolution[0][0]-dx) + fabs(_Resolution[0][1]-dy) + fabs(_Resolution[0][2]-dz);

  if (level > 0 || temp > 0.000001) {
    cout << "Resampling image ... "; cout.flush();
    // Create resampling filter
    irtkResampling<irtkGreyPixel> resample(_Resolution[level][0], _Resolution[level][1], _Resolution[level][2]);
    resample.SetInput (_image);
    resample.SetOutput(_image);
    resample.Run();
    cout << "done" << endl;
  }

  // Find out the min and max values in target image, ignoring padding
  max = MIN_GREY;
  min = MAX_GREY;
  for (t = 0; t < _image->GetT(); t++) {
    for (k = 0; k < _image->GetZ(); k++) {
      for (j = 0; j < _image->GetY(); j++) {
        for (i = 0; i < _image->GetX(); i++) {
          if (_mask->Get(i, j, k) > 0) {
            if (_image->Get(i, j, k, t) > max)
              max = _image->Get(i, j, k, t);
            if (_image->Get(i, j, k, t) < min)
              min = _image->Get(i, j, k, t);
          }
        }
      }
    }
  }

  // Check whether dynamic range of data is not to large
  if (max - min > MAX_GREY) {
    cerr << this->NameOfClass()
         << "::Initialize: Dynamic range of target is too large" << endl;
    exit(1);
  } else {
    for (t = 0; t < _image->GetT(); t++) {
      for (k = 0; k < _image->GetZ(); k++) {
        for (j = 0; j < _image->GetY(); j++) {
          for (i = 0; i < _image->GetX(); i++) {
            if (_mask->Get(i, j, k) > 0) {
              _image->Put(i, j, k, t, _image->Get(i, j, k, t) - min);
            } else {
              _image->Put(i, j, k, t, -1);
            }
          }
        }
      }
    }
  }

  // Allocate memory for metric
  switch (_SimilarityMeasure) {
  case SSD:
    _metric = new irtkSSDSimilarityMetric;
    break;
  case CC:
    // Rescale images by an integer factor if necessary
    _metric = new irtkCrossCorrelationSimilarityMetric;
    break;
  case JE:
    // Rescale images by an integer factor if necessary
    nbins = irtkCalculateNumberOfBins(_image, _NumberOfBins, min, max);
    _metric = new irtkJointEntropySimilarityMetric(nbins, nbins);
    break;
  case MI:
    // Rescale images by an integer factor if necessary
    nbins = irtkCalculateNumberOfBins(_image, _NumberOfBins, min, max);
    _metric = new irtkMutualInformationSimilarityMetric(nbins, nbins);
    break;
  case NMI:
    // Rescale images by an integer factor if necessary
    nbins = irtkCalculateNumberOfBins(_image, _NumberOfBins, min, max);
    _metric = new irtkNormalisedMutualInformationSimilarityMetric(nbins, nbins);
    break;
  case CR_XY:
    // Rescale images by an integer factor if necessary
    nbins = irtkCalculateNumberOfBins(_image, _NumberOfBins, min, max);
    _metric = new irtkCorrelationRatioXYSimilarityMetric(nbins, nbins);
    break;
  case CR_YX:
    // Rescale images by an integer factor if necessary
    nbins = irtkCalculateNumberOfBins(_image, _NumberOfBins, min, max);
    _metric = new irtkCorrelationRatioYXSimilarityMetric(nbins, nbins);
    break;
  case LC:
    _metric = new irtkLabelConsistencySimilarityMetric;
    break;
  case K:
    // Rescale images by an integer factor if necessary
    nbins = irtkCalculateNumberOfBins(_image, _NumberOfBins, min, max);
    _metric = new irtkKappaSimilarityMetric(nbins, nbins);
    break;
  default:
    cerr << "irtkMotionTracking: Unknown metric" << endl;
    exit(1);
  }

  // Setup the interpolator
  _interpolator = irtkInterpolateImageFunction::New(_InterpolationMode, _image);

  // Setup interpolation for the source image
  _interpolator->SetInput(_image);
  _interpolator->Initialize();

  // Calculate the source image domain in which we can interpolate
  _interpolator->Inside(_x1,_y1, _z1, _x2, _y2, _z2);

  // Setup the optimizer
  switch (_OptimizationMethod) {
  case DownhillDescent:
    _optimizer = new irtkDownhillDescentOptimizer;
    break;
  case GradientDescent:
    _optimizer = new irtkGradientDescentOptimizer;
    break;
  case SteepestGradientDescent:
    _optimizer = new irtkSteepestGradientDescentOptimizer;
    break;
  case ConjugateGradientDescent:
    _optimizer = new irtkConjugateGradientDescentOptimizer;
    break;
  default:
    cerr << "irtkMotionTracking: Unkown optimizer" << endl;
    exit(1);
  }
  _optimizer->SetTransformation(_transformation);
  _optimizer->SetRegistration(this);

  // Print some debugging information
  cout << "Image sequence" << endl;
  _image->Print();
  cout << "Range is from " << min << " to " << max << endl;

  // Print initial transformation
  cout << "Initial transformation for level = " << level+1 << endl;;
  _transformation->Print();

}

void irtkMotionTracking::Finalize()
{}

void irtkMotionTracking::Finalize(int level)
{
  // Print final transformation
  cout << "Final transformation for level = " << level+1 << endl;;
  _transformation->Print();

  // Swap image back with temp space copies (see Initialize)
  swap(tmp_image, _image);

  delete tmp_image;
  delete _metric;
  delete _optimizer;
  delete _interpolator;
}

void irtkMotionTracking::Run()
{
  int i, j, level;
  char buffer[256];
  double step, epsilon, delta, maxChange;

  // Print debugging information
  this->Debug("irtkMotionTracking::Run");

  if (_image == NULL) {
    cerr << "Registration::Run: Filter has no image input" << endl;
    exit(1);
  }

  if (_mask == NULL) {
    cerr << "Registration::Run: Filter has no mask input" << endl;
    exit(1);
  }

  if (_transformation == NULL) {
    cerr << "irtkMotionTracking::Run: Filter has no transformation output" << endl;
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

    // Initialize for this level
    this->Initialize(level);

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
          this->Print();
        } else {
          sprintf(buffer, "log_%.3d_%.3d_%.3d.dof", level, i+1, j+1);
          this->Print();
          break;
        }
      }
      step = step / 2;
      delta = delta / 2.0;
    }

    // Do the final cleaning up for this level
    this->Finalize(level);

#ifdef HISTORY
    history->Print();
#endif

  }

  // Do the final cleaning up for all levels
  this->Finalize();
}

double irtkMotionTracking::EvaluateGradient(float step, float *dx)
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

bool irtkMotionTracking::Read(char *buffer1, char *buffer2, int &level)
{
  int i, n, ok = false;
  double dx, dy, dz;

  // Resolution level
  if (strstr(buffer1, "Resolution level") != NULL) {
    level = atoi(buffer2)-1;
    ok = true;
  }

  // Image blurring
  if (strstr(buffer1, "Image blurring (in mm)") != NULL) {
    if (level == -1) {
      for (i = 0; i < MAX_NO_RESOLUTIONS; i++) {
        this->_Blurring[i] = pow(2.0, double(i)) * atof(buffer2);
      }
    } else {
      this->_Blurring[level] = atof(buffer2);
    }
    ok = true;
  }
  // Image resolution
  if (strstr(buffer1, "Image resolution (in mm)") != NULL) {
    _image->GetPixelSize(&dx, &dy, &dz);
    if (level == -1) {
      n = sscanf(buffer2, "%lf %lf %lf", &(this->_Resolution[0][0]),  &(this->_Resolution[0][1]),  &(this->_Resolution[0][2]));
      if (n == 1) {
        this->_Resolution[0][1] = this->_Resolution[0][0];
        this->_Resolution[0][2] = this->_Resolution[0][0];
      }
      if (this->_Resolution[0][0] == 0) this->_Resolution[0][0] = dx;
      if (this->_Resolution[0][1] == 0) this->_Resolution[0][1] = dy;
      if (this->_Resolution[0][2] == 0) this->_Resolution[0][2] = dz;
      for (i = 0; i < MAX_NO_RESOLUTIONS; i++) {
        this->_Resolution[i][0] = pow(2.0, double(i)) * this->_Resolution[0][0];
        this->_Resolution[i][1] = pow(2.0, double(i)) * this->_Resolution[0][1];
        this->_Resolution[i][2] = pow(2.0, double(i)) * this->_Resolution[0][2];
      }
    } else {
      n = sscanf(buffer2, "%lf %lf %lf", &(this->_Resolution[level][0]),  &(this->_Resolution[level][1]),  &(this->_Resolution[level][2]));
      if (n == 1) {
        this->_Resolution[level][1] = this->_Resolution[level][0];
        this->_Resolution[level][2] = this->_Resolution[level][0];
      }
      if (this->_Resolution[level][0] == 0) this->_Resolution[level][0] = dx;
      if (this->_Resolution[level][1] == 0) this->_Resolution[level][1] = dy;
      if (this->_Resolution[level][2] == 0) this->_Resolution[level][2] = dz;
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
    cerr << "irtkMotionTracking::Read: Can't parse line " << buffer1 << endl;
    exit(1);
  }

  return ok;
}

void irtkMotionTracking::Write(ostream &to)
{
  int i;

  to << "\n#\n# Registration parameters\n#\n\n";
  to << "No. of resolution levels          = " << this->_NumberOfLevels << endl;
  to << "No. of bins                       = " << this->_NumberOfBins << endl;
  to << "Epsilon                           = " << this->_Epsilon << endl;
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
  	cerr << "irtkMotionTracking::Write: Interpolation mode not supported" << endl;
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
    to << "Image blurring (in mm)            = " << this->_Blurring[i] << endl;
    to << "Image resolution (in mm)          = " << this->_Resolution[i][0] << " " << this->_Resolution[i][1] << " " << this->_Resolution[i][2] << endl;
    to << "No. of iterations                 = " << this->_NumberOfIterations[i] << endl;
    to << "No. of steps                      = " << this->_NumberOfSteps[i] << endl;
    to << "Length of steps                   = " << this->_LengthOfSteps[i] << endl;
    to << "Delta                             = " << this->_Delta[i] << endl;
  }
}

void irtkMotionTracking::Read(char *filename)
{
  int level;
  char buffer1[255], *buffer2;

  ifstream from(filename);

  if (!from) {
    cerr << "irtkMotionTracking::Read: Can't open file " << filename
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

void irtkMotionTracking::Write(char *filename)
{
  ofstream to(filename);

  if (!to) {
    cerr << "irtkMotionTracking::Write: Can't open file " << filename;
    exit(1);
  }

  this->Write(to);
}

void irtkMotionTracking::UpdateLUT()
{
  int i, j, k, l;
  double x, y, z, t;

  // Print debugging information
  this->Debug("irtkImageFreeFormRegistration::UpdateLUT");

  float *ptr2affd = _affdLookupTable;
  float *ptr2mffd = _mffdLookupTable;
  for (l = 0; l < _image->GetT(); l++) {
    t = _image->ImageToTime(l);
    for (k = 0; k < _image->GetZ(); k++) {
      for (j = 0; j < _image->GetY(); j++) {
        for (i = 0; i < _image->GetX(); i++) {
          x = i;
          y = j;
          z = k;
          _image->ImageToWorld(x, y, z);
          *ptr2mffd = x;
          ptr2mffd++;
          *ptr2mffd = y;
          ptr2mffd++;
          *ptr2mffd = z;
          ptr2mffd++;
          _affd->LocalDisplacement(x, y, z, t);
          *ptr2affd = x;
          ptr2affd++;
          *ptr2affd = y;
          ptr2affd++;
          *ptr2affd = z;
          ptr2mffd++;
        }
      }
    }
  }
}

double irtkMotionTracking::Evaluate()
{
  // Image coordinates
  int i, j, k, l;
  // World coordinates
  double x, y, z, t;
  // Pointer to reference data
  irtkGreyPixel *ptr2image;
  irtkGreyPixel *ptr2mask;
  irtkGreyPixel *ptr2tmp;

  // Print debugging information
  this->Debug("irtkMotionTracking::Evaluate");

  // Initialize metric
  _metric->Reset();

  // Loop over all voxels in the target (reference) volume
  ptr2tmp   = _tmpImage->GetPointerToVoxels();
  for (l = 0; l < _image->GetT(); l++) {
    ptr2image = _image->GetPointerToVoxels();
    ptr2mask  = _mask ->GetPointerToVoxels();
    t = _image->ImageToTime(l);
    for (k = 0; k < _image->GetZ(); k++) {
      for (j = 0; j < _image->GetY(); j++) {
        for (i = 0; i < _image->GetX(); i++) {
          // Check whether reference point is valid
          if (*ptr2mask >= 0) {
            x = i;
            y = j;
            z = k;
            _image->ImageToWorld(x, y, z);
            _affd->LocalDisplacement(x, y, z, t);
            _image->WorldToImage(x, y, z);
            // Check whether transformed point is inside volume
            if ((x > _x1) && (x < _x2) && (y > _y1) && (y < _y2) && (z > _z1) && (z < _z2)) {
              // Add sample to metric
              *ptr2tmp =  round(_interpolator->EvaluateInside(x, y, z, l));
              _metric->Add(*ptr2image, *ptr2tmp);
            } else {
              *ptr2tmp = -1;
            }
          }
          // Increment pointers to next voxel
          ptr2tmp++;
          ptr2mask++;
          ptr2image++;
        }
      }
    }
  }

  // Evaluate similarity measure
  double similarity = _metric->Evaluate();

  // Return similarity measure + penalty terms
  return similarity;
}

double irtkMotionTracking::EvaluateDerivative(int index, double step)
{
  float *ptr;
  irtkPoint p1, p2;
  double bi, bj, bk, dx, dy, dz, p[3];
  int i, j, k, t, i1, i2, j1, j2, k1, k2, t1, t2, dim;
  irtkGreyPixel *ptr2image, *ptr2mask, *ptr2tmp;
  irtkSimilarityMetric *tmpMetricA, *tmpMetricB;

  // Print debugging information
  this->Debug("irtkMotionTracking::EvaluateDerivative(int, double)");

  tmpMetricA = _tmpMetricA;
  tmpMetricB = _tmpMetricB;

  // Initialize metrics for forward and backward derivative steps
  tmpMetricA->Reset(_metric);
  tmpMetricB->Reset(_metric);

  // Calculate bounding box of control point in world coordinates
  _affd->BoundingBox(index, p1, p2);
  _image->WorldToImage(p1);
  _image->WorldToImage(p2);

  // Calculate bounding box of control point in image coordinates
  _affd->BoundingBox(_image, index, i1, j1, k1, t1, i2, j2, k2, t2, 1.0 / _SpeedupFactor);

  // Calculate incremental changes in lattice coordinates when looping
  // over target
  dx = (FFDLOOKUPTABLESIZE-1)/(p2._x-p1._x);
  dy = (FFDLOOKUPTABLESIZE-1)/(p2._y-p1._y);
  dz = (FFDLOOKUPTABLESIZE-1)/(p2._z-p1._z);

  // Calculate whether this DOF corresponds to x, y or z-displacement
  dim = int(index / (_affd->GetX()*_affd->GetY()*_affd->GetZ()*_affd->GetT()));

  // Loop over all voxels in the target (reference) volume
  for (t = 0; t < _image->GetT(); t++) {
    for (k = k1; k <= k2; k++) {
      bk = step * _localLookupTable[round(dz*(k-p1._z))];
      for (j = j1; j <= j2; j++) {
        ptr2image = _image->GetPointerToVoxels(i1, j, k);
        ptr2mask  = _mask ->GetPointerToVoxels(i1, j, k);
        ptr        = &(_affdLookupTable[3*_image->VoxelToIndex(i1, j, k)]);
        bj = bk * _localLookupTable[round(dy*(j-p1._y))];
        ptr2tmp  = _tmpImage->GetPointerToVoxels(i1, j, k, t);
        for (i = i1; i <= i2; i++) {

          // Check whether reference point is valid
          if (*ptr2mask >= 0) {
            bi = bj * _localLookupTable[round(dx*(i-p1._x))];

            // Delete old samples from both metrics
            if (*ptr2tmp != -1) {
              tmpMetricA->Delete(*ptr2image, *ptr2tmp);
              tmpMetricB->Delete(*ptr2image, *ptr2tmp);
            }

            p[0] = ptr[0];
            p[1] = ptr[1];
            p[2] = ptr[2];
            p[dim] += bi;

            // Convert transformed point to image coordinates
            _image->WorldToImage(p[0], p[1], p[2]);

            // Check whether transformed point is inside volume
            if ((p[0] > _x1) && (p[0] < _x2) && (p[1] > _y1) && (p[1] < _y2) && (p[2] > _z1) && (p[2] < _z2)) {

              // Add sample to metric
              tmpMetricA->Add(*ptr2image, round(_interpolator->EvaluateInside(p[0], p[1], p[2], t)));
            }

            p[0] = ptr[0];
            p[1] = ptr[1];
            p[2] = ptr[2];
            p[dim] -= bi;

            // Convert transformed point to image coordinates
            _image->WorldToImage(p[0], p[1], p[2]);

            // Check whether transformed point is inside volume
            if ((p[0] > _x1) && (p[0] < _x2) && (p[1] > _y1) && (p[1] < _y2) && (p[2] > _z1) && (p[2] < _z2)) {

              // Add sample to metric
              tmpMetricB->Add(*ptr2image, round(_interpolator->EvaluateInside(p[0], p[1], p[2], t)));
            }
          }

          // Increment pointers to next voxel
          ptr2image++;
          ptr2mask++;
          ptr2tmp++;
        }
      }
    }
  }

  // Evaluate similarity measure
  double similarityA = tmpMetricA->Evaluate();

  // Evaluate similarity measure
  double similarityB = tmpMetricB->Evaluate();

  return similarityA - similarityB;
}

