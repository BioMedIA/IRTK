/*=========================================================================

  Library   : Image Registration Toolkit (IRTK)
  Module    : $Id$
  Copyright : Imperial College, Department of Computing
              Visual Information Processing (VIP), 2009 onwards
  Date      : $Date$
  Version   : $Revision$
  Changes   : $Author$

=========================================================================*/

#include <irtkRegistration.h>

#include <irtkResampling.h>

#include <irtkGradientImage.h>

#include <irtkGaussianBlurring.h>

#define HISTORY

irtkGreyImage *tmp_image;

#ifdef HAS_VTK

irtkModelRegistration::irtkModelRegistration()
{
  int i;

  for (i = 0; i < MAX_NO_RESOLUTIONS; i++) {
    // Default parameters for image
    _ImageBlurring[i]      = 0;
    _ImageResolution[i][0] = 0;
    _ImageResolution[i][1] = 0;
    _ImageResolution[i][2] = 0;

    // Default parameters for optimization
    _NumberOfIterations[i] = 20;
    _NumberOfSteps[i]      = 5;
    _LengthOfSteps[i]      = 2;
    _Delta[i]              = 0;
  }

  // Default parameters for registration
  _NumberOfLevels     = 1;

  // Default parameters for optimization
  _OptimizationMethod = DownhillDescent;
  _Epsilon            = 0;

  // Default parameters for debugging
  _DebugFlag = false;

  // Set inputs
  _image = NULL;
  _model = NULL;

  // Set output
  _transformation = NULL;

  // Allocate optimizer object
  _optimizer = NULL;

#ifdef HISTORY
  history = new irtkHistory;
#endif
}

irtkModelRegistration::~irtkModelRegistration()
{
#ifdef HISTORY
  delete history;
#endif
}

void irtkModelRegistration::Initialize()
{
}

void irtkModelRegistration::Initialize(int level)
{
  double dx, dy, dz, temp;

  // Copy image to temp space
  tmp_image = new irtkGreyImage(*_image);

  // Swap source and target with temp space copies
  swap(tmp_image, _image);

  // Blur image if necessary
  if (_ImageBlurring[level] > 0) {
    cout << "Blurring image ... "; cout.flush();
    irtkGaussianBlurring<irtkGreyPixel> blurring(_ImageBlurring[level]);
    blurring.SetInput (_image);
    blurring.SetOutput(_image);
    blurring.Run();
    cout << "done" << endl;
  }

  _image->GetPixelSize(&dx, &dy, &dz);
  temp = fabs(_ImageResolution[0][0]-dx) + fabs(_ImageResolution[0][1]-dy) + fabs(_ImageResolution[0][2]-dz);

  if (level > 0 || temp > 0.000001) {
    cout << "Resampling image ... "; cout.flush();
    // Create resampling filter
    irtkResampling<irtkGreyPixel> resample(_ImageResolution[level][0], _ImageResolution[level][1], _ImageResolution[level][2]);
    resample.SetInput (_image);
    resample.SetOutput(_image);
    
    // Create interpolator for resampling
    irtkInterpolateImageFunction *interpolator = new irtkLinearInterpolateImageFunction;
    resample.SetInterpolator(interpolator);
    resample.Run();
    delete interpolator;
    
    cout << "done" << endl;
  }

  // Setup metric
  _model->GetPointData()->SetActiveScalars("IntensityProfile");
	if (_model->GetPointData()->GetScalars() == NULL){
		cout << "Using gradient similarity metric" << endl;
		_metric = new irtkModelGradientSimilarityMetric(_image);
	} else {
		cout << "Using correlation similarity metric" << endl;
		_metric = new irtkModelCorrelationSimilarityMetric(_image);
	}

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
    cerr << "Unknown optimizer" << endl;
    exit(1);
  }
  _optimizer->SetTransformation(_transformation);
  _optimizer->SetRegistration(this);

  // Print some debugging information
  _image->Print();

  // Print initial transformation
  cout << "Initial transformation for level = " << level+1 << endl;;
  _transformation->Print();

}

void irtkModelRegistration::Finalize()
{}

void irtkModelRegistration::Finalize(int level)
{
  // Print final transformation
  cout << "Final transformation for level = " << level+1 << endl;;
  _transformation->Print();

  // Swap source and target back with temp space copies (see Initialize)
  swap(tmp_image, _image);

  delete tmp_image;
  delete _optimizer;
}

void irtkModelRegistration::Run()
{
  int i, j, level;
  char buffer[256];
  double step, epsilon, delta, maxChange;

  // Print debugging information
  this->Debug("irtkModelRegistration::Run");

  if (_model == NULL) {
    cerr << "Registration::Run: Filter has no model input" << endl;
    exit(1);
  }

  if (_image == NULL) {
    cerr << "Registration::Run: Filter has no image input" << endl;
    exit(1);
  }

  if (_transformation == NULL) {
    cerr << "irtkModelRegistration::Run: Filter has no transformation output" << endl;
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

    // Do the final cleaning up for this level
    this->Finalize(level);

#ifdef HISTORY
    history->Print();
#endif

  }

  // Do the final cleaning up for all levels
  this->Finalize();
}

double irtkModelRegistration::EvaluateGradient(float step, float *dx)
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

bool irtkModelRegistration::Read(char *buffer1, char *buffer2, int &level)
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
        this->_ImageBlurring[i] = pow(2.0, double(i)) * atof(buffer2);
      }
    } else {
      this->_ImageBlurring[level] = atof(buffer2);
    }
    ok = true;
  }
  // Image resolution
  if (strstr(buffer1, "Image resolution (in mm)") != NULL) {
    _image->GetPixelSize(&dx, &dy, &dz);
    if (level == -1) {
      n = sscanf(buffer2, "%lf %lf %lf", &(this->_ImageResolution[0][0]),  &(this->_ImageResolution[0][1]),  &(this->_ImageResolution[0][2]));
      if (n == 1) {
        this->_ImageResolution[0][1] = this->_ImageResolution[0][0];
        this->_ImageResolution[0][2] = this->_ImageResolution[0][0];
      }
      if (this->_ImageResolution[0][0] == 0) this->_ImageResolution[0][0] = dx;
      if (this->_ImageResolution[0][1] == 0) this->_ImageResolution[0][1] = dy;
      if (this->_ImageResolution[0][2] == 0) this->_ImageResolution[0][2] = dz;
      for (i = 0; i < MAX_NO_RESOLUTIONS; i++) {
        this->_ImageResolution[i][0] = pow(2.0, double(i)) * this->_ImageResolution[0][0];
        this->_ImageResolution[i][1] = pow(2.0, double(i)) * this->_ImageResolution[0][1];
        this->_ImageResolution[i][2] = pow(2.0, double(i)) * this->_ImageResolution[0][2];
      }
    } else {
      n = sscanf(buffer2, "%lf %lf %lf", &(this->_ImageResolution[level][0]),  &(this->_ImageResolution[level][1]),  &(this->_ImageResolution[level][2]));
      if (n == 1) {
        this->_ImageResolution[level][1] = this->_ImageResolution[level][0];
        this->_ImageResolution[level][2] = this->_ImageResolution[level][0];
      }
      if (this->_ImageResolution[level][0] == 0) this->_ImageResolution[level][0] = dx;
      if (this->_ImageResolution[level][1] == 0) this->_ImageResolution[level][1] = dy;
      if (this->_ImageResolution[level][2] == 0) this->_ImageResolution[level][2] = dz;
    }
    ok = true;
  }
  if (strstr(buffer1, "No. of resolution levels") != NULL) {
    this->_NumberOfLevels = atoi(buffer2);
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
  
  /*
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
*/
  
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
    cerr << "irtkModelRegistration::Read: Can't parse line " << buffer1 << endl;
    exit(1);
  }

  return ok;
}

void irtkModelRegistration::Write(ostream &to)
{
  int i;

  to << "\n#\n# Registration parameters\n#\n\n";
  to << "No. of resolution levels          = " << this->_NumberOfLevels << endl;  
  to << "Epsilon                           = " << this->_Epsilon << endl;

  /*
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
  */

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
    to << "Image blurring (in mm)            = " << this->_ImageBlurring[i] << endl;
    to << "Image resolution (in mm)          = " << this->_ImageResolution[i][0] << " " << this->_ImageResolution[i][1] << " " << this->_ImageResolution[i][2] << endl;
    to << "No. of iterations                 = " << this->_NumberOfIterations[i] << endl;
    to << "No. of steps                      = " << this->_NumberOfSteps[i] << endl;
    to << "Length of steps                   = " << this->_LengthOfSteps[i] << endl;
    to << "Delta                             = " << this->_Delta[i] << endl;
  }
}

void irtkModelRegistration::Read(char *filename)
{
  int level;
  char buffer1[255], *buffer2;

  ifstream from(filename);

  if (!from) {
    cerr << "irtkModelRegistration::Read: Can't open file " << filename
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

void irtkModelRegistration::Write(char *filename)
{
  ofstream to(filename);

  if (!to) {
    cerr << "irtkModelRegistration::Write: Can't open file " << filename;
    exit(1);
  }

  this->Write(to);
}

#endif
