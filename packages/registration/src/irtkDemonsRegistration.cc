/*=========================================================================

  Library   : Image Registration Toolkit (IRTK)
  Module    : $Id$
  Copyright : Imperial College, Department of Computing
              Visual Information Processing (VIP), 2008 onwards
  Date      : $Date$
  Version   : $Revision$
  Changes   : $Author$

=========================================================================*/

#include <irtkDemonsRegistration.h>

#include <irtkResamplingWithPadding.h>

#include <irtkGaussianBlurring.h>

#include <irtkGradientImageFilter.h>

#define EPSILON 0.001

int regrid = False;

irtkRealImage *tmp_target_dem, *tmp_source_dem;

irtkDemonsRegistration::irtkDemonsRegistration()
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

    // Default parameters for smoothing
    _Smoothing[i]           = 0;

  }

  // Default parameters for registration
  _NumberOfLevels     = 1;
  _NumberOfIterations = 200;
  _StepSize           = 0.005;
  _Epsilon            = 0.001;
  _Regridding         = 10;
  _Symmetric          = False;

  _TargetPadding = MIN_GREY;
  _SourcePadding = MIN_GREY;

  // Set inputs
  _target = NULL;
  _source = NULL;

  // Set output
  _transformation1 = NULL;
  _transformation2 = NULL;

  // Set ffds
  _ffd1 = NULL;
  _ffd2 = NULL;

  // Debug flag
  _DebugFlag = False;

  // Allocate interpolation object
  _interpolator1 = new irtkLinearInterpolateImageFunction;
  _interpolator2 = new irtkLinearInterpolateImageFunction;

  // Allocate test transformation
  _transhist1 = new irtkFluidFreeFormTransformation;
  _transhist2 = new irtkFluidFreeFormTransformation;
}

irtkDemonsRegistration::~irtkDemonsRegistration()
{
  delete _interpolator1;
  delete _interpolator2;
}

void irtkDemonsRegistration::GuessParameter()
{
  if ((_target == NULL) || (_source == NULL)) {
    cerr << "irtkImageRigidRegistration::GuessParameter: Target and source image not found" << endl;
    exit(1);
  }
}

void irtkDemonsRegistration::Initialize()
{
  if (_source == NULL) {
    cerr << "Registration::Run: Filter has no source input" << endl;
    exit(1);
  }

  if (_target == NULL) {
    cerr << "Registration::Run: Filter has no target input" << endl;
    exit(1);
  }

  if ((_transformation1 == NULL) || (_transformation2 == NULL)) {
    cerr << "irtkDemonsRegistration::Run: Filter has no transformation output" << endl;
    exit(1);
  }
}

void irtkDemonsRegistration::Initialize(int level)
{
  irtkImageAttributes attr;

  // Copy source and target to temp space
  tmp_target_dem = new irtkRealImage(*_target);
  tmp_source_dem = new irtkRealImage(*_source);

  // Swap source and target with temp space copies
  swap(tmp_target_dem, _target);
  swap(tmp_source_dem, _source);

  // Blur images if necessary
  if (_TargetBlurring[level] > 0) {
    cout << "Blurring target ... "; cout.flush();
    irtkGaussianBlurringWithPadding<irtkRealPixel> blurring(_TargetBlurring[level],
        _TargetPadding);
    blurring.SetInput (_target);
    blurring.SetOutput(_target);
    blurring.Run();
    cout << "done" << endl;
  }

  if (_SourceBlurring[level] > 0) {
    cout << "Blurring source ... "; cout.flush();
    irtkGaussianBlurringWithPadding<irtkRealPixel> blurring(_SourceBlurring[level],
        _SourcePadding);
    blurring.SetInput (_source);
    blurring.SetOutput(_source);
    blurring.Run();
    cout << "done" << endl;
  }

  // Resample images if necessary
  if (_TargetResolution > 0) {
    cout << "Resampling target ... "; cout.flush();
    // Create resampling filter
    irtkResamplingWithPadding<irtkRealPixel> resample(_TargetResolution[level][0],
        _TargetResolution[level][1],
        _TargetResolution[level][2],
        _TargetPadding);
    resample.SetInput (_target);
    resample.SetOutput(_target);
    resample.Run();
    cout << "done" << endl;
  }
  if (_SourceResolution > 0) {
    cout << "Resampling source ... "; cout.flush();
    // Create resampling filter
    irtkResamplingWithPadding<irtkRealPixel> resample(_SourceResolution[level][0],
        _SourceResolution[level][1],
        _SourceResolution[level][2],
        _SourcePadding);
    resample.SetInput (_source);
    resample.SetOutput(_source);
    resample.Run();
    cout << "done" << endl;
  }

  if (_ffd1 == NULL) {
    // We are at the highest level of the pyramid, so create initial free-form deformations
    _ffd1 = new irtkLinearFreeFormTransformation(*_target, _target->GetXSize(), _target->GetYSize(), _target->GetZSize());
    _ffd2 = new irtkLinearFreeFormTransformation(*_source, _source->GetXSize(), _source->GetYSize(), _source->GetZSize());

    // Push local transformation back on transformation stack
    _transformation1->PushLocalTransformation(_ffd1);
    _transformation2->PushLocalTransformation(_ffd2);
  } else {
    // Create linear free-form deformations
    _ffd1 = new irtkLinearFreeFormTransformation(*_target, _target->GetXSize(), _target->GetYSize(), _target->GetZSize());
    _ffd2 = new irtkLinearFreeFormTransformation(*_source, _source->GetXSize(), _source->GetYSize(), _source->GetZSize());

    // Compose exisiting free-formations with new ones
    _ffd1->Compose(_transformation1->PopLocalTransformation());
    _ffd2->Compose(_transformation2->PopLocalTransformation());

    // Push local transformation back on transformation stack
    _transformation1->PushLocalTransformation(_ffd1);
    _transformation2->PushLocalTransformation(_ffd2);
  }

  // Extract image attributes
  attr = _target->GetImageAttributes();
  attr._t = 3;

  // Create local displacement field with same size as target
  _local1.Initialize(attr);
  _local2.Initialize(attr);

  // Setup interpolation for the source image
  _interpolator1->SetInput(_source);
  _interpolator1->Initialize();
  _interpolator2->SetInput(_target);
  _interpolator2->Initialize();

  // Calculate the source image domain in which we can interpolate
  _interpolator1->Inside(_source_x1, _source_y1, _source_z1,
                         _source_x2, _source_y2, _source_z2);
  _interpolator2->Inside(_target_x1, _target_y1, _target_z1,
                         _target_x2, _target_y2, _target_z2);

  // Temporary image equals target image
  _sourceTmp = *_target;
  _targetTmp = *_target;

  // Generate transformed tmp images
  _imagetransformation1.SetInput (_source, _transformation1);
  _imagetransformation1.SetOutput(&_sourceTmp);
  _imagetransformation1.PutInterpolator(_interpolator1);
  _imagetransformation1.Run();

  _imagetransformation2.SetInput (_target, _transformation2);
  _imagetransformation2.SetOutput(&_targetTmp);
  _imagetransformation2.PutInterpolator(_interpolator2);
  _imagetransformation2.Run();

  // Extract image attributes
  attr = _target->GetImageAttributes();

  // Compute gradient of target image
  irtkGradientImageFilter<irtkRealPixel> gradient1(irtkGradientImageFilter<irtkRealPixel>::GRADIENT_VECTOR);
  gradient1.SetInput (_target);
  gradient1.SetOutput(&_targetGradient);
  gradient1.Run();

  // Compute gradient of source image
  irtkGradientImageFilter<irtkRealPixel> gradient2(irtkGradientImageFilter<irtkRealPixel>::GRADIENT_VECTOR);
  gradient2.SetInput (_source);
  gradient2.SetOutput(&_sourceGradient);
  gradient2.Run();

  cout << "Target image (reference)" << endl;
  _target->Print();

  cout << "Source image (transform)" << endl;
  _source->Print();

  // Print level
  cout << "Starting level " << level+1 << endl;;
}

void irtkDemonsRegistration::Finalize()
{}

void irtkDemonsRegistration::Finalize(int level)
{
  // Print level
  cout << "Finished level " << level+1 << endl;;

  // Swap source and target back with temp space copies (see Initialize)
  swap(tmp_target_dem, _target);
  swap(tmp_source_dem, _source);

  delete tmp_target_dem;
  delete tmp_source_dem;
}

void irtkDemonsRegistration::Smooth(double sigma)
{
  if (sigma > 0) {
    irtkGaussianBlurring<double> blurring(sigma);

    // Smooth displacement
    blurring.SetInput(&_local1);
    blurring.SetOutput(&_local1);
    blurring.Run();

    // Smooth displacement
    if (_Symmetric) {
      blurring.SetInput(&_local2);
      blurring.SetOutput(&_local2);
      blurring.Run();
    }
  }

  if (_StepSize > 0) {
    _local1 *= _StepSize;
    _local2 *= _StepSize;
  }
}

void irtkDemonsRegistration::Update()
{
  irtkLinearFreeFormTransformation *ffd1 = new irtkLinearFreeFormTransformation(_local1);
  irtkLinearFreeFormTransformation *ffd2 = new irtkLinearFreeFormTransformation(_local2);

  // Concatenate displacement fields
  _ffd1->Compose(ffd1);
  _ffd2->Compose(ffd2);

  delete ffd1;
  delete ffd2;

  // Update source image
  _imagetransformation1.Run();

  // Compute gradient of source image
  irtkGradientImageFilter<irtkRealPixel> gradient1(irtkGradientImageFilter<irtkRealPixel>::GRADIENT_VECTOR);
  gradient1.SetInput (&_sourceTmp);
  gradient1.SetOutput(&_sourceGradient);
  gradient1.Run();

  // If symmetric version:
  if (_Symmetric) {

    // Update target image
    _imagetransformation2.Run();

    // Compute gradient of target image
    irtkGradientImageFilter<irtkRealPixel> gradient2(irtkGradientImageFilter<irtkRealPixel>::GRADIENT_VECTOR);
    gradient2.SetInput (&_targetTmp);
    gradient2.SetOutput(&_targetGradient);
    gradient2.Run();
  }
}

void irtkDemonsRegistration::Run()
{
  int iter, level;

  // Do the initial set up
  this->Initialize();

  // Save pre-processed images if we are debugging
  if (_DebugFlag == True) {
    char buffer[255];
    sprintf(buffer, "source_%d.nii.gz", 1);
    _source->Write(buffer);
    sprintf(buffer, "target_%d.nii.gz", 1);
    _target->Write(buffer);
  }

  // Loop over levels
  for (level = _NumberOfLevels-1; level >= 0; level--) {

    // Print resolution level
    cout << "Resolution level no. " << level+1 << endl;

    // Initialize for this level
    this->Initialize(level);

    // Run the registration filter
    for (iter = 0; iter < _NumberOfIterations; iter++) {

      // Print out iteration number
      cout << "Iteration = " << iter + 1 << " (out of " << _NumberOfIterations << ") "<< endl;

      // Compute force
      if (Force() < _Epsilon) {
        break;
      }

      // Smooth force
      Smooth(_Smoothing[level]);

      // Automatically compute step size
      /*
      if (iter == 0) {
        double min, max;
        _local1.GetMinMax(&min, &max);
        _StepSize = 0.5 / (fabs(min) > fabs(max) ? fabs(min) : fabs(max));
        cout << "Automatically computing step size as " << _StepSize << endl;
        _local1 *= _StepSize;
        _local2 *= _StepSize;
      }
*/

      if (_DebugFlag == True) {
        char buffer[255];

        sprintf(buffer, "local_%.3d.nii.gz", iter);
        _local1.Write(buffer);
        sprintf(buffer, "tmp_source_%.3d.nii.gz", iter);
        _sourceTmp.Write(buffer);
        sprintf(buffer, "tmp_target_%.3d.nii.gz", iter);
        _targetTmp.Write(buffer);
      }

      // Regrid
      if ((iter+1) % _Regridding == 0) {
        regrid = True;
        cout << "regridding at iteration = " << iter+1 << endl;
      } else {
        regrid = False;
      }
      regrid = False;

      // Update images and displacement fields
      Update();
    }

    // Do the final cleaning up
    this->Finalize(level);
  }

  // Do the final cleaning up
  this->Finalize();
}

double irtkDemonsRegistration::Force()
{
  int i, j, k;
  double ssd, mag, rms;
  irtkMatrix m(3, 3);

  m(0, 0) = _target->GetImageAttributes()._xaxis[0];
  m(1, 0) = _target->GetImageAttributes()._xaxis[1];
  m(2, 0) = _target->GetImageAttributes()._xaxis[2];
  m(0, 1) = _target->GetImageAttributes()._yaxis[0];
  m(1, 1) = _target->GetImageAttributes()._yaxis[1];
  m(2, 1) = _target->GetImageAttributes()._yaxis[2];
  m(0, 2) = _target->GetImageAttributes()._zaxis[0];
  m(1, 2) = _target->GetImageAttributes()._zaxis[1];
  m(2, 2) = _target->GetImageAttributes()._zaxis[2];

  rms = 0;
  for (k = 0; k < _target->GetZ(); k++) {
    for (j = 0; j < _target->GetY(); j++) {
      for (i = 0; i < _target->GetX(); i++) {
        mag = _sourceGradient(i, j, k, 0)*_sourceGradient(i, j, k, 0) + _sourceGradient(i, j, k, 1)*_sourceGradient(i, j, k, 1) +
              _sourceGradient(i, j, k, 2)*_sourceGradient(i, j, k, 2);
        ssd = (_targetTmp(i, j, k) - _sourceTmp(i, j, k)) / (mag + 0.0001 + (_targetTmp(i, j, k) - _sourceTmp(i, j, k)) * (_targetTmp(i, j, k) - _sourceTmp(i, j, k)));
        _local1(i, j, k, 0) = ssd * (m(0, 0) * _sourceGradient(i, j, k, 0) + m(0, 1) * _sourceGradient(i, j, k, 1) + m(0, 2) * _sourceGradient(i, j, k, 2));
        _local1(i, j, k, 1) = ssd * (m(1, 0) * _sourceGradient(i, j, k, 0) + m(1, 1) * _sourceGradient(i, j, k, 1) + m(1, 2) * _sourceGradient(i, j, k, 2));
        _local1(i, j, k, 2) = ssd * (m(2, 0) * _sourceGradient(i, j, k, 0) + m(2, 1) * _sourceGradient(i, j, k, 1) + m(2, 2) * _sourceGradient(i, j, k, 2));
        rms += fabs(_targetTmp(i, j, k) - _sourceTmp(i, j, k));
      }
    }
  }

  if (_Symmetric) {
    for (k = 0; k < _target->GetZ(); k++) {
      for (j = 0; j < _target->GetY(); j++) {
        for (i = 0; i < _target->GetX(); i++) {
          ssd  = _sourceTmp(i, j, k) - _targetTmp(i, j, k);
          _local2(i, j, k, 0) = ssd * (m(0, 0) * _targetGradient(i, j, k, 0) + m(0, 1) * _targetGradient(i, j, k, 1) + m(0, 2) * _targetGradient(i, j, k, 2));
          _local2(i, j, k, 1) = ssd * (m(1, 0) * _targetGradient(i, j, k, 0) + m(1, 1) * _targetGradient(i, j, k, 1) + m(1, 2) * _targetGradient(i, j, k, 2));
          _local2(i, j, k, 2) = ssd * (m(2, 0) * _targetGradient(i, j, k, 0) + m(2, 1) * _targetGradient(i, j, k, 1) + m(2, 2) * _targetGradient(i, j, k, 2));
          rms += _local2(i, j, k, 0)*_local2(i, j, k, 0) + _local2(i, j, k, 1)*_local2(i, j, k, 1) + _local2(i, j, k, 2)*_local2(i, j, k, 2);
        }
      }
    }
  }
  cout << "SSD Metric = " << sqrt(rms) / _target->GetNumberOfVoxels() << endl;
  return sqrt(rms) / _target->GetNumberOfVoxels();
}

double irtkDemonsRegistration::Force2()
{
  double rms;
  int i, j, k, l, m, n, nsample, size = 3;

  // Create images for gradient
  irtkGenericImage<double> Imean(_target->GetImageAttributes());
  irtkGenericImage<double> Jmean(_target->GetImageAttributes());
  irtkGenericImage<double> A(_target->GetImageAttributes());
  irtkGenericImage<double> B(_target->GetImageAttributes());
  irtkGenericImage<double> C(_target->GetImageAttributes());

  for (k = 0; k < _target->GetZ(); k++) {
    for (j = 0; j < _target->GetY(); j++) {
      for (i = 0; i < _target->GetX(); i++) {
        nsample = 0;
        Imean(i, j, k) = 0;
        Jmean(i, j, k) = 0;
        for (n = -size; n <= size; n++) {
          if ((k+n >= 0) && (k+n < _target->GetZ())) {
            for (m = -size; m <= size; m++) {
              if ((j+m >= 0) && (j+m < _target->GetY())) {
                for (l = -size; l <= size; l++) {
                  if ((i+l >= 0) && (i+l < _target->GetX())) {
                    Imean(i, j, k) += _targetTmp(i+l, j+m, k+n);
                    Jmean(i, j, k) += _sourceTmp(i+l, j+m, k+n);
                    nsample++;
                  }
                }
              }
            }
          }
        }
        Imean(i, j, k) /= (double)nsample;
        Jmean(i, j, k) /= (double)nsample;
      }
    }
  }

  for (k = 0; k < _target->GetZ(); k++) {
    for (j = 0; j < _target->GetY(); j++) {
      for (i = 0; i < _target->GetX(); i++) {
        nsample = 0;
        for (n = -size; n <= size; n++) {
          if ((k+n >= 0) && (k+n < _target->GetZ())) {
            for (m = -size; m <= size; m++) {
              if ((j+m >= 0) && (j+m < _target->GetY())) {
                for (l = -size; l <= size; l++) {
                  if ((i+l >= 0) && (i+l < _target->GetX())) {
                    A(i, j, k) += (_targetTmp(i+l, j+m, k+n) - Imean(i+l, j+m, k+n)) *
                                  (_sourceTmp(i+l, j+m, k+n) - Jmean(i+l, j+m, k+n));
                    B(i, j, k) += (_targetTmp(i+l, j+m, k+n) - Imean(i+l, j+m, k+n)) *
                                  (_targetTmp(i+l, j+m, k+n) - Imean(i+l, j+m, k+n));
                    C(i, j, k) += (_sourceTmp(i+l, j+m, k+n) - Jmean(i+l, j+m, k+n)) *
                                  (_sourceTmp(i+l, j+m, k+n) - Jmean(i+l, j+m, k+n));
                    nsample++;
                  }
                }
              }
            }
          }
        }
        A(i, j, k) /= (double)nsample;
        B(i, j, k) /= (double)nsample;
        C(i, j, k) /= (double)nsample;
      }
    }
  }

  rms = 0;
  for (k = 0; k < _target->GetZ(); k++) {
    for (j = 0; j < _target->GetY(); j++) {
      for (i = 0; i < _target->GetX(); i++) {
        if ((B(i, j, k) > EPSILON) && (C(i, j, k) > EPSILON)) {
          double x = -2 * A(i, j, k) / (B(i, j, k) * C(i, j, k)) * (Jmean(i, j, k) - A(i, j, k) / B(i, j, k) * Imean(i, j, k));
          if ((fabs(Imean(i, j, k) - _targetTmp(i, j, k)) > EPSILON) && (fabs(Jmean(i, j, k) - _sourceTmp(i, j, k)) > EPSILON)) {
            _local1(i, j, k, 0) = x * _sourceGradient(i, j, k, 0);
            _local1(i, j, k, 1) = x * _sourceGradient(i, j, k, 1);
            _local1(i, j, k, 2) = x * _sourceGradient(i, j, k, 2);
          } else {
            _local1(i, j, k, 0) = 0;
            _local1(i, j, k, 1) = 0;
            _local1(i, j, k, 2) = 0;
          }
        } else {
          _local1(i, j, k, 0) = 0;
          _local1(i, j, k, 1) = 0;
          _local1(i, j, k, 2) = 0;
        }
      }
    }
  }

  if (_Symmetric) {
    for (k = 0; k < _target->GetZ(); k++) {
      for (j = 0; j < _target->GetY(); j++) {
        for (i = 0; i < _target->GetX(); i++) {
          _local2(i, j, k, 0) = 0; //x * _targetGradient(i, j, k, 0) * _StepSize;
          _local2(i, j, k, 1) = 0; //x * _targetGradient(i, j, k, 1) * _StepSize;
          _local2(i, j, k, 2) = 0; //x * _targetGradient(i, j, k, 2) * _StepSize;
          //        }
        }
      }
    }
  }
  cout << "NCC Metric = " << sqrt(rms) / _target->GetNumberOfVoxels() << endl;
  return sqrt(rms) / _target->GetNumberOfVoxels();
}

Bool irtkDemonsRegistration::Read(char *buffer1, char *buffer2, int &level)
{
  int i, n, ok = False;
  double dx, dy, dz;

  // Resolution level
  if (strstr(buffer1, "Resolution level") != NULL) {
    level = atoi(buffer2)-1;
    ok = True;
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
    ok = True;
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
    ok = True;
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
    ok = True;
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
    ok = True;
  }
  if (strstr(buffer1, "No. of resolution levels") != NULL) {
    this->_NumberOfLevels = atoi(buffer2);
    ok = True;
  }
  if (strstr(buffer1, "No. of iterations") != NULL) {
    this->_NumberOfIterations = atoi(buffer2);
    ok = True;
  }
  if (strstr(buffer1, "Regridding") != NULL) {
    this->_Regridding = atoi(buffer2);
    ok = True;
  }
  // Source blurring
  if (strstr(buffer1, "Smoothing") != NULL) {
    if (level == -1) {
      for (i = 0; i < MAX_NO_RESOLUTIONS; i++) {
        this->_Smoothing[i] = pow(2.0, double(i)) * atof(buffer2);
      }
    } else {
      this->_Smoothing[level] = atof(buffer2);
    }
    ok = True;
  }
  if (strstr(buffer1, "Step size") != NULL) {
    this->_StepSize = atoi(buffer2);
    ok = True;
  }
  if (strstr(buffer1, "Epsilon") != NULL) {
    this->_Epsilon = atof(buffer2);
    ok = True;
  }
  if (strstr(buffer1, "Padding value") != NULL) {
    this->_TargetPadding = atoi(buffer2);
    ok = True;
  }

  if (ok == False) {
    cerr << "irtkDemonsRegistration::Read: Can't parse line " << buffer1 << endl;
    exit(1);
  }

  return ok;
}

void irtkDemonsRegistration::Write(ostream &to)
{
  int i;

  to << "\n#\n# Registration parameters\n#\n\n";
  to << "No. of resolution levels          = " << this->_NumberOfLevels << endl;
  to << "Epsilon                           = " << this->_Epsilon << endl;
  to << "Padding value                     = " << this->_TargetPadding << endl;
  to << "No. of iterations                 = " << this->_NumberOfIterations << endl;
  to << "Step size                         = " << this->_StepSize << endl;
  to << "Regridding                        = " << this->_Regridding << endl;

  for (i = 0; i < this->_NumberOfLevels; i++) {
    to << "\n#\n# Registration parameters for resolution level " << i+1 << "\n#\n\n";
    to << "Resolution level                  = " << i+1 << endl;
    to << "Target blurring (in mm)           = " << this->_TargetBlurring[i] << endl;
    to << "Target resolution (in mm)         = " << this->_TargetResolution[i][0] << " " << this->_TargetResolution[i][1] << " " << this->_TargetResolution[i][2] << endl;
    to << "Source blurring (in mm)           = " << this->_SourceBlurring[i] << endl;
    to << "Source resolution (in mm)         = " << this->_SourceResolution[i][0] << " " << this->_SourceResolution[i][1] << " " << this->_SourceResolution[i][2] << endl;
    to << "Smoothing                         = " << this->_Smoothing[i] << endl;
  }
}

void irtkDemonsRegistration::Read(char *filename)
{
  int level;
  char buffer1[255], *buffer2;

  ifstream from(filename);

  if (!from) {
    cerr << "irtkDemonsRegistration::Read: Can't open file " << filename
    << endl;
    exit(1);
  }

  level = -1;
  while (from.eof() != True) {
    if (read_line(from, buffer1, buffer2) != 0) {
      if (this->Read(buffer1, buffer2, level) == False) {
        cerr << "Couldn't parse line: " << buffer1 << endl;
      }
    }
  }
}

void irtkDemonsRegistration::Write(char *filename)
{
  ofstream to(filename);

  if (!to) {
    cerr << "irtkDemonsRegistration::Write: Can't open file " << filename;
    exit(1);
  }

  this->Write(to);
}
