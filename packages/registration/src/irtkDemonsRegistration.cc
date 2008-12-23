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

irtkDemonsRegistration::irtkDemonsRegistration()
{
  // Default parameters for target image
  _TargetBlurring     = 0;
  _TargetResolution   = 0;

  // Default parameters for source image
  _SourceBlurring     = 0;
  _SourceResolution   = 0;

  // Default parameters for registration
  _NumberOfLevels     = 1;
  _NumberOfIterations = 10;
  _ReductionFactor    = 2;

  // Default parameters for smoothing
  _Smoothing          = 4;

  // Default parameters for optimization
  _InterpolationMode  = Interpolation_Linear;

  _TargetPadding = MIN_GREY;
  _SourcePadding = MIN_GREY;

  // Set inputs
  _target = NULL;
  _source = NULL;

  // Set output
  _transformation = NULL;

  // Allocate interpolation object
  _interpolator = new irtkLinearInterpolateImageFunction<irtkGreyPixel>;
}

irtkDemonsRegistration::~irtkDemonsRegistration()
{
  delete _interpolator;
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

  if (_transformation == NULL) {
    cerr << "irtkDemonsRegistration::Run: Filter has no transformation output" << endl;
    exit(1);
  }

  // Push local transformation back on transformation stack
  _ffd = (irtkLinearFreeFormTransformation *)_transformation->PopLocalTransformation();

  // Invert transformation
  _transformation->Invert();

  // Blur images if necessary
  if (_TargetBlurring > 0) {
    cout << "Blurring target ... "; cout.flush();
    irtkGaussianBlurringWithPadding<irtkGreyPixel> blurring(_TargetBlurring,
        _TargetPadding);
    blurring.SetInput (_target);
    blurring.SetOutput(_target);
    blurring.Run();
    cout << "done" << endl;
  }

  if (_SourceBlurring > 0) {
    cout << "Blurring source ... "; cout.flush();
    irtkGaussianBlurringWithPadding<irtkGreyPixel> blurring(_SourceBlurring,
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
    irtkResamplingWithPadding<irtkGreyPixel> resample(_TargetResolution,
        _TargetResolution,
        _TargetResolution,
        _TargetPadding);
    resample.SetInput (_target);
    resample.SetOutput(_target);
    resample.Run();
    cout << "done" << endl;
  }
  if (_SourceResolution > 0) {
    cout << "Resampling source ... "; cout.flush();
    // Create resampling filter
    irtkResamplingWithPadding<irtkGreyPixel> resample(_SourceResolution,
        _SourceResolution,
        _SourceResolution,
        _SourcePadding);
    resample.SetInput (_source);
    resample.SetOutput(_source);
    resample.Run();
    cout << "done" << endl;
  }

  cout << "Target image (reference)" << endl;
  _target->Print();

  cout << "Source image (transform)" << endl;
  _source->Print();
}

void irtkDemonsRegistration::Initialize(int level)
{
  int i, j, k;
  double x, y, z;

  // Push local transformation back on transformation stack
  _transformation->PushLocalTransformation(_ffd);

  // Create global displacement field with same size as target
  _globalDX.Initialize(*_target);
  _globalDY.Initialize(*_target);
  _globalDZ.Initialize(*_target);

  // Calculate global displacement field
  for (k = 0; k < _target->GetZ(); k++) {
    for (j = 0; j < _target->GetY(); j++) {
      for (i = 0; i < _target->GetX(); i++) {
        // Current point
        x = i;
        y = j;
        z = k;
        // Transform point into world coordinates
        _target->ImageToWorld(x, y, z);
        // Calculate displacement
        _transformation->Displacement(x, y, z);
        // Store displacement
        _globalDX(i, j, k) = x;
        _globalDY(i, j, k) = y;
        _globalDZ(i, j, k) = z;
      }
    }
  }

  // Push local transformation back on transformation stack
  _ffd = (irtkLinearFreeFormTransformation *)_transformation->PopLocalTransformation();

  // Create global displacement field with same size as target
  _localDX.Initialize(*_target);
  _localDY.Initialize(*_target);
  _localDZ.Initialize(*_target);

  // Setup interpolation for the source image
  _interpolator->SetInput(_source);
  _interpolator->Initialize();

  // Calculate the source image domain in which we can interpolate
  _interpolator->Inside(_source_x1, _source_y1, _source_z1,
                        _source_x2, _source_y2, _source_z2);

  // Temporary image equals source image
  _tmp = *_target;

  // Update tmp image
  for (k = 0; k < _target->GetZ(); k++) {
    for (j = 0; j < _target->GetY(); j++) {
      for (i = 0; i < _target->GetX(); i++) {
        x = i;
        y = j;
        z = k;
        // Transform point into world coordinates
        _target->ImageToWorld(x, y, z);

        // Transform point
        x += _globalDX(i, j, k);
        y += _globalDY(i, j, k);
        z += _globalDZ(i, j, k);

        // Transform point into image coordinates
        _source->WorldToImage(x, y, z);

        // Interpolate intensity
        if ((x > _source_x1) && (x < _source_x2) &&
            (y > _source_y1) && (y < _source_y2) &&
            (z > _source_z1) && (z < _source_z2)) {
          _tmp(i, j, k) = round(_interpolator->EvaluateInside(x, y, z));
        }
      }
    }
  }

  // Create images for gradient
  _targetGradientX.Initialize(*_target);
  _targetGradientY.Initialize(*_target);
  _targetGradientZ.Initialize(*_target);

  // Calculate gradient
  for (k = 1; k < _target->GetZ()-1; k++) {
    for (j = 1; j < _target->GetY()-1; j++) {
      for (i = 1; i < _target->GetX()-1; i++) {
        _targetGradientX(i, j, k) = _target->Get(i+1, j, k) - _target->Get(i-1, j, k);
        _targetGradientY(i, j, k) = _target->Get(i, j+1, k) - _target->Get(i, j-1, k);
        _targetGradientZ(i, j, k) = _target->Get(i, j, k+1) - _target->Get(i, j, k-1);
      }
    }
  }
}

void irtkDemonsRegistration::Finalize()
{
  // Push local transformation back on transformation stack
  _transformation->PushLocalTransformation(_ffd);

  // Invert transformation
  _transformation->Invert();
}

void irtkDemonsRegistration::Finalize(int level)
{
  int i, j, k;
  double x, y, z, x1, y1, z1, x2, y2, z2;

  // Create interpolators for displacement fields
  irtkLinearInterpolateImageFunction<irtkRealPixel> *interpolatorDX =
    new irtkLinearInterpolateImageFunction<irtkRealPixel>;
  interpolatorDX->SetInput(&_globalDX);
  interpolatorDX->Initialize();
  irtkLinearInterpolateImageFunction<irtkRealPixel> *interpolatorDY =
    new irtkLinearInterpolateImageFunction<irtkRealPixel>;
  interpolatorDY->SetInput(&_globalDY);
  interpolatorDY->Initialize();
  irtkLinearInterpolateImageFunction<irtkRealPixel> *interpolatorDZ =
    new irtkLinearInterpolateImageFunction<irtkRealPixel>;
  interpolatorDZ->SetInput(&_globalDZ);
  interpolatorDZ->Initialize();

  for (k = 0; k < _ffd->GetZ(); k++) {
    for (j = 0; j < _ffd->GetY(); j++) {
      for (i = 0; i < _ffd->GetX(); i++) {
        x = i;
        y = j;
        z = k;
        // Current point into world coordinates
        _ffd->LatticeToWorld(x, y, z);

        // Calculate previous displacement
        x1 = x;
        y1 = y;
        z1 = z;
        _transformation->Displacement(x1, y1, z1);

        // Calculate current displacement
        _target->WorldToImage(x, y, z);
        x2 =  interpolatorDX->Evaluate(x, y, z);
        y2 =  interpolatorDY->Evaluate(x, y, z);
        z2 =  interpolatorDZ->Evaluate(x, y, z);

        _ffd->Put(i, j, k, x2 - x1, y2 - y1, z2 - z1);
      }
    }
  }

  // Delete interpolators for displacement fields
  delete interpolatorDX;
  delete interpolatorDY;
  delete interpolatorDZ;
}

void irtkDemonsRegistration::Run(irtkGreyImage target, irtkGreyImage source, int level)
{
  int i, j, k, iter;
  double x, y, z, xsize, ysize, zsize, ssd, grad;
  irtkGreyImage *ptr2target, *ptr2source;

  // Resample target image
  target.GetPixelSize(&xsize, &ysize, &zsize);
  irtkResamplingWithPadding<irtkGreyPixel> resample1(_ReductionFactor * xsize,
      _ReductionFactor * ysize,
      _ReductionFactor * zsize,
      -1);
  resample1.SetInput (&target);
  resample1.SetOutput(&target);
  resample1.Run();

  // Resample source image
  source.GetPixelSize(&xsize, &ysize, &zsize);
  irtkResamplingWithPadding<irtkGreyPixel> resample2(_ReductionFactor * xsize,
      _ReductionFactor * ysize,
      _ReductionFactor * zsize,
      -1);
  resample2.SetInput (&source);
  resample2.SetOutput(&source);
  resample2.Run();

  // Here appears the magic
  ptr2target = _target;
  ptr2source = _source;
  _target = &target;
  _source = &source;

  // Recursion
  if (level < _NumberOfLevels) {
    this->Run(target, source, level + 1);
  }

  // Initialize for this level
  this->Initialize(level);

  // Print resolution level
  cout << "Resolution level no. " << level << endl;

  // Run the registration filter
  for (iter = 0; iter < _NumberOfIterations; iter++) {

    // Print out iteration number
    cout << "Iteration = " << iter + 1 << " (out of " << _NumberOfIterations << ") "<< endl;

    for (k = 0; k < _target->GetZ(); k++) {
      for (j = 0; j < _target->GetY(); j++) {
        for (i = 0; i < _target->GetX(); i++) {
          ssd  = _target->Get(i, j, k) - _tmp(i, j, k);
          grad = _targetGradientX(i, j, k) * _targetGradientX(i, j, k) +
                 _targetGradientY(i, j, k) * _targetGradientY(i, j, k) +
                 _targetGradientZ(i, j, k) * _targetGradientZ(i, j, k);
          if ((fabs(ssd) > 0) && (grad > 0)) {
            _localDX(i, j, k) = ssd / (grad + ssd * ssd) * _targetGradientX(i, j, k);
            _localDY(i, j, k) = ssd / (grad + ssd * ssd) * _targetGradientY(i, j, k);
            _localDZ(i, j, k) = ssd / (grad + ssd * ssd) * _targetGradientZ(i, j, k);
          } else {
            _localDX(i, j, k) = 0;
            _localDY(i, j, k) = 0;
            _localDZ(i, j, k) = 0;
          }
        }
      }
    }

    irtkGaussianBlurring<irtkRealPixel> blurring(_Smoothing);

    // Smooth X-displacement
    blurring.SetInput(&_localDX);
    blurring.SetOutput(&_localDX);
    blurring.Run();

    // Smooth Y-displacement
    blurring.SetInput(&_localDY);
    blurring.SetOutput(&_localDY);
    blurring.Run();

    // Smooth Z-displacement
    blurring.SetInput(&_localDZ);
    blurring.SetOutput(&_localDZ);
    blurring.Run();

    // Add displacement field
    _globalDX += _localDX;
    _globalDY += _localDY;
    _globalDZ += _localDZ;

    // Update tmp image
    for (k = 0; k < _target->GetZ(); k++) {
      for (j = 0; j < _target->GetY(); j++) {
        for (i = 0; i < _target->GetX(); i++) {
          x = i;
          y = j;
          z = k;
          // Transform point into world coordinates
          _target->ImageToWorld(x, y, z);

          // Transform point
          x += _globalDX(i, j, k);
          y += _globalDY(i, j, k);
          z += _globalDZ(i, j, k);

          // Transform point into image coordinates
          _source->WorldToImage(x, y, z);

          // Interpolate intensity
          if ((x > _source_x1) && (x < _source_x2) &&
              (y > _source_y1) && (y < _source_y2) &&
              (z > _source_z1) && (z < _source_z2)) {
            _tmp(i, j, k) = round(_interpolator->EvaluateInside(x, y, z));
          }
        }
      }
    }
  }

  // Do the final cleaning up for this level
  this->Finalize(level);

  // Here disappers the magic
  _target = ptr2target;
  _source = ptr2source;
}

void irtkDemonsRegistration::Run()
{
  int i, j, k, iter;
  double x, y, z, ssd, grad;

  // Do the initial set up for all levels
  this->Initialize();

  // Run the registration filter at high resolutions
  if (_NumberOfLevels > 1) {
    this->Run(*_target, *_source, 2);
  }

  // Initialize for this level
  this->Initialize(1);

  // Start on this level
  cout << "Resolution level no. 1" << endl;

  // Run the registration filter
  for (iter = 0; iter < _NumberOfIterations; iter++) {

    // Print out iteration number
    cout << "Iteration = " << iter + 1 << " (out of " << _NumberOfIterations << ") "<< endl;

    for (k = 0; k < _target->GetZ(); k++) {
      for (j = 0; j < _target->GetY(); j++) {
        for (i = 0; i < _target->GetX(); i++) {
          ssd  = _target->Get(i, j, k) - _tmp(i, j, k);
          grad = _targetGradientX(i, j, k) * _targetGradientX(i, j, k) +
                 _targetGradientY(i, j, k) * _targetGradientY(i, j, k) +
                 _targetGradientZ(i, j, k) * _targetGradientZ(i, j, k);
          if ((fabs(ssd) > 0) && (grad > 0)) {
            _localDX(i, j, k) = ssd / (grad + ssd * ssd) * _targetGradientX(i, j, k);
            _localDY(i, j, k) = ssd / (grad + ssd * ssd) * _targetGradientY(i, j, k);
            _localDZ(i, j, k) = ssd / (grad + ssd * ssd) * _targetGradientZ(i, j, k);
          } else {
            _localDX(i, j, k) = 0;
            _localDY(i, j, k) = 0;
            _localDZ(i, j, k) = 0;
          }
        }
      }
    }

    irtkGaussianBlurring<irtkRealPixel> blurring(_Smoothing);

    // Smooth X-displacement
    blurring.SetInput(&_localDX);
    blurring.SetOutput(&_localDX);
    blurring.Run();

    // Smooth Y-displacement
    blurring.SetInput(&_localDY);
    blurring.SetOutput(&_localDY);
    blurring.Run();

    // Smooth Z-displacement
    blurring.SetInput(&_localDZ);
    blurring.SetOutput(&_localDZ);
    blurring.Run();

    // Add displacement field
    _globalDX += _localDX;
    _globalDY += _localDY;
    _globalDZ += _localDZ;

    // Update tmp image
    for (k = 0; k < _target->GetZ(); k++) {
      for (j = 0; j < _target->GetY(); j++) {
        for (i = 0; i < _target->GetX(); i++) {
          x = i;
          y = j;
          z = k;
          // Transform point into world coordinates
          _target->ImageToWorld(x, y, z);

          // Transform point
          x += _globalDX(i, j, k);
          y += _globalDY(i, j, k);
          z += _globalDZ(i, j, k);

          // Transform point into image coordinates
          _source->WorldToImage(x, y, z);

          // Interpolate intensity
          if ((x > _source_x1) && (x < _source_x2) &&
              (y > _source_y1) && (y < _source_y2) &&
              (z > _source_z1) && (z < _source_z2)) {
            _tmp(i, j, k) = round(_interpolator->EvaluateInside(x, y, z));
          }
        }
      }
    }
  }

  // Do the final cleaning up for this level
  this->Finalize(1);

  // Do the final cleaning up for all levels
  this->Finalize();
}

void irtkDemonsRegistration::SetParameter(const irtkDemonsRegistration *r)
{
  // Default parameters for target image
  _TargetBlurring     = r->_TargetBlurring;
  _TargetResolution   = r->_TargetResolution;

  // Default parameters for source image
  _SourceBlurring     = r->_SourceBlurring;
  _SourceResolution   = r->_SourceResolution;

  // Default parameters for registration
  _NumberOfLevels     = r->_NumberOfLevels;
  _NumberOfIterations = r->_NumberOfIterations;
  _ReductionFactor    = r->_ReductionFactor;

  // Default parameters for smoothing
  _Smoothing          = r->_Smoothing;

  // Default parameters for optimization
  _InterpolationMode  = r->_InterpolationMode;
}

void irtkDemonsRegistration::Read(char *filename)
{
  ifstream from(filename);
  if (!from) {
    cerr << "irtkDemonsRegistration::Read: Can't open file " << filename
         << endl;
    exit(1);
  }

  // Read out base registration parameters
  from >> this;
}

void irtkDemonsRegistration::Write(char *filename)
{
  ofstream to(filename);
  if (!to) {
    cerr << "irtkRegistration::Write: Can't open file " << filename;
    exit(1);
  }

  // Write out base registration parameters
  to << this;
}

istream& operator>> (istream &from, irtkDemonsRegistration *r)
{
  int ok;
  char buffer1[255], *buffer2;

  while (from.eof() != True) {

    ok = False;
    if (read_line(from, buffer1, buffer2) == 0) break;

    // Target blurring
    if (strstr(buffer1, "Target blurring (in mm)") != NULL) {
      r->_TargetBlurring = atof(buffer2);
      cout << "Target blurring is ... " << r->_TargetBlurring << endl;
      ok = True;
    }
    // Target resolution
    if (strstr(buffer1, "Target resolution (in mm)") != NULL) {
      r->_TargetResolution = atof(buffer2);
      cout << "Target resolution is ... " << r->_TargetResolution << endl;
      ok = True;
    }
    // Source blurring
    if (strstr(buffer1, "Source blurring (in mm)") != NULL) {
      r->_SourceBlurring = atof(buffer2);
      cout << "Source blurring is ... " << r->_SourceBlurring << endl;
      ok = True;
    }
    // Source resolution
    if (strstr(buffer1, "Source resolution (in mm)") != NULL) {
      r->_SourceResolution = atof(buffer2);
      cout << "Source resolution is ... " << r->_SourceResolution << endl;
      ok = True;
    }
    if (strstr(buffer1, "No. of resolution levels") != NULL) {
      r->_NumberOfLevels = atoi(buffer2);
      cout << "No. of resolution levels is ... " << r->_NumberOfLevels << endl;
      ok = True;
    }
    if (strstr(buffer1, "No. of iterations") != NULL) {
      r->_NumberOfIterations = atoi(buffer2);
      cout << "No. of iterations is ... " << r->_NumberOfIterations << endl;
      ok = True;
    }
    if (strstr(buffer1, "Reduction factor") != NULL) {
      r->_ReductionFactor = atof(buffer2);
      cout << "Reduction factor for images is ... " << r->_ReductionFactor
           << endl;
      ok = True;
    }
    if (strstr(buffer1, "Smoothing factor") != NULL) {
      r->_Smoothing = atof(buffer2);
      cout << "Smoothing factor for deformations is ... " << r->_Smoothing
           << endl;
      ok = True;
    }
    // Check if we parse every line
    if (ok != True) {
      cerr << "istream& operator>> irtkRegistration: Ignoring line "
           << buffer2 << endl;
    }
  }

  return from;
}

ostream& operator<< (ostream &to, const irtkDemonsRegistration *r)
{
  to << "#\n# Target image parameters\n#\n\n";
  to << "Target blurring (in mm)           = " << r->_TargetBlurring << endl;
  to << "Target resolution (in mm)         = " << r->_TargetResolution << endl;

  to << "\n#\n# Source image parameters\n#\n\n";
  to << "Source blurring (in mm)           = " << r->_SourceBlurring << endl;
  to << "Source resolution (in mm)         = " << r->_SourceResolution << endl;

  to << "\n#\n# Registration parameters\n#\n\n";
  to << "No. of resolution levels          = " << r->_NumberOfLevels << endl;
  to << "No. of iterations                 = " << r->_NumberOfIterations << endl;
  to << "Reduction factor                  = " << r->_ReductionFactor << endl;
  to << "Smoothing factor                  = " << r->_Smoothing << endl;

  return to;
}

