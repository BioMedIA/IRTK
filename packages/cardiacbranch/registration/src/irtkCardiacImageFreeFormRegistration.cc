#include <irtkRegistration.h>

#include <irtkGaussianBlurring.h>

// Used as temporary memory for transformed intensities
irtkGreyImage *_tmpImage;

// Used as temporary memory for transformed intensities
irtkGreyImage *_ttmpImage;

// Used as temporary memory for weight on untagged image
irtkGreyImage *_weight;

// The original target and source images
extern irtkGreyImage *tmp_target, *tmp_source;

irtkGreyImage *tmp_ttarget, *tmp_tsource;

#include <irtkMultiThreadedImageFreeFormRegistration.h>

irtkCardiacImageFreeFormRegistration::irtkCardiacImageFreeFormRegistration()
{
  // Print debugging information
  this->Debug("irtkCardiacImageFreeFormRegistration::irtkCardiacImageFreeFormRegistration");

  // Allocate interpolation object
  _tinterpolator = NULL;

  // Default optimization
  _OptimizationMethod = GradientDescent;

  // Default speedup factor
  _SpeedupFactor = 1;

  // Default parameters for non-rigid registration
  _Lambda1     = 0;
  _Lambda2     = 0;
  _Lambda3     = 0;
  _Weight      = 1;
  _Tweight     = 1;
  _Thresholdmax= 100;
  _Thresholdmin= 10;
  _DX          = 20;
  _DY          = 20;
  _DZ          = 20;
  _Subdivision = True;
  _Mode        = RegisterXYZ;
}

void irtkCardiacImageFreeFormRegistration::GuessParameter()
{
  int i;
  double xsize, ysize, zsize, spacing;

  if ((_target == NULL) || (_source == NULL) ||(_ttarget == NULL) || (_tsource == NULL)  ) {
    cerr << "irtkCardiacImageFreeFormRegistration::GuessParameter: one of the tagged or untagged image not found" << endl;
    exit(1);
  }

  // Default parameters for registration
  _NumberOfLevels     = 3;
  _NumberOfBins       = 64;

  // Default parameters for optimization
  _SimilarityMeasure  = NMI;
  _OptimizationMethod = GradientDescent;
  _Epsilon            = 0.0001;

  // Read target pixel size
  _target->GetPixelSize(&xsize, &ysize, &zsize);

  // Use xsize as spacing
  spacing = xsize;

  // Default target parameters
  _TargetBlurring[0]      = GuessResolution(xsize, ysize, zsize) / 2.0;
  _TargetResolution[0][0] = GuessResolution(xsize, ysize, zsize);
  _TargetResolution[0][1] = GuessResolution(xsize, ysize, zsize);
  _TargetResolution[0][2] = GuessResolution(xsize, ysize, zsize);

  for (i = 1; i < _NumberOfLevels; i++) {
    _TargetBlurring[i]      = _TargetBlurring[i-1] * 2;
    _TargetResolution[i][0] = _TargetResolution[i-1][0] * 2;
    _TargetResolution[i][1] = _TargetResolution[i-1][1] * 2;
    _TargetResolution[i][2] = _TargetResolution[i-1][2] * 2;
  }

  // Read source pixel size
  _source->GetPixelSize(&xsize, &ysize, &zsize);

  // Default source parameters
  _SourceBlurring[0]      = GuessResolution(xsize, ysize, zsize) / 2.0;
  _SourceResolution[0][0] = GuessResolution(xsize, ysize, zsize);
  _SourceResolution[0][1] = GuessResolution(xsize, ysize, zsize);
  _SourceResolution[0][2] = GuessResolution(xsize, ysize, zsize);

  for (i = 1; i < _NumberOfLevels; i++) {
    _SourceBlurring[i]      = _SourceBlurring[i-1] * 2;
    _SourceResolution[i][0] = _SourceResolution[i-1][0] * 2;
    _SourceResolution[i][1] = _SourceResolution[i-1][1] * 2;
    _SourceResolution[i][2] = _SourceResolution[i-1][2] * 2;
  }

  // Default parameters for non-rigid registration
  _Lambda1            = 0;
  _Lambda2            = 0;
  _Lambda3            = 0;
  _Weight             = 1;
  _Tweight            = 1;
  _Thresholdmax       = 110;
  _Thresholdmin       = 10;
  _DX                 =_target->GetX() * spacing / 10.0;
  _DY                 =_target->GetX() * spacing / 10.0;
  _DZ                 =_target->GetX() * spacing / 10.0;
  _Subdivision        = True;

  // Remaining parameters
  for (i = 0; i < _NumberOfLevels; i++) {
    _NumberOfIterations[i] = 10;
    _NumberOfSteps[i]      = 4;
    _LengthOfSteps[i]      = _DX / 8.0 * pow(2.0, i);
  }

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

void irtkCardiacImageFreeFormRegistration::Initialize()
{
  // Print debugging information
  this->Debug("irtkCardiacImageFreeFormRegistration::Initialize");

  // Initialize base class
  this->irtkImageRegistration::Initialize();
  
  // Pointer to multi-level FFD
  _mffd = (irtkMultiLevelFreeFormTransformation *)_transformation;

  // Create FFD
  if (_mffd->NumberOfLevels() == 0) {
    _affd = new irtkBSplineFreeFormTransformation(*_target, this->_DX, this->_DY, this->_DZ);
  } else {
    _affd = (irtkBSplineFreeFormTransformation *)_mffd->PopLocalTransformation();
  }

  // Initialize pointers
  _tmpImage         = NULL;
  _ttmpImage        = NULL;
  _weight           = NULL;
  _affdLookupTable  = NULL;
  _mffdLookupTable  = NULL;
  _localLookupTable = new float [FFDLOOKUPTABLESIZE];

}

void irtkCardiacImageFreeFormRegistration::Initialize(int level)
{
  int i, j, k, n, t;
  double u, x, y, z;
  double dx, dy, dz, temp;
  float *ptr;
  irtkGreyPixel ttarget_min, ttarget_max, ttarget_nbins;
  irtkGreyPixel tsource_min, tsource_max, tsource_nbins;

  // Print debugging information
  this->Debug("irtkCardiacImageFreeFormRegistration::Initialize(int)");

  // Initialize base class
  this->irtkImageRegistration::Initialize(level);

   // Copy source and target to temp space
  tmp_ttarget = new irtkGreyImage(*_ttarget);
  tmp_tsource = new irtkGreyImage(*_tsource);

  // Swap source and target with temp space copies
  swap(tmp_ttarget, _ttarget);
  swap(tmp_tsource, _tsource);

  // Blur images if necessary
  if (_TargetBlurring[level] > 0) {
    cout << "Blurring taggedtarget ... "; cout.flush();
    irtkGaussianBlurringWithPadding<irtkGreyPixel> blurring(_TargetBlurring[level], _TargetPadding);
    blurring.SetInput (_ttarget);
    blurring.SetOutput(_ttarget);
    blurring.Run();
    cout << "done" << endl;
  }

  if (_SourceBlurring[level] > 0) {
    cout << "Blurring taggedsource ... "; cout.flush();
    irtkGaussianBlurring<irtkGreyPixel> blurring(_SourceBlurring[level]);
    blurring.SetInput (_tsource);
    blurring.SetOutput(_tsource);
    blurring.Run();
    cout << "done" << endl;
  }

  _ttarget->GetPixelSize(&dx, &dy, &dz);
  temp = fabs(_TargetResolution[0][0]-dx) + fabs(_TargetResolution[0][1]-dy) + fabs(_TargetResolution[0][2]-dz);

  if (level > 0 || temp > 0.000001) {
    cout << "Resampling tagged target ... "; cout.flush();
    // Create resampling filter
    irtkResamplingWithPadding<irtkGreyPixel> resample(_TargetResolution[level][0],
        _TargetResolution[level][1],
        _TargetResolution[level][2],
        _TargetPadding);
    resample.SetInput (_ttarget);
    resample.SetOutput(_ttarget);
    resample.Run();
    cout << "done" << endl;
  }

  _tsource->GetPixelSize(&dx, &dy, &dz);
  temp = fabs(_SourceResolution[0][0]-dx) + fabs(_SourceResolution[0][1]-dy) + fabs(_SourceResolution[0][2]-dz);

  if (level > 0 || temp > 0.000001) {
    cout << "Resampling tagged source ... "; cout.flush();
    // Create resampling filter
    irtkResamplingWithPadding<irtkGreyPixel> resample(_SourceResolution[level][0],
        _SourceResolution[level][1],
        _SourceResolution[level][2], MIN_GREY);

    resample.SetInput (_tsource);
    resample.SetOutput(_tsource);
    resample.Run();
    cout << "done" << endl;
  }
  
  CorrespondTaggedImages(_target,_source);
  // Find out the min and max values in target image, ignoring padding
  ttarget_max = MIN_GREY;
  ttarget_min = MAX_GREY;
  for (t = 0; t < _ttarget->GetT(); t++) {
    for (k = 0; k < _ttarget->GetZ(); k++) {
      for (j = 0; j < _ttarget->GetY(); j++) {
        for (i = 0; i < _ttarget->GetX(); i++) {
          if (_ttarget->Get(i, j, k, t) > _TargetPadding) {
            if (_ttarget->Get(i, j, k, t) > ttarget_max)
              ttarget_max = _ttarget->Get(i, j, k, t);
            if (_ttarget->Get(i, j, k, t) < ttarget_min)
              ttarget_min = _ttarget->Get(i, j, k, t);
          } else {
            _ttarget->Put(i, j, k, t, _TargetPadding);
          }
        }
      }
    }
  }

  // Find out the min and max values in source image, ignoring padding
  tsource_max = MIN_GREY;
  tsource_min = MAX_GREY;
  for (t = 0; t < _tsource->GetT(); t++) {
    for (k = 0; k < _tsource->GetZ(); k++) {
      for (j = 0; j < _tsource->GetY(); j++) {
        for (i = 0; i < _tsource->GetX(); i++) {
          if (_tsource->Get(i, j, k, t) > tsource_max)
            tsource_max = _tsource->Get(i, j, k, t);
          if (_tsource->Get(i, j, k, t) < tsource_min)
            tsource_min = _tsource->Get(i, j, k, t);
        }
      }
    }
  }

  // Check whether dynamic range of data is not to large
  if (ttarget_max - ttarget_min > MAX_GREY) {
    cerr << this->NameOfClass()
         << "::Initialize: Dynamic range of taggedtarget is too large" << endl;
    exit(1);
  } else {
    for (t = 0; t < _ttarget->GetT(); t++) {
      for (k = 0; k < _ttarget->GetZ(); k++) {
        for (j = 0; j < _ttarget->GetY(); j++) {
          for (i = 0; i < _ttarget->GetX(); i++) {
            if (_ttarget->Get(i, j, k, t) > _TargetPadding) {
              _ttarget->Put(i, j, k, t, _ttarget->Get(i, j, k, t) - ttarget_min);
            } else {
              _ttarget->Put(i, j, k, t, -1);
            }
          }
        }
      }
    }
  }

  if ((_SimilarityMeasure == SSD) || (_SimilarityMeasure == CC) ||
      (_SimilarityMeasure == LC)  || (_SimilarityMeasure == K) || (_SimilarityMeasure == ML)) {
    if (tsource_max - ttarget_min > MAX_GREY) {
      cerr << this->NameOfClass()
           << "::Initialize: Dynamic range of tagged source is too large" << endl;
      exit(1);
    } else {
      for (t = 0; t < _tsource->GetT(); t++) {
        for (k = 0; k < _tsource->GetZ(); k++) {
          for (j = 0; j < _tsource->GetY(); j++) {
            for (i = 0; i < _tsource->GetX(); i++) {
              _tsource->Put(i, j, k, t, _tsource->Get(i, j, k, t) - ttarget_min);
            }
          }
        }
      }
    }
  } else {
    if (tsource_max - tsource_min > MAX_GREY) {
      cerr << this->NameOfClass()
           << "::Initialize: Dynamic range of tagged source is too large" << endl;
      exit(1);
    } else {
      for (t = 0; t < _tsource->GetT(); t++) {
        for (k = 0; k < _tsource->GetZ(); k++) {
          for (j = 0; j < _tsource->GetY(); j++) {
            for (i = 0; i < _tsource->GetX(); i++) {
              _tsource->Put(i, j, k, t, _tsource->Get(i, j, k, t) - tsource_min);
            }
          }
        }
      }
    }
  }
  // Setup the interpolator
  _tinterpolator = irtkInterpolateImageFunction::New(_InterpolationMode, _tsource);

  // Setup interpolation for the source image
  _tinterpolator->SetInput(_tsource);
  _tinterpolator->Initialize();

  // Tell optimizer which transformation to optimize
  _optimizer->SetTransformation(_affd);

  // Allocate memory for metric
  _tmpMetricA = irtkSimilarityMetric::New(_metric);
  _tmpMetricB = irtkSimilarityMetric::New(_metric);

  _tmpImage = new irtkGreyImage(_target->GetX(),
                                _target->GetY(),
                                _target->GetZ(),
                                _target->GetT());
  _ttmpImage = new irtkGreyImage(_ttarget->GetX(),
                                _ttarget->GetY(),
                                _ttarget->GetZ(),
                                _ttarget->GetT());
  _weight = new irtkGreyImage(_ttarget->GetX(),
                                _ttarget->GetY(),
                                _ttarget->GetZ(),
                                _ttarget->GetT());

  this->EvaluateWeight(_weight,_target);
  char buffer[255];
  sprintf(buffer, "d:\\%d.gipl", level);
  _weight->Write(buffer);

  n = _target->GetNumberOfVoxels() * 3 / _target->GetT();

  // Allocate memory for lookup table for single-level FFD
  _affdLookupTable  = new float[n];

  // Allocate memory for lookup table for multi-level FFD
  _mffdLookupTable  = new float[n];

  // Allocate memory for local lookuptable
  _localLookupTable = new float [FFDLOOKUPTABLESIZE];

  for (i = 0; i < FFDLOOKUPTABLESIZE; i++) {
    u = i / (FFDLOOKUPTABLESIZE / 4.0);
    j = (int)floor(u);
    u = u - j;
    _localLookupTable[i] = _affd->B(j, 1-u);
  }

  // Initialize lookup table for multi-level FFD (this is done only once)
  ptr = _mffdLookupTable;
  for (k = 0; k < _target->GetZ(); k++) {
    for (j = 0; j < _target->GetY(); j++) {
      for (i = 0; i < _target->GetX(); i++) {
        x = i;
        y = j;
        z = k;
        _target->ImageToWorld(x, y, z);
        _mffd->Transform(x, y, z);
        ptr[0] = x;
        ptr[1] = y;
        ptr[2] = z;
        ptr += 3;
      }
    }
  }

  // Padding of FFD
  irtkPadding(*tmp_target, this->_TargetPadding, _affd);
  // Padding of FFD
  irtkPadding(*tmp_ttarget, this->_TargetPadding, _affd);

  // Register in the x-direction only
  if (_Mode == RegisterX) {
    for (i = 0; i < _affd->GetX(); i++) {
      for (j = 0; j < _affd->GetY(); j++) {
        for (k = 0; k < _affd->GetZ(); k++) {
          _Status sx, sy, sz;
          _affd->GetStatus(i, j, k, sx, sy, sz);
          _affd->PutStatus(i, j, k, sx, _Passive, _Passive);
        }
      }
    }
  }

  // Register in the y-direction only
  if (_Mode == RegisterY) {
    for (i = 0; i < _affd->GetX(); i++) {
      for (j = 0; j < _affd->GetY(); j++) {
        for (k = 0; k < _affd->GetZ(); k++) {
          _Status sx, sy, sz;
          _affd->GetStatus(i, j, k, sx, sy, sz);
          _affd->PutStatus(i, j, k, _Passive, sy, _Passive);
        }
      }
    }
  }

  // Register in the x- and y-direction only
  if (_Mode == RegisterXY) {
    for (i = 0; i < _affd->GetX(); i++) {
      for (j = 0; j < _affd->GetY(); j++) {
        for (k = 0; k < _affd->GetZ(); k++) {
          _Status sx, sy, sz;
          _affd->GetStatus(i, j, k, sx, sy, sz);
          _affd->PutStatus(i, j, k, sx, sy, _Passive);
        }
      }
    }
  }
   // Print some debugging information
  cout << "Tagged Target image (reference)" << endl;
  _ttarget->Print();
  cout << "Range is from " << ttarget_min << " to " << ttarget_max << endl;

  cout << "Tagged Source image (transform)" << endl;
  _tsource->Print();
  cout << "Range is from " << tsource_min << " to " << tsource_max << endl;


}

void irtkCardiacImageFreeFormRegistration::Finalize()
{
  // Print debugging information
  this->Debug("irtkCardiacImageFreeFormRegistration::Finalize");

  // Push local transformation back on transformation stack
  _mffd->PushLocalTransformation(_affd);

  // Finalize base class
  this->irtkImageRegistration::Finalize();

  delete []_localLookupTable;
}

void irtkCardiacImageFreeFormRegistration::Finalize(int level)
{
  // Print debugging information
  this->Debug("irtkCardiacImageFreeFormRegistration::Finalize(int)");

  // Finalize base class
  this->irtkImageRegistration::Finalize(level);

  // Check if we are not at the lowest level of resolution
  if (level != 0) {
    if (this->_Subdivision == True) {
      _affd->Subdivide();
    } else {
      // Push local transformation back on transformation stack
      _mffd->PushLocalTransformation(_affd);

      // Create new FFD
      _affd = new irtkBSplineFreeFormTransformation(*_target,
          this->_DX / pow(2.0, this->_NumberOfLevels-level),
          this->_DY / pow(2.0, this->_NumberOfLevels-level),
          this->_DZ / pow(2.0, this->_NumberOfLevels-level));
    }
  }

  // Swap source and target back with temp space copies (see Initialize)
  swap(tmp_ttarget, _ttarget);
  swap(tmp_tsource, _tsource);

  delete tmp_ttarget;
  delete tmp_tsource;
  delete _tmpImage;
  delete _ttmpImage;
  delete _weight;
  delete _tmpMetricA;
  delete _tmpMetricB;
  delete []_affdLookupTable;
  delete []_mffdLookupTable;
  delete _tinterpolator;
}

void irtkCardiacImageFreeFormRegistration::UpdateLUT()
{
  int i, j, k;
  double x, y, z;
  float *ptr2mffd;
  float *ptr2affd;

  // Print debugging information
  this->Debug("irtkCardiacImageFreeFormRegistration::UpdateLUT");

  ptr2affd = _affdLookupTable;
  ptr2mffd = _mffdLookupTable;
  for (k = 0; k < _target->GetZ(); k++) {
    for (j = 0; j < _target->GetY(); j++) {
      for (i = 0; i < _target->GetX(); i++) {
        x = i;
        y = j;
        z = k;
        _target->ImageToWorld(x, y, z);
        _affd->LocalDisplacement(x, y, z);
        ptr2affd[0] = x + ptr2mffd[0];
        ptr2affd[1] = y + ptr2mffd[1];
        ptr2affd[2] = z + ptr2mffd[2];
        ptr2mffd += 3;
        ptr2affd++;
        ptr2affd++;
        ptr2affd++;
      }
    }
  }
}

double irtkCardiacImageFreeFormRegistration::SmoothnessPenalty()
{
  int i, j, k;
  double x, y, z, penalty;

  penalty = 0;
  for (k = 0; k < _affd->GetZ(); k++) {
    for (j = 0; j < _affd->GetY(); j++) {
      for (i = 0; i < _affd->GetX(); i++) {
        x = i;
        y = j;
        z = k;
        _affd->LatticeToWorld(x, y, z);
        penalty += _affd->Bending(x, y, z);
      }
    }
  }
  return -penalty / _affd->NumberOfDOFs();
}

double irtkCardiacImageFreeFormRegistration::SmoothnessPenalty(int index)
{
  int i, j, k;
  double x, y, z;

  _affd->IndexToLattice(index, i, j, k);
  x = i;
  y = j;
  z = k;
  _affd->LatticeToWorld(x, y, z);
  return -_affd->Bending(x, y, z);
}

double irtkCardiacImageFreeFormRegistration::VolumePreservationPenalty()
{
  int i, j, k;
  double x, y, z, penalty;

  penalty = 0;
  for (k = 0; k < _affd->GetZ(); k++) {
    for (j = 0; j < _affd->GetY(); j++) {
      for (i = 0; i < _affd->GetZ(); i++) {
        x = i;
        y = j;
        z = k;
        _affd->LatticeToWorld(x, y, z);
        // Torsten Rohlfing et al. MICCAI'01 (w/o scaling correction):
        penalty += fabs(log(_affd->irtkTransformation::Jacobian(x, y, z)));
      }
    }
  }

  // Normalize sum by number of DOFs
  return penalty / (double) _affd->NumberOfDOFs();
}

double irtkCardiacImageFreeFormRegistration::VolumePreservationPenalty(int index)
{
  int i, j, k;
  double x, y, z;

  _affd->IndexToLattice(index, i, j, k);
  x = i;
  y = j;
  z = k;
  _affd->LatticeToWorld(x, y, z);
  return fabs(log(_affd->irtkTransformation::Jacobian(x, y, z)));
}


double irtkCardiacImageFreeFormRegistration::TopologyPreservationPenalty()
{
  int i, j, k;
  double x, y, z, jac, penalty;

  penalty = 0;
  for (k = 0; k < _affd->GetZ()-1; k++) {
    for (j = 0; j < _affd->GetY()-1; j++) {
      for (i = 0; i < _affd->GetZ()-1; i++) {
        x = i+0.5;
        y = j+0.5;
        z = k+0.5;
        _affd->LatticeToWorld(x, y, z);
        jac = _affd->irtkTransformation::Jacobian(x, y, z);
        if (jac < 0.3) {
          penalty += 10*jac*jac + 0.1/(jac*jac) - 2.0;
        }
      }
    }
  }
  return -penalty;
}

double irtkCardiacImageFreeFormRegistration::TopologyPreservationPenalty(int index)
{
  int i, j, k, l, m, n;
  double x, y, z, jac, penalty;

  penalty = 0;
  for (l = 0; l <= 1; l++) {
    for (m = 0; m <= 1; m++) {
      for (n = 0; n <= 1; n++) {
        _affd->IndexToLattice(index, i, j, k);
        x = i+l-0.5;
        y = j+m-0.5;
        z = k+n-0.5;
        _affd->LatticeToWorld(x, y, z);
        jac = _affd->irtkTransformation::Jacobian(x, y, z);
        if (jac < 0.3) {
          penalty += 10*jac*jac + 0.1/(jac*jac) - 2.0;
        }
      }
    }
  }
  return -penalty;
}

double irtkCardiacImageFreeFormRegistration::Evaluate()
{
#ifndef HAS_TBB
  // Image coordinates
  int i, j, k, t;
  // World coordinates
  double x, y, z;
  // Pointer to reference data
  irtkGreyPixel *ptr2target;
  irtkGreyPixel *ptr2ttarget;
  irtkGreyPixel *ptr2tmp;
  irtkGreyPixel *ptr2ttmp;
  irtkGreyPixel *ptr2weight;
  float *ptr;
#endif
     /*used to test if the correspond programe is correct or not.
	  *_target->Write("D:\\target.gipl");
	  *_source->Write("D:\\source.gipl");
	  *_ttarget->Write("D:\\ttarget.gipl");
	  *_tsource->Write("D:\\tsource.gipl");*/
  // Print debugging information
  this->Debug("irtkCardiacImageFreeFormRegistration::Evaluate");

  // Initialize metric
  _metric->Reset();

#ifdef HAS_TBB
  irtkMultiThreadedImageFreeFormRegistrationEvaluate evaluate(this);
  parallel_reduce(blocked_range<int>(0, _target->GetZ(), 1), evaluate);
#else

  // Loop over all voxels in the target (reference) volume
  ptr2target = _target->GetPointerToVoxels();
  ptr2ttarget= _ttarget->GetPointerToVoxels();
  ptr2tmp    = _tmpImage->GetPointerToVoxels();
  ptr2ttmp   = _ttmpImage->GetPointerToVoxels();
  ptr2weight = _weight->GetPointerToVoxels();
  for (t = 0; t < _target->GetT(); t++) {
    ptr        = _mffdLookupTable;
    for (k = 0; k < _target->GetZ(); k++) {
      for (j = 0; j < _target->GetY(); j++) {
        for (i = 0; i < _target->GetX(); i++) {
          // Check whether reference point is valid
          if (*ptr2target >= 0 || *ptr2ttarget >= 0) {
            x = i;
            y = j;
            z = k;
            _target->ImageToWorld(x, y, z);
            _affd->LocalDisplacement(x, y, z);
            x += ptr[0];
            y += ptr[1];
            z += ptr[2];
            _source->WorldToImage(x, y, z);
            // Check whether transformed point is inside volume
            if ((x > _source_x1) && (x < _source_x2) &&
                (y > _source_y1) && (y < _source_y2) &&
                (z > _source_z1) && (z < _source_z2)) {
              // Add sample to metric
				if (*ptr2target >= 0) {
				  *ptr2tmp =  round(_interpolator->EvaluateInside(x, y, z, t));
				  _metric->Add(*ptr2target, *ptr2tmp, *ptr2weight);
			    }else
					*ptr2tmp = -1;
				if (*ptr2ttarget >= 0){
			      *ptr2ttmp =  round(_tinterpolator->EvaluateInside(x, y, z, t));
				  _metric->Add(*ptr2ttarget, *ptr2ttmp,(MaxWeight - *ptr2weight));
			    }
			    else
				  *ptr2ttmp = -1;
            } else {
              *ptr2tmp = -1;
			  *ptr2ttmp = -1;
            }
          }
          // Increment pointers to next voxel
          ptr2tmp++;
          ptr2target++;
		  ptr2ttmp++;
          ptr2ttarget++;
		  ptr2weight++;
          ptr += 3;
        }
      }
    }
  }

#endif
  double similarity = _metric->Evaluate();

  // Add penalty for smoothness
  if (this->_Lambda1 > 0) {
    similarity += this->_Lambda1*this->SmoothnessPenalty();
  }
  // Add penalty for volume preservation
  if (this->_Lambda2 > 0) {
    similarity += this->_Lambda2*this->VolumePreservationPenalty();
  }
  // Add penalty for topology preservation
  if (this->_Lambda3 > 0) {
    similarity += this->_Lambda3*this->TopologyPreservationPenalty();
  }

  // Return similarity measure + penalty terms
  return similarity;
}

double irtkCardiacImageFreeFormRegistration::EvaluateDerivative(int index, double step)
{
  float *ptr;
  irtkPoint p1, p2;
  double bi, bj, bk, dx, dy, dz, p[3];
  int i, j, k, i1, i2, j1, j2, k1, k2, dim, t;
  irtkGreyPixel *ptr2target, *ptr2tmp;
  irtkGreyPixel *ptr2ttarget, *ptr2ttmp, *ptr2weight;
  irtkSimilarityMetric *tmpMetricA, *tmpMetricB;

  // Print debugging information
  this->Debug("irtkImageFreeFormRegistration::EvaluateDerivative(int, double)");

#ifdef HAS_TBB
  // Create similarity metric if necessary
  if (queue.pop_if_present(tmpMetricA) == false) {
    tmpMetricA = irtkSimilarityMetric::New(_metric);
  }
  // Create similarity metric if necessary
  if (queue.pop_if_present(tmpMetricB) == false) {
    tmpMetricB = irtkSimilarityMetric::New(_metric);
  }
#else
  tmpMetricA = _tmpMetricA;
  tmpMetricB = _tmpMetricB;
#endif

  // Initialize metrics for forward and backward derivative steps
  tmpMetricA->Reset(_metric);
  tmpMetricB->Reset(_metric);

  // Calculate bounding box of control point in world coordinates
  _affd->BoundingBox(index, p1, p2);
  _target->WorldToImage(p1);
  _target->WorldToImage(p2);

  // Calculate bounding box of control point in image coordinates
  _affd->BoundingBox(_target, index, i1, j1, k1, i2, j2, k2, 1.0 / _SpeedupFactor);

  // Calculate incremental changes in lattice coordinates when looping
  // over target
  dx = (FFDLOOKUPTABLESIZE-1)/(p2._x-p1._x);
  dy = (FFDLOOKUPTABLESIZE-1)/(p2._y-p1._y);
  dz = (FFDLOOKUPTABLESIZE-1)/(p2._z-p1._z);

  // Calculate whether this DOF corresponds to x, y or z-displacement
  dim = int(index / (_affd->GetX()*_affd->GetY()*_affd->GetZ()));

  // Loop over all voxels in the target (reference) volume
  for (t = 0; t < _target->GetT(); t++) {
    for (k = k1; k <= k2; k++) {
      bk = step * _localLookupTable[round(dz*(k-p1._z))];
      for (j = j1; j <= j2; j++) {
        bj = bk * _localLookupTable[round(dy*(j-p1._y))];
		ptr        = &(_affdLookupTable[3*_target->VoxelToIndex(i1, j, k)]); 
		ptr2target = _target->GetPointerToVoxels(i1, j, k, t);               
        ptr2tmp  = _tmpImage->GetPointerToVoxels(i1, j, k, t);
		ptr2ttarget = _ttarget->GetPointerToVoxels(i1, j, k, t);        
        ptr2ttmp  = _ttmpImage->GetPointerToVoxels(i1, j, k, t);
		ptr2weight = _weight->GetPointerToVoxels(i1, j, k, t);
        for (i = i1; i <= i2; i++) {

          // Check whether reference point is valid
          if (*ptr2target >= 0 || *ptr2ttarget >= 0) {
            bi = bj * _localLookupTable[round(dx*(i-p1._x))];

            // Delete old samples from both metrics
            if (*ptr2tmp != -1 && *ptr2target >= 0) {
              tmpMetricA->Delete(*ptr2target, *ptr2tmp, *ptr2weight);
              tmpMetricB->Delete(*ptr2target, *ptr2tmp, *ptr2weight);
            }

			// Delete old samples from both metrics
            if (*ptr2ttmp != -1 && *ptr2ttarget >= 0) {
              tmpMetricA->Delete(*ptr2ttarget, *ptr2ttmp, (MaxWeight - *ptr2weight));
              tmpMetricB->Delete(*ptr2ttarget, *ptr2ttmp, (MaxWeight - *ptr2weight));
            }

            p[0] = ptr[0];
            p[1] = ptr[1];
            p[2] = ptr[2];
            p[dim] += bi;

            // Convert transformed point to image coordinates
            _source->WorldToImage(p[0], p[1], p[2]);

            // Check whether transformed point is inside volume
            if ((p[0] > _source_x1) && (p[0] < _source_x2) &&
                (p[1] > _source_y1) && (p[1] < _source_y2) &&
                (p[2] > _source_z1) && (p[2] < _source_z2)) {
              // Add sample to metric
			  if (*ptr2target >= 0)
              tmpMetricA->Add(*ptr2target, round(_interpolator->EvaluateInside(p[0], p[1], p[2], t)), *ptr2weight);
			  if (*ptr2ttarget >= 0)
              // Add sample to metric
              tmpMetricA->Add(*ptr2ttarget, round(_tinterpolator->EvaluateInside(p[0], p[1], p[2], t)), (MaxWeight - *ptr2weight));
            }
            p[0] = ptr[0];
            p[1] = ptr[1];
            p[2] = ptr[2];
            p[dim] -= bi;

            // Convert transformed point to image coordinates
            _source->WorldToImage(p[0], p[1], p[2]);

            // Check whether transformed point is inside volume
            if ((p[0] > _source_x1) && (p[0] < _source_x2) &&
                (p[1] > _source_y1) && (p[1] < _source_y2) &&
                (p[2] > _source_z1) && (p[2] < _source_z2)) {

              if (*ptr2target >= 0)
              tmpMetricB->Add(*ptr2target, round(_interpolator->EvaluateInside(p[0], p[1], p[2], t)), *ptr2weight);
			  if (*ptr2ttarget >= 0)
              // Add sample to metric
              tmpMetricB->Add(*ptr2ttarget, round(_tinterpolator->EvaluateInside(p[0], p[1], p[2], t)), (MaxWeight - *ptr2weight));
            }
          }

          // Increment pointers to next voxel
          ptr2target++;
          ptr2tmp++;
		  ptr2ttarget++;
          ptr2ttmp++;
          ptr += 3;
        }
      }
    }
  }

  // Save value of DOF for which we calculate the derivative
  double dof = _affd->Get(index);

  // Evaluate similarity measure
  double similarityA = tmpMetricA->Evaluate();

  // Add penalties
  _affd->Put(index, dof + step);

  // Smoothness
  if (this->_Lambda1 > 0) {
    similarityA += this->_Lambda1*this->SmoothnessPenalty(index);
  }
  // Volume preservation
  if (this->_Lambda2 > 0) {
    similarityA += this->_Lambda2*this->VolumePreservationPenalty(index);
  }
  // Topology preservation
  if (this->_Lambda3 > 0) {
    similarityA += this->_Lambda3*this->TopologyPreservationPenalty(index);
  }

  // Evaluate similarity measure
  double similarityB = tmpMetricB->Evaluate();
  // Add penalties
  _affd->Put(index, dof - step);

  // Smoothness
  if (this->_Lambda1 > 0) {
    similarityB += this->_Lambda1*this->SmoothnessPenalty(index);
  }
  // Volume preservation
  if (this->_Lambda2 > 0) {
    similarityB += this->_Lambda2*this->VolumePreservationPenalty(index);
  }
  // Topology preservation
  if (this->_Lambda3 > 0) {
    similarityB += this->_Lambda3*this->TopologyPreservationPenalty(index);
  }

  // Restore value of DOF for which we calculate the derivative
  _affd->Put(index, dof);

#ifdef HAS_TBB
  queue.push(tmpMetricA);
  queue.push(tmpMetricB);
#endif

  return similarityA - similarityB;
}

double irtkCardiacImageFreeFormRegistration::EvaluateGradient(float step, float *dx)
{
  int i;
  double norm;

  // Update lookup table
  this->UpdateLUT();

#ifdef HAS_TBB
  parallel_for(blocked_range<int>(0, _affd->NumberOfDOFs(), 1), irtkMultiThreadedImageFreeFormRegistrationEvaluateGradient(this, dx, step));
#else
  for (i = 0; i < _affd->NumberOfDOFs(); i++) {
    if (_affd->irtkTransformation::GetStatus(i) == _Active) {
      dx[i] = this->EvaluateDerivative(i, step);
    } else {
      dx[i] = 0;
    }
  }
#endif

  // Calculate norm of vector
  norm = 0;
  for (i = 0; i < _affd->NumberOfDOFs(); i++) {
    norm += dx[i] * dx[i];
  }

  // Normalize vector
  norm = sqrt(norm);
  if (norm > 0) {
    for (i = 0; i < _affd->NumberOfDOFs(); i++) {
      dx[i] /= norm;
    }
  } else {
    for (i = 0; i < _affd->NumberOfDOFs(); i++) {
      dx[i] = 0;
    }
  }

  return norm;
}

Bool irtkCardiacImageFreeFormRegistration::Read(char *buffer1, char *buffer2, int &level)
{
  int ok = False;

  if (strstr(buffer1, "Speedup factor ") != NULL) {
    this->_SpeedupFactor = atof(buffer2);
    cout << "Speedup factor is ... " << this->_SpeedupFactor << endl;
    ok = True;
  }
  if ((strstr(buffer1, "Lambda ") != NULL) ||
      (strstr(buffer1, "Lambda1") != NULL)) {
    this->_Lambda1 = atof(buffer2);
    cout << "Lambda 1 is ... " << this->_Lambda1 << endl;
    ok = True;
  }
  if (strstr(buffer1, "Lambda2") != NULL) {
    this->_Lambda2 = atof(buffer2);
    cout << "Lambda 2 is ... " << this->_Lambda2 << endl;
    ok = True;
  }
  if (strstr(buffer1, "Lambda3") != NULL) {
    this->_Lambda3 = atof(buffer2);
    cout << "Lambda 3 is ... " << this->_Lambda3 << endl;
    ok = True;
  }
  if (strstr(buffer1, "UntaggedWeight") != NULL) {
    this->_Weight = atof(buffer2);
    cout << "UntaggedWeight is ... " << this->_Weight << endl;
    ok = True;
  }
  if (strstr(buffer1, "TaggedWeight") != NULL) {
    this->_Tweight = atof(buffer2);
    cout << "TaggedWeight is ... " << this->_Tweight << endl;
    ok = True;
  }
  if (strstr(buffer1, "Thresholdmax") != NULL) {
    this->_Thresholdmax = atof(buffer2);
    cout << "Thresholdmax is ... " << this->_Thresholdmax << endl;
    ok = True;
  }
  if (strstr(buffer1, "Thresholdmin") != NULL) {
    this->_Thresholdmin = atof(buffer2);
    cout << "Thresholdmin is ... " << this->_Thresholdmin << endl;
    ok = True;
  }
  if (strstr(buffer1, "Control point spacing in X") != NULL) {
    this->_DX = atof(buffer2);
    cout << "Control point spacing in X is ... " << this->_DX << endl;
    ok = True;
  }
  if (strstr(buffer1, "Control point spacing in Y") != NULL) {
    this->_DY = atof(buffer2);
    cout << "Control point spacing in Y is ... " << this->_DY << endl;
    ok = True;
  }
  if (strstr(buffer1, "Control point spacing in Z") != NULL) {
    this->_DZ = atof(buffer2);
    cout << "Control point spacing in Z is ... " << this->_DZ << endl;
    ok = True;
  }
  if (strstr(buffer1, "Subdivision") != NULL) {
    if ((strcmp(buffer2, "False") == 0) || (strcmp(buffer2, "No") == 0)) {
      this->_Subdivision = False;
      cout << "Subdivision is ... false" << endl;
    } else {
      if ((strcmp(buffer2, "True") == 0) || (strcmp(buffer2, "Yes") == 0)) {
        this->_Subdivision = True;
        cout << "Subdivision is ... true" << endl;
      } else {
        cerr << "Can't read boolean value = " << buffer2 << endl;
        exit(1);
      }
    }
    ok = True;
  }

  if (ok == False) {
    return this->irtkImageRegistration::Read(buffer1, buffer2, level);
  } else {
    return ok;
  }
}

void irtkCardiacImageFreeFormRegistration::SetInput(irtkGreyImage *target, irtkGreyImage *source, irtkGreyImage *ttarget, irtkGreyImage *tsource)
{
  _ttarget = ttarget;
  _tsource = tsource;
  this->irtkImageRegistration::SetInput(target, source);
}

void irtkCardiacImageFreeFormRegistration::CorrespondTaggedImages (irtkGreyImage *target, irtkGreyImage *source)
{
  int i,j,k,t;
  double x,y,z;
  irtkGreyPixel *ptr2tmp;
  // Copy source and target to temp space
  irtkGreyImage *tmptarget = new irtkGreyImage(*target);
  irtkGreyImage *tmpsource = new irtkGreyImage(*source);

  // Setup the interpolator
  irtkInterpolateImageFunction *tinterpolator = irtkInterpolateImageFunction::New(_InterpolationMode, _ttarget);

  // Setup interpolation for the source image
  tinterpolator->SetInput(_ttarget);
  tinterpolator->Initialize();

  // Setup the interpolator
  irtkInterpolateImageFunction *sinterpolator = irtkInterpolateImageFunction::New(_InterpolationMode, _tsource);

  // Setup interpolation for the source image
  sinterpolator->SetInput(_tsource);
  sinterpolator->Initialize();

  
  int ttarget_x,ttarget_y,ttarget_z;
  int tsource_x,tsource_y,tsource_z;
  // Compute domain on which the linear interpolation is defined
  ttarget_x = _ttarget->GetX() - 1;
  ttarget_y = _ttarget->GetY() - 1;
  ttarget_z = _ttarget->GetZ() - 1;

  tsource_x = _tsource->GetX() - 1;
  tsource_y = _tsource->GetY() - 1;
  tsource_z = _tsource->GetZ() - 1;

  for (t = 0; t < target->GetT(); t++) {
    for (k = 0; k < target->GetZ(); k++) {
      for (j = 0; j < target->GetY(); j++) {
        for (i = 0; i < target->GetX(); i++) {
            ptr2tmp = tmptarget->GetPointerToVoxels(i, j, k, t);
			x = i;
            y = j;
            z = k;
            target->ImageToWorld(x, y, z);
            _ttarget->WorldToImage(x, y, z);
            // Check whether transformed point is inside volume
            if ((x > 0) && (x < ttarget_x) &&
                (y > 0) && (y < ttarget_y) &&
                (z > 0) && (z < ttarget_z)) {
              // Add sample to metric
			  
				*ptr2tmp =  round(tinterpolator->EvaluateInside(x, y, z, t));
            } else {
              *ptr2tmp = -1;
            }
          }
        }
      }
    }

  for (t = 0; t < source->GetT(); t++) {
    for (k = 0; k < source->GetZ(); k++) {
      for (j = 0; j < source->GetY(); j++) {
        for (i = 0; i < source->GetX(); i++) {
            ptr2tmp = tmpsource->GetPointerToVoxels(i, j, k, t);
			x = i;
            y = j;
            z = k;
            source->ImageToWorld(x, y, z);
            _tsource->WorldToImage(x, y, z);
            // Check whether transformed point is inside volume
            if ((x > 0) && (x < ttarget_x) &&
                (y > 0) && (y < ttarget_y) &&
                (z > 0) && (z < ttarget_z)) {
              // Add sample to metric
			  
				*ptr2tmp =  round(sinterpolator->EvaluateInside(x, y, z, t));
            } else {
              *ptr2tmp = -1;
            }
          }
        }
      }
    }


  // Swap source and target with temp space copies
 swap(tmptarget, _ttarget);
 swap(tmpsource, _tsource);
}

int irtkCardiacImageFreeFormRegistration::WeightFunction (double value){

	double weight = 0;
	if (value >= _Thresholdmin && value <= _Thresholdmax){
      weight = MaxWeight - _Tweight*MaxWeight/(_Tweight+_Weight);
	}else{
      weight = MaxWeight;
	}
	return round(weight);
}

void irtkCardiacImageFreeFormRegistration::EvaluateWeight (irtkGreyImage *weight, irtkGreyImage *target)
{
  int i,j,k,t;
  double x,y,z;
  irtkGreyPixel *ptr2weight,*ptr2target;
  //evaluate weighting function

  //caculate the weight
  for (t = 0; t < target->GetT(); t++) {
    for (k = 0; k < target->GetZ(); k++) {
      for (j = 0; j < target->GetY(); j++) {
        for (i = 0; i < target->GetX(); i++) {
            ptr2weight = weight->GetPointerToVoxels(i, j, k, t);
			ptr2target = target->GetPointerToVoxels(i, j, k, t);
			if (*ptr2target >=0)
				*ptr2weight = WeightFunction(*ptr2target);
			else
				*ptr2weight = -1;
			
          }
        }
      }
    }
}

void irtkCardiacImageFreeFormRegistration::Write(ostream &to)
{
  to << "\n#\n# Non-rigid registration parameters\n#\n\n";
  to << "Lambda1                           = " << this->_Lambda1 << endl;
  to << "Lambda2                           = " << this->_Lambda2 << endl;
  to << "Lambda3                           = " << this->_Lambda3 << endl;
  to << "UntaggedWeight                    = " << this->_Weight << endl;
  to << "TaggedWeight                      = " << this->_Tweight << endl;
  to << "Thresholdmax                      = " << this->_Thresholdmax << endl;
  to << "Thresholdmin                      = " << this->_Thresholdmin << endl;
  to << "Control point spacing in X        = " << this->_DX << endl;
  to << "Control point spacing in Y        = " << this->_DY << endl;
  to << "Control point spacing in Z        = " << this->_DZ << endl;
  if (_Subdivision == True) {
    to << "Subdivision                       = True" << endl;
  } else {
    to << "Subdivision                       = False" << endl;
  }
  to << "Speedup factor                    = " << this->_SpeedupFactor << endl;

  this->irtkImageRegistration::Write(to);
}
