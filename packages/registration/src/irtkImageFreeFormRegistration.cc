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

// Used as temporary memory for transformed intensities
irtkGreyImage *_tmpImage;

// The original target and source images
extern irtkGreyImage *tmp_target, *tmp_source;

#include <irtkMultiThreadedImageFreeFormRegistration.h>

irtkImageFreeFormRegistration::irtkImageFreeFormRegistration()
{
  // Print debugging information
  this->Debug("irtkImageFreeFormRegistration::irtkImageFreeFormRegistration");

  // Default optimization
  _OptimizationMethod = GradientDescent;

  // Default speedup factor
  _SpeedupFactor = 1;

  // Default parameters for non-rigid registration
  _Lambda1     = 0;
  _Lambda2     = 0;
  _Lambda3     = 0;
  _DX          = 20;
  _DY          = 20;
  _DZ          = 20;
  _Subdivision = true;
  _Mode        = RegisterXYZ;
  _MFFDMode    = true;
}

void irtkImageFreeFormRegistration::GuessParameter()
{
  int i;
  double xsize, ysize, zsize, spacing;

  if ((_target == NULL) || (_source == NULL)) {
    cerr << "irtkImageFreeFormRegistration::GuessParameter: Target and source image not found" << endl;
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
  _DX                 =_target->GetX() * spacing / 10.0;
  _DY                 =_target->GetX() * spacing / 10.0;
  _DZ                 =_target->GetX() * spacing / 10.0;
  _Subdivision        = true;

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

void irtkImageFreeFormRegistration::Initialize()
{
  // Print debugging information
  this->Debug("irtkImageFreeFormRegistration::Initialize");

  // Initialize base class
  this->irtkImageRegistration::Initialize();
  int i, j;
  double u;

  // Pointer to multi-level FFD
  _mffd = (irtkMultiLevelFreeFormTransformation *)_transformation;

  // Create FFD
  if (_MFFDMode == false) {
    if (_mffd->NumberOfLevels() == 0) {
      _affd = new irtkBSplineFreeFormTransformation(*_target, this->_DX, this->_DY, this->_DZ);
    } else {
      _affd = (irtkBSplineFreeFormTransformation *)_mffd->PopLocalTransformation();
    }
  } else {
    _affd = new irtkBSplineFreeFormTransformation(*_target, this->_DX, this->_DY, this->_DZ);
  }

  // Initialize pointers
  _tmpImage         = NULL;
  _affdLookupTable  = NULL;
  _mffdLookupTable  = NULL;
    // Allocate memory for local lookuptable
  _localLookupTable = new float [FFDLOOKUPTABLESIZE];

  for (i = 0; i < FFDLOOKUPTABLESIZE; i++) {
    u = i / (FFDLOOKUPTABLESIZE / 4.0);
    j = (int)floor(u);
    u = u - j;
    _localLookupTable[i] = _affd->B(j, 1-u);
  }

}

void irtkImageFreeFormRegistration::Initialize(int level)
{
  int i, j, k, n;
  double x, y, z;
  float *ptr;

  // Print debugging information
  this->Debug("irtkImageFreeFormRegistration::Initialize(int)");

  // Initialize base class
  this->irtkImageRegistration::Initialize(level);

  // Tell optimizer which transformation to optimize
  _optimizer->SetTransformation(_affd);

  // Allocate memory for metric
  _tmpMetricA = irtkSimilarityMetric::New(_metric);
  _tmpMetricB = irtkSimilarityMetric::New(_metric);

  _tmpImage = new irtkGreyImage(_target->GetX(),
                                _target->GetY(),
                                _target->GetZ(),
                                _target->GetT());

  n = _target->GetNumberOfVoxels() * 3 / _target->GetT();

  // Allocate memory for lookup table for single-level FFD
  _affdLookupTable  = new float[n];

  // Allocate memory for lookup table for multi-level FFD
  _mffdLookupTable  = new float[n];

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

}

void irtkImageFreeFormRegistration::Finalize()
{
  // Print debugging information
  this->Debug("irtkImageFreeFormRegistration::Finalize");

  // Push local transformation back on transformation stack
  _mffd->PushLocalTransformation(_affd);

  // Finalize base class
  this->irtkImageRegistration::Finalize();

  delete []_localLookupTable;
}

void irtkImageFreeFormRegistration::Finalize(int level)
{
  // Print debugging information
  this->Debug("irtkImageFreeFormRegistration::Finalize(int)");

  // Finalize base class
  this->irtkImageRegistration::Finalize(level);

  // Check if we are not at the lowest level of resolution
  if (level != 0) {
    if (this->_Subdivision == true) {
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

  delete _tmpImage;
  delete _tmpMetricA;
  delete _tmpMetricB;
  delete []_affdLookupTable;
  delete []_mffdLookupTable;
}

void irtkImageFreeFormRegistration::UpdateLUT()
{
  int i, j, k;
  double x, y, z;
  float *ptr2mffd;
  float *ptr2affd;

  // Print debugging information
  this->Debug("irtkImageFreeFormRegistration::UpdateLUT");

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
        ptr2affd += 3;
      }
    }
  }
}

double irtkImageFreeFormRegistration::SmoothnessPenalty()
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

double irtkImageFreeFormRegistration::SmoothnessPenalty(int index)
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

double irtkImageFreeFormRegistration::VolumePreservationPenalty()
{
  int i, j, k;
  double x, y, z, penalty, jacobian;
  irtkMatrix jac,tmp_jac;

  penalty = 0;
  for (k = 0; k < _affd->GetZ(); k++) {
    for (j = 0; j < _affd->GetY(); j++) {
      for (i = 0; i < _affd->GetX(); i++) {
        x = i;
        y = j;
        z = k;
        _affd->LatticeToWorld(x, y, z);
        // Torsten Rohlfing et al. MICCAI'01 (w/o scaling correction):
		_affd->Jacobian(tmp_jac,x,y,z);
		// Calculate jacobian
		jac.Initialize(3, 3);
		_mffd->LocalJacobian(jac, x, y, z);

		// Subtract identity matrix
		tmp_jac(0, 0) = tmp_jac(0, 0) - 1;
		tmp_jac(1, 1) = tmp_jac(1, 1) - 1;
		tmp_jac(2, 2) = tmp_jac(2, 2) - 1;

		// Add jacobian
		jac += tmp_jac;
		// Determinant of Jacobian of deformation derivatives
		jacobian = (jac(0, 0)*jac(1, 1)*jac(2, 2) + jac(0, 1)*jac(1, 2)*jac(2, 0) +
			jac(0, 2)*jac(1, 0)*jac(2, 1) - jac(0, 2)*jac(1, 1)*jac(2, 0) -
			jac(0, 0)*jac(1, 2)*jac(2, 1) - jac(0, 1)*jac(1, 0)*jac(2, 2));
		if(jacobian < 0.0001) jacobian = 0.0001;
		penalty += fabs(log(jacobian));
      }
    }
  }

  // Normalize sum by number of DOFs
  return -penalty / (double) _affd->NumberOfDOFs();
}

double irtkImageFreeFormRegistration::VolumePreservationPenalty(int index)
{
    int i, j, k, i1, j1, k1, i2, j2, k2, count;
    double x, y, z, jacobian, penalty;

    _affd->IndexToLattice(index, i, j, k);
    penalty = 0;
    count = 0;
    k1 = (k-1)>0?(k-1):0;
    j1 = (j-1)>0?(j-1):0;
    i1 = (i-1)>0?(i-1):0;
    k2 = (k+2) < _affd->GetZ()? (k+2) : _affd->GetZ();
    j2 = (j+2) < _affd->GetY()? (j+2) : _affd->GetY();
    i2 = (i+2) < _affd->GetX()? (i+2) : _affd->GetX();
    for (k = k1; k < k2; k++) {
        for (j = j1; j < j2; j++) {
            for (i = i1; i < i2; i++) {
                x = i;
                y = j;
                z = k;
                _affd->LatticeToWorld(x, y, z);
                // Torsten Rohlfing et al. MICCAI'01 (w/o scaling correction):
                irtkMatrix jac,tmp_jac;
                _affd->Jacobian(tmp_jac,x,y,z);
                // Calculate jacobian
                jac.Initialize(3, 3);
                _mffd->LocalJacobian(jac, x, y, z);

                // Subtract identity matrix
                tmp_jac(0, 0) = tmp_jac(0, 0) - 1;
                tmp_jac(1, 1) = tmp_jac(1, 1) - 1;
                tmp_jac(2, 2) = tmp_jac(2, 2) - 1;

                // Add jacobian
                jac += tmp_jac;
                // Determinant of Jacobian of deformation derivatives
                jacobian = (jac(0, 0)*jac(1, 1)*jac(2, 2) + jac(0, 1)*jac(1, 2)*jac(2, 0) +
                    jac(0, 2)*jac(1, 0)*jac(2, 1) - jac(0, 2)*jac(1, 1)*jac(2, 0) -
                    jac(0, 0)*jac(1, 2)*jac(2, 1) - jac(0, 1)*jac(1, 0)*jac(2, 2));
                if(jacobian < 0.0001) jacobian = 0.0001;
                penalty += fabs(log(jacobian));
                count ++;
            }
        }
    }
    return -penalty/count;
}


double irtkImageFreeFormRegistration::TopologyPreservationPenalty()
{
  int i, j, k;
  double x, y, z, jac, penalty;

  penalty = 0;
  for (k = -1; k < _affd->GetZ(); k++) {
    for (j = -1; j < _affd->GetY(); j++) {
      for (i = -1; i < _affd->GetX(); i++) {
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

double irtkImageFreeFormRegistration::TopologyPreservationPenalty(int index)
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

double irtkImageFreeFormRegistration::Evaluate()
{
#ifndef HAS_TBB
  // Image coordinates
  int i, j, k, t;
  // World coordinates
  double x, y, z;
  // Pointer to reference data
  irtkGreyPixel *ptr2target;
  irtkGreyPixel *ptr2tmp;
  float *ptr;
#endif

  // Print debugging information
  this->Debug("irtkImageFreeFormRegistration::Evaluate");

  // Initialize metric
  _metric->Reset();

#ifdef HAS_TBB
  irtkMultiThreadedImageFreeFormRegistrationEvaluate evaluate(this);
  parallel_reduce(blocked_range<int>(0, _target->GetZ(), 1), evaluate);
#else

  // Loop over all voxels in the target (reference) volume
  ptr2target = _target->GetPointerToVoxels();
  ptr2tmp    = _tmpImage->GetPointerToVoxels();
  for (t = 0; t < _target->GetT(); t++) {
    ptr        = _mffdLookupTable;
    for (k = 0; k < _target->GetZ(); k++) {
      for (j = 0; j < _target->GetY(); j++) {
        for (i = 0; i < _target->GetX(); i++) {
          // Check whether reference point is valid
          if (*ptr2target >= 0) {
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
              *ptr2tmp =  round(_interpolator->EvaluateInside(x, y, z, t));
              _metric->Add(*ptr2target, *ptr2tmp);
            } else {
              *ptr2tmp = -1;
            }
          }
          // Increment pointers to next voxel
          ptr2tmp++;
          ptr2target++;
          ptr += 3;
        }
      }
    }
  }

#endif

  // Evaluate similarity measure
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

double irtkImageFreeFormRegistration::EvaluateDerivative(int index, double step)
{
  float *ptr;
  irtkPoint p1, p2;
  double bi, bj, bk, dx, dy, dz, p[3];
  int i, j, k, i1, i2, j1, j2, k1, k2, dim, t;
  irtkGreyPixel *ptr2target, *ptr2tmp;
  irtkSimilarityMetric *tmpMetricA, *tmpMetricB;

  // Print debugging information
  this->Debug("irtkImageFreeFormRegistration::EvaluateDerivative(int, double)");

#ifdef HAS_TBB
  // Create similarity metric if necessary
  if (sim_queue.pop_if_present(tmpMetricA) == false) {
    tmpMetricA = irtkSimilarityMetric::New(_metric);
  }
  // Create similarity metric if necessary
  if (sim_queue.pop_if_present(tmpMetricB) == false) {
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
        ptr2target = _target->GetPointerToVoxels(i1, j, k, t);
        ptr        = &(_affdLookupTable[3*_target->VoxelToIndex(i1, j, k)]);
        bj = bk * _localLookupTable[round(dy*(j-p1._y))];
        ptr2tmp  = _tmpImage->GetPointerToVoxels(i1, j, k, t);
        for (i = i1; i <= i2; i++) {

          // Check whether reference point is valid
          if (*ptr2target >= 0) {
            bi = bj * _localLookupTable[round(dx*(i-p1._x))];

            // Delete old samples from both metrics
            if (*ptr2tmp != -1) {
              tmpMetricA->Delete(*ptr2target, *ptr2tmp);
              tmpMetricB->Delete(*ptr2target, *ptr2tmp);
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
              tmpMetricA->Add(*ptr2target, round(_interpolator->EvaluateInside(p[0], p[1], p[2], t)));
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

              // Add sample to metric
              tmpMetricB->Add(*ptr2target, round(_interpolator->EvaluateInside(p[0], p[1], p[2], t)));
            }
          }

          // Increment pointers to next voxel
          ptr2target++;
          ptr2tmp++;
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
  sim_queue.push(tmpMetricA);
  sim_queue.push(tmpMetricB);
#endif

  return similarityA - similarityB;
}

double irtkImageFreeFormRegistration::EvaluateGradient(float step, float *dx)
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

bool irtkImageFreeFormRegistration::Read(char *buffer1, char *buffer2, int &level)
{
  int ok = false;

  if (strstr(buffer1, "Speedup factor ") != NULL) {
    this->_SpeedupFactor = atof(buffer2);
    cout << "Speedup factor is ... " << this->_SpeedupFactor << endl;
    ok = true;
  }
  if ((strstr(buffer1, "Lambda ") != NULL) ||
      (strstr(buffer1, "Lambda1") != NULL)) {
    this->_Lambda1 = atof(buffer2);
    cout << "Lambda 1 is ... " << this->_Lambda1 << endl;
    ok = true;
  }
  if (strstr(buffer1, "Lambda2") != NULL) {
    this->_Lambda2 = atof(buffer2);
    cout << "Lambda 2 is ... " << this->_Lambda2 << endl;
    ok = true;
  }
  if (strstr(buffer1, "Lambda3") != NULL) {
    this->_Lambda3 = atof(buffer2);
    cout << "Lambda 3 is ... " << this->_Lambda3 << endl;
    ok = true;
  }
  if (strstr(buffer1, "Control point spacing in X") != NULL) {
    this->_DX = atof(buffer2);
    cout << "Control point spacing in X is ... " << this->_DX << endl;
    ok = true;
  }
  if (strstr(buffer1, "Control point spacing in Y") != NULL) {
    this->_DY = atof(buffer2);
    cout << "Control point spacing in Y is ... " << this->_DY << endl;
    ok = true;
  }
  if (strstr(buffer1, "Control point spacing in Z") != NULL) {
    this->_DZ = atof(buffer2);
    cout << "Control point spacing in Z is ... " << this->_DZ << endl;
    ok = true;
  }
  if (strstr(buffer1, "Subdivision") != NULL) {
    if ((strcmp(buffer2, "False") == 0) || (strcmp(buffer2, "No") == 0)) {
      this->_Subdivision = false;
      cout << "Subdivision is ... false" << endl;
    } else {
      if ((strcmp(buffer2, "True") == 0) || (strcmp(buffer2, "Yes") == 0)) {
        this->_Subdivision = true;
        cout << "Subdivision is ... true" << endl;
      } else {
        cerr << "Can't read boolean value = " << buffer2 << endl;
        exit(1);
      }
    }
    ok = true;
  }
  if (strstr(buffer1, "MFFDMode") != NULL) {
    if ((strcmp(buffer2, "False") == 0) || (strcmp(buffer2, "No") == 0)) {
      this->_MFFDMode = false;
      cout << "MFFDMode is ... false" << endl;
    } else {
      if ((strcmp(buffer2, "True") == 0) || (strcmp(buffer2, "Yes") == 0)) {
        this->_MFFDMode = true;
        cout << "MFFDMode is ... true" << endl;
      } else {
        cerr << "Can't read boolean value = " << buffer2 << endl;
        exit(1);
      }
    }
    ok = true;
  }
  if (ok == false) {
    return this->irtkImageRegistration::Read(buffer1, buffer2, level);
  } else {
    return ok;
  }
}

void irtkImageFreeFormRegistration::Write(ostream &to)
{
  to << "\n#\n# Non-rigid registration parameters\n#\n\n";
  to << "Lambda1                           = " << this->_Lambda1 << endl;
  to << "Lambda2                           = " << this->_Lambda2 << endl;
  to << "Lambda3                           = " << this->_Lambda3 << endl;
  to << "Control point spacing in X        = " << this->_DX << endl;
  to << "Control point spacing in Y        = " << this->_DY << endl;
  to << "Control point spacing in Z        = " << this->_DZ << endl;
  if (_Subdivision == true) {
    to << "Subdivision                       = True" << endl;
  } else {
    to << "Subdivision                       = False" << endl;
  }
  to << "Speedup factor                    = " << this->_SpeedupFactor << endl;
  if (_MFFDMode == true) {
    to << "MFFDMode                          = True" << endl;
  } else {
    to << "MFFDMode                          = False" << endl;
  }

  this->irtkImageRegistration::Write(to);
}
