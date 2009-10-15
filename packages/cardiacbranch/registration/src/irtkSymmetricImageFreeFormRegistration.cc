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
irtkGreyImage *_tmpSourceImage, *_tmpTargetImage;

// The original target and source images
extern irtkGreyImage *tmp_target, *tmp_source;

irtkSymmetricImageFreeFormRegistration::irtkSymmetricImageFreeFormRegistration()
{
  // Print debugging information
  this->Debug("irtkSymmetricImageFreeFormRegistration::irtkSymmetricImageFreeFormRegistration");

  // Default optimization
  _OptimizationMethod = GradientDescent;

  // Default speedup factor
  _SpeedupFactor = 1;

  // Default parameters for non-rigid registration
  _Lambda1     = 0;
  _Lambda2     = 0;
  _Lambda3     = 0;
  _Lambda4     = 0;
  _DX          = 20;
  _DY          = 20;
  _DZ          = 20;
  _Subdivision = True;
}

void irtkSymmetricImageFreeFormRegistration::GuessParameter()
{
  int i;
  double xsize, ysize, zsize, spacing;

  if ((_target == NULL) || (_source == NULL)) {
    cerr << "irtkSymmetricImageFreeFormRegistration::GuessParameter: Target and source image not found" << endl;
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
  _Subdivision        = True;

  // Remaining parameters
  for (i = 0; i < _NumberOfLevels; i++) {
    _NumberOfIterations[i] = 10;
    _NumberOfSteps[i]      = 4;
    _LengthOfSteps[i]      = _DX / 8.0 * pow(2.0, i);
  }

  // Try to guess padding by looking at voxel values in all eight corners of the volume:
  // If all values are the same we assume that they correspond to the padding value
  _Padding = MIN_GREY;
  if ((_target->Get(_target->GetX()-1, 0, 0)                                 == _target->Get(0, 0, 0)) &&
      (_target->Get(0, _target->GetY()-1, 0)                                 == _target->Get(0, 0, 0)) &&
      (_target->Get(0, 0, _target->GetZ()-1)                                 == _target->Get(0, 0, 0)) &&
      (_target->Get(_target->GetX()-1, _target->GetY()-1, 0)                 == _target->Get(0, 0, 0)) &&
      (_target->Get(0, _target->GetY()-1, _target->GetZ()-1)                 == _target->Get(0, 0, 0)) &&
      (_target->Get(_target->GetX()-1, 0, _target->GetZ()-1)                 == _target->Get(0, 0, 0)) &&
      (_target->Get(_target->GetX()-1, _target->GetY()-1, _target->GetZ()-1) == _target->Get(0, 0, 0))) {
    _Padding = _target->Get(0, 0, 0);
  }
}

void irtkSymmetricImageFreeFormRegistration::Initialize()
{
  // Print debugging information
  this->Debug("irtkSymmetricImageFreeFormRegistration::Initialize");

  // Initialize base class
  this->irtkSymmetricImageRegistration::Initialize();

  // Pointer to multi-level FFD
  _mffd1 = (irtkMultiLevelFreeFormTransformation *)_transformation1;
  _mffd2 = (irtkMultiLevelFreeFormTransformation *)_transformation2;

  // Create first FFD
  if (_mffd1->NumberOfLevels() == 0) {
    _affd1 = new irtkBSplineFreeFormTransformation(*_target, this->_DX, this->_DY, this->_DZ);
  } else {
    _affd1 = (irtkBSplineFreeFormTransformation *)_mffd1->PopLocalTransformation();
  }

  // Create second FFD
  if (_mffd2->NumberOfLevels() == 0) {
    _affd2 = new irtkBSplineFreeFormTransformation(*_source, this->_DX, this->_DY, this->_DZ);
  } else {
    _affd2 = (irtkBSplineFreeFormTransformation *)_mffd2->PopLocalTransformation();
  }

  // Initialize pointers
  _tmpTargetImage    = NULL;
  _tmpSourceImage    = NULL;
  _affdLookupTable1  = NULL;
  _affdLookupTable2  = NULL;
  _mffdLookupTable1  = NULL;
  _mffdLookupTable2  = NULL;
  _localLookupTable = new float [FFDLOOKUPTABLESIZE];

}

void irtkSymmetricImageFreeFormRegistration::Initialize(int level)
{
  int i, j, k, n1, n2;
  double u, x, y, z;
  float *ptr;

  // Print debugging information
  this->Debug("irtkSymmetricImageFreeFormRegistration::Initialize(int)");

  // Initialize base class
  this->irtkSymmetricImageRegistration::Initialize(level);

  // Tell optimizer which transformation to optimize
  _optimizer->SetTransformation(_affd1, _affd2);

  // Allocate memory for metric
  _tmp1MetricA = irtkSimilarityMetric::New(_metric1);
  _tmp1MetricB = irtkSimilarityMetric::New(_metric1);
  _tmp2MetricA = irtkSimilarityMetric::New(_metric2);
  _tmp2MetricB = irtkSimilarityMetric::New(_metric2);

  _tmpTargetImage = new irtkGreyImage(_target->GetX(),
                                      _target->GetY(),
                                      _target->GetZ(),
                                      _target->GetT());
  _tmpSourceImage = new irtkGreyImage(_source->GetX(),
                                      _source->GetY(),
                                      _source->GetZ(),
                                      _source->GetT());

  n1 = _target->GetNumberOfVoxels() * 3 / _target->GetT();
  n2 = _source->GetNumberOfVoxels() * 3 / _source->GetT();

  // Allocate memory for lookup table for single-level FFD
  _affdLookupTable1  = new float[n1];
  _affdLookupTable2  = new float[n2];

  // Allocate memory for lookup table for multi-level FFD
  _mffdLookupTable1  = new float[n1];
  _mffdLookupTable2  = new float[n2];

  // Allocate memory for local lookuptable
  _localLookupTable = new float [FFDLOOKUPTABLESIZE];

  for (i = 0; i < FFDLOOKUPTABLESIZE; i++) {
    u = i / (FFDLOOKUPTABLESIZE / 4.0);
    j = (int)floor(u);
    u = u - j;
    _localLookupTable[i] = _affd1->B(j, 1-u);
  }

  // Initialize lookup table for multi-level FFD (this is done only once)
  ptr = _mffdLookupTable1;
  for (k = 0; k < _target->GetZ(); k++) {
    for (j = 0; j < _target->GetY(); j++) {
      for (i = 0; i < _target->GetX(); i++) {
        x = i;
        y = j;
        z = k;
        _target->ImageToWorld(x, y, z);
        _mffd1->Transform(x, y, z);
        ptr[0] = x;
        ptr[1] = y;
        ptr[2] = z;
        ptr += 3;
      }
    }
  }

  // Initialize lookup table for multi-level FFD (this is done only once)
  ptr = _mffdLookupTable2;
  for (k = 0; k < _source->GetZ(); k++) {
    for (j = 0; j < _source->GetY(); j++) {
      for (i = 0; i < _source->GetX(); i++) {
        x = i;
        y = j;
        z = k;
        _source->ImageToWorld(x, y, z);
        _mffd2->Transform(x, y, z);
        ptr[0] = x;
        ptr[1] = y;
        ptr[2] = z;
        ptr += 3;
      }
    }
  }

  // Padding of FFD
  irtkPadding(*tmp_target, this->_Padding, _affd1);
  irtkPadding(*tmp_source, this->_Padding, _affd2);

}

void irtkSymmetricImageFreeFormRegistration::Finalize()
{
  // Print debugging information
  this->Debug("irtkSymmetricImageFreeFormRegistration::Finalize");

  // Push local transformation back on transformation stack
  _mffd1->PushLocalTransformation(_affd1);
  _mffd2->PushLocalTransformation(_affd2);

  // Finalize base class
  this->irtkSymmetricImageRegistration::Finalize();

  delete []_localLookupTable;
}

void irtkSymmetricImageFreeFormRegistration::Finalize(int level)
{
  // Print debugging information
  this->Debug("irtkSymmetricImageFreeFormRegistration::Finalize(int)");

  // Finalize base class
  this->irtkSymmetricImageRegistration::Finalize(level);

  // Check if we are not at the lowest level of resolution
  if (level != 0) {
    if (this->_Subdivision == True) {
      _affd1->Subdivide();
      _affd2->Subdivide();
    } else {
      // Push local transformation back on transformation stack
      _mffd1->PushLocalTransformation(_affd1);
      _mffd2->PushLocalTransformation(_affd2);

      // Create new FFD
      _affd1 = new irtkBSplineFreeFormTransformation(*_target,
          this->_DX / pow(2.0, this->_NumberOfLevels-level),
          this->_DY / pow(2.0, this->_NumberOfLevels-level),
          this->_DZ / pow(2.0, this->_NumberOfLevels-level));
      // Create new FFD
      _affd2 = new irtkBSplineFreeFormTransformation(*_source,
          this->_DX / pow(2.0, this->_NumberOfLevels-level),
          this->_DY / pow(2.0, this->_NumberOfLevels-level),
          this->_DZ / pow(2.0, this->_NumberOfLevels-level));
    }
  }

  delete _tmpTargetImage;
  delete _tmpSourceImage;
  delete _tmp1MetricA;
  delete _tmp2MetricA;
  delete _tmp1MetricB;
  delete _tmp2MetricB;
  delete []_affdLookupTable1;
  delete []_affdLookupTable2;
  delete []_mffdLookupTable1;
  delete []_mffdLookupTable2;
}

void irtkSymmetricImageFreeFormRegistration::UpdateLUT()
{
  int i, j, k;
  double x, y, z;
  float *ptr2mffd;
  float *ptr2affd;

  // Print debugging information
  this->Debug("irtkSymmetricImageFreeFormRegistration::UpdateLUT");

  ptr2affd = _affdLookupTable1;
  ptr2mffd = _mffdLookupTable1;
  for (k = 0; k < _target->GetZ(); k++) {
    for (j = 0; j < _target->GetY(); j++) {
      for (i = 0; i < _target->GetX(); i++) {
        x = i;
        y = j;
        z = k;
        _target->ImageToWorld(x, y, z);
        _affd1->LocalDisplacement(x, y, z);
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

  ptr2affd = _affdLookupTable2;
  ptr2mffd = _mffdLookupTable2;
  for (k = 0; k < _source->GetZ(); k++) {
    for (j = 0; j < _source->GetY(); j++) {
      for (i = 0; i < _source->GetX(); i++) {
        x = i;
        y = j;
        z = k;
        _source->ImageToWorld(x, y, z);
        _affd2->LocalDisplacement(x, y, z);
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

double irtkSymmetricImageFreeFormRegistration::SmoothnessPenalty()
{
  int i, j, k;
  double x, y, z, penalty;

  penalty = 0;
  for (k = 0; k < _affd1->GetZ(); k++) {
    for (j = 0; j < _affd1->GetY(); j++) {
      for (i = 0; i < _affd1->GetX(); i++) {
        x = i;
        y = j;
        z = k;
        _affd1->LatticeToWorld(x, y, z);
        penalty += _affd1->Bending(x, y, z);
      }
    }
  }
  return -penalty / _affd1->NumberOfDOFs();
}

double irtkSymmetricImageFreeFormRegistration::SmoothnessPenalty(int index)
{
  int i, j, k;
  double x, y, z;

  _affd1->IndexToLattice(index, i, j, k);
  x = i;
  y = j;
  z = k;
  _affd1->LatticeToWorld(x, y, z);
  return -_affd1->Bending(x, y, z);
}

double irtkSymmetricImageFreeFormRegistration::VolumePreservationPenalty()
{
  int i, j, k;
  double x, y, z, penalty;

  penalty = 0;
  for (k = 0; k < _affd1->GetZ(); k++) {
    for (j = 0; j < _affd1->GetY(); j++) {
      for (i = 0; i < _affd1->GetZ(); i++) {
        x = i;
        y = j;
        z = k;
        _affd1->LatticeToWorld(x, y, z);
        // Torsten Rohlfing et al. MICCAI'01 (w/o scaling correction):
        penalty += fabs(log(_affd1->irtkTransformation::Jacobian(x, y, z)));
      }
    }
  }

  // Normalize sum by number of DOFs
  return penalty / (double) _affd1->NumberOfDOFs();
}

double irtkSymmetricImageFreeFormRegistration::VolumePreservationPenalty(int index)
{
  int i, j, k;
  double x, y, z;

  _affd1->IndexToLattice(index, i, j, k);
  x = i;
  y = j;
  z = k;
  _affd1->LatticeToWorld(x, y, z);
  return fabs(log(_affd1->irtkTransformation::Jacobian(x, y, z)));
}

double irtkSymmetricImageFreeFormRegistration::TopologyPreservationPenalty()
{
  int i, j, k;
  double x, y, z, jac, penalty;

  penalty = 0;
  for (k = 0; k < _affd1->GetZ()-1; k++) {
    for (j = 0; j < _affd1->GetY()-1; j++) {
      for (i = 0; i < _affd1->GetZ()-1; i++) {
        x = i+0.5;
        y = j+0.5;
        z = k+0.5;
        _affd1->LatticeToWorld(x, y, z);
        jac = _affd1->irtkTransformation::Jacobian(x, y, z);
        if (jac < 0.3) {
          penalty += 10*jac*jac + 0.1/(jac*jac) - 2.0;
        }
      }
    }
  }
  return -penalty;
}

double irtkSymmetricImageFreeFormRegistration::TopologyPreservationPenalty(int index)
{
  int i, j, k, l, m, n;
  double x, y, z, jac, penalty;

  penalty = 0;
  for (l = 0; l <= 1; l++) {
    for (m = 0; m <= 1; m++) {
      for (n = 0; n <= 1; n++) {
        _affd1->IndexToLattice(index, i, j, k);
        x = i+l-0.5;
        y = j+m-0.5;
        z = k+n-0.5;
        _affd1->LatticeToWorld(x, y, z);
        jac = _affd1->irtkTransformation::Jacobian(x, y, z);
        if (jac < 0.3) {
          penalty += 10*jac*jac + 0.1/(jac*jac) - 2.0;
        }
      }
    }
  }
  return -penalty;
}

double irtkSymmetricImageFreeFormRegistration::InverseConsistencyPenalty()
{
  int i, j, k;
  double x1, y1, z1, x2, y2, z2, penalty;

  penalty = 0;
  for (k = 0; k < _affd1->GetZ(); k++) {
    for (j = 0; j < _affd1->GetY(); j++) {
      for (i = 0; i < _affd1->GetZ(); i++) {
        x1 = i;
        y1 = j;
        z1 = k;
        _affd1->LatticeToWorld(x1, y1, z1);
        x2 = x1;
        y2 = y1;
        z2 = z1;
        _affd1->Transform(x2, y2, z2);
        _affd2->Transform(x2, y2, z2);
        penalty += (x1-x2)*(x1-x2)+(y1-y2)*(y1-y2)+(z1-z2)*(z1-z2);
      }
    }
  }
  for (k = 0; k < _affd2->GetZ(); k++) {
    for (j = 0; j < _affd2->GetY(); j++) {
      for (i = 0; i < _affd2->GetZ(); i++) {
        x1 = i;
        y1 = j;
        z1 = k;
        _affd2->LatticeToWorld(x1, y1, z1);
        x2 = x1;
        y2 = y1;
        z2 = z1;
        _affd2->Transform(x2, y2, z2);
        _affd1->Transform(x2, y2, z2);
        penalty += (x1-x2)*(x1-x2)+(y1-y2)*(y1-y2)+(z1-z2)*(z1-z2);
      }
    }
  }

  // Normalize sum by number of DOFs
  return -penalty / (double) _affd1->NumberOfDOFs();
}

double irtkSymmetricImageFreeFormRegistration::InverseConsistencyPenalty(int index)
{
  int i, j, k;
  double ice, x1, y1, z1, x2, y2, z2;

  ice = 0;
  _affd1->IndexToLattice(index, i, j, k);
  x1 = i;
  y1 = j;
  z1 = k;
  _affd1->LatticeToWorld(x1, y1, z1);
  x2 = x1;
  y2 = y1;
  z2 = z1;
  _affd1->Transform(x2, y2, z2);
  _affd2->Transform(x2, y2, z2);
  ice += (x1-x2)*(x1-x2)+(y1-y2)*(y1-y2)+(z1-z2)*(z1-z2);
  _affd2->IndexToLattice(index, i, j, k);
  x1 = i;
  y1 = j;
  z1 = k;
  _affd2->LatticeToWorld(x1, y1, z1);
  x2 = x1;
  y2 = y1;
  z2 = z1;
  _affd2->Transform(x2, y2, z2);
  _affd1->Transform(x2, y2, z2);
  ice += (x1-x2)*(x1-x2)+(y1-y2)*(y1-y2)+(z1-z2)*(z1-z2);
  return -ice;
}

double irtkSymmetricImageFreeFormRegistration::Evaluate()
{
  // Image coordinates
  int i, j, k, t;
  // World coordinates
  double x, y, z;
  // Pointer to reference data
  irtkGreyPixel *ptr2target;
  irtkGreyPixel *ptr2source;
  irtkGreyPixel *ptr2tmp;
  float *ptr;

  // Print debugging information
  this->Debug("irtkSymmetricImageFreeFormRegistration::Evaluate");

  // Initialize metric
  _metric1->Reset();
  _metric2->Reset();

  // Loop over all voxels in the target (reference) volume
  ptr2target = _target->GetPointerToVoxels();
  ptr2tmp    = _tmpSourceImage->GetPointerToVoxels();
  for (t = 0; t < _target->GetT(); t++) {
    ptr        = _mffdLookupTable1;
    for (k = 0; k < _target->GetZ(); k++) {
      for (j = 0; j < _target->GetY(); j++) {
        for (i = 0; i < _target->GetX(); i++) {
          // Check whether reference point is valid
          if (*ptr2target >= 0) {
            x = i;
            y = j;
            z = k;
            _target->ImageToWorld(x, y, z);
            _affd1->LocalDisplacement(x, y, z);
            x += ptr[0];
            y += ptr[1];
            z += ptr[2];
            _source->WorldToImage(x, y, z);
            // Check whether transformed point is inside volume
            if ((x > _source_x1) && (x < _source_x2) &&
                (y > _source_y1) && (y < _source_y2) &&
                (z > _source_z1) && (z < _source_z2)) {
              // Add sample to metric
              *ptr2tmp =  round(_interpolator1->EvaluateInside(x, y, z, t));
              _metric1->Add(*ptr2target, *ptr2tmp);
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

  // Loop over all voxels in the source (reference) volume
  ptr2source = _source->GetPointerToVoxels();
  ptr2tmp    = _tmpTargetImage->GetPointerToVoxels();
  for (t = 0; t < _source->GetT(); t++) {
    ptr        = _mffdLookupTable2;
    for (k = 0; k < _source->GetZ(); k++) {
      for (j = 0; j < _source->GetY(); j++) {
        for (i = 0; i < _source->GetX(); i++) {
          // Check whether reference point is valid
          if (*ptr2source >= 0) {
            x = i;
            y = j;
            z = k;
            _source->ImageToWorld(x, y, z);
            _affd2->LocalDisplacement(x, y, z);
            x += ptr[0];
            y += ptr[1];
            z += ptr[2];
            _target->WorldToImage(x, y, z);
            // Check whether transformed point is inside volume
            if ((x > _target_x1) && (x < _target_x2) &&
                (y > _target_y1) && (y < _target_y2) &&
                (z > _target_z1) && (z < _target_z2)) {
              // Add sample to metric
              *ptr2tmp =  round(_interpolator2->EvaluateInside(x, y, z, t));
              _metric2->Add(*ptr2source, *ptr2tmp);
            } else {
              *ptr2tmp = -1;
            }
          }
          // Increment pointers to next voxel
          ptr2tmp++;
          ptr2source++;
          ptr += 3;
        }
      }
    }
  }

  // Evaluate similarity measure
  double similarity = (_metric1->Evaluate() + _metric2->Evaluate()) / 2.0;

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
  // Add inverse consistency
  if (this->_Lambda4 > 0) {
    similarity += this->_Lambda4*this->InverseConsistencyPenalty();
  }

  // Return similarity measure + penalty terms
  return similarity;
}

double irtkSymmetricImageFreeFormRegistration::EvaluateDerivative1(int index, double step)
{
  float *ptr;
  irtkPoint p1, p2;
  double bi, bj, bk, dx, dy, dz, p[3];
  int i, j, k, i1, i2, j1, j2, k1, k2, dim, t;
  irtkGreyPixel *ptr2target, *ptr2tmp;
  irtkSimilarityMetric *tmp1MetricA, *tmp1MetricB;

  // Print debugging information
  this->Debug("irtkSymmetricImageFreeFormRegistration::EvaluateDerivative1(int, double)");

  tmp1MetricA = _tmp1MetricA;
  tmp1MetricB = _tmp1MetricB;

  // Initialize metrics for forward and backward derivative steps
  tmp1MetricA->Reset(_metric1);
  tmp1MetricB->Reset(_metric1);

  // Calculate bounding box of control point in world coordinates
  _affd1->BoundingBox(index, p1, p2);
  _target->WorldToImage(p1);
  _target->WorldToImage(p2);

  // Calculate bounding box of control point in image coordinates
  _affd1->BoundingBox(_target, index, i1, j1, k1, i2, j2, k2, 1.0 / _SpeedupFactor);

  // Calculate incremental changes in lattice coordinates when looping
  // over target
  dx = (FFDLOOKUPTABLESIZE-1)/(p2._x-p1._x);
  dy = (FFDLOOKUPTABLESIZE-1)/(p2._y-p1._y);
  dz = (FFDLOOKUPTABLESIZE-1)/(p2._z-p1._z);

  // Calculate whether this DOF corresponds to x, y or z-displacement
  dim = int(index / (_affd1->GetX()*_affd1->GetY()*_affd1->GetZ()));

  // Loop over all voxels in the target (reference) volume
  for (t = 0; t < _target->GetT(); t++) {
    for (k = k1; k <= k2; k++) {
      bk = step * _localLookupTable[round(dz*(k-p1._z))];
      for (j = j1; j <= j2; j++) {
        ptr2target = _target->GetPointerToVoxels(i1, j, k, t);
        ptr        = &(_affdLookupTable1[3*_target->VoxelToIndex(i1, j, k)]);
        bj = bk * _localLookupTable[round(dy*(j-p1._y))];
        ptr2tmp  = _tmpSourceImage->GetPointerToVoxels(i1, j, k, t);
        for (i = i1; i <= i2; i++) {

          // Check whether reference point is valid
          if (*ptr2target >= 0) {
            bi = bj * _localLookupTable[round(dx*(i-p1._x))];

            // Delete old samples from both metrics
            if (*ptr2tmp != -1) {
              tmp1MetricA->Delete(*ptr2target, *ptr2tmp);
              tmp1MetricB->Delete(*ptr2target, *ptr2tmp);
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
              tmp1MetricA->Add(*ptr2target, round(_interpolator1->EvaluateInside(p[0], p[1], p[2], t)));
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
              tmp1MetricB->Add(*ptr2target, round(_interpolator1->EvaluateInside(p[0], p[1], p[2], t)));
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
  double dof = _affd1->Get(index);

  // Evaluate similarity measure
  double similarityA = tmp1MetricA->Evaluate();

  // Add penalties
  _affd1->Put(index, dof + step);

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
  // Inverse consistency
  if (this->_Lambda4 > 0) {
    similarityA += this->_Lambda4*this->InverseConsistencyPenalty(index);
  }

  // Evaluate similarity measure
  double similarityB = tmp1MetricB->Evaluate();

  // Add penalties
  _affd1->Put(index, dof - step);

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
  // Inverse consistency
  if (this->_Lambda4 > 0) {
    similarityB += this->_Lambda4*this->InverseConsistencyPenalty(index);
  }

  // Restore value of DOF for which we calculate the derivative
  _affd1->Put(index, dof);

  return similarityA - similarityB;
}

double irtkSymmetricImageFreeFormRegistration::EvaluateDerivative2(int index, double step)
{
  float *ptr;
  irtkPoint p1, p2;
  double bi, bj, bk, dx, dy, dz, p[3];
  int i, j, k, i1, i2, j1, j2, k1, k2, dim, t;
  irtkGreyPixel *ptr2source, *ptr2tmp;
  irtkSimilarityMetric *tmp2MetricA, *tmp2MetricB;

  // Print debugging information
  this->Debug("irtkSymmetricImageFreeFormRegistration::EvaluateDerivative2(int, double)");

  tmp2MetricA = _tmp2MetricA;
  tmp2MetricB = _tmp2MetricB;

  // Initialize metrics for forward and backward derivative steps
  tmp2MetricA->Reset(_metric2);
  tmp2MetricB->Reset(_metric2);

  // Calculate bounding box of control point in world coordinates
  _affd2->BoundingBox(index, p1, p2);
  _source->WorldToImage(p1);
  _source->WorldToImage(p2);

  // Calculate bounding box of control point in image coordinates
  _affd2->BoundingBox(_source, index, i1, j1, k1, i2, j2, k2, 1.0 / _SpeedupFactor);

  // Calculate incremental changes in lattice coordinates when looping
  // over target
  dx = (FFDLOOKUPTABLESIZE-1)/(p2._x-p1._x);
  dy = (FFDLOOKUPTABLESIZE-1)/(p2._y-p1._y);
  dz = (FFDLOOKUPTABLESIZE-1)/(p2._z-p1._z);

  // Calculate whether this DOF corresponds to x, y or z-displacement
  dim = int(index / (_affd2->GetX()*_affd2->GetY()*_affd2->GetZ()));

  // Loop over all voxels in the target (reference) volume
  for (t = 0; t < _source->GetT(); t++) {
    for (k = k1; k <= k2; k++) {
      bk = step * _localLookupTable[round(dz*(k-p1._z))];
      for (j = j1; j <= j2; j++) {
        ptr2source = _source->GetPointerToVoxels(i1, j, k, t);
        ptr        = &(_affdLookupTable2[3*_source->VoxelToIndex(i1, j, k)]);
        bj = bk * _localLookupTable[round(dy*(j-p1._y))];
        ptr2tmp  = _tmpTargetImage->GetPointerToVoxels(i1, j, k, t);
        for (i = i1; i <= i2; i++) {

          // Check whether reference point is valid
          if (*ptr2source >= 0) {
            bi = bj * _localLookupTable[round(dx*(i-p1._x))];

            // Delete old samples from both metrics
            if (*ptr2tmp != -1) {
              tmp2MetricA->Delete(*ptr2source, *ptr2tmp);
              tmp2MetricB->Delete(*ptr2source, *ptr2tmp);
            }

            p[0] = ptr[0];
            p[1] = ptr[1];
            p[2] = ptr[2];
            p[dim] += bi;

            // Convert transformed point to image coordinates
            _target->WorldToImage(p[0], p[1], p[2]);

            // Check whether transformed point is inside volume
            if ((p[0] > _target_x1) && (p[0] < _target_x2) &&
                (p[1] > _target_y1) && (p[1] < _target_y2) &&
                (p[2] > _target_z1) && (p[2] < _target_z2)) {

              // Add sample to metric
              tmp2MetricA->Add(*ptr2source, round(_interpolator2->EvaluateInside(p[0], p[1], p[2], t)));
            }

            p[0] = ptr[0];
            p[1] = ptr[1];
            p[2] = ptr[2];
            p[dim] -= bi;

            // Convert transformed point to image coordinates
            _target->WorldToImage(p[0], p[1], p[2]);

            // Check whether transformed point is inside volume
            if ((p[0] > _target_x1) && (p[0] < _target_x2) &&
                (p[1] > _target_y1) && (p[1] < _target_y2) &&
                (p[2] > _target_z1) && (p[2] < _target_z2)) {

              // Add sample to metric
              tmp2MetricB->Add(*ptr2source, round(_interpolator2->EvaluateInside(p[0], p[1], p[2], t)));
            }
          }

          // Increment pointers to next voxel
          ptr2source++;
          ptr2tmp++;
          ptr += 3;
        }
      }
    }
  }

  // Save value of DOF for which we calculate the derivative
  double dof = _affd2->Get(index);

  // Evaluate similarity measure
  double similarityA = tmp2MetricA->Evaluate();

  // Add penalties
  _affd2->Put(index, dof + step);

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
  // Inverse consistency
  if (this->_Lambda4 > 0) {
    similarityA += this->_Lambda4*this->InverseConsistencyPenalty(index);
  }

  // Evaluate similarity measure
  double similarityB = tmp2MetricB->Evaluate();

  // Add penalties
  _affd2->Put(index, dof - step);

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
  // Inverse consistency
  if (this->_Lambda4 > 0) {
    similarityB += this->_Lambda4*this->InverseConsistencyPenalty(index);
  }

  // Restore value of DOF for which we calculate the derivative
  _affd2->Put(index, dof);

  return similarityA - similarityB;
}


double irtkSymmetricImageFreeFormRegistration::EvaluateGradient(float step, float *dx)
{
  int i;
  double norm;

  // Update lookup table
  this->UpdateLUT();

  for (i = 0; i < _affd1->NumberOfDOFs(); i++) {
    if (_affd1->irtkTransformation::GetStatus(i) == _Active) {
      dx[i] = this->EvaluateDerivative1(i, step);
    } else {
      dx[i] = 0;
    }
  }
  for (i = 0; i < _affd2->NumberOfDOFs(); i++) {
    if (_affd2->irtkTransformation::GetStatus(i) == _Active) {
      dx[i+_affd1->NumberOfDOFs()] = this->EvaluateDerivative2(i, step);
    } else {
      dx[i+_affd1->NumberOfDOFs()] = 0;
    }
  }

  // Calculate norm of vector
  norm = 0;
  for (i = 0; i < _affd1->NumberOfDOFs() + _affd2->NumberOfDOFs(); i++) {
    norm += dx[i] * dx[i];
  }

  // Normalize vector
  norm = sqrt(norm);
  if (norm > 0) {
    for (i = 0; i < _affd1->NumberOfDOFs() + _affd2->NumberOfDOFs(); i++) {
      dx[i] /= norm;
    }
  } else {
    for (i = 0; i < _affd1->NumberOfDOFs() + _affd2->NumberOfDOFs(); i++) {
      dx[i] = 0;
    }
  }

  return norm;
}

Bool irtkSymmetricImageFreeFormRegistration::Read(char *buffer1, char *buffer2, int &level)
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
  if (strstr(buffer1, "Lambda4") != NULL) {
    this->_Lambda4 = atof(buffer2);
    cout << "Lambda 4 is ... " << this->_Lambda3 << endl;
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
    return this->irtkSymmetricImageRegistration::Read(buffer1, buffer2, level);
  } else {
    return ok;
  }
}

void irtkSymmetricImageFreeFormRegistration::Write(ostream &to)
{
  to << "\n#\n# Non-rigid registration parameters\n#\n\n";
  to << "Lambda1                           = " << this->_Lambda1 << endl;
  to << "Lambda2                           = " << this->_Lambda2 << endl;
  to << "Lambda3                           = " << this->_Lambda3 << endl;
  to << "Lambda4                           = " << this->_Lambda4 << endl;
  to << "Control point spacing in X        = " << this->_DX << endl;
  to << "Control point spacing in Y        = " << this->_DY << endl;
  to << "Control point spacing in Z        = " << this->_DZ << endl;
  if (_Subdivision == True) {
    to << "Subdivision                       = True" << endl;
  } else {
    to << "Subdivision                       = False" << endl;
  }
  to << "Speedup factor                    = " << this->_SpeedupFactor << endl;

  this->irtkSymmetricImageRegistration::Write(to);
}
