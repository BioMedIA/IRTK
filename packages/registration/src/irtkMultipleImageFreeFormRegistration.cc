/*=========================================================================

Library   : Image Registration Toolkit (IRTK)
Module    : $Id$
Copyright : Imperial College, Department of Computing
Visual Information Processing (VIP), 2008 onwards
Date      : $Date$
Version   : $Revision$
Changes   : $Author$

=========================================================================*/

#ifdef HAS_VTK

#include <irtkRegistration.h>

#undef HAS_TBB

// Used as temporary memory for transformed intensities
irtkGreyImage **_mtmpImage;

vtkPolyData *_tmpptarget;

// The original target and source images
extern irtkGreyImage **tmp_mtarget, **tmp_msource;

#include <irtkMultiThreadedImageFreeFormRegistration.h>

irtkMultipleImageFreeFormRegistration::irtkMultipleImageFreeFormRegistration()
{
  // Print debugging information
  this->Debug("irtkMultipleImageFreeFormRegistration::irtkMultipleImageFreeFormRegistration");

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
  _MFFDMode    = true;

  // weight
  _weight = NULL;

  _level = 0;
}

void irtkMultipleImageFreeFormRegistration::GuessParameter()
{
  int i;
  double xsize, ysize, zsize, spacing;

  if ((_target == NULL) || (_source == NULL)) {
    cerr << "irtkMultipleImageFreeFormRegistration::GuessParameter: Target and source image not found" << endl;
    exit(1);
  }

  // Default parameters for registration
  _NumberOfLevels     = 3;
  _NumberOfBins       = 64;

  // Default parameters for optimization
  _SimilarityMeasure  = NMI;
  _OptimizationMethod = GradientDescent;
  _Epsilon            = 0.0001;

#ifdef HAS_VTK

  if (_ptarget != NULL && _psource != NULL) {
    _Lregu          = 0.001;
  } else {
    _Lregu          = 0;
  }

#endif

  // Read target pixel size
  _target[0]->GetPixelSize(&xsize, &ysize, &zsize);

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
  _source[0]->GetPixelSize(&xsize, &ysize, &zsize);

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
  _DX                 =_target[0]->GetX() * spacing / 10.0;
  _DY                 =_target[0]->GetX() * spacing / 10.0;
  _DZ                 =_target[0]->GetX() * spacing / 10.0;
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
  if ((_target[0]->Get(_target[0]->GetX()-1, 0, 0)                                 == _target[0]->Get(0, 0, 0)) &&
      (_target[0]->Get(0, _target[0]->GetY()-1, 0)                                 == _target[0]->Get(0, 0, 0)) &&
      (_target[0]->Get(0, 0, _target[0]->GetZ()-1)                                 == _target[0]->Get(0, 0, 0)) &&
      (_target[0]->Get(_target[0]->GetX()-1, _target[0]->GetY()-1, 0)                 == _target[0]->Get(0, 0, 0)) &&
      (_target[0]->Get(0, _target[0]->GetY()-1, _target[0]->GetZ()-1)                 == _target[0]->Get(0, 0, 0)) &&
      (_target[0]->Get(_target[0]->GetX()-1, 0, _target[0]->GetZ()-1)                 == _target[0]->Get(0, 0, 0)) &&
      (_target[0]->Get(_target[0]->GetX()-1, _target[0]->GetY()-1, _target[0]->GetZ()-1) == _target[0]->Get(0, 0, 0))) {
    _TargetPadding = _target[0]->Get(0, 0, 0);
  }
}

void irtkMultipleImageFreeFormRegistration::Initialize()
{
  int i, j;
  double u;

  // Print debugging information
  this->Debug("irtkMultipleImageFreeFormRegistration::Initialize");

  // Initialize base class
  this->irtkMultipleImageRegistration::Initialize();

  // Pointer to multi-level FFD
  _mffd = (irtkMultiLevelFreeFormTransformation *)_transformation;

  // Create FFD
  if (_MFFDMode == false) {
    if (_mffd->NumberOfLevels() == 0) {
      _affd = new irtkBSplineFreeFormTransformation(*_target[0], this->_DX, this->_DY, this->_DZ);
    } else {
      _affd = (irtkBSplineFreeFormTransformation *)_mffd->PopLocalTransformation();
    }
  } else {
    _affd = new irtkBSplineFreeFormTransformation(*_target[0], this->_DX, this->_DY, this->_DZ);
  }

  // Initialize pointers
  _mtmpImage         = NULL;
  _tmpptarget        = NULL;
  _affdLookupTable  = NULL;
  _mffdLookupTable  = NULL;
  _localLookupTable = new float [FFDLOOKUPTABLESIZE];

  for (i = 0; i < FFDLOOKUPTABLESIZE; i++) {
    u = i / (FFDLOOKUPTABLESIZE / 4.0);
    j = (int)floor(u);
    u = u - j;
    _localLookupTable[i] = _affd->B(j, 1-u);
  }

}

void irtkMultipleImageFreeFormRegistration::Initialize(int level)
{
  int i, j, k, l, n;
  double x, y, z;
  float *ptr;

  // Print debugging information
  this->Debug("irtkMultipleImageFreeFormRegistration::Initialize(int)");

  // Initialize base class
  this->irtkMultipleImageRegistration::Initialize(level);

  // Set level
  _level = level;

  // Tell optimizer which transformation to optimize
  _optimizer->SetTransformation(_affd);

  // Allocate memory for metric
  _tmpMetricA = new irtkSimilarityMetric*[_numberOfImages];
  _tmpMetricB = new irtkSimilarityMetric*[_numberOfImages];

  // Allocate memory for temp image
  _mtmpImage = new irtkGreyImage*[_numberOfImages];

  // Allocate memory for lookup tables
  _affdLookupTable = new float*[_numberOfImages];
  _mffdLookupTable = new float*[_numberOfImages];

  for (l = 0; l < _numberOfImages; l++) {
    if(_weight != NULL) {
      cout << "image " << l << "th's external weight is " << _weight[_level][l] << endl;
    }
    _tmpMetricA[l] = irtkSimilarityMetric::New(_metric[l]);
    _tmpMetricB[l] = irtkSimilarityMetric::New(_metric[l]);
    _mtmpImage[l] = new irtkGreyImage(_target[l]->GetX(),
                                      _target[l]->GetY(),
                                      _target[l]->GetZ(),
                                      _target[l]->GetT());

    n = _target[l]->GetNumberOfVoxels() * 3 / _target[l]->GetT();

    // Allocate memory for lookup table for single-level FFD
    _affdLookupTable[l]  = new float[n];

    // Allocate memory for lookup table for multi-level FFD
    _mffdLookupTable[l]  = new float[n];

    // Initialize lookup table for multi-level FFD (this is done only once)
    ptr = _mffdLookupTable[l];
    for (k = 0; k < _target[l]->GetZ(); k++) {
      for (j = 0; j < _target[l]->GetY(); j++) {
        for (i = 0; i < _target[l]->GetX(); i++) {
          x = i;
          y = j;
          z = k;
          _target[l]->ImageToWorld(x, y, z);
          _mffd->Transform(x, y, z);
          ptr[0] = x;
          ptr[1] = y;
          ptr[2] = z;
          ptr += 3;
        }
      }
    }
  }
  #ifdef HAS_VTK
  if(_ptarget != NULL){
      _tmpptarget = vtkPolyData::New();
      _tmpptarget->DeepCopy(_ptarget);
      double p[3];
      for (i = 0; i < _ptarget->GetNumberOfPoints(); i++) {
          _ptarget->GetPoints()->GetPoint(i,p);
          _mffd->Transform(p[0],p[1],p[2]);
          _tmpptarget->GetPoints()->SetPoint(i,p);
      }
  }
  #endif
}

void irtkMultipleImageFreeFormRegistration::Finalize()
{
  // Print debugging information
  this->Debug("irtkMultipleImageFreeFormRegistration::Finalize");

  // Push local transformation back on transformation stack
  _mffd->PushLocalTransformation(_affd);

  // Finalize base class
  this->irtkMultipleImageRegistration::Finalize();

  delete []_localLookupTable;
}

void irtkMultipleImageFreeFormRegistration::Finalize(int level)
{
  // Print debugging information
  this->Debug("irtkMultipleImageFreeFormRegistration::Finalize(int)");

  // Finalize base class
  this->irtkMultipleImageRegistration::Finalize(level);

  // Check if we are not at the lowest level of resolution
  if (level != 0) {
    if (this->_Subdivision == true) {
      _affd->Subdivide();
    } else {
      // Push local transformation back on transformation stack
      _mffd->PushLocalTransformation(_affd);

      // Create new FFD
      _affd = new irtkBSplineFreeFormTransformation(*_target[0],
              this->_DX / pow(2.0, this->_NumberOfLevels-level),
              this->_DY / pow(2.0, this->_NumberOfLevels-level),
              this->_DZ / pow(2.0, this->_NumberOfLevels-level));
    }
  }
#ifdef HAS_VTK
  if(_ptarget != NULL){
      _tmpptarget->Delete();
  }
#endif
  for (int n = 0; n < _numberOfImages; n++) {
    delete _mtmpImage[n];
    delete _affdLookupTable[n];
    delete _mffdLookupTable[n];
    delete _tmpMetricA[n];
    delete _tmpMetricB[n];
  }
  delete []_mtmpImage;
  delete []_tmpMetricA;
  delete []_tmpMetricB;
  delete []_affdLookupTable;
  delete []_mffdLookupTable;
}

void irtkMultipleImageFreeFormRegistration::UpdateLUT()
{
  int i, j, k, n;
  double x, y, z;
  float *ptr2mffd;
  float *ptr2affd;

  // Print debugging information
  this->Debug("irtkMultipleImageFreeFormRegistration::UpdateLUT");

  for (n = 0; n < _numberOfImages; n++) {
    ptr2affd = _affdLookupTable[n];
    ptr2mffd = _mffdLookupTable[n];
    for (k = 0; k < _target[n]->GetZ(); k++) {
      for (j = 0; j < _target[n]->GetY(); j++) {
        for (i = 0; i < _target[n]->GetX(); i++) {
          x = i;
          y = j;
          z = k;
          _target[n]->ImageToWorld(x, y, z);
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
}

double irtkMultipleImageFreeFormRegistration::LandMarkPenalty(int index)
{
  int i,k;
#ifdef HAS_VTK
  vtkIdType j;
  double dx = 0, dy = 0, dz = 0, min, max, d = 0, distance = 0 , p[3], q[3];
  irtkPoint p1, p2, pt;

  if(index != -1) {
    _affd->BoundingBox(index, p1, p2);
    _target[0]->WorldToImage(p1);
    _target[0]->WorldToImage(p2);
    dx = (FFDLOOKUPTABLESIZE-1)/(p2._x-p1._x);
    dy = (FFDLOOKUPTABLESIZE-1)/(p2._y-p1._y);
    dz = (FFDLOOKUPTABLESIZE-1)/(p2._z-p1._z);

    min = 0;
    max = FFDLOOKUPTABLESIZE-1;
  }

  if (_ptarget == NULL || _psource == NULL) {
      return 0;
  } else if(_ptarget->GetNumberOfPoints() == 0 || _psource->GetNumberOfPoints() == 0) {
      return 0;
  }

  vtkCellLocator *locator = vtkCellLocator::New(); 
  locator->SetDataSet(_psource); // data represents the surface 
  locator->BuildLocator(); 

  vtkGenericCell *tmp = vtkGenericCell::New();

  int valid = 0;
  if(index != -1) {
      for (i = 0; i < _ptarget->GetNumberOfPoints(); i++) {
          _ptarget->GetPoints()->GetPoint(i,p);
          pt._x = p[0];
          pt._y = p[1];
          pt._z = p[2];
          _target[0]->WorldToImage(pt);
          if(round(dz*(pt._z-p1._z))<=max && round(dz*(pt._z-p1._z))>=min
              &&round(dy*(pt._y-p1._y))<=max && round(dy*(pt._y-p1._y))>=min
              &&round(dx*(pt._x-p1._x))<=max && round(dx*(pt._x-p1._x))>=min) {
                  valid = 1;
                  break;
          }
      }
  }else{
      valid = 1;
  }
  if(valid == 1) {
      for (i = 0; i < _tmpptarget->GetNumberOfPoints(); i++) {
          _tmpptarget->GetPoints()->GetPoint(i,p);
          q[0] = p[0]; q[1] = p[1]; q[2] = p[2];
          _affd->LocalDisplacement(q[0],q[1],q[2]);
          p[0] += q[0];
          p[1] += q[1];
          p[2] += q[2];
          locator->FindClosestPoint(p,q,tmp,j,k,d);
          distance += d;
      }
  }else{
      distance = 0;
  }
  tmp->Delete();
  locator->Delete();
  return -(distance/double(_tmpptarget->GetNumberOfPoints()));
#else
  return 0;
#endif

}

double irtkMultipleImageFreeFormRegistration::SmoothnessPenalty()
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
        //penalty += _mffd->Bending(x,y,z);
        penalty += _affd->Bending(x, y, z);
      }
    }
  }
  return -penalty / _affd->NumberOfDOFs();
}

double irtkMultipleImageFreeFormRegistration::SmoothnessPenalty(int index)
{
  int i, j, k;
  double x, y, z, penalty;

  penalty = 0;

  _affd->IndexToLattice(index, i, j, k);
  x = i;
  y = j;
  z = k;
  _affd->LatticeToWorld(x, y, z);
  //penalty += _mffd->Bending(x,y,z);
  penalty += _affd->Bending(x, y, z);
  return -penalty;
}

double irtkMultipleImageFreeFormRegistration::VolumePreservationPenalty()
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
        //jacobian = _affd->irtkTransformation::Jacobian(x, y, z);
        //if (jacobian < 0.001)
        //jacobian = 0.001;
        // Torsten Rohlfing et al. MICCAI'01 (w/o scaling correction):
        penalty += fabs(log(jacobian));
      }
    }
  }

  // Normalize sum by number of DOFs
  return -penalty / (double) _affd->NumberOfDOFs();
}

double irtkMultipleImageFreeFormRegistration::VolumePreservationPenalty(int index)
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


double irtkMultipleImageFreeFormRegistration::TopologyPreservationPenalty()
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

double irtkMultipleImageFreeFormRegistration::TopologyPreservationPenalty(int index)
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

double irtkMultipleImageFreeFormRegistration::Evaluate()
{
#ifndef HAS_TBB
  // Image coordinates
  int i, j, k, t, n;
  // World coordinates
  double x, y, z;
  double *sweight = new double[_numberOfImages];
  // Pointer to reference data
  irtkGreyPixel *ptr2target;
  irtkGreyPixel *ptr2tmp;
  float *ptr;
#endif

  // Print debugging information
  this->Debug("irtkMultipleImageFreeFormRegistration::Evaluate");

#ifdef HAS_TBB
  irtkMultiThreadedImageFreeFormRegistrationEvaluate evaluate(this);
  parallel_reduce(blocked_range<int>(0, _target->GetZ(), 1), evaluate);
#else

  for (n = 0; n < _numberOfImages; n++) {
    // Initialize metric
    _metric[n]->Reset();
    sweight[n] = 0;
    // Loop over all voxels in the target (reference) volume
    ptr2target = _target[n]->GetPointerToVoxels();
    ptr2tmp    = _mtmpImage[n]->GetPointerToVoxels();
    for (t = 0; t < _target[n]->GetT(); t++) {
      ptr        = _mffdLookupTable[n];
      for (k = 0; k < _target[n]->GetZ(); k++) {
        for (j = 0; j < _target[n]->GetY(); j++) {
          for (i = 0; i < _target[n]->GetX(); i++) {
            // Check whether reference point is valid
            if (*ptr2target >= 0) {
              x = i;
              y = j;
              z = k;
              _target[n]->ImageToWorld(x, y, z);
              _affd->LocalDisplacement(x, y, z);
              x += ptr[0];
              y += ptr[1];
              z += ptr[2];
              _source[n]->WorldToImage(x, y, z);
              // Check whether transformed point is inside volume
              if ((x > _source_x1[n]) && (x < _source_x2[n]) &&
                  (y > _source_y1[n]) && (y < _source_y2[n]) &&
                  ((_target[n]->GetZ() == 1 && round(z) == 0)
                   ||( _target[n]->GetZ() != 1
                       && z > _source_z1[n] && z < _source_z2[n]))) {
                // Add sample to metric
                *ptr2tmp =  round(_interpolator[n]->EvaluateInside(x, y, z, t));
                _metric[n]->Add(*ptr2target, *ptr2tmp);
                if(_weight == NULL)
                  sweight[n]++;
                else
                  sweight[n] += _weight[_level][n];
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
  }
#endif

  // Evaluate similarity measure
  double similarity = combine_mysimilarity(_metric,sweight,_numberOfImages);

  // Add penalty for landmark regulation
  if (this->_Lregu > 0) {
    similarity += this->_Lregu * this->LandMarkPenalty();
  }
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

  delete []sweight;

  // Return similarity measure + penalty terms
  return similarity;
}

double irtkMultipleImageFreeFormRegistration::EvaluateDerivative(int index, double step)
{
  float *ptr;
  irtkPoint p1, p2;
  double bi, bj, bk, dx, dy, dz, p[3],x,y,z;
  int i, j, k, i1, i2, j1, j2, k1, k2, dim, t, n,min,max;
  irtkGreyPixel *ptr2target, *ptr2tmp;
  irtkSimilarityMetric **tmpMetricA, **tmpMetricB;
  double *weight = new double[_numberOfImages];

  // Print debugging information
  this->Debug("irtkMultipleImageFreeFormRegistration::EvaluateDerivative(int, double)");

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

  // Calculate whether this DOF corresponds to x, y or z-displacement
  dim = int(index / (_affd->GetX()*_affd->GetY()*_affd->GetZ()));

  for (n = 0; n < _numberOfImages; n++) {
    // Initialize metrics for forward and backward derivative steps
    tmpMetricA[n]->Reset(_metric[n]);
    tmpMetricB[n]->Reset(_metric[n]);
    weight[n] = 0;

    // Calculate bounding box of control point in world coordinates
    _affd->BoundingBox(index, p1, p2);
    _target[0]->WorldToImage(p1);
    _target[0]->WorldToImage(p2);

    // Calculate bounding box of control point in image coordinates
    _affd->MultiBoundingBox(_target[n], index, i1, j1, k1, i2, j2, k2, 1.0 / _SpeedupFactor);

    // Calculate incremental changes in lattice coordinates when looping
    // over target
    dx = (FFDLOOKUPTABLESIZE-1)/(p2._x-p1._x);
    dy = (FFDLOOKUPTABLESIZE-1)/(p2._y-p1._y);
    dz = (FFDLOOKUPTABLESIZE-1)/(p2._z-p1._z);

    min = round((FFDLOOKUPTABLESIZE-1)*(0.5 - 0.5/_SpeedupFactor));
    max = round((FFDLOOKUPTABLESIZE-1)*(0.5 + 0.5/_SpeedupFactor));

    // Loop over all voxels in the target (reference) volume
    for (t = 0; t < _target[n]->GetT(); t++) {
      for (k = k1; k <= k2; k++) {
        for (j = j1; j <= j2; j++) {
          ptr2target = _target[n]->GetPointerToVoxels(i1, j, k, t);
          ptr        = &(_affdLookupTable[n][3*_target[n]->VoxelToIndex(i1, j, k)]);
          ptr2tmp  = _mtmpImage[n]->GetPointerToVoxels(i1, j, k, t);
          for (i = i1; i <= i2; i++) {
            x = i; y = j; z = k;
            _target[n]->ImageToWorld(x,y,z);
            _target[0]->WorldToImage(x,y,z);
            if(round(dz*(z-p1._z))<=max && round(dz*(z-p1._z))>=min) {
              bk = step * _localLookupTable[round(dz*(z-p1._z))];
              if(round(dy*(y-p1._y))<=max && round(dy*(y-p1._y))>=min) {
                bj = bk * _localLookupTable[round(dy*(y-p1._y))];
                // Check whether reference point is valid
                if (*ptr2target >= 0 && round(dx*(x-p1._x))<=max && round(dx*(x-p1._x))>=min) {
                  bi = bj * _localLookupTable[round(dx*(x-p1._x))];
                  // Delete old samples from both metrics
                  if (*ptr2tmp != -1) {
                    tmpMetricA[n]->Delete(*ptr2target, *ptr2tmp);
                    tmpMetricB[n]->Delete(*ptr2target, *ptr2tmp);
                  }

                  p[0] = ptr[0];
                  p[1] = ptr[1];
                  p[2] = ptr[2];
                  p[dim] += bi;

                  // Convert transformed point to image coordinates
                  _source[n]->WorldToImage(p[0], p[1], p[2]);

                  // Check whether transformed point is inside volume
                  if ((p[0] > _source_x1[n]) && (p[0] < _source_x2[n]) &&
                      (p[1] > _source_y1[n]) && (p[1] < _source_y2[n]) &&
                      ((_target[n]->GetZ() == 1 && round(p[2]) == 0)
                       ||( _target[n]->GetZ() != 1
                           && p[2] > _source_z1[n] && p[2] < _source_z2[n]))) {

                    // Add sample to metric
                    if(_weight == NULL)
                      weight[n]++;
                    else
                      weight[n] += _weight[_level][n];
                    tmpMetricA[n]->Add(*ptr2target, round(_interpolator[n]->EvaluateInside(p[0], p[1], p[2], t)));
                  }

                  p[0] = ptr[0];
                  p[1] = ptr[1];
                  p[2] = ptr[2];
                  p[dim] -= bi;

                  // Convert transformed point to image coordinates
                  _source[n]->WorldToImage(p[0], p[1], p[2]);

                  // Check whether transformed point is inside volume
                  if ((p[0] > _source_x1[n]) && (p[0] < _source_x2[n]) &&
                      (p[1] > _source_y1[n]) && (p[1] < _source_y2[n]) &&
                      ((_target[n]->GetZ() == 1 && round(p[2]) == 0)
                       ||( _target[n]->GetZ() != 1
                           && p[2] > _source_z1[n] && p[2] < _source_z2[n]))) {

                    // Add sample to metric
                    if(_weight == NULL)
                      weight[n]++;
                    else
                      weight[n] += _weight[_level][n];
                    tmpMetricB[n]->Add(*ptr2target, round(_interpolator[n]->EvaluateInside(p[0], p[1], p[2], t)));
                  }
                }
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
  }

  // Save value of DOF for which we calculate the derivative
  double dof = _affd->Get(index);

  // Evaluate similarity measure
  double similarityA = combine_mysimilarity(tmpMetricA,weight,_numberOfImages);

  // Add penalties
  _affd->Put(index, dof + step);

  // Smoothness
  if (this->_Lregu > 0) {
    similarityA += this->_Lregu * this->LandMarkPenalty(index);
  }
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
  double similarityB = combine_mysimilarity(tmpMetricB,weight,_numberOfImages);

  // Add penalties
  _affd->Put(index, dof - step);

  // Smoothness
  if (this->_Lregu > 0) {
    similarityB += this->_Lregu * this->LandMarkPenalty(index);
  }
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

  delete []weight;

  return similarityA - similarityB;
}

double irtkMultipleImageFreeFormRegistration::EvaluateGradient(float step, float *dx)
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

bool irtkMultipleImageFreeFormRegistration::Read(char *buffer1, char *buffer2, int &level)
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
    return this->irtkMultipleImageRegistration::Read(buffer1, buffer2, level);
  } else {
    return ok;
  }
}

void irtkMultipleImageFreeFormRegistration::Write(ostream &to)
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
    to << "MFFDMode                       = True" << endl;
  } else {
    to << "MFFDMode                       = False" << endl;
  }

  this->irtkMultipleImageRegistration::Write(to);
}

#endif // HAS_VTK

