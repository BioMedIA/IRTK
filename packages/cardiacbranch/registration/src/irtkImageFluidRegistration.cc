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

#include <irtkGradientDescentConstrainedOptimizer.h>

// Used as temporary memory for transformed intensities
irtkGreyImage *_tmpImage;

// Used as temporary memory for transformed intensities
irtkGreyImage *_tmpSource;

// Used as temporary memory for metric
irtkSimilarityMetric *_tmpMetricA, *_tmpMetricB;

// Used as lookup table for transformed coordinates including level n. This
// lookup table needs to be calculated each time a control point has been
// modified.
float *_affdLookupTable;

// Used as lookup table for the contribution of each control point. This
// lookup table needs to be calculated only once.
float *_localLookupTable;

irtkImageFluidRegistration::irtkImageFluidRegistration()
{
  // Print debugging information
  this->Debug("irtkImageFluidRegistration::irtkImageFluidRegistration");

  // Default optimization
  _OptimizationMethod = GradientDescentConstrained;

  // Default speedup factor
  _SpeedupFactor = 1;

  // Default parameters for non-rigid registration
  _DX          = 20;
  _DY          = 20;
  _DZ          = 20;
  _NumberOfTimeSteps = 10;
}

void irtkImageFluidRegistration::GuessParameter()
{
  cerr << "irtkImageFluidRegistration::GuessParameter: Not yet implemented" << endl;
  exit(1);
}

void irtkImageFluidRegistration::Initialize()
{
  // Print debugging information
  this->Debug("irtkImageFluidRegistration::Initialize");

  // Initialize base class
  this->irtkImageRegistration::Initialize();

  // Pointer to fluid FFD
  _mffd = (irtkFluidFreeFormTransformation *)_transformation;

  // Initialize pointers
  _tmpImage         = NULL;
  _affdLookupTable  = NULL;
  _localLookupTable = new float [FFDLOOKUPTABLESIZE];
}

void irtkImageFluidRegistration::Initialize(int level)
{
  int i, j;
  double u;

  // Print debugging information
  this->Debug("irtkImageFluidRegistration::Initialize(int)");

  // Initialize base class
  this->irtkImageRegistration::Initialize(level);

  // Tell optimizer what the maximum allowed deformation is
  irtkGradientDescentConstrainedOptimizer *o = dynamic_cast<irtkGradientDescentConstrainedOptimizer *>(_optimizer);
  if (o == NULL) {
    cerr << "irtkImageFluidRegistration::Initialize: Optimizer must be of type irtkGradientDescentConstrainedOptimizer" << endl;
    exit(1);
  }
  o->SetLimits(this->_DX * 0.4);

  // Allocate memory for metric
  _tmpMetricA = irtkSimilarityMetric::New(_metric);
  _tmpMetricB = irtkSimilarityMetric::New(_metric);

  // Allocate tmp image
  _tmpImage  = new irtkGreyImage(*_target);

  // Allocate tmp source image
  _tmpSource = new irtkGreyImage(*_source);

  // Allocate memory for lookup table for single-level FFD
  _affdLookupTable  = new float[_target->GetNumberOfVoxels()*3];

  for (i = 0; i < FFDLOOKUPTABLESIZE; i++) {
    u = i / (FFDLOOKUPTABLESIZE / 4.0);
    j = (int)floor(u);
    u = u - j;
    _localLookupTable[i] = _affd->B(j, 1-u);
  }
}

void irtkImageFluidRegistration::Finalize()
{
  // Print debugging information
  this->Debug("irtkImageFluidRegistration::Finalize");

  // Push local transformation back on transformation stack
  _mffd->PushLocalTransformation(_affd);

  // Finalize base class
  this->irtkImageRegistration::Finalize();

  delete []_localLookupTable;
}

void irtkImageFluidRegistration::Finalize(int level)
{
  // Print debugging information
  this->Debug("irtkImageFluidRegistration::Finalize(int)");

  // Finalize base class
  this->irtkImageRegistration::Finalize(level);

  delete _tmpImage;
  delete _tmpSource;
  delete _tmpMetricA;
  delete _tmpMetricB;
  delete []_affdLookupTable;
}

void irtkImageFluidRegistration::UpdateLUT()
{
  int i, j, k;
  double x, y, z;
  float *ptr2affd;

  // Print debugging information
  this->Debug("irtkImageFluidRegistration::UpdateLUT");

  ptr2affd = _affdLookupTable;
  for (k = 0; k < _target->GetZ(); k++) {
    for (j = 0; j < _target->GetY(); j++) {
      for (i = 0; i < _target->GetX(); i++) {
        x = i;
        y = j;
        z = k;
        _target->ImageToWorld(x, y, z);
        _affd->Transform(x, y, z);
        ptr2affd[0] = x;
        ptr2affd[1] = y;
        ptr2affd[2] = z;
        ptr2affd++;
        ptr2affd++;
        ptr2affd++;
      }
    }
  }
}

double irtkImageFluidRegistration::Evaluate()
{
  // Image coordinates
  int i, j, k;
  // World coordinates
  double x, y, z;
  // Pointer to reference data
  irtkGreyPixel *ptr2target;
  irtkGreyPixel *ptr2tmp;

  // Print debugging information
  this->Debug("irtkImageFluidRegistration::Evaluate");

  // Initialize metric
  _metric->Reset();

  // Loop over all voxels in the target (reference) volume
  ptr2target = _target->GetPointerToVoxels();
  ptr2tmp    = _tmpImage->GetPointerToVoxels();
  for (k = 0; k < _target->GetZ(); k++) {
    for (j = 0; j < _target->GetY(); j++) {
      for (i = 0; i < _target->GetX(); i++) {
        // Check whether reference point is valid
        if (*ptr2target >= 0) {
          x = i;
          y = j;
          z = k;
          _target->ImageToWorld(x, y, z);
          _affd->Transform(x, y, z);
          _source->WorldToImage(x, y, z);
          // Check whether transformed point is inside volume
          if ((x > _source_x1) && (x < _source_x2) &&
              (y > _source_y1) && (y < _source_y2) &&
              (z > _source_z1) && (z < _source_z2)) {
            // Add sample to metric
            *ptr2tmp =  round(_interpolator->EvaluateInside(x, y, z));
            _metric->Add(*ptr2target, *ptr2tmp);
          } else {
            *ptr2tmp = -1;
          }
        }
        // Increment pointers to next voxel
        ptr2tmp++;
        ptr2target++;
      }
    }
  }

  // Evaluate similarity measure
  double similarity = _metric->Evaluate();

  // Return similarity measure + penalty terms
  return similarity;
}

double irtkImageFluidRegistration::EvaluateDerivative(int index, double step)
{
  float *ptr;
  double bi, bj, bk, dx, dy, dz, p[3];
  int i, j, k, i1, i2, j1, j2, k1, k2, dim;
  irtkPoint p1, p2;
  irtkGreyPixel *ptr2target, *ptr2tmp;

  // Print debugging information
  this->Debug("irtkImageFluidRegistration::EvaluateDerivative(int, double)");

  // Initialize metrics for forward and backward derivative steps
  _tmpMetricA->Reset(_metric);
  _tmpMetricB->Reset(_metric);

  // Calculate bounding box of control point in world coordinates
  _affd->BoundingBox(index, p1, p2);
  _target->WorldToImage(p1);
  _target->WorldToImage(p2);

  // Calculate bounding box of control point in image coordinates
  _affd->BoundingBox(_target, index, i1, j1, k1, i2, j2, k2);

  // Calculate incremental changes in lattice coordinates when looping
  // over target
  dx = (FFDLOOKUPTABLESIZE-1)/(p2._x-p1._x);
  dy = (FFDLOOKUPTABLESIZE-1)/(p2._y-p1._y);
  dz = (FFDLOOKUPTABLESIZE-1)/(p2._z-p1._z);

  // Calculate whether this DOF corresponds to x, y or z-displacement
  dim = int(index / (_affd->GetX()*_affd->GetY()*_affd->GetZ()));

  // Loop over all voxels in the target (reference) volume
  for (k = k1; k <= k2; k++) {
    bk = step * _localLookupTable[round(dz*(k-p1._z))];
    for (j = j1; j <= j2; j++) {
      ptr2target = _target->GetPointerToVoxels(i1, j, k);
      ptr        = &(_affdLookupTable[3*_target->VoxelToIndex(i1, j, k)]);
      bj = bk * _localLookupTable[round(dy*(j-p1._y))];
      ptr2tmp  = _tmpImage->GetPointerToVoxels(i1, j, k);
      for (i = i1; i <= i2; i++) {

        // Check whether reference point is valid
        if (*ptr2target >= 0) {
          bi = bj * _localLookupTable[round(dx*(i-p1._x))];

          // Delete old samples from both metrics
          if (*ptr2tmp != -1) {
            _tmpMetricA->Delete(*ptr2target, *ptr2tmp);
            _tmpMetricB->Delete(*ptr2target, *ptr2tmp);
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
            _tmpMetricA->Add(*ptr2target, round(_interpolator->EvaluateInside(p[0], p[1], p[2])));
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
            _tmpMetricB->Add(*ptr2target, round(_interpolator->EvaluateInside(p[0], p[1], p[2])));
          }
        }

        // Increment pointers to next voxel
        ptr2target++;
        ptr2tmp++;
        ptr += 3;
      }
    }
  }

  // Evaluate similarity measure
  double similarityA = _tmpMetricA->Evaluate();

  // Evaluate similarity measure
  double similarityB = _tmpMetricB->Evaluate();

  return similarityA - similarityB;
}

double irtkImageFluidRegistration::EvaluateGradient(float step, float *dx)
{
  int i;
  double norm;

  // Update lookup table
  this->UpdateLUT();

  for (i = 0; i < _affd->NumberOfDOFs(); i++) {
    if (_affd->irtkTransformation::GetStatus(i) == _Active) {
      dx[i] = this->EvaluateDerivative(i, step);
    } else {
      dx[i] = 0;
    }
  }

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

void irtkImageFluidRegistration::Run()
{
  int i, j, k, level;
  char buffer[256];
  double step, epsilon;

  // Print debugging information
  this->Debug("irtkImageFluidRegistration::Run");

  if (_source == NULL) {
    cerr << "irtkImageFluidRegistration::Run: Filter has no source input" << endl;
    exit(1);
  }

  if (_target == NULL) {
    cerr << "irtkImageFluidRegistration::Run: Filter has no target input" << endl;
    exit(1);
  }

  if (_transformation == NULL) {
    cerr << "irtkImageFluidRegistration::Run: Filter has no transformation output" << endl;
    exit(1);
  }

  // Do the initial set up for all levels
  this->Initialize();

  this->Write(cout);

  for (level = _NumberOfLevels-1; level >= 0; level--) {

    // Initial step size
    step = _LengthOfSteps[level];

    // Print resolution level
    cout << "Resolution level no. " << level+1 << " (step sizes ";
    cout << step << " to " << step / pow(2.0, static_cast<double>(_NumberOfSteps[level]-1)) << ")\n";

#ifdef HISTORY
    history->Clear();
#endif

    // Initialize for this level
    this->Initialize(level);

    // Save pre-processed images if we are debugging
    sprintf(buffer, "source_%d.gipl", level);
    if (_DebugFlag == True) _source->Write(buffer);
    sprintf(buffer, "target_%d.gipl", level);
    if (_DebugFlag == True) _target->Write(buffer);

    // Run the registration filter
    for (k = 0; k < this->_NumberOfTimeSteps; k++) {

      // Initial step size
      step = _LengthOfSteps[level];

      // Create new FFD
      _affd = new irtkBSplineFreeFormTransformation(*_target, this->_DX, this->_DY, this->_DZ);

      // Tell optimizer which transformation to optimize
      _optimizer->SetTransformation(_affd);

      // Create image transformation
      irtkImageTransformation *imagetransformation = new irtkImageTransformation;
      imagetransformation->SetInput (_source, _mffd);
      imagetransformation->SetOutput(_tmpSource);
      imagetransformation->PutInterpolator(_interpolator);

      // Transform image
      imagetransformation->Run();

      // Delete image transformation
      delete imagetransformation;

      // Setup interpolation for the source image
      _interpolator->SetInput(_tmpSource);
      _interpolator->Initialize();

      // Calculate the source image domain in which we can interpolate
      _interpolator->Inside(_source_x1, _source_y1, _source_z1,
                            _source_x2, _source_y2, _source_z2);

      for (i = 0; i < _NumberOfSteps[level]; i++) {
        _optimizer->SetStepSize(step);
        _optimizer->SetEpsilon(_Epsilon);
        for (j = 0; j < _NumberOfIterations[level]; j++) {

          cout << "Time Step = " << k+1 << ", iteration = " << j + 1 << " (out of " << this->_NumberOfIterations[level];
          cout << "), step size = " << step << endl;

          // Optimize at lowest level of resolution
          epsilon = _optimizer->Run();

          // Check whether we made any improvement or not
          if (epsilon > _Epsilon) {
            sprintf(buffer, "log_%d_%d_%d.dof", level, i+1, j+1);
            if (_DebugFlag == True) _transformation->Write(buffer);
            this->Print();
          } else {
            sprintf(buffer, "log_%d_%d_%d.dof", level, i+1, j+1);
            if (_DebugFlag == True) _transformation->Write(buffer);
            this->Print();
            break;
          }
        }
        step = step / 2;
      }
      // Add local transformation
      _mffd->PushLocalTransformation(_affd);
    }
    // Do the final cleaning up for this level
    this->Finalize(level);
  }

#ifdef HISTORY
  history->Print();
#endif

  // Do the final cleaning up for all levels
  this->Finalize();
}

Bool irtkImageFluidRegistration::Read(char *buffer1, char *buffer2, int &level)
{
  int ok = False;

  if (strstr(buffer1, "No. of time steps") != NULL) {
    this->_NumberOfTimeSteps = atoi(buffer2);
    cout << "No. of time steps are ... " << this->_NumberOfTimeSteps << endl;
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

  if (ok == False) {
    return this->irtkImageRegistration::Read(buffer1, buffer2, level);
  } else {
    return ok;
  }
}

void irtkImageFluidRegistration::Write(ostream &to)
{
  to << "\n#\n# Non-rigid fluid registration parameters\n#\n\n";
  to << "No. of time steps                 = " << this->_NumberOfTimeSteps << endl;
  to << "Control point spacing in X        = " << this->_DX << endl;
  to << "Control point spacing in Y        = " << this->_DY << endl;
  to << "Control point spacing in Z        = " << this->_DZ << endl;

  this->irtkImageRegistration::Write(to);
}
