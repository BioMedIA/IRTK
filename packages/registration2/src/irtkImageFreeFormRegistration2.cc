/*=========================================================================

  Library   : Image Registration Toolkit (IRTK)
  Module    : $Id$
  Copyright : Imperial College, Department of Computing
              Visual Information Processing (VIP), 2009 onwards
  Date      : $Date$
  Version   : $Revision$
  Changes   : $Author$

=========================================================================*/

#include <irtkRegistration2.h>

#ifdef WIN32
#include <time.h>
#else
#include <sys/resource.h>
#endif

irtkImageFreeFormRegistration2::irtkImageFreeFormRegistration2()
{
  // Print debugging information
  this->Debug("irtkImageFreeFormRegistration2::irtkImageFreeFormRegistration2");

  // Default optimization
  _OptimizationMethod = GradientDescent;

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

void irtkImageFreeFormRegistration2::GuessParameter()
{
  int i;
  double xsize, ysize, zsize, spacing;

  if ((_target == NULL) || (_source == NULL)) {
    cerr << "irtkImageFreeFormRegistration2::GuessParameter: Target and source image not found" << endl;
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

void irtkImageFreeFormRegistration2::Initialize()
{
  // Print debugging information
  this->Debug("irtkImageFreeFormRegistration2::Initialize");

  // Initialize base class
  this->irtkImageRegistration2::Initialize();

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
  _mffd->PushLocalTransformation(_affd);
}

void irtkImageFreeFormRegistration2::Initialize(int level)
{
  int i, j, k;

  // Print debugging information
  this->Debug("irtkImageFreeFormRegistration2::Initialize(int)");

  // Initialize base class
  this->irtkImageRegistration2::Initialize(level);

  // Padding of FFD
  irtkPadding(*_target, this->_TargetPadding, _affd);

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

void irtkImageFreeFormRegistration2::Finalize()
{
  // Print debugging information
  this->Debug("irtkImageFreeFormRegistration2::Finalize");

  // Finalize base class
  this->irtkImageRegistration2::Finalize();
}

void irtkImageFreeFormRegistration2::Finalize(int level)
{
  // Print debugging information
  this->Debug("irtkImageFreeFormRegistration2::Finalize(int)");

  // Finalize base class
  this->irtkImageRegistration2::Finalize(level);

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
}

double irtkImageFreeFormRegistration2::EvaluateGradient(double *gradient)
{
  double basis, pos[3], norm;
  int i, j, k, i1, i2, j1, j2, k1, k2, x, y, z, index, index2, index3;

  // Compute gradient with respect to displacements
  this->irtkImageRegistration2::EvaluateGradient(gradient);

  // Start timing
  clock_t start, end;
  double cpu_time_used;
  start = clock();

  // Initialize gradient to 0
  for (i = 0; i < _affd->NumberOfDOFs(); i++) {
    gradient[i] = 0;
  }

  // Loop over control points
  for (z = 0; z < _affd->GetZ(); z++) {
    for (y = 0; y < _affd->GetY(); y++) {
      for (x = 0; x < _affd->GetX(); x++) {

        // Compute DoFs corresponding to the control point
        index  = _affd->LatticeToIndex(x, y, z);
        index2 = index+_affd->GetX()*_affd->GetY()*_affd->GetZ();
        index3 = index+2*_affd->GetX()*_affd->GetY()*_affd->GetZ();

        // Check if any DoF corresponding to the control point is active
        if ((_affd->irtkTransformation::GetStatus(index) == _Active) || (_affd->irtkTransformation::GetStatus(index2) == _Active) || (_affd->irtkTransformation::GetStatus(index3) == _Active)) {

        	// If so, calculate bounding box of control point in image coordinates
          _affd->BoundingBox(_target, index, i1, j1, k1, i2, j2, k2, 1.0);

          // Loop over all voxels in the target (reference) volume
          for (k = k1; k <= k2; k++) {
            for (j = j1; j <= j2; j++) {
              for (i = i1; i <= i2; i++) {

                // Check whether reference point is valid
                if ((_target->Get(i, j, k) >= 0) && (_transformedSource(i, j, k) >= 0)) {

                  // Convert position from voxel coordinates to world coordinates
                  pos[0] = i;
                  pos[1] = j;
                  pos[2] = k;
                  _target->ImageToWorld(pos[0], pos[1], pos[2]);

                  // Convert world coordinates into lattice coordinates
                  _affd->WorldToLattice(pos[0], pos[1], pos[2]);

                  // Compute B-spline tensor product at pos
                  basis = _affd->B(pos[0] - x) * _affd->B(pos[1] - y) * _affd->B(pos[2] - z);

                  // Convert voxel-based gradient into gradient with respect to parameters (chain rule)
                  //
                  // NOTE: This currently assumes that the control points displacements are aligned with the world coordinate displacements
                  //
                  gradient[index]  += basis * _similarityGradient(i, j, k, 0);
                  gradient[index2] += basis * _similarityGradient(i, j, k, 1);
                  gradient[index3] += basis * _similarityGradient(i, j, k, 2);
                }
              }
            }
          }
        }
      }
    }
  }

  // Calculate norm of vector
  norm = 0;
  for (i = 0; i < _affd->NumberOfDOFs(); i++) {
    norm += gradient[i] * gradient[i];
  }

  // Normalize vector
  norm = sqrt(norm);
  if (norm > 0) {
    for (i = 0; i < _affd->NumberOfDOFs(); i++) {
      if (_affd->irtkTransformation::GetStatus(i) == _Active) {
        gradient[i] /= norm;
      } else {
        gradient[i] = 0;
      }
    }
  } else {
    for (i = 0; i < _affd->NumberOfDOFs(); i++) {
      gradient[i] = 0;
    }
  }

  // Stop timing
  end = clock();
  cpu_time_used = ((double) (end - start)) / CLOCKS_PER_SEC;
  //cout << "CPU time for irtkImageFreeFormRegistration2::EvaluateGradient() = " << cpu_time_used << endl;

  return norm;
}

void irtkImageFreeFormRegistration2::Run()
{
  int i, j, k, level, update, updateGradient;
  char buffer[256];
  double *gradient, step, delta, similarity, new_similarity, old_similarity;

  // Print debugging information
  this->Debug("irtkImageRegistration2::Run");

  if (_source == NULL) {
    cerr << "Registration::Run: Filter has no source input" << endl;
    exit(1);
  }

  if (_target == NULL) {
    cerr << "Registration::Run: Filter has no target input" << endl;
    exit(1);
  }

  if (_transformation == NULL) {
    cerr << "irtkImageRegistration2::Run: Filter has no transformation output" << endl;
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
    cout << "Delta: " << delta << " to ";
    cout << delta / pow(2.0, static_cast<double>(_NumberOfSteps[level]-1)) << "\n";

    // Initialize for this level
    this->Initialize(level);

    // Save pre-processed images if we are debugging
    sprintf(buffer, "source_%d.nii.gz", level);
    if (_DebugFlag == true) _source->Write(buffer);
    sprintf(buffer, "target_%d.nii.gz", level);
    if (_DebugFlag == true) _target->Write(buffer);

    // Allocate memory for gradient vector
    gradient = new double[_affd->NumberOfDOFs()];

    // Update image
    update = true;

    // Update gradient
    updateGradient = true;

    // Run the registration filter at this resolution
    for (i = 0; i < _NumberOfSteps[level]; i++) {
      for (j = 0; j < _NumberOfIterations[level]; j++) {

        cout << "Iteration = " << j + 1 << " (out of " << _NumberOfIterations[level];
        cout << "), step size = " << step << endl;

        // Update source image
        if (update == true) {
          this->Update();
          update = false;
        }

        // Compute current metric value
        old_similarity = new_similarity = similarity = this->Evaluate();
        cout << "Current metric value is " << similarity << endl;

        // Compute gradient of similarity metric
        if (updateGradient == true) {
          this->EvaluateGradient(gradient);
          updateGradient = false;
        }

        // Step along gradient direction until no further improvement is necessary
        do {
          new_similarity = similarity;
          for (k = 0; k < _affd->NumberOfDOFs(); k++) {
            _affd->Put(k, _affd->Get(k) + step * gradient[k]);
          }

          // We have just changed the transformation parameters, so definitely need to update
          this->Update();
          update = false;

          // Compute new similarity
          similarity = this->Evaluate();

          if (similarity > new_similarity + _Epsilon) {
            cout << "New metric value is " << similarity << endl;
            updateGradient = true;
          } else {
            // Last step was no improvement, so back track
            for (k = 0; k < _affd->NumberOfDOFs(); k++) {
              _affd->Put(k, _affd->Get(k) - step * gradient[k]);
            }
            update = true;
          }
        } while (similarity > new_similarity + _Epsilon);

        // Check whether we made any improvement or not
        if (new_similarity - old_similarity > _Epsilon) {
          sprintf(buffer, "log_%.3d_%.3d_%.3d.dof", level, i+1, j+1);
          if (_DebugFlag == true) _affd->irtkTransformation::Write(buffer);
        } else {
          sprintf(buffer, "log_%.3d_%.3d_%.3d.dof", level, i+1, j+1);
          if (_DebugFlag == true) _affd->irtkTransformation::Write(buffer);
          break;
        }
      }
      step = step / 2;
      delta = delta / 2.0;
    }

    // Delete gradient
    delete gradient;

    // Do the final cleaning up for this level
    this->Finalize(level);
  }

  // Do the final cleaning up for all levels
  this->Finalize();
}

bool irtkImageFreeFormRegistration2::Read(char *buffer1, char *buffer2, int &level)
{
  int ok = false;

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
    if ((strcmp(buffer2, "false") == 0) || (strcmp(buffer2, "No") == 0)) {
      this->_Subdivision = false;
      cout << "Subdivision is ... false" << endl;
    } else {
      if ((strcmp(buffer2, "true") == 0) || (strcmp(buffer2, "Yes") == 0)) {
        this->_Subdivision = true;
        cout << "Subdivision is ... true" << endl;
      } else {
        cerr << "Can't read boolean value = " << buffer2 << endl;
        exit(1);
      }
    }
    ok = true;
  }

  if (ok == false) {
    return this->irtkImageRegistration2::Read(buffer1, buffer2, level);
  } else {
    return ok;
  }
}

void irtkImageFreeFormRegistration2::Write(ostream &to)
{
  to << "\n#\n# Non-rigid registration parameters\n#\n\n";
  to << "Lambda1                           = " << this->_Lambda1 << endl;
  to << "Lambda2                           = " << this->_Lambda2 << endl;
  to << "Lambda3                           = " << this->_Lambda3 << endl;
  to << "Control point spacing in X        = " << this->_DX << endl;
  to << "Control point spacing in Y        = " << this->_DY << endl;
  to << "Control point spacing in Z        = " << this->_DZ << endl;
  if (_Subdivision == true) {
    to << "Subdivision                       = true" << endl;
  } else {
    to << "Subdivision                       = false" << endl;
  }

  this->irtkImageRegistration2::Write(to);
}
