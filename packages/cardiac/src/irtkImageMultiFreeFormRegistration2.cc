/*=========================================================================

  Library   : Image Registration Toolkit (IRTK)
  Module    : $Id$
  Copyright : Imperial College, Department of Computing
              Visual Information Processing (VIP), 2009 onwards
  Date      : $Date$
  Version   : $Revision$
  Changes   : $Author$

=========================================================================*/

#include <irtkCardiac.h>

#include <irtkGradientImageFilter.h>

#ifdef WIN32
#include <time.h>
#else
#include <sys/resource.h>
#endif

irtkGreyImage *tmp_segmentation;

irtkImageMultiFreeFormRegistration2::irtkImageMultiFreeFormRegistration2()
{
  // Print debugging information
  this->Debug("irtkImageMultiFreeFormRegistration2::irtkImageMultiFreeFormRegistration2");

  // Default optimization
  _OptimizationMethod = GradientDescent;

  // Default parameters for non-rigid registration
  _NumberOfLabels = 0;
  _Lambda1     = 0;
  _Lambda2     = 0;
  _Lambda3     = 0;
  _LambdaC     = NULL;
  _DX          = 20;
  _DY          = 20;
  _DZ          = 20;
  _Subdivision = true;
  _Mode        = RegisterXYZ;
  _MFFDMode    = true;
  _adjugate     = NULL;
  _determine   = NULL;
}

void irtkImageMultiFreeFormRegistration2::SetInput(irtkGreyImage *target, irtkGreyImage *source, irtkGreyImage *segmentation)
{

    //check dimension and resolution
    if(target->GetX() != segmentation->GetX()
        || target->GetY() != segmentation->GetY()
        || target->GetZ() != segmentation->GetZ()
        || target->GetT() != segmentation->GetT()){
            cerr << "target and segmentation size does not correspond!" <<endl;
            exit(0);
    }
    _target = target;
    _source = source;
    _segmentation = segmentation;

    irtkGreyPixel current,min,max,*ptr;
    int i,sum,*map;
    _segmentation->GetMinMax(&min,&max);

    map = new int[max-min+1];

    for (current = min; current <= max; current ++){
        ptr = _segmentation->GetPointerToVoxels();
        sum = 0;
        for (i = 0; i < _segmentation->GetNumberOfVoxels(); i++){
            if(*ptr == current) sum++;
            ptr++;
        }
        if(sum > 0){
            _NumberOfLabels++;
            map[current - min] = _NumberOfLabels - 1;
        }else{
            map[current - min] = -1;
        }
    }
    ptr = _segmentation->GetPointerToVoxels();
    for (i = 0; i < _segmentation->GetNumberOfVoxels(); i++){
        *ptr = map[*ptr];
        ptr++;
    }
    delete []map;
}

void irtkImageMultiFreeFormRegistration2::SetOutput(irtkTransformation **mtransformation)
{
    // Print debugging information
    this->Debug("irtkImageFreeFormRegistration::SetOutput");

    if(_segmentation == NULL){
        cerr << "please set input first" <<endl;
        exit(0);
    }

    for(int i = 0; i < _NumberOfLabels; i++)
    if (strcmp(mtransformation[i]->NameOfClass(),
        "irtkMultiLevelFreeFormTransformation") != 0) {
            cerr << "irtkImageFreeFormRegistration::SetOutput: Transformation must be "
                << "irtkMultiLevelFreeFormTransformation" << endl;
            exit(0);
    }
    _mtransformation = mtransformation;
    _transformation = mtransformation[0];
}

void irtkImageMultiFreeFormRegistration2::GuessParameter()
{
  int i;
  double xsize, ysize, zsize, spacing;

  if ((_target == NULL) || (_source == NULL)) {
    cerr << "irtkImageMultiFreeFormRegistration2::GuessParameter: Target and source image not found" << endl;
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

void irtkImageMultiFreeFormRegistration2::Initialize()
{
  // Print debugging information
  this->Debug("irtkImageMultiFreeFormRegistration2::Initialize");

  // Initialize base class
  this->irtkImageRegistration2::Initialize();

  _mffd = new irtkMultiLevelFreeFormTransformation*[_NumberOfLabels];
  _affd = new irtkBSplineFreeFormTransformation*[_NumberOfLabels];
  _adjugate = new irtkMatrix*[_NumberOfLabels];
  _determine = new double*[_NumberOfLabels];

  for(int i = 0; i < _NumberOfLabels; i++){

  // Pointer to multi-level FFD
  _mffd[i] = (irtkMultiLevelFreeFormTransformation *)_mtransformation[i];

  // Create FFD
  if (_MFFDMode == false) {
    if (_mffd[i]->NumberOfLevels() == 0) {
      _affd[i] = new irtkBSplineFreeFormTransformation(*_target, this->_DX, this->_DY, this->_DZ);
    } else {
      _affd[i] = (irtkBSplineFreeFormTransformation *)_mffd[i]->PopLocalTransformation();
    }
  } else {
    _affd[i] = new irtkBSplineFreeFormTransformation(*_target, this->_DX, this->_DY, this->_DZ);
  }
  _mffd[i]->PushLocalTransformation(_affd[i]);
  _adjugate[i] = new irtkMatrix[_affd[i]->NumberOfDOFs()/3];
  _determine[i] = new double[_affd[i]->NumberOfDOFs()/3];
  }
}

void irtkImageMultiFreeFormRegistration2::Initialize(int level)
{
  int i, j, k, l;
  double dx, dy, dz, temp;

  // Print debugging information
  this->Debug("irtkImageMultiFreeFormRegistration2::Initialize(int)");

  // Copy segmentation to temp space
  tmp_segmentation = new irtkGreyImage(*_segmentation);
  swap(tmp_segmentation, _segmentation);
  // Resample segmentation if necessary
  _segmentation->GetPixelSize(&dx, &dy, &dz);
  temp = fabs(_TargetResolution[0][0]-dx) + fabs(_TargetResolution[0][1]-dy) + fabs(_TargetResolution[0][2]-dz);

  if (level > 0 || temp > 0.000001) {
      cout << "Resampling segmentation ... "; cout.flush();
      irtkNearestNeighborInterpolateImageFunction sinter;
      // Create resampling filter
      irtkResampling<irtkGreyPixel> resample(_TargetResolution[level][0],
          _TargetResolution[level][1],
          _TargetResolution[level][2]);
      resample.SetInput (_segmentation);
      resample.SetOutput(_segmentation);
      resample.SetInterpolator(&sinter);
      resample.Run();
      cout << "done" << endl;
  }

  _segmentation->Write("testseg.nii.gz");

  // Initialize base class
  this->irtkImageRegistration2::Initialize(level);

  for(l = 0; l < _NumberOfLabels; l ++){

  // Padding of FFD
  irtkPadding(*_target, this->_TargetPadding, _affd[l]);
  irtkGreyImage paddingimage(*_segmentation);
  irtkGreyPixel *paddingptr = paddingimage.GetPointerToVoxels();
  for (i = 0; i < paddingimage.GetNumberOfVoxels(); i ++){
      if(*paddingptr != l)
          *paddingptr = -1;
      paddingptr ++;
  }
  irtkPadding(paddingimage, -1, _affd[l]);


  // Register in the x-direction only
  if (_Mode == RegisterX) {
    for (i = 0; i < _affd[l]->GetX(); i++) {
      for (j = 0; j < _affd[l]->GetY(); j++) {
        for (k = 0; k < _affd[l]->GetZ(); k++) {
          _Status sx, sy, sz;
          _affd[l]->GetStatus(i, j, k, sx, sy, sz);
          _affd[l]->PutStatus(i, j, k, sx, _Passive, _Passive);
        }
      }
    }
  }

  // Register in the y-direction only
  if (_Mode == RegisterY) {
    for (i = 0; i < _affd[l]->GetX(); i++) {
      for (j = 0; j < _affd[l]->GetY(); j++) {
        for (k = 0; k < _affd[l]->GetZ(); k++) {
          _Status sx, sy, sz;
          _affd[l]->GetStatus(i, j, k, sx, sy, sz);
          _affd[l]->PutStatus(i, j, k, _Passive, sy, _Passive);
        }
      }
    }
  }

  // Register in the x- and y-direction only
  if (_Mode == RegisterXY) {
    for (i = 0; i < _affd[l]->GetX(); i++) {
      for (j = 0; j < _affd[l]->GetY(); j++) {
        for (k = 0; k < _affd[l]->GetZ(); k++) {
          _Status sx, sy, sz;
          _affd[l]->GetStatus(i, j, k, sx, sy, sz);
          _affd[l]->PutStatus(i, j, k, sx, sy, _Passive);
        }
      }
    }
  }
  }
}

void irtkImageMultiFreeFormRegistration2::Finalize()
{
  // Print debugging information
  this->Debug("irtkImageMultiFreeFormRegistration2::Finalize");

  // Finalize base class
  this->irtkImageRegistration2::Finalize();

  delete []_mffd;
  delete []_affd;
  for(int i = 0; i < _NumberOfLabels; i++){
      delete []_adjugate[i];
      delete []_determine[i];
  }
  delete []_determine;
  delete []_adjugate;
}

void irtkImageMultiFreeFormRegistration2::Finalize(int level)
{
  // Print debugging information
  this->Debug("irtkImageMultiFreeFormRegistration2::Finalize(int)");

  // Finalize base class
  this->irtkImageRegistration2::Finalize(level);

  swap(tmp_segmentation, _segmentation);

  delete tmp_segmentation;

  // Check if we are not at the lowest level of resolution
  if (level != 0) {
    for(int i = 0; i < _NumberOfLabels; i++){
    if (this->_Subdivision == true) {
      _affd[i]->Subdivide();
    } else {
      // Create new FFD
      _affd[i] = new irtkBSplineFreeFormTransformation(*_target,
              this->_DX / pow(2.0, this->_NumberOfLevels-level),
              this->_DY / pow(2.0, this->_NumberOfLevels-level),
              this->_DZ / pow(2.0, this->_NumberOfLevels-level));
      // Push local transformation back on transformation stack
      _mffd[i]->PushLocalTransformation(_affd[i]);
    }
    delete []_adjugate[i];
    delete []_determine[i];
    _adjugate[i] = new irtkMatrix[_affd[i]->NumberOfDOFs()/3];
    _determine[i] = new double[_affd[i]->NumberOfDOFs()/3];
    }
  }
}

double irtkImageMultiFreeFormRegistration2::VolumePreservationPenalty()
{
    int k, l;
    double penalty, jacobian, lambdac, count;
    irtkMatrix jac;
    penalty = 0;
    count = 0;
    for (l = 0; l < _NumberOfLabels; l++){
        if(_LambdaC == NULL){
            lambdac = 1;
        }else{
            lambdac = _LambdaC[l];
        }
        for (k = 0; k < _affd[l]->NumberOfDOFs()/3; k++) {
            if(_affd[l]->GetStatus(k) == _Active){
                // Determinant of Jacobian of deformation derivatives
                jacobian = _determine[l][k];
                penalty += lambdac*pow(log(jacobian),2);
                count += lambdac;
            }
        }
    }

    // Normalize sum by number of DOFs
    return -penalty / count;
}

void irtkImageMultiFreeFormRegistration2::VolumePreservationPenalty(int l,int index, double *drv)
{
    int i, j, k, m, n, o, i1, j1, k1, i2, j2, k2, count, index1;
    double jacobian;
    irtkMatrix jac,det_drv[3];

    _affd[l]->IndexToLattice(index, i, j, k);
    o = i;
    m = j;
    n = k;
    count = 0;
    k1 = (k-1)>0?(k-1):0;
    j1 = (j-1)>0?(j-1):0;
    i1 = (i-1)>0?(i-1):0;
    k2 = (k+2) < _affd[l]->GetZ()? (k+2) : _affd[l]->GetZ();
    j2 = (j+2) < _affd[l]->GetY()? (j+2) : _affd[l]->GetY();
    i2 = (i+2) < _affd[l]->GetX()? (i+2) : _affd[l]->GetX();

    for(i=0;i<3;i++)
        drv[i] = 0;

    for (k = k1; k < k2; k++) {
        for (j = j1; j < j2; j++) {
            for (i = i1; i < i2; i++) {
                if(((k == n && j == m) 
                    || (k == n && i == o)
                    || (j == m && i == o))
                    &&(k != n || j != m || i != o)){
                        index1 = _affd[l]->LatticeToIndex(i,j,k);
                        // Torsten Rohlfing et al. MICCAI'01 (w/o scaling correction):
                        // Calculate jacobian
                        irtkMatrix jac = _adjugate[l][index1];

                        // find jacobian derivatives
                        jacobian = _determine[l][index1];

                        // if jacobian < 0
                        jacobian = (2.0*log(jacobian))/jacobian;

                        _affd[l]->JacobianDetDerivative(det_drv,i-o,j-m,k-n);

                        double tmpdrv[3];

                        for(int p = 0; p < 3; p++){
                            // trace * adjugate * derivative
                            tmpdrv[p] = jac(0,p)*det_drv[p](p,0) + jac(1,p)*det_drv[p](p,1) + jac(2,p)*det_drv[p](p,2);
                            // * rest pf the vplume preservatipn derivative
                            drv[p] += (jacobian*tmpdrv[p]);
                        }
                        count ++;
                }
            }
        }
    }

    for(o=0;o<3;o++)
        drv[l] = -drv[l]/count;

    return;
}

double irtkImageMultiFreeFormRegistration2::Evaluate()
{
    double similarity;

    // Evaluate similarity
    similarity = this->irtkImageRegistration2::Evaluate();

    // Add penalty for smoothness
    if (this->_Lambda1 > 0) {
    }
    // Add penalty for volume preservation
    if (this->_Lambda2 > 0) {
        similarity += this->_Lambda2*this->VolumePreservationPenalty();
    }
    // Add penalty for topology preservation
    if (this->_Lambda3 > 0) {
    }

    //Return similarity measure + penalty terms
    return similarity;
}

double irtkImageMultiFreeFormRegistration2::EvaluateGradient(int l, double *gradient)
{
  double basis, pos[3], norm;
  int i, j, k, i1, i2, j1, j2, k1, k2, x, y, z, index, index2, index3;

  // Initialize gradient to 0
  for (i = 0; i < _affd[l]->NumberOfDOFs(); i++) {
    gradient[i] = 0;
  }

  // Loop over control points
  for (z = 0; z < _affd[l]->GetZ(); z++) {
    for (y = 0; y < _affd[l]->GetY(); y++) {
      for (x = 0; x < _affd[l]->GetX(); x++) {

        // Compute DoFs corresponding to the control point
        index  = _affd[l]->LatticeToIndex(x, y, z);
        index2 = index+_affd[l]->GetX()*_affd[l]->GetY()*_affd[l]->GetZ();
        index3 = index+2*_affd[l]->GetX()*_affd[l]->GetY()*_affd[l]->GetZ();

        // Check if any DoF corresponding to the control point is active
        if ((_affd[l]->GetStatus(index) == _Active) 
            || (_affd[l]->GetStatus(index2) == _Active) 
            || (_affd[l]->GetStatus(index3) == _Active)) {

            // If so, calculate bounding box of control point in image coordinates
            _affd[l]->BoundingBox(_target, index, i1, j1, k1, i2, j2, k2, 1.0);

            // Loop over all voxels in the target (reference) volume
            for (k = k1; k <= k2; k++) {
                for (j = j1; j <= j2; j++) {
                    for (i = i1; i <= i2; i++) {

                        // Check whether reference point is valid
                        if ((_target->Get(i, j, k) >= 0) && (_transformedSource(i, j, k) >= 0) && (_segmentation->Get(i,j,k) == l)) {

                            // Convert position from voxel coordinates to world coordinates
                            pos[0] = i;
                            pos[1] = j;
                            pos[2] = k;
                            _target->ImageToWorld(pos[0], pos[1], pos[2]);

                            // Convert world coordinates into lattice coordinates
                            _affd[l]->WorldToLattice(pos[0], pos[1], pos[2]);

                            // Compute B-spline tensor product at pos
                            basis = _affd[l]->B(pos[0] - x) * _affd[l]->B(pos[1] - y) * _affd[l]->B(pos[2] - z);

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
            // Add penalty for smoothness
            if (this->_Lambda1 > 0) {

            }
            // Add penalty for volume preservation
            if (this->_Lambda2 > 0) {
                // Get inverted jacobian
                double det_dev[3];
                this->VolumePreservationPenalty(l,index,det_dev);
                if(_LambdaC == NULL){
                    gradient[index]  += this->_Lambda2 * det_dev[0];
                    gradient[index2] += this->_Lambda2 * det_dev[1];
                    gradient[index3] += this->_Lambda2 * det_dev[2];
                }else{
                    gradient[index]  += this->_Lambda2 * this->_LambdaC[l] * det_dev[0];
                    gradient[index2] += this->_Lambda2 * this->_LambdaC[l] * det_dev[1];
                    gradient[index3] += this->_Lambda2 * this->_LambdaC[l] * det_dev[2];
                }
            }
            // Add penalty for topology preservation
            if (this->_Lambda3 > 0) {

            }
        }
      }
    }
  }

  // Calculate norm of vector
  norm = 0;
  for (i = 0; i < _affd[l]->NumberOfDOFs(); i++) {
    norm += gradient[i] * gradient[i];
  }

  // Normalize vector
  norm = sqrt(norm);
  if (norm > 0) {
    for (i = 0; i < _affd[l]->NumberOfDOFs(); i++) {
      if (_affd[l]->irtkTransformation::GetStatus(i) == _Active) {
        gradient[i] /= norm;
      } else {
        gradient[i] = 0;
      }
    }
  } else {
    for (i = 0; i < _affd[l]->NumberOfDOFs(); i++) {
      gradient[i] = 0;
    }
  }

  return norm;
}

void irtkImageMultiFreeFormRegistration2::Update()
{

  // Image transformation
  irtkImageTransformation _imagetransformation;

  // Generate transformed tmp image
  _transformedSource = *_target;

  int i, j, k;
  double x, y, z;
  for (k = 0; k < _target->GetZ(); k++) {
    for (j = 0; j < _target->GetY(); j++) {
      for (i = 0; i < _target->GetX(); i++) {
        x = i;
        y = j;
        z = k;
        _target->ImageToWorld(x, y, z);
        _mtransformation[_segmentation->Get(i,j,k)]->Transform(x, y, z);
        _source->WorldToImage(x, y, z);
        // Check whether transformed point is inside volume
        if ((x > _source_x1) && (x < _source_x2) &&
            (y > _source_y1) && (y < _source_y2) &&
            (z > _source_z1) && (z < _source_z2)) {
          _transformedSource(i, j, k) = round(_interpolator->EvaluateInside(x, y, z, 0));
        } else {
        	_transformedSource(i, j, k) = -1;
        }
      }
    }
  }

  // Compute gradient of source image
  irtkGradientImageFilter<double> gradient(irtkGradientImageFilter<double>::GRADIENT_VECTOR);
  gradient.SetInput (&_transformedSource);
  gradient.SetOutput(&_transformedSourceGradient);
  gradient.Run();

  if(_Lambda2 > 0){
      int index,l;
      double jacobian;
      for(l = 0; l < this->_NumberOfLabels; l++){
          // Update Jacobian and jacobian determine;
          for (index = 0; index < _affd[l]->NumberOfDOFs()/3; index++) {
              _affd[l]->IndexToLattice(index,i,j,k);
              x = i;
              y = j;
              z = k;
              _affd[l]->LatticeToWorld(x,y,z);
              irtkMatrix jac;
              jac.Initialize(3,3);
              _mffd[l]->LocalJacobian(jac,x,y,z);
              jac.Adjugate(jacobian);
              if(jacobian < 0.00001) jacobian = 0.00001;
              _determine[l][index] = jacobian;
              _adjugate[l][index] = jac;
          }
      }
  }

}

void irtkImageMultiFreeFormRegistration2::Run()
{
  int i, j, k, l, level, update, updateGradient;
  char buffer[256];
  double **gradient, step, delta, similarity, new_similarity, old_similarity;

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

  gradient = new double*[_NumberOfLabels];

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

    for(l = 0; l < _NumberOfLabels; l++){
        // Allocate memory for gradient vector
        gradient[l] = new double[_affd[l]->NumberOfDOFs()];
    }

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
            // Evaluate similarity gradient
            this->irtkImageRegistration2::EvaluateGradient(gradient[0]);
            for(l = 0; l < _NumberOfLabels; l++){
                this->EvaluateGradient(l,gradient[l]);
            }
          updateGradient = false;
        }

        // Step along gradient direction until no further improvement is necessary
        do {
          new_similarity = similarity;
          for(l = 0; l < _NumberOfLabels; l++){
          for (k = 0; k < _affd[l]->NumberOfDOFs(); k++) {
            _affd[l]->Put(k, _affd[l]->Get(k) + step * gradient[l][k]);
          }
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
              for(l = 0; l < _NumberOfLabels; l++){
                  for (k = 0; k < _affd[l]->NumberOfDOFs(); k++) {
                      _affd[l]->Put(k, _affd[l]->Get(k) - step * gradient[l][k]);
                  }
              }
            update = true;
          }
        } while (similarity > new_similarity + _Epsilon);

        // Check whether we made any improvement or not
        if (new_similarity - old_similarity > _Epsilon) {
            if (_DebugFlag == true) 
            {for(l = 0; l < _NumberOfLabels; l++){
                sprintf(buffer, "log_%.3d_%.3d_%.3d%d.dof", level, i+1, j+1, l);
            _affd[l]->irtkTransformation::Write(buffer);
            }
            }
        } else {
            if (_DebugFlag == true) 
            {for(l = 0; l < _NumberOfLabels; l++){
                sprintf(buffer, "log_%.3d_%.3d_%.3d%d.dof", level, i+1, j+1, l);
            _affd[l]->irtkTransformation::Write(buffer);
            }
            }
            break;
        }
      }
      step = step / 2;
      delta = delta / 2.0;
    }

    // Delete gradient
    for(l = 0; l < _NumberOfLabels; l++){
        delete []gradient[l];
    }

    // Do the final cleaning up for this level
    this->Finalize(level);
  }

  delete []gradient;

  // Do the final cleaning up for all levels
  this->Finalize();
}

bool irtkImageMultiFreeFormRegistration2::Read(char *buffer1, char *buffer2, int &level)
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

  if (ok == false) {
    return this->irtkImageRegistration2::Read(buffer1, buffer2, level);
  } else {
    return ok;
  }
}

void irtkImageMultiFreeFormRegistration2::Write(ostream &to)
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
  if (_MFFDMode == true) {
      to << "MFFDMode                       = True" << endl;
  } else {
      to << "MFFDMode                       = False" << endl;
  }

  this->irtkImageRegistration2::Write(to);
}
