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

#define MAX_NO_LINE_ITERATIONS 12

irtkDiscontinuousFreeFormRegistration::irtkDiscontinuousFreeFormRegistration()
{
  // Print debugging information
  this->Debug("irtkDiscontinuousFreeFormRegistration::irtkDiscontinuousFreeFormRegistration");

  // Default optimization
  _OptimizationMethod = GradientDescent;

  // Default parameters for non-rigid registration
  _Lambda1     = 0;
  _Mode        = RegisterXYZ;
  _adjugate    = NULL;
  _determinant = NULL;
  _bssm = new irtkMultiLevelFreeFormTransformation;
  _bsfd = NULL;
  _mffd = NULL;
  _affd = NULL;
  _LevelOfStart = 0;
  _NumberOfModels = 0;
  _currentgradient = NULL;
  _previousFFD = 0;
  _compressrate = 0.5;
}

irtkDiscontinuousFreeFormRegistration::~irtkDiscontinuousFreeFormRegistration()
{
  delete _bssm;
}

void irtkDiscontinuousFreeFormRegistration::GuessParameter()
{
  int i;
  double xsize, ysize, zsize, spacing;
  irtkGreyPixel min,max;

  if ((_target == NULL) || (_source == NULL)) {
    cerr << "irtkDiscontinuousFreeFormRegistration::GuessParameter: Target and source image not found" << endl;
    exit(1);
  }

  // Default parameters for registration
  _NumberOfLevels     = 2;
  // Initial guess of bins with range
  _target->GetMinMax(&min,&max);
  _NumberOfBins       = round((max - min)/5.0);
  // Initial guess of bins with number of voxels
  if(_NumberOfBins > round(_target->GetNumberOfVoxels() / 1000.0)){
      _NumberOfBins = round(_target->GetNumberOfVoxels() / 1000.0);
  }
  // Don't be too small
  if(_NumberOfBins < 16){
      _NumberOfBins = 16;
  }

  // Default parameters for optimization
  _SimilarityMeasure  = NMI;
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
  _Lambda1            = 0; //recommended value 1

  // Remaining parameters
  for (i = 0; i < _NumberOfLevels; i++) {
    _NumberOfIterations[i] = 100;
    _MinStep[i]            = 0.01;
    _MaxStep[i]            = 4.0;
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

void irtkDiscontinuousFreeFormRegistration::Initialize()
{
    int i,j,k;
  // Print debugging information
  this->Debug("irtkDiscontinuousFreeFormRegistration::Initialize");

  // Initialize base class
  this->irtkImageRegistration2::Initialize();

  // Pointer to multi-level FFD
  _mffd = (irtkMultiLevelFreeFormTransformation *)_transformation;

  _LevelOfStart = _mffd->NumberOfLevels();

  // Intialize number of models
  double dx,dy,dz;
  dx = _target->GetXSize()*_target->GetX()/4;
  dy = _target->GetYSize()*_target->GetY()/4;
  if(_target->GetZ() > 1)
      dz = _target->GetZSize()*_target->GetZ()/4;
  else
      dz = 1;
  _NumberOfModels = 0;
  while(dx*2.0 > _target->GetXSize() && dy*2.0 > _target->GetYSize() 
      && (dz*2.0 > _target->GetZSize() || _target->GetZ() < 2)){
          _bsfd = new irtkBSplineFreeFormTransformation3D(*_target, dx, dy, dz);
          _bssm->PushLocalTransformation(_bsfd);

          _affd = new irtkBSplineFreeFormTransformation3D(*_target, dx, dy, dz);
          _mffd->PushLocalTransformation(_affd);
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

          dx /= 2; dy /= 2; 
          if(_target->GetZ() > 1) dz /= 2;
          _NumberOfModels++;
  }
  cout << "Number of models used: " << _NumberOfModels << endl;

  _adjugate = NULL;
  _determinant = NULL;
  _sparsityindex = new double[_NumberOfModels];
}

void irtkDiscontinuousFreeFormRegistration::Initialize(int level)
{

  // Print debugging information
  this->Debug("irtkDiscontinuousFreeFormRegistration::Initialize(int)");

  // Initialize base class
  this->irtkImageRegistration2::Initialize(level);

  // Set up _currentgradient residual
  irtkImageAttributes attr = _target->GetImageAttributes();

  attr._t = 3;
  _gradientResidual.Initialize(attr);

}

void irtkDiscontinuousFreeFormRegistration::Finalize()
{
  // Print debugging information
  this->Debug("irtkDiscontinuousFreeFormRegistration::Finalize");

  // Finalize base class
  this->irtkImageRegistration2::Finalize();

  delete []_sparsityindex;
}

void irtkDiscontinuousFreeFormRegistration::Finalize(int level)
{
  // Print debugging information
  this->Debug("irtkDiscontinuousFreeFormRegistration::Finalize(int)");

  // Finalize base class
  this->irtkImageRegistration2::Finalize(level);

  if(_currentgradient != NULL){
      delete []_currentgradient;
      _currentgradient = NULL;
  }
}

void irtkDiscontinuousFreeFormRegistration::UpdateVolume(bool affdchange)
{
    if(_Lambda1 > 0) {
        // Print debugging information
        this->Debug("irtkDiscontinuousFreeFormRegistration::UpdateVolume(bool affdchange)");
        if(affdchange){
            if(_adjugate != NULL) delete []_adjugate;
            if(_determinant != NULL) delete []_determinant;
            _adjugate = new irtkMatrix[_affd->NumberOfDOFs()/3];
            _determinant = new double[_affd->NumberOfDOFs()/3];
        }

        int index, index2, index3, i, j, k;
        double x,y,z,jacobian;
        // Update Jacobian and jacobian determinants;
        for (index = 0; index < _affd->NumberOfDOFs()/3; index++) {
            index2 = _affd->NumberOfDOFs()/3 + index;
            index3 = _affd->NumberOfDOFs()/3*2 + index;
            _affd->IndexToLattice(index,i,j,k);
            x = i;
            y = j;
            z = k;
            _affd->LatticeToWorld(x,y,z);
            irtkMatrix jac;
            jac.Initialize(3,3);
            _mffd->LocalJacobian(jac,x,y,z);
            jac.Adjugate(jacobian);
            if(jacobian < 0.0000001) jacobian = 0.0000001;
            _determinant[index] = jacobian;
            _adjugate[index] = jac;
        }
    }
}

void irtkDiscontinuousFreeFormRegistration::Update(bool updateGradient)
{
  // Print debugging information
  this->Debug("irtkDiscontinuousFreeFormRegistration::Update()");

  // Finalize base class
  this->irtkImageRegistration2::Update(updateGradient);

}

double irtkDiscontinuousFreeFormRegistration::VolumePreservationPenalty()
{
  int k;
  double penalty, jacobian;

  penalty = 0;
  for (k = 0; k < _affd->NumberOfDOFs()/3; k++) {
    // Determinant of Jacobian of deformation derivatives
    jacobian = _determinant[k];
    penalty += pow(log(jacobian), 2.0);
  }

  // Normalize sum by number of DOFs
  return -penalty / (double) _affd->NumberOfDOFs();
}

void irtkDiscontinuousFreeFormRegistration::VolumePreservationPenaltyGradient()
{
  int i, j, k, l, m, n, o, i1, j1, k1, i2, j2, k2, x, y, z, count, index, index1, index2, index3;
  double jacobian, drv[3];
  irtkMatrix jac,det_drv[3];

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


                  _affd->IndexToLattice(index, i, j, k);
                  l = i;
                  m = j;
                  n = k;
                  count = 0;
                  k1 = (k-1)>0?(k-1):0;
                  j1 = (j-1)>0?(j-1):0;
                  i1 = (i-1)>0?(i-1):0;
                  k2 = (k+2) < _affd->GetZ()? (k+2) : _affd->GetZ();
                  j2 = (j+2) < _affd->GetY()? (j+2) : _affd->GetY();
                  i2 = (i+2) < _affd->GetX()? (i+2) : _affd->GetX();

                  for(i=0;i<3;i++)
                      drv[i] = 0;

                  for (k = k1; k < k2; k++) {
                      for (j = j1; j < j2; j++) {
                          for (i = i1; i < i2; i++) {
                              if(k != n || j != m || i != l) {
                                  index1 = _affd->LatticeToIndex(i,j,k);
                                  // Torsten Rohlfing et al. MICCAI'01 (w/o scaling correction):
                                  // Calculate jacobian
                                  irtkMatrix jac = _adjugate[index1];

                                  // find jacobian derivatives
                                  jacobian = _determinant[index1];

                                  // if jacobian < 0
                                  jacobian = (2.0*log(jacobian))/jacobian;

                                  _affd->JacobianDetDerivative(det_drv,i-l,j-m,k-n);

                                  double tmpdrv[3];

                                  for(o = 0; o < 3; o++) {
                                      // trace * adjugate * derivative
                                      tmpdrv[o] = jac(0,o)*det_drv[o](o,0) + jac(1,o)*det_drv[o](o,1) + jac(2,o)*det_drv[o](o,2);
                                      // * rest of the volume preservation derivative
                                      drv[o] += (jacobian*tmpdrv[o]);
                                  }
                                  count ++;
                              }
                          }
                      }
                  }

                  for(l=0;l<3;l++)
                      drv[l] = -drv[l]/count;

                  _currentgradient[index]  += this->_Lambda1  * drv[0];
                  _currentgradient[index2] += this->_Lambda1  * drv[1];
                  _currentgradient[index3] += this->_Lambda1  * drv[2];
              }
          }
      }
  }

  return;
}

double irtkDiscontinuousFreeFormRegistration::Evaluate()
{
  double tmp, similarity;

  // Evaluate similarity
  similarity = this->irtkImageRegistration2::Evaluate();

  // Add penalty for volume preservation
  if (this->_Lambda1 > 0) {
    tmp = this->_Lambda1*this->VolumePreservationPenalty();
    cout << "\t Volume = " << tmp;
    similarity += tmp;
    cout << endl;
  }

  //Return similarity measure + penalty terms
  return similarity;
}

void irtkDiscontinuousFreeFormRegistration::ChooseLevel()
{
   int i,index;
   double maxsparsity;
       double temp = 0;
       maxsparsity = 0;
       index = _NumberOfModels - 1;
       for(i = 0; i < _NumberOfModels; i++){
           maxsparsity += _sparsityindex[i];
       }
       maxsparsity = maxsparsity*_compressrate;
       for(i = 0; i < _NumberOfModels; i++){
           temp += _sparsityindex[i];
           if(temp > maxsparsity){
               index = i;
               break;
           }
       }
       if(index+_LevelOfStart > _currentffdlevel){
           _currentffdlevel = index+_LevelOfStart;
       }else{
           maxsparsity = 0;
           for(i = 0; i < _NumberOfModels; i++){
               maxsparsity += _sparsityindex[i];
           }
           index = _currentffdlevel - _LevelOfStart;
           temp = 0;
           for(i = 0; i < index + 1; i++){
               temp += _sparsityindex[i];
           }
           _compressrate = temp / maxsparsity;
       }
       _affd = (irtkBSplineFreeFormTransformation *)_mffd->GetLocalTransformation(_currentffdlevel);
       cout << "Chosen level based on sparsity is: " << _currentffdlevel << " compress rate is: " << 1 - _compressrate << endl;
}

void irtkDiscontinuousFreeFormRegistration::DomainDecomposition()
{
    this->Debug("irtkDiscontinuousFreeFormRegistration::DomainDecomposition()");

    int i,j,k,t,index,index1,index2,index3;
    double max_x,max_y,max_z,cthresholdx,cthresholdy,cthresholdz;
    // x y z position
    double *x_pos;
    double *y_pos;
    double *z_pos;

    // x y z displacement
    double *x_dis;
    double *y_dis;
    double *z_dis;

    double energyl1[3];

    //debug
    if (_DebugFlag == true)
      _similarityGradient.Write("similaritygradient1.nii.gz");

    max_x = 0;
    max_y = 0;
    max_z = 0;

    // put similarity to residual
    for(i = 0; i < _similarityGradient.GetX(); i++){
        for(j = 0; j < _similarityGradient.GetY(); j++){
            for(k = 0; k < _similarityGradient.GetZ(); k++){
                _gradientResidual(i, j, k, 0) = _similarityGradient(i,j,k,0);
                if(fabs(_similarityGradient(i,j,k,0))>max_x) 
                    max_x = fabs(_similarityGradient(i,j,k,0));
                _gradientResidual(i, j, k, 1) = _similarityGradient(i,j,k,1);
                if(fabs(_similarityGradient(i,j,k,1))>max_y) 
                    max_y = fabs(_similarityGradient(i,j,k,1));
                _gradientResidual(i, j, k, 2) = _similarityGradient(i,j,k,2);
                if(fabs(_similarityGradient(i,j,k,2))>max_z) 
                    max_z = fabs(_similarityGradient(i,j,k,2));
            }
        }
    }

    irtkHistogram_1D<double> xhistogram(1000);
    irtkHistogram_1D<double> yhistogram(1000);
    irtkHistogram_1D<double> zhistogram(1000);

    xhistogram.PutMin(0);
    yhistogram.PutMin(0);
    zhistogram.PutMin(0);

    xhistogram.PutMax(max_x+0.0001);
    yhistogram.PutMax(max_y+0.0001);
    zhistogram.PutMax(max_z+0.0001);

    for (t = 0; t < _NumberOfModels; t++){

        // generate displacement to be approximated
        x_pos = new double[_gradientResidual.GetNumberOfVoxels()/3];
        y_pos = new double[_gradientResidual.GetNumberOfVoxels()/3];
        z_pos = new double[_gradientResidual.GetNumberOfVoxels()/3];

        x_dis = new double[_gradientResidual.GetNumberOfVoxels()/3];
        y_dis = new double[_gradientResidual.GetNumberOfVoxels()/3];
        z_dis = new double[_gradientResidual.GetNumberOfVoxels()/3];

        index = 0;

        for(i = 0; i < _gradientResidual.GetX(); i++){
            for(j = 0; j < _gradientResidual.GetY(); j++){
                for(k = 0; k < _gradientResidual.GetZ(); k++){
                    x_pos[index] = i;
                    y_pos[index] = j;
                    z_pos[index] = k;
                    _gradientResidual.ImageToWorld(x_pos[index],y_pos[index],z_pos[index]);
                    x_dis[index] = _gradientResidual(i, j, k, 0);
                    y_dis[index] = _gradientResidual(i, j, k, 1);
                    if(_gradientResidual.GetZ()>1)
                        z_dis[index] = _gradientResidual(i, j, k, 2);
                    else
                        z_dis[index] = 0;
                    index ++;
                }
            }
        }

        // get b spline model
        _bsfd = (irtkBSplineFreeFormTransformation *)_bssm->GetLocalTransformation(t);

        // approximate using finest b spline model.
        _bsfd->ApproximateAsNew(x_pos,y_pos,z_pos,x_dis,y_dis,z_dis,_gradientResidual.GetNumberOfVoxels()/3);

        index = 0;

        // put residual back to _currentgradient residual
        for(i = 0; i < _gradientResidual.GetX(); i++){
            for(j = 0; j < _gradientResidual.GetY(); j++){
                for(k = 0; k < _gradientResidual.GetZ(); k++){
                    _gradientResidual(i, j, k, 0) = x_dis[index];
                    _gradientResidual(i, j, k, 1) = y_dis[index];
                    if(_gradientResidual.GetZ() > 1){
                        _gradientResidual(i, j, k, 2) = z_dis[index];
                    }
                    index ++;
                }
            }
        }

        //calculate energy
        irtkRealImage debug(_bsfd->GetX(),_bsfd->GetY(),_bsfd->GetZ(),3);
        energyl1[0] = 0;
        energyl1[1] = 0;
        energyl1[2] = 0;
        index = 0;
        for(i = 0; i < _bsfd->GetX(); i++){
            for(j = 0; j < _bsfd->GetY(); j++){
                for(k = 0; k < _bsfd->GetZ(); k++){
                    index1  = _bsfd->LatticeToIndex(i,j,k);
                    index2 = index1+_bsfd->GetX()*_bsfd->GetY()*_bsfd->GetZ();
                    if (_DebugFlag == true){
                        debug(i, j, k, 0) = _bsfd->Get(index1);
                        debug(i, j, k, 1) = _bsfd->Get(index2); 
                    }
                    energyl1[0] += _bsfd->Get(index1);
                    xhistogram.AddSample(fabs(_bsfd->Get(index1)));
                    energyl1[1] += _bsfd->Get(index2);
                    yhistogram.AddSample(fabs(_bsfd->Get(index2)));
                    if(_bsfd->GetZ() > 1){
                        index3 = index1+2*_bsfd->GetX()*_bsfd->GetY()*_bsfd->GetZ();
                        if (_DebugFlag == true){
                            debug(i, j, k, 2) = _bsfd->Get(index3);
                        }
                        energyl1[2] += _bsfd->Get(index3);
                        zhistogram.AddSample(fabs(_bsfd->Get(index3)));
                    }                  
                    index ++;
                }
            }
        }

        for(i = 0; i < 3; i++){
            energyl1[i] = energyl1[i]/index;
        }

        cout << "level " << t << " L1 norm energy = " << energyl1[0] << " " << energyl1[1] << " " << energyl1[2] 
        << " sum = " << fabs(energyl1[0])+fabs(energyl1[1])+fabs(energyl1[2]) << endl;

        _sparsityindex[t] = fabs(energyl1[0])+fabs(energyl1[1])+fabs(energyl1[2]);

        // debug information
        if (_DebugFlag == true){
            char buffer[255];
            sprintf(buffer,"level%d.nii.gz",t);
            debug.Write(buffer);
        }

        delete []x_pos;
        delete []y_pos;
        delete []z_pos;
        delete []x_dis;
        delete []y_dis;
        delete []z_dis;
    }
    //debug information
    
    if (_DebugFlag == true){
        char buffer[255];
        sprintf(buffer,"level%d.nii.gz",t);
        _gradientResidual.Write(buffer);
        _bssm->irtkTransformation::Write("model.dof.gz");
    }

    double ratio = 0.99;
    cthresholdx = xhistogram.CDFToVal(ratio);
    cthresholdy = yhistogram.CDFToVal(ratio);
    if(_gradientResidual.GetZ() > 1)
        cthresholdz = zhistogram.CDFToVal(ratio);
    else
        cthresholdz = 0;

    // compress using the coefficients and project back to similarity gradient
    for(i = 0; i < _similarityGradient.GetX(); i++){
        for(j = 0; j < _similarityGradient.GetY(); j++){
            for(k = 0; k < _similarityGradient.GetZ(); k++){
                if(fabs(_gradientResidual(i,j,k,0)) < cthresholdx)
                    _gradientResidual(i,j,k,0) = 0;
                if(fabs(_gradientResidual(i,j,k,1)) < cthresholdy)
                    _gradientResidual(i,j,k,1) = 0;
                if(fabs(_gradientResidual(i,j,k,2)) < cthresholdz)
                    _gradientResidual(i,j,k,2) = 0;
            }
        }
    }
    // compress
    for (t = 0; t < _NumberOfModels; t++){
        _bsfd = (irtkBSplineFreeFormTransformation *)_bssm->GetLocalTransformation(t);
        for(i = 0; i < _bsfd->GetX(); i++){
            for(j = 0; j < _bsfd->GetY(); j++){
                for(k = 0; k < _bsfd->GetZ(); k++){
                    index1  = _bsfd->LatticeToIndex(i,j,k);
                    index2 = index1+_bsfd->GetX()*_bsfd->GetY()*_bsfd->GetZ();
                    if(fabs(_bsfd->Get(index1)) < cthresholdx)
                        _bsfd->Put(index1,0);
                    if(fabs(_bsfd->Get(index2)) < cthresholdy)
                        _bsfd->Put(index2,0);
                    if(_bsfd->GetZ() > 1){
                        index3 = index1+2*_bsfd->GetX()*_bsfd->GetY()*_bsfd->GetZ();
                        if(fabs(_bsfd->Get(index3)) < cthresholdz)
                            _bsfd->Put(index3,0);
                    }                  
                }
            }
        }
    }
    // project
    double p1[3],p2[3];
    for (k = 0; k < _similarityGradient.GetZ(); k++) {
        for (j = 0; j < _similarityGradient.GetY(); j++) {
            for (i = 0; i < _similarityGradient.GetX(); i++) {
                p1[0] = i;
                p1[1] = j;
                p1[2] = k;
                _similarityGradient.ImageToWorld(p1[0], p1[1], p1[2]);
                p2[0] = p1[0];
                p2[1] = p1[1];
                p2[2] = p1[2];
                _bssm->Transform(p2[0], p2[1], p2[2]);
                p2[0] -= p1[0];
                p2[1] -= p1[1];
                p2[2] -= p1[2];
                _similarityGradient(i, j, k, 0) = p2[0] + _gradientResidual(i, j, k, 0);
                _similarityGradient(i, j, k, 1) = p2[1] + _gradientResidual(i, j, k, 1);
                _similarityGradient(i, j, k, 2) = p2[2] + _gradientResidual(i, j, k, 2);
            }
        }
    }

    //debug
    if (_DebugFlag == true)
        _similarityGradient.Write("similaritygradient2.nii.gz");
}

void irtkDiscontinuousFreeFormRegistration::EvaluateGradient2D()
{
  double basis, pos[3];
  int i, j, i1, i2, j1, j2, k1, k2, x, y, index, index2, index3;

  // Initialize _currentgradient to zero
  for (i = 0; i < _affd->NumberOfDOFs(); i++) {
    _currentgradient[i] = 0;
  }

  // Loop over control points
  for (y = 0; y < _affd->GetY(); y++) {
    for (x = 0; x < _affd->GetX(); x++) {

      // Compute DoFs corresponding to the control point
      index  = _affd->LatticeToIndex(x, y, 0);
      index2 = index+_affd->GetX()*_affd->GetY()*_affd->GetZ();
      index3 = index+2*_affd->GetX()*_affd->GetY()*_affd->GetZ();

      // Check if any DoF corresponding to the control point is active
      if ((_affd->irtkTransformation::GetStatus(index) == _Active) || (_affd->irtkTransformation::GetStatus(index2) == _Active) || (_affd->irtkTransformation::GetStatus(index3) == _Active)) {

        // If so, calculate bounding box of control point in image coordinates
        _affd->BoundingBox(_target, index, i1, j1, k1, i2, j2, k2, 1.0);

        // Loop over all voxels in the target (reference) volume
        for (j = j1; j <= j2; j++) {
          for (i = i1; i <= i2; i++) {

            // Check whether reference point is valid
            if ((_target->Get(i, j, 0) >= 0) && (_transformedSource(i, j, 0) >= 0)) {

              // Convert position from voxel coordinates to world coordinates
              pos[0] = i;
              pos[1] = j;
              pos[2] = 0;
              _target->ImageToWorld(pos[0], pos[1], pos[2]);

              // Convert world coordinates into lattice coordinates
              _affd->WorldToLattice(pos[0], pos[1], pos[2]);

              // Compute B-spline tensor product at pos
              basis = _affd->B(pos[0] - x) * _affd->B(pos[1] - y);

              // Convert voxel-based _currentgradient into _currentgradient with respect to parameters (chain rule)
              //
              // NOTE: This currently assumes that the control points displacements are aligned with the world coordinate displacements
              //
              _currentgradient[index]  += basis * _similarityGradient(i, j, 0, 0);
              _currentgradient[index2] += basis * _similarityGradient(i, j, 0, 1);
              _currentgradient[index3] += 0;
            }

          }
        }
      }
    }
  }
}

void irtkDiscontinuousFreeFormRegistration::EvaluateGradient3D()
{
  double basis, pos[3];
  int i, j, k, i1, i2, j1, j2, k1, k2, x, y, z, index, index2, index3;

  // Initialize _currentgradient to zero
  for (i = 0; i < _affd->NumberOfDOFs(); i++) {
    _currentgradient[i] = 0;
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
          _affd->BoundingBox(_target, index, i1, j1, k1, i2, j2, k2, 1);

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

                  // Convert voxel-based _currentgradient into _currentgradient with respect to parameters (chain rule)
                  //
                  // NOTE: This currently assumes that the control points displacements are aligned with the world coordinate displacements
                  //
                  _currentgradient[index]  += basis * _similarityGradient(i, j, k, 0);
                  _currentgradient[index2] += basis * _similarityGradient(i, j, k, 1);
                  _currentgradient[index3] += basis * _similarityGradient(i, j, k, 2);
                }
              }
            }
          }
        }
      }
    }
  }
}

double irtkDiscontinuousFreeFormRegistration::EvaluateGradient(double *gradient)
{
  double norm, max_length;
  int i, x, y, z, index, index2, index3;
  static double *g = NULL, *h = NULL, gg, dgg, gamma;

  // Compute _currentgradient with respect to displacements
  this->irtkImageRegistration2::EvaluateGradient(gradient);

  if(*gradient > 0){
      this->DomainDecomposition();
  }
  this->ChooseLevel();
  this->UpdateVolume(true);

  // Allocate memory for _currentgradient vector
  if(_currentgradient != NULL){
    delete []_currentgradient;
  }
  _currentgradient = new double[_affd->NumberOfDOFs()];

  if (_affd->GetZ() == 1) {
    this->EvaluateGradient2D();
  } else {
    this->EvaluateGradient3D();
  }

  // Update _currentgradient
  if (_CurrentIteration == 0 || _previousFFD != _affd->NumberOfDOFs()) {
    // First iteration, so let's initialize
    _previousFFD = _affd->NumberOfDOFs();
    if (g != NULL) delete []g;
    g = new double [_previousFFD];
    if (h != NULL) delete []h;
    h = new double [_previousFFD];
    for (i = 0; i < _previousFFD; i++) {
      g[i] = -_currentgradient[i];
      h[i] = g[i];
    }
  } else {
    // Update _currentgradient direction to be conjugate
    gg = 0;
    dgg = 0;
    for (i = 0; i < _previousFFD; i++) {
      gg  += g[i]*h[i];
      dgg += (_currentgradient[i]+g[i])*_currentgradient[i];
    }
    gamma = dgg/gg;
    for (i = 0; i < _previousFFD; i++) {
      g[i] = -_currentgradient[i];
      h[i] = g[i] + gamma*h[i];
      _currentgradient[i] = -h[i];
    }
  }

  if (this->_Lambda1 > 0){
    this->VolumePreservationPenaltyGradient();
  }

  // Calculate maximum of _currentgradient vector
  max_length = 0;
  for (z = 0; z < _affd->GetZ(); z++) {
    for (y = 0; y < _affd->GetY(); y++) {
      for (x = 0; x < _affd->GetX(); x++) {
        index  = _affd->LatticeToIndex(x, y, z);
        index2 = index+_affd->GetX()*_affd->GetY()*_affd->GetZ();
        index3 = index+2*_affd->GetX()*_affd->GetY()*_affd->GetZ();
        norm = sqrt(_currentgradient[index] * _currentgradient[index] + _currentgradient[index2] * _currentgradient[index2] + _currentgradient[index3] * _currentgradient[index3]);
        if (norm > max_length) max_length = norm;
      }
    }
  }

  // Deal with active and passive control points
  for (i = 0; i < _affd->NumberOfDOFs(); i++) {
    if (_affd->irtkTransformation::GetStatus(i) == _Passive) {
      _currentgradient[i] = 0;
    }
  }

  return max_length;
}

void irtkDiscontinuousFreeFormRegistration::Run()
{
  int i, k;
  char buffer[256];
  double gradient, delta, step, min_step, max_step, max_length, best_similarity, new_similarity, old_similarity;

  // Print debugging information
  this->Debug("irtkDiscontinuousFreeFormRegistration::Run");

  if (_source == NULL) {
    cerr << "irtkDiscontinuousFreeFormRegistration::Run: Filter has no source input" << endl;
    exit(1);
  }

  if (_target == NULL) {
    cerr << "irtkDiscontinuousFreeFormRegistration::Run: Filter has no target input" << endl;
    exit(1);
  }

  if (_transformation == NULL) {
    cerr << "irtkDiscontinuousFreeFormRegistration::Run: Filter has no transformation output" << endl;
    exit(1);
  }

  // Do the initial set up for all levels
  this->Initialize();

  // Loop over levels
  for (_CurrentLevel = _NumberOfLevels-1; _CurrentLevel >= 0; _CurrentLevel--) {
    
    gradient = 1.0;

    // Initial step size
    min_step = _MinStep[_CurrentLevel];
    max_step = _MaxStep[_CurrentLevel];

    // Print resolution level
    cout << "Resolution level no. " << _CurrentLevel+1 << " (step sizes " << min_step << " to " << max_step  << ")\n";

    // Initialize for this level
    this->Initialize(_CurrentLevel);

    // Save pre-processed images if we are debugging
    sprintf(buffer, "source_%d.nii.gz", _CurrentLevel);
    if (_DebugFlag == true) _source->Write(buffer);
    sprintf(buffer, "target_%d.nii.gz", _CurrentLevel);
    if (_DebugFlag == true) _target->Write(buffer);

    // Run the registration filter at this resolution
    _CurrentIteration = 0;
    _compressrate = 0.6;
    _currentffdlevel = _LevelOfStart;
    while (_CurrentIteration < _NumberOfIterations[_CurrentLevel]) {
      cout << "Iteration = " << _CurrentIteration + 1 << " (out of " << _NumberOfIterations[_CurrentLevel] << ")"<< endl;

      // Update source image
      this->Update(true);

      // Compute current metric value
      best_similarity = old_similarity = this->Evaluate();
      cout << "Current best metric value is " << best_similarity << endl;

      // Compute _currentgradient of similarity metric. The function EvaluateGradient() returns the maximum control point length in the _currentgradient
      max_length = this->EvaluateGradient(&gradient);

      // Step along _currentgradient direction until no further improvement is necessary
      i = 0;
      delta = 0;
      step = max_step;
      do {
        double current = step / max_length;

        // Move along _currentgradient direction
        for (k = 0; k < _affd->NumberOfDOFs(); k++) {
          _affd->Put(k, _affd->Get(k) + current * _currentgradient[k]);
        }

        // We have just changed the transformation parameters, so we need to update
        this->Update(false);
        this->UpdateVolume(false);

        // Compute new similarity
        new_similarity = this->Evaluate();

        if (new_similarity > best_similarity + _Epsilon) {
          cout << "New metric value is " << new_similarity << "; step = " << step << endl;
          best_similarity = new_similarity;
          delta += step;
          step = step * 1.1;
          if (step > max_step) step = max_step;

        } else {
          // Last step was no improvement, so back track
          //cout << "Rejected metric value is " << new_similarity << "; step = " << step << endl;
          for (k = 0; k < _affd->NumberOfDOFs(); k++) {
            _affd->Put(k, _affd->Get(k) - current * _currentgradient[k]);
          }
          step = step * 0.5;
        }
        i++;
      } while ((i < MAX_NO_LINE_ITERATIONS) && (step > min_step));

      _CurrentIteration++;

      // Check for convergence
      if (delta == 0) {
          if(_currentffdlevel >= _mffd->NumberOfLevels() - 1){
              break;
          }else{
              _compressrate += 0.1;
              _currentffdlevel ++;
              gradient = -1.0;
          }
      }else{
          _compressrate = 0.6;
          _currentffdlevel = floor(_currentffdlevel/2.0);
          gradient = 1.0;
      }
    }

    // Do the final cleaning up for this level
    this->Finalize(_CurrentLevel);
  }

  // Do the final cleaning up for all levels
  this->Finalize();
}

bool irtkDiscontinuousFreeFormRegistration::Read(char *buffer1, char *buffer2, int &level)
{
  int ok = false;

  if ((strstr(buffer1, "Volume preserve") != NULL)) {
    this->_Lambda1 = atof(buffer2);
    cout << "Volume preservation penalty is ... " << this->_Lambda1 << endl;
    ok = true;
  }

  if (ok == false) {
    return this->irtkImageRegistration2::Read(buffer1, buffer2, level);
  } else {
    return ok;
  }
}

void irtkDiscontinuousFreeFormRegistration::Write(ostream &to)
{
  to << "\n#\n# Sparse non-rigid registration parameters\n#\n\n";
  to << "Volume preserve                     = " << this->_Lambda1 << endl;

  this->irtkImageRegistration2::Write(to);
}
