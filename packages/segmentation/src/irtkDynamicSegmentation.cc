/*=========================================================================

  Library   : Image Registration Toolkit (IRTK)
  Module    : $Id: irtkDynamicSegmentation.cc 101 2009-11-04 16:27:54Z dr $
  Copyright : Imperial College, Department of Computing
              Visual Information Processing (VIP), 2008 onwards
  Date      : $Date: 2009-11-04 16:27:54 +0000 (‰∏? 04 ÂçÅ‰∏ÄÊú?2009) $
  Version   : $Revision: 101 $
  Changes   : $Author: dr $

=========================================================================*/

#define GATLAS 2.0

#include <irtkSegmentationFunction.h>

template <class VoxelType> irtkDynamicSegmentation<VoxelType>::irtkDynamicSegmentation()
{
  // Set in- and outputs
  _input  = NULL;
  _output = NULL;
  /// Input transformation
  _mffd = NULL;

  /// Gaussian distributions
  _gaussian = NULL;

  /// Gaussian normal
  _normal = NULL;

  /// motion maginitude
  _motion = NULL;

  /// jacobian determine
  _jacobian = NULL;

  /// pointer
  _multiptr = NULL;

  /// atlas
  _atlas = NULL;

  /// graphcut
  _graphcut = NULL;

  /// NumberofLabels
  _NumberofLabels = 0;

  _mode = 1;
}

template <class VoxelType> irtkDynamicSegmentation<VoxelType>::~irtkDynamicSegmentation()
{
    //check input output
    if(_input == NULL || _atlas == NULL 
        || _mffd == NULL){
    }else{
        // Do the final cleaning up
        this->Finalize();
    }
    // Set in- and outputs
    _input  = NULL;
    _output = NULL;
    /// Input transformation
    _mffd = NULL;

    /// 3D Gaussian image
    _gaussian = NULL;

    /// Gaussian normal
    _normal = NULL;

    /// motion maginitude
    _motion = NULL;

    /// jacobian determine
    _jacobian = NULL;

    /// pointer
    _multiptr = NULL;

    /// atlas
    _atlas = NULL;

    /// graphcut
    _graphcut = NULL;

    /// NumberofLabels
    _NumberofLabels = 0;
}

template <class VoxelType> void irtkDynamicSegmentation<VoxelType>::SetInput(irtkGenericImage<VoxelType> *input,irtkMultiLevelFreeFormTransformation *mffd, irtkRealImage **atlas, int numberoflabels)
{
    int i,j;
    //check input output
    if(_input == NULL || _atlas == NULL 
        || _mffd == NULL){
    }else{
        // Do the final cleaning up
        this->Finalize();
    }

    if (input != NULL) {
        _input = input;
    } else {
        cerr << "irtkDynamicSegmentation::SetInput: Input is not an image\n";
        exit(1);
    }
    if (atlas != NULL) {
        _atlas = atlas;
    } else {
        cerr << "irtkDynamicSegmentation::SetInput: Segmentation is not an image\n";
        exit(1);
    }

    _NumberofLabels = numberoflabels+1;

    //padding
    irtkRealPixel **aptr = new irtkRealPixel*[_NumberofLabels];
    _atlas = new irtkRealImage*[_NumberofLabels];
    _atlas[0] = new irtkRealImage(atlas[0]->GetImageAttributes());
    aptr[0] = _atlas[0]->GetPointerToVoxels();
    for (i = 1; i < _NumberofLabels; i++){
        _atlas[i] = new irtkRealImage(*atlas[i-1]);
        aptr[i] = _atlas[i]->GetPointerToVoxels();
    }
    for(j = 0; j < _input->GetNumberOfVoxels(); j ++){
        *aptr[0] = 255.0;
        for (i = 1; i < _NumberofLabels; i++){
            *aptr[0] -= *aptr[i];
            aptr[i]++;
        }
        if(*aptr[0] < 0) *aptr[0] = 0;
        aptr[0]++;
    }
    delete []aptr;

    for (i = 0; i < _NumberofLabels - 1; i++){
        if(_atlas[i]->GetX() != _input->GetX()
            || _atlas[i]->GetY() != _input->GetY()
            || _atlas[i]->GetZ() != _input->GetZ()){
                cerr << "irtkDynamicSegmentation::SetInput: Atlas size does not equal to Input size\n";
                exit(1);
        }
    }

    if (mffd != NULL &&
        strcmp(mffd->NameOfClass(), "irtkMultiLevelFreeFormTransformation") == 0) {
        _mffd = mffd;
    } else {
        cerr << "irtkDynamicSegmentation::SetInput: Transformation is not irtkMultiLevelFreeFormTransformation\n";
        exit(1);
    }

    // Do the initial set up
    this->Initialize();
}



template <class VoxelType> void irtkDynamicSegmentation<VoxelType>::SetOutput(irtkGenericImage<VoxelType> *image)
{
  if (image != NULL) {
    _output = image;
  } else {
    cerr << "irtkDynamicSegmentation::SetOutput: Output is not an image\n";
    exit(1);
  }

  if(_output->GetX() != _input->GetX()
      || _output->GetY() != _input->GetY()
      || _output->GetZ() != _input->GetZ()){
          cerr << "irtkDynamicSegmentation::SetInput: Output size does not equal to Input size\n";
          exit(1);
  }
}

template <class VoxelType> void irtkDynamicSegmentation<VoxelType>::UpdateMFFD(int t)
{
    int i,j,k;
    double p1[3],p2[3],jac;
    //build motion
    for (k = 0; k < _motion->GetZ(); k++) {
        for (j = 0; j < _motion->GetY(); j++) {
            for (i = 0; i < _motion->GetX(); i++) {
                p1[0] = i; p1[1] = j; p1[2] = k;
                _motion->ImageToWorld(p1[0],p1[1],p1[2]);
                p2[0] = p1[0];
                p2[1] = p1[1];
                p2[2] = p1[2];
                _mffd->Transform(p2[0],p2[1],p2[2]);
                p2[0] -= p1[0];
                p2[1] -= p1[1];
                p2[2] -= p1[2];
                /*if(sqrt(p2[0]*p2[0] + p2[1]*p2[1] + p2[2]*p2[2])*10.0 > _motion->GetAsDouble(i,j,k)){
                    _motion->PutAsDouble(i,j,k,sqrt(p2[0]*p2[0] + p2[1]*p2[1] + p2[2]*p2[2])*10.0);
                }*/
                _motion->PutAsDouble(i,j,k,(_motion->GetAsDouble(i,j,k)*double(t)
                    + sqrt(p2[0]*p2[0] + p2[1]*p2[1] + p2[2]*p2[2])*10.0)/double(t+1));
            }
        }
    }
    _motion->Write("motion.nii.gz");
    //build jacobian
    for (k = 0; k < _jacobian->GetZ(); k++) {
        for (j = 0; j < _jacobian->GetY(); j++) {
            for (i = 0; i < _jacobian->GetX(); i++) {
                p1[0] = i; p1[1] = j; p1[2] = k;
                _jacobian->ImageToWorld(p1[0],p1[1],p1[2]);
                jac = _mffd->irtkTransformation::Jacobian(p1[0],p1[1],p1[2]);
                if (jac < 0.1) jac = 0.1;
                /*if(fabs(log(jac)*100) > fabs(_jacobian->GetAsDouble(i,j,k))){
                    _jacobian->PutAsDouble(i,j,k,log(jac)*100);
                }*/
                _jacobian->PutAsDouble(i,j,k,(_jacobian->GetAsDouble(i,j,k)*double(t)
                    + log(jac)*100)/double(t+1));
            }
        }
    }
    _jacobian->Write("jacobian.nii.gz");
}

template <class VoxelType> void irtkDynamicSegmentation<VoxelType>::Initialize()
{
    int i,j,k;
    double p1[3],p2[3],jac;

    //build motion
    _motion = new irtkGenericImage<VoxelType>(_input->GetImageAttributes());
    for (k = 0; k < _motion->GetZ(); k++) {
        for (j = 0; j < _motion->GetY(); j++) {
            for (i = 0; i < _motion->GetX(); i++) {
                p1[0] = i; p1[1] = j; p1[2] = k;
                _motion->ImageToWorld(p1[0],p1[1],p1[2]);
                p2[0] = p1[0];
                p2[1] = p1[1];
                p2[2] = p1[2];
                _mffd->Transform(p2[0],p2[1],p2[2]);
                p2[0] -= p1[0];
                p2[1] -= p1[1];
                p2[2] -= p1[2];
                _motion->PutAsDouble(i,j,k,sqrt(p2[0]*p2[0] + p2[1]*p2[1] + p2[2]*p2[2])*10);
            }
        }
    }
    _motion->Write("motion.nii.gz");
    //build jacobian
    _jacobian = new irtkGenericImage<VoxelType>(_input->GetImageAttributes());
    for (k = 0; k < _jacobian->GetZ(); k++) {
        for (j = 0; j < _jacobian->GetY(); j++) {
            for (i = 0; i < _jacobian->GetX(); i++) {
                p1[0] = i; p1[1] = j; p1[2] = k;
                _jacobian->ImageToWorld(p1[0],p1[1],p1[2]);
                jac = _mffd->irtkTransformation::Jacobian(p1[0],p1[1],p1[2]);
                if (jac < 0.1) jac = 0.1;
                _jacobian->PutAsDouble(i,j,k,log(jac)*100);
            }
        }
    }
    _jacobian->Write("jacobian.nii.gz");
    //build pointer
    _multiptr = new irtkGenericImage<VoxelType>*[3];
    _multiptr[0] = _input;
    _multiptr[1] = _motion;
    _multiptr[2] = _jacobian;

    //build probability
    _probability.AddProbabilityMaps(_NumberofLabels,_atlas);
    
    //graphcut
    _graphcut = new irtkImageGraphCut<VoxelType>;
    _graphcut->SetMode(3);
}

template <class VoxelType> void irtkDynamicSegmentation<VoxelType>::Finalize()
{
    int i;
    delete _motion;
    delete _jacobian;
    delete []_multiptr;
    for (i = 0; i < 3; i++){
        delete []_gaussian[i];
        delete []_normal[i];
    }
    for (i = 0; i < _NumberofLabels; i++){
        delete _atlas[i];
    }
    delete []_atlas;
    delete []_normal;
    delete []_gaussian;
    delete _graphcut;
    _NumberofLabels = 0;
}

template <class VoxelType> void irtkDynamicSegmentation<VoxelType>::Run()
{

  //check input output
  if(_input == NULL || _atlas == NULL 
      || _output == NULL || _mffd == NULL){
      cerr << "input or output not set yet" << endl;
      exit(1);
  }

  int k,j,i;
  //build Gaussian
  _gaussian = new irtkGaussian*[3];
  _normal = new double*[3];
  double mean, sdm;
  for (k = 0; k < 3; k++){
      _gaussian[k] = new irtkGaussian[_NumberofLabels];
      _normal[k] = new double[_NumberofLabels];
      VoxelType *ptr;
      switch(k){
      case 0:
          ptr = _input->GetPointerToVoxels();
          break;
      case 1:
          ptr = _motion->GetPointerToVoxels();
          break;
      case 2:
          ptr = _jacobian->GetPointerToVoxels();
          break;
      default:
          ptr = NULL;
          break;
      }
      mean = 0; sdm = 0;
      //evaluate means
      for(j = 0; j < _input->GetNumberOfVoxels(); j ++){
          mean += *ptr;
          ptr++;
      }  
      mean = mean/_input->GetNumberOfVoxels();
      //evaluate stds
      switch(k){
      case 0:
          ptr = _input->GetPointerToVoxels();
          break;
      case 1:
          ptr = _motion->GetPointerToVoxels();
          break;
      case 2:
          ptr = _jacobian->GetPointerToVoxels();
          break;
      default:
          break;
      }
      for(j = 0; j < _input->GetNumberOfVoxels(); j ++){
          sdm += pow((*ptr - mean),2);
          ptr++;
      }
      sdm = sdm/_input->GetNumberOfVoxels();
      if(sdm > 0){
          //em segmentation
          irtkEMClassification *classification;
          classification= new irtkEMClassification(_NumberofLabels, _atlas);

          if(k == 0){
              classification->SetPadding(0);
              classification->SetInput(*_input);
          }else if(k == 1){
              classification->SetPadding(-1);
              classification->SetInput(*_motion);
          }else if(k == 2){
              VoxelType min,max;
              _jacobian->GetMinMax(&min,&max);
              classification->SetPadding(min-1);
              classification->SetInput(*_jacobian);
          }
          classification->Initialise();

          double rel_diff,previous;
          double treshold = 0.0001;
          int iterations = 20;
          cout << "Segmenting "<<_NumberofLabels<<" tissues"<<endl;
          i=0; rel_diff = 100;
          do {
              previous = rel_diff;
              cout << "Iteration = " << i+1 << " / " << iterations << endl;
              rel_diff = classification->Iterate(i);
              i++;
          } while ((rel_diff>treshold)&&(i<iterations)&&(rel_diff < previous));

          double *cmean,*cvariance;
          cmean = new double[_NumberofLabels];
          cvariance = new double[_NumberofLabels];
          classification->GetMean(cmean);
          classification->GetVariance(cvariance);

          //weights for graph cut
          for(j = 0; j < _NumberofLabels; j ++) {
              _gaussian[k][j].Initialise(cmean[j],cvariance[j]*cvariance[j]);
              _normal[k][j] = 1.0/_gaussian[k][j].Getnorm();
          }
          delete []cvariance;
          delete []cmean;
          delete classification;
      }else{
          for(j = 0; j < _NumberofLabels; j++){
              _normal[k][j] = 0;
              _gaussian[k][j].Initialise(1,1);
          }
          _mode = 0;
      }
  }

  // Evaluate probability
  VoxelType *ptr_input,*ptr_motion,*ptr_jacobian;
  ptr_input = _input->GetPointerToVoxels();
  ptr_motion = _motion->GetPointerToVoxels();
  ptr_jacobian = _jacobian->GetPointerToVoxels();
  double *value = new double[_NumberofLabels];
  double norm;
  _probability.First();
  for(int j = 0; j < _input->GetNumberOfVoxels(); j ++){
      norm = 0;
      for(int i = 0; i < _NumberofLabels; i++){
          value[i] = 0;
          if(_normal[0][i] > 0)
              value[i] += _gaussian[0][i].Evaluate(*ptr_input);
          if(_normal[1][i] > 0)
              value[i] += _gaussian[1][i].Evaluate(*ptr_motion);
          if(_normal[2][i] > 0)
              value[i] += _gaussian[2][i].Evaluate(*ptr_jacobian);
          value[i] *= _probability.GetValue(i);
          norm+=value[i];
      }
      if(norm == 0){
          value[0] = 1.0;
      }else{
          for(int i = 0; i < _NumberofLabels; i++){
              value[i] = value[i] / norm;
              _probability.SetValue(i,value[i]);
          }
      }
      _probability.Next();
      ptr_input++;
      ptr_motion++;
      ptr_jacobian++;
  }
  delete []value;

  // Get probability
  irtkRealImage **test = new irtkRealImage*[_NumberofLabels];
  for(int n=0;n<_NumberofLabels;n++) {
      test[n] = new irtkRealImage(_input->GetImageAttributes());
      *test[n] = _probability.GetImage(n);
  }

  // Do graphcut
  if(_mode == 1){
      _graphcut->SetInput(3,_multiptr,_NumberofLabels,test);
  }else{
      _graphcut->SetInput(_input,_NumberofLabels,test);
  }
  _graphcut->SetOutput(_output);
  _graphcut->Run(1,1);
  for(int n=0;n<_NumberofLabels;n++) {
    delete test[n];
  }
  delete []test;
}

template class irtkDynamicSegmentation<unsigned char>;
template class irtkDynamicSegmentation<short>;
template class irtkDynamicSegmentation<unsigned short>;
template class irtkDynamicSegmentation<float>;
template class irtkDynamicSegmentation<double>;
