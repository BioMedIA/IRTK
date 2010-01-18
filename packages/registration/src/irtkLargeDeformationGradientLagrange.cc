/*=========================================================================

  Date      : $Date: 29.06.2009$
  Changes   : $Authors: Laurent Risser, Francois-Xavier Vialard$

=========================================================================*/

#include <irtkLargeDeformationGradientLagrange.h>

///++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
///                                       GENERAL CLASS FUNCTIONS
///++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
template <class VoxelType> LargeDefGradLagrange<VoxelType>::LargeDefGradLagrange(){
  //default parameters
  epsilon=0.1;
  iteration_nb=10;
  NbTimeSubdiv=10;
  MaxVelocityUpdate=2.;  //rem: Delta Voxels = 1
  sigmaX1=-1.; sigmaY1=-1.; sigmaZ1=-1.;
  sigmaX2=-1.; sigmaY2=-1.; sigmaZ2=-1.;
  sigmaX3=-1.; sigmaY3=-1.; sigmaZ3=-1.;
  sigmaX4=-1.; sigmaY4=-1.; sigmaZ4=-1.;
  NbKernels=1;
  Margin=0;
  alpha=0.01;
  gamma=0.1;
  WghtVelField=1.;
  reparametrization = 200;
  strcpy(PrefixInputVF,"Null");
  strcpy(PrefixOutputVF,"Null");
}

template <class VoxelType> LargeDefGradLagrange<VoxelType>::~LargeDefGradLagrange(void)
{}

template <class VoxelType> Bool LargeDefGradLagrange<VoxelType>::RequiresBuffering(void)
{
  return True;
}

template <class VoxelType> const char *LargeDefGradLagrange<VoxelType>::NameOfClass()
{
  return "LargeDefGradLagrange";
}


///++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
///                        SUB-FUNCTIONS TO PERFORM THE REGISTRATION  (1: main subfunctions)
///++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++



template <class VoxelType> void LargeDefGradLagrange<VoxelType>::HighPassFilter(float* Image, float SigmaFilter){   //added
  int x,y,z;   //added
  irtkGaussianBlurring<float> gaussianBlurring(SigmaFilter);   //added
     //added
  for (z = 0; z < this->NZ; z++) for (y = 0; y < this->NY; y++) for (x = 0; x < this->NX; x++)   //added
        Image3DTemp.Put(x, y, z,Image[ptSF(x,y,z)]);   //added
     //added
  gaussianBlurring.SetInput (&this->Image3DTemp);   //added
  gaussianBlurring.SetOutput(&this->Image3DTemp);   //added
  gaussianBlurring.Run();   //added
     //added
  for (z = 0; z < this->NZ; z++) for (y = 0; y < this->NY; y++) for (x = 0; x < this->NX; x++)   //added
        Image[ptSF(x,y,z)] -= Image3DTemp.Get(x, y, z);   //added
}   //added





///allocate all variables used for the gradient descent (Beg 2005) of the current 3D image from the treated 4D time sequence.
///Compute also the dimension of the scalar and vector fields in use
template <class VoxelType> void LargeDefGradLagrange<VoxelType>::AllocateAllVariables(void){
  int i;
  
  //time step between two subdivision
  this->DeltaTimeSubdiv=1./(static_cast<float>(NbTimeSubdiv-1));
  
  //local variables calculations (mostly to improve the readability)
  this->NX=this->_input->GetX();
  this->NY=this->_input->GetY();
  this->NZ=this->_input->GetZ();
  this->NT=this->_input->GetT();
  this->NXtY=this->NX*this->NY;
  this->NXtYtZ=this->NXtY*this->NZ;
  this->NXt3=this->NX*3;
  this->NXtYt3=this->NXt3*this->NY;
  this->NXtYtZt3=this->NXtYt3*this->NZ;
  
  cout << "Image size: " << this->NX <<  " , "  <<  this->NY  <<  " , "  << this->NZ  <<  " , " << this->NT  << "\n";
  
  //... temporary casts of the template and target images
  //    -->  ImTemplate[ptSF(x,y,z)]= gray level of the template image at (x,y,z)
  //    -->  ImTarget[ptSF(x,y,z)]= gray level of the target image at (x,y,z)
  this->ImTemplate = new float [this->NXtYtZ];
  this->ImTarget = new float [this->NXtYtZ];
  
  
  //... velocity field and VelocityFieldTemp (only for the time reparameterization)
  //    -->  VelocityField[i][ptVF(x,y,z,0)]= direction ex of the vector at (x,y,z)
  //    -->  VelocityField[i][ptVF(x,y,z,1)]= direction ey of the vector at (x,y,z)
  //    -->  VelocityField[i][ptVF(x,y,z,2)]= direction ez of the vector at (x,y,z)
  //    -->  where n is the id of the velocity field
  this->VelocityField= new float* [this->NbTimeSubdiv];
  for (i=0;i<this->NbTimeSubdiv;i++) this->VelocityField[i]=new float [this->NXtYtZt3];
  
  this->VelocityFieldTemp= new float* [this->NbTimeSubdiv];
  for (i=0;i<this->NbTimeSubdiv;i++) this->VelocityFieldTemp[i]=new float [this->NXtYtZt3];

  
  //... forward mapping
  //    -->  ForwardMapping[i][ptVF(x,y,z,0)]= coordinate x at time i corresponding to (x,y,z) at virtual time 0
  //    -->  ForwardMapping[i][ptVF(x,y,z,1)]= coordinate y at time i corresponding to (x,y,z) at virtual time 0
  //    -->  ForwardMapping[i][ptVF(x,y,z,2)]= coordinate z at time i corresponding to (x,y,z) at virtual time 0
  this->ForwardMapping= new float* [this->NbTimeSubdiv];
  for (i=0;i<this->NbTimeSubdiv;i++) this->ForwardMapping[i]=new float [this->NXtYtZt3];
  
  //... backward mapping
  //    -->  BackwardMapping[i][ptVF(x,y,z,0)]= coordinate x at time i corresponding to (x,y,z) at virtual time 1
  //    -->  BackwardMapping[i][ptVF(x,y,z,1)]= coordinate y at time i corresponding to (x,y,z) at virtual time 1
  //    -->  BackwardMapping[i][ptVF(x,y,z,2)]= coordinate z at time i corresponding to (x,y,z) at virtual time 1
  this->BackwardMapping= new float* [this->NbTimeSubdiv];
  for (i=0;i<this->NbTimeSubdiv;i++) this->BackwardMapping[i]=new float [this->NXtYtZt3];
  
  //... temporary image transformed using the forward mapping from time 0
  //    -->  J0[ptSF(x,y,z)]= gray level of the transformed image J0 at (x,y,z)
  this->J0 = new float [this->NXtYtZ];
  
  //... temporary image transformed using the backward mapping from time 1
  //    -->  J1[ptSF(x,y,z)]= gray level of the transformed image J1 at (x,y,z)
  this->J1 = new float [this->NXtYtZ];
  
  //... gradient of J0
  //    -->  GradJ0[ptVF(x,y,z,0)]= gradient of J0 in direction ex at (x,y,z)
  //    -->  GradJ0[ptVF(x,y,z,1)]= gradient of J0 in direction ey at (x,y,z)
  //    -->  GradJ0[ptVF(x,y,z,2)]= gradient of J0 in direction ez at (x,y,z)
  this->GradJ0 = new float [this->NXtYtZt3];
  
  //... determinent of the Jacobians  
  //    -->  Jacobians[ptSF(x,y,z)]= jacobian at (x,y,z)
  this->DetJacobians = new float [this->NXtYtZ];
  
  //... Energy Gradient
  //    -->  GradE[i][ptVF(x,y,z,0)]= Energy gradient at time i in direction ex at (x,y,z)
  //    -->  GradE[i][ptVF(x,y,z,1)]= Energy gradient at time i in direction ey at (x,y,z)
  //    -->  GradE[i][ptVF(x,y,z,2)]= Energy gradient at time i in direction ez at (x,y,z)
  this->GradE= new float* [this->NbTimeSubdiv];
  for (i=0;i<this->NbTimeSubdiv;i++)   this->GradE[i] = new float [this->NXtYtZt3];
  
  
  //Jacobian Matrix
  //     --> JacobianMatrix[3*i+j]: value of the matrix in (i,j)
  this->JacobianMatrix = new float[9];
  
  //Norm of the velocity field
  //     --> norm[i]: norm of the velicity field at subdivision time i
  this->norm = new float [this->NbTimeSubdiv];
  
  //temporary image to use when irtk functions are called
  this->Image3DTemp = irtkGenericImage<float>(this->NX, this->NY, this->NZ,1);
  
  //images to perform the FFT
  this->NXfft=(int)(pow(2.,floor((log((double)this->NX)/log(2.))+0.99999))+0.00001); //smaller size higher than 'this->NX' and being a power of 2
  this->NYfft=(int)(pow(2.,floor((log((double)this->NY)/log(2.))+0.99999))+0.00001); // ... 'this->NY' ...
  this->NZfft=(int)(pow(2.,floor((log((double)this->NZ)/log(2.))+0.99999))+0.00001); // ... 'this->NZ' ...
  this->SXfft=(NXfft-this->NX)/2;   // first X coordinate where the treated images are copied in the images for the fft
  this->SYfft=(NYfft-this->NY)/2;   // first X coordinate where the treated images are copied in the images for the fft
  this->SZfft=(NZfft-this->NZ)/2;   // first X coordinate where the treated images are copied in the images for the fft
  
  cout << "Images to perform FFTs: " << this->NXfft << " , " << this->NYfft  << " , " << this->NZfft  << " (" << this->SXfft  << "," <<  this->SYfft << "," << this->SZfft  << ")\n";
  
  this->RealSignalForFFT = irtkGenericImage<float>(this->NXfft, this->NYfft, this->NZfft,1); //image  - real part
  this->ImagSignalForFFT = irtkGenericImage<float>(this->NXfft, this->NYfft, this->NZfft,1); //image  - imaginary part
  this->RealFilterForFFT = irtkGenericImage<float>(this->NXfft, this->NYfft, this->NZfft,1); //filter - real part
  this->ImagFilterForFFT = irtkGenericImage<float>(this->NXfft, this->NYfft, this->NZfft,1); //filter - imaginary part

  
}


///initiate the gradient descent (Beg 2005) for the current 3D image of the 4D time sequence
template <class VoxelType> void LargeDefGradLagrange<VoxelType>::InitiateGradientDescent(int TimeLoc){
  int subivId, x, y, z;
  int DistClosestEdge;

  //1) cast the values of the template and target images at time t in float
  for (z = 0; z < this->NZ; z++)  for (y = 0; y < this->NY; y++) for (x = 0; x < this->NX; x++)
        this->ImTemplate[ptSF(x,y,z)]=static_cast<float>(this->_input->Get(x, y, z, TimeLoc));
  
  for (z = 0; z < this->NZ; z++)  for (y = 0; y < this->NY; y++) for (x = 0; x < this->NX; x++)
        this->ImTarget[ptSF(x,y,z)]=static_cast<float>(this->target_image.Get(x, y, z, TimeLoc));
  
  //2) take into accout the margin
    for (z = 0; z < this->NZ; z++)  for (y = 0; y < this->NY; y++) for (x = 0; x < this->NX; x++){
    DistClosestEdge=z+1;
    if (y+1<DistClosestEdge) DistClosestEdge=y+1;
    if (x+1<DistClosestEdge) DistClosestEdge=x+1;
    if (this->NZ-z<DistClosestEdge) DistClosestEdge=this->NZ-z;
    if (this->NY-y<DistClosestEdge) DistClosestEdge=this->NY-y;
    if (this->NX-x<DistClosestEdge) DistClosestEdge=this->NX-x;
    if (DistClosestEdge<=this->Margin){
      this->ImTemplate[ptSF(x,y,z)]=this->ImTemplate[ptSF(x,y,z)]*((DistClosestEdge-1.)/this->Margin)*((DistClosestEdge-1.)/this->Margin);
      this->ImTarget[ptSF(x,y,z)]=this->ImTarget[ptSF(x,y,z)]*((DistClosestEdge-1.)/this->Margin)*((DistClosestEdge-1.)/this->Margin);
    }
  }
  
  //3) Gray level regularisation  (the gray levels of the target and the template are sampled between 0 and 100)
  this->MinGrayLevelImTemplate=this->ImTemplate[ptSF(1,1,1)];
  for (z = 0; z < this->NZ; z++)  for (y = 0; y < this->NY; y++) for (x = 0; x < this->NX; x++)
        if (this->ImTemplate[ptSF(x,y,z)]<this->MinGrayLevelImTemplate) this->MinGrayLevelImTemplate=this->ImTemplate[ptSF(x,y,z)];
  
  this->MaxGrayLevelImTemplate=this->ImTemplate[ptSF(1,1,1)];
  for (z = 0; z < this->NZ; z++)  for (y = 0; y < this->NY; y++) for (x = 0; x < this->NX; x++)
        if (this->ImTemplate[ptSF(x,y,z)]>this->MaxGrayLevelImTemplate) this->MaxGrayLevelImTemplate=this->ImTemplate[ptSF(x,y,z)];

  for (z = 0; z < this->NZ; z++)  for (y = 0; y < this->NY; y++) for (x = 0; x < this->NX; x++)
        this->ImTemplate[ptSF(x,y,z)]=100.*(this->ImTemplate[ptSF(x,y,z)]-this->MinGrayLevelImTemplate)/(this->MaxGrayLevelImTemplate-this->MinGrayLevelImTemplate);
  
  this->MinGrayLevelImTarget=this->ImTarget[ptSF(1,1,1)];
  for (z = 0; z < this->NZ; z++)  for (y = 0; y < this->NY; y++) for (x = 0; x < this->NX; x++)
        if (this->ImTarget[ptSF(x,y,z)]<this->MinGrayLevelImTarget) this->MinGrayLevelImTarget=this->ImTarget[ptSF(x,y,z)];
  
  this->MaxGrayLevelImTarget=this->ImTarget[ptSF(1,1,1)];
  for (z = 0; z < this->NZ; z++)  for (y = 0; y < this->NY; y++) for (x = 0; x < this->NX; x++)
        if (this->ImTarget[ptSF(x,y,z)]>this->MaxGrayLevelImTarget) this->MaxGrayLevelImTarget=this->ImTarget[ptSF(x,y,z)];

  for (z = 0; z < this->NZ; z++)  for (y = 0; y < this->NY; y++) for (x = 0; x < this->NX; x++)
        this->ImTarget[ptSF(x,y,z)]=100.*(this->ImTarget[ptSF(x,y,z)]-this->MinGrayLevelImTarget)/(this->MaxGrayLevelImTarget-this->MinGrayLevelImTarget);

  //3.added) high pass filter  //added
  if (alpha<0.){  //added
    this->HighPassFilter(this->ImTemplate,5);  //added
    this->HighPassFilter(this->ImTarget,5);  //added
  }  //added
  
//   for (z = 0; z < this->NZ; z++) for (y = 0; y < this->NY; y++) for (x = 0; x < this->NX; x++){  //added
//     this->Image3DTemp.Put(x, y, z, 0, this->ImTemplate[ptSF(x,y,z)]);  //added
//   }  //added
//   this->Image3DTemp.Write("Template.nii");  //added
//     
//   for (z = 0; z < this->NZ; z++) for (y = 0; y < this->NY; y++) for (x = 0; x < this->NX; x++){  //added
//     this->Image3DTemp.Put(x, y, z, 0, this->ImTarget[ptSF(x,y,z)]);  //added
//   }  //added
//   this->Image3DTemp.Write("Target.nii");  //added

  
  
  //4) initiate the velocity field
  if (strcmp(PrefixInputVF,"Null")!=0){
    this->LoadVelocityFields(PrefixInputVF);
  }
  else{
    for (subivId=0;subivId<this->NbTimeSubdiv;subivId++) for (z=0;z<this->NZ;z++) for (y=0;y<this->NY;y++) for (x=0;x<this->NX;x++){
      this->VelocityField[subivId][ptVF(x,y,z,0)]=0;
      this->VelocityField[subivId][ptVF(x,y,z,1)]=0;
      this->VelocityField[subivId][ptVF(x,y,z,2)]=0;
    }
  }

  //5) initiate the values of the forward mapping at subdivision time 0
  for (z=0;z<this->NZ;z++) for (y=0;y<this->NY;y++) for (x=0;x<this->NX;x++){
    this->ForwardMapping[0][ptVF(x,y,z,0)]=0.;
    this->ForwardMapping[0][ptVF(x,y,z,1)]=0.;
    this->ForwardMapping[0][ptVF(x,y,z,2)]=0.;
  }

  //6) initiate the values of the backward mapping at subdivision time 1
  for (z=0;z<this->NZ;z++) for (y=0;y<this->NY;y++) for (x=0;x<this->NX;x++){
    this->BackwardMapping[this->NbTimeSubdiv-1][ptVF(x,y,z,0)]=0.;
    this->BackwardMapping[this->NbTimeSubdiv-1][ptVF(x,y,z,1)]=0.;
    this->BackwardMapping[this->NbTimeSubdiv-1][ptVF(x,y,z,2)]=0.;
  }

  //7) initiate GradE
  for (subivId=0;subivId<this->NbTimeSubdiv;subivId++) for (z=0;z<this->NZ;z++) for (y=0;y<this->NY;y++) for (x=0;x<this->NX;x++){
    this->GradE[subivId][ptVF(x,y,z,0)]=0;
    this->GradE[subivId][ptVF(x,y,z,1)]=0;
    this->GradE[subivId][ptVF(x,y,z,2)]=0;
  }

  //8)Initiate the images for the fft + Draw and transform in Fourier spaces the kernel of the filter
  //8.1) Initiate all images for fft at 0
  for (z=0;z<this->NZfft;z++) for (y=0;y<this->NYfft;y++) for (x=0;x<this->NXfft;x++) this->RealSignalForFFT.Put(x,y,z,0,0);
  for (z=0;z<this->NZfft;z++) for (y=0;y<this->NYfft;y++) for (x=0;x<this->NXfft;x++) this->ImagSignalForFFT.Put(x,y,z,0,0);
  for (z=0;z<this->NZfft;z++) for (y=0;y<this->NYfft;y++) for (x=0;x<this->NXfft;x++) this->RealFilterForFFT.Put(x,y,z,0,0);
  for (z=0;z<this->NZfft;z++) for (y=0;y<this->NYfft;y++) for (x=0;x<this->NXfft;x++) this->ImagFilterForFFT.Put(x,y,z,0,0);

  //8.2) Draw the filter
  //8.2.1) defaulf filter -original filter of Beg's paper 2005
  this->RealFilterForFFT.Put(0,0,0,0,(float)(2*this->alpha*3.+this->gamma));
  this->RealFilterForFFT.Put(1,0,0,0,(float)(this->alpha));
  this->RealFilterForFFT.Put(0,1,0,0,(float)(this->alpha));
  if (this->NZfft>1) this->RealFilterForFFT.Put(0,0,1,0,(float)(this->alpha));
  this->RealFilterForFFT.Put(this->NXfft-1,0,0,0,(float)(this->alpha));
  this->RealFilterForFFT.Put(0,this->NYfft-1,0,0,(float)(this->alpha));
  if (this->NZfft>1) this->RealFilterForFFT.Put(0,0,this->NZfft-1,0,(float)(this->alpha));

  
  //8.2.2) isotropic or anisotropic gaussian filter
  if ((this->sigmaX1>0.)&&(this->sigmaX2<0.))
    MakeAnisotropicGaussianFilter(this->weight1,this->sigmaX1,this->sigmaY1,this->sigmaZ1,&this->RealFilterForFFT,&this->ImagFilterForFFT);

  //8.2.3) sum of anisotropic filters
  if ((this->sigmaX2>0.)&&(this->sigmaX3<0.))
    MakeSumOf2AnisotropicGaussianFilters(this->weight1,this->sigmaX1,this->sigmaY1,this->sigmaZ1,this->weight2,this->sigmaX2,this->sigmaY2,this->sigmaZ2,&this->RealFilterForFFT,&this->ImagFilterForFFT);
  
  if ((this->sigmaX3>0.)&&(this->sigmaX4<0.))
    MakeSumOf3AnisotropicGaussianFilters(this->weight1,this->sigmaX1,this->sigmaY1,this->sigmaZ1,this->weight2,this->sigmaX2,this->sigmaY2,this->sigmaZ2,this->weight3,this->sigmaX3,this->sigmaY3,this->sigmaZ3,&this->RealFilterForFFT,&this->ImagFilterForFFT);
  
  if (this->sigmaX4>0.)
    MakeSumOf4AnisotropicGaussianFilters(this->weight1,this->sigmaX1,this->sigmaY1,this->sigmaZ1,this->weight2,this->sigmaX2,this->sigmaY2,this->sigmaZ2,this->weight3,this->sigmaX3,this->sigmaY3,this->sigmaZ3,this->weight4,this->sigmaX4,this->sigmaY4,this->sigmaZ4,&this->RealFilterForFFT,&this->ImagFilterForFFT);
  
  //this->RealFilterForFFT.Write("Filter.nii");

  //8.2.4) chain of anisotropic gaussian filters
  // --- NOT TREATED HERE: filters are defined in the main loop of the program ---
  
  //8.3) transformation of the filter in Fourier spaces (once for all)
  DirectFFT(&this->RealFilterForFFT,&this->ImagFilterForFFT);

  

  //9) initiate norm, Length and Energy
  for (subivId=0;subivId<this->NbTimeSubdiv;subivId++) this->norm[subivId]=0;
  this->length=0;
  this->EnergyTot=0;
  this->EnergyVF=0;
  this->EnergyDI=0;
}


///compute the forward mapping
template <class VoxelType> void LargeDefGradLagrange<VoxelType>::ComputeForwardMapping(void){
  float VecTemp[3];  //velocity vector that targets the current treated voxel
  float VecTemp2[3]; //temporary vector to propagate the mapping
  int ConvergenceSteps;
  int i,x,y,z;
  int timeSubdiv;
      
  //initialisation
  ConvergenceSteps=3;
  
  //JO at the first time subdivision
  for (z = 0; z < this->NZ; z++) for (y = 0; y < this->NY; y++) for (x = 0; x < this->NX; x++){
    this->ForwardMapping[0][ptVF(x,y,z,0)]=static_cast<float>(x);
    this->ForwardMapping[0][ptVF(x,y,z,1)]=static_cast<float>(y);
    this->ForwardMapping[0][ptVF(x,y,z,2)]=static_cast<float>(z);
  }
  
  //JO at the other time subdivisions
  for (timeSubdiv=1;timeSubdiv<this->NbTimeSubdiv;timeSubdiv++){
    for (z = 0; z < this->NZ; z++) for (y = 0; y < this->NY; y++) for (x = 0; x < this->NX; x++){
      //init
      VecTemp[0]=0.; 
      VecTemp[1]=0.;
      VecTemp[2]=0.;
      
      //convergence
      for (i=0;i<ConvergenceSteps;i++){
        this->GetVectorFromVelocityField(timeSubdiv-1,x-VecTemp[0],y-VecTemp[1],z-VecTemp[2],VecTemp);
        VecTemp[0]=(VecTemp[0]+this->VelocityField[timeSubdiv][ptVF(x,y,z,0)])*(this->DeltaTimeSubdiv)/2;
        VecTemp[1]=(VecTemp[1]+this->VelocityField[timeSubdiv][ptVF(x,y,z,1)])*(this->DeltaTimeSubdiv)/2;
        VecTemp[2]=(VecTemp[2]+this->VelocityField[timeSubdiv][ptVF(x,y,z,2)])*(this->DeltaTimeSubdiv)/2;
      }
      
      //find the original coordinates
      this->GetCoordFromForwardMapping(timeSubdiv-1,x-VecTemp[0],y-VecTemp[1],z-VecTemp[2],VecTemp2);
      
      this->ForwardMapping[timeSubdiv][ptVF(x,y,z,0)]=VecTemp2[0];
      this->ForwardMapping[timeSubdiv][ptVF(x,y,z,1)]=VecTemp2[1];
      this->ForwardMapping[timeSubdiv][ptVF(x,y,z,2)]=VecTemp2[2];
    }
  }
}




///compute the backward mapping
template <class VoxelType> void LargeDefGradLagrange<VoxelType>::ComputeBackwardMapping(void){
  float VecTemp[3];  //velocity vector that targets the current treated voxel
  float VecTemp2[3]; //temporary vector to propagate the mapping
  int ConvergenceSteps;
  int i,x,y,z;
  int timeSubdiv;
  int LastTimSub;

  //initialisation
  ConvergenceSteps=3;
  LastTimSub=this->NbTimeSubdiv-1;
  
  //JO at the last time subdivision
  for (z = 0; z < this->NZ; z++) for (y = 0; y < this->NY; y++) for (x = 0; x < this->NX; x++){
    this->BackwardMapping[LastTimSub][ptVF(x,y,z,0)]=static_cast<float>(x);
    this->BackwardMapping[LastTimSub][ptVF(x,y,z,1)]=static_cast<float>(y);
    this->BackwardMapping[LastTimSub][ptVF(x,y,z,2)]=static_cast<float>(z);
  }
  
  //JO at the other time subdivisions
  for (timeSubdiv=LastTimSub-1;timeSubdiv>=0;timeSubdiv--){
    for (z = 0; z < this->NZ; z++) for (y = 0; y < this->NY; y++) for (x = 0; x < this->NX; x++){
      //init
      VecTemp[0]=0.; 
      VecTemp[1]=0.;
      VecTemp[2]=0.;
      
      //convergence
      for (i=0;i<ConvergenceSteps;i++){
        this->GetVectorFromVelocityField(timeSubdiv+1,x+VecTemp[0],y+VecTemp[1],z+VecTemp[2],VecTemp);
        VecTemp[0]=(VecTemp[0]+this->VelocityField[timeSubdiv][ptVF(x,y,z,0)])*(this->DeltaTimeSubdiv)/2;
        VecTemp[1]=(VecTemp[1]+this->VelocityField[timeSubdiv][ptVF(x,y,z,1)])*(this->DeltaTimeSubdiv)/2;
        VecTemp[2]=(VecTemp[2]+this->VelocityField[timeSubdiv][ptVF(x,y,z,2)])*(this->DeltaTimeSubdiv)/2;
      }
      
      //new value
      this->GetCoordFromBackwardMapping(timeSubdiv+1,x+VecTemp[0],y+VecTemp[1],z+VecTemp[2],VecTemp2);
      
      this->BackwardMapping[timeSubdiv][ptVF(x,y,z,0)]=VecTemp2[0];
      this->BackwardMapping[timeSubdiv][ptVF(x,y,z,1)]=VecTemp2[1];
      this->BackwardMapping[timeSubdiv][ptVF(x,y,z,2)]=VecTemp2[2];

    }
  }
}



///compute J0 using the forward mapping   (J0 is the transformed template image at the time subdivision timeSubdiv)
template <class VoxelType> void LargeDefGradLagrange<VoxelType>::ComputeJ0(int timeSubdiv){
  int x,y,z;
  
  for (z = 0; z < this->NZ; z++) for (y = 0; y < this->NY; y++) for (x = 0; x < this->NX; x++){
    this->J0[ptSF(x,y,z)]=this->GetGreyLevelFromTemplate(this->ForwardMapping[timeSubdiv][ptVF(x,y,z,0)], this->ForwardMapping[timeSubdiv][ptVF(x,y,z,1)], this->ForwardMapping[timeSubdiv][ptVF(x,y,z,2)]);
  }
}


///compute J1 using the backward mapping   (J1 is the transformed target image at the time subdivision timeSubdiv)
template <class VoxelType> void LargeDefGradLagrange<VoxelType>::ComputeJ1(int timeSubdiv){
  int x,y,z;
  
  for (z = 0; z < this->NZ; z++) for (y = 0; y < this->NY; y++) for (x = 0; x < this->NX; x++){
    this->J1[ptSF(x,y,z)]=this->GetGreyLevelFromTarget(this->BackwardMapping[timeSubdiv][ptVF(x,y,z,0)],this->BackwardMapping[timeSubdiv][ptVF(x,y,z,1)],this->BackwardMapping[timeSubdiv][ptVF(x,y,z,2)]);
  }
}


///compute the gradient of J0
template <class VoxelType> void LargeDefGradLagrange<VoxelType>::ComputeGradientJ0(void){
  int x,y,z;
  
  //Energy gradient in direction x, y, z
  for (z = 1; z < this->NZ-2; z++) for (y = 1; y < this->NY-2; y++) for (x = 1; x < this->NX-2; x++){
    this->GradJ0[ptVF(x,y,z,0)]=(this->J0[ptSF(x+1,y,z)]-this->J0[ptSF(x-1,y,z)])/2.;
    this->GradJ0[ptVF(x,y,z,1)]=(this->J0[ptSF(x,y+1,z)]-this->J0[ptSF(x,y-1,z)])/2.;
    this->GradJ0[ptVF(x,y,z,2)]=(this->J0[ptSF(x,y,z+1)]-this->J0[ptSF(x,y,z-1)])/2.;
  }
  
  //boundary condition z=0
  for (y = 1; y < this->NY-1; y++) for (x = 1; x < this->NX-1; x++){
    this->GradJ0[ptVF(x,y,0,0)]=0;
    this->GradJ0[ptVF(x,y,0,1)]=0;
    this->GradJ0[ptVF(x,y,0,2)]=0;
  }
  
  //boundary condition z=this->NZ-1
  for (y = 1; y < this->NY-1; y++) for (x = 1; x < this->NX-1; x++){
    this->GradJ0[ptVF(x,y,this->NZ-1,0)]=0;
    this->GradJ0[ptVF(x,y,this->NZ-1,1)]=0;
    this->GradJ0[ptVF(x,y,this->NZ-1,2)]=0;
  }
  
  //2D images (this->NZ==1)
  for (y = 1; y < this->NY-2; y++) for (x = 1; x < this->NX-2; x++){
    this->GradJ0[ptVF(x,y,0,0)]=(this->J0[ptSF(x+1,y,0)]-this->J0[ptSF(x-1,y,0)])/2.;
    this->GradJ0[ptVF(x,y,0,1)]=(this->J0[ptSF(x,y+1,0)]-this->J0[ptSF(x,y-1,0)])/2.;
    this->GradJ0[ptVF(x,y,0,2)]=0.;
  }
  
  //boundary condition y=0
  for (z = 0; z < this->NZ; z++) for (x = 1; x < this->NX-1; x++){
    this->GradJ0[ptVF(x,0,z,0)]=0;
    this->GradJ0[ptVF(x,0,z,1)]=0;
    this->GradJ0[ptVF(x,0,z,2)]=0;
  }
  
  //boundary condition y=this->NY-1
  for (z = 0; z < this->NZ; z++) for (x = 1; x < this->NX-1; x++){
    this->GradJ0[ptVF(x,this->NY-1,z,0)]=0;
    this->GradJ0[ptVF(x,this->NY-1,z,1)]=0;
    this->GradJ0[ptVF(x,this->NY-1,z,2)]=0;
  }

  //boundary condition x=0
  for (z = 0; z < this->NZ; z++) for (y = 0; y < this->NY; y++){
    this->GradJ0[ptVF(0,y,z,0)]=0;
    this->GradJ0[ptVF(0,y,z,1)]=0;
    this->GradJ0[ptVF(0,y,z,2)]=0;
  }
  
  //boundary condition x=this->NX-1
  for (z = 0; z < this->NZ; z++) for (y = 0; y < this->NY; y++){
    this->GradJ0[ptVF(this->NX-1,y,z,0)]=this->GradJ0[ptVF(this->NX-2,y,z,0)];
    this->GradJ0[ptVF(this->NX-1,y,z,1)]=this->GradJ0[ptVF(this->NX-2,y,z,1)];
    this->GradJ0[ptVF(this->NX-1,y,z,2)]=this->GradJ0[ptVF(this->NX-2,y,z,2)];
  }
}

///Compute the determinant of the jacobian
template <class VoxelType> void LargeDefGradLagrange<VoxelType>::ComputeJacobianDeterminant(int timeSubdiv){
  int x,y,z,i,j;
  
  
  //default values (for the boundary conditions)
  for (z = 0; z < this->NZ; z++) for (y = 0; y < this->NY; y++) for (x = 0; x < this->NX; x++)
        this->DetJacobians[ptSF(x,y,z)] =1.;
  
  //3D image
  for (z = 1; z < this->NZ-2; z++) for (y = 1; y < this->NY-2; y++) for (x = 1; x < this->NX-2; x++){
    for(i=0;i<3;i++)for(j=0;j<3;j++){
      this->JacobianMatrix[3*i+j] = (this->BackwardMapping[timeSubdiv][ptVF(x + (int)(0==i),y + (int)(1==i),z + (int)(2==i),j)] - this->BackwardMapping[timeSubdiv][ptVF(x - (int)(0==i),y - (int)(1==i),z - (int)(2==i),j)])/2.;
    }
      this->DetJacobians[ptSF(x,y,z)] = this->CptDetJac3d(this->JacobianMatrix);
  }
  
  
  //2D image where this->NZ==1
  if (this->NZ==1){
    for (y = 1; y < this->NY-2; y++) for (x = 1; x < this->NX-2; x++){
      for(i=0;i<3;i++)for(j=0;j<3;j++) this->JacobianMatrix[3*i+j] =0;
      this->JacobianMatrix[3*2+2]=1.;
      for(i=0;i<2;i++)for(j=0;j<2;j++){
        this->JacobianMatrix[3*i+j] = (this->BackwardMapping[timeSubdiv][ptVF(x + (int)(0==i),y + (int)(1==i),0,j)] - this->BackwardMapping[timeSubdiv][ptVF(x - (int)(0==i),y - (int)(1==i),0,j)])/2.;
      }
      this->DetJacobians[ptSF(x,y,z)] = this->CptDetJac3d(this->JacobianMatrix);
    }
  }
}



///Compute the energy gradients
template <class VoxelType> void LargeDefGradLagrange<VoxelType>::ComputeEnergyGradient(int timeSubdiv){
  int x,y,z,i;
  float temp;
  
  //loop on the vector directions (x,y,z)
  for (i=0;i<3;i++){
    //compute the scalar field (one dimension out of the vector field) to smooth
    for (z=0;z<this->NZfft;z++) for (y=0;y<this->NYfft;y++) for (x=0;x<this->NXfft;x++) this->RealSignalForFFT.Put(x,y,z,0,0);
    for (z=0;z<this->NZfft;z++) for (y=0;y<this->NYfft;y++) for (x=0;x<this->NXfft;x++) this->ImagSignalForFFT.Put(x,y,z,0,0);
    
    for (z = 0; z < this->NZ; z++) for (y = 0; y < this->NY; y++) for (x = 0; x < this->NX; x++){
      temp=(this->J0[ptSF(x,y,z)] - this->J1[ptSF(x,y,z)]) * this->DetJacobians[ptSF(x,y,z)] * this->GradJ0[ptVF(x,y,z,i)];
      this->RealSignalForFFT.Put(SXfft+x,SYfft+y,SZfft+z,0,static_cast<float>(temp));
    }
    
     //smooth the scalar field
    ConvolutionInFourierNoFilterTransfo(&this->RealSignalForFFT,&this->ImagSignalForFFT,&this->RealFilterForFFT,&this->ImagFilterForFFT);
     
    //save the smoothed scalar field
    for (z = 0; z < this->NZ; z++) for (y = 0; y < this->NY; y++) for (x = 0; x < this->NX; x++){
      this->GradE[timeSubdiv][ptVF(x,y,z,i)] = this->WghtVelField*2*this->VelocityField[timeSubdiv][ptVF(x,y,z,i)] - 2*this->RealSignalForFFT.Get(SXfft+x,SYfft+y,SZfft+z,0);
    }
  }
}


///Update VelocityField with with the energy gradients
template <class VoxelType> void LargeDefGradLagrange<VoxelType>::UpdateVelocityField(int IterationsNb){
  int x, y, z, i;
  float MaxGrad,MultFactor,GradSVG;
  double LocGrad;
  
  //0) advice on the weigth to give to the kernels
  if ((IterationsNb==0)&&(this->NbKernels>1)){
    cout << "To have kernel of similar influence, put the following weights:\n";
    GradSVG=1;
    for (i=0;i<this->NbTimeSubdiv;i++){
      MaxGrad=0;
      //3D image
      for (z = 1; z < this->NZ-2; z++) for (y = 1; y < this->NY-2; y++) for (x = 1; x < this->NX-2; x++){
        LocGrad=sqrt(pow((double)this->GradE[i][ptVF(x,y,z,0)],2.0)+pow((double)this->GradE[i][ptVF(x,y,z,1)],2.0)+pow((double)this->GradE[i][ptVF(x,y,z,2)],2.0));
        if (MaxGrad<LocGrad) MaxGrad=(float)LocGrad;
      }
      //2D image
      if (this->NZ==1) for (y = 1; y < this->NY-2; y++) for (x = 1; x < this->NX-2; x++){
        LocGrad=sqrt(pow((double)this->GradE[i][ptVF(x,y,0,0)],2.0)+pow((double)this->GradE[i][ptVF(x,y,0,1)],2.0));
        if (MaxGrad<LocGrad) MaxGrad=(float)LocGrad;
      }
      
      //cout << i << " -> " << MaxGrad << "\n";
      
      if (i==0){
        GradSVG=MaxGrad;
        cout << "Kernel of time subdivision " << i << ": cst*" << 1. << "/[entered weight]\n";
      }
      else
        cout << "Kernel of time subdivision " << i << ": cst*" << GradSVG/MaxGrad << "/[entered weight]\n";
    }
  }
  
  //1) compute the multiplication factor...
  //...3D images
  MaxGrad=0;
  for (i=0;i<this->NbTimeSubdiv;i++) for (z = 1; z < this->NZ-2; z++) for (y = 1; y < this->NY-2; y++) for (x = 1; x < this->NX-2; x++){
    LocGrad=sqrt(pow((double)this->GradE[i][ptVF(x,y,z,0)],2.0)+pow((double)this->GradE[i][ptVF(x,y,z,1)],2.0)+pow((double)this->GradE[i][ptVF(x,y,z,2)],2.0));
    if (MaxGrad<LocGrad) MaxGrad=(float)LocGrad;
  }
  
  //...2D images
  if (this->NZ==1){
    for (i=0;i<this->NbTimeSubdiv;i++) for (y = 1; y < this->NY-2; y++) for (x = 1; x < this->NX-2; x++){
      LocGrad=sqrt(pow((double)this->GradE[i][ptVF(x,y,0,0)],2.0)+pow((double)this->GradE[i][ptVF(x,y,0,1)],2.0));
      if (MaxGrad<LocGrad) MaxGrad=(float)LocGrad;
    }
  }
  
  //2) maximum update control
  if (MaxGrad>this->MaxVelocityUpdate) MultFactor=this->MaxVelocityUpdate/MaxGrad;
  else MultFactor=0.1;
  if (MultFactor>0.1) MultFactor=0.1;
  
  cout << "MultFactor = " << MultFactor << "\n";
    
  //3) update the vector field...
  //...3D images
  for (i=0;i<this->NbTimeSubdiv;i++) for (z = 1; z < this->NZ-2; z++) for (y = 1; y < this->NY-2; y++) for (x = 1; x < this->NX-2; x++){
    this->VelocityField[i][ptVF(x,y,z,0)]=this->VelocityField[i][ptVF(x,y,z,0)]-this->GradE[i][ptVF(x,y,z,0)]*MultFactor;
    this->VelocityField[i][ptVF(x,y,z,1)]=this->VelocityField[i][ptVF(x,y,z,1)]-this->GradE[i][ptVF(x,y,z,1)]*MultFactor;
    this->VelocityField[i][ptVF(x,y,z,2)]=this->VelocityField[i][ptVF(x,y,z,2)]-this->GradE[i][ptVF(x,y,z,2)]*MultFactor;
  }
  
  //...2D images
  if (this->NZ==1){
    for (i=0;i<this->NbTimeSubdiv;i++) for (y = 1; y < this->NY-2; y++) for (x = 1; x < this->NX-2; x++){
      this->VelocityField[i][ptVF(x,y,0,0)]=this->VelocityField[i][ptVF(x,y,0,0)]-this->GradE[i][ptVF(x,y,0,0)]*MultFactor;
      this->VelocityField[i][ptVF(x,y,0,1)]=this->VelocityField[i][ptVF(x,y,0,1)]-this->GradE[i][ptVF(x,y,0,1)]*MultFactor;
    }
  }
}


///save the result of the gradient descent (Beg 2005) for the current 3D image of the 4D time sequence
template <class VoxelType> void LargeDefGradLagrange<VoxelType>::SaveResultGradientDescent(int TimeLoc){
  int x, y, z;
  
  
  //init -> compute the forward mapping and import the original input template (non pre-treated)
  this->ComputeForwardMapping();
  
  for (z = 0; z < this->NZ; z++)  for (y = 0; y < this->NY; y++) for (x = 0; x < this->NX; x++)
        this->ImTemplate[ptSF(x,y,z)]=static_cast<float>(this->_input->Get(x, y, z, TimeLoc));
  
  /*VIRE  VIRE  VIRE  VIRE  VIRE  VIRE  VIRE  VIRE  VIRE  VIRE*/
/*  for (z = 0; z < this->NZ; z++)  for (y = 0; y < this->NY; y++) for (x = 0; x < this->NX; x++)
        this->ImTemplate[ptSF(x,y,z)]=static_cast<float>(0);
  
  for (z = 0; z < this->NZ; z++)  for (y = 0; y < this->NY; y++) for (x = 0; x < this->NX; x++){
    //if ((z>0)&&(z<this->NZ-1)&&(x>=20+70)&&(x<=29+70)&&(y>=20+70)&&(y<=29+70)) this->ImTemplate[ptSF(x,y,z)]=static_cast<float>(100);
    if ((z>1)&&(z<this->NZ-2)&&(x>=49)&&(x<=60)&&(y>=54)&&(y<=65)) this->ImTemplate[ptSF(x,y,z)]=static_cast<float>(100);
    if ((z>1)&&(z<this->NZ-2)&&(x>=49)&&(x<=75)&&(y>=66)&&(y<=71)) this->ImTemplate[ptSF(x,y,z)]=static_cast<float>(100);
    if ((z>1)&&(z<this->NZ-2)&&(x>=64)&&(x<=75)&&(y>=54)&&(y<=65)) this->ImTemplate[ptSF(x,y,z)]=static_cast<float>(100);
  }*/
  /*VIRE  VIRE  VIRE  VIRE  VIRE  VIRE  VIRE  VIRE  VIRE  VIRE*/
  
  //deform the original template
  this->ComputeJ0(this->NbTimeSubdiv-1);
  
  //save the deformed template in the output of the class
  for (z = 0; z < this->NZ; z++) for (y = 0; y < this->NY; y++) for (x = 0; x < this->NX; x++){
      this->_output->Put(x, y, z, TimeLoc, static_cast<VoxelType>(this->J0[ptSF(x,y,z)]));
  }
  
  //if required, save the velocity field and the template deformations across the time subdivisions
  if (strcmp(this->PrefixOutputVF,"Null")!=0){
    this->SaveVelocityFields(this->PrefixOutputVF);
    this->SaveDeformations(this->PrefixOutputVF);
    this->SaveGridDeformations(this->PrefixOutputVF);
    //this->SaveForwardMapping(this->PrefixOutputVF);
  }
  
}


///save the deformations in time subdivisions (not the convergence)
template <class VoxelType> void LargeDefGradLagrange<VoxelType>::SaveDeformations(char Prefix[256]){
  int TimeLoc,x, y, z;
  irtkGenericImage<VoxelType> Temp4DField;
  char FileName[256];
  char Deformations[256];
  
  //intialisation
  Temp4DField = irtkGenericImage<float>(this->NX, this->NY, this->NZ,this->NbTimeSubdiv);
  strcpy(Deformations,"_Deformations.nii");
  
  //save the deformations
  for (TimeLoc=0;TimeLoc<this->NbTimeSubdiv;TimeLoc++){
    this->ComputeJ0(TimeLoc);
    
    for (z = 0; z < this->NZ; z++) for (y = 0; y < this->NY; y++) for (x = 0; x < this->NX; x++){
      Temp4DField.Put(x, y, z, TimeLoc, this->J0[ptSF(x,y,z)]);
    }
  }
  strcpy(FileName,Prefix);
  strcat(FileName,Deformations);
  Temp4DField.Write(FileName);
}



///save the deformations of a grid in time subdivisions
template <class VoxelType> void LargeDefGradLagrange<VoxelType>::SaveGridDeformations(char Prefix[256]){
  int TimeLoc,x, y, z;
  irtkGenericImage<VoxelType> Temp4DField;
  char FileName[256];
  char Deformations[256];
  double epsilon;
  
  epsilon=0.000001;
  
  //intialisation of the field to save
  Temp4DField = irtkGenericImage<float>(this->NX, this->NY, this->NZ,this->NbTimeSubdiv);
  
  //initialisation of the grid
  for (z = 0; z < this->NZ; z++)  for (y = 0; y < this->NY; y++) for (x = 0; x < this->NX; x++)
        this->ImTemplate[ptSF(x,y,z)]=static_cast<float>(0);
  
  if (this->NZ>1) for (z = 0; z < this->NZ; z++)  for (y = 0; y < this->NY; y++) for (x = 0; x < this->NX; x++)
          if ( ( abs((double)(z/5)-(((double)z)/5.)) < epsilon ) || ( abs((double)(y/5)-(((double)y)/5.)) < epsilon )|| ( abs((double)(x/5)-(((double)x)/5.)) < epsilon ) )
            this->ImTemplate[ptSF(x,y,z)]=static_cast<float>(100);
  
  if (this->NZ==1) for (y = 0; y < this->NY; y++) for (x = 0; x < this->NX; x++)
          if ( ( abs((double)(y/5)-(((double)y)/5.)) < epsilon )|| ( abs((double)(x/5)-(((double)x)/5.)) < epsilon ) )
            this->ImTemplate[ptSF(x,y,0)]=static_cast<float>(100);
  
  //compute the deformations of the grid
  for (TimeLoc=0;TimeLoc<this->NbTimeSubdiv;TimeLoc++){
    this->ComputeJ0(TimeLoc);
    
    for (z = 0; z < this->NZ; z++) for (y = 0; y < this->NY; y++) for (x = 0; x < this->NX; x++){
      Temp4DField.Put(x, y, z, TimeLoc, this->J0[ptSF(x,y,z)]);
    }
  }
  
  //save the deformations of the grid
  strcpy(Deformations,"_GridDef.nii");
  strcpy(FileName,Prefix);
  strcat(FileName,Deformations);
  Temp4DField.Write(FileName);
  
  //recompute the original template (without preprocessing)
  for (z = 0; z < this->NZ; z++)  for (y = 0; y < this->NY; y++) for (x = 0; x < this->NX; x++)
    this->ImTemplate[ptSF(x,y,z)]=static_cast<float>(this->_input->Get(x, y, z, 0));
}



///save the velocity fields
template <class VoxelType> void LargeDefGradLagrange<VoxelType>::SaveVelocityFields(char Prefix[256]){
  int TimeLoc,x, y, z;
  irtkGenericImage<VoxelType> Temp4DField;
  char FileName[256];
  char VelocityField_X[256];
  char VelocityField_Y[256];
  char VelocityField_Z[256];
  
  //intialisation
  Temp4DField = irtkGenericImage<float>(this->NX, this->NY, this->NZ,this->NbTimeSubdiv);
  strcpy(VelocityField_X,"_VelocityField_X.nii");
  strcpy(VelocityField_Y,"_VelocityField_Y.nii");
  strcpy(VelocityField_Z,"_VelocityField_Z.nii");
  
  //save the velocity field in direction X
  for (TimeLoc=0;TimeLoc<this->NbTimeSubdiv;TimeLoc++)  for (z = 0; z < this->NZ; z++) for (y = 0; y < this->NY; y++) for (x = 0; x < this->NX; x++)
          Temp4DField.Put(x, y, z, TimeLoc, this->VelocityField[TimeLoc][ptVF(x,y,z,0)]*100);
  strcpy(FileName,Prefix);
  strcat(FileName,VelocityField_X);
  Temp4DField.Write(FileName);

  //save the velocity field in direction Y
  for (TimeLoc=0;TimeLoc<this->NbTimeSubdiv;TimeLoc++)  for (z = 0; z < this->NZ; z++) for (y = 0; y < this->NY; y++) for (x = 0; x < this->NX; x++)
          Temp4DField.Put(x, y, z, TimeLoc, this->VelocityField[TimeLoc][ptVF(x,y,z,1)]*100);
  strcpy(FileName,Prefix);
  strcat(FileName,VelocityField_Y);
  Temp4DField.Write(FileName);

  //save the velocity field in direction Z
  for (TimeLoc=0;TimeLoc<this->NbTimeSubdiv;TimeLoc++)  for (z = 0; z < this->NZ; z++) for (y = 0; y < this->NY; y++) for (x = 0; x < this->NX; x++)
          Temp4DField.Put(x, y, z, TimeLoc, this->VelocityField[TimeLoc][ptVF(x,y,z,2)]*100);
  strcpy(FileName,Prefix);
  strcat(FileName,VelocityField_Z);
  Temp4DField.Write(FileName);
  
  
  
  cout << "REMARK: THE VELOCITIES ARE SAVED AS INTEGERS AND THEN ARE MULTIPLIED BY 100\n";
}


///save the velocity fields
template <class VoxelType> void LargeDefGradLagrange<VoxelType>::SaveForwardMapping(char Prefix[256]){
  int TimeLoc,x, y, z;
  irtkGenericImage<VoxelType> Temp4DField;
  char FileName[256];
  char ForwardMapping_X[256];
  char ForwardMapping_Y[256];
  char ForwardMapping_Z[256];
  char ForwardMapping_N[256];
  float tempFloat;
  
  //intialisation
  Temp4DField = irtkGenericImage<float>(this->NX, this->NY, this->NZ,this->NbTimeSubdiv);
  strcpy(ForwardMapping_X,"_ForwardMapping_X.nii");
  strcpy(ForwardMapping_Y,"_ForwardMapping_Y.nii");
  strcpy(ForwardMapping_Z,"_ForwardMapping_Z.nii");
  strcpy(ForwardMapping_N,"_ForwardMapping_Norm.nii");
  
  //save the forward mapping in direction X
  for (TimeLoc=0;TimeLoc<this->NbTimeSubdiv;TimeLoc++)  for (z = 0; z < this->NZ; z++) for (y = 0; y < this->NY; y++) for (x = 0; x < this->NX; x++)
          Temp4DField.Put(x, y, z, TimeLoc, (this->ForwardMapping[TimeLoc][ptVF(x,y,z,0)]-this->ForwardMapping[0][ptVF(x,y,z,0)])*100);
  strcpy(FileName,Prefix);
  strcat(FileName,ForwardMapping_X);
  Temp4DField.Write(FileName);

  //save the forward mapping in direction Y
  for (TimeLoc=0;TimeLoc<this->NbTimeSubdiv;TimeLoc++)  for (z = 0; z < this->NZ; z++) for (y = 0; y < this->NY; y++) for (x = 0; x < this->NX; x++)
          Temp4DField.Put(x, y, z, TimeLoc, (this->ForwardMapping[TimeLoc][ptVF(x,y,z,1)]-this->ForwardMapping[0][ptVF(x,y,z,1)])*100);
  strcpy(FileName,Prefix);
  strcat(FileName,ForwardMapping_Y);
  Temp4DField.Write(FileName);

  //save the forward mapping in direction Z
  for (TimeLoc=0;TimeLoc<this->NbTimeSubdiv;TimeLoc++)  for (z = 0; z < this->NZ; z++) for (y = 0; y < this->NY; y++) for (x = 0; x < this->NX; x++)
          Temp4DField.Put(x, y, z, TimeLoc, (this->ForwardMapping[TimeLoc][ptVF(x,y,z,2)]-this->ForwardMapping[0][ptVF(x,y,z,2)])*100);
  strcpy(FileName,Prefix);
  strcat(FileName,ForwardMapping_Z);
  Temp4DField.Write(FileName);
  
  //save the norm of the forward mapping
  for (TimeLoc=0;TimeLoc<this->NbTimeSubdiv;TimeLoc++)  for (z = 0; z < this->NZ; z++) for (y = 0; y < this->NY; y++) for (x = 0; x < this->NX; x++){
    tempFloat= (this->ForwardMapping[TimeLoc][ptVF(x,y,z,0)]-this->ForwardMapping[0][ptVF(x,y,z,0)])*(this->ForwardMapping[TimeLoc][ptVF(x,y,z,0)]-this->ForwardMapping[0][ptVF(x,y,z,0)]);
    tempFloat+=(this->ForwardMapping[TimeLoc][ptVF(x,y,z,1)]-this->ForwardMapping[0][ptVF(x,y,z,1)])*(this->ForwardMapping[TimeLoc][ptVF(x,y,z,1)]-this->ForwardMapping[0][ptVF(x,y,z,1)]);
    tempFloat+=(this->ForwardMapping[TimeLoc][ptVF(x,y,z,2)]-this->ForwardMapping[0][ptVF(x,y,z,2)])*(this->ForwardMapping[TimeLoc][ptVF(x,y,z,2)]-this->ForwardMapping[0][ptVF(x,y,z,2)]);
    tempFloat=(float)sqrt((double)tempFloat);
    Temp4DField.Put(x, y, z, TimeLoc, tempFloat*100);
  }
  
  strcpy(FileName,Prefix);
  strcat(FileName,ForwardMapping_N);
  Temp4DField.Write(FileName);
  
  cout << "REMARK: THE MAPPINGS ARE SAVED AS INTEGERS AND THEN ARE MULTIPLIED BY 100\n";
}


///load the velocity fields
template <class VoxelType> void LargeDefGradLagrange<VoxelType>::LoadVelocityFields(char Prefix[256]){
  int TimeLoc,x, y, z;
  irtkGenericImage<float> Temp4DField;
  char FileName[256];
  char VelocityField_X[256];
  char VelocityField_Y[256];
  char VelocityField_Z[256];
  float Xfactor_x,Yfactor_x,Zfactor_x,Tfactor_x;
  float Xfactor_y,Yfactor_y,Zfactor_y,Tfactor_y;
  float Xfactor_z,Yfactor_z,Zfactor_z,Tfactor_z;
  
  //1) intialisation
  //Temp4DField = irtkGenericImage<float>(this->NX, this->NY, this->NZ,this->NbTimeSubdiv);
  strcpy(VelocityField_X,"_VelocityField_X.nii");
  strcpy(VelocityField_Y,"_VelocityField_Y.nii");
  strcpy(VelocityField_Z,"_VelocityField_Z.nii");
  
  
  //2) load the velocity field in direction X
  strcpy(FileName,Prefix);
  strcat(FileName,VelocityField_X);
  Temp4DField.Read(FileName);
  
  //multiplication factors
  Xfactor_x=0.000001+static_cast<float>(Temp4DField.GetX())/static_cast<float>(this->NX);
  Yfactor_x=0.000001+static_cast<float>(Temp4DField.GetY())/static_cast<float>(this->NY);
  Zfactor_x=0.000001+static_cast<float>(Temp4DField.GetZ())/static_cast<float>(this->NZ);
  Tfactor_x=0.000001+static_cast<float>(Temp4DField.GetT())/static_cast<float>(this->NbTimeSubdiv);
  
  //copy the values
  if ((this->NX==Temp4DField.GetX())&&(this->NY==Temp4DField.GetY())&&(this->NZ==Temp4DField.GetZ())&&(this->NbTimeSubdiv==Temp4DField.GetT())){
    for (TimeLoc=0;TimeLoc<this->NbTimeSubdiv;TimeLoc++)  for (z = 0; z < this->NZ; z++) for (y = 0; y < this->NY; y++) for (x = 0; x < this->NX; x++)
            this->VelocityField[TimeLoc][ptVF(x,y,z,0)]=static_cast<float>(Temp4DField.Get(x,y,z,TimeLoc))/100;  //TO REMOVE "/100"
  }
  else{
    cout << "Interpolation on " << FileName << "\n";
    for (TimeLoc=0;TimeLoc<this->NbTimeSubdiv;TimeLoc++)  for (z = 0; z < this->NZ; z++) for (y = 0; y < this->NY; y++) for (x = 0; x < this->NX; x++)
            //this->VelocityField[TimeLoc][ptVF(x,y,z,0)]=static_cast<float>(Temp4DField.Get((int)(x*Xfactor_x),(int)(y*Yfactor_x),(int)(z*Zfactor_x),(int)(TimeLoc*Tfactor_x)))/Xfactor_x;
            this->VelocityField[TimeLoc][ptVF(x,y,z,0)]=GetGreyLevelFromFloatGenericImage(&Temp4DField,x*Xfactor_x, y*Yfactor_x, z*Zfactor_x, TimeLoc*Tfactor_x)/(Xfactor_x*100);  //TO REMOVE "/100"
  }
  
  
  //3) load the velocity field in direction Y
  strcpy(FileName,Prefix);
  strcat(FileName,VelocityField_Y);
  Temp4DField.Read(FileName);
  
  //multiplication factors
  Xfactor_y=0.000001+static_cast<float>(Temp4DField.GetX())/static_cast<float>(this->NX);
  Yfactor_y=0.000001+static_cast<float>(Temp4DField.GetY())/static_cast<float>(this->NY);
  Zfactor_y=0.000001+static_cast<float>(Temp4DField.GetZ())/static_cast<float>(this->NZ);
  Tfactor_y=0.000001+static_cast<float>(Temp4DField.GetT())/static_cast<float>(this->NbTimeSubdiv);
  
  //copy the values
  if ((this->NX==Temp4DField.GetX())&&(this->NY==Temp4DField.GetY())&&(this->NZ==Temp4DField.GetZ())&&(this->NbTimeSubdiv==Temp4DField.GetT())){
    for (TimeLoc=0;TimeLoc<this->NbTimeSubdiv;TimeLoc++)  for (z = 0; z < this->NZ; z++) for (y = 0; y < this->NY; y++) for (x = 0; x < this->NX; x++)
            this->VelocityField[TimeLoc][ptVF(x,y,z,1)]=static_cast<float>(Temp4DField.Get(x,y,z,TimeLoc))/100;  //TO REMOVE "/100"
  }
  else{
    cout << "Interpolation on " << FileName << "\n";
    for (TimeLoc=0;TimeLoc<this->NbTimeSubdiv;TimeLoc++)  for (z = 0; z < this->NZ; z++) for (y = 0; y < this->NY; y++) for (x = 0; x < this->NX; x++)
            //this->VelocityField[TimeLoc][ptVF(x,y,z,1)]=static_cast<float>(Temp4DField.Get((int)(x*Xfactor_y),(int)(y*Yfactor_y),(int)(z*Zfactor_y),(int)(TimeLoc*Tfactor_y)))/Yfactor_y;
            this->VelocityField[TimeLoc][ptVF(x,y,z,1)]=GetGreyLevelFromFloatGenericImage(&Temp4DField,x*Xfactor_y, y*Yfactor_y, z*Zfactor_y, TimeLoc*Tfactor_y)/(Yfactor_y*100);  //TO REMOVE "/100"
  }
  
  
  //4) load the velocity field in direction Z
  strcpy(FileName,Prefix);
  strcat(FileName,VelocityField_Z);
  Temp4DField.Read(FileName);
  
  //multiplication factors
  Xfactor_z=0.000001+static_cast<float>(Temp4DField.GetX())/static_cast<float>(this->NX);
  Yfactor_z=0.000001+static_cast<float>(Temp4DField.GetY())/static_cast<float>(this->NY);
  Zfactor_z=0.000001+static_cast<float>(Temp4DField.GetZ())/static_cast<float>(this->NZ);
  Tfactor_z=0.000001+static_cast<float>(Temp4DField.GetT())/static_cast<float>(this->NbTimeSubdiv);
  
  //copy the values
  if ((this->NX==Temp4DField.GetX())&&(this->NY==Temp4DField.GetY())&&(this->NZ==Temp4DField.GetZ())&&(this->NbTimeSubdiv==Temp4DField.GetT())){
    for (TimeLoc=0;TimeLoc<this->NbTimeSubdiv;TimeLoc++)  for (z = 0; z < this->NZ; z++) for (y = 0; y < this->NY; y++) for (x = 0; x < this->NX; x++)
            this->VelocityField[TimeLoc][ptVF(x,y,z,2)]=static_cast<float>(Temp4DField.Get(x,y,z,TimeLoc))/100.;  //TO REMOVE "/100"
  }
  else{
    cout << "Interpolation on " << FileName << "\n";
    for (TimeLoc=0;TimeLoc<this->NbTimeSubdiv;TimeLoc++)  for (z = 0; z < this->NZ; z++) for (y = 0; y < this->NY; y++) for (x = 0; x < this->NX; x++)
            //this->VelocityField[TimeLoc][ptVF(x,y,z,2)]=static_cast<float>(Temp4DField.Get((int)(x*Xfactor_z),(int)(y*Yfactor_z),(int)(z*Zfactor_z),(int)(TimeLoc*Tfactor_z)))/Zfactor_z;
            this->VelocityField[TimeLoc][ptVF(x,y,z,2)]=GetGreyLevelFromFloatGenericImage(&Temp4DField,x*Xfactor_z, y*Yfactor_z, z*Zfactor_z, TimeLoc*Tfactor_z)/(Zfactor_z*100);  //TO REMOVE "/100"
  }
  
  //5) warnings
  
  if ((Xfactor_x!=Xfactor_y)||(Xfactor_x!=Xfactor_z)||(Yfactor_x!=Yfactor_y)||(Yfactor_x!=Yfactor_z)||(Zfactor_x!=Zfactor_y)||(Zfactor_x!=Zfactor_z)||(Tfactor_x!=Tfactor_y)||(Tfactor_x!=Tfactor_z))
    cout << "Warning: all input velocity fields have not the same size!\n";

  if ((abs(Xfactor_x/Yfactor_x)>1.01)||(abs(Xfactor_x/Zfactor_x)>1.01)||(abs(Xfactor_x/Yfactor_x)<0.99)||(abs(Xfactor_x/Zfactor_x)<0.99))
    cout << "Warning: Anisotropic resampling in x direction!\n";

  if ((abs(Xfactor_y/Yfactor_y)>1.01)||(abs(Xfactor_y/Zfactor_y)>1.01)||(abs(Xfactor_y/Yfactor_y)<0.99)||(abs(Xfactor_y/Zfactor_y)<0.99))
    cout << "Warning: Anisotropic resampling in y direction!\n";

  if ((abs(Xfactor_z/Yfactor_z)>1.01)||(abs(Xfactor_z/Zfactor_z)>1.01)||(abs(Xfactor_z/Yfactor_z)<0.99)||(abs(Xfactor_z/Zfactor_z)<0.99))
    cout << "Warning: Anisotropic resampling in z direction!\n";

  //cout << "X factor=" << Xfactor_x << ", Y factor=" << Yfactor_x << ", Z factor=" << Zfactor_x << ", T factor=" << Tfactor_x << "\n";

}

///++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
///     SUB-FUNCTIONS TO PERFORM THE REGISTRATION (2: Energies and velocity field reparameterization)
///++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

/// compute the norm of the vector field at each time: this is the L^2 product of the vector field.
/// Also compute the total length
template <class VoxelType> void LargeDefGradLagrange<VoxelType>::ComputeNormsAndLength(){
  int t,z,y,x,d;
  double normTmp;
  
  //1) compute norm
  for (t=0;t<this->NbTimeSubdiv;t++){
    //energy at the time subdivision t...
    normTmp=0;
    
    for(d=0;d<3;d++){
      for (z=0;z<this->NZfft;z++) for (y=0;y<this->NYfft;y++) for (x=0;x<this->NXfft;x++) this->RealSignalForFFT.Put(x,y,z,0,0);
      for (z=0;z<this->NZfft;z++) for (y=0;y<this->NYfft;y++) for (x=0;x<this->NXfft;x++) this->ImagSignalForFFT.Put(x,y,z,0,0);
      
      for (z = 0; z < this->NZ; z++) for (y = 0; y < this->NY; y++) for (x = 0; x < this->NX; x++)
            this->RealSignalForFFT.Put(this->SXfft+x,this->SYfft+y,this->SZfft+z,0,static_cast<float>(this->VelocityField[t][ptVF(x,y,z,d)]));
      
      ConvolutionInFourierNoFilterTransfo(&this->RealSignalForFFT,&this->ImagSignalForFFT,&this->RealFilterForFFT,&this->ImagFilterForFFT);
      
      for (z = 0; z < this->NZ; z++) for (y = 0; y < this->NY; y++) for (x = 0; x < this->NX; x++)
            normTmp+=(double)(this->RealSignalForFFT.Get(this->SXfft+x,this->SYfft+y,this->SZfft+z,0)*this->VelocityField[t][ptVF(x,y,z,d)]);
    }
    
    this->norm[t] = (float)sqrt(normTmp);
  }

  //2) compute norm
  this->length = 0;
  for (t=0;t<this->NbTimeSubdiv;t++) this->length += this->norm[t]/this->NbTimeSubdiv;
}




///Compute the objective function that is minimized in the gradient descent:
/// !!! Requires to have an up to date this->norm (launch ComputeNormsAndLength) !!!
template <class VoxelType> void LargeDefGradLagrange<VoxelType>::ComputeCostFunction(){
  int k,z,y,x;
  
  // reinitialisation of the energy variables
  this->EnergyVF = 0;
  this->EnergyDI = 0;

  //cost from the velocity field
  for (k=0;k<this->NbTimeSubdiv;k++){this->EnergyVF+=this->norm[k]*this->norm[k];}
  this->EnergyVF /=  this->NbTimeSubdiv;
  
  
  //cost of the matching between the deformed template and the target
  this->ComputeJ0(this->NbTimeSubdiv-1);
  for (z = 0; z < this->NZ; z++) for (y = 0; y < this->NY; y++) for (x = 0; x < this->NX; x++){
    this->EnergyDI += (this->J0[ptSF(x,y,z)] - this->ImTarget[ptSF(x,y,z)]) *(this->J0[ptSF(x,y,z)] - this->ImTarget[ptSF(x,y,z)]);
  }
  
  this->EnergyDI/=this->NX*this->NY*this->NZ;
  
  //total cost
  this->EnergyTot = this->WghtVelField *this->EnergyVF + this->EnergyDI;
}



/// Compute the time reparameterization by a linear approximation to be used in the speed reparameterization
template <class VoxelType> float* LargeDefGradLagrange<VoxelType>::ComputeLinearReparameterizationMapping(float time){
  float temp=0;
  float temp2=0;
  int i=-1;
  float* result=new float[3];  // result[0] is the inf time to interpolate
                              // result[1] is the sup time to interpolate
                             // result[2] is the barycenter to interpolate result[0] and result[2] for the linear approximation
  
  while ((temp<this->length*(time+1))&&(i<this->NbTimeSubdiv - 1)){
    temp2 = temp;
    i +=1;
    temp += this->norm[i];
  }
  if (temp==temp2){
    result[0] = (i>=0)?1:0 * i;
    result[1] = (i>=0)?1:0 * i;
    result[2] = 0;
  }
  else{
    result[0] = i-1;
    result[1] = i;
    result[2] = (-temp2 + this->length*(time+1))/(temp-temp2);
  }
  return result;
}




/// Reparameterization of the velocity field to have a constant speed velocity field with (almost) equal final tranformation
/// Assuming that ComputeNorm and ComputeLength are already computed.
template <class VoxelType> void LargeDefGradLagrange<VoxelType>::SpeedReparameterization(){
  int t,x,y,z,d,t0,t1;
  float barycenter,rescaleFactor;
          //reparamHelper is the result of ComputeLinearReparameterizationMapping needed for the linear time interpolation.
  float* reparamHelper;
  
  for (t=0;t<this->NbTimeSubdiv;t++){
    reparamHelper = this->ComputeLinearReparameterizationMapping(t);
    t0 = reparamHelper[0];
    t1 = reparamHelper[1];
    barycenter = reparamHelper[2];
    //cout << "t0: "<< t0 <<"/ t1: "<< t1 << "/ barycentre: " << barycenter << "/ Length: " << this->length << "\n";
    if(t0==-1){
      for (z = 0; z < this->NZ; z++) for (y = 0; y < this->NY; y++) for (x = 0; x < this->NX; x++) for (d = 0; d < 3; d++){
        this->VelocityFieldTemp[t][ptVF(x,y,z,d)]  = barycenter * this->VelocityField[t1][ptVF(x,y,z,d)];
      }
    }
    else{
      for (z = 0; z < this->NZ; z++) for (y = 0; y < this->NY; y++) for (x = 0; x < this->NX; x++) for (d = 0; d < 3; d++){
        this->VelocityFieldTemp[t][ptVF(x,y,z,d)]  = (1-barycenter) * this->VelocityField[t0][ptVF(x,y,z,d)] + barycenter * this->VelocityField[t1][ptVF(x,y,z,d)];
      }
    }
  }
  for (t=0;t<this->NbTimeSubdiv;t++){
    for (z = 0; z < this->NZ; z++) for (y = 0; y < this->NY; y++) for (x = 0; x < this->NX; x++) for (d = 0; d < 3; d++){
      this->VelocityField[t][ptVF(x,y,z,d)] = this->VelocityFieldTemp[t][ptVF(x,y,z,d)];
    }
  }
  
  // +++ Important to call this function here +++
  this->ComputeNormsAndLength();

  for (t=0;t<this->NbTimeSubdiv;t++){
    if (this->norm[t]>0){
      rescaleFactor = this->length/this->norm[t];
      for (z = 0; z < this->NZ; z++) for (y = 0; y < this->NY; y++) for (x = 0; x < this->NX; x++) for (d = 0; d < 3; d++){
        this->VelocityField[t][ptVF(x,y,z,d)] *= rescaleFactor;
      }
    }
  }
  //cout << "Speed reparameterization done"<< "\n";
}



///++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
///     SUB-FUNCTIONS TO PERFORM THE REGISTRATION (3: Copies and multi-kernel stuffs)
///++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

template <class VoxelType> void LargeDefGradLagrange<VoxelType>::KernelManager(int TimeSubdiv){
  int CycleID;
  double LocTime;
  
  //timesubdivision nomalized between 0 and 1
  LocTime=((double)TimeSubdiv)/((double)this->NbTimeSubdiv-1);
  
  //choice of the kernel
  if (this->NbKernels==1) CycleID=1;
  if (this->NbKernels==2){
    if (LocTime<0.5) CycleID=1;
    else CycleID=2;
  }
  if (this->NbKernels==3){
    if (LocTime<0.3333) CycleID=1;
    else if (LocTime<0.6666) CycleID=2;
    else CycleID=3;
  }
  if (this->NbKernels==4){
    if (LocTime<0.25) CycleID=1;
    else if (LocTime<0.5) CycleID=2;
    else if (LocTime<0.75) CycleID=3;
    else CycleID=4;
  }
  
  //compute the kernel
  if (CycleID==1){
    //cout << "CycleID==1\n";
    MakeAnisotropicGaussianFilter(this->weight1,this->sigmaX1,this->sigmaY1,this->sigmaZ1,&this->RealFilterForFFT,&this->ImagFilterForFFT);
    DirectFFT(&this->RealFilterForFFT,&this->ImagFilterForFFT);
  }
  if (CycleID==2){
    //cout << "CycleID==2\n";
    MakeAnisotropicGaussianFilter(this->weight2,this->sigmaX2,this->sigmaY2,this->sigmaZ2,&this->RealFilterForFFT,&this->ImagFilterForFFT);
    DirectFFT(&this->RealFilterForFFT,&this->ImagFilterForFFT);
  }
  if (CycleID==3){
    //cout << "CycleID==3\n";
    MakeAnisotropicGaussianFilter(this->weight3,this->sigmaX3,this->sigmaY3,this->sigmaZ3,&this->RealFilterForFFT,&this->ImagFilterForFFT);
    DirectFFT(&this->RealFilterForFFT,&this->ImagFilterForFFT);
  }
  if (CycleID==4){
    //cout << "CycleID==4\n";
    MakeAnisotropicGaussianFilter(this->weight4,this->sigmaX4,this->sigmaY4,this->sigmaZ4,&this->RealFilterForFFT,&this->ImagFilterForFFT);
    DirectFFT(&this->RealFilterForFFT,&this->ImagFilterForFFT);
  }
}

///++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
///            SUB-FUNCTIONS TO PERFORM THE REGISTRATION  (4: interpolations)
///++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


//Perform a trilinear interplolation and returns the velocity field at non integer coordinates.
//++++ The output is in VecTemp ++++ 
template <class VoxelType> void LargeDefGradLagrange<VoxelType>::GetVectorFromVelocityField(int subivId,float x, float y, float z,float VecTemp[3]){
  int xi,yi,zi;
  float xwm,ywm,zwm,xwp,ywp,zwp;
  float wmmm,wmmp,wmpm,wmpp,wpmm,wpmp,wppm,wppp;
  
  xi=static_cast<int>(x);  xwm=1-(x-static_cast<float>(xi));  xwp=x-static_cast<float>(xi);
  yi=static_cast<int>(y);  ywm=1-(y-static_cast<float>(yi));  ywp=y-static_cast<float>(yi);
  zi=static_cast<int>(z);  zwm=1-(z-static_cast<float>(zi));  zwp=z-static_cast<float>(zi);
  
  if ((xi<0.)||(xi>=this->NX-1)||(yi<0.)||(yi>=this->NY-1)||(zi<0.)||(zi>=this->NZ-1)){
    VecTemp[0]=0;
    VecTemp[1]=0;
    VecTemp[2]=0;
    if (this->NZ==1) GetVectorFrom2DVelocityField(subivId,x,y,z,VecTemp); //2D image
  }
  else{
  wmmm=xwm*ywm*zwm;
  wmmp=xwm*ywm*zwp;
  wmpm=xwm*ywp*zwm;
  wmpp=xwm*ywp*zwp;
  wpmm=xwp*ywm*zwm;
  wpmp=xwp*ywm*zwp;
  wppm=xwp*ywp*zwm;
  wppp=xwp*ywp*zwp;
  
  //cout << "sum = " << wmmm+wmmp+wmpm+wmpp+wpmm+wpmp+wppm+wppp << "\n";
  
  VecTemp[0]= wmmm*this->VelocityField[subivId][ptVF(xi,yi,zi,0)];
  VecTemp[0]+=wmmp*this->VelocityField[subivId][ptVF(xi,yi,zi+1,0)];
  VecTemp[0]+=wmpm*this->VelocityField[subivId][ptVF(xi,yi+1,zi,0)];
  VecTemp[0]+=wmpp*this->VelocityField[subivId][ptVF(xi,yi+1,zi+1,0)];
  VecTemp[0]+=wpmm*this->VelocityField[subivId][ptVF(xi+1,yi,zi,0)];
  VecTemp[0]+=wpmp*this->VelocityField[subivId][ptVF(xi+1,yi,zi+1,0)];
  VecTemp[0]+=wppm*this->VelocityField[subivId][ptVF(xi+1,yi+1,zi,0)];
  VecTemp[0]+=wppp*this->VelocityField[subivId][ptVF(xi+1,yi+1,zi+1,0)];
  
  VecTemp[1]= wmmm*this->VelocityField[subivId][ptVF(xi,yi,zi,1)];
  VecTemp[1]+=wmmp*this->VelocityField[subivId][ptVF(xi,yi,zi+1,1)];
  VecTemp[1]+=wmpm*this->VelocityField[subivId][ptVF(xi,yi+1,zi,1)];
  VecTemp[1]+=wmpp*this->VelocityField[subivId][ptVF(xi,yi+1,zi+1,1)];
  VecTemp[1]+=wpmm*this->VelocityField[subivId][ptVF(xi+1,yi,zi,1)];
  VecTemp[1]+=wpmp*this->VelocityField[subivId][ptVF(xi+1,yi,zi+1,1)];
  VecTemp[1]+=wppm*this->VelocityField[subivId][ptVF(xi+1,yi+1,zi,1)];
  VecTemp[1]+=wppp*this->VelocityField[subivId][ptVF(xi+1,yi+1,zi+1,1)];
  
  VecTemp[2]= wmmm*this->VelocityField[subivId][ptVF(xi,yi,zi,2)];
  VecTemp[2]+=wmmp*this->VelocityField[subivId][ptVF(xi,yi,zi+1,2)];
  VecTemp[2]+=wmpm*this->VelocityField[subivId][ptVF(xi,yi+1,zi,2)];
  VecTemp[2]+=wmpp*this->VelocityField[subivId][ptVF(xi,yi+1,zi+1,2)];
  VecTemp[2]+=wpmm*this->VelocityField[subivId][ptVF(xi+1,yi,zi,2)];
  VecTemp[2]+=wpmp*this->VelocityField[subivId][ptVF(xi+1,yi,zi+1,2)];
  VecTemp[2]+=wppm*this->VelocityField[subivId][ptVF(xi+1,yi+1,zi,2)];
  VecTemp[2]+=wppp*this->VelocityField[subivId][ptVF(xi+1,yi+1,zi+1,2)];
  }
}

//Perform a bilinear interplolation in a 3D image where NZ==1 and returns the velocity field at non integer coordinates. (here, VecTemp[2]=0)
//++++ The output is in VecTemp ++++ 
template <class VoxelType> void LargeDefGradLagrange<VoxelType>::GetVectorFrom2DVelocityField(int subivId,float x, float y, float z,float VecTemp[3]){
  int xi,yi;
  float xwm,ywm,xwp,ywp;
  float wmm,wmp,wpm,wpp;
  
  xi=static_cast<int>(x);  xwm=1-(x-static_cast<float>(xi));  xwp=x-static_cast<float>(xi);
  yi=static_cast<int>(y);  ywm=1-(y-static_cast<float>(yi));  ywp=y-static_cast<float>(yi);
  
  if ((xi<0.)||(xi>=this->NX-1)||(yi<0.)||(yi>=this->NY-1)){
    VecTemp[0]=0;
    VecTemp[1]=0;
    VecTemp[2]=0;
  }
  else{
    wmm=xwm*ywm;
    wmp=xwm*ywp;
    wpm=xwp*ywm;
    wpp=xwp*ywp;
  
    VecTemp[0]= wmm*this->VelocityField[subivId][ptVF(xi,yi,0,0)];
    VecTemp[0]+=wmp*this->VelocityField[subivId][ptVF(xi,yi+1,0,0)];
    VecTemp[0]+=wpm*this->VelocityField[subivId][ptVF(xi+1,yi,0,0)];
    VecTemp[0]+=wpp*this->VelocityField[subivId][ptVF(xi+1,yi+1,0,0)];
  
    VecTemp[1]= wmm*this->VelocityField[subivId][ptVF(xi,yi,0,1)];
    VecTemp[1]+=wmp*this->VelocityField[subivId][ptVF(xi,yi+1,0,1)];
    VecTemp[1]+=wpm*this->VelocityField[subivId][ptVF(xi+1,yi,0,1)];
    VecTemp[1]+=wpp*this->VelocityField[subivId][ptVF(xi+1,yi+1,0,1)];
  
    VecTemp[2]=0;
  }
}





//Perform a trilinear interplolation and returns the forward mapping at non integer coordinates.
//++++ The output is in VecTemp ++++ 
template <class VoxelType> void LargeDefGradLagrange<VoxelType>::GetCoordFromForwardMapping(int subivId,float x, float y, float z,float VecTemp[3]){
  int xi,yi,zi;
  float xwm,ywm,zwm,xwp,ywp,zwp;
  float wmmm,wmmp,wmpm,wmpp,wpmm,wpmp,wppm,wppp;
  
  xi=static_cast<int>(x);  xwm=1-(x-static_cast<float>(xi));  xwp=x-static_cast<float>(xi);
  yi=static_cast<int>(y);  ywm=1-(y-static_cast<float>(yi));  ywp=y-static_cast<float>(yi);
  zi=static_cast<int>(z);  zwm=1-(z-static_cast<float>(zi));  zwp=z-static_cast<float>(zi);
  
  if (xi<0.) xi=0.0001;
  if (xi>=this->NX-1) xi=this->NX-1.0001;
  if (yi<0.) yi=0.0001;
  if (yi>=this->NY-1) yi=this->NY-1.0001;
  if (zi<0.) zi=0.0001;
  if (zi>=this->NZ-1) zi=this->NZ-1.0001;
  
  if (this->NZ==1){ //2D image
    GetCoordFrom2DForwardMapping(subivId,x,y,z,VecTemp);
    return;
  }
  
  wmmm=xwm*ywm*zwm;
  wmmp=xwm*ywm*zwp;
  wmpm=xwm*ywp*zwm;
  wmpp=xwm*ywp*zwp;
  wpmm=xwp*ywm*zwm;
  wpmp=xwp*ywm*zwp;
  wppm=xwp*ywp*zwm;
  wppp=xwp*ywp*zwp;

  //cout << "sum = " << wmmm+wmmp+wmpm+wmpp+wpmm+wpmp+wppm+wppp << "\n";

  VecTemp[0]= wmmm*this->ForwardMapping[subivId][ptVF(xi,yi,zi,0)];
  VecTemp[0]+=wmmp*this->ForwardMapping[subivId][ptVF(xi,yi,zi+1,0)];
  VecTemp[0]+=wmpm*this->ForwardMapping[subivId][ptVF(xi,yi+1,zi,0)];
  VecTemp[0]+=wmpp*this->ForwardMapping[subivId][ptVF(xi,yi+1,zi+1,0)];
  VecTemp[0]+=wpmm*this->ForwardMapping[subivId][ptVF(xi+1,yi,zi,0)];
  VecTemp[0]+=wpmp*this->ForwardMapping[subivId][ptVF(xi+1,yi,zi+1,0)];
  VecTemp[0]+=wppm*this->ForwardMapping[subivId][ptVF(xi+1,yi+1,zi,0)];
  VecTemp[0]+=wppp*this->ForwardMapping[subivId][ptVF(xi+1,yi+1,zi+1,0)];

  VecTemp[1]= wmmm*this->ForwardMapping[subivId][ptVF(xi,yi,zi,1)];
  VecTemp[1]+=wmmp*this->ForwardMapping[subivId][ptVF(xi,yi,zi+1,1)];
  VecTemp[1]+=wmpm*this->ForwardMapping[subivId][ptVF(xi,yi+1,zi,1)];
  VecTemp[1]+=wmpp*this->ForwardMapping[subivId][ptVF(xi,yi+1,zi+1,1)];
  VecTemp[1]+=wpmm*this->ForwardMapping[subivId][ptVF(xi+1,yi,zi,1)];
  VecTemp[1]+=wpmp*this->ForwardMapping[subivId][ptVF(xi+1,yi,zi+1,1)];
  VecTemp[1]+=wppm*this->ForwardMapping[subivId][ptVF(xi+1,yi+1,zi,1)];
  VecTemp[1]+=wppp*this->ForwardMapping[subivId][ptVF(xi+1,yi+1,zi+1,1)];

  VecTemp[2]= wmmm*this->ForwardMapping[subivId][ptVF(xi,yi,zi,2)];
  VecTemp[2]+=wmmp*this->ForwardMapping[subivId][ptVF(xi,yi,zi+1,2)];
  VecTemp[2]+=wmpm*this->ForwardMapping[subivId][ptVF(xi,yi+1,zi,2)];
  VecTemp[2]+=wmpp*this->ForwardMapping[subivId][ptVF(xi,yi+1,zi+1,2)];
  VecTemp[2]+=wpmm*this->ForwardMapping[subivId][ptVF(xi+1,yi,zi,2)];
  VecTemp[2]+=wpmp*this->ForwardMapping[subivId][ptVF(xi+1,yi,zi+1,2)];
  VecTemp[2]+=wppm*this->ForwardMapping[subivId][ptVF(xi+1,yi+1,zi,2)];
  VecTemp[2]+=wppp*this->ForwardMapping[subivId][ptVF(xi+1,yi+1,zi+1,2)];
}


//Perform a bilinear interplolation in a 3D image where NZ==1 and returns the forward mapping at non integer coordinates.
//++++ The output is in VecTemp ++++ 
template <class VoxelType> void LargeDefGradLagrange<VoxelType>::GetCoordFrom2DForwardMapping(int subivId,float x, float y, float z,float VecTemp[3]){
  int xi,yi;
  float xwm,ywm,xwp,ywp;
  float wmm,wmp,wpm,wpp;
  
  xi=static_cast<int>(x);  xwm=1-(x-static_cast<float>(xi));  xwp=x-static_cast<float>(xi);
  yi=static_cast<int>(y);  ywm=1-(y-static_cast<float>(yi));  ywp=y-static_cast<float>(yi);
  
  if (xi<0.) xi=0.0001;
  if (xi>=this->NX-1) xi=this->NX-1.0001;
  if (yi<0.) yi=0.0001;
  if (yi>=this->NY-1) yi=this->NY-1.0001;
  
  wmm=xwm*ywm;
  wmp=xwm*ywp;
  wpm=xwp*ywm;
  wpp=xwp*ywp;


  VecTemp[0]= wmm*this->ForwardMapping[subivId][ptVF(xi,yi,0,0)];
  VecTemp[0]+=wmp*this->ForwardMapping[subivId][ptVF(xi,yi+1,0,0)];
  VecTemp[0]+=wpm*this->ForwardMapping[subivId][ptVF(xi+1,yi,0,0)];
  VecTemp[0]+=wpp*this->ForwardMapping[subivId][ptVF(xi+1,yi+1,0,0)];

  VecTemp[1]= wmm*this->ForwardMapping[subivId][ptVF(xi,yi,0,1)];
  VecTemp[1]+=wmp*this->ForwardMapping[subivId][ptVF(xi,yi+1,0,1)];
  VecTemp[1]+=wpm*this->ForwardMapping[subivId][ptVF(xi+1,yi,0,1)];
  VecTemp[1]+=wpp*this->ForwardMapping[subivId][ptVF(xi+1,yi+1,0,1)];

  VecTemp[2]= 0;
}


//Perform a trilinear interplolation and returns the backward mapping at non integer coordinates.
//++++ The output is in VecTemp ++++
template <class VoxelType> void LargeDefGradLagrange<VoxelType>::GetCoordFromBackwardMapping(int subivId,float x, float y, float z,float VecTemp[3]){
  int xi,yi,zi;
  float xwm,ywm,zwm,xwp,ywp,zwp;
  float wmmm,wmmp,wmpm,wmpp,wpmm,wpmp,wppm,wppp;
  
  xi=static_cast<int>(x);  xwm=1-(x-static_cast<float>(xi));  xwp=x-static_cast<float>(xi);
  yi=static_cast<int>(y);  ywm=1-(y-static_cast<float>(yi));  ywp=y-static_cast<float>(yi);
  zi=static_cast<int>(z);  zwm=1-(z-static_cast<float>(zi));  zwp=z-static_cast<float>(zi);
  
  if (xi<0.) xi=0.0001;
  if (xi>=this->NX-1) xi=this->NX-1.0001;
  if (yi<0.) yi=0.0001;
  if (yi>=this->NY-1) yi=this->NY-1.0001;
  if (zi<0.) zi=0.0001;
  if (zi>=this->NZ-1) zi=this->NZ-1.0001;
  
  if (this->NZ==1){ //2D image
    GetCoordFrom2DBackwardMapping(subivId,x,y,z,VecTemp);
    return;
  }
  
  wmmm=xwm*ywm*zwm;
  wmmp=xwm*ywm*zwp;
  wmpm=xwm*ywp*zwm;
  wmpp=xwm*ywp*zwp;
  wpmm=xwp*ywm*zwm;
  wpmp=xwp*ywm*zwp;
  wppm=xwp*ywp*zwm;
  wppp=xwp*ywp*zwp;

  //cout << "sum = " << wmmm+wmmp+wmpm+wmpp+wpmm+wpmp+wppm+wppp << "\n";

  VecTemp[0]= wmmm*this->BackwardMapping[subivId][ptVF(xi,yi,zi,0)];
  VecTemp[0]+=wmmp*this->BackwardMapping[subivId][ptVF(xi,yi,zi+1,0)];
  VecTemp[0]+=wmpm*this->BackwardMapping[subivId][ptVF(xi,yi+1,zi,0)];
  VecTemp[0]+=wmpp*this->BackwardMapping[subivId][ptVF(xi,yi+1,zi+1,0)];
  VecTemp[0]+=wpmm*this->BackwardMapping[subivId][ptVF(xi+1,yi,zi,0)];
  VecTemp[0]+=wpmp*this->BackwardMapping[subivId][ptVF(xi+1,yi,zi+1,0)];
  VecTemp[0]+=wppm*this->BackwardMapping[subivId][ptVF(xi+1,yi+1,zi,0)];
  VecTemp[0]+=wppp*this->BackwardMapping[subivId][ptVF(xi+1,yi+1,zi+1,0)];

  VecTemp[1]= wmmm*this->BackwardMapping[subivId][ptVF(xi,yi,zi,1)];
  VecTemp[1]+=wmmp*this->BackwardMapping[subivId][ptVF(xi,yi,zi+1,1)];
  VecTemp[1]+=wmpm*this->BackwardMapping[subivId][ptVF(xi,yi+1,zi,1)];
  VecTemp[1]+=wmpp*this->BackwardMapping[subivId][ptVF(xi,yi+1,zi+1,1)];
  VecTemp[1]+=wpmm*this->BackwardMapping[subivId][ptVF(xi+1,yi,zi,1)];
  VecTemp[1]+=wpmp*this->BackwardMapping[subivId][ptVF(xi+1,yi,zi+1,1)];
  VecTemp[1]+=wppm*this->BackwardMapping[subivId][ptVF(xi+1,yi+1,zi,1)];
  VecTemp[1]+=wppp*this->BackwardMapping[subivId][ptVF(xi+1,yi+1,zi+1,1)];

  VecTemp[2]= wmmm*this->BackwardMapping[subivId][ptVF(xi,yi,zi,2)];
  VecTemp[2]+=wmmp*this->BackwardMapping[subivId][ptVF(xi,yi,zi+1,2)];
  VecTemp[2]+=wmpm*this->BackwardMapping[subivId][ptVF(xi,yi+1,zi,2)];
  VecTemp[2]+=wmpp*this->BackwardMapping[subivId][ptVF(xi,yi+1,zi+1,2)];
  VecTemp[2]+=wpmm*this->BackwardMapping[subivId][ptVF(xi+1,yi,zi,2)];
  VecTemp[2]+=wpmp*this->BackwardMapping[subivId][ptVF(xi+1,yi,zi+1,2)];
  VecTemp[2]+=wppm*this->BackwardMapping[subivId][ptVF(xi+1,yi+1,zi,2)];
  VecTemp[2]+=wppp*this->BackwardMapping[subivId][ptVF(xi+1,yi+1,zi+1,2)];
}


//Perform a bilinear interplolation in a 3D image where NZ==1 and returns the backward mapping at non integer coordinates.
//++++ The output is in VecTemp ++++ 
template <class VoxelType> void LargeDefGradLagrange<VoxelType>::GetCoordFrom2DBackwardMapping(int subivId,float x, float y, float z,float VecTemp[3]){
  int xi,yi;
  float xwm,ywm,xwp,ywp;
  float wmm,wmp,wpm,wpp;
  
  xi=static_cast<int>(x);  xwm=1-(x-static_cast<float>(xi));  xwp=x-static_cast<float>(xi);
  yi=static_cast<int>(y);  ywm=1-(y-static_cast<float>(yi));  ywp=y-static_cast<float>(yi);
  
  if (xi<0.) xi=0.0001;
  if (xi>=this->NX-1) xi=this->NX-1.0001;
  if (yi<0.) yi=0.0001;
  if (yi>=this->NY-1) yi=this->NY-1.0001;
  
  wmm=xwm*ywm;
  wmp=xwm*ywp;
  wpm=xwp*ywm;
  wpp=xwp*ywp;


  VecTemp[0]= wmm*this->BackwardMapping[subivId][ptVF(xi,yi,0,0)];
  VecTemp[0]+=wmp*this->BackwardMapping[subivId][ptVF(xi,yi+1,0,0)];
  VecTemp[0]+=wpm*this->BackwardMapping[subivId][ptVF(xi+1,yi,0,0)];
  VecTemp[0]+=wpp*this->BackwardMapping[subivId][ptVF(xi+1,yi+1,0,0)];

  VecTemp[1]= wmm*this->BackwardMapping[subivId][ptVF(xi,yi,0,1)];
  VecTemp[1]+=wmp*this->BackwardMapping[subivId][ptVF(xi,yi+1,0,1)];
  VecTemp[1]+=wpm*this->BackwardMapping[subivId][ptVF(xi+1,yi,0,1)];
  VecTemp[1]+=wpp*this->BackwardMapping[subivId][ptVF(xi+1,yi+1,0,1)];

  VecTemp[2]= 0;
}




//Perform a trilinear interplolation and returns the grey level of the template image at non integer coordinates.
template <class VoxelType> float LargeDefGradLagrange<VoxelType>::GetGreyLevelFromTemplate(float x, float y, float z){
  float InterpoGreyLevel;
  int xi,yi,zi;
  float xwm,ywm,zwm,xwp,ywp,zwp;
  float wmmm,wmmp,wmpm,wmpp,wpmm,wpmp,wppm,wppp;

  xi=static_cast<int>(x);  xwm=1-(x-static_cast<float>(xi));  xwp=x-static_cast<float>(xi);
  yi=static_cast<int>(y);  ywm=1-(y-static_cast<float>(yi));  ywp=y-static_cast<float>(yi);
  zi=static_cast<int>(z);  zwm=1-(z-static_cast<float>(zi));  zwp=z-static_cast<float>(zi);

  if ((xi<0.)||(xi>=this->NX-1)||(yi<0.)||(yi>=this->NY-1)||(zi<0.)||(zi>=this->NZ-1)){
    InterpoGreyLevel=0;
    if (this->NZ==1) InterpoGreyLevel=GetGreyLevelFrom2DTemplate(x,y,z);
  }
  else{
    wmmm=xwm*ywm*zwm; wmmp=xwm*ywm*zwp; wmpm=xwm*ywp*zwm; wmpp=xwm*ywp*zwp;
    wpmm=xwp*ywm*zwm; wpmp=xwp*ywm*zwp; wppm=xwp*ywp*zwm; wppp=xwp*ywp*zwp;

  //cout << "sum = " << wmmm+wmmp+wmpm+wmpp+wpmm+wpmp+wppm+wppp << "\n";
  
    InterpoGreyLevel= wmmm*this->ImTemplate[ptSF(xi,yi,zi)];
    InterpoGreyLevel+=wmmp*this->ImTemplate[ptSF(xi,yi,zi+1)];
    InterpoGreyLevel+=wmpm*this->ImTemplate[ptSF(xi,yi+1,zi)];
    InterpoGreyLevel+=wmpp*this->ImTemplate[ptSF(xi,yi+1,zi+1)];
    InterpoGreyLevel+=wpmm*this->ImTemplate[ptSF(xi+1,yi,zi)];
    InterpoGreyLevel+=wpmp*this->ImTemplate[ptSF(xi+1,yi,zi+1)];
    InterpoGreyLevel+=wppm*this->ImTemplate[ptSF(xi+1,yi+1,zi)];
    InterpoGreyLevel+=wppp*this->ImTemplate[ptSF(xi+1,yi+1,zi+1)];
  }

  return InterpoGreyLevel;
}


//Perform a bilinear interplolation in a 3D image where NZ==1 and returns
template <class VoxelType> float LargeDefGradLagrange<VoxelType>::GetGreyLevelFrom2DTemplate(float x, float y, float z){
  float InterpoGreyLevel;
  int xi,yi;
  float xwm,ywm,xwp,ywp;
  float wmm,wmp,wpm,wpp;
  
  xi=static_cast<int>(x);  xwm=1-(x-static_cast<float>(xi));  xwp=x-static_cast<float>(xi);
  yi=static_cast<int>(y);  ywm=1-(y-static_cast<float>(yi));  ywp=y-static_cast<float>(yi);
  
  if ((xi<0.)||(xi>=this->NX-1)||(yi<0.)||(yi>=this->NY-1)){
    InterpoGreyLevel=0;
  }
  else{
    wmm=xwm*ywm;
    wmp=xwm*ywp;
    wpm=xwp*ywm;
    wpp=xwp*ywp;
  
    InterpoGreyLevel= wmm*this->ImTemplate[ptSF(xi,yi,0)];
    InterpoGreyLevel+=wmp*this->ImTemplate[ptSF(xi,yi+1,0)];
    InterpoGreyLevel+=wpm*this->ImTemplate[ptSF(xi+1,yi,0)];
    InterpoGreyLevel+=wpp*this->ImTemplate[ptSF(xi+1,yi+1,0)];
  }
  
  return InterpoGreyLevel;
}




//Perform a trilinear interplolation and returns the grey level of the target image at non integer coordinates.
template <class VoxelType> float LargeDefGradLagrange<VoxelType>::GetGreyLevelFromTarget(float x, float y, float z){
  float InterpoGreyLevel;
  int xi,yi,zi;
  float xwm,ywm,zwm,xwp,ywp,zwp;
  float wmmm,wmmp,wmpm,wmpp,wpmm,wpmp,wppm,wppp;

  xi=static_cast<int>(x);  xwm=1-(x-static_cast<float>(xi));  xwp=x-static_cast<float>(xi);
  yi=static_cast<int>(y);  ywm=1-(y-static_cast<float>(yi));  ywp=y-static_cast<float>(yi);
  zi=static_cast<int>(z);  zwm=1-(z-static_cast<float>(zi));  zwp=z-static_cast<float>(zi);

  if ((xi<0.)||(xi>=this->NX-1)||(yi<0.)||(yi>=this->NY-1)||(zi<0.)||(zi>=this->NZ-1)){
    InterpoGreyLevel=0;
    if (this->NZ==1) InterpoGreyLevel=GetGreyLevelFrom2DTarget(x,y,z);
  }
  else{
    wmmm=xwm*ywm*zwm; wmmp=xwm*ywm*zwp; wmpm=xwm*ywp*zwm; wmpp=xwm*ywp*zwp;
    wpmm=xwp*ywm*zwm; wpmp=xwp*ywm*zwp; wppm=xwp*ywp*zwm; wppp=xwp*ywp*zwp;

  //cout << "sum = " << wmmm+wmmp+wmpm+wmpp+wpmm+wpmp+wppm+wppp << "\n";
  
    InterpoGreyLevel= wmmm*this->ImTarget[ptSF(xi,yi,zi)];
    InterpoGreyLevel+=wmmp*this->ImTarget[ptSF(xi,yi,zi+1)];
    InterpoGreyLevel+=wmpm*this->ImTarget[ptSF(xi,yi+1,zi)];
    InterpoGreyLevel+=wmpp*this->ImTarget[ptSF(xi,yi+1,zi+1)];
    InterpoGreyLevel+=wpmm*this->ImTarget[ptSF(xi+1,yi,zi)];
    InterpoGreyLevel+=wpmp*this->ImTarget[ptSF(xi+1,yi,zi+1)];
    InterpoGreyLevel+=wppm*this->ImTarget[ptSF(xi+1,yi+1,zi)];
    InterpoGreyLevel+=wppp*this->ImTarget[ptSF(xi+1,yi+1,zi+1)];
  }

  return InterpoGreyLevel;
}

//Perform a bilinear interplolation in a 3D image where NZ==1 and returns the grey level
template <class VoxelType> float LargeDefGradLagrange<VoxelType>::GetGreyLevelFrom2DTarget(float x, float y, float z){
  float InterpoGreyLevel;
  int xi,yi;
  float xwm,ywm,xwp,ywp;
  float wmm,wmp,wpm,wpp;
  
  xi=static_cast<int>(x);  xwm=1-(x-static_cast<float>(xi));  xwp=x-static_cast<float>(xi);
  yi=static_cast<int>(y);  ywm=1-(y-static_cast<float>(yi));  ywp=y-static_cast<float>(yi);
  
  if ((xi<0.)||(xi>=this->NX-1)||(yi<0.)||(yi>=this->NY-1)){
    InterpoGreyLevel=0;
  }
  else{
    wmm=xwm*ywm;
    wmp=xwm*ywp;
    wpm=xwp*ywm;
    wpp=xwp*ywp;
  
    InterpoGreyLevel= wmm*this->ImTarget[ptSF(xi,yi,0)];
    InterpoGreyLevel+=wmp*this->ImTarget[ptSF(xi,yi+1,0)];
    InterpoGreyLevel+=wpm*this->ImTarget[ptSF(xi+1,yi,0)];
    InterpoGreyLevel+=wpp*this->ImTarget[ptSF(xi+1,yi+1,0)];
  }
  
  return InterpoGreyLevel;
}



//Perform a trilinear interplolation and returns the grey level of the float generic image at non integer coordinates.
template <class VoxelType> float LargeDefGradLagrange<VoxelType>::GetGreyLevelFromFloatGenericImage(irtkGenericImage<float>* Temp4DField, float x, float y, float z, float t){
  float InterpoGreyLevel;
  int xi,yi,zi,ti;
  float xwm,ywm,zwm,twm,xwp,ywp,zwp,twp;
  float wmmmm,wmmpm,wmpmm,wmppm,wpmmm,wpmpm,wppmm,wpppm;
  float wmmmp,wmmpp,wmpmp,wmppp,wpmmp,wpmpp,wppmp,wpppp;
  
  xi=static_cast<int>(x);  xwm=1-(x-static_cast<float>(xi));  xwp=x-static_cast<float>(xi);
  yi=static_cast<int>(y);  ywm=1-(y-static_cast<float>(yi));  ywp=y-static_cast<float>(yi);
  zi=static_cast<int>(z);  zwm=1-(z-static_cast<float>(zi));  zwp=z-static_cast<float>(zi);
  ti=static_cast<int>(t);  twm=1-(t-static_cast<float>(ti));  twp=t-static_cast<float>(ti);

  if ((xi<0.)||(xi>=Temp4DField->GetX()-1)||(yi<0.)||(yi>=Temp4DField->GetY()-1)||(zi<0.)||(zi>=Temp4DField->GetZ()-1)||(ti<0.)||(ti>=Temp4DField->GetT()-1)){
    InterpoGreyLevel=0;
    if (ti==Temp4DField->GetT()-1)InterpoGreyLevel=Temp4DField->Get(xi,yi,zi,ti);
    if (Temp4DField->GetZ()==1) InterpoGreyLevel=GetGreyLevelFrom2DFloatGenericImage(Temp4DField,x, y, z, t);
  }
  else{
    wmmmm=xwm*ywm*zwm*twm; wmmpm=xwm*ywm*zwp*twm; wmpmm=xwm*ywp*zwm*twm; wmppm=xwm*ywp*zwp*twm;
    wpmmm=xwp*ywm*zwm*twm; wpmpm=xwp*ywm*zwp*twm; wppmm=xwp*ywp*zwm*twm; wpppm=xwp*ywp*zwp*twm;
    wmmmp=xwm*ywm*zwm*twp; wmmpp=xwm*ywm*zwp*twp; wmpmp=xwm*ywp*zwm*twp; wmppp=xwm*ywp*zwp*twp;
    wpmmp=xwp*ywm*zwm*twp; wpmpp=xwp*ywm*zwp*twp; wppmp=xwp*ywp*zwm*twp; wpppp=xwp*ywp*zwp*twp;
    
    //cout << "sum m = " << wmmmm+wmmpm+wmpmm+wmppm+wpmmm+wpmpm+wppmm+wpppm << "\n";
    //cout << "sum p = " << wmmmp+wmmpp+wmpmp+wmppp+wpmmp+wpmpp+wppmp+wpppp << "\n";
    //cout << xi << " " << yi << " " << zi << " " << ti << " " << Temp4DField->Get(xi,yi,zi,ti) << "\n";
    
    InterpoGreyLevel= wmmmm*static_cast<float>(Temp4DField->Get(xi,yi,zi,ti));
    InterpoGreyLevel+=wmmpm*static_cast<float>(Temp4DField->Get(xi,yi,zi+1,ti));
    InterpoGreyLevel+=wmpmm*static_cast<float>(Temp4DField->Get(xi,yi+1,zi,ti));
    InterpoGreyLevel+=wmppm*static_cast<float>(Temp4DField->Get(xi,yi+1,zi+1,ti));
    InterpoGreyLevel+=wpmmm*static_cast<float>(Temp4DField->Get(xi+1,yi,zi,ti));
    InterpoGreyLevel+=wpmpm*static_cast<float>(Temp4DField->Get(xi+1,yi,zi+1,ti));
    InterpoGreyLevel+=wppmm*static_cast<float>(Temp4DField->Get(xi+1,yi+1,zi,ti));
    InterpoGreyLevel+=wpppm*static_cast<float>(Temp4DField->Get(xi+1,yi+1,zi+1,ti));
    InterpoGreyLevel+=wmmmp*static_cast<float>(Temp4DField->Get(xi,yi,zi,ti+1));
    InterpoGreyLevel+=wmmpp*static_cast<float>(Temp4DField->Get(xi,yi,zi+1,ti+1));
    InterpoGreyLevel+=wmpmp*static_cast<float>(Temp4DField->Get(xi,yi+1,zi,ti+1));
    InterpoGreyLevel+=wmppp*static_cast<float>(Temp4DField->Get(xi,yi+1,zi+1,ti+1));
    InterpoGreyLevel+=wpmmp*static_cast<float>(Temp4DField->Get(xi+1,yi,zi,ti+1));
    InterpoGreyLevel+=wpmpp*static_cast<float>(Temp4DField->Get(xi+1,yi,zi+1,ti+1));
    InterpoGreyLevel+=wppmp*static_cast<float>(Temp4DField->Get(xi+1,yi+1,zi,ti+1));
    InterpoGreyLevel+=wpppp*static_cast<float>(Temp4DField->Get(xi+1,yi+1,zi+1,ti+1));
  }
  
  return InterpoGreyLevel;
}


//Perform a trilinear interplolation and returns the grey level of the float generic image at non integer coordinates.
template <class VoxelType> float LargeDefGradLagrange<VoxelType>::GetGreyLevelFrom2DFloatGenericImage(irtkGenericImage<float>* Temp4DField, float x, float y, float z, float t){
  float InterpoGreyLevel;
  int xi,yi,ti;
  float xwm,ywm,twm,xwp,ywp,twp;
  float wmmm,wmpm,wpmm,wppm;
  float wmmp,wmpp,wpmp,wppp;
  
  xi=static_cast<int>(x);  xwm=1-(x-static_cast<float>(xi));  xwp=x-static_cast<float>(xi);
  yi=static_cast<int>(y);  ywm=1-(y-static_cast<float>(yi));  ywp=y-static_cast<float>(yi);
  ti=static_cast<int>(t);  twm=1-(t-static_cast<float>(ti));  twp=t-static_cast<float>(ti);

  if ((xi<0.)||(xi>=Temp4DField->GetX()-1)||(yi<0.)||(yi>=Temp4DField->GetY()-1)||(ti<0.)||(ti>=Temp4DField->GetT()-1)){
    InterpoGreyLevel=0;
    if (ti==Temp4DField->GetT()-1)InterpoGreyLevel=Temp4DField->Get(xi,yi,0,ti);
  }
  else{
    wmmm=xwm*ywm*twm; wmpm=xwm*ywp*twm;
    wpmm=xwp*ywm*twm; wppm=xwp*ywp*twm;
    wmmp=xwm*ywm*twp; wmpp=xwm*ywp*twp;
    wpmp=xwp*ywm*twp; wppp=xwp*ywp*twp;
    
    
    InterpoGreyLevel= wmmm*static_cast<float>(Temp4DField->Get(xi,yi,0,ti));
    InterpoGreyLevel+=wmpm*static_cast<float>(Temp4DField->Get(xi,yi+1,0,ti));
    InterpoGreyLevel+=wpmm*static_cast<float>(Temp4DField->Get(xi+1,yi,0,ti));
    InterpoGreyLevel+=wppm*static_cast<float>(Temp4DField->Get(xi+1,yi+1,0,ti));
    InterpoGreyLevel+=wmmp*static_cast<float>(Temp4DField->Get(xi,yi,0,ti+1));
    InterpoGreyLevel+=wmpp*static_cast<float>(Temp4DField->Get(xi,yi+1,0,ti+1));
    InterpoGreyLevel+=wpmp*static_cast<float>(Temp4DField->Get(xi+1,yi,0,ti+1));
    InterpoGreyLevel+=wppp*static_cast<float>(Temp4DField->Get(xi+1,yi+1,0,ti+1));
  }
  
  return InterpoGreyLevel;
}



///++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
///                                 FUNCTIONS TO PERFORM THE REGISTRATION
///++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++



///Function to solve the registration using the gradient descent algorithm of Beg 05
template <class VoxelType> void LargeDefGradLagrange<VoxelType>::Run_Default(void){  
  int t;
  int IterationsNb;
  int TimeSubdiv;
  int ReparamCount;
  
  //1) INITIALISATION
  
  //1.1) Initial set up
  this->Initialize();
  
  //1.2) allocations...
  this->AllocateAllVariables();
  
  //2) LOOP ON THE TIME FRAMES
  for (t = 0; t < this->NT; t++) {
    cout << "Image " << t+1 << " / " << this->NT << "\n";
    
    //2.1 ) INITIALISATION
    this->InitiateGradientDescent(t);
    
    //2.2) LOOP OF THE GRADIENT DESCENT
    IterationsNb=0;
    ReparamCount=0;
    while (IterationsNb<this->iteration_nb){
      cout << "Iteration Number " << IterationsNb+1 << " / " << this->iteration_nb << "\n";
      
      //2.2.1) compute the forward mapping on space
      //cout << "Compute the spatial forward Mapping\n";
      this->ComputeForwardMapping();
      
      //2.2.2) compute the backward mapping on space
      //cout << "Compute the spatial Backward Mapping\n";
      this->ComputeBackwardMapping();
      
      //2.2.4) LOOP ON THE TIME SUBDIVISIONS BETWEEN TIME=0 (template) AND TIME=1 (target)
      for (TimeSubdiv=0;TimeSubdiv<this->NbTimeSubdiv;TimeSubdiv++){
        //cout << "Subdivision Number " << TimeSubdiv+1 << " / " << this->NbTimeSubdiv << "\n";
        
        //2.2.4.0) manage the current filter in case of multikernel registration
        if (this->NbKernels>1) KernelManager(TimeSubdiv);
        
        //2.2.4.1) compute the temporary image transformed using the forward mapping from time 0 -> J0
        //cout << "Computing J0 at time step t=" << TimeSubdiv << "\n";
        this->ComputeJ0(TimeSubdiv);
        
        //2.2.4.3) compute the temporary image transformed using the backward mapping from time 1 -> J1
        //cout << "Computing J1 at time step t=" << TimeSubdiv << "\n";
        this->ComputeJ1(TimeSubdiv);
        
        //2.2.4.4) compute gradient of J0
        //cout << "Compute Grad J0\n";
        this->ComputeGradientJ0();
        
        //2.2.4.5) compute the determinant of the jacobian of the transformation
        //cout << "Computing Det Jacobian\n";
        this->ComputeJacobianDeterminant(TimeSubdiv);
        
        //2.2.4.6) compute the gradient of energy
        //cout << "ComputeEnergyGradient\n";
        this->ComputeEnergyGradient(TimeSubdiv);
      }
      
      //2.2.5) update the velocity fields
      //cout << "Update Vector Field\n";
      this->UpdateVelocityField(IterationsNb);

      
      /*BEGIN TO ADD AGAIN   -   BEGIN TO ADD AGAIN   -   BEGIN TO ADD AGAIN   -   BEGIN TO ADD AGAIN   -   BEGIN TO ADD AGAIN*/
/*      
      //2.2.6) Compute the norms and length of the velocity field plus the total energy
      //cout << "Compute Norms, Length and Energy\n";
      this->ComputeNormsAndLength();
      this->ComputeCostFunction();
      
      //2.2.7) update convergence test parameters
      cout << "\nIteration " << IterationsNb << ":\n";
      cout << "Norm = ";
      for (TimeSubdiv=0;TimeSubdiv<this->NbTimeSubdiv;TimeSubdiv++) cout << this->norm[TimeSubdiv] << " ";
      cout << "\n";
      cout << "Length = " << this->length << "\n";
      cout << "Energy / Velocity Field = " << this->EnergyVF << "\n";
      cout << "Energy / Difference Images = " << this->EnergyDI << "\n";
      cout << "Energy / Total  = " << this->EnergyTot << "\n \n";
      */
      
      /*END TO ADD AGAIN   -   END TO ADD AGAIN   -   END TO ADD AGAIN   -   END TO ADD AGAIN   -   END TO ADD AGAIN*/
      
      IterationsNb++;
      
      //2.2.8) velocity field reparametrization
      ReparamCount++;
      if (ReparamCount==this->reparametrization){
        cout << "Speed reparametrization\n";
        this->ComputeNormsAndLength();   //to remove when TO ADD AGAIN will be added again
        this->ComputeCostFunction();     //to remove when TO ADD AGAIN will be added again
        this->SpeedReparameterization();
        ReparamCount=0;
      }
    }
    
    //2.3) save the filtered temporary 3D image in VoxelType in the  output image at time t
    //cout << "Save the result\n";
    this->SaveResultGradientDescent(t);
  }
  
  //3) END OF THE FUNCTION
  // Do the final cleaning up
  this->Finalize();
}





///++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
///                                        MAIN RUN FUNCTION 
///++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


///run function
template <class VoxelType> void LargeDefGradLagrange<VoxelType>::Run(void)
{
  this->Run_Default();
}





///++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
///                                        TEMPLATE TYPES
///++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


template class LargeDefGradLagrange<irtkBytePixel>;
template class LargeDefGradLagrange<irtkGreyPixel>;
template class LargeDefGradLagrange<irtkRealPixel>;


