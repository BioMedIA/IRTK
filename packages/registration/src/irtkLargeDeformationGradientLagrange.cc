/*=========================================================================

  Date      : $Date: 29.06.2009$
  Changes   : $Authors: Laurent Risser, Francois-Xavier Vialard$

=========================================================================*/

#include <irtkImage.h>
#include <irtkLargeDeformationGradientLagrange.h>

///++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
///                                       GENERAL CLASS FUNCTIONS
///++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
template <class VoxelType> LargeDefGradLagrange<VoxelType>::LargeDefGradLagrange(){
  //default parameters
  epsilon=0.1;
  iteration_nb=10;
  NbTimeSubdiv=10;
  MaxVelocityUpdate=3.;
  sigma=3.;
  Margin=0;
  UNDETERMINED_VALUE=-300000000; //has to be smaller than all computed determinants of jacobians and momentum (that means possibly very small)
  DeltaVox=1.;
  alpha=0.01;
  gamma=0.1;
  reparametrization = 10;
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
///                               SUB-FUNCTIONS TO PERFORM THE REGISTRATION  (level 1)
///++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

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
  
  
  //... velocity field
  //    -->  VelocityField[i][ptVF(x,y,z,0)]= direction ex of the vector at (x,y,z)
  //    -->  VelocityField[i][ptVF(x,y,z,1)]= direction ey of the vector at (x,y,z)
  //    -->  VelocityField[i][ptVF(x,y,z,2)]= direction ez of the vector at (x,y,z)
  this->VelocityField= new float* [this->NbTimeSubdiv];
  for (i=0;i<this->NbTimeSubdiv;i++) this->VelocityField[i]=new float [this->NXtYtZt3];
  
  //... direct mapping
  //    -->  DirectMapping[i][ptVF(x,y,z,0)]= coordinate x at time i corresponding to (x,y,z) at virtual time 0
  //    -->  DirectMapping[i][ptVF(x,y,z,1)]= coordinate y at time i corresponding to (x,y,z) at virtual time 0
  //    -->  DirectMapping[i][ptVF(x,y,z,2)]= coordinate z at time i corresponding to (x,y,z) at virtual time 0
  this->DirectMapping= new float* [this->NbTimeSubdiv];
  for (i=0;i<this->NbTimeSubdiv;i++) this->DirectMapping[i]=new float [this->NXtYtZt3];
  
  //... inverse mapping
  //    -->  InverseMapping[i][ptVF(x,y,z,0)]= coordinate x at time i corresponding to (x,y,z) at virtual time 1
  //    -->  InverseMapping[i][ptVF(x,y,z,1)]= coordinate y at time i corresponding to (x,y,z) at virtual time 1
  //    -->  InverseMapping[i][ptVF(x,y,z,2)]= coordinate z at time i corresponding to (x,y,z) at virtual time 1
  this->InverseMapping= new float* [this->NbTimeSubdiv];
  for (i=0;i<this->NbTimeSubdiv;i++) this->InverseMapping[i]=new float [this->NXtYtZt3];
  
  //... temporary image transformed using the direct mapping from time 0
  //    -->  J0[ptSF(x,y,z)]= gray level of the transformed image J0 at (x,y,z)
  this->J0 = new float [this->NXtYtZ];
  
  //... temporary image transformed using the inverse mapping from time 1
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
}


///initiate the gradient descent (Beg 2005) for the current 3D image of the 4D time sequence
template <class VoxelType> void LargeDefGradLagrange<VoxelType>::InitiateGradientDescent(int TimeLoc){
  int subivId, x, y, z;
  float MinGrayLevelImTemplate,MaxGrayLevelImTemplate,MinGrayLevelImTarget,MaxGrayLevelImTarget;
  int DistClosestEdge;
  
  //1) cast the values of the template and target images at time t in float
  for (z = 0; z < this->NZ; z++)  for (y = 0; y < this->NY; y++) for (x = 0; x < this->NX; x++)
        this->ImTemplate[ptSF(x,y,z)]=static_cast<float>(this->_input->Get(x, y, z, TimeLoc));
  
  for (z = 0; z < this->NZ; z++)  for (y = 0; y < this->NY; y++) for (x = 0; x < this->NX; x++)
        this->ImTarget[ptSF(x,y,z)]=static_cast<float>(this->target_image.Get(x, y, z, TimeLoc));
  
  //2) take into accout the margin and normalize the gray levels of ImTemplate and ImTarget between 0 and 100
  //2.1) Margins
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
  
  //2.2)Gray level regularisation
  MinGrayLevelImTemplate=this->ImTemplate[ptSF(1,1,1)];
  for (z = 1; z < this->NZ-1; z++)  for (y = 1; y < this->NY-1; y++) for (x = 1; x < this->NX-1; x++)
        if (this->ImTemplate[ptSF(x,y,z)]<MinGrayLevelImTemplate) MinGrayLevelImTemplate=this->ImTemplate[ptSF(x,y,z)];
  
  MaxGrayLevelImTemplate=this->ImTemplate[ptSF(1,1,1)];
  for (z = 1; z < this->NZ-1; z++)  for (y = 1; y < this->NY-1; y++) for (x = 1; x < this->NX-1; x++)
        if (this->ImTemplate[ptSF(x,y,z)]>MaxGrayLevelImTemplate) MaxGrayLevelImTemplate=this->ImTemplate[ptSF(x,y,z)];

  for (z = 0; z < this->NZ; z++)  for (y = 0; y < this->NY; y++) for (x = 0; x < this->NX; x++)
        this->ImTemplate[ptSF(x,y,z)]=100.*(this->ImTemplate[ptSF(x,y,z)]-MinGrayLevelImTemplate)/(MaxGrayLevelImTemplate-MinGrayLevelImTemplate);
  
  MinGrayLevelImTarget=this->ImTarget[ptSF(1,1,1)];
  for (z = 1; z < this->NZ-1; z++)  for (y = 1; y < this->NY-1; y++) for (x = 1; x < this->NX-1; x++)
        if (this->ImTarget[ptSF(x,y,z)]<MinGrayLevelImTarget) MinGrayLevelImTarget=this->ImTarget[ptSF(x,y,z)];
  
  MaxGrayLevelImTarget=this->ImTarget[ptSF(1,1,1)];
  for (z = 1; z < this->NZ-1; z++)  for (y = 1; y < this->NY-1; y++) for (x = 1; x < this->NX-1; x++)
        if (this->ImTarget[ptSF(x,y,z)]>MaxGrayLevelImTarget) MaxGrayLevelImTarget=this->ImTarget[ptSF(x,y,z)];

  for (z = 0; z < this->NZ; z++)  for (y = 0; y < this->NY; y++) for (x = 0; x < this->NX; x++)
        this->ImTarget[ptSF(x,y,z)]=100.*(this->ImTarget[ptSF(x,y,z)]-MinGrayLevelImTarget)/(MaxGrayLevelImTarget-MinGrayLevelImTarget);

  //3) initiate the vector field
  if (strcmp(PrefixInputVF,"Null")!=0){
    this->LoadVelocityFields(PrefixInputVF);
  }
  else{
    for (subivId=0;subivId<this->NbTimeSubdiv;subivId++) for (z=0;z<this->NZ;z++) for (y=0;y<this->NY;y++) for (x=0;x<this->NX;x++){
      this->VelocityField[subivId][ptVF(x,y,z,0)]=0; //(static_cast<float>(pow(static_cast<double>(x-y)/20,2.)))*(static_cast<float>(subivId)/5);
      this->VelocityField[subivId][ptVF(x,y,z,1)]=0; //(static_cast<float>(pow(static_cast<double>(x+y-2*64)/20,2.)))*(static_cast<float>(subivId)/5);
      this->VelocityField[subivId][ptVF(x,y,z,2)]=0;
    }
  }
  
  //4) initiate the values of the direct mapping at subdivision time 0
  for (z=0;z<this->NZ;z++) for (y=0;y<this->NY;y++) for (x=0;x<this->NX;x++){
    this->DirectMapping[0][ptVF(x,y,z,0)]=0.;
    this->DirectMapping[0][ptVF(x,y,z,1)]=0.;
    this->DirectMapping[0][ptVF(x,y,z,2)]=0.;
  }
  
  //5) initiate the values of the inverse mapping at subdivision time 1
  for (z=0;z<this->NZ;z++) for (y=0;y<this->NY;y++) for (x=0;x<this->NX;x++){
    this->InverseMapping[this->NbTimeSubdiv-1][ptVF(x,y,z,0)]=0.;
    this->InverseMapping[this->NbTimeSubdiv-1][ptVF(x,y,z,1)]=0.;
    this->InverseMapping[this->NbTimeSubdiv-1][ptVF(x,y,z,2)]=0.;
  }
  
  //6) initiate GradE
  for (subivId=0;subivId<this->NbTimeSubdiv;subivId++) for (z=0;z<this->NZ;z++) for (y=0;y<this->NY;y++) for (x=0;x<this->NX;x++){
    this->GradE[subivId][ptVF(x,y,z,0)]=0;
    this->GradE[subivId][ptVF(x,y,z,1)]=0;
    this->GradE[subivId][ptVF(x,y,z,2)]=0;
  }
  
  //7) initiate norm, Length and Energy
  for (subivId=0;subivId<this->NbTimeSubdiv;subivId++) this->norm[subivId]=0;
  this->length=0;
  this->EnergyTot=0;
  this->EnergyVF=0;
  this->EnergyDI=0;
}


///compute the direct mapping
template <class VoxelType> void LargeDefGradLagrange<VoxelType>::ComputeDirectMapping(void){
  float VecTemp[3];  //velocity vector that targets the current treated voxel
  float VecTemp2[3]; //temporary vector to propagate the mapping
  int ConvergenceSteps;
  int i,x,y,z;
  int timeSubdiv;
      
  //initialisation
  ConvergenceSteps=3;
  
  //JO at the first time subdivision
  for (z = 0; z < this->NZ; z++) for (y = 0; y < this->NY; y++) for (x = 0; x < this->NX; x++){
    this->DirectMapping[0][ptVF(x,y,z,0)]=static_cast<float>(x);
    this->DirectMapping[0][ptVF(x,y,z,1)]=static_cast<float>(y);
    this->DirectMapping[0][ptVF(x,y,z,2)]=static_cast<float>(z);
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
        VecTemp[0]=(VecTemp[0]+this->VelocityField[timeSubdiv][ptVF(x,y,z,0)])*(this->DeltaTimeSubdiv/this->DeltaVox)/2;
        VecTemp[1]=(VecTemp[1]+this->VelocityField[timeSubdiv][ptVF(x,y,z,1)])*(this->DeltaTimeSubdiv/this->DeltaVox)/2;
        VecTemp[2]=(VecTemp[2]+this->VelocityField[timeSubdiv][ptVF(x,y,z,2)])*(this->DeltaTimeSubdiv/this->DeltaVox)/2;
      }
      
      //find the original coordinates
      this->GetCoordFromDirectMapping(timeSubdiv-1,x-VecTemp[0],y-VecTemp[1],z-VecTemp[2],VecTemp2);
      
      this->DirectMapping[timeSubdiv][ptVF(x,y,z,0)]=VecTemp2[0];
      this->DirectMapping[timeSubdiv][ptVF(x,y,z,1)]=VecTemp2[1];
      this->DirectMapping[timeSubdiv][ptVF(x,y,z,2)]=VecTemp2[2];
    }
  }
}


///compute the inverse mapping
template <class VoxelType> void LargeDefGradLagrange<VoxelType>::ComputeInverseMapping(void){
  float VecTemp[3];  //velocity vector that targets the current treated voxel
  float VecTemp2[3]; //temporary vector to propagate the mapping
  int ConvergenceSteps;
  int i,x,y,z;
  int timeSubdiv;
  int LastTimSub;

  //initialisation
  ConvergenceSteps=3;
  LastTimSub=this->NbTimeSubdiv-1;
  
  //JO at the first time subdivision
  for (z = 0; z < this->NZ; z++) for (y = 0; y < this->NY; y++) for (x = 0; x < this->NX; x++){
    this->InverseMapping[LastTimSub][ptVF(x,y,z,0)]=static_cast<float>(x);
    this->InverseMapping[LastTimSub][ptVF(x,y,z,1)]=static_cast<float>(y);
    this->InverseMapping[LastTimSub][ptVF(x,y,z,2)]=static_cast<float>(z);
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
        VecTemp[0]=(VecTemp[0]+this->VelocityField[timeSubdiv][ptVF(x,y,z,0)])*(this->DeltaTimeSubdiv/this->DeltaVox)/2;
        VecTemp[1]=(VecTemp[1]+this->VelocityField[timeSubdiv][ptVF(x,y,z,1)])*(this->DeltaTimeSubdiv/this->DeltaVox)/2;
        VecTemp[2]=(VecTemp[2]+this->VelocityField[timeSubdiv][ptVF(x,y,z,2)])*(this->DeltaTimeSubdiv/this->DeltaVox)/2;
      }
      
      //new value
      this->GetCoordFromInverseMapping(timeSubdiv+1,x+VecTemp[0],y+VecTemp[1],z+VecTemp[2],VecTemp2);
      
      this->InverseMapping[timeSubdiv][ptVF(x,y,z,0)]=VecTemp2[0];
      this->InverseMapping[timeSubdiv][ptVF(x,y,z,1)]=VecTemp2[1];
      this->InverseMapping[timeSubdiv][ptVF(x,y,z,2)]=VecTemp2[2];

    }
  }
}


///compute J0 using the direct mapping
template <class VoxelType> void LargeDefGradLagrange<VoxelType>::ComputeJ0(int timeSubdiv){
  int x,y,z;
  
  for (z = 0; z < this->NZ; z++) for (y = 0; y < this->NY; y++) for (x = 0; x < this->NX; x++){
    this->J0[ptSF(x,y,z)]=this->GetGreyLevelFromTemplate(this->DirectMapping[timeSubdiv][ptVF(x,y,z,0)], this->DirectMapping[timeSubdiv][ptVF(x,y,z,1)], this->DirectMapping[timeSubdiv][ptVF(x,y,z,2)]);
  }
}


///compute J1 using the inverse mapping
template <class VoxelType> void LargeDefGradLagrange<VoxelType>::ComputeJ1(int timeSubdiv){
  int x,y,z;
  
  for (z = 0; z < this->NZ; z++) for (y = 0; y < this->NY; y++) for (x = 0; x < this->NX; x++){
    this->J1[ptSF(x,y,z)]=this->GetGreyLevelFromTarget(this->InverseMapping[timeSubdiv][ptVF(x,y,z,0)],this->InverseMapping[timeSubdiv][ptVF(x,y,z,1)],this->InverseMapping[timeSubdiv][ptVF(x,y,z,2)]);
  }
}



///compute the gradient of J0  (may be improved at the boundaries)
template <class VoxelType> void LargeDefGradLagrange<VoxelType>::ComputeGradientJ0(void){
  int x,y,z;
  
  //Energy gradient in direction x, y, z
  for (z = 1; z < this->NZ-2; z++) for (y = 1; y < this->NY-2; y++) for (x = 1; x < this->NX-2; x++){
    if ((this->J0[ptSF(x+1,y,z)]>this->UNDETERMINED_VALUE)&&(this->J0[ptSF(x-1,y,z)]>this->UNDETERMINED_VALUE)&&
      (this->J0[ptSF(x,y+1,z)]>this->UNDETERMINED_VALUE)&&(this->J0[ptSF(x,y-1,z)]>this->UNDETERMINED_VALUE)&&
        (this->J0[ptSF(x,y,z+1)]>this->UNDETERMINED_VALUE)&&(this->J0[ptSF(x,y,z-1)]>this->UNDETERMINED_VALUE)){
          this->GradJ0[ptVF(x,y,z,0)]=(this->J0[ptSF(x+1,y,z)]-this->J0[ptSF(x-1,y,z)])/(2.*this->DeltaVox);
          this->GradJ0[ptVF(x,y,z,1)]=(this->J0[ptSF(x,y+1,z)]-this->J0[ptSF(x,y-1,z)])/(2.*this->DeltaVox);
          this->GradJ0[ptVF(x,y,z,2)]=(this->J0[ptSF(x,y,z+1)]-this->J0[ptSF(x,y,z-1)])/(2.*this->DeltaVox);
        }
        else{
          this->GradJ0[ptVF(x,y,z,0)]=this->UNDETERMINED_VALUE;
          this->GradJ0[ptVF(x,y,z,1)]=this->UNDETERMINED_VALUE;
          this->GradJ0[ptVF(x,y,z,2)]=this->UNDETERMINED_VALUE;
        }
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
  int test;
  
  for (z = 1; z < this->NZ-1; z++) for (y = 1; y < this->NY-1; y++) for (x = 1; x < this->NX-1; x++){
    test=0;
    for(i=0;i<3;i++)for(j=0;j<3;j++){
      this->JacobianMatrix[3*i+j] = (this->InverseMapping[timeSubdiv][ptVF(x + (int)(0==i),y + (int)(1==i),z + (int)(2==i),j)] - this->InverseMapping[timeSubdiv][ptVF(x - (int)(0==i),y - (int)(1==i),z - (int)(2==i),j)])/(2.*this->DeltaVox);
      if ((this->InverseMapping[timeSubdiv][ptVF(x + (int)(0==i),y + (int)(1==i),z + (int)(2==i),j)]<this->UNDETERMINED_VALUE+1)||(this->InverseMapping[timeSubdiv][ptVF(x - (int)(0==i),y - (int)(1==i),z - (int)(2==i),j)]<this->UNDETERMINED_VALUE+1))
        test=1;
    }
    if (test==0){
      this->DetJacobians[ptSF(x,y,z)] = this->CptDetJac3d(this->JacobianMatrix);
      if (this->DetJacobians[ptSF(x,y,z)]<this->UNDETERMINED_VALUE) cout << "WARNING: TOO HIGH -(DETERMINANT OF THE JACOBIAN) IN " << x << " " <<  y << " " << z << "\n";
    }
    else
      this->DetJacobians[ptSF(x,y,z)] = this->UNDETERMINED_VALUE;
  }
  
  //boundary conditions: be aware that the boundary of the image does not move!
  for (y = 0; y < this->NY; y++) for (x = 0; x < this->NX; x++) this->DetJacobians[ptSF(x,y,this->NZ-1)] = 1/(this->DeltaVox*this->DeltaVox*this->DeltaVox);
  for (y = 0; y < this->NY; y++) for (x = 0; x < this->NX; x++) this->DetJacobians[ptSF(x,y,0)] = 1/(this->DeltaVox*this->DeltaVox*this->DeltaVox);	
  for (z = 1; z < this->NZ-1; z++) for (y = 1; y < this->NY-1; y++) this->DetJacobians[ptSF(0,y,z)] = 1/(this->DeltaVox*this->DeltaVox*this->DeltaVox);
  for (z = 1; z < this->NZ-1; z++) for (y = 1; y < this->NY-1; y++) this->DetJacobians[ptSF(this->NX - 1,y,z)] = 1/(this->DeltaVox*this->DeltaVox*this->DeltaVox);
  for (z = 1; z < this->NZ-1; z++) for (x = 1; x < this->NX-1; x++) this->DetJacobians[ptSF(x,0,z)] = 1/(this->DeltaVox*this->DeltaVox*this->DeltaVox);
  for (z = 1; z < this->NZ-1; z++) for (x = 1; x < this->NX-1; x++) this->DetJacobians[ptSF(x,this->NY-1,z)] = 1/(this->DeltaVox*this->DeltaVox*this->DeltaVox);
}



///Compute the energy gradients
template <class VoxelType> void LargeDefGradLagrange<VoxelType>::ComputeEnergyGradient(int timeSubdiv){
  int x,y,z,i;
  float temp;
  irtkRealPixel paddingValue;
  
  //initialisation
  paddingValue=this->UNDETERMINED_VALUE;
  irtkGaussianBlurringWithPadding<float> gaussianBlurring(this->sigma, paddingValue);

  
  //loop on the vector directions (x,y,z)
  for (i=0;i<3;i++){
    //compute the scalar field (one dimension out of the vector field) to smooth
    for (z = 0; z < this->NZ; z++) for (y = 0; y < this->NY; y++) for (x = 0; x < this->NX; x++){
      if ((this->J0[ptSF(x,y,z)]>paddingValue)&&(this->J1[ptSF(x,y,z)]>paddingValue)&&(this->DetJacobians[ptSF(x,y,z)]>paddingValue)&&(this->GradJ0[ptVF(x,y,z,i)]>paddingValue)){
        temp=(this->J0[ptSF(x,y,z)] - this->J1[ptSF(x,y,z)]) * this->DetJacobians[ptSF(x,y,z)] * this->GradJ0[ptVF(x,y,z,i)];
        this->Image3DTemp.Put(x, y, z, 0, static_cast<float>(temp));
        if (static_cast<float>(temp)<this->UNDETERMINED_VALUE) cout << "WARNING: TOO HIGH -(MOMENTUM) IN " << x << " " <<  y << " " << z << "\n";
      }
      else
        this->Image3DTemp.Put(x, y, z, 0, this->UNDETERMINED_VALUE-1);
    }
    
    //smooth the scalar field
    gaussianBlurring.SetInput (&this->Image3DTemp);
    gaussianBlurring.SetOutput(&this->Image3DTemp);
    gaussianBlurring.Run();
    
    //save the smoothed scalar field
    for (z = 0; z < this->NZ; z++) for (y = 0; y < this->NY; y++) for (x = 0; x < this->NX; x++){
      this->GradE[timeSubdiv][ptVF(x,y,z,i)] = 2*this->VelocityField[i][ptVF(x,y,z,i)] - 2*this->Image3DTemp.Get(x, y, z, 0)/(this->sigma*this->sigma);
    }
  }
}

///Update VelocityField with with the energy gradients
template <class VoxelType> void LargeDefGradLagrange<VoxelType>::UpdateVelocityField(void){
  int x, y, z, i;
  float MaxGrad,MultFactor;
  double LocGrad;
  
  //1) compute the multiplication factor
  MaxGrad=0;
  for (i=0;i<this->NbTimeSubdiv;i++) for (z = 1; z < this->NZ-2; z++) for (y = 1; y < this->NY-2; y++) for (x = 1; x < this->NX-2; x++){
    if ((this->J0[ptSF(x,y,z)]>this->UNDETERMINED_VALUE)&&(this->J1[ptSF(x,y,z)]>this->UNDETERMINED_VALUE)&&(this->DetJacobians[ptSF(x,y,z)]>this->UNDETERMINED_VALUE)&&
         (this->GradJ0[ptVF(x,y,z,0)]>this->UNDETERMINED_VALUE)&&(this->GradJ0[ptVF(x,y,z,1)]>this->UNDETERMINED_VALUE)&&(this->GradJ0[ptVF(x,y,z,2)]>this->UNDETERMINED_VALUE)&&
         (this->GradE[i][ptVF(x,y,z,0)]>this->UNDETERMINED_VALUE)&&(this->GradE[i][ptVF(x,y,z,1)]>this->UNDETERMINED_VALUE)&&(this->GradE[i][ptVF(x,y,z,2)]>this->UNDETERMINED_VALUE)){
           LocGrad=sqrt(pow((double)this->GradE[i][ptVF(x,y,z,0)],2.0)+pow((double)this->GradE[i][ptVF(x,y,z,1)],2.0)+pow((double)this->GradE[i][ptVF(x,y,z,2)],2.0));
           if (MaxGrad<LocGrad) MaxGrad=(float)LocGrad;
         }
  }
  
  if (MaxGrad>this->MaxVelocityUpdate) MultFactor=this->MaxVelocityUpdate/MaxGrad;
  else MultFactor=0.1;
  if (MultFactor>0.1) MultFactor=0.1;
    
  //2) update the vector field
  for (i=0;i<this->NbTimeSubdiv;i++) for (z = 1; z < this->NZ-2; z++) for (y = 1; y < this->NY-2; y++) for (x = 1; x < this->NX-2; x++){
    if ((this->J0[ptSF(x,y,z)]>this->UNDETERMINED_VALUE)&&(this->J1[ptSF(x,y,z)]>this->UNDETERMINED_VALUE)&&(this->DetJacobians[ptSF(x,y,z)]>this->UNDETERMINED_VALUE)&&
         (this->GradJ0[ptVF(x,y,z,0)]>this->UNDETERMINED_VALUE)&&(this->GradJ0[ptVF(x,y,z,1)]>this->UNDETERMINED_VALUE)&&(this->GradJ0[ptVF(x,y,z,2)]>this->UNDETERMINED_VALUE)&&
         (this->GradE[i][ptVF(x,y,z,0)]>this->UNDETERMINED_VALUE)&&(this->GradE[i][ptVF(x,y,z,1)]>this->UNDETERMINED_VALUE)&&(this->GradE[i][ptVF(x,y,z,2)]>this->UNDETERMINED_VALUE)){
      this->VelocityField[i][ptVF(x,y,z,0)]=this->VelocityField[i][ptVF(x,y,z,0)]-this->GradE[i][ptVF(x,y,z,0)]*MultFactor;
      this->VelocityField[i][ptVF(x,y,z,1)]=this->VelocityField[i][ptVF(x,y,z,1)]-this->GradE[i][ptVF(x,y,z,1)]*MultFactor;
      this->VelocityField[i][ptVF(x,y,z,2)]=this->VelocityField[i][ptVF(x,y,z,2)]-this->GradE[i][ptVF(x,y,z,2)]*MultFactor;
    }
    else{
      this->VelocityField[i][ptVF(x,y,z,0)]=this->UNDETERMINED_VALUE;
      this->VelocityField[i][ptVF(x,y,z,1)]=this->UNDETERMINED_VALUE;
      this->VelocityField[i][ptVF(x,y,z,2)]=this->UNDETERMINED_VALUE;
    }
  }

}


///compute the norm of the velocity fields at all subdivision times and the length of all paths
///WARNING: "J0" MUST BE COMPUTED AT THE TIME SUBDIVISION "this->NbTimeSubdiv-1"
template <class VoxelType> void LargeDefGradLagrange<VoxelType>::ComputeNormsAndLengthAndEnergy(){
  int t,z,y,x,d;
  double SumLoc;
  double ValLoc;
  float deltaT;
  
  //0) parameter to compute the norm
  deltaT=1.0/this->NbTimeSubdiv;

  //1) compute the norm of the vector field
  for (t=0;t<this->NbTimeSubdiv;t++){
    SumLoc=0;
    for (z = 1; z < this->NZ-1; z++) for (y = 1; y < this->NY-1; y++) for (x = 1; x < this->NX-1; x++) for (d = 0; d < 3; d++)
            if ((this->VelocityField[t][ptVF(x,y,z,d)]>this->UNDETERMINED_VALUE)&&
                 (this->VelocityField[t][ptVF(x+1,y,z,d)]>this->UNDETERMINED_VALUE)&&(this->VelocityField[t][ptVF(x-1,y,z,d)]>this->UNDETERMINED_VALUE)&&
                 (this->VelocityField[t][ptVF(x,y+1,z,d)]>this->UNDETERMINED_VALUE)&&(this->VelocityField[t][ptVF(x,y-1,z,d)]>this->UNDETERMINED_VALUE)&&
                 (this->VelocityField[t][ptVF(x,y,z+1,d)]>this->UNDETERMINED_VALUE)&&(this->VelocityField[t][ptVF(x,y,z-1,d)]>this->UNDETERMINED_VALUE)){
              ValLoc=static_cast<double>(this->gamma*this->VelocityField[t][ptVF(x,y,z,d)]);
              ValLoc-=static_cast<double>(this->alpha*(this->VelocityField[t][ptVF(x+1,y,z,d)]-2*this->VelocityField[t][ptVF(x,y,z,d)]+this->VelocityField[t][ptVF(x-1,y,z,d)])/(this->DeltaVox*this->DeltaVox));
              ValLoc-=static_cast<double>(this->alpha*(this->VelocityField[t][ptVF(x,y+1,z,d)]-2*this->VelocityField[t][ptVF(x,y,z,d)]+this->VelocityField[t][ptVF(x,y-1,z,d)])/(this->DeltaVox*this->DeltaVox));
              ValLoc-=static_cast<double>(this->alpha*(this->VelocityField[t][ptVF(x,y,z+1,d)]-2*this->VelocityField[t][ptVF(x,y,z,d)]+this->VelocityField[t][ptVF(x,y,z-1,d)])/(this->DeltaVox*this->DeltaVox));
              SumLoc+=ValLoc*ValLoc;
                 }
                 this->norm[t] = static_cast<float>(sqrt(SumLoc));
    //cout << "norm[" << t << "] = " << this->norm[t] << "\n";
  }

  //2) compute the length
  this->length=0;
  for (t=0;t<this->NbTimeSubdiv;t++) this->length+=this->norm[t]*deltaT;
  //cout << "length = " << length << "\n";
  
  
  //3) compute the energy
  this->EnergyVF=0;
  
  //3.1) Energy from the velocity field
  for (t=0;t<this->NbTimeSubdiv;t++) this->EnergyVF+=this->norm[t]*this->norm[t]*deltaT;
  
  //3.2) Energy from the difference between the deformed template and the target
  //this->ComputeJ0(this->NbTimeSubdiv-1);   //TO ADD IF "J0" IS NOT COMPUTED AT "this->NbTimeSubdiv-1"
  ValLoc=0;
  for (z = 1; z < this->NZ-1; z++) for (y = 1; y < this->NY-1; y++) for (x = 1; x < this->NX-1; x++)
        if ((this->J0[ptSF(x,y,z)]>this->UNDETERMINED_VALUE)&&(this->J1[ptSF(x,y,z)]>this->UNDETERMINED_VALUE))
          ValLoc+=static_cast<double>((this->J0[ptSF(x,y,z)]-this->ImTarget[ptSF(x,y,z)])*(this->J0[ptSF(x,y,z)]-this->ImTarget[ptSF(x,y,z)]));
  
  this->EnergyDI=static_cast<float>(ValLoc/(this->NX*this->NY*this->NZ));
  
  //3.3) Total Energy
  this->EnergyTot=EnergyVF+EnergyDI;
  
}


///Reparametrize the velocity field to have a constant speed.
///REMARK: THE VARIABLES Norms AND Length MUST JUST HAVE BEEN COMPUTED
template <class VoxelType> void LargeDefGradLagrange<VoxelType>::SpeedReparametrization(){
  int t,z,y,x,d;
  float* s;  //integral of norm
  float* h;  //inverse of s
  float deltaT;
  float step;
  float ht,sht,ht_wm,ht_wp;
  int ht_i;
  float h_derivative;
  
  deltaT=1.0/this->NbTimeSubdiv;
  
  //allocations
  s = new float [this->NbTimeSubdiv+1];
  h = new float [this->NbTimeSubdiv+1];
  

  //compute s
  s[0]=0;   //we consider we start at subdivision time 0 and finish at 1
  for (t=1;t<this->NbTimeSubdiv+1;t++){
    s[t]=s[t-1]+this->norm[t-1]*deltaT/this->length;
  }
  
  //for (t=0;t<this->NbTimeSubdiv+1;t++) cout << "s[" << t*deltaT << "] = " << s[t] << "\n";
  
  //estimate h, the inverse of s
  h[0]=0;
  for (t=1;t<this->NbTimeSubdiv;t++){
    //initial guess and research step
    h[t]=t*deltaT;
    step=0.01;
    
    //optimal result search
    while (step>0.0001){ //do the linear stuff
      ht=h[t]/deltaT; ht_i=static_cast<int>(ht);  ht_wm=1-(ht-static_cast<float>(ht_i));  ht_wp=ht-static_cast<float>(ht_i);
      if (ht_i<0){ht_i=0; ht_wm=1; ht_wp=0;}
      if (ht_i>this->NbTimeSubdiv-1) {ht_i=this->NbTimeSubdiv-1; ht_wm=0; ht_wp=1;}
      sht=ht_wm*s[ht_i]+ht_wp*s[ht_i+1];
      
      if (sht>t*deltaT){
        while (sht>t*deltaT){
          h[t]-=step;
          ht=h[t]/deltaT; ht_i=static_cast<int>(ht);  ht_wm=1-(ht-static_cast<float>(ht_i));  ht_wp=ht-static_cast<float>(ht_i);
          if (ht_i<0){ht_i=0; ht_wm=1; ht_wp=0;}
          if (ht_i>this->NbTimeSubdiv-1) {ht_i=this->NbTimeSubdiv-1; ht_wm=0; ht_wp=1;}
          sht=ht_wm*s[ht_i]+ht_wp*s[ht_i+1];
        }
      }
      else{
        while (sht<t*deltaT){
          h[t]+=step;
          ht=h[t]/deltaT; ht_i=static_cast<int>(ht);  ht_wm=1-(ht-static_cast<float>(ht_i));  ht_wp=ht-static_cast<float>(ht_i);
          if (ht_i<0){ht_i=0; ht_wm=1; ht_wp=0;}
          if (ht_i>this->NbTimeSubdiv-1) {ht_i=this->NbTimeSubdiv-1; ht_wm=0; ht_wp=1;}
          sht=ht_wm*s[ht_i]+ht_wp*s[ht_i+1];
        }
        step/=2;
      }
    }
  }
  h[this->NbTimeSubdiv]=1;
  
  //for (t=0;t<this->NbTimeSubdiv+1;t++) cout << "h[" << t*deltaT << "] = " <<  h[t] << "\n";
  
  //compute the derivative of h and reparametrize the vector field
  for (t=0;t<this->NbTimeSubdiv;t++){
    h_derivative=(h[t+1]-h[t])/deltaT;
    //cout << "h'[" << t << "] = " <<  h_derivative << "\n";
    
    for (z = 0; z < this->NZ; z++) for (y = 0; y < this->NY; y++) for (x = 0; x < this->NX; x++) for (d = 0; d < 3; d++)
            if (this->VelocityField[t][ptVF(x,y,z,d)]>this->UNDETERMINED_VALUE)
              this->VelocityField[t][ptVF(x,y,z,d)]=this->VelocityField[t][ptVF(x,y,z,d)]*h_derivative;

  }
  
  //deallocations
  delete s;
  delete h;
}



///save the result of the gradient descent (Beg 2005) for the current 3D image of the 4D time sequence
template <class VoxelType> void LargeDefGradLagrange<VoxelType>::SaveResultGradientDescent(int TimeLoc){
  int x, y, z;
  
  //init -> compute the direct mapping and import the original input template (non pre-treated)
  this->ComputeDirectMapping();
  
  for (z = 0; z < this->NZ; z++)  for (y = 0; y < this->NY; y++) for (x = 0; x < this->NX; x++)
        this->ImTemplate[ptSF(x,y,z)]=static_cast<float>(this->_input->Get(x, y, z, TimeLoc));
  
  //deform the original template
  this->ComputeJ0(this->NbTimeSubdiv-1);
  
  //save the deformed template in the output of the class
  for (z = 0; z < this->NZ; z++) for (y = 0; y < this->NY; y++) for (x = 0; x < this->NX; x++){
    if (this->J0[ptSF(x,y,z)]>this->UNDETERMINED_VALUE)
      this->_output->Put(x, y, z, TimeLoc, static_cast<VoxelType>(this->J0[ptSF(x,y,z)]));
    else
      this->_output->Put(x, y, z, TimeLoc, static_cast<VoxelType>(-1));
  }
  
  //if required, save the velocity field and the template deformations across the time subdivisions
  if (strcmp(this->PrefixOutputVF,"Null")!=0){
    this->SaveVelocityFields(this->PrefixOutputVF);
    this->SaveDeformations(this->PrefixOutputVF);
  }
  
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
          Temp4DField.Put(x, y, z, TimeLoc, this->VelocityField[TimeLoc][ptVF(x,y,z,0)]*100);  //TO REMOVE "*100"
  strcpy(FileName,Prefix);
  strcat(FileName,VelocityField_X);
  Temp4DField.Write(FileName);

  //save the velocity field in direction Y
  for (TimeLoc=0;TimeLoc<this->NbTimeSubdiv;TimeLoc++)  for (z = 0; z < this->NZ; z++) for (y = 0; y < this->NY; y++) for (x = 0; x < this->NX; x++)
          Temp4DField.Put(x, y, z, TimeLoc, this->VelocityField[TimeLoc][ptVF(x,y,z,1)]*100);  //TO REMOVE "*100"
  strcpy(FileName,Prefix);
  strcat(FileName,VelocityField_Y);
  Temp4DField.Write(FileName);

  //save the velocity field in direction Z
  for (TimeLoc=0;TimeLoc<this->NbTimeSubdiv;TimeLoc++)  for (z = 0; z < this->NZ; z++) for (y = 0; y < this->NY; y++) for (x = 0; x < this->NX; x++)
          Temp4DField.Put(x, y, z, TimeLoc, this->VelocityField[TimeLoc][ptVF(x,y,z,2)]*100);  //TO REMOVE "*100"
  strcpy(FileName,Prefix);
  strcat(FileName,VelocityField_Z);
  Temp4DField.Write(FileName);
  
  cout << "REMARK: THE VELOCITIES ARE SAVED AS INTEGERS AND THEN ARE MULTIPLIED BY 100\n";
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
          if (this->J0[ptSF(x,y,z)]>this->UNDETERMINED_VALUE)
            Temp4DField.Put(x, y, z, TimeLoc, this->J0[ptSF(x,y,z)]);
          else
            Temp4DField.Put(x, y, z, TimeLoc, -1);
    }
  }
  strcpy(FileName,Prefix);
  strcat(FileName,Deformations);
  Temp4DField.Write(FileName);
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




///++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
///                               SUB-FUNCTIONS TO PERFORM THE REGISTRATION  (level 2)
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
    VecTemp[0]=this->UNDETERMINED_VALUE;
    VecTemp[1]=this->UNDETERMINED_VALUE;
    VecTemp[2]=this->UNDETERMINED_VALUE;
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

//Perform a trilinear interplolation and returns the direct mapping at non integer coordinates.
//++++ The output is in VecTemp ++++ 
template <class VoxelType> void LargeDefGradLagrange<VoxelType>::GetCoordFromDirectMapping(int subivId,float x, float y, float z,float VecTemp[3]){
  int xi,yi,zi;
  float xwm,ywm,zwm,xwp,ywp,zwp;
  float wmmm,wmmp,wmpm,wmpp,wpmm,wpmp,wppm,wppp;
  
  xi=static_cast<int>(x);  xwm=1-(x-static_cast<float>(xi));  xwp=x-static_cast<float>(xi);
  yi=static_cast<int>(y);  ywm=1-(y-static_cast<float>(yi));  ywp=y-static_cast<float>(yi);
  zi=static_cast<int>(z);  zwm=1-(z-static_cast<float>(zi));  zwp=z-static_cast<float>(zi);
  
  if ((xi<0.)||(xi>=this->NX-1)||(yi<0.)||(yi>=this->NY-1)||(zi<0.)||(zi>=this->NZ-1)){
    VecTemp[0]=this->UNDETERMINED_VALUE;
    VecTemp[1]=this->UNDETERMINED_VALUE;
    VecTemp[2]=this->UNDETERMINED_VALUE;
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
  
    VecTemp[0]= wmmm*this->DirectMapping[subivId][ptVF(xi,yi,zi,0)];
    VecTemp[0]+=wmmp*this->DirectMapping[subivId][ptVF(xi,yi,zi+1,0)];
    VecTemp[0]+=wmpm*this->DirectMapping[subivId][ptVF(xi,yi+1,zi,0)];
    VecTemp[0]+=wmpp*this->DirectMapping[subivId][ptVF(xi,yi+1,zi+1,0)];
    VecTemp[0]+=wpmm*this->DirectMapping[subivId][ptVF(xi+1,yi,zi,0)];
    VecTemp[0]+=wpmp*this->DirectMapping[subivId][ptVF(xi+1,yi,zi+1,0)];
    VecTemp[0]+=wppm*this->DirectMapping[subivId][ptVF(xi+1,yi+1,zi,0)];
    VecTemp[0]+=wppp*this->DirectMapping[subivId][ptVF(xi+1,yi+1,zi+1,0)];
  
    VecTemp[1]= wmmm*this->DirectMapping[subivId][ptVF(xi,yi,zi,1)];
    VecTemp[1]+=wmmp*this->DirectMapping[subivId][ptVF(xi,yi,zi+1,1)];
    VecTemp[1]+=wmpm*this->DirectMapping[subivId][ptVF(xi,yi+1,zi,1)];
    VecTemp[1]+=wmpp*this->DirectMapping[subivId][ptVF(xi,yi+1,zi+1,1)];
    VecTemp[1]+=wpmm*this->DirectMapping[subivId][ptVF(xi+1,yi,zi,1)];
    VecTemp[1]+=wpmp*this->DirectMapping[subivId][ptVF(xi+1,yi,zi+1,1)];
    VecTemp[1]+=wppm*this->DirectMapping[subivId][ptVF(xi+1,yi+1,zi,1)];
    VecTemp[1]+=wppp*this->DirectMapping[subivId][ptVF(xi+1,yi+1,zi+1,1)];
  
    VecTemp[2]= wmmm*this->DirectMapping[subivId][ptVF(xi,yi,zi,2)];
    VecTemp[2]+=wmmp*this->DirectMapping[subivId][ptVF(xi,yi,zi+1,2)];
    VecTemp[2]+=wmpm*this->DirectMapping[subivId][ptVF(xi,yi+1,zi,2)];
    VecTemp[2]+=wmpp*this->DirectMapping[subivId][ptVF(xi,yi+1,zi+1,2)];
    VecTemp[2]+=wpmm*this->DirectMapping[subivId][ptVF(xi+1,yi,zi,2)];
    VecTemp[2]+=wpmp*this->DirectMapping[subivId][ptVF(xi+1,yi,zi+1,2)];
    VecTemp[2]+=wppm*this->DirectMapping[subivId][ptVF(xi+1,yi+1,zi,2)];
    VecTemp[2]+=wppp*this->DirectMapping[subivId][ptVF(xi+1,yi+1,zi+1,2)];
  }
}


//Perform a trilinear interplolation and returns the inverse mapping at non integer coordinates.
//++++ The output is in VecTemp ++++ 
template <class VoxelType> void LargeDefGradLagrange<VoxelType>::GetCoordFromInverseMapping(int subivId,float x, float y, float z,float VecTemp[3]){
  int xi,yi,zi;
  float xwm,ywm,zwm,xwp,ywp,zwp;
  float wmmm,wmmp,wmpm,wmpp,wpmm,wpmp,wppm,wppp;
  
  xi=static_cast<int>(x);  xwm=1-(x-static_cast<float>(xi));  xwp=x-static_cast<float>(xi);
  yi=static_cast<int>(y);  ywm=1-(y-static_cast<float>(yi));  ywp=y-static_cast<float>(yi);
  zi=static_cast<int>(z);  zwm=1-(z-static_cast<float>(zi));  zwp=z-static_cast<float>(zi);
  
  if ((xi<0.)||(xi>=this->NX-1)||(yi<0.)||(yi>=this->NY-1)||(zi<0.)||(zi>=this->NZ-1)){
    VecTemp[0]=this->UNDETERMINED_VALUE;
    VecTemp[1]=this->UNDETERMINED_VALUE;
    VecTemp[2]=this->UNDETERMINED_VALUE;
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
  
    VecTemp[0]= wmmm*this->InverseMapping[subivId][ptVF(xi,yi,zi,0)];
    VecTemp[0]+=wmmp*this->InverseMapping[subivId][ptVF(xi,yi,zi+1,0)];
    VecTemp[0]+=wmpm*this->InverseMapping[subivId][ptVF(xi,yi+1,zi,0)];
    VecTemp[0]+=wmpp*this->InverseMapping[subivId][ptVF(xi,yi+1,zi+1,0)];
    VecTemp[0]+=wpmm*this->InverseMapping[subivId][ptVF(xi+1,yi,zi,0)];
    VecTemp[0]+=wpmp*this->InverseMapping[subivId][ptVF(xi+1,yi,zi+1,0)];
    VecTemp[0]+=wppm*this->InverseMapping[subivId][ptVF(xi+1,yi+1,zi,0)];
    VecTemp[0]+=wppp*this->InverseMapping[subivId][ptVF(xi+1,yi+1,zi+1,0)];
  
    VecTemp[1]= wmmm*this->InverseMapping[subivId][ptVF(xi,yi,zi,1)];
    VecTemp[1]+=wmmp*this->InverseMapping[subivId][ptVF(xi,yi,zi+1,1)];
    VecTemp[1]+=wmpm*this->InverseMapping[subivId][ptVF(xi,yi+1,zi,1)];
    VecTemp[1]+=wmpp*this->InverseMapping[subivId][ptVF(xi,yi+1,zi+1,1)];
    VecTemp[1]+=wpmm*this->InverseMapping[subivId][ptVF(xi+1,yi,zi,1)];
    VecTemp[1]+=wpmp*this->InverseMapping[subivId][ptVF(xi+1,yi,zi+1,1)];
    VecTemp[1]+=wppm*this->InverseMapping[subivId][ptVF(xi+1,yi+1,zi,1)];
    VecTemp[1]+=wppp*this->InverseMapping[subivId][ptVF(xi+1,yi+1,zi+1,1)];
  
    VecTemp[2]= wmmm*this->InverseMapping[subivId][ptVF(xi,yi,zi,2)];
    VecTemp[2]+=wmmp*this->InverseMapping[subivId][ptVF(xi,yi,zi+1,2)];
    VecTemp[2]+=wmpm*this->InverseMapping[subivId][ptVF(xi,yi+1,zi,2)];
    VecTemp[2]+=wmpp*this->InverseMapping[subivId][ptVF(xi,yi+1,zi+1,2)];
    VecTemp[2]+=wpmm*this->InverseMapping[subivId][ptVF(xi+1,yi,zi,2)];
    VecTemp[2]+=wpmp*this->InverseMapping[subivId][ptVF(xi+1,yi,zi+1,2)];
    VecTemp[2]+=wppm*this->InverseMapping[subivId][ptVF(xi+1,yi+1,zi,2)];
    VecTemp[2]+=wppp*this->InverseMapping[subivId][ptVF(xi+1,yi+1,zi+1,2)];
  }
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
    InterpoGreyLevel=this->UNDETERMINED_VALUE;
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
    InterpoGreyLevel=this->UNDETERMINED_VALUE;
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
  
  //2) BIG LOOP ON THE TIME FRAMES
  for (t = 0; t < this->NT; t++) {
    cout << "Image " << t+1 << " / " << this->NT << "\n";
    
    //2.1 ) INITIALISATION
    this->InitiateGradientDescent(t);
    
    //2.2) SMALLER LOOP OF THE GRADIENT DESCENT
    IterationsNb=0;
    ReparamCount=0;
    while (IterationsNb<this->iteration_nb){
      cout << "Iteration Number " << IterationsNb+1 << " / " << this->iteration_nb << "\n";
      
      //2.2.1) compute the temporary image transformed using the direct mapping from time 0 -> J0
      //cout << "Compute direct Mapping\n";
      this->ComputeDirectMapping();
      
      //2.2.2) compute the temporary image transformed using the inverse mapping from time 1 -> J1
      //cout << "Compute Inverse Mapping\n";
      this->ComputeInverseMapping();
      
      //2.2.3) EVEN SMALLER LOOP ON THE TIME SUBDIVISIONS BETWEEN TIME=0 (template) AND TIME=1 (target)
      for (TimeSubdiv=0;TimeSubdiv<this->NbTimeSubdiv;TimeSubdiv++){
        //cout << "Subdivision Number " << TimeSubdiv+1 << " / " << this->NbTimeSubdiv << "\n";
        
        //2.2.3.1) compute the temporary image transformed using the direct mapping from time 0 -> J0
        //cout << "Computing J0 at time step t=" << TimeSubdiv << "\n";
        this->ComputeJ0(TimeSubdiv);
        
        //2.2.3.2) compute the temporary image transformed using the inverse mapping from time 1 -> J1
        //cout << "Computing J1 at time step t=" << TimeSubdiv << "\n";
        this->ComputeJ1(TimeSubdiv);

        //2.2.3.3) compute gradient of J0
        //cout << "Compute Grad J0\n";
        this->ComputeGradientJ0();
        
        //2.2.3.4) compute the determinant of the jacobian of the transformation
        //cout << "Computing Det Jacobian\n";
        this->ComputeJacobianDeterminant(TimeSubdiv);
        
        //2.2.3.5) compute the gradient of energy
        //cout << "ComputeEnergyGradient\n";
        this->ComputeEnergyGradient(TimeSubdiv);
      }
      
      //2.2.4) update the velocity field
      //cout << "Update Vector Field\n";
      this->UpdateVelocityField();
      
      //2.2.5) Compute the norms and length of the velocity field plus the total energy
      //cout << "Compute Norms, Length and Energy\n";
      this->ComputeNormsAndLengthAndEnergy();
      
      //2.2.6) update convergence test parameters
      IterationsNb++;
      
      cout << "\nIteration " << IterationsNb << ":\n";
      cout << "Norm = ";
      for (TimeSubdiv=0;TimeSubdiv<this->NbTimeSubdiv;TimeSubdiv++) cout << this->norm[TimeSubdiv] << " ";
      cout << "\n";
      cout << "Length = " << this->length << "\n";
      cout << "Energy / Velocity Field = " << this->EnergyVF << "\n";
      cout << "Energy / Difference Images = " << this->EnergyDI << "\n";
      cout << "Energy / Total  = " << this->EnergyTot << "\n \n";
      
      //2.2.7) velocity field reparametrization 
      ReparamCount++;
      if (ReparamCount==this->reparametrization){
        cout << "Speed reparametrization\n";
        this->SpeedReparametrization();
        ReparamCount=0;
      }

    }
    
    //2.3) save the filtered temporary 3D image in VoxelType in the  output image at time t
    //cout << "Save the result\n";
    SaveResultGradientDescent(t);
    
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


