/*=========================================================================
  Date      : $Date: 31.03.2010$
=========================================================================*/

#include <irtkLargeDeformationShooting.h>
#include <fstream>

///constructor
EulerianShooting::EulerianShooting()
{
  //default parameters
  strcpy(this->SourceImageName,"Null");
  strcpy(this->TargetImageName,"Null");
  this->NbTimes=10;
  this->NbIter=10;
  this->weight1 = 1.0; this->sigmaX1 = 3.0; this->sigmaY1 = 3.0; this->sigmaZ1 = 3.0;
  this->weight2 = 0.0; this->sigmaX2 = 3.0; this->sigmaY2 = 3.0; this->sigmaZ2 = 3.0;
  this->weight3 = 0.0; this->sigmaX3 = 3.0; this->sigmaY3 = 3.0; this->sigmaZ3 = 3.0;
  this->weight4 = 0.0; this->sigmaX4 = 3.0; this->sigmaY4 = 3.0; this->sigmaZ4 = 3.0;
  this->Margin = 3;
  this->GreyLevAlign = 0;
  this->GLA_Padding_Src=-1.;
  this->GLA_Padding_Trg=-1.;
  this->alpha=0.001;
  this->MaxUpdate = 0.5;
  this->indicatorInitialMomentum = 0;
  this->OutIniMoTxt = 0;
  this->OutVeloField = 0;
  this->OutDistEnSim = 0;
  this->OutDeformation = 0;
  this->indicatorRungeKutta = 0;
  strcpy(PrefixOutputs,"Outputs");
  // Default indicator is MinMod
  this->indicatorLimiter = 1;
  this->DeltaX = 1.0;      // To compute the gradient step on the space

}

///destructor
EulerianShooting::~EulerianShooting(void)
{
}

///Read and treat input images
void EulerianShooting::ReadAndTreatInputImages(void){
  int x, y, z;
  int DistClosestEdge;
  double mean1,mean2,std1,std2;
  float PaddingValue;
  int NbVoxelsOK;

  //1) BASIC READING STUFFS
  
  //1.1) read files
  ImTemplate.Read(this->SourceImageName);
  ImTarget.Read(this->TargetImageName);
  if (this->indicatorInitialMomentum>0){
    if (this->indicatorInitialMomentum==1){
      InputInitialMomentum.Read_and_Interpolate(this->InputInitialMomentumName,ImTemplate.NX,ImTemplate.NY,ImTemplate.NZ);
    }
    else{
      //put the input image in the initial momentum map to initiate its headers and size
      //InputInitialMomentum.Read(this->SourceImageName);
      InputInitialMomentum.CreateVoidField(ImTemplate.NX,ImTemplate.NY,ImTemplate.NZ,1);
      //fill the scalar field with the values in the vectorized image
      float tempFl;
      FILE * MyFile;
      
      MyFile=fopen(this->InputInitialMomentumName,"r");
      
      for(z=0;z<ImTemplate.NZ;z++) for(y=0;y<ImTemplate.NY;y++) for(x=0;x<ImTemplate.NX;x++){
        fscanf(MyFile,"%f\n",&tempFl);
        InputInitialMomentum.P(tempFl,x, y, z);
      }
      
      fclose(MyFile);
    }
  }
  
  //1.2) check whether  3D or 2D images are opened
  if (ImTemplate.NT>1) cout << "Source image depends on time!!!";
  if (ImTarget.NT>1) cout << "Target image depends on time!!!";

  //1.3) check whether source and target images have the same size
  if ((ImTemplate.NX!=ImTarget.NX)) cout << "Source and target images do not have the same size!!!";
  if ((ImTemplate.NY!=ImTarget.NY)) cout << "Source and target images do not have the same size!!!";
  if ((ImTemplate.NZ!=ImTarget.NZ)) cout << "Source and target images do not have the same size!!!";
  
  //1.4) variables containing the size of the image
  this->NX=ImTemplate.NX;
  this->NY=ImTemplate.NY;
  this->NZ=ImTemplate.NZ;
  this->NT=1;
  
  cout << "Image size: " << this->NX <<  " , "  <<  this->NY  <<  " , "  << this->NZ  << "\n";
  
  
  //for  (z = 0; z < this->NZ; z++)  for (y = 0; y < this->NY; y++) for (x = 0; x < this->NX; x++) InputInitialMomentum.P(-InputInitialMomentum.G(x,y,z)*50.,x,y,z);
  
  
  //2) COMPUTE THE MARGINS
  if (this->NZ > 1)
  {
    for (z = 0; z < this->NZ; z++)  for (y = 0; y < this->NY; y++) for (x = 0; x < this->NX; x++)
    {
      DistClosestEdge=z+1;
      if (y+1<DistClosestEdge) DistClosestEdge=y+1;
      if (x+1<DistClosestEdge) DistClosestEdge=x+1;
      if (this->NZ-z<DistClosestEdge) DistClosestEdge=this->NZ-z;
      if (this->NY-y<DistClosestEdge) DistClosestEdge=this->NY-y;
      if (this->NX-x<DistClosestEdge) DistClosestEdge=this->NX-x;
      if (DistClosestEdge<=this->Margin)
      {
        this->ImTemplate.P(0,x,y,z);
        this->ImTarget.P(0,x,y,z);
      }
    }
  }
  else
  {
    for (y = 0; y < this->NY; y++) for (x = 0; x < this->NX; x++)
    {
      DistClosestEdge=y+1;
      if (x+1<DistClosestEdge) DistClosestEdge=x+1;
      if (this->NY-y<DistClosestEdge) DistClosestEdge=this->NY-y;
      if (this->NX-x<DistClosestEdge) DistClosestEdge=this->NX-x;
      if (DistClosestEdge<=this->Margin)
      {
        this->ImTemplate.P(0,x,y,0);
        this->ImTarget.P(0,x,y,0);
      }
    }
  }

  
  //3) LINEAR ALIGNMENT OF THE GREY LEVELS OF ImTarget ON THOSE OF ImTemplate
  PaddingValue=10;
  
  if (GreyLevAlign!=0){
    //compute mean and std dev of the source and target images
    mean1=0.;
    NbVoxelsOK=0;
    for (z = 0; z < this->NZ; z++)  for (y = 0; y < this->NY; y++) for (x = 0; x < this->NX; x++) if (this->ImTemplate.G(x,y,z)>GLA_Padding_Src){
      mean1+=(double)this->ImTemplate.G(x,y,z);
      NbVoxelsOK++;
    }
    mean1/=(double)(NbVoxelsOK);
    
    mean2=0.;
    NbVoxelsOK=0;
    for (z = 0; z < this->NZ; z++)  for (y = 0; y < this->NY; y++) for (x = 0; x < this->NX; x++) if (this->ImTarget.G(x,y,z)>GLA_Padding_Trg){
      mean2+=(double)this->ImTarget.G(x,y,z);
      NbVoxelsOK++;
    }
    mean2/=(double)(NbVoxelsOK);
    
    std1=0.;
    NbVoxelsOK=0;
    for (z = 0; z < this->NZ; z++)  for (y = 0; y < this->NY; y++) for (x = 0; x < this->NX; x++) if (this->ImTemplate.G(x,y,z)>GLA_Padding_Src){
      std1+=pow((double)this->ImTemplate.G(x,y,z)-mean1,2.);
      NbVoxelsOK++;
    }
    std1/=(double)(NbVoxelsOK);
    std1=sqrt(std1);
    
    std2=0.;
    NbVoxelsOK=0;
    for (z = 0; z < this->NZ; z++)  for (y = 0; y < this->NY; y++) for (x = 0; x < this->NX; x++) if (this->ImTarget.G(x,y,z)>GLA_Padding_Trg){
      std2+=pow((double)this->ImTarget.G(x,y,z)-mean2,2.);
      NbVoxelsOK++;
    }
    std2/=(double)(NbVoxelsOK);
    std2=sqrt(std2);
    
    cout << "Template: mean=" << mean1 << ", stddev=" << std1 << ".    Target: mean=" << mean2 << ", stddev=" << std2 << "\n";
    
    
    for (z = 0; z < this->NZ; z++)  for (y = 0; y < this->NY; y++) for (x = 0; x < this->NX; x++)
      this->ImTarget.P((this->ImTarget.G(x,y,z)-(float)mean2)*(((float)std1)/((float)std2))+(float)mean1,x,y,z);
    
    
    for (z = 0; z < this->NZ; z++)  for (y = 0; y < this->NY; y++) for (x = 0; x < this->NX; x++)
      if ((this->ImTarget.G(x,y,z)<(GLA_Padding_Trg-(float)mean2)*(((float)std1)/((float)std2))+(float)mean1)||(this->ImTarget.G(x,y,z)<GLA_Padding_Src))
        this->ImTarget.P(0.,x,y,z);
    
    //this->ImTarget.Write("TrgRegistration.nii",this->SourceImageName);
    //this->ImTemplate.Write("SrcRegistration.nii",this->SourceImageName);
  }
  
}

///allocation of all variables used in the program and tune some options
void EulerianShooting::AllocateVariablesShooting(void)
{
  // Iteration Number in the flow calculation.
  this->IterationNumber = NbTimes-1;
  
  // distance between two times: there exist NbTime times.
  this->DeltaTimeSubdiv=1./(static_cast<float>(this->IterationNumber));
  
  // InvDiffeo (or Diffeo) is the list of InvDiffeos (or Diffeo) indexed by the time
  this->InvDiffeo.CreateVoidField(this->NX,this->NY,this->NZ,this->NbTimes);
  this->Diffeo.CreateVoidField(this->NX,this->NY,this->NZ,this->NbTimes);
  
  // velocity field 
  this->VelocityField.CreateVoidField(this->NX,this->NY,this->NZ,1);
  
  // This is the adjoint of the velocity field
  this->AdjointVectorField.CreateVoidField(this->NX,this->NY,this->NZ,1);
  
  // Adjoint variable associated with the image
  this->AdjointImage.CreateVoidField(this->NX,this->NY,this->NZ);
  
  // Variable to store the optimal momentum
  this->OptimizedMomentum.CreateVoidField(this->NX,this->NY,this->NZ,1);
  
  // Momentum used in the algorithm
  this->Momentum.CreateVoidField(this->NX,this->NY,this->NZ,1);
  
  // The Initial Momentum is at the core of the Shooting method: it is given by the user otherwise the initial momentum is null.
  this->InitialMomentum.CreateVoidField(this->NX,this->NY,this->NZ,1,0.0);
  if (this->indicatorInitialMomentum>0) {DeepCopy(&this->InputInitialMomentum,&this->InitialMomentum,0);}
  
  // Adjoint momentum
  this->AdjointMomentum.CreateVoidField(this->NX,this->NY,this->NZ);
  
  // Current image in the Shooting method
  this->Image.CreateVoidField(this->NX,this->NY,this->NZ,1);
  
  // Gradient of the image
  this->NablaI.CreateVoidField(this->NX,this->NY,this->NZ,1);
  
  // Gradient of the adjoint of the momentum
  this->NablaAdM.CreateVoidField(this->NX,this->NY,this->NZ,1);
  
  // Temporary variables
  this->TempAdImage.CreateVoidField(this->NX,this->NY,this->NZ,1);
  this->TempAdMomentum.CreateVoidField(this->NX,this->NY,this->NZ,1);
  this->TempScalarField.CreateVoidField(this->NX,this->NY,this->NZ);
  this->TempScalarField3.CreateVoidField(this->NX,this->NY,this->NZ);
  this->TempVectorField.CreateVoidField(this->NX,this->NY,this->NZ);
  this->TempDiffeo.CreateVoidField(this->NX,this->NY,this->NZ,1);
  this->TempInvDiffeo.CreateVoidField(this->NX,this->NY,this->NZ,1);
  
  // Prepare the temporary variables for the different schemes
  if (this->indicatorRungeKutta==1)
  {
    this->TempDiffeoLocalTable = new VectorField[5];
    this->TempInvDiffeoLocalTable = new VectorField[5];
    int h;
    for (h=0;h<5;h++) 
    {
      this->TempDiffeoLocalTable[h].CreateVoidField(this->NX,this->NY,this->NZ,1);
      this->TempInvDiffeoLocalTable[h].CreateVoidField(this->NX,this->NY,this->NZ,1);
    }
  }
  else 
  {
    this->TempInvDiffeoLocal.CreateVoidField(this->NX,this->NY,this->NZ,1);
    this->TempDiffeoLocal.CreateVoidField(this->NX,this->NY,this->NZ,1);
  }
  
  // Gradient of the functional 
  this->GradientMomentum.CreateVoidField(this->NX,this->NY,this->NZ);
  
  //???
  int x,y,z;
  for (z=0;z<this->NZ;z++) for (y=0;y<this->NY;y++) for (x=0;x<this->NX;x++)
  {
    TempInvDiffeo.P(static_cast<float>(x),0,x,y,z,0);
    TempDiffeo.P(static_cast<float>(x),0,x,y,z,0);
    InvDiffeo.P(static_cast<float>(x),0,x,y,z,0);
    Diffeo.P(static_cast<float>(x),0,x,y,z,0);
    TempInvDiffeo.P(static_cast<float>(y),1,x,y,z,0);
    TempDiffeo.P(static_cast<float>(y),1,x,y,z,0);
    InvDiffeo.P(static_cast<float>(y),1,x,y,z,0);
    Diffeo.P(static_cast<float>(y),1,x,y,z,0);
    TempInvDiffeo.P(static_cast<float>(z),2,x,y,z,0);
    TempDiffeo.P(static_cast<float>(z),2,x,y,z,0);
    InvDiffeo.P(static_cast<float>(z),2,x,y,z,0);
    Diffeo.P(static_cast<float>(z),2,x,y,z,0);
  }
  
  // Initiate the class to filter the velocity field
  FFTconvolver.InitiateConvolver(this->NX,this->NY,this->NZ,this->weight1,this->sigmaX1,this->sigmaY1,this->sigmaZ1,this->weight2,this->sigmaX2,this->sigmaY2,this->sigmaZ2,this->weight3,this->sigmaX3,this->sigmaY3,this->sigmaZ3,this->weight4,this->sigmaX4,this->sigmaY4,this->sigmaZ4);
  
  // choose the spatial scheme: default: UpWind, 1: MinMod, 2: SuperBee
  switch(this->indicatorLimiter)
  {
    case 0 :
      this->Limiter = this->UpWindLimiter;
      cout << "Shooting using UpWindLimiter" <<"\n";
      break;
    case 2 :
      this->Limiter = this->SuperBeeLimiter;
      cout << "Shooting using SuperBeeLimiter" <<"\n";
      break;
    default :
      this->Limiter = this->MinModLimiter;
      cout << "Shooting using MinModLimiter" <<"\n";
      break;
  }
}


/// Implements the advection scheme on the inverse of the diffeomorphism and the Euler scheme on the diffeomorphism
void EulerianShooting::SchemeStep(void)
{
	TransportMomentum(&this->InitialMomentum, &this->TempInvDiffeo, &this->Momentum, this->DeltaX);
	TransportImage(&this->ImTemplate, &this->TempInvDiffeo, &this->Image);
	Cpt_Grad_ScalarField(&this->Image,&this->NablaI,0,this->DeltaX);
	ComputeVelocityField();
	int j,x,y,z;
	float temp,deltaBB,deltaB,deltaF,deltaFF,eta;
	float maxEta = 0.0;
        this->TempInvDiffeoLocal.PutToAllVoxels(0.0);
        this->TempDiffeoLocal.PutToAllVoxels(0.0);
	for (j=0;j<3;j++)
	{
		for (z = 2; z < this->NZ-2; z++) for (y = 2; y < this->NY-2; y++) for (x = 2; x < this->NX-2; x++)
		{
			temp=0.0;
			/// Computation on the first dimension of the scheme.
			eta = this->DeltaTimeSubdiv / this->DeltaX * this->VelocityField.G(0,x,y,z);
			if (abs(eta)>maxEta){maxEta = abs(eta);}
			deltaB = (this->TempInvDiffeo.G(j,x,y,z) - this->TempInvDiffeo.G(j,x-1,y,z));
			deltaF = (this->TempInvDiffeo.G(j,x+1,y,z) - this->TempInvDiffeo.G(j,x,y,z));
			if (eta>=0.0)
			{
				deltaBB = (this->TempInvDiffeo.G(j,x-1,y,z) - this->TempInvDiffeo.G(j,x-2,y,z));
				temp += -eta * (deltaB + 0.5 * (1.0 - eta) * (this->Limiter(deltaB,deltaF) - this->Limiter(deltaBB,deltaB)));
			}
			else
			{
				deltaFF = (this->TempInvDiffeo.G(j,x+2,y,z) - this->TempInvDiffeo.G(j,x+1,y,z));				
				temp += eta * (-deltaF + 0.5 * (1.0 + eta) * (this->Limiter(deltaFF,deltaF) - this->Limiter(deltaF,deltaB)));
			}
			/// Computation on the second dimension of the scheme.
			eta = this->DeltaTimeSubdiv / this->DeltaX * this->VelocityField.G(1,x,y,z);
			if (abs(eta)>maxEta){maxEta = abs(eta);}
			deltaB = (this->TempInvDiffeo.G(j,x,y,z) - this->TempInvDiffeo.G(j,x,y-1,z));
			deltaF = (this->TempInvDiffeo.G(j,x,y+1,z) - this->TempInvDiffeo.G(j,x,y,z));
			if (eta>=0.0)
			{
				deltaBB = (this->TempInvDiffeo.G(j,x,y-1,z) - this->TempInvDiffeo.G(j,x,y-2,z));
				temp += -eta * (deltaB + 0.5 * (1.0 - eta) * (this->Limiter(deltaB,deltaF) - this->Limiter(deltaBB,deltaB)));
			}
			else
			{
				deltaFF = (this->TempInvDiffeo.G(j,x,y+2,z) - this->TempInvDiffeo.G(j,x,y+1,z));				
				temp += eta * (-deltaF + 0.5 * (1.0 + eta) * (this->Limiter(deltaFF,deltaF) - this->Limiter(deltaF,deltaB)));
			}
			/// Computation on the third dimension of the scheme.
			eta = this->DeltaTimeSubdiv / this->DeltaX * this->VelocityField.G(2,x,y,z);
			if (abs(eta)>maxEta){maxEta = abs(eta);}
			deltaB = (this->TempInvDiffeo.G(j,x,y,z) - this->TempInvDiffeo.G(j,x,y,z-1));
			deltaF = (this->TempInvDiffeo.G(j,x,y,z+1) - this->TempInvDiffeo.G(j,x,y,z));
			if (eta>=0.0)
			{
				deltaBB = (this->TempInvDiffeo.G(j,x,y,z-1) - this->TempInvDiffeo.G(j,x,y,z-2));
				temp += -eta * (deltaB + 0.5 * (1.0 - eta) * (this->Limiter(deltaB,deltaF) - this->Limiter(deltaBB,deltaB)));
			}
			else
			{
				deltaFF = (this->TempInvDiffeo.G(j,x,y,z+2) - this->TempInvDiffeo.G(j,x,y,z+1));				
				temp += eta * (-deltaF + 0.5 * (1.0 + eta) * (this->Limiter(deltaFF,deltaF) - this->Limiter(deltaF,deltaB)));
			}
			this->TempInvDiffeoLocal.P(temp,j,x,y,z);
		}
	}
	if (this->NZ==1)
	{
		z=0;
		for (j=0;j<2;j++)
		{
			for (y = 2; y < this->NY-2; y++) for (x = 2; x < this->NX-2; x++)
			{
				temp=0.0;
				/// Computation on the first dimension of the scheme.
				eta = this->DeltaTimeSubdiv / this->DeltaX * this->VelocityField.G(0,x,y,z);
				if (abs(eta)>maxEta){maxEta = abs(eta);}
				deltaB = (this->TempInvDiffeo.G(j,x,y,z) - this->TempInvDiffeo.G(j,x-1,y,z));
				deltaF = (this->TempInvDiffeo.G(j,x+1,y,z) - this->TempInvDiffeo.G(j,x,y,z));
				if (eta>=0.0)
				{
					deltaBB = (this->TempInvDiffeo.G(j,x-1,y,z) - this->TempInvDiffeo.G(j,x-2,y,z));
					temp += -eta * (deltaB + 0.5 * (1.0 - eta) * (this->Limiter(deltaB,deltaF) - this->Limiter(deltaBB,deltaB)));
				}
				else
				{
					deltaFF = (this->TempInvDiffeo.G(j,x+2,y,z) - this->TempInvDiffeo.G(j,x+1,y,z));				
					temp += eta * (-deltaF + 0.5 * (1.0 + eta) * (this->Limiter(deltaFF,deltaF) - this->Limiter(deltaF,deltaB)));
				}
				/// Computation on the second dimension of the scheme.
				eta = this->DeltaTimeSubdiv / this->DeltaX * this->VelocityField.G(1,x,y,z);
				if (abs(eta)>maxEta){maxEta = abs(eta);}
				deltaB = (this->TempInvDiffeo.G(j,x,y,z) - this->TempInvDiffeo.G(j,x,y-1,z));
				deltaF = (this->TempInvDiffeo.G(j,x,y+1,z) - this->TempInvDiffeo.G(j,x,y,z));
				if (eta>=0.0)
				{
					deltaBB = (this->TempInvDiffeo.G(j,x,y-1,z) - this->TempInvDiffeo.G(j,x,y-2,z));
					temp += -eta * (deltaB + 0.5 * (1.0 - eta) * (this->Limiter(deltaB,deltaF) - this->Limiter(deltaBB,deltaB)));
				}
				else
				{
					deltaFF = (this->TempInvDiffeo.G(j,x,y+2,z) - this->TempInvDiffeo.G(j,x,y+1,z));				
					temp += eta * (-deltaF + 0.5 * (1.0 + eta) * (this->Limiter(deltaFF,deltaF) - this->Limiter(deltaF,deltaB)));
				}
				this->TempInvDiffeoLocal.P(temp,j,x,y,z);
			}
		}
	}
	for (j=0;j<3;j++) for (z = 0; z < this->NZ; z++) for (y = 0; y < this->NY; y++) for (x = 0; x < this->NX; x++)
	{
		this->TempInvDiffeo.Add(this->TempInvDiffeoLocal.G(j,x,y,z),j,x,y,z);
		this->TempDiffeoLocal.P(this->VelocityField.G(j,this->TempDiffeo.G(0,x,y,z),this->TempDiffeo.G(1,x,y,z),this->TempDiffeo.G(2,x,y,z)),j,x,y,z);	
	}
	for (j=0;j<3;j++) for (z = 0; z < this->NZ; z++) for (y = 0; y < this->NY; y++) for (x = 0; x < this->NX; x++)
	{
		this->TempDiffeo.Add(TempDiffeoLocal.G(j,x,y,z) * this->DeltaTimeSubdiv,j,x,y,z);
	}
	if (maxEta>1){ cout << " CFL condition not respected  :   " << maxEta <<" > 1" <<"\n"; }
}

/// Implements the advection scheme on the inverse of the diffeomorphism and the Euler scheme on the diffeomorphism with the corresponding arguments
void EulerianShooting::SchemeStep(VectorField *TempInvDiffeoLoc, VectorField *, VectorField *Output1, VectorField *Output2,int t1, int t2)
{
	TransportMomentum(&this->InitialMomentum, TempInvDiffeoLoc, &this->Momentum,this->DeltaX);
	TransportImage(&this->ImTemplate, TempInvDiffeoLoc, &this->Image);
	Cpt_Grad_ScalarField(&this->Image,&this->NablaI,0,this->DeltaX);
	ComputeVelocityField();
	int j,x,y,z;
	float temp,deltaBB,deltaB,deltaF,deltaFF,eta;
	float maxEta = 0.0;
	for (j=0;j<3;j++)
	{
		for (z = 2; z < this->NZ-2; z++) for (y = 2; y < this->NY-2; y++) for (x = 2; x < this->NX-2; x++)
		{
			temp=0.0;
			/// Computation on the first dimension of the scheme.
			eta = this->DeltaTimeSubdiv / this->DeltaX * this->VelocityField.G(0,x,y,z);
			if (abs(eta)>maxEta){maxEta = abs(eta);}
			deltaB = (TempInvDiffeoLoc->G(j,x,y,z) - TempInvDiffeoLoc->G(j,x-1,y,z));
			deltaF = (TempInvDiffeoLoc->G(j,x+1,y,z) - TempInvDiffeoLoc->G(j,x,y,z));
			if (eta>=0.0)
			{
				deltaBB = (TempInvDiffeoLoc->G(j,x-1,y,z) - TempInvDiffeoLoc->G(j,x-2,y,z));
				temp += -eta * (deltaB + 0.5 * (1.0 - eta) * (this->Limiter(deltaB,deltaF) - this->Limiter(deltaBB,deltaB)));
			}
			else
			{
				deltaFF = (TempInvDiffeoLoc->G(j,x+2,y,z) - TempInvDiffeoLoc->G(j,x+1,y,z));				
				temp += eta * (-deltaF + 0.5 * (1.0 + eta) * (this->Limiter(deltaFF,deltaF) - this->Limiter(deltaF,deltaB)));
			}
			/// Computation on the second dimension of the scheme.
			eta = this->DeltaTimeSubdiv / this->DeltaX * this->VelocityField.G(1,x,y,z);
			if (abs(eta)>maxEta){maxEta = abs(eta);}
			deltaB = (TempInvDiffeoLoc->G(j,x,y,z) - TempInvDiffeoLoc->G(j,x,y-1,z));
			deltaF = (TempInvDiffeoLoc->G(j,x,y+1,z) - TempInvDiffeoLoc->G(j,x,y,z));
			if (eta>=0.0)
			{
				deltaBB = (TempInvDiffeoLoc->G(j,x,y-1,z) - TempInvDiffeoLoc->G(j,x,y-2,z));
				temp += -eta * (deltaB + 0.5 * (1.0 - eta) * (this->Limiter(deltaB,deltaF) - this->Limiter(deltaBB,deltaB)));
			}
			else
			{
				deltaFF = (TempInvDiffeoLoc->G(j,x,y+2,z) - TempInvDiffeoLoc->G(j,x,y+1,z));				
				temp += eta * (-deltaF + 0.5 * (1.0 + eta) * (this->Limiter(deltaFF,deltaF) - this->Limiter(deltaF,deltaB)));
			}
			/// Computation on the third dimension of the scheme.
			eta = this->DeltaTimeSubdiv / this->DeltaX * this->VelocityField.G(2,x,y,z);
			if (abs(eta)>maxEta){maxEta = abs(eta);}
			deltaB = (TempInvDiffeoLoc->G(j,x,y,z) - TempInvDiffeoLoc->G(j,x,y,z-1));
			deltaF = (TempInvDiffeoLoc->G(j,x,y,z+1) - TempInvDiffeoLoc->G(j,x,y,z));
			if (eta>=0.0)
			{
				deltaBB = (TempInvDiffeoLoc->G(j,x,y,z-1) - TempInvDiffeoLoc->G(j,x,y,z-2));
				temp += -eta * (deltaB + 0.5 * (1.0 - eta) * (this->Limiter(deltaB,deltaF) - this->Limiter(deltaBB,deltaB)));
			}
			else
			{
				deltaFF = (TempInvDiffeoLoc->G(j,x,y,z+2) - TempInvDiffeoLoc->G(j,x,y,z+1));				
				temp += eta * (-deltaF + 0.5 * (1.0 + eta) * (this->Limiter(deltaFF,deltaF) - this->Limiter(deltaF,deltaB)));
			}
			Output1->P(temp,j,x,y,z,t1);
		}
	}
	if (this->NZ==1)
	{
		z=0;
		for (j=0;j<2;j++)
		{
			for (y = 2; y < this->NY-2; y++) for (x = 2; x < this->NX-2; x++)
			{
				temp=0.0;
				/// Computation on the first dimension of the scheme.
				eta = this->DeltaTimeSubdiv / this->DeltaX * this->VelocityField.G(0,x,y,z);
				if (abs(eta)>maxEta){maxEta = abs(eta);}
				deltaB = (TempInvDiffeoLoc->G(j,x,y,z) - TempInvDiffeoLoc->G(j,x-1,y,z));
				deltaF = (TempInvDiffeoLoc->G(j,x+1,y,z) - TempInvDiffeoLoc->G(j,x,y,z));
				if (eta>=0.0)
				{
					deltaBB = (TempInvDiffeoLoc->G(j,x-1,y,z) - TempInvDiffeoLoc->G(j,x-2,y,z));
					temp += -eta * (deltaB + 0.5 * (1.0 - eta) * (this->Limiter(deltaB,deltaF) - this->Limiter(deltaBB,deltaB)));
				}
				else
				{
					deltaFF = (TempInvDiffeoLoc->G(j,x+2,y,z) - TempInvDiffeoLoc->G(j,x+1,y,z));				
					temp += eta * (-deltaF + 0.5 * (1.0 + eta) * (this->Limiter(deltaFF,deltaF) - this->Limiter(deltaF,deltaB)));
				}
				/// Computation on the second dimension of the scheme.
				eta = this->DeltaTimeSubdiv / this->DeltaX * this->VelocityField.G(1,x,y,z);
				if (abs(eta)>maxEta){maxEta = abs(eta);}
				deltaB = (TempInvDiffeoLoc->G(j,x,y,z) - TempInvDiffeoLoc->G(j,x,y-1,z));
				deltaF = (TempInvDiffeoLoc->G(j,x,y+1,z) - TempInvDiffeoLoc->G(j,x,y,z));
				if (eta>=0.0)
				{
					deltaBB = (TempInvDiffeoLoc->G(j,x,y-1,z) - TempInvDiffeoLoc->G(j,x,y-2,z));
					temp += -eta * (deltaB + 0.5 * (1.0 - eta) * (this->Limiter(deltaB,deltaF) - this->Limiter(deltaBB,deltaB)));
				}
				else
				{
					deltaFF = (TempInvDiffeoLoc->G(j,x,y+2,z) - TempInvDiffeoLoc->G(j,x,y+1,z));				
					temp += eta * (-deltaF + 0.5 * (1.0 + eta) * (this->Limiter(deltaFF,deltaF) - this->Limiter(deltaF,deltaB)));
				}
				Output1->P(temp,j,x,y,z,t1);
			}
		}
	}
	for (j=0;j<3;j++) for (z = 0; z < this->NZ; z++) for (y = 0; y < this->NY; y++) for (x = 0; x < this->NX; x++)
	{
		Output2->P(this->DeltaTimeSubdiv*this->VelocityField.G(j,this->TempDiffeo.G(0,x,y,z),this->DeltaTimeSubdiv*this->TempDiffeo.G(1,x,y,z),this->DeltaTimeSubdiv*this->TempDiffeo.G(2,x,y,z)),j,x,y,z,t2);	
	}
	if (maxEta>1){ cout << " CFL condition not respected  :   " << maxEta <<" > 1" <<"\n"; }
}

/// Runge Kutta scheme
void EulerianShooting::RungeKutta(void)
{
	this->SchemeStep(&this->TempInvDiffeo,&this->TempDiffeo,&this->TempInvDiffeoLocalTable[0],&this->TempDiffeoLocalTable[0]);
	SumVectorField(&this->TempInvDiffeo,&this->TempInvDiffeoLocalTable[0],&this->TempInvDiffeoLocalTable[1]);
	SumVectorField(&this->TempDiffeo,&this->TempDiffeoLocalTable[0],&this->TempDiffeoLocalTable[1]);

	this->SchemeStep(&this->TempInvDiffeoLocalTable[1],&this->TempDiffeoLocalTable[1],&this->TempInvDiffeoLocalTable[2],&this->TempDiffeoLocalTable[2]);
	SumVectorField(&this->TempInvDiffeo,&this->TempInvDiffeoLocalTable[1],&this->TempInvDiffeoLocalTable[3],0,0,0,0.75,0.25);
	SumVectorField(&this->TempDiffeo,&this->TempDiffeoLocalTable[1],&this->TempDiffeoLocalTable[3],0,0,0,0.75,0.25);
	AddVectorField(&this->TempInvDiffeoLocalTable[2],&this->TempInvDiffeoLocalTable[3],0.25);
	AddVectorField(&this->TempDiffeoLocalTable[2],&this->TempDiffeoLocalTable[3],0.25);

	this->SchemeStep(&this->TempInvDiffeoLocalTable[3],&this->TempDiffeoLocalTable[3],&this->TempInvDiffeoLocalTable[0],&this->TempDiffeoLocalTable[0]);
	SumVectorField(&this->TempInvDiffeo,&this->TempInvDiffeoLocalTable[3],&this->TempInvDiffeoLocalTable[4],0,0,0,1.0/3.0,2.0/3.0);
	SumVectorField(&this->TempDiffeo,&this->TempDiffeoLocalTable[3],&this->TempDiffeoLocalTable[4],0,0,0,1.0/3.0,2.0/3.0);
	SumVectorField(&this->TempInvDiffeoLocalTable[4],&this->TempInvDiffeoLocalTable[0],&this->TempInvDiffeo,0,0,0,1.0,2.0/3.0);
	SumVectorField(&this->TempDiffeoLocalTable[4],&this->TempDiffeoLocalTable[0],&this->TempDiffeo,0,0,0,1.0,2.0/3.0);
}

/// Computes the velocity field and requires NablaI and Momentum computed
void EulerianShooting::ComputeVelocityField(void)
{
	float temp;
	int i,x,y,z;
	for (i=0;i<3;i++)
	{
		for (z = 0; z < this->NZ; z++) for (y = 0; y < this->NY; y++) for (x = 0; x < this->NX; x++)
		{
			temp = - this->Momentum.G(x,y,z) * this->NablaI.G(i,x,y,z);
			this->FFTconvolver.P(static_cast<float>(temp),x,y,z);
		}
		FFTconvolver.Convolution();
		for (z = 0; z < this->NZ; z++) for (y = 0; y < this->NY; y++) for (x = 0; x < this->NX; x++)
		{
			this->VelocityField.P(FFTconvolver.G(x,y,z),i,x,y,z);
		}
	}
}

/// Compute the velocity field with arguments
void EulerianShooting::ComputeVelocityField(ScalarField *Momentum,VectorField * NablaI)
{	
	float temp;
	int i,x,y,z;
	for (i=0;i<3;i++)
	{
		for (z = 0; z < this->NZ; z++) for (y = 0; y < this->NY; y++) for (x = 0; x < this->NX; x++)
		{
			temp = - Momentum->G(x,y,z) * NablaI->G(i,x,y,z);
			this->FFTconvolver.P(static_cast<float>(temp),x,y,z);
		}
		FFTconvolver.Convolution();
		for (z = 0; z < this->NZ; z++) for (y = 0; y < this->NY; y++) for (x = 0; x < this->NX; x++)
		{
			this->VelocityField.P(FFTconvolver.G(x,y,z),i,x,y,z);
		}
	}
}

/// Update NablaAdM and compute the adjoint vector field 
void EulerianShooting::ComputeAdjointVectorField(void)
{
	Cpt_Grad_ScalarField(&this->AdjointMomentum,&this->NablaAdM,0,this->DeltaX);
	float temp;
	int i,x,y,z;
	for (i=0;i<3;i++)
	{
		for (z = 0; z < this->NZ; z++) for (y = 0; y < this->NY; y++) for (x = 0; x < this->NX; x++)
		{
			temp =  this->Momentum.G(x,y,z) * this->NablaAdM.G(i,x,y,z) - this->AdjointImage.G(x,y,z) * this->NablaI.G(i,x,y,z);
			this->FFTconvolver.P(static_cast<float>(temp),x,y,z);
		}
		FFTconvolver.Convolution();
		for (z = 0; z < this->NZ; z++) for (y = 0; y < this->NY; y++) for (x = 0; x < this->NX; x++)
		{
			this->AdjointVectorField.P(FFTconvolver.G(x,y,z),i,x,y,z);
		}
	}
}

/// Implementation of SuperBee Limiter
float EulerianShooting::SuperBeeLimiter(float a, float b)
{
	float temp1,temp2,result;
	temp1 = MinModLimiter(a,2*b);
	temp2 = MinModLimiter(2*a,b);
	if (a<b){result=b;}
	else {result=a;}
	return result;
}

/// Run the shooting from an initial Momentum and compute the gradient of the energy of the velocity field w.r.t the initial momentum
void EulerianShooting::Shooting(void)
{   
	this->InitializeVariables();
	this->Cost=0.0;
	this->Energy =0.0;
	int k;
	this->Scheme();	
	/// Compute the gradient of the norm w.r.t the initial momentum 
	ScalarProduct(&this->NablaI,&this->VelocityField,&this->GradientMomentum,0,-this->alpha);
	this->Energy += 0.5 * DotProduct(&this->GradientMomentum,&this->InitialMomentum);
	cout << this->Energy << " Energy of the vector field "<<"\n";
	DeepCopy(&this->TempInvDiffeo,&this->InvDiffeo,1);
	DeepCopy(&this->TempDiffeo,&this->Diffeo,1);
	for (k=1;k<this->IterationNumber;k++)
	{
		this->Scheme();
		DeepCopy(&this->TempInvDiffeo,&this->InvDiffeo,k+1);
		DeepCopy(&this->TempDiffeo,&this->Diffeo,k+1);
	}
	TransportImage(&this->ImTemplate, &this->TempInvDiffeo, &this->Image);
	///add the similiraty measure to the cost 
	this->Cost += this->SimilarityMeasure();
	cout << this->Cost  << " Similarity Measure "<<"\n";
	this->Cost += this->Energy;
	cout << this->Cost  << " Global Cost "<<"\n";
}

/// Perform the gradient descent
void EulerianShooting::GradientDescent(int Niter, float gradientStep)
{	
	float temp;
	int localCounter = 0;
        float optimizedCost = (this->ImTemplate.GetMaxAbsVal() + this->ImTarget.GetMaxAbsVal());
        optimizedCost*=(optimizedCost *this->ImTemplate.GetNbVoxels());
	float currentCost = optimizedCost;
	int i; 
	for (i=0;i<Niter;i++)
	{
		cout <<" "<<"\n";
		cout <<"		Gradient Iteration Number "<<i+1<<"\n";
		this->Gradient();
		if (this->Cost<currentCost)
		{
			if (this->Cost < optimizedCost)
			{
				cout <<"Global Cost Decreasing "<<i+1<<"\n";
				DeepCopy(&this->InitialMomentum,&this->OptimizedMomentum,0);
				optimizedCost = this->Cost;
				localCounter = 0;
			}
		}
		if(this->Cost > currentCost) 
		{
			cout <<"Global Cost Increasing "<<i+1<<"\n";
			localCounter++;
			if (localCounter==2) 
			{
				gradientStep *= 0.8;
				localCounter = 0;
			}
			cout << "  Local Counter " << localCounter << " and Gradient Step "<< gradientStep <<"\n";
		}
		currentCost = this->Cost;
	    Cpt_Grad_ScalarField(&this->ImTemplate,&this->NablaI,0,this->DeltaX);
		this->ComputeVelocityField(&this->GradientMomentum,&this->NablaI);
                this->MaxVectorField = this->VelocityField.GetMaxAbsVal();
		if(this->MaxVectorField>this->MaxUpdate){temp = this->MaxUpdate/this->MaxVectorField;}
		else {temp=1.0;}
		AddScalarField(&this->GradientMomentum,&this->InitialMomentum,-temp*gradientStep);
	}
}

/// Compute the gradient of the shooting w.r.t. the initial momentum and stores the result in GradientMomentum
void EulerianShooting::Gradient(void)
{
	int k;
	this->Shooting();
	this->InitializeAdjointVariables();
	TransportMomentum(&this->AdjointImage,&this->Diffeo,&this->TempAdImage,this->DeltaX,this->IterationNumber);
	for (k=this->IterationNumber;k>0;k--)
	{
		TransportMomentum(&this->TempAdImage,&this->InvDiffeo,&this->AdjointImage,this->DeltaX,k);
		TransportImage(&this->TempAdMomentum, &this->InvDiffeo, &this->AdjointMomentum,k);
		TransportImage(&this->ImTemplate, &this->InvDiffeo, &this->Image,k);
		Cpt_Grad_ScalarField(&this->Image,&this->NablaI,0,this->DeltaX);
		TransportMomentum(&this->InitialMomentum,&this->InvDiffeo,&this->Momentum,this->DeltaX,k);
		this->ComputeAdjointVectorField();
		
		// Compute the increment for TempAdImage
		Product(&this->Momentum,&this->AdjointVectorField,&this->TempVectorField);
		Cpt_Grad_Scal_VectorField(&this->TempVectorField,&this->TempScalarField,0,this->DeltaX);
		TransportMomentum(&this->TempScalarField,&this->Diffeo,&this->TempScalarField3,this->DeltaX,k);
		AddScalarField(&this->TempScalarField3,&this->TempAdImage,this->DeltaTimeSubdiv);

		// Compute the increment for TempAdMomentum
		ScalarProduct(&this->NablaI,&this->AdjointVectorField,&this->TempScalarField);
		AddTransportImage(&this->TempScalarField,&this->Diffeo,&this->TempAdMomentum,-this->DeltaTimeSubdiv,k);
	}
	AddScalarField(&this->TempAdMomentum,&this->GradientMomentum,-1.0);
}

/// Initialize the adjoint variables to do when computing the gradient
void EulerianShooting::InitializeAdjointVariables(void)
{
	int x,y,z;
	for (z = 0; z < this->NZ; z++) for (y = 0; y < this->NY; y++) for (x = 0; x < this->NX; x++)
	{
		this->AdjointImage.P(this->ImTarget.G(x,y,z) - this->Image.G(x,y,z),x,y,z);
		this->TempAdMomentum.P(0.0,x,y,z);
	}
}

/// Initialize the temporary diffeomorphisms for the shooting
void EulerianShooting::InitializeVariables(void)
{
	int x,y,z;
	for (z = 0; z < this->NZ; z++) for (y = 0; y < this->NY; y++) for (x = 0; x < this->NX; x++)
	{
		TempInvDiffeo.P(static_cast<float>(x),0,x,y,z,0);
		TempDiffeo.P(static_cast<float>(x),0,x,y,z,0);
		TempInvDiffeo.P(static_cast<float>(y),1,x,y,z,0);
		TempDiffeo.P(static_cast<float>(y),1,x,y,z,0);
		TempInvDiffeo.P(static_cast<float>(z),2,x,y,z,0);
		TempDiffeo.P(static_cast<float>(z),2,x,y,z,0);
	}
}

/// Choice of the scheme !!!!!  TO BE MODIFIED
void EulerianShooting::Scheme(void)
{	
	if (this->indicatorRungeKutta==1) {	this->RungeKutta(); }
	else { this->SchemeStep(); }
}

/// Compute the sum of square of the difference (can be replaced with CalcSqrtSumOfSquaredDif )
float EulerianShooting::SimilarityMeasure()
{
	float result = 0.0;
	int x,y,z;
	for (z = 0; z < this->NZ; z++) for (y = 0; y < this->NY; y++) for (x = 0; x < this->NX; x++)
	{
		result += pow((double)(this->Image.G(x,y,z) - this->ImTarget.G(x,y,z)),2);
	}
	return 0.5*result;
}

/// Run the registration
void EulerianShooting::Run(void)
{
	this->ReadAndTreatInputImages();
	this->AllocateVariablesShooting();
    this->GradientDescent(this->NbIter,0.5);
	this->SaveResult();
}

/// Save the results
void EulerianShooting::SaveResult(void)
{
	char Output_Momentum[256];
	char Output_Initial_VectorFieldX[256];
	char Output_Initial_VectorFieldY[256];
	char Output_Initial_VectorFieldZ[256];

	//save the initial momentum in a nifti file
	strcpy(Output_Momentum,this->PrefixOutputs);
	strcat(Output_Momentum,"_InitMomentum.nii");
	this->ShootingShow();  //<- ???
	this->InitialMomentum.Write(Output_Momentum,this->SourceImageName);
	
	//save the initial momentum in a vectorized text file
        if (this->OutIniMoTxt==1){
          strcpy(Output_Momentum,this->PrefixOutputs);
	  strcat(Output_Momentum,"_InitMomentum.txt");
	  this->ShootingShow();  //<- ???
	  this->InitialMomentum.WriteInAscii(Output_Momentum);
        }
	
	//save the correponding velocity vector field
        if (this->OutVeloField==1){
          Cpt_Grad_ScalarField(&this->ImTemplate,&this->NablaI,0,this->DeltaX);
          this->ComputeVelocityField(&this->InitialMomentum,&this->NablaI);
          strcpy(Output_Initial_VectorFieldX,this->PrefixOutputs);
          strcat(Output_Initial_VectorFieldX,"_Initial_VectorFieldX.nii");
          strcpy(Output_Initial_VectorFieldY,this->PrefixOutputs);
          strcat(Output_Initial_VectorFieldY,"_Initial_VectorFieldY.nii");
          strcpy(Output_Initial_VectorFieldZ,this->PrefixOutputs);
          strcat(Output_Initial_VectorFieldZ,"_Initial_VectorFieldZ.nii");
  
          this->VelocityField.Write(Output_Initial_VectorFieldX,Output_Initial_VectorFieldY,Output_Initial_VectorFieldZ);
        }
	// Save the distance squared and the SSD
        if (this->OutDistEnSim==1){
          char FileName[256];
          strcpy(FileName,this->PrefixOutputs);
          strcat(FileName,"_Distance.txt");
          cout <<FileName<<"\n";
          ofstream MyFile(FileName, ios::out);
          MyFile << "Energy: " << this->Energy << " -- Similarity Measure: " << this->Cost - this->Energy << endl;
          MyFile.close();
        }
}

///Show the deformation of the source image
void EulerianShooting::ShootingShow(void)
{
  char FileName[256];
  char FileName1[256];
  char FileName2[256];
  char FileName3[256];
  char FileName4[256];
  int k;
  ScalarField OutputImage;
  ScalarField OutputMomentum;
  VectorField OutputVectorField;
  //initialisation
  this->InitializeVariables();
  OutputImage.CreateVoidField(this->NX,this->NY,this->NZ,this->NbTimes);
  OutputMomentum.CreateVoidField(this->NX,this->NY,this->NZ,this->NbTimes);
  OutputVectorField.CreateVoidField(this->NX,this->NY,this->NZ,this->NbTimes);
  DeepCopy(&this->ImTemplate,&OutputImage,0);
  DeepCopy(&this->InitialMomentum,&OutputMomentum,0);
  
  //compute the deformations
	Cpt_Grad_ScalarField(&this->ImTemplate,&this->NablaI,0,this->DeltaX);
	this->ComputeVelocityField(&this->InitialMomentum,&this->NablaI);
	DeepCopy(&this->VelocityField,&OutputVectorField,0);
  for (k=0;k<this->IterationNumber;k++)
  {	
  this->Scheme();
	  TransportMomentum(&this->InitialMomentum, &this->TempInvDiffeo, &this->Momentum, this->DeltaX);
	  TransportImage(&this->ImTemplate, &this->TempInvDiffeo, &this->Image);
	  Cpt_Grad_ScalarField(&this->Image,&this->NablaI,0,this->DeltaX);
	  ComputeVelocityField();
    DeepCopy(&this->TempInvDiffeo,&this->InvDiffeo,k+1);
    DeepCopy(&this->TempDiffeo,&this->Diffeo,k+1);
    DeepCopy(&this->Momentum,&OutputMomentum,k+1);
    DeepCopy(&this->Image,&OutputImage,k+1);
	  DeepCopy(&this->VelocityField,&OutputVectorField,k+1);
  }
  
  //save the deformations
  strcpy(FileName,this->PrefixOutputs);
  strcpy(FileName1,this->PrefixOutputs);
  strcpy(FileName2,this->PrefixOutputs);
  strcpy(FileName3,this->PrefixOutputs);
  strcpy(FileName4,this->PrefixOutputs);
  strcat(FileName,"_Deformation.nii");
  strcat(FileName1,"_VelocityField_X.nii");
  strcat(FileName2,"_VelocityField_Y.nii");
  strcat(FileName3,"_VelocityField_Z.nii");
  strcat(FileName4,"_MomentumEvolution.nii");
  if (OutDeformation==1){
    OutputImage.Write(FileName,this->SourceImageName);
    OutputMomentum.Write(FileName4,this->SourceImageName);
  }
  if (OutVeloField==1){
    OutputVectorField.Write(FileName1,FileName2,FileName3);
  }
}
