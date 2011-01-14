/*=========================================================================
 
 Date      : $Date: 29.04.2010$
 Changes   : $Authors: Laurent Risser, Francois-Xavier Vialard$
 
 =========================================================================*/

#include <irtkLargeDeformationGradientLagrange.h>

///++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
///                                   CONSTRUCTOR AND DESTRUCTOR
///++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
LargeDefGradLagrange::LargeDefGradLagrange(void){
	int i;
	
	//default parameters
	epsilon=0.2;
	iteration_nb=10;
	NbTimeSubdiv=10;
	MaxVelocityUpdate=0.4;  //rem: Delta Voxels = 1
	weight1=100.; sigmaX1=1.;  sigmaY1=1.;  sigmaZ1=1.;
	weight2=0.;   sigmaX2=-1.; sigmaY2=-1.; sigmaZ2=-1.;
	weight3=0.;   sigmaX3=-1.; sigmaY3=-1.; sigmaZ3=-1.;
	weight4=0.;   sigmaX4=-1.; sigmaY4=-1.; sigmaZ4=-1.;
	weight5=0.;   sigmaX5=-1.; sigmaY5=-1.; sigmaZ5=-1.;
	weight6=0.;   sigmaX6=-1.; sigmaY6=-1.; sigmaZ6=-1.;
	weight7=0.;   sigmaX7=-1.; sigmaY7=-1.; sigmaZ7=-1.;
	NbKernels=1;
	SplitKernels=0;
	symmetric=0;
	Margin=0;
	WghtVelField=0.000001; //previously 1
	RefMaxGrad=-1.;
	GreyLevAlign=0;
	GLA_Padding_Src=-1.;
	GLA_Padding_Trg=-1.;
	FlowLength=0;
	DetJacobian=0;
	FinalDefVec=0;
	FinalDefInvVec=0;
	CptInitMomentum=0;
	ShowSSD=0;
	strcpy(PrefixInputs,"Null");
	strcpy(PrefixOutputs,"Outputs");
	strcpy(MaskFile,"Null");
	strcpy(MappingSrc,"Null");
	strcpy(MappingTrg,"Null");
	for (i=0;i<100;i++) strcpy(SourceFiles[i],"Null");
	for (i=0;i<100;i++) strcpy(TargetFiles[i],"Null");
	for (i=0;i<100;i++) weightChannel[i]=1.;
	NbChannels=0;
	MeasureTypicAmp=0;
}

LargeDefGradLagrange::~LargeDefGradLagrange(void){}


///++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
///                        SUB-FUNCTIONS TO PERFORM THE REGISTRATION
///++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


///initiate the gradient descent (Beg 2005) for the current 3D image of the 4D time sequence
void LargeDefGradLagrange::ReadAndTreatInputImages(void){
	int x, y, z;
	int DistClosestEdge;
	int i;
	double mean1,mean2,std1,std2;
	ScalarField Mask;
	
	//1) CREATE THE SOURCE AND TARGET IMAGES 3D *[Nb Channels] FOR THE CALCULATIONS
	//    -->  ImTemplate[c].G(x,y,z) = gray level of the c'th template channel at (x,y,z)
	//    -->  ImTarget[c].G(x,y,z)  = gray level of the c'th target channel at (x,y,z)
	this->ImTemplate = new ScalarField [this->NbChannels];
	this->ImTarget = new ScalarField [this->NbChannels];
	
	//2) READ INPUTS
	for (i=0;i<this->NbChannels;i++){
		//2.1) read files
		ImTemplate[i].Read(this->SourceFiles[i]);
		ImTarget[i].Read(this->TargetFiles[i]);
		
		//2.2) check whether  3D or 2D images are opened
		if (ImTemplate[i].NT>1) cout << "Source image " << i << " depends on time!!!";
		if (ImTarget[i].NT>1) cout << "Target image " << i << " depends on time!!!";
		
		//2.3) check whether source and target images have the same size
		if ((ImTemplate[i].NX!=ImTarget[i].NX)) cout << "Source and target images " << i << " do not have the same size!!!";
		if ((ImTemplate[i].NY!=ImTarget[i].NY)) cout << "Source and target images " << i << " do not have the same size!!!";
		if ((ImTemplate[i].NZ!=ImTarget[i].NZ)) cout << "Source and target images " << i << " do not have the same size!!!";
		
		//2.4) check whether the channels have the same size
		if (i>0){
			if ((ImTemplate[i].NX!=ImTemplate[i-1].NX)) cout << "Images " << i << " and " << i-1 << " do not have the same size!!!";
			if ((ImTemplate[i].NY!=ImTemplate[i-1].NY)) cout << "Images " << i << " and " << i-1 << " do not have the same size!!!";
			if ((ImTemplate[i].NZ!=ImTemplate[i-1].NZ)) cout << "Images " << i << " and " << i-1 << " do not have the same size!!!";
		}
	}
	
	//2.5) variables containing the size of the image
	this->NX=ImTemplate[0].NX;
	this->NY=ImTemplate[0].NY;
	this->NZ=ImTemplate[0].NZ;
	this->NT=1;
	
	cout << "Image size: " << this->NX <<  " , "  <<  this->NY  <<  " , "  << this->NZ  << "\n";
	
	//3) CREATE THE MASK
	if (strcmp(this->MaskFile,"Null")!=0){
		//read the mask
		Mask.Read(this->MaskFile);
		
		if ((ImTemplate[0].NX!=Mask.NX)) cout << "The image(s) and the mask do not have the same size!!!";
		if ((ImTemplate[0].NY!=Mask.NY)) cout << "The image(s) and the mask do not have the same size!!!";
		if ((ImTemplate[0].NZ!=Mask.NZ)) cout << "The image(s) and the mask do not have the same size!!!";
		
		//mask the image
		for (i=0;i<this->NbChannels;i++){
			for(z=0;z<this->NZ;z++) for(y=0;y<this->NY;y++) for(x=0;x<this->NX;x++) if (Mask.G(x,y,z)<0.001)
				this->ImTemplate[i].P(0.,x,y,z);
			
			for(z=0;z<this->NZ;z++) for(y=0;y<this->NY;y++) for(x=0;x<this->NX;x++) if (Mask.G(x,y,z)<0.001)
				this->ImTarget[i].P(0.,x,y,z);
		}
	}  
	//4) COMPUTE THE MARGINS
	for (i=0;i<this->NbChannels;i++){
		for (z = 0; z < this->NZ; z++)  for (y = 0; y < this->NY; y++) for (x = 0; x < this->NX; x++){
			DistClosestEdge=z+1;
			if (y+1<DistClosestEdge) DistClosestEdge=y+1;
			if (x+1<DistClosestEdge) DistClosestEdge=x+1;
			if (this->NZ-z<DistClosestEdge) DistClosestEdge=this->NZ-z;
			if (this->NY-y<DistClosestEdge) DistClosestEdge=this->NY-y;
			if (this->NX-x<DistClosestEdge) DistClosestEdge=this->NX-x;
			if (DistClosestEdge<=this->Margin){
				//this->ImTemplate[i].P(this->ImTemplate[i].G(x,y,z)*((DistClosestEdge-1.)/this->Margin)*((DistClosestEdge-1.)/this->Margin),x,y,z);
				//this->ImTarget[i].P(this->ImTarget[i].G(x,y,z)*((DistClosestEdge-1.)/this->Margin)*((DistClosestEdge-1.)/this->Margin),x,y,z);
				this->ImTemplate[i].P(0,x,y,z);
				this->ImTarget[i].P(0,x,y,z);
			}
		}
	}
	
	//5) LINEAR ALIGNMENT OF THE GREY LEVELS OF ImTarget ON THOSE OF ImTemplate
	float PaddingValue;
	int NbVoxelsOK;
	PaddingValue=10;
	
	if (GreyLevAlign!=0) for (i=0;i<this->NbChannels;i++){
		//compute mean and std dev of the source and target images
		mean1=0.;
		NbVoxelsOK=0;
		for (z = 0; z < this->NZ; z++)  for (y = 0; y < this->NY; y++) for (x = 0; x < this->NX; x++) if (this->ImTemplate[i].G(x,y,z)>GLA_Padding_Src){
			mean1+=(double)this->ImTemplate[i].G(x,y,z);
			NbVoxelsOK++;
		}
		mean1/=(double)(NbVoxelsOK);
		
		mean2=0.;
		NbVoxelsOK=0;
		for (z = 0; z < this->NZ; z++)  for (y = 0; y < this->NY; y++) for (x = 0; x < this->NX; x++) if (this->ImTarget[i].G(x,y,z)>GLA_Padding_Trg){
			mean2+=(double)this->ImTarget[i].G(x,y,z);
			NbVoxelsOK++;
		}
		mean2/=(double)(NbVoxelsOK);
		
		std1=0.;
		NbVoxelsOK=0;
		for (z = 0; z < this->NZ; z++)  for (y = 0; y < this->NY; y++) for (x = 0; x < this->NX; x++) if (this->ImTemplate[i].G(x,y,z)>GLA_Padding_Src){
			std1+=pow((double)this->ImTemplate[i].G(x,y,z)-mean1,2.);
			NbVoxelsOK++;
		}
		std1/=(double)(NbVoxelsOK);
		std1=sqrt(std1);
		
		std2=0.;
		NbVoxelsOK=0;
		for (z = 0; z < this->NZ; z++)  for (y = 0; y < this->NY; y++) for (x = 0; x < this->NX; x++) if (this->ImTarget[i].G(x,y,z)>GLA_Padding_Trg){
			std2+=pow((double)this->ImTarget[i].G(x,y,z)-mean2,2.);
			NbVoxelsOK++;
		}
		std2/=(double)(NbVoxelsOK);
		std2=sqrt(std2);
		
		cout << "Template: mean=" << mean1 << ", stddev=" << std1 << ".    Target: mean=" << mean2 << ", stddev=" << std2 << "\n";
		
		
		for (z = 0; z < this->NZ; z++)  for (y = 0; y < this->NY; y++) for (x = 0; x < this->NX; x++)
			this->ImTarget[i].P((this->ImTarget[i].G(x,y,z)-(float)mean2)*(((float)std1)/((float)std2))+(float)mean1,x,y,z);
		
		
		for (z = 0; z < this->NZ; z++)  for (y = 0; y < this->NY; y++) for (x = 0; x < this->NX; x++)
			if ((this->ImTarget[i].G(x,y,z)<(GLA_Padding_Trg-(float)mean2)*(((float)std1)/((float)std2))+(float)mean1)||(this->ImTarget[i].G(x,y,z)<GLA_Padding_Src))
				this->ImTarget[i].P(0.,x,y,z);
	}
	
	//this->ImTarget[0].Write("TrgNew.nii",SourceFiles[0]);
	//this->ImTemplate[0].Write("SrcNew.nii",SourceFiles[0]);
	
	
	
	//6) MAPPING OF THE SOURCE IMAGE
	this->MappingSrcImag.CreateVoidField(this->NX,this->NY,this->NZ);
	
	if (strcmp(this->MappingSrc,"Null")!=0){
		Mask.Read(this->MappingSrc);
		for (z = 0; z < this->NZ; z++)  for (y = 0; y < this->NY; y++) for (x = 0; x < this->NX; x++){
			this->MappingSrcImag.P(Mask.G(x,y,z,0),0,x,y,z);
			this->MappingSrcImag.P(Mask.G(x,y,z,1),1,x,y,z);
			this->MappingSrcImag.P(Mask.G(x,y,z,2),2,x,y,z);
		}
	}
	else{
		for (z = 0; z < this->NZ; z++)  for (y = 0; y < this->NY; y++) for (x = 0; x < this->NX; x++){
			this->MappingSrcImag.P(static_cast<float>(x),0,x,y,z);
			this->MappingSrcImag.P(static_cast<float>(y),1,x,y,z);
			this->MappingSrcImag.P(static_cast<float>(z),2,x,y,z);
		}
	}
	
	//7) MAPPING OF THE TARGET IMAGE
	this->MappingTrgImag.CreateVoidField(this->NX,this->NY,this->NZ);
	
	if (strcmp(this->MappingTrg,"Null")!=0){
		Mask.Read(this->MappingTrg);
		for (z = 0; z < this->NZ; z++)  for (y = 0; y < this->NY; y++) for (x = 0; x < this->NX; x++){
			this->MappingTrgImag.P(Mask.G(x,y,z,0),0,x,y,z);
			this->MappingTrgImag.P(Mask.G(x,y,z,1),1,x,y,z);
			this->MappingTrgImag.P(Mask.G(x,y,z,2),2,x,y,z);
		}
	}
	else{
		for (z = 0; z < this->NZ; z++)  for (y = 0; y < this->NY; y++) for (x = 0; x < this->NX; x++){
			this->MappingTrgImag.P(static_cast<float>(x),0,x,y,z);
			this->MappingTrgImag.P(static_cast<float>(y),1,x,y,z);
			this->MappingTrgImag.P(static_cast<float>(z),2,x,y,z);
		}
	}
	
	//8) Identity MAPPING
	this->MappingId.CreateVoidField(this->NX,this->NY,this->NZ);
	for (z = 0; z < this->NZ; z++)  for (y = 0; y < this->NY; y++) for (x = 0; x < this->NX; x++){
		this->MappingId.P(static_cast<float>(x),0,x,y,z);
		this->MappingId.P(static_cast<float>(y),1,x,y,z);
		this->MappingId.P(static_cast<float>(z),2,x,y,z);
	}
	
}


///allocate all variables used for the gradient descent (Beg 2005) of the current 3D image from the treated 4D time sequence.
///Compute also the dimension of the scalar and vector fields in use
void LargeDefGradLagrange::AllocateAllVariables(void){
	int i;
	
	//time step between two subdivision
	this->DeltaTimeSubdiv=1./(static_cast<float>(NbTimeSubdiv-1));
	
	//4) initiate the velocity field
	//... velocity field
	//    -->  VelocityField.G(0,x,y,z,i)= direction ex of the vector at (x,y,z)
	//    -->  VelocityField.G(1,x,y,z,i)= direction ey of the vector at (x,y,z)
	//    -->  VelocityField.G(2,x,y,z,i)= direction ez of the vector at (x,y,z)
	//    -->  where n is the id of the velocity field
	if (strcmp(PrefixInputs,"Null")!=0)
		this->LoadVelocityFields(PrefixInputs);  //NbTimeSubdiv should be checked
	else
		this->VelocityField.CreateVoidField(this->NX,this->NY,this->NZ,this->NbTimeSubdiv);
	
	
	if (SplitKernels!=0){  //contribution of each kernel on the velocity field evaluated
		this->SplittedVelocityField=new VectorField [this->NbKernels];
		
		if (strcmp(PrefixInputs,"Null")!=0){
			this->LoadSplittedVelocityFields(PrefixInputs);  //NbTimeSubdiv should be checked
		}
		else{
			for (i=0;i<this->NbKernels;i++ )
				this->SplittedVelocityField[i].CreateVoidField(this->NX,this->NY,this->NZ,this->NbTimeSubdiv);
		}
	}
	
	//... forward mapping
	//    -->  ForwardMapping.G(0,x,y,z,i) = coordinate x at time i corresponding to (x,y,z) at time 0
	//    -->  ForwardMapping.G(1,x,y,z,i) = coordinate y at time i corresponding to (x,y,z) at time 0
	//    -->  ForwardMapping.G(2,x,y,z,i) = coordinate z at time i corresponding to (x,y,z) at time 0
	this->ForwardMapping.CreateVoidField(this->NX,this->NY,this->NZ,this->NbTimeSubdiv);
	
	//... backward mapping
	//    -->  BackwardMapping.G(0,x,y,z,i) = coordinate x at time i corresponding to (x,y,z) at time 1
	//    -->  BackwardMapping.G(1,x,y,z,i) = coordinate y at time i corresponding to (x,y,z) at time 1
	//    -->  BackwardMapping.G(2,x,y,z,i) = coordinate z at time i corresponding to (x,y,z) at time 1
	this->BackwardMapping.CreateVoidField(this->NX,this->NY,this->NZ,this->NbTimeSubdiv);
	
	//... temporary image transformed using the forward mapping from time 0
	//    -->  J0.G(x,y,z) = gray level of the transformed image J0 at (x,y,z)
	this->J0.CreateVoidField(this->NX,this->NY,this->NZ);
	
	//... temporary image transformed using the backward mapping from time 1
	//    -->  J1.G(x,y,z) = gray level of the transformed image J1 at (x,y,z)
	this->J1.CreateVoidField(this->NX,this->NY,this->NZ);
	
	//... gradient of J
	//    -->  GradJ.G(0,x,y,z)= gradient of J0 in direction ex at (x,y,z)
	//    -->  GradJ.G(1,x,y,z)= gradient of J0 in direction ey at (x,y,z)
	//    -->  GradJ.G(2,x,y,z)= gradient of J0 in direction ez at (x,y,z)
	this->GradJ.CreateVoidField(this->NX,this->NY,this->NZ);
	
	//... determinent of the Jacobians  
	//    -->  Jacobians.G(x,y,z)= determinant of the jacobian at (x,y,z)
	this->DetJacobians.CreateVoidField(this->NX,this->NY,this->NZ);
	
	//... Energy Gradient
	//    -->  GradE.G(0,i,x,y,z) = Energy gradient at time i in direction ex at (x,y,z)
	//    -->  GradE.G(1,i,x,y,z) = Energy gradient at time i in direction ey at (x,y,z)
	//    -->  GradE.G(2,i,x,y,z) = Energy gradient at time i in direction ez at (x,y,z)
	this->GradE.CreateVoidField(this->NX,this->NY,this->NZ,this->NbTimeSubdiv);
	
	//contribution of each kernel in the energy gradient
	if (this->SplitKernels!=0){
		this->SplittedGradE=new VectorField [this->NbKernels];
		for (i=0;i<this->NbKernels;i++ )
			this->SplittedGradE[i].CreateVoidField(this->NX,this->NY,this->NZ,this->NbTimeSubdiv);
	}
	
	//Initiate the initial momentum
	if (CptInitMomentum!=0){
		if ((this->SplitKernels!=0)||(this->NbChannels!=1)){ 
			//we can only compute the initial momentum when we have one channel and no splitted kernel (should be extended)
			cout << "Sorry, we were too lazy to program the estimation of the initial momentum with several channels or splitted kernels!\n";
			CptInitMomentum=0;
		}
		else{
			GradInitialMomentum.CreateVoidField(this->NX,this->NY,this->NZ);
			
			if (strcmp(PrefixInputs,"Null")==0)
				InitialMomentum.CreateVoidField(this->NX,this->NY,this->NZ);
			else{
				char ImName[256];
				char Suffix[256];
				strcpy(ImName,PrefixInputs);
				strcpy(Suffix,"_InitMomentum.nii");
				strcat(ImName,Suffix);
				InitialMomentum.Read_and_Interpolate(ImName,this->NX,this->NY,this->NZ);
			}
		}
	}
}



///Compute the energy gradients
void LargeDefGradLagrange::ComputeEnergyGradient(int timeSubdiv,int IdChannel){
	int x,y,z,i,k;
	float temp;
	
	
	if (SplitKernels==0){  //CASE 1: NON-SPLITTED KERNEL
		//loop on the vector directions (x,y,z)
		for (i=0;i<3;i++){
			//compute the scalar field (one dimension out of the vector field) to smooth
			for (z = 0; z < this->NZ; z++) for (y = 0; y < this->NY; y++) for (x = 0; x < this->NX; x++){
				temp=(this->J0.G(x,y,z) - this->J1.G(x,y,z)) * this->DetJacobians.G(x,y,z) * this->GradJ.G(i,x,y,z);
				this->FFTconvolver.P(static_cast<float>(temp),x,y,z);
			}
			
			//smooth the scalar field
			FFTconvolver.Convolution();
			
			//set the gradient of Energy...
			if (IdChannel==0){//...first channel
				for (z = 0; z < this->NZ; z++) for (y = 0; y < this->NY; y++) for (x = 0; x < this->NX; x++){
					this->GradE.P(this->WghtVelField*2*this->VelocityField.G(i,x,y,z,timeSubdiv) - 2*FFTconvolver.G(x,y,z),i,x,y,z,timeSubdiv);
				}
			}
			else{//...other channels -> just update the values computed before
				for (z = 0; z < this->NZ; z++) for (y = 0; y < this->NY; y++) for (x = 0; x < this->NX; x++){
					this->GradE.P(this->GradE.G(i,x,y,z,timeSubdiv)+weightChannel[IdChannel]*(- 2*FFTconvolver.G(x,y,z)),i,x,y,z,timeSubdiv);
				}
			}
		}
		
		//computation of the initial momentum if asked
		if (CptInitMomentum!=0) if (timeSubdiv==0){
			for (z = 0; z < this->NZ; z++) for (y = 0; y < this->NY; y++) for (x = 0; x < this->NX; x++){
				temp=(this->J0.G(x,y,z) - this->J1.G(x,y,z)) * this->DetJacobians.G(x,y,z);
				this->GradInitialMomentum.P(this->WghtVelField*2*this->InitialMomentum.G(x,y,z)-2*temp,x,y,z);
				//if (this->J0.G(x,y,z)!=this->J1.G(x,y,z)) cout << x << " " << y <<  " " << this->GradInitialMomentum.G(x,y,z) <<"\n";
			}
		}
		
		
	}
	else{  //CASE 2: SPLITTED KERNEL
		for (k=0;k<this->NbKernels;k++){
			//2.1) define the kernel
			if (k==0) this->FFTconvolver.ChangeKernel(this->weight1,this->sigmaX1,this->sigmaY1,this->sigmaZ1);
			if (k==1) this->FFTconvolver.ChangeKernel(this->weight2,this->sigmaX2,this->sigmaY2,this->sigmaZ2);
			if (k==2) this->FFTconvolver.ChangeKernel(this->weight3,this->sigmaX3,this->sigmaY3,this->sigmaZ3);
			if (k==3) this->FFTconvolver.ChangeKernel(this->weight4,this->sigmaX4,this->sigmaY4,this->sigmaZ4);
			if (k==4) this->FFTconvolver.ChangeKernel(this->weight5,this->sigmaX5,this->sigmaY5,this->sigmaZ5);
			if (k==5) this->FFTconvolver.ChangeKernel(this->weight6,this->sigmaX6,this->sigmaY6,this->sigmaZ6);
			if (k==6) this->FFTconvolver.ChangeKernel(this->weight7,this->sigmaX7,this->sigmaY7,this->sigmaZ7);
			
			//2.2) do the work like in case 1 for each kernel
			//loop on the vector directions (x,y,z)
			for (i=0;i<3;i++){
				//compute the scalar field (one dimension out of the vector field) to smooth
				for (z = 0; z < this->NZ; z++) for (y = 0; y < this->NY; y++) for (x = 0; x < this->NX; x++){
					temp=(this->J0.G(x,y,z) - this->J1.G(x,y,z)) * this->DetJacobians.G(x,y,z) * this->GradJ.G(i,x,y,z);
					this->FFTconvolver.P(static_cast<float>(temp),x,y,z);
				}
				
				//smooth the scalar field
				FFTconvolver.Convolution();
				
				//set the gradient of Energy...
				if (IdChannel==0){//...first channel  - rem: weightChannel[0]=1, so not expressed here
					//splitted grad E
					for (z = 0; z < this->NZ; z++) for (y = 0; y < this->NY; y++) for (x = 0; x < this->NX; x++)
						this->SplittedGradE[k].P(this->WghtVelField*2*this->SplittedVelocityField[k].G(i,x,y,z,timeSubdiv) - 2*FFTconvolver.G(x,y,z),i,x,y,z,timeSubdiv);
					
					//contribution to grad E
					if (k==0){
						for (z = 0; z < this->NZ; z++) for (y = 0; y < this->NY; y++) for (x = 0; x < this->NX; x++)
							this->GradE.P(this->WghtVelField*2*this->SplittedVelocityField[k].G(i,x,y,z,timeSubdiv) - 2*FFTconvolver.G(x,y,z),i,x,y,z,timeSubdiv);
					}
					else{
						for (z = 0; z < this->NZ; z++) for (y = 0; y < this->NY; y++) for (x = 0; x < this->NX; x++)
							this->GradE.P(this->GradE.G(i,x,y,z,timeSubdiv)- 2*FFTconvolver.G(x,y,z),i,x,y,z,timeSubdiv);
					}
				}
				else{//...other channels -> just update the values computed before
					//splitted grad E
					for (z = 0; z < this->NZ; z++) for (y = 0; y < this->NY; y++) for (x = 0; x < this->NX; x++)
						this->SplittedGradE[k].P(this->SplittedGradE[k].G(i,x,y,z,timeSubdiv)+weightChannel[IdChannel]*(-2*FFTconvolver.G(x,y,z)),i,x,y,z,timeSubdiv);
					
					//contribution to grad E
					for (z = 0; z < this->NZ; z++) for (y = 0; y < this->NY; y++) for (x = 0; x < this->NX; x++)
						this->GradE.P(this->GradE.G(i,x,y,z,timeSubdiv)+weightChannel[IdChannel]*(-2*FFTconvolver.G(x,y,z)),i,x,y,z,timeSubdiv);
				}
			}
		}
	}
}


///Update VelocityField with with the energy gradients.
///Return MaxGrad/this->RefMaxGrad
float LargeDefGradLagrange::UpdateVelocityField(int IterationNb){
	int x, y, z, i, k;
	float MaxGrad,MultFactor;
	double LocGrad;
	
	//1) Compute the maximum of gradient in all time frames...
	//...3D images
	MaxGrad=0;
	for (i=0;i<this->NbTimeSubdiv;i++) for (z = 1; z < this->NZ-2; z++) for (y = 1; y < this->NY-2; y++) for (x = 1; x < this->NX-2; x++){
		LocGrad=sqrt(pow((double)this->GradE.G(0,x,y,z,i),2.)+pow((double)this->GradE.G(1,x,y,z,i),2.)+pow((double)this->GradE.G(2,x,y,z,i),2.));
		if (MaxGrad<(float)LocGrad) MaxGrad=(float)LocGrad;
	}
	
	//...2D images
	if (this->NZ==1){
		for (i=0;i<this->NbTimeSubdiv;i++) for (y = 1; y < this->NY-2; y++) for (x = 1; x < this->NX-2; x++){
			LocGrad=sqrt(pow((double)this->GradE.G(0,x,y,0,i),2.)+pow((double)this->GradE.G(1,x,y,0,i),2.));
			if (MaxGrad<(float)LocGrad) MaxGrad=(float)LocGrad;
		}
	}
	
	//2) maximum update control at the first iteration
	if ((IterationNb==0)&&(RefMaxGrad<0.)) {
		this->RefMaxGrad=MaxGrad;
		if (this->RefMaxGrad==0){
			cout << "It seems that the registered images are identical\n";
			this->RefMaxGrad=1;
		}
		cout << "\n\nRefMaxGrad is set to " << this->RefMaxGrad << ". Keep this value if you continue these\n";
		cout << "computations (using -PrefixInputs) to manage well the convergence.\n\n";
	}
	
	//3) compute the MultFactor
	if (MaxGrad>this->RefMaxGrad) MultFactor=this->MaxVelocityUpdate/MaxGrad;
	else MultFactor=this->MaxVelocityUpdate/(this->RefMaxGrad);
	
	//4) Messages
	if ((IterationNb==0)&&(MultFactor>0.01)) cout << "\nThe weight on the kernels is perhaps too low!!!\n \n";
	
	cout << " -> MaxGrad/RefMaxGrad=" << MaxGrad/this->RefMaxGrad  << "\n";
	
	//5) update the vector field...
	//...3D images
	for (i=0;i<this->NbTimeSubdiv;i++) for (z = 1; z < this->NZ-2; z++) for (y = 1; y < this->NY-2; y++) for (x = 1; x < this->NX-2; x++){
		this->VelocityField.P(this->VelocityField.G(0,x,y,z,i)-this->GradE.G(0,x,y,z,i)*MultFactor,0,x,y,z,i);
		this->VelocityField.P(this->VelocityField.G(1,x,y,z,i)-this->GradE.G(1,x,y,z,i)*MultFactor,1,x,y,z,i);
		this->VelocityField.P(this->VelocityField.G(2,x,y,z,i)-this->GradE.G(2,x,y,z,i)*MultFactor,2,x,y,z,i);
	}
	
	//...2D images
	if (this->NZ==1){
		for (i=0;i<this->NbTimeSubdiv;i++) for (y = 1; y < this->NY-2; y++) for (x = 1; x < this->NX-2; x++){
			this->VelocityField.P(this->VelocityField.G(0,x,y,0,i)-this->GradE.G(0,x,y,0,i)*MultFactor,0,x,y,0,i);
			this->VelocityField.P(this->VelocityField.G(1,x,y,0,i)-this->GradE.G(1,x,y,0,i)*MultFactor,1,x,y,0,i);
		}
	}
	
	//6) IF WE WANT TO MEASURE THE INITIAL MOMENTUM
	if (CptInitMomentum!=0){
		//...3D image
		for (z = 1; z < this->NZ-2; z++) for (y = 1; y < this->NY-2; y++) for (x = 1; x < this->NX-2; x++){
			this->InitialMomentum.P(this->InitialMomentum.G(x,y,z)+this->GradInitialMomentum.G(x,y,z)*MultFactor,x,y,z);
		}
		
		//...2D images
		if (this->NZ==1){
			for (y = 1; y < this->NY-2; y++) for (x = 1; x < this->NX-2; x++){
				this->InitialMomentum.P(this->InitialMomentum.G(x,y,0)+this->GradInitialMomentum.G(x,y,0)*MultFactor,x,y,0);
			}
		}
	}
	
	
	//7) IF SPLITTED KERNEL: update the contribution of each kernel in the velocity field
	if (SplitKernels!=0){
		for (k=0;k<this->NbKernels;k++){
			//...3D images
			for (i=0;i<this->NbTimeSubdiv;i++) for (z = 1; z < this->NZ-2; z++) for (y = 1; y < this->NY-2; y++) for (x = 1; x < this->NX-2; x++){
				this->SplittedVelocityField[k].P(this->SplittedVelocityField[k].G(0,x,y,z,i)-this->SplittedGradE[k].G(0,x,y,z,i)*MultFactor,0,x,y,z,i);
				this->SplittedVelocityField[k].P(this->SplittedVelocityField[k].G(1,x,y,z,i)-this->SplittedGradE[k].G(1,x,y,z,i)*MultFactor,1,x,y,z,i);
				this->SplittedVelocityField[k].P(this->SplittedVelocityField[k].G(2,x,y,z,i)-this->SplittedGradE[k].G(2,x,y,z,i)*MultFactor,2,x,y,z,i);
			}
			
			//...2D images
			if (this->NZ==1){
				for (i=0;i<this->NbTimeSubdiv;i++) for (y = 1; y < this->NY-2; y++) for (x = 1; x < this->NX-2; x++){
					this->SplittedVelocityField[k].P(this->SplittedVelocityField[k].G(0,x,y,0,i)-this->SplittedGradE[k].G(0,x,y,0,i)*MultFactor,0,x,y,0,i);
					this->SplittedVelocityField[k].P(this->SplittedVelocityField[k].G(1,x,y,0,i)-this->SplittedGradE[k].G(1,x,y,0,i)*MultFactor,1,x,y,0,i);
				}
			}
		}
	}
	
	return MaxGrad/this->RefMaxGrad;
}


///save the result of the gradient descent (Beg 2005) for the current 3D image of the 4D time sequence
void LargeDefGradLagrange::SaveResultGradientDescent(void){
	//init -> compute the forward mapping and import the original input template (non pre-treated)
	CptMappingFromVeloField(0,&this->MappingSrcImag,&this->VelocityField,&this->ForwardMapping);
	
	//whole transformations
	this->SaveVelocityFields(&this->VelocityField,this->PrefixOutputs);
	this->SaveDeformations(this->PrefixOutputs);
	if (this->FlowLength==1) this->SaveGlobalFlowLength(this->PrefixOutputs);
	if (this->DetJacobian==1) this->SaveDetJacobian(this->PrefixOutputs);
	if (this->FinalDefVec==1) this->SaveVecDeformation(this->PrefixOutputs);
	if (this->FinalDefInvVec==1) this->SaveInvVecDeformation(this->PrefixOutputs);
	if (this->CptInitMomentum==1) this->SaveInitMomentum(this->PrefixOutputs);
	
	//transformations due to a kernel when using the sum of kernels
	if (SplitKernels!=0){
		this->SaveSplittedVelocityFields(this->PrefixOutputs);
		this->SaveSplittedDeformations(this->PrefixOutputs);
		if (this->FlowLength==1) this->SaveSplittedFlowLength(this->PrefixOutputs);
		if (this->DetJacobian==1) this->SaveSplittedDetJacobian(this->PrefixOutputs);
		if (this->FinalDefVec==1) this->SaveSplittedVecDeformation(this->PrefixOutputs);
	}
}


///++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
///                          FUNCTIONS TO SAVE AND LOAD THE VARIOUS STRUCTURES
///++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


///load the velocity fields
void LargeDefGradLagrange::LoadVelocityFields(char Prefix[256]){
	char FileNameX[256];
	char FileNameY[256];
	char FileNameZ[256];
	char VelocityField_X[256];
	char VelocityField_Y[256];
	char VelocityField_Z[256];
	
	//1) intialisation
	strcpy(FileNameX,Prefix);
	strcpy(VelocityField_X,"_VelocityField_X.nii");
	strcat(FileNameX,VelocityField_X);
	strcpy(FileNameY,Prefix);
	strcpy(VelocityField_Y,"_VelocityField_Y.nii");
	strcat(FileNameY,VelocityField_Y);
	strcpy(FileNameZ,Prefix);
	strcpy(VelocityField_Z,"_VelocityField_Z.nii");
	strcat(FileNameZ,VelocityField_Z);
	
	this->VelocityField.Read_and_Interpolate(FileNameX,FileNameY,FileNameZ,this->NX,this->NY,this->NZ,1);
}


///load the velocity fields
void LargeDefGradLagrange::LoadSplittedVelocityFields(char Prefix[256]){
	char FileNameX[256];
	char FileNameY[256];
	char FileNameZ[256];
	char VelocityField_X1[256];  char VelocityField_Y1[256];  char VelocityField_Z1[256];
	char VelocityField_X2[256];  char VelocityField_Y2[256];  char VelocityField_Z2[256];
	char VelocityField_X3[256];  char VelocityField_Y3[256];  char VelocityField_Z3[256];
	char VelocityField_X4[256];  char VelocityField_Y4[256];  char VelocityField_Z4[256];
	char VelocityField_X5[256];  char VelocityField_Y5[256];  char VelocityField_Z5[256];
	char VelocityField_X6[256];  char VelocityField_Y6[256];  char VelocityField_Z6[256];
	char VelocityField_X7[256];  char VelocityField_Y7[256];  char VelocityField_Z7[256];
	
	//intialisation
	strcpy(VelocityField_X1,"_VelocityField_X1.nii");  strcpy(VelocityField_Y1,"_VelocityField_Y1.nii");  strcpy(VelocityField_Z1,"_VelocityField_Z1.nii");
	strcpy(VelocityField_X2,"_VelocityField_X2.nii");  strcpy(VelocityField_Y2,"_VelocityField_Y2.nii");  strcpy(VelocityField_Z2,"_VelocityField_Z2.nii");
	strcpy(VelocityField_X3,"_VelocityField_X3.nii");  strcpy(VelocityField_Y3,"_VelocityField_Y3.nii");  strcpy(VelocityField_Z3,"_VelocityField_Z3.nii");
	strcpy(VelocityField_X4,"_VelocityField_X4.nii");  strcpy(VelocityField_Y4,"_VelocityField_Y4.nii");  strcpy(VelocityField_Z4,"_VelocityField_Z4.nii");
	strcpy(VelocityField_X5,"_VelocityField_X5.nii");  strcpy(VelocityField_Y5,"_VelocityField_Y5.nii");  strcpy(VelocityField_Z5,"_VelocityField_Z5.nii");
	strcpy(VelocityField_X6,"_VelocityField_X6.nii");  strcpy(VelocityField_Y6,"_VelocityField_Y6.nii");  strcpy(VelocityField_Z6,"_VelocityField_Z6.nii");
	strcpy(VelocityField_X7,"_VelocityField_X7.nii");  strcpy(VelocityField_Y7,"_VelocityField_Y7.nii");  strcpy(VelocityField_Z7,"_VelocityField_Z7.nii");
	
	//velocity field 1
	if (SplitKernels!=0){
		if (this->NbKernels>0){
			strcpy(FileNameX,Prefix);
			strcat(FileNameX,VelocityField_X1);
			strcpy(FileNameY,Prefix);
			strcat(FileNameY,VelocityField_Y1);
			strcpy(FileNameZ,Prefix);
			strcat(FileNameZ,VelocityField_Z1);
			this->SplittedVelocityField[0].Read_and_Interpolate(FileNameX,FileNameY,FileNameZ,this->NX,this->NY,this->NZ,1);
		}
		
		//velocity field 2
		if (this->NbKernels>1){
			strcpy(FileNameX,Prefix);
			strcat(FileNameX,VelocityField_X2);
			strcpy(FileNameY,Prefix);
			strcat(FileNameY,VelocityField_Y2);
			strcpy(FileNameZ,Prefix);
			strcat(FileNameZ,VelocityField_Z2);
			this->SplittedVelocityField[1].Read_and_Interpolate(FileNameX,FileNameY,FileNameZ,this->NX,this->NY,this->NZ,1);
		}
		
		//velocity field 3
		if (this->NbKernels>2){
			strcpy(FileNameX,Prefix);
			strcat(FileNameX,VelocityField_X3);
			strcpy(FileNameY,Prefix);
			strcat(FileNameY,VelocityField_Y3);
			strcpy(FileNameZ,Prefix);
			strcat(FileNameZ,VelocityField_Z3);
			this->SplittedVelocityField[2].Read_and_Interpolate(FileNameX,FileNameY,FileNameZ,this->NX,this->NY,this->NZ,1);
		}
		
		//velocity field 4
		if (this->NbKernels>3){
			strcpy(FileNameX,Prefix);
			strcat(FileNameX,VelocityField_X4);
			strcpy(FileNameY,Prefix);
			strcat(FileNameY,VelocityField_Y4);
			strcpy(FileNameZ,Prefix);
			strcat(FileNameZ,VelocityField_Z4);
			this->SplittedVelocityField[3].Read_and_Interpolate(FileNameX,FileNameY,FileNameZ,this->NX,this->NY,this->NZ,1);
		}
		
		//velocity field 5
		if (this->NbKernels>4){
			strcpy(FileNameX,Prefix);
			strcat(FileNameX,VelocityField_X5);
			strcpy(FileNameY,Prefix);
			strcat(FileNameY,VelocityField_Y5);
			strcpy(FileNameZ,Prefix);
			strcat(FileNameZ,VelocityField_Z5);
			this->SplittedVelocityField[4].Read_and_Interpolate(FileNameX,FileNameY,FileNameZ,this->NX,this->NY,this->NZ,1);
		}
		
		//velocity field 6
		if (this->NbKernels>5){
			strcpy(FileNameX,Prefix);
			strcat(FileNameX,VelocityField_X6);
			strcpy(FileNameY,Prefix);
			strcat(FileNameY,VelocityField_Y6);
			strcpy(FileNameZ,Prefix);
			strcat(FileNameZ,VelocityField_Z6);
			this->SplittedVelocityField[5].Read_and_Interpolate(FileNameX,FileNameY,FileNameZ,this->NX,this->NY,this->NZ,1);
		}
		
		//velocity field 7
		if (this->NbKernels>6){
			strcpy(FileNameX,Prefix);
			strcat(FileNameX,VelocityField_X7);
			strcpy(FileNameY,Prefix);
			strcat(FileNameY,VelocityField_Y7);
			strcpy(FileNameZ,Prefix);
			strcat(FileNameZ,VelocityField_Z7);
			this->SplittedVelocityField[6].Read_and_Interpolate(FileNameX,FileNameY,FileNameZ,this->NX,this->NY,this->NZ,1);
		}
	}
}



///save the velocity fields
void LargeDefGradLagrange::SaveVelocityFields(VectorField * VelocityFieldLoc,char Prefix[256]){
	char FileNameX[256];
	char FileNameY[256];
	char FileNameZ[256];
	char VelocityField_X[256];
	char VelocityField_Y[256];
	char VelocityField_Z[256];
	
	//intialisation
	strcpy(FileNameX,Prefix);
	strcpy(VelocityField_X,"_VelocityField_X.nii");
	strcat(FileNameX,VelocityField_X);
	strcpy(FileNameY,Prefix);
	strcpy(VelocityField_Y,"_VelocityField_Y.nii");
	strcat(FileNameY,VelocityField_Y);
	strcpy(FileNameZ,Prefix);
	strcpy(VelocityField_Z,"_VelocityField_Z.nii");
	strcat(FileNameZ,VelocityField_Z);
	
	//save the velocity field
	VelocityFieldLoc->Write(FileNameX,FileNameY,FileNameZ);
}

///save the velocity fields
void LargeDefGradLagrange::SaveSplittedVelocityFields(char Prefix[256]){
	char FileNameX[256];
	char FileNameY[256];
	char FileNameZ[256];
	char VelocityField_X1[256];  char VelocityField_Y1[256];  char VelocityField_Z1[256];
	char VelocityField_X2[256];  char VelocityField_Y2[256];  char VelocityField_Z2[256];
	char VelocityField_X3[256];  char VelocityField_Y3[256];  char VelocityField_Z3[256];
	char VelocityField_X4[256];  char VelocityField_Y4[256];  char VelocityField_Z4[256];
	char VelocityField_X5[256];  char VelocityField_Y5[256];  char VelocityField_Z5[256];
	char VelocityField_X6[256];  char VelocityField_Y6[256];  char VelocityField_Z6[256];
	char VelocityField_X7[256];  char VelocityField_Y7[256];  char VelocityField_Z7[256];
	
	//intialisation
	strcpy(VelocityField_X1,"_VelocityField_X1.nii");  strcpy(VelocityField_Y1,"_VelocityField_Y1.nii");  strcpy(VelocityField_Z1,"_VelocityField_Z1.nii");
	strcpy(VelocityField_X2,"_VelocityField_X2.nii");  strcpy(VelocityField_Y2,"_VelocityField_Y2.nii");  strcpy(VelocityField_Z2,"_VelocityField_Z2.nii");
	strcpy(VelocityField_X3,"_VelocityField_X3.nii");  strcpy(VelocityField_Y3,"_VelocityField_Y3.nii");  strcpy(VelocityField_Z3,"_VelocityField_Z3.nii");
	strcpy(VelocityField_X4,"_VelocityField_X4.nii");  strcpy(VelocityField_Y4,"_VelocityField_Y4.nii");  strcpy(VelocityField_Z4,"_VelocityField_Z4.nii");
	strcpy(VelocityField_X5,"_VelocityField_X5.nii");  strcpy(VelocityField_Y5,"_VelocityField_Y5.nii");  strcpy(VelocityField_Z5,"_VelocityField_Z5.nii");
	strcpy(VelocityField_X6,"_VelocityField_X6.nii");  strcpy(VelocityField_Y6,"_VelocityField_Y6.nii");  strcpy(VelocityField_Z6,"_VelocityField_Z6.nii");
	strcpy(VelocityField_X7,"_VelocityField_X7.nii");  strcpy(VelocityField_Y7,"_VelocityField_Y7.nii");  strcpy(VelocityField_Z7,"_VelocityField_Z7.nii");
	
	//velocity field 1
	if (SplitKernels!=0){
		if (this->NbKernels>0){
			strcpy(FileNameX,Prefix);
			strcat(FileNameX,VelocityField_X1);
			strcpy(FileNameY,Prefix);
			strcat(FileNameY,VelocityField_Y1);
			strcpy(FileNameZ,Prefix);
			strcat(FileNameZ,VelocityField_Z1);
			this->SplittedVelocityField[0].Write(FileNameX,FileNameY,FileNameZ);
		}
		
		//velocity field 2
		if (this->NbKernels>1){
			strcpy(FileNameX,Prefix);
			strcat(FileNameX,VelocityField_X2);
			strcpy(FileNameY,Prefix);
			strcat(FileNameY,VelocityField_Y2);
			strcpy(FileNameZ,Prefix);
			strcat(FileNameZ,VelocityField_Z2);
			this->SplittedVelocityField[1].Write(FileNameX,FileNameY,FileNameZ);
		}
		
		//velocity field 3
		if (this->NbKernels>2){
			strcpy(FileNameX,Prefix);
			strcat(FileNameX,VelocityField_X3);
			strcpy(FileNameY,Prefix);
			strcat(FileNameY,VelocityField_Y3);
			strcpy(FileNameZ,Prefix);
			strcat(FileNameZ,VelocityField_Z3);
			this->SplittedVelocityField[2].Write(FileNameX,FileNameY,FileNameZ);
		}
		
		//velocity field 4
		if (this->NbKernels>3){
			strcpy(FileNameX,Prefix);
			strcat(FileNameX,VelocityField_X4);
			strcpy(FileNameY,Prefix);
			strcat(FileNameY,VelocityField_Y4);
			strcpy(FileNameZ,Prefix);
			strcat(FileNameZ,VelocityField_Z4);
			this->SplittedVelocityField[3].Write(FileNameX,FileNameY,FileNameZ);
		}
		
		//velocity field 5
		if (this->NbKernels>4){
			strcpy(FileNameX,Prefix);
			strcat(FileNameX,VelocityField_X5);
			strcpy(FileNameY,Prefix);
			strcat(FileNameY,VelocityField_Y5);
			strcpy(FileNameZ,Prefix);
			strcat(FileNameZ,VelocityField_Z5);
			this->SplittedVelocityField[4].Write(FileNameX,FileNameY,FileNameZ);
		}
		
		//velocity field 6
		if (this->NbKernels>5){
			strcpy(FileNameX,Prefix);
			strcat(FileNameX,VelocityField_X6);
			strcpy(FileNameY,Prefix);
			strcat(FileNameY,VelocityField_Y6);
			strcpy(FileNameZ,Prefix);
			strcat(FileNameZ,VelocityField_Z6);
			this->SplittedVelocityField[5].Write(FileNameX,FileNameY,FileNameZ);
		}
		
		//velocity field 7
		if (this->NbKernels>6){
			strcpy(FileNameX,Prefix);
			strcat(FileNameX,VelocityField_X7);
			strcpy(FileNameY,Prefix);
			strcat(FileNameY,VelocityField_Y7);
			strcpy(FileNameZ,Prefix);
			strcat(FileNameZ,VelocityField_Z7);
			this->SplittedVelocityField[6].Write(FileNameX,FileNameY,FileNameZ);
		}
	}
}



///save the deformations in time subdivisions (not the convergence)
void LargeDefGradLagrange::SaveDeformations(char Prefix[256]){
	int TimeLoc,x, y, z;
	ScalarField Temp4DField;
	ScalarField Temp3DField;
	char FileName[256];
	char Deformations[256];
	char FinalDef[256];
	ScalarField source_image;
	
	//read the original input image of the 1st channel (with no treatments)
	source_image.Read(this->SourceFiles[0]);
	
	//intialisation
	Temp4DField.CreateVoidField(this->NX, this->NY, this->NZ,this->NbTimeSubdiv);
	strcpy(Deformations,"_Deformation.nii");
	
	Temp3DField.CreateVoidField(this->NX, this->NY, this->NZ);
	strcpy(FinalDef,"_FinalDefSrc.nii");
	
	//save the deformations
	for (TimeLoc=0;TimeLoc<this->NbTimeSubdiv;TimeLoc++){
		Project3Dimage(&source_image,&this->ForwardMapping,&this->J0,TimeLoc);
		
		for (z = 0; z < this->NZ; z++) for (y = 0; y < this->NY; y++) for (x = 0; x < this->NX; x++)
			Temp4DField.P(this->J0.G(x,y,z),x, y, z, TimeLoc);
		
		if (TimeLoc==this->NbTimeSubdiv-1)
			for (z = 0; z < this->NZ; z++) for (y = 0; y < this->NY; y++) for (x = 0; x < this->NX; x++)
				Temp3DField.P(this->J0.G(x,y,z),x, y, z);
	}
	
	
	
	/*TO REMOVE   TO REMOVE   TO REMOVE   TO REMOVE   TO REMOVE   TO REMOVE   TO REMOVE   TO REMOVE   TO REMOVE*/
	
	//NEW 2D IMAGES:
	/*for (y = 0; y < this->NY; y++) for (x = 0; x < this->NX; x++) Temp3DField.P(0.,x, y,0);
	 
	 double dx,dy;
	 for (y = 1; y < this->NY-1; y++) for (x = 1; x < this->NX-1; x++) for (dy=-0.45;dy<0.46;dy+=0.05) for (dx=-0.45;dx<0.46;dx+=0.05){
	 
	 //if ((x+dx>=100.5-10)&&(x+dx<=100.5+10)&&(y+dy>=80.5-10)&&(y+dy<=80.5+10)) Temp3DField.Add(1.,x, y, 0);
	 
	 if ((x+dx>=95.5-10)&&(x+dx<=95.5+10)&&(y+dy>=85.5-10)&&(y+dy<=85.5+10)) 
	 if (!((x+dx>=95.5-1)&&(x+dx<=95.5+1)&&(y+dy>=95.5-4)&&(y+dy<=95.5-0))) 
	 Temp3DField.Add(1.,x, y, 0);
	 
	 if ((x+dx>=100.5-2)&&(x+dx<=100.5+2)&&(y+dy>=111.5-2)&&(y+dy<=111.5+2)) Temp3DField.Add(1.,x, y, 0);
	 
	 
	 }*/
	
	
	//all tests:
	/*  for (z = 0; z < this->NZ; z++) for (y = 0; y < this->NY; y++) for (x = 0; x < this->NX; x++)
	 Temp3DField.P(0.,x, y, z);*/
	
	//3D grid
	/*  for (z = 0; z < this->NZ; z++) for (y = 0; y < this->NY; y++) for (x = 0; x < this->NX; x++)
	 if (((double)(z/10)==((double)z)/10.)||((double)(y/10)==((double)y)/10.)||((double)(x/10)==((double)x)/10.))
	 Temp3DField.P(100.,x, y, z);*/
	
	//double circle
	/*  for (z = 0; z < this->NZ; z++) for (y = 0; y < this->NY; y++) for (x = 0; x < this->NX; x++)
	 if (sqrt(pow((float)(y)-80,2)+pow((float)(x)-100,2))<20) Temp3DField.P(100.,x, y, z);
	 
	 for (z = 0; z < this->NZ; z++) for (y = 0; y < this->NY; y++) for (x = 0; x < this->NX; x++)
	 if ((x>=99)&&(x<=101)&&(y>=113+4)&&(y<=115+4)) Temp3DField.P(100.,x, y, z);*/
	
	
	/*  for (z = 0; z < this->NZ; z++) for (y = 0; y < this->NY; y++) for (x = 0; x < this->NX; x++)
	 if (sqrt(pow((float)(y)-80,2)+pow((float)(x)-78,2))<20) Temp3DField.P(100.,x, y, z);
	 
	 for (z = 0; z < this->NZ; z++) for (y = 0; y < this->NY; y++) for (x = 0; x < this->NX; x++)
	 if ((x>=99)&&(x<=103)&&(y>=113+4)&&(y<=114+4)) Temp3DField.P(100.,x, y, z);*/
	
	//complex circle
    //for (z = 0; z < this->NZ; z++) for (y = 0; y < this->NY; y++) for (x = 0; x < this->NX; x++)
	//      if (sqrt(pow((float)(y)-64,2)+pow((float)(x)-64,2))<20) Temp3DField.P(100.,x, y, z);
	
	/* for (z = 0; z < this->NZ; z++) for (y = 0; y < this->NY; y++) for (x = 0; x < this->NX; x++)
	 if (sqrt(pow((float)(y)-64,2)+pow((float)(x)-64,2))<25) Temp3DField.P(100.,x, y, z);
	 
	 for (z = 0; z < this->NZ; z++) for (y = 0; y < this->NY; y++) for (x = 0; x < this->NX; x++)
	 if (sqrt(pow((float)(y)-64,2)+pow((float)(x)-39,2))<2.3) Temp3DField.P(0.,x, y, z);
	 
	 for (z = 0; z < this->NZ; z++) for (y = 0; y < this->NY; y++) for (x = 0; x < this->NX; x++)
	 if (sqrt(pow((float)(y)-64,2)+pow((float)(x)-41,2))<2.) Temp3DField.P(0.,x, y, z);
	 
	 for (z = 0; z < this->NZ; z++) for (y = 0; y < this->NY; y++) for (x = 0; x < this->NX; x++)
	 if (sqrt(pow((float)(y)-64,2)+pow((float)(x)-43,2))<2.) Temp3DField.P(0.,x, y, z);
	 
	 for (z = 0; z < this->NZ; z++) for (y = 0; y < this->NY; y++) for (x = 0; x < this->NX; x++)
	 if (sqrt(pow((float)(y)-64,2)+pow((float)(x)-45,2))<1) Temp3DField.P(0.,x, y, z);
	 */ 
	/*TO REMOVE   TO REMOVE   TO REMOVE   TO REMOVE   TO REMOVE   TO REMOVE   TO REMOVE   TO REMOVE   TO REMOVE*/
	
	strcpy(FileName,Prefix);
	strcat(FileName,Deformations);
	Temp4DField.Write(FileName,SourceFiles[0]);
	
	
	strcpy(FileName,Prefix);
	strcat(FileName,FinalDef);
	Temp3DField.Write(FileName,SourceFiles[0]);
	//Temp3DField.Write(FileName,"./RefImages/ComplexCircleSrc.nii.gz");
}

///save the deformations due to each kernel in the context of sum of kernels
void LargeDefGradLagrange::SaveSplittedDeformations(char Prefix[256]){
	int TimeLoc,x, y, z;
	ScalarField Temp4DField;
	char FileName[256];
	char DeformationsK1[256];
	char DeformationsK2[256];
	char DeformationsK3[256];
	char DeformationsK4[256];
	char DeformationsK5[256];
	char DeformationsK6[256];
	char DeformationsK7[256];
	ScalarField source_image;
	int k;
	
	//read the original input image of the 1st channel (with no treatments)
	source_image.Read(this->SourceFiles[0]);
	
	//intialisation
	Temp4DField.CreateVoidField(this->NX, this->NY, this->NZ,this->NbTimeSubdiv);
	strcpy(DeformationsK1,"_DeformationsK1.nii");
	strcpy(DeformationsK2,"_DeformationsK2.nii");
	strcpy(DeformationsK3,"_DeformationsK3.nii");
	strcpy(DeformationsK4,"_DeformationsK4.nii");
	strcpy(DeformationsK5,"_DeformationsK5.nii");
	strcpy(DeformationsK6,"_DeformationsK6.nii");
	strcpy(DeformationsK7,"_DeformationsK7.nii");
	
	//save the deformations
	for (k=0;k<this->NbKernels;k++){
		//compute the current forward mapping
		if (k==0) CptPartialMappingFromVeloFields(0,&this->MappingSrcImag,&this->MappingId,&this->VelocityField,&this->SplittedVelocityField[0],&this->ForwardMapping);
		if (k==1) CptPartialMappingFromVeloFields(0,&this->MappingSrcImag,&this->MappingId,&this->VelocityField,&this->SplittedVelocityField[1],&this->ForwardMapping);
		if (k==2) CptPartialMappingFromVeloFields(0,&this->MappingSrcImag,&this->MappingId,&this->VelocityField,&this->SplittedVelocityField[2],&this->ForwardMapping);
		if (k==3) CptPartialMappingFromVeloFields(0,&this->MappingSrcImag,&this->MappingId,&this->VelocityField,&this->SplittedVelocityField[3],&this->ForwardMapping);
		if (k==4) CptPartialMappingFromVeloFields(0,&this->MappingSrcImag,&this->MappingId,&this->VelocityField,&this->SplittedVelocityField[4],&this->ForwardMapping);
		if (k==5) CptPartialMappingFromVeloFields(0,&this->MappingSrcImag,&this->MappingId,&this->VelocityField,&this->SplittedVelocityField[5],&this->ForwardMapping);
		if (k==6) CptPartialMappingFromVeloFields(0,&this->MappingSrcImag,&this->MappingId,&this->VelocityField,&this->SplittedVelocityField[6],&this->ForwardMapping);
		
		//compute the deformations
		for (TimeLoc=0;TimeLoc<this->NbTimeSubdiv;TimeLoc++){
			Project3Dimage(&source_image,&this->ForwardMapping,&this->J0,TimeLoc);
			
			for (z = 0; z < this->NZ; z++) for (y = 0; y < this->NY; y++) for (x = 0; x < this->NX; x++)
				Temp4DField.P(this->J0.G(x,y,z),x, y, z, TimeLoc);
		}
		
		//save the deformation
		strcpy(FileName,Prefix);
		if (k==0) strcat(FileName,DeformationsK1);
		if (k==1) strcat(FileName,DeformationsK2);
		if (k==2) strcat(FileName,DeformationsK3);
		if (k==3) strcat(FileName,DeformationsK4);
		if (k==4) strcat(FileName,DeformationsK5);
		if (k==5) strcat(FileName,DeformationsK6);
		if (k==6) strcat(FileName,DeformationsK7);
		Temp4DField.Write(FileName,SourceFiles[0]);
	}
	
	//recompute to entire ForwardMapping
	CptMappingFromVeloField(0,&this->MappingSrcImag,&this->VelocityField,&this->ForwardMapping);
}


///save the vector field that transforms [source] into [target]
void LargeDefGradLagrange::SaveVecDeformation(char Prefix[256]){
	int x, y, z;
	ScalarField Temp3DField;
	char FileName[256];
	char VecDef_X[256];
	char VecDef_Y[256];
	char VecDef_Z[256];
	
	//intialisation
	CptMappingFromVeloField(this->VelocityField.NT-1,&this->MappingTrgImag,&this->VelocityField,&this->BackwardMapping);
	
	Temp3DField.CreateVoidField(this->NX, this->NY, this->NZ);
	
	strcpy(VecDef_X,"_VecDef_X.nii");
	strcpy(VecDef_Y,"_VecDef_Y.nii");
	strcpy(VecDef_Z,"_VecDef_Z.nii");
	
	
	//save the forward mapping in direction X
	for (z = 0; z < this->NZ; z++) for (y = 0; y < this->NY; y++) for (x = 0; x < this->NX; x++)
        Temp3DField.P(this->BackwardMapping.G(0,x,y,z)-(float)x,x, y, z);
	strcpy(FileName,Prefix);
	strcat(FileName,VecDef_X);
	Temp3DField.Write(FileName,SourceFiles[0]);
	
	//save the forward mapping in direction Y
	for (z = 0; z < this->NZ; z++) for (y = 0; y < this->NY; y++) for (x = 0; x < this->NX; x++)
        Temp3DField.P(this->BackwardMapping.G(1,x,y,z)-(float)y,x, y, z);
	strcpy(FileName,Prefix);
	strcat(FileName,VecDef_Y);
	Temp3DField.Write(FileName,SourceFiles[0]);
	
	//save the forward mapping in direction Z
	for (z = 0; z < this->NZ; z++) for (y = 0; y < this->NY; y++) for (x = 0; x < this->NX; x++)
        Temp3DField.P(this->BackwardMapping.G(2,x,y,z)-(float)z,x, y, z);
	strcpy(FileName,Prefix);
	strcat(FileName,VecDef_Z);
	Temp3DField.Write(FileName,SourceFiles[0]);
}



///save the vector field that transforms [source] into [target]
void LargeDefGradLagrange::SaveSplittedVecDeformation(char Prefix[256]){
	int x, y, z,k;
	ScalarField Temp3DField;
	VectorField Temp3DVecField;
	char FileName[256];
	char VecDef_X_K1[256];  char VecDef_Y_K1[256];  char VecDef_Z_K1[256];
	char VecDef_X_K2[256];  char VecDef_Y_K2[256];  char VecDef_Z_K2[256];
	char VecDef_X_K3[256];  char VecDef_Y_K3[256];  char VecDef_Z_K3[256];
	char VecDef_X_K4[256];  char VecDef_Y_K4[256];  char VecDef_Z_K4[256];
	char VecDef_X_K5[256];  char VecDef_Y_K5[256];  char VecDef_Z_K5[256];
	char VecDef_X_K6[256];  char VecDef_Y_K6[256];  char VecDef_Z_K6[256];
	char VecDef_X_K7[256];  char VecDef_Y_K7[256];  char VecDef_Z_K7[256];
	
	//intialisation
	Temp3DField.CreateVoidField(this->NX, this->NY, this->NZ);
	Temp3DVecField.CreateVoidField(this->NX, this->NY, this->NZ);
	
	strcpy(VecDef_X_K1,"_VecDef_X_K1.nii");  strcpy(VecDef_Y_K1,"_VecDef_Y_K1.nii");  strcpy(VecDef_Z_K1,"_VecDef_Z_K1.nii");
	strcpy(VecDef_X_K2,"_VecDef_X_K2.nii");  strcpy(VecDef_Y_K2,"_VecDef_Y_K2.nii");  strcpy(VecDef_Z_K2,"_VecDef_Z_K2.nii");
	strcpy(VecDef_X_K3,"_VecDef_X_K3.nii");  strcpy(VecDef_Y_K3,"_VecDef_Y_K3.nii");  strcpy(VecDef_Z_K3,"_VecDef_Z_K3.nii");
	strcpy(VecDef_X_K4,"_VecDef_X_K4.nii");  strcpy(VecDef_Y_K4,"_VecDef_Y_K4.nii");  strcpy(VecDef_Z_K4,"_VecDef_Z_K4.nii");
	strcpy(VecDef_X_K4,"_VecDef_X_K5.nii");  strcpy(VecDef_Y_K4,"_VecDef_Y_K5.nii");  strcpy(VecDef_Z_K4,"_VecDef_Z_K5.nii");
	strcpy(VecDef_X_K4,"_VecDef_X_K6.nii");  strcpy(VecDef_Y_K4,"_VecDef_Y_K6.nii");  strcpy(VecDef_Z_K4,"_VecDef_Z_K6.nii");
	strcpy(VecDef_X_K4,"_VecDef_X_K7.nii");  strcpy(VecDef_Y_K4,"_VecDef_Y_K7.nii");  strcpy(VecDef_Z_K4,"_VecDef_Z_K7.nii");
	
	for (k=0;k<this->NbKernels;k++){
		
		ComputeLagrangianPartialMapping(0,this->NbTimeSubdiv,&this->VelocityField,&this->SplittedVelocityField[k],&Temp3DVecField);
		
		//save the forward mapping in direction X
		for (z = 0; z < this->NZ; z++) for (y = 0; y < this->NY; y++) for (x = 0; x < this->NX; x++)
			Temp3DField.P(Temp3DVecField.G(0,x,y,z)-(float)x,x, y, z);
		strcpy(FileName,Prefix);
		if (k==0) strcat(FileName,VecDef_X_K1);
		if (k==1) strcat(FileName,VecDef_X_K2);
		if (k==2) strcat(FileName,VecDef_X_K3);
		if (k==3) strcat(FileName,VecDef_X_K4);
		if (k==4) strcat(FileName,VecDef_X_K5);
		if (k==5) strcat(FileName,VecDef_X_K6);
		if (k==6) strcat(FileName,VecDef_X_K7);
		Temp3DField.Write(FileName,SourceFiles[0]);
		
		//save the forward mapping in direction Y
		for (z = 0; z < this->NZ; z++) for (y = 0; y < this->NY; y++) for (x = 0; x < this->NX; x++)
			Temp3DField.P(Temp3DVecField.G(1,x,y,z)-(float)y,x, y, z);
		strcpy(FileName,Prefix);
		if (k==0) strcat(FileName,VecDef_Y_K1);
		if (k==1) strcat(FileName,VecDef_Y_K2);
		if (k==2) strcat(FileName,VecDef_Y_K3);
		if (k==3) strcat(FileName,VecDef_Y_K4);
		if (k==4) strcat(FileName,VecDef_Y_K5);
		if (k==5) strcat(FileName,VecDef_Y_K6);
		if (k==6) strcat(FileName,VecDef_Y_K7);
		Temp3DField.Write(FileName,SourceFiles[0]);
		
		//save the forward mapping in direction Z
		for (z = 0; z < this->NZ; z++) for (y = 0; y < this->NY; y++) for (x = 0; x < this->NX; x++)
			Temp3DField.P(Temp3DVecField.G(2,x,y,z)-(float)z,x, y, z);
		strcpy(FileName,Prefix);
		if (k==0) strcat(FileName,VecDef_Z_K1);
		if (k==1) strcat(FileName,VecDef_Z_K2);
		if (k==2) strcat(FileName,VecDef_Z_K3);
		if (k==3) strcat(FileName,VecDef_Z_K4);
		if (k==4) strcat(FileName,VecDef_Z_K5);
		if (k==5) strcat(FileName,VecDef_Z_K6);
		if (k==6) strcat(FileName,VecDef_Z_K7);
		Temp3DField.Write(FileName,SourceFiles[0]);
	}
}



///save the vector field that transforms [target] into [source]
void LargeDefGradLagrange::SaveInvVecDeformation(char Prefix[256]){
	int x, y, z;
	ScalarField Temp3DField;
	char FileName[256];
	char VecInvDef_X[256];
	char VecInvDef_Y[256];
	char VecInvDef_Z[256];
	
	//intialisation
	CptMappingFromVeloField(0,&this->MappingSrcImag,&this->VelocityField,&this->ForwardMapping);
	
	Temp3DField.CreateVoidField(this->NX, this->NY, this->NZ);
	
	strcpy(VecInvDef_X,"_VecInvDef_X.nii");
	strcpy(VecInvDef_Y,"_VecInvDef_Y.nii");
	strcpy(VecInvDef_Z,"_VecInvDef_Z.nii");
	
	
	//save the forward mapping in direction X
	for (z = 0; z < this->NZ; z++) for (y = 0; y < this->NY; y++) for (x = 0; x < this->NX; x++)
        Temp3DField.P(this->ForwardMapping.G(0,x,y,z,this->NbTimeSubdiv-1)-(float)x,x, y, z);
	strcpy(FileName,Prefix);
	strcat(FileName,VecInvDef_X);
	Temp3DField.Write(FileName,SourceFiles[0]);
	
	//save the forward mapping in direction Y
	for (z = 0; z < this->NZ; z++) for (y = 0; y < this->NY; y++) for (x = 0; x < this->NX; x++)
        Temp3DField.P(this->ForwardMapping.G(1,x,y,z,this->NbTimeSubdiv-1)-(float)y,x, y, z);
	strcpy(FileName,Prefix);
	strcat(FileName,VecInvDef_Y);
	Temp3DField.Write(FileName,SourceFiles[0]);
	
	//save the forward mapping in direction Z
	for (z = 0; z < this->NZ; z++) for (y = 0; y < this->NY; y++) for (x = 0; x < this->NX; x++)
        Temp3DField.P(this->ForwardMapping.G(2,x,y,z,this->NbTimeSubdiv-1)-(float)z,x, y, z);
	strcpy(FileName,Prefix);
	strcat(FileName,VecInvDef_Z);
	Temp3DField.Write(FileName,SourceFiles[0]);
}




///save the  initial momentum that transforms [target] into [source]
void LargeDefGradLagrange::SaveInitMomentum(char Prefix[256]){
	char FileName[256];
	char InitM[256];
	
	if (this->CptInitMomentum!=0){
		//1) Save the niftii image
		strcpy(InitM,"_InitMomentum.nii");
		
		strcpy(FileName,Prefix);
		strcat(FileName,InitM);
		InitialMomentum.Write(FileName,SourceFiles[0]);
		
		//2) save the corresponding vectorized image in a text file
		strcpy(InitM,"_InitMomentum.txt");
		
		strcpy(FileName,Prefix);
		strcat(FileName,InitM);
		InitialMomentum.WriteInAscii(FileName);
	}
}



///save the total length of the flow of deformation from each voxel of the image
void LargeDefGradLagrange::SaveGlobalFlowLength(char Prefix[256]){
	char VeloLength[256];
	char EvoVeloLength[256];
	
	strcpy(VeloLength,"_TotalAOD.nii");
	this->SaveFlowLength(&this->VelocityField,&this->VelocityField,this->PrefixOutputs,VeloLength);
	
	strcpy(EvoVeloLength,"EvoAOD.nii");
	this->SaveEvoFlowLength(&this->VelocityField,&this->VelocityField,this->PrefixOutputs,EvoVeloLength);
}

///save the splitted length of the flow of deformation from each voxel of the image
void LargeDefGradLagrange::SaveSplittedFlowLength(char Prefix[256]){
	char VeloLengthK1[256];
	char VeloLengthK2[256];
	char VeloLengthK3[256];
	char VeloLengthK4[256];
	char VeloLengthK5[256];
	char VeloLengthK6[256];
	char VeloLengthK7[256];
	
	strcpy(VeloLengthK1,"_TotalAOD_K1.nii");
	strcpy(VeloLengthK2,"_TotalAOD_K2.nii");
	strcpy(VeloLengthK3,"_TotalAOD_K3.nii");
	strcpy(VeloLengthK4,"_TotalAOD_K4.nii");
	strcpy(VeloLengthK5,"_TotalAOD_K5.nii");
	strcpy(VeloLengthK6,"_TotalAOD_K6.nii");
	strcpy(VeloLengthK7,"_TotalAOD_K7.nii");
	
	if (this->NbKernels>0) this->SaveFlowLength(&this->VelocityField,&this->SplittedVelocityField[0],this->PrefixOutputs,VeloLengthK1);
	if (this->NbKernels>1) this->SaveFlowLength(&this->VelocityField,&this->SplittedVelocityField[1],this->PrefixOutputs,VeloLengthK2);
	if (this->NbKernels>2) this->SaveFlowLength(&this->VelocityField,&this->SplittedVelocityField[2],this->PrefixOutputs,VeloLengthK3);
	if (this->NbKernels>3) this->SaveFlowLength(&this->VelocityField,&this->SplittedVelocityField[3],this->PrefixOutputs,VeloLengthK4);
	if (this->NbKernels>4) this->SaveFlowLength(&this->VelocityField,&this->SplittedVelocityField[4],this->PrefixOutputs,VeloLengthK5);
	if (this->NbKernels>5) this->SaveFlowLength(&this->VelocityField,&this->SplittedVelocityField[5],this->PrefixOutputs,VeloLengthK6);
	if (this->NbKernels>6) this->SaveFlowLength(&this->VelocityField,&this->SplittedVelocityField[6],this->PrefixOutputs,VeloLengthK7);
	
	
	strcpy(VeloLengthK1,"_EvoAOD_K1.nii");
	strcpy(VeloLengthK2,"_EvoAOD_K2.nii");
	strcpy(VeloLengthK3,"_EvoAOD_K3.nii");
	strcpy(VeloLengthK4,"_EvoAOD_K4.nii");
	strcpy(VeloLengthK5,"_EvoAOD_K5.nii");
	strcpy(VeloLengthK6,"_EvoAOD_K6.nii");
	strcpy(VeloLengthK7,"_EvoAOD_K7.nii");
	
	if (this->NbKernels>0) this->SaveEvoFlowLength(&this->VelocityField,&this->SplittedVelocityField[0],this->PrefixOutputs,VeloLengthK1);
	if (this->NbKernels>1) this->SaveEvoFlowLength(&this->VelocityField,&this->SplittedVelocityField[1],this->PrefixOutputs,VeloLengthK2);
	if (this->NbKernels>2) this->SaveEvoFlowLength(&this->VelocityField,&this->SplittedVelocityField[2],this->PrefixOutputs,VeloLengthK3);
	if (this->NbKernels>3) this->SaveEvoFlowLength(&this->VelocityField,&this->SplittedVelocityField[3],this->PrefixOutputs,VeloLengthK4);
	if (this->NbKernels>4) this->SaveEvoFlowLength(&this->VelocityField,&this->SplittedVelocityField[4],this->PrefixOutputs,VeloLengthK5);
	if (this->NbKernels>5) this->SaveEvoFlowLength(&this->VelocityField,&this->SplittedVelocityField[5],this->PrefixOutputs,VeloLengthK6);
	if (this->NbKernels>6) this->SaveEvoFlowLength(&this->VelocityField,&this->SplittedVelocityField[6],this->PrefixOutputs,VeloLengthK7);
	
	
}




///By following the flow defined by the velocity field 'VeloField4Flow' PROJECT AT T=0 the contribution of
///'VeloField4Measure' in the total length of the flow from each point of the field.
/// * 'VeloField4Measure' is assumed to be part of a linear decomposition of 'VeloField4Flow'.
/// * If 'VeloField4Measure'=='VeloField4Flow' then the length of the flow defined by 'VeloField4Flow'
///   is computed.
void LargeDefGradLagrange::SaveFlowLength(VectorField * VeloField4Flow,VectorField * VeloField4Measure,char Prefix[256],char Suffix[256]){
	ScalarField LengthOfFlow;
	char FlowLength[256];
	char FileName[256];
	
	CptLengthOfFlow(VeloField4Flow,VeloField4Measure,&LengthOfFlow);
	
	
	strcpy(FlowLength,Suffix);
	strcpy(FileName,Prefix);
	strcat(FileName,FlowLength);
	LengthOfFlow.Write(FileName,SourceFiles[0]);
}

///By following the flow defined by the velocity field 'VeloField4Flow' FOLLOW IN TIME the contribution of
///'VeloField4Measure' in the length of the flow from each point of the field.
/// * 'VeloField4Measure' is assumed to be part of a linear decomposition of 'VeloField4Flow'.
/// * If 'VeloField4Measure'=='VeloField4Flow' then the length of the flow defined by 'VeloField4Flow'
///   is computed.
void LargeDefGradLagrange::SaveEvoFlowLength(VectorField * VeloField4Flow,VectorField * VeloField4Measure,char Prefix[256],char Suffix[256]){
	ScalarField LengthOfFlow;
	char FlowLength[256];
	char FileName[256];
	
	CptEvoLengthOfFlow(VeloField4Flow,VeloField4Measure,&LengthOfFlow);
	
	
	strcpy(FlowLength,Suffix);
	strcpy(FileName,Prefix);
	strcat(FileName,FlowLength);
	LengthOfFlow.Write(FileName,SourceFiles[0]);
}




///save the map of the determinant of Jacobians
void LargeDefGradLagrange::SaveDetJacobian(char Prefix[256]){
	char FileName[256];
	char StrDetJacobians[256];
	
	//compute the determinant of jacobian
	CptMappingFromVeloField(this->VelocityField.NT-1,&this->MappingTrgImag,&this->VelocityField,&this->BackwardMapping);
	
	Cpt_JacobianDeterminant(&BackwardMapping,&DetJacobians,0);
	
	strcpy(StrDetJacobians,"_DetJacobian.nii");
	
	strcpy(FileName,Prefix);
	strcat(FileName,StrDetJacobians);
	DetJacobians.Write(FileName,SourceFiles[0]);
}


///save the map of the determinant of Jacobians for each kernel in the context of sum of kernels
void LargeDefGradLagrange::SaveSplittedDetJacobian(char Prefix[256]){
	VectorField Temp3DVecField;
	int k;
	char FileName[256];
	char StrDetJacobiansK1[256];
	char StrDetJacobiansK2[256];
	char StrDetJacobiansK3[256];
	char StrDetJacobiansK4[256];
	char StrDetJacobiansK5[256];
	char StrDetJacobiansK6[256];
	char StrDetJacobiansK7[256];
	
	strcpy(StrDetJacobiansK1,"_DetJacobianK1.nii");
	strcpy(StrDetJacobiansK2,"_DetJacobianK2.nii");
	strcpy(StrDetJacobiansK3,"_DetJacobianK3.nii");
	strcpy(StrDetJacobiansK4,"_DetJacobianK4.nii");
	strcpy(StrDetJacobiansK5,"_DetJacobianK5.nii");
	strcpy(StrDetJacobiansK6,"_DetJacobianK6.nii");
	strcpy(StrDetJacobiansK7,"_DetJacobianK7.nii");
	Temp3DVecField.CreateVoidField(this->NX, this->NY, this->NZ);
	
	//compute the determinant of jacobian
	for (k=0;k<this->NbKernels;k++){
		ComputeLagrangianPartialMapping(0,this->NbTimeSubdiv,&this->VelocityField,&this->SplittedVelocityField[k],&Temp3DVecField);
		Cpt_JacobianDeterminant(&Temp3DVecField,&DetJacobians,0);
		
		strcpy(FileName,Prefix);
		if (k==0) strcat(FileName,StrDetJacobiansK1);
		if (k==1) strcat(FileName,StrDetJacobiansK2);
		if (k==2) strcat(FileName,StrDetJacobiansK3);
		if (k==3) strcat(FileName,StrDetJacobiansK4);
		if (k==4) strcat(FileName,StrDetJacobiansK5);
		if (k==5) strcat(FileName,StrDetJacobiansK6);
		if (k==6) strcat(FileName,StrDetJacobiansK7);
		DetJacobians.Write(FileName,SourceFiles[0]);
	}
}




///++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
///                                      RUN FUNCTIONS
///++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


///Measure of the inverse typical amplitude of the deformations for given source
///and target images plus a simple Gaussian Kernel with a weight = 1.
float LargeDefGradLagrange::Run_MeasureTypicalAmplitude(void){
	float MaxGrad;
	double LocGrad;
	int x,y,z,i;
	
	//1) INITIALISATION
	
	//1.1) make sure that only a simple kernel with mono-channel images and no input velocity field will be used
	this->weight1=1.;
	this->weight2=0.;   this->sigmaX2=-1.; this->sigmaY2=-1.; this->sigmaZ2=-1.;
	this->weight3=0.;   this->sigmaX3=-1.; this->sigmaY3=-1.; this->sigmaZ3=-1.;
	this->weight4=0.;   this->sigmaX4=-1.; this->sigmaY4=-1.; this->sigmaZ4=-1.;
	this->weight5=0.;   this->sigmaX5=-1.; this->sigmaY5=-1.; this->sigmaZ5=-1.;
	this->weight6=0.;   this->sigmaX6=-1.; this->sigmaY6=-1.; this->sigmaZ6=-1.;
	this->weight7=0.;   this->sigmaX7=-1.; this->sigmaY7=-1.; this->sigmaZ7=-1.;
	this->NbKernels=1;
	this->SplitKernels=0;
	this->NbChannels=1;
	for (i=1;i<100;i++) strcpy(this->SourceFiles[i],"Null");
	for (i=1;i<100;i++) strcpy(this->TargetFiles[i],"Null");
	strcpy(this->PrefixInputs,"Null");
	
	//1.2) Pre-treatment of the inuput images (grey level alignment + margins)
	this->ReadAndTreatInputImages();
	
	//1.3) Allocations of the scalar and vector fields + definition of global parameters
	this->AllocateAllVariables();
	
	//1.4) Initiate the class to smooth the images
	FFTconvolver.InitiateConvolver(this->NX,this->NY,this->NZ,1.,this->sigmaX1,this->sigmaY1,this->sigmaZ1);
	
	//2) MEASURE OF THE TYPICAL AMPLITUDE
	
	//2.1) compute the forward mapping on space
	CptMappingFromVeloField(0,&this->MappingSrcImag,&this->VelocityField,&this->ForwardMapping);
	
	//2.2) compute the backward mapping on space
	CptMappingFromVeloField(this->VelocityField.NT-1,&this->MappingTrgImag,&this->VelocityField,&this->BackwardMapping);
	
	//2.3) LOOP ON THE TIME SUBDIVISIONS
	//2.3.1) compute the determinant of the jacobian of the transformation
	Cpt_JacobianDeterminant(&BackwardMapping,&DetJacobians,0);
	
	//2.3.2) compute the temporary image transformed using the forward mapping from time 0 -> J0
	Project3Dimage(&this->ImTemplate[0],&this->ForwardMapping,&this->J0,0);
	
	//2.3.3) compute the temporary image transformed using the backward mapping from time 1 -> J1
	Project3Dimage(&this->ImTarget[0],&this->BackwardMapping,&this->J1,0);
	
	//2.3.4) compute gradient of J0
	Cpt_Grad_ScalarField(&this->J0,&this->GradJ);
	
	//2.3.5) compute the gradient of energy
	this->ComputeEnergyGradient(0,0);
	
	//2.4) Compute the maximum of gradient in all time frames...
	//...3D images
	MaxGrad=0;
	for (z = 1; z < this->NZ-2; z++) for (y = 1; y < this->NY-2; y++) for (x = 1; x < this->NX-2; x++){
		LocGrad=sqrt(pow((double)this->GradE.G(0,x,y,z,0),2.)+pow((double)this->GradE.G(1,x,y,z,0),2.)+pow((double)this->GradE.G(2,x,y,z,0),2.));
		if (MaxGrad<(float)LocGrad) MaxGrad=(float)LocGrad;
	}
	
	//...2D images
	if (this->NZ==1){
		for (y = 1; y < this->NY-2; y++) for (x = 1; x < this->NX-2; x++){
			LocGrad=sqrt(pow((double)this->GradE.G(0,x,y,0,0),2.)+pow((double)this->GradE.G(1,x,y,0,0),2.));
			if (MaxGrad<(float)LocGrad) MaxGrad=(float)LocGrad;
		}
	}
	
	
	//3) RETURN THE INVERSE TYPICAL AMPLITUDE
	// the *10000 is only here for readability convenience
	
	//cout << "Typical update = " << MaxGrad/10000. << "\n";
	cout << "Inverse typical update = " << 10000./MaxGrad << "\n";
	cout << "-> Use this value as the weight of a Gaussian kernel to have the same apparent weights everywhere.\n";
	
	return 1./MaxGrad;
}

///Function to solve the registration using the gradient descent algorithm of Beg 05
void LargeDefGradLagrange::Run_Default(void){
	int IterationStopper;
	int IterationNb;
	int TimeSubdiv;
	int IdChannel;
	float SqrtSSD;
	float NormaMaxGrad;  //[maximum gradient at the current iteration] / [maximum gradient at the first iteration]
	float PreviousNormaMaxGrad[7];
	int i;
	
	//1) INITIALISATION
	
	//1.1) Pre-treatment of the inuput images (grey level alignment + margins)
	this->ReadAndTreatInputImages();
	
	//1.2) Allocations of the scalar and vector fields + definition of global parameters
	this->AllocateAllVariables();
	
	//1.3) Initiate the class to smooth the images
	FFTconvolver.InitiateConvolver(this->NX,this->NY,this->NZ,this->weight1,this->sigmaX1,this->sigmaY1,this->sigmaZ1,this->weight2,this->sigmaX2,this->sigmaY2,this->sigmaZ2,this->weight3,this->sigmaX3,this->sigmaY3,this->sigmaZ3,this->weight4,this->sigmaX4,this->sigmaY4,this->sigmaZ4,this->weight5,this->sigmaX5,this->sigmaY5,this->sigmaZ5,this->weight6,this->sigmaX6,this->sigmaY6,this->sigmaZ6,this->weight7,this->sigmaX7,this->sigmaY7,this->sigmaZ7);
	
	
	//2) GRADIENT DESCENT
	if (this->iteration_nb==0) IterationStopper=1;
	else IterationStopper=0;
	IterationNb=0;
	for (i=0;i<7;i++) PreviousNormaMaxGrad[i]=1.;
	while (IterationStopper==0){
		cout << "Iteration Number " << IterationNb+1 << " / " << this->iteration_nb << "\n";
		
		//2.1) compute the forward mapping on space
		CptMappingFromVeloField(0,&this->MappingSrcImag,&this->VelocityField,&this->ForwardMapping);
		
		//2.2) compute the backward mapping on space
		CptMappingFromVeloField(this->VelocityField.NT-1,&this->MappingTrgImag,&this->VelocityField,&this->BackwardMapping);
		
		//2.3) eventually also computes the mapping from t=0.5
		if (this->symmetric==1) CptMappingFromVeloField((this->VelocityField.NT-1)/2,&this->MappingId,&this->VelocityField,&this->MappingFromT05);
		
		//2.3) LOOP ON THE TIME SUBDIVISIONS AND THE CHANNELS
		for (TimeSubdiv=0;TimeSubdiv<this->NbTimeSubdiv;TimeSubdiv++){//LOOP ON THE TIME SUBDIVISIONS
			
			//2.3.1) compute the determinant of the jacobian of the transformation
			if (this->symmetric==0) Cpt_JacobianDeterminant(&BackwardMapping,&DetJacobians,TimeSubdiv);
			else  Cpt_JacobianDeterminant(&MappingFromT05,&DetJacobians,TimeSubdiv);
			
			for (IdChannel=0;IdChannel<this->NbChannels;IdChannel++){//LOOP ON THE CHANNELS
				//2.3.2) compute the temporary image transformed using the forward mapping from time 0 -> J0
				Project3Dimage(&this->ImTemplate[IdChannel],&this->ForwardMapping,&this->J0,TimeSubdiv);
				
				//2.3.3) compute the temporary image transformed using the backward mapping from time 1 -> J1
				Project3Dimage(&this->ImTarget[IdChannel],&this->BackwardMapping,&this->J1,TimeSubdiv);
				
				//2.3.4) compute gradient of J
				if (this->symmetric==0)  Cpt_Grad_ScalarField(&this->J0,&this->GradJ);
				else{
					if(TimeSubdiv<=(this->NbTimeSubdiv-1)/2)  Cpt_Grad_ScalarField(&this->J0,&this->GradJ);
					else  Cpt_Grad_ScalarField(&this->J1,&this->GradJ);
				}
				
				//2.3.5) compute the gradient of energy
				this->ComputeEnergyGradient(TimeSubdiv,IdChannel);
				
				//2.3.6) Measure the convergence of the similarity between the source and target image
				if (ShowSSD==1) if (TimeSubdiv==this->NbTimeSubdiv-1){
					SqrtSSD=CalcSqrtSumOfSquaredDif(&this->J0,&this->J1);
					cout << "Sqrt SSD = " << SqrtSSD << " (Channel " << IdChannel << ")\n";
				}
			}
		}
		
		//2.4) update the velocity fields
		NormaMaxGrad=this->UpdateVelocityField(IterationNb);
		
		//2.5) end of the convergence...
		//...controled by the number of iterations
		IterationNb++;
		if (IterationNb>=this->iteration_nb) IterationStopper=1;
		
		//...controled by the maximum gradient
		if ((NormaMaxGrad<epsilon)&&(NormaMaxGrad>PreviousNormaMaxGrad[0])&&(NormaMaxGrad>PreviousNormaMaxGrad[2])&&(NormaMaxGrad>PreviousNormaMaxGrad[4])&&(NormaMaxGrad>PreviousNormaMaxGrad[6]))  IterationStopper=1;
		PreviousNormaMaxGrad[6]=PreviousNormaMaxGrad[5];
		PreviousNormaMaxGrad[5]=PreviousNormaMaxGrad[4];
		PreviousNormaMaxGrad[4]=PreviousNormaMaxGrad[3];
		PreviousNormaMaxGrad[3]=PreviousNormaMaxGrad[2];
		PreviousNormaMaxGrad[2]=PreviousNormaMaxGrad[1];
		PreviousNormaMaxGrad[1]=PreviousNormaMaxGrad[0];
		PreviousNormaMaxGrad[0]=NormaMaxGrad;
	}
	
	//3) SAVE THE RESULTS
	this->SaveResultGradientDescent();
}



///++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
///                                        MAIN RUN FUNCTION 
///++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


///run function
void LargeDefGradLagrange::Run(void)
{
	if (this->MeasureTypicAmp!=1)
		this->Run_Default();
	else
		this->Run_MeasureTypicalAmplitude();
}

