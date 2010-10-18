/*=========================================================================

  Date      : $Date: 03.03.2009$
  Changes   : $Author: Laurent Risser $

=========================================================================*/

#include <irtkImage.h>
#include <irtkAnisoDiffusion.h>


template <class VoxelType> anisoDiffusion<VoxelType>::anisoDiffusion(){
  //default parameters
  ax=3;
  ay=3;
  az=3;
  at=3;
  dx=1;
  dy=1;
  dz=1;
  dt=1;
  dTau=1;
  ITERATIONS_NB=5;
  TimeDependent=false;
  SemiImplicit=true;
}

template <class VoxelType> anisoDiffusion<VoxelType>::~anisoDiffusion(void)
{}

template <class VoxelType> bool anisoDiffusion<VoxelType>::RequiresBuffering(void)
{
  return true;
}

template <class VoxelType> const char *anisoDiffusion<VoxelType>::NameOfClass()
{
  return "anisoDiffusion";
}


template <class VoxelType> void anisoDiffusion<VoxelType>::Run_3D_semiImplicit(){
	int i, j, x, y, z, t;
	double  ax,ay,az,at,dx,dy,dz,dt;
	float dTau;
	float*** imageE;
	float*** imageO;
	int NBX,NBY,NBZ,NBT;
	double dIdx,dIdy,dIdz;
	float DivDgradI;
	float *Va;
	float *Vb; 
	float *Vc;
	float *Vd;
	float *Vx;
	int n;
	int ITERATIONS_NB;
	float Dxx_div_dxSq,Dyy_div_dySq,Dzz_div_dzSq;
	int iteration;
	float DivPowDxSqu,DivPowDySqu,DivPowDzSqu,DivPowDtSqu;
	
	//1) INITIALISATION
	
	// Do the initial set up
	this->Initialize();
	
	//variables definition
	ax=this->ax;
	ay=this->ay;
	az=this->az;
	at=this->at;
	dx=this->dx;
	dy=this->dy;
	dz=this->dz;
	dt=this->dt;
	dTau=this->dTau;
	ITERATIONS_NB=this->ITERATIONS_NB;
	NBX=this->_input->GetX()+2;  //for boundary effects
	NBY=this->_input->GetY()+2;  //for boundary effects
	NBZ=this->_input->GetZ()+2;  //for boundary effects
	NBT=this->_input->GetT();
	cout << "Image size: " << (NBX-2) <<  " , "  <<  (NBY-2)  <<  " , "  << (NBZ-2)  <<  " , " << NBT  << " + boundaries \n";
	
	//temporary input and output images
	imageE= (float***) malloc (NBZ*sizeof(float**));
	for (i=0;i<NBZ;i++) imageE[i]= (float**) malloc (NBY*sizeof(float*));
	for (i=0;i<NBZ;i++) for (j=0;j<NBY;j++) imageE[i][j]= (float*) malloc (NBX*sizeof(float));
	
	imageO= (float***) malloc (NBZ*sizeof(float**));
	for (i=0;i<NBZ;i++) imageO[i]= (float**) malloc (NBY*sizeof(float*));
	for (i=0;i<NBZ;i++) for (j=0;j<NBY;j++) imageO[i][j]= (float*) malloc (NBX*sizeof(float));
	
	//precomputed values
	DivPowDxSqu=1./pow(dx,2);
	DivPowDySqu=1./pow(dy,2);
	DivPowDzSqu=1./pow(dz,2);
	DivPowDtSqu=1./pow(dt,2);
	
	//temporary variables dedicated to the semi implicit scheme (+4 is to avoid boundary effects)
	n=max(max(max(NBX,NBY),NBZ),NBT)+4; //for boundary effects
	Va=(float*)malloc(n*sizeof(float));
	Vb=(float*)malloc(n*sizeof(float));
	Vc=(float*)malloc(n*sizeof(float));
	Vd=(float*)malloc(n*sizeof(float));
	Vx=(float*)malloc(n*sizeof(float));

	//2) ANISOTROPIC DIFFUSION
	for (t = 0; t < NBT; t++) {
		cout << "Image " << t+1 << " / " << NBT << "\n";
		//2.1) convert the values of the input image at time t in double in a temporary 3D image
		for (z = 0; z < NBZ-2; z++)  for (y = 0; y < NBY-2; y++) for (x = 0; x < NBX-2; x++)
			imageE[z+1][y+1][x+1]=static_cast<float>(this->_input->Get(x, y, z, t));
		
		for (z = 0; z < NBZ-2; z++)  for (y = 0; y < NBY-2; y++) for (x = 0; x < NBX-2; x++)
			imageO[z+1][y+1][x+1]=static_cast<float>(this->_input->Get(x, y, z, t));
		
		//image extension to avoid boundary effects
		for (z = 0; z < NBZ; z++) for (y = 0; y < NBY; y++) imageE[z][y][0]=imageE[z][y][1];
		for (z = 0; z < NBZ; z++) for (y = 0; y < NBY; y++) imageE[z][y][NBX-1]=imageE[z][y][NBX-2];
		for (z = 0; z < NBZ; z++) for (y = 0; y < NBY; y++) imageO[z][y][0]=imageO[z][y][1];
		for (z = 0; z < NBZ; z++) for (y = 0; y < NBY; y++) imageO[z][y][NBX-1]=imageO[z][y][NBX-2];
		
		for (z = 0; z < NBZ; z++) for (x = 0; x < NBX; x++) imageE[z][0][x]=imageE[z][1][x];
		for (z = 0; z < NBZ; z++) for (x = 0; x < NBX; x++) imageE[z][NBY-1][x]=imageE[z][NBY-2][x];
		for (z = 0; z < NBZ; z++) for (x = 0; x < NBX; x++) imageO[z][0][x]=imageO[z][1][x];
		for (z = 0; z < NBZ; z++) for (x = 0; x < NBX; x++) imageO[z][NBY-1][x]=imageO[z][NBY-2][x];
		
		for (y = 0; y < NBY; y++) for (x = 0; x < NBX; x++) imageE[0][y][x]=imageE[1][y][x];
		for (y = 0; y < NBY; y++) for (x = 0; x < NBX; x++) imageE[NBZ-1][y][x]=imageE[NBZ-2][y][x];
		for (y = 0; y < NBY; y++) for (x = 0; x < NBX; x++) imageO[0][y][x]=imageO[1][y][x];
		for (y = 0; y < NBY; y++) for (x = 0; x < NBX; x++) imageO[NBZ-1][y][x]=imageO[NBZ-2][y][x];
		
		//2.2) diffusion in the temporary 3D image - ADI semi implicit scheme
		for (iteration=0 ; iteration<ITERATIONS_NB; iteration++){
			cout << "| Iteration " << iteration+1 << " / " << ITERATIONS_NB << "\n";
			
			//2.2.2) diffusion - x implicit / y,z explicit
			//2.2.2.1) explicit part
			for (z = 1; z < NBZ-1; z++) for (y = 1; y < NBY-1; y++) for (x = 1; x < NBX-1; x++) {
				dIdy=(imageE[z][y+1][x]-imageE[z][y-1][x])/(2*dy);
				Dyy_div_dySq=static_cast<float>((1-exp(-3.314/pow((dIdy/ay),4)))*DivPowDySqu);
				dIdz=(imageE[z+1][y][x]-imageE[z-1][y][x])/(2*dz);
				Dzz_div_dzSq=static_cast<float>((1-exp(-3.314/pow((dIdz/az),4)))*DivPowDzSqu);
				//new value of the voxel
				DivDgradI=(imageE[z][y+1][x]-2*imageE[z][y][x]+imageE[z][y-1][x])*Dyy_div_dySq+
					(imageE[z+1][y][x]-2*imageE[z][y][x]+imageE[z-1][y][x])*Dzz_div_dzSq;
				
				imageO[z][y][x]=imageE[z][y][x]+(dTau/3.)*DivDgradI;
			}
				
			//2.2.2.2) implicit part
			for (z = 1; z < NBZ-1; z++) for (y = 1; y < NBY-1; y++) {
				for (x = 1; x < NBX-1; x++){
					dIdx=(imageE[z][y][x+1]-imageE[z][y][x-1])/(2*dx);
					Dxx_div_dxSq=static_cast<float>((1-exp(-3.314/pow((dIdx/ax),4)))*DivPowDxSqu);
					Va[x+1]=(dTau/3.)*Dxx_div_dxSq;
					Vb[x+1]=-1-2*(dTau/3.)*Dxx_div_dxSq;
					Vc[x+1]=(dTau/3.)*Dxx_div_dxSq;
                                        Vd[x+1]=imageE[z][y][x]; //why not imageO ???
				}
				Va[1]=Va[3]; Va[0]=Va[4]; Va[NBX]=Va[NBX-2]; Va[NBX+1]=Va[NBX-3]; //to avoid boundary effects
				Vb[1]=Vb[3]; Vb[0]=Vb[4]; Vb[NBX]=Vb[NBX-2]; Vb[NBX+1]=Vb[NBX-3]; //to avoid boundary effects
				Vc[1]=Vc[3]; Vc[0]=Vc[4]; Vc[NBX]=Vc[NBX-2]; Vc[NBX+1]=Vc[NBX-3]; //to avoid boundary effects
				Vd[1]=Vd[3]; Vd[0]=Vd[4]; Vd[NBX]=Vd[NBX-2]; Vd[NBX+1]=Vd[NBX-3]; //to avoid boundary effects
				TridiagonalSolveFloat(Va,Vb,Vc,Vd,Vx,NBX+2);
				for (x = 1; x < NBX-1; x++) imageO[z][y][x]=-Vx[x+1];
			}
			
			//2.2.3) diffusion - y implicit / x,z explicit
			//2.2.3.1) explicit part
			for (z = 1; z < NBZ-1; z++) for (y = 1; y < NBY-1; y++) for (x = 1; x < NBX-1; x++) {
				dIdx=(imageO[z][y][x+1]-imageO[z][y][x-1])/(2*dx);
				Dxx_div_dxSq=static_cast<float>((1-exp(-3.314/pow((dIdx/ax),4)))*DivPowDxSqu);
				dIdz=(imageO[z+1][y][x]-imageO[z-1][y][x])/(2*dz);
				Dzz_div_dzSq=static_cast<float>((1-exp(-3.314/pow((dIdz/az),4)))*DivPowDzSqu);
				//new value of the voxel
				DivDgradI=(imageO[z][y][x+1]-2*imageO[z][y][x]+imageO[z][y][x-1])*Dxx_div_dxSq+
						(imageO[z+1][y][x]-2*imageO[z][y][x]+imageO[z-1][y][x])*Dzz_div_dzSq;
				
				imageE[z][y][x]=imageO[z][y][x]+(dTau/3.)*DivDgradI;
			}
			
			//2.2.3.2) implicit part
			for (z = 1; z < NBZ-1; z++) for (x = 1; x < NBX-1; x++){
				for (y = 1; y < NBY-1; y++){
					dIdy=(imageO[z][y+1][x]-imageO[z][y-1][x])/(2*dy);
					Dyy_div_dySq=static_cast<float>((1-exp(-3.314/pow((dIdy/ay),4)))*DivPowDySqu);
					Va[y+1]=(dTau/3.)*Dyy_div_dySq;
					Vb[y+1]=-1-2*(dTau/3.)*Dyy_div_dySq;
					Vc[y+1]=(dTau/3.)*Dyy_div_dySq;
					Vd[y+1]=imageO[z][y][x];
				}
				Va[1]=Va[3]; Va[0]=Va[4]; Va[NBY]=Va[NBY-2]; Va[NBY+1]=Va[NBY-3]; //to avoid boundary effects
				Vb[1]=Vb[3]; Vb[0]=Vb[4]; Vb[NBY]=Vb[NBY-2]; Vb[NBY+1]=Vb[NBY-3]; //to avoid boundary effects
				Vc[1]=Vc[3]; Vc[0]=Vc[4]; Vc[NBY]=Vc[NBY-2]; Vc[NBY+1]=Vc[NBY-3]; //to avoid boundary effects
				Vd[1]=Vd[3]; Vd[0]=Vd[4]; Vd[NBY]=Vd[NBY-2]; Vd[NBY+1]=Vd[NBY-3]; //to avoid boundary effects
				TridiagonalSolveFloat(Va,Vb,Vc,Vd,Vx,NBY+2);
				for (y = 1; y < NBY-1; y++) imageE[z][y][x]=-Vx[y+1];
			}
		
			//2.2.4) diffusion - z implicit / x,y explicit
			//2.2.4.1) explicit part
			for (z = 1; z < NBZ-1; z++) for (y = 1; y < NBY-1; y++) for (x = 1; x < NBX-1; x++) {
				dIdx=(imageE[z][y][x+1]-imageE[z][y][x-1])/(2*dx);
				Dxx_div_dxSq=static_cast<float>((1-exp(-3.314/pow((dIdx/ax),4)))*DivPowDxSqu);
				dIdy=(imageE[z][y+1][x]-imageE[z][y-1][x])/(2*dy);
				Dyy_div_dySq=static_cast<float>((1-exp(-3.314/pow((dIdy/ay),4)))*DivPowDySqu);
				//new value of the voxel
				DivDgradI=(imageE[z][y][x+1]-2*imageE[z][y][x]+imageE[z][y][x-1])*Dxx_div_dxSq+
					(imageE[z][y+1][x]-2*imageE[z][y][x]+imageE[z][y-1][x])*Dyy_div_dySq;
					
				imageO[z][y][x]=imageE[z][y][x]+(dTau/3.)*DivDgradI;
			}
			
			//2.2.4.2) implicit part
			for (y = 1; y < NBY-1; y++) for (x = 1; x < NBX-1; x++) {
				for (z = 1; z < NBZ-1; z++){
					dIdz=(imageE[z+1][y][x]-imageE[z-1][y][x])/(2*dz);
					Dzz_div_dzSq=static_cast<float>((1-exp(-3.314/pow((dIdz/az),4)))*DivPowDzSqu);
					Va[z+1]=(dTau/3.)*Dzz_div_dzSq;
					Vb[z+1]=-1-2*(dTau/3.)*Dzz_div_dzSq;
					Vc[z+1]=(dTau/3.)*Dzz_div_dzSq;
					Vd[z+1]=imageE[z][y][x];
				}
				Va[1]=Va[3]; Va[0]=Va[4]; Va[NBZ]=Va[NBZ-2]; Va[NBZ+1]=Va[NBZ-3]; //to avoid boundary effects
				Vb[1]=Vb[3]; Vb[0]=Vb[4]; Vb[NBZ]=Vb[NBZ-2]; Vb[NBZ+1]=Vb[NBZ-3]; //to avoid boundary effects
				Vc[1]=Vc[3]; Vc[0]=Vc[4]; Vc[NBZ]=Vc[NBZ-2]; Vc[NBZ+1]=Vc[NBZ-3]; //to avoid boundary effects
				Vd[1]=Vd[3]; Vd[0]=Vd[4]; Vd[NBZ]=Vd[NBZ-2]; Vd[NBZ+1]=Vd[NBZ-3]; //to avoid boundary effects
				TridiagonalSolveFloat(Va,Vb,Vc,Vd,Vx,NBZ+2);
				for (z = 1; z < NBZ-1; z++) imageO[z][y][x]=-Vx[z+1];
			}
			
			/*for (y = 1; y < NBY-1; y++) for (x = 1; x < NBX-1; x++) {
				for (z = 1; z < NBZ-1; z++){
					dIdz=(imageE[z+1][y][x]-imageE[z-1][y][x])/(2*dz);
					Dzz_div_dzSq=static_cast<float>((1-exp(-3.314/pow((dIdz/az),4)))*DivPowDzSqu);
					Va[z-1]=(dTau/3.)*Dzz_div_dzSq;
					Vb[z-1]=-1-2*(dTau/3.)*Dzz_div_dzSq;
					Vc[z-1]=(dTau/3.)*Dzz_div_dzSq;
					Vd[z-1]=imageE[z][y][x];
				}
				TridiagonalSolveFloat(Va,Vb,Vc,Vd,Vx,NBZ-2);
				for (z = 1; z < NBZ-1; z++) imageO[z][y][x]=-Vx[z-1];
			}*/
		
			//2.2.5) temporary output image is reinjected in temporary input image
			for (z = 0; z < NBZ; z++) for (y = 0; y < NBY; y++) for (x = 0; x < NBX; x++)
				imageE[z][y][x]=imageO[z][y][x];
				
		}
		//2.3) save the filtered temporary 3D image in VoxelType in the  output image at time t
		for (z = 0; z < NBZ-2; z++) for (y = 0; y < NBY-2; y++) for (x = 0; x < NBX-2; x++) 
			this->_output->Put(x, y, z, t, static_cast<VoxelType>(imageE[z+1][y+1][x+1]));
	}
	
	//3) END OF THE FUNCTION
	// Do the final cleaning up
	this->Finalize();
}



template <class VoxelType> void anisoDiffusion<VoxelType>::Run_4D_semiImplicit(){
	int i, j, x, y, z, t;
	double  ax,ay,az,at,dx,dy,dz,dt;
	float dTau;
	float**** imageE;
	float**** imageO;
	int NBX,NBY,NBZ,NBT;
	double dIdx,dIdy,dIdz,dIdt;
	float DivDgradI;
	float *Va;
	float *Vb; 
	float *Vc;
	float *Vd;
	float *Vx;
	int n;
	int ITERATIONS_NB;
	float Dxx_div_dxSq,Dyy_div_dySq,Dzz_div_dzSq,Dtt_div_dtSq;
	int iteration;
	float DivPowDxSqu,DivPowDySqu,DivPowDzSqu,DivPowDtSqu;
	
	//1) INITIALISATION
	
	// Do the initial set up
	this->Initialize();
	
	//variables definition
	ax=this->ax;
	ay=this->ay;
	az=this->az;
	at=this->at;
	dx=this->dx;
	dy=this->dy;
	dz=this->dz;
	dt=this->dt;
	dTau=this->dTau;
	ITERATIONS_NB=this->ITERATIONS_NB;
	NBX=this->_input->GetX()+2;
	NBY=this->_input->GetY()+2;
	NBZ=this->_input->GetZ()+2;
	NBT=this->_input->GetT()+2;
	cout << "Image size: " << (NBX-2) <<  " , "  <<  (NBY-2)  <<  " , "  << (NBZ-2)  <<  " , " << (NBT-2)  << " + boundaries \n";
	
	//temporary input and output images
	imageE= (float****) malloc (NBT*sizeof(float***));
	for (t=0;t<NBT;t++) imageE[t]= (float***) malloc (NBZ*sizeof(float**));
	for (t=0;t<NBT;t++) for (i=0;i<NBZ;i++) imageE[t][i]= (float**) malloc (NBY*sizeof(float*));
	for (t=0;t<NBT;t++) for (i=0;i<NBZ;i++) for (j=0;j<NBY;j++) imageE[t][i][j]= (float*) malloc (NBX*sizeof(float));
	
	imageO= (float****) malloc (NBT*sizeof(float***));
	for (t=0;t<NBT;t++) imageO[t]= (float***) malloc (NBZ*sizeof(float**));
	for (t=0;t<NBT;t++) for (i=0;i<NBZ;i++) imageO[t][i]= (float**) malloc (NBY*sizeof(float*));
	for (t=0;t<NBT;t++) for (i=0;i<NBZ;i++) for (j=0;j<NBY;j++) imageO[t][i][j]= (float*) malloc (NBX*sizeof(float));
	
	//precomputed values
	DivPowDxSqu=1./pow(dx,2);
	DivPowDySqu=1./pow(dy,2);
	DivPowDzSqu=1./pow(dz,2);
	DivPowDtSqu=1./pow(dt,2);
	
	//temporary variables dedicated to the semi implicit scheme
	n=max(max(max(NBX,NBY),NBZ),NBT)+4; //for boundary effects
	Va=(float*)malloc(n*sizeof(float));
	Vb=(float*)malloc(n*sizeof(float));
	Vc=(float*)malloc(n*sizeof(float));
	Vd=(float*)malloc(n*sizeof(float));
	Vx=(float*)malloc(n*sizeof(float));

	//2) ANISOTROPIC DIFFUSION
		
	//2.1) convert the values of the input image at time t in double in a temporary 3D image
	for (t = 0; t < NBT-2; t++) for (z = 0; z < NBZ-2; z++)  for (y = 0; y < NBY-2; y++) for (x = 0; x < NBX-2; x++)
		imageE[t+1][z+1][y+1][x+1]=static_cast<float>(this->_input->Get(x, y, z, t));
	
	for (t = 0; t < NBT-2; t++) for (z = 0; z < NBZ-2; z++)  for (y = 0; y < NBY-2; y++) for (x = 0; x < NBX-2; x++)
		imageO[t+1][z+1][y+1][x+1]=static_cast<float>(this->_input->Get(x, y, z, t));
	
	for (t = 0; t < NBT; t++)  for (z = 0; z < NBZ; z++) for (y = 0; y < NBY; y++) imageE[t][z][y][0]=imageE[t][z][y][1];
	for (t = 0; t < NBT; t++)  for (z = 0; z < NBZ; z++) for (y = 0; y < NBY; y++) imageE[t][z][y][NBX-1]=imageE[t][z][y][NBX-2];
	for (t = 0; t < NBT; t++)  for (z = 0; z < NBZ; z++) for (y = 0; y < NBY; y++) imageO[t][z][y][0]=imageO[t][z][y][1];
	for (t = 0; t < NBT; t++)  for (z = 0; z < NBZ; z++) for (y = 0; y < NBY; y++) imageO[t][z][y][NBX-1]=imageO[t][z][y][NBX-2];
	
	for (t = 0; t < NBT; t++)  for (z = 0; z < NBZ; z++) for (x = 0; x < NBX; x++) imageE[t][z][0][x]=imageE[t][z][1][x];
	for (t = 0; t < NBT; t++)  for (z = 0; z < NBZ; z++) for (x = 0; x < NBX; x++) imageE[t][z][NBY-1][x]=imageE[t][z][NBY-2][x];
	for (t = 0; t < NBT; t++)  for (z = 0; z < NBZ; z++) for (x = 0; x < NBX; x++) imageO[t][z][0][x]=imageO[t][z][1][x];
	for (t = 0; t < NBT; t++)  for (z = 0; z < NBZ; z++) for (x = 0; x < NBX; x++) imageO[t][z][NBY-1][x]=imageO[t][z][NBY-2][x];
	
	
	for (t = 0; t < NBT; t++)  for (y = 0; y < NBY; y++) for (x = 0; x < NBX; x++) imageE[t][0][y][x]=imageE[t][1][y][x];
	for (t = 0; t < NBT; t++)  for (y = 0; y < NBY; y++) for (x = 0; x < NBX; x++) imageE[t][NBZ-1][y][x]=imageE[t][NBZ-2][y][x];
	for (t = 0; t < NBT; t++)  for (y = 0; y < NBY; y++) for (x = 0; x < NBX; x++) imageO[t][0][y][x]=imageO[t][1][y][x];
	for (t = 0; t < NBT; t++)  for (y = 0; y < NBY; y++) for (x = 0; x < NBX; x++) imageO[t][NBZ-1][y][x]=imageO[t][NBZ-2][y][x];

	for (z = 0; z < NBZ; z++)  for (y = 0; y < NBY; y++) for (x = 0; x < NBX; x++) imageE[0][z][y][x]=imageE[0][z][y][x];
	for (z = 0; z < NBZ; z++)  for (y = 0; y < NBY; y++) for (x = 0; x < NBX; x++) imageE[NBT-1][z][y][x]=imageE[NBT-1][z][y][x];
	for (z = 0; z < NBZ; z++)  for (y = 0; y < NBY; y++) for (x = 0; x < NBX; x++) imageO[0][z][y][x]=imageO[0][z][y][x];
	for (z = 0; z < NBZ; z++)  for (y = 0; y < NBY; y++) for (x = 0; x < NBX; x++) imageO[NBT-1][z][y][x]=imageO[NBT-1][z][y][x];

	
	//2.2) DIFFUSION
	//2.2) diffusion in the temporary 4D image - ADI semi implicit scheme
	for (iteration=0 ; iteration<ITERATIONS_NB; iteration++){
		cout << "| Iteration " << iteration+1 << " / " << ITERATIONS_NB << "\n";
		
		//2.2.1) diffusion - x implicit / y,z,t explicit
		//2.2.1.1) explicit part
		for (t = 1; t < NBT-1; t++) for (z = 1; z < NBZ-1; z++) for (y = 1; y < NBY-1; y++) for (x = 1; x < NBX-1; x++) {
			dIdy=(imageE[t][z][y+1][x]-imageE[t][z][y-1][x])/(2*dy);
			Dyy_div_dySq=static_cast<float>((1-exp(-3.314/pow((dIdy/ay),4)))*DivPowDySqu);
			dIdz=(imageE[t][z+1][y][x]-imageE[t][z-1][y][x])/(2*dz);
			Dzz_div_dzSq=static_cast<float>((1-exp(-3.314/pow((dIdz/az),4)))*DivPowDzSqu);
			dIdt=(imageE[t+1][z][y][x]-imageE[t-1][z][y][x])/(2*dt);
			Dtt_div_dtSq=static_cast<float>((1-exp(-3.314/pow((dIdt/at),4)))*DivPowDzSqu);
			//new value of the voxel
			DivDgradI=(imageE[t][z][y+1][x]-2*imageE[t][z][y][x]+imageE[t][z][y-1][x])*Dyy_div_dySq+
					(imageE[t][z+1][y][x]-2*imageE[t][z][y][x]+imageE[t][z-1][y][x])*Dzz_div_dzSq+
					(imageE[t+1][z][y][x]-2*imageE[t][z][y][x]+imageE[t-1][z][y][x])*Dtt_div_dtSq;

			imageO[t][z][y][x]=imageE[t][z][y][x]+(dTau/4.)*DivDgradI;
		}
		
		//2.2.1.2) implicit part
		for (t = 1; t < NBT-1; t++) for (z = 1; z < NBZ-1; z++) for (y = 1; y < NBY-1; y++) {
			for (x = 1; x < NBX-1; x++){
				dIdx=(imageE[t][z][y][x+1]-imageE[t][z][y][x-1])/(2*dx);
				Dxx_div_dxSq=static_cast<float>((1-exp(-3.314/pow((dIdx/ax),4)))*DivPowDxSqu);
				Va[x+1]=(dTau/4.)*Dxx_div_dxSq;
				Vb[x+1]=-1-2*(dTau/4.)*Dxx_div_dxSq;
				Vc[x+1]=(dTau/4.)*Dxx_div_dxSq;
				Vd[x+1]=imageE[t][z][y][x];
			}
			Va[1]=Va[3]; Va[0]=Va[4]; Va[NBX]=Va[NBX-2]; Va[NBX+1]=Va[NBX-3]; //to avoid boundary effects
			Vb[1]=Vb[3]; Vb[0]=Vb[4]; Vb[NBX]=Vb[NBX-2]; Vb[NBX+1]=Vb[NBX-3]; //to avoid boundary effects
			Vc[1]=Vc[3]; Vc[0]=Vc[4]; Vc[NBX]=Vc[NBX-2]; Vc[NBX+1]=Vc[NBX-3]; //to avoid boundary effects
			Vd[1]=Vd[3]; Vd[0]=Vd[4]; Vd[NBX]=Vd[NBX-2]; Vd[NBX+1]=Vd[NBX-3]; //to avoid boundary effects
			TridiagonalSolveFloat(Va,Vb,Vc,Vd,Vx,NBX+2);
			for (x = 1; x < NBX-1; x++) imageO[t][z][y][x]=-Vx[x+1];
		}
		
		//2.2.2) diffusion - y implicit / x,z,t explicit
		//2.2.2.1) explicit part
		for (t = 1; t < NBT-1; t++) for (z = 1; z < NBZ-1; z++) for (y = 1; y < NBY-1; y++) for (x = 1; x < NBX-1; x++) {
			dIdx=(imageO[t][z][y][x+1]-imageO[t][z][y][x-1])/(2*dx);
			Dxx_div_dxSq=static_cast<float>((1-exp(-3.314/pow((dIdx/ax),4)))*DivPowDxSqu);
			dIdz=(imageO[t][z+1][y][x]-imageO[t][z-1][y][x])/(2*dz);
			Dzz_div_dzSq=static_cast<float>((1-exp(-3.314/pow((dIdz/az),4)))*DivPowDzSqu);
			dIdt=(imageO[t+1][z][y][x]-imageO[t-1][z][y][x])/(2*dt);
			Dtt_div_dtSq=static_cast<float>((1-exp(-3.314/pow((dIdt/at),4)))*DivPowDzSqu);
			//new value of the voxel
			DivDgradI=(imageO[t][z][y][x+1]-2*imageO[t][z][y][x]+imageO[t][z][y][x-1])*Dxx_div_dxSq+
					(imageO[t][z+1][y][x]-2*imageO[t][z][y][x]+imageO[t][z-1][y][x])*Dzz_div_dzSq+
					(imageO[t+1][z][y][x]-2*imageO[t][z][y][x]+imageO[t-1][z][y][x])*Dtt_div_dtSq;
			
			imageE[t][z][y][x]=imageO[t][z][y][x]+(dTau/4.)*DivDgradI;
		}
		
		//2.2.2.2) implicit part
		for (t = 1; t < NBT-1; t++) for (z = 1; z < NBZ-1; z++) for (x = 1; x < NBX-1; x++){
			for (y = 1; y < NBY-1; y++){
				dIdy=(imageO[t][z][y+1][x]-imageO[t][z][y-1][x])/(2*dy);
				Dyy_div_dySq=static_cast<float>((1-exp(-3.314/pow((dIdy/ay),4)))*DivPowDySqu);
				Va[y+1]=(dTau/4.)*Dyy_div_dySq;
				Vb[y+1]=-1-2*(dTau/4.)*Dyy_div_dySq;
				Vc[y+1]=(dTau/4.)*Dyy_div_dySq;
				Vd[y+1]=imageO[t][z][y][x];
			}
			Va[1]=Va[3]; Va[0]=Va[4]; Va[NBY]=Va[NBY-2]; Va[NBY+1]=Va[NBY-3]; //to avoid boundary effects
			Vb[1]=Vb[3]; Vb[0]=Vb[4]; Vb[NBY]=Vb[NBY-2]; Vb[NBY+1]=Vb[NBY-3]; //to avoid boundary effects
			Vc[1]=Vc[3]; Vc[0]=Vc[4]; Vc[NBY]=Vc[NBY-2]; Vc[NBY+1]=Vc[NBY-3]; //to avoid boundary effects
			Vd[1]=Vd[3]; Vd[0]=Vd[4]; Vd[NBY]=Vd[NBY-2]; Vd[NBY+1]=Vd[NBY-3]; //to avoid boundary effects
			TridiagonalSolveFloat(Va,Vb,Vc,Vd,Vx,NBY+2);
			for (y = 1; y < NBY-1; y++) imageE[t][z][y][x]=-Vx[y+1];
		}
	
		//2.2.3) diffusion - z implicit / x,y,t explicit
		//2.2.3.1) explicit part
		for (t = 1; t < NBT-1; t++) for (z = 1; z < NBZ-1; z++) for (y = 1; y < NBY-1; y++) for (x = 1; x < NBX-1; x++) {
			dIdx=(imageE[t][z][y][x+1]-imageE[t][z][y][x-1])/(2*dx);
			Dxx_div_dxSq=static_cast<float>((1-exp(-3.314/pow((dIdx/ax),4)))*DivPowDxSqu);
			dIdy=(imageE[t][z][y+1][x]-imageE[t][z][y-1][x])/(2*dy);
			Dyy_div_dySq=static_cast<float>((1-exp(-3.314/pow((dIdy/ay),4)))*DivPowDySqu);
			dIdt=(imageE[t+1][z][y][x]-imageE[t-1][z][y][x])/(2*dt);
			Dtt_div_dtSq=static_cast<float>((1-exp(-3.314/pow((dIdt/at),4)))*DivPowDzSqu);
			//new value of the voxel
			DivDgradI=(imageE[t][z][y][x+1]-2*imageE[t][z][y][x]+imageE[t][z][y][x-1])*Dxx_div_dxSq+
					(imageE[t][z][y+1][x]-2*imageE[t][z][y][x]+imageE[t][z][y-1][x])*Dyy_div_dySq+
					(imageE[t+1][z][y][x]-2*imageE[t][z][y][x]+imageE[t-1][z][y][x])*Dtt_div_dtSq;
	
			imageO[t][z][y][x]=imageE[t][z][y][x]+(dTau/4.)*DivDgradI;
		}
		
		//2.2.3.2) implicit part
		for (t = 1; t < NBT-1; t++) for (y = 1; y < NBY-1; y++) for (x = 1; x < NBX-1; x++) {
			for (z = 1; z < NBZ-1; z++){
				dIdz=(imageE[t][z+1][y][x]-imageE[t][z-1][y][x])/(2*dz);
				Dzz_div_dzSq=static_cast<float>((1-exp(-3.314/pow((dIdz/az),4)))*DivPowDzSqu);
				Va[z+1]=(dTau/4.)*Dzz_div_dzSq;
				Vb[z+1]=-1-2*(dTau/4.)*Dzz_div_dzSq;
				Vc[z+1]=(dTau/4.)*Dzz_div_dzSq;
				Vd[z+1]=imageE[t][z][y][x];
			}
			Va[1]=Va[3]; Va[0]=Va[4]; Va[NBZ]=Va[NBZ-2]; Va[NBZ+1]=Va[NBZ-3]; //to avoid boundary effects
			Vb[1]=Vb[3]; Vb[0]=Vb[4]; Vb[NBZ]=Vb[NBZ-2]; Vb[NBZ+1]=Vb[NBZ-3]; //to avoid boundary effects
			Vc[1]=Vc[3]; Vc[0]=Vc[4]; Vc[NBZ]=Vc[NBZ-2]; Vc[NBZ+1]=Vc[NBZ-3]; //to avoid boundary effects
			Vd[1]=Vd[3]; Vd[0]=Vd[4]; Vd[NBZ]=Vd[NBZ-2]; Vd[NBZ+1]=Vd[NBZ-3]; //to avoid boundary effects
			TridiagonalSolveFloat(Va,Vb,Vc,Vd,Vx,NBZ+2);
			for (z = 1; z < NBZ-1; z++) imageO[t][z][y][x]=-Vx[z+1];
		}
		
		//2.2.4) diffusion - t implicit / x,y,z explicit
		//2.2.4.1) explicit part
		for (t = 1; t < NBT-1; t++) for (z = 1; z < NBZ-1; z++) for (y = 1; y < NBY-1; y++) for (x = 1; x < NBX-1; x++) {
			dIdx=(imageO[t][z][y][x+1]-imageO[t][z][y][x-1])/(2*dx);
			Dxx_div_dxSq=static_cast<float>((1-exp(-3.314/pow((dIdx/ax),4)))*DivPowDxSqu);
			dIdy=(imageO[t][z][y+1][x]-imageO[t][z][y-1][x])/(2*dy);
			Dyy_div_dySq=static_cast<float>((1-exp(-3.314/pow((dIdy/ay),4)))*DivPowDySqu);
			dIdz=(imageO[t][z+1][y][x]-imageO[t][z-1][y][x])/(2*dz);
			Dzz_div_dzSq=static_cast<float>((1-exp(-3.314/pow((dIdz/az),4)))*DivPowDzSqu);
			//new value of the voxel
			DivDgradI=(imageO[t][z][y][x+1]-2*imageO[t][z][y][x]+imageO[t][z][y][x-1])*Dxx_div_dxSq+
					(imageO[t][z][y+1][x]-2*imageO[t][z][y][x]+imageO[t][z][y-1][x])*Dyy_div_dySq+
					(imageO[t][z+1][y][x]-2*imageO[t][z][y][x]+imageO[t][z-1][y][x])*Dzz_div_dzSq;
	
			imageE[t][z][y][x]=imageO[t][z][y][x]+(dTau/4.)*DivDgradI;
		}
		
		//2.2.4.2) implicit part
		for (z = 1; z < NBZ-1; z++) for (y = 1; y < NBY-1; y++) for (x = 1; x < NBX-1; x++) {
			for (t = 1; t < NBT-1; t++){
				dIdt=(imageO[t+1][z][y][x]-imageO[t-1][z][y][x])/(2*dt);
				Dtt_div_dtSq=static_cast<float>((1-exp(-3.314/pow((dIdt/at),4)))*DivPowDzSqu);
				Va[t+1]=(dTau/4.)*Dtt_div_dtSq;
				Vb[t+1]=-1-2*(dTau/4.)*Dtt_div_dtSq;
				Vc[t+1]=(dTau/4.)*Dtt_div_dtSq;
				Vd[t+1]=imageO[t][z][y][x];
			}
			Va[1]=Va[3]; Va[0]=Va[4]; Va[NBT]=Va[NBT-2]; Va[NBT+1]=Va[NBT-3]; //to avoid boundary effects
			Vb[1]=Vb[3]; Vb[0]=Vb[4]; Vb[NBT]=Vb[NBT-2]; Vb[NBT+1]=Vb[NBT-3]; //to avoid boundary effects
			Vc[1]=Vc[3]; Vc[0]=Vc[4]; Vc[NBT]=Vc[NBT-2]; Vc[NBT+1]=Vc[NBT-3]; //to avoid boundary effects
			Vd[1]=Vd[3]; Vd[0]=Vd[4]; Vd[NBT]=Vd[NBT-2]; Vd[NBT+1]=Vd[NBT-3]; //to avoid boundary effects
			TridiagonalSolveFloat(Va,Vb,Vc,Vd,Vx,NBT+2);
			for (t = 1; t < NBT-1; t++) imageE[t][z][y][x]=-Vx[t+1];
		}
	}
	
	//2.3) save the filtered temporary 3D image in VoxelType in the  output image at time t
	for (t = 0; t < NBT-2; t++) for (z = 0; z < NBZ-2; z++) for (y = 0; y < NBY-2; y++) for (x = 0; x < NBX-2; x++) 
		this->_output->Put(x, y, z, t, static_cast<VoxelType>(imageE[t+1][z+1][y+1][x+1]));
	
	
	//3) END OF THE FUNCTION
	// Do the final cleaning up
	this->Finalize();	
}


template <class VoxelType> void anisoDiffusion<VoxelType>::Run_4D_Explicit(){
	int i, j, x, y, z, t;
	double  ax,ay,az,at,dx,dy,dz,dt;
	float dTau;
	float**** imageE;
	float**** imageO;
	int NBX,NBY,NBZ,NBT;
	double dIdx,dIdy,dIdz,dIdt;
	float DivDgradI;
	int ITERATIONS_NB;
	float Dxx_div_dxSq,Dyy_div_dySq,Dzz_div_dzSq,Dtt_div_dtSq;
	int iteration;
	float DivPowDxSqu,DivPowDySqu,DivPowDzSqu,DivPowDtSqu;
	
	//1) INITIALISATION
	
	// Do the initial set up
	this->Initialize();
	
	//variables definition
	ax=this->ax;
	ay=this->ay;
	az=this->az;
	at=this->at;
	dx=this->dx;
	dy=this->dy;
	dz=this->dz;
	dt=this->dt;
	dTau=this->dTau;
	ITERATIONS_NB=this->ITERATIONS_NB;
	NBX=this->_input->GetX();
	NBY=this->_input->GetY();
	NBZ=this->_input->GetZ();
	NBT=this->_input->GetT();
	cout << "Image size: " << NBX <<  " , "  <<  NBY  <<  " , "  << NBZ  <<  " , " << NBT  << "\n";
	
	//precomputed values
	DivPowDxSqu=1./pow(dx,2);
	DivPowDySqu=1./pow(dy,2);
	DivPowDzSqu=1./pow(dz,2);
	DivPowDtSqu=1./pow(dt,2);
	
	//temporary input and output images and diffusion tensor field
	imageE= (float****) malloc (NBT*sizeof(float***));
	for (t=0;t<NBT;t++) imageE[t]= (float***) malloc (NBZ*sizeof(float**));
	for (t=0;t<NBT;t++) for (i=0;i<NBZ;i++) imageE[t][i]= (float**) malloc (NBY*sizeof(float*));
	for (t=0;t<NBT;t++) for (i=0;i<NBZ;i++) for (j=0;j<NBY;j++) imageE[t][i][j]= (float*) malloc (NBX*sizeof(float));
	
	imageO= (float****) malloc (NBT*sizeof(float***));
	for (t=0;t<NBT;t++) imageO[t]= (float***) malloc (NBZ*sizeof(float**));
	for (t=0;t<NBT;t++) for (i=0;i<NBZ;i++) imageO[t][i]= (float**) malloc (NBY*sizeof(float*));
	for (t=0;t<NBT;t++) for (i=0;i<NBZ;i++) for (j=0;j<NBY;j++) imageO[t][i][j]= (float*) malloc (NBX*sizeof(float));
	
	//2) ANISOTROPIC DIFFUSION
	
	//2.1) convert the values of the input image at time t in double in a temporary 3D image
	for (t = 0; t < NBT; t++) for (z = 0; z < NBZ; z++)  for (y = 0; y < NBY; y++) for (x = 0; x < NBX; x++)
				imageE[t][z][y][x]=static_cast<float>(this->_input->Get(x, y, z, t));
	
	for (t = 0; t < NBT; t++) for (z = 0; z < NBZ; z++) for (y = 0; y < NBY; y++){ 
		imageO[t][z][y][0]=imageE[t][z][y][0];
		imageO[t][z][y][NBX-1]=imageE[t][z][y][NBX-1];
	}
	
	for (t = 0; t < NBT; t++) for (z = 0; z < NBZ; z++) for (x = 0; x < NBX; x++){ 
		imageO[t][z][0][x]=imageE[t][z][0][x];
		imageO[t][z][NBY-1][x]=imageE[t][z][NBY-1][x];
	}
	
	for (t = 0; t < NBT; t++) for (y = 0; y < NBY; y++) for (x = 0; x < NBX; x++){ 
		imageO[t][0][y][x]=imageE[t][0][y][x];
		imageO[t][NBZ-1][y][x]=imageE[t][NBZ-1][y][x];
	}
	
	for (z = 0; z < NBZ; z++) for (y = 0; y < NBY; y++) for (x = 0; x < NBX; x++){ 
		imageO[0][z][y][x]=imageE[0][z][y][x];
		imageO[NBT-1][z][y][x]=imageE[NBT-1][z][y][x];
	}

	
	//2.2) diffusion in the temporary 3D image - explicit scheme
	for (iteration=0 ; iteration<ITERATIONS_NB; iteration++){
		cout << "| Iteration " << iteration+1 << " / " << ITERATIONS_NB << "\n";
		
		for (t = 1; t < NBT-1; t++) for (z = 1; z < NBZ-1; z++) for (y = 1; y < NBY-1; y++) for (x = 1; x < NBX-1; x++) {
			dIdx=(imageE[t][z][y][x+1]-imageE[t][z][y][x-1])/(2*dx);
			Dxx_div_dxSq=static_cast<float>((1-exp(-3.314/pow((dIdx/ax),4)))*DivPowDxSqu);
			dIdy=(imageE[t][z][y+1][x]-imageE[t][z][y-1][x])/(2*dy);
			Dyy_div_dySq=static_cast<float>((1-exp(-3.314/pow((dIdy/ay),4)))*DivPowDySqu);
			dIdz=(imageE[t][z+1][y][x]-imageE[t][z-1][y][x])/(2*dz);
			Dzz_div_dzSq=static_cast<float>((1-exp(-3.314/pow((dIdz/az),4)))*DivPowDzSqu);
			dIdt=(imageE[t+1][z][y][x]-imageE[t-1][z][y][x])/(2*dt);
			Dtt_div_dtSq=static_cast<float>((1-exp(-3.314/pow((dIdt/at),4)))*DivPowDzSqu);
			
			//new value of the voxel
			DivDgradI=(imageE[t][z][y][x+1]-2*imageE[t][z][y][x]+imageE[t][z][y][x-1])*Dxx_div_dxSq+
					(imageE[t][z][y+1][x]-2*imageE[t][z][y][x]+imageE[t][z][y-1][x])*Dyy_div_dySq+
					(imageE[t][z+1][y][x]-2*imageE[t][z][y][x]+imageE[t][z-1][y][x])*Dzz_div_dzSq+
					(imageE[t+1][z][y][x]-2*imageE[t][z][y][x]+imageE[t-1][z][y][x])*Dzz_div_dzSq;
			
			imageO[t][z][y][x]=imageE[t][z][y][x]+(dTau)*DivDgradI;
		}
		
	
		//2.2.5) temporary output image is reinjected in temporary input image
		for (t = 0; t < NBT; t++) for (z = 0; z < NBZ; z++) for (y = 0; y < NBY; y++) for (x = 0; x < NBX; x++)
			imageE[t][z][y][x]=imageO[t][z][y][x];

	}
	
	//2.3) save the filtered temporary 3D image in VoxelType in the  output image at time t
	for (t = 0; t < NBT; t++) for (z = 0; z < NBZ; z++) for (y = 0; y < NBY; y++) for (x = 0; x < NBX; x++) 
		this->_output->Put(x, y, z, t, static_cast<VoxelType>(imageE[t][z][y][x]));

	
	//3) END OF THE FUNCTION
	// Do the final cleaning up
	this->Finalize();
}

///+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
///+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
///+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
/*
template <class VoxelType> void anisoDiffusion<VoxelType>::Run_3D_Explicit(){
	int i, j, x, y, z, t;
	double  ax,ay,az,at,dx,dy,dz,dt;
	float dTau;
	float*** imageE;
	float*** imageO;
	int NBX,NBY,NBZ,NBT;
	double dIdx,dIdy,dIdz;
	float DivDgradI;
	int ITERATIONS_NB;
	float Dxx_div_dxSq,Dyy_div_dySq,Dzz_div_dzSq;
	int iteration;
	float DivPowDxSqu,DivPowDySqu,DivPowDzSqu,DivPowDtSqu;
	
	//1) INITIALISATION
	
	// Do the initial set up
	this->Initialize();
	
	//variables definition
	ax=this->ax;
	ay=this->ay;
	az=this->az;
	at=this->at;
	dx=this->dx;
	dy=this->dy;
	dz=this->dz;
	dt=this->dt;
	dTau=this->dTau;
	ITERATIONS_NB=this->ITERATIONS_NB;
	NBX=this->_input->GetX();
	NBY=this->_input->GetY();
	NBZ=this->_input->GetZ();
	NBT=this->_input->GetT();
	cout << "Image size: " << NBX <<  " , "  <<  NBY  <<  " , "  << NBZ  <<  " , " << NBT  << "\n";
        cout << "TOTO ";
	//precomputed values
	DivPowDxSqu=1./pow(dx,2);
	DivPowDySqu=1./pow(dy,2);
	DivPowDzSqu=1./pow(dz,2);
	DivPowDtSqu=1./pow(dt,2);
	
	//temporary input and output images
	imageE= (float***) malloc (NBZ*sizeof(float**));
	for (i=0;i<NBZ;i++) imageE[i]= (float**) malloc (NBY*sizeof(float*));
	for (i=0;i<NBZ;i++) for (j=0;j<NBY;j++) imageE[i][j]= (float*) malloc (NBX*sizeof(float));
	
	imageO= (float***) malloc (NBZ*sizeof(float**));
	for (i=0;i<NBZ;i++) imageO[i]= (float**) malloc (NBY*sizeof(float*));
	for (i=0;i<NBZ;i++) for (j=0;j<NBY;j++) imageO[i][j]= (float*) malloc (NBX*sizeof(float));
	

	//2) ANISOTROPIC DIFFUSION
	for (t = 0; t < NBT; t++) {
		cout << "Image " << t+1 << " / " << NBT << "\n";
		
		//2.1) convert the values of the input image at time t in double in a temporary 3D image
		for (z = 0; z < NBZ; z++)  for (y = 0; y < NBY; y++) for (x = 0; x < NBX; x++)
			imageE[z][y][x]=static_cast<float>(this->_input->Get(x, y, z, t));
		
		for (z = 0; z < NBZ; z++)  for (y = 0; y < NBY; y++) for (x = 0; x < NBX; x++)
			imageO[z][y][x]=static_cast<float>(this->_input->Get(x, y, z, t));
		
		
		//2.2) diffusion in the temporary 3D image - explicit scheme
		for (iteration=0 ; iteration<ITERATIONS_NB; iteration++){
			cout << "| Iteration " << iteration+1 << " / " << ITERATIONS_NB << "\n";
			
			for (z = 1; z < NBZ-1; z++) for (y = 1; y < NBY-1; y++) for (x = 1; x < NBX-1; x++) {
				dIdx=(imageE[z][y][x+1]-imageE[z][y][x-1])/(2*dx);
				Dxx_div_dxSq=static_cast<float>((1-exp(-3.314/pow((dIdx/ax),4)))*DivPowDxSqu);
				dIdy=(imageE[z][y+1][x]-imageE[z][y-1][x])/(2*dy);
				Dyy_div_dySq=static_cast<float>((1-exp(-3.314/pow((dIdy/ay),4)))*DivPowDySqu);
				dIdz=(imageE[z+1][y][x]-imageE[z-1][y][x])/(2*dz);
				Dzz_div_dzSq=static_cast<float>((1-exp(-3.314/pow((dIdz/az),4)))*DivPowDzSqu);
				
				//new value of the voxel
				DivDgradI=(imageE[z][y][x+1]-2*imageE[z][y][x]+imageE[z][y][x-1])*Dxx_div_dxSq+
					(imageE[z][y+1][x]-2*imageE[z][y][x]+imageE[z][y-1][x])*Dyy_div_dySq+
					(imageE[z+1][y][x]-2*imageE[z][y][x]+imageE[z-1][y][x])*Dzz_div_dzSq;
	
				imageO[z][y][x]=imageE[z][y][x]+(dTau)*DivDgradI;
			}
				
		
			//2.2.5) temporary output image is reinjected in temporary input image
			for (z = 0; z < NBZ; z++) for (y = 0; y < NBY; y++) for (x = 0; x < NBX; x++)
						imageE[z][y][x]=imageO[z][y][x];
	
		}
		
		//2.3) save the filtered temporary 3D image in VoxelType in the  output image at time t
		for (z = 0; z < NBZ; z++) for (y = 0; y < NBY; y++) for (x = 0; x < NBX; x++) 
			this->_output->Put(x, y, z, t, static_cast<VoxelType>(imageE[z][y][x]));
	}
	
	//3) END OF THE FUNCTION
	// Do the final cleaning up
	this->Finalize();
}
*/







///+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
///+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
///+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
///                                          BEGIN TENSOR VOTING PROJECT
///+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
///+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
///+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

//Stucture utilisee pour stoquer une image 3D de float.
//remarque : en utilisant des float, l'image est 4 fois plus grosse que l'image initiale en unsigned char
typedef struct
{
  float ***image;  // l'image elle meme
  int NBZ;                 // |
  int NBY;                 // |-> dimensions de l'image
  int NBX;                 // |
} Image3Dfloat;

// Field[0][0] = sum v_x^2   / Field[0][1] = sum v_x v_y / Field[0][2] = sum v_x v_z
// Field[1][0] = sum v_y v_x / Field[1][1] = sum v_y^2   / Field[1][2] = sum v_y v_z
// Field[2][0] = sum v_z v_x / Field[2][1] = sum v_z v_y / Field[2][2] = sum v_z^2
typedef struct
{
  Image3Dfloat Field[3][3];
  int NBZ;                 // |
  int NBY;                 // |-> dimensions des images
  int NBX;                 // |
} TensorField;




//cree une image 3d codee en float avec toutes ses valeurs a zero.
extern void CreateImage3DFloat(Image3Dfloat * img3d,int NZ,int NY,int NX){
  int i, j, k;

  img3d->NBZ=NZ;
  img3d->NBY=NY;
  img3d->NBX=NX;

  img3d->image = (float***)malloc((img3d->NBZ)*sizeof(float**));
  for (i=0;i<img3d->NBZ;i++)
    img3d->image[i]=(float**)malloc((img3d->NBY)*sizeof(float*));
  for (i=0;i<img3d->NBZ;i++)
    for (j=0;j<img3d->NBY;j++)
      img3d->image[i][j]=(float*)malloc((img3d->NBX)*sizeof(float));

  for (i=0;i<img3d->NBZ;i++)
    for (j=0;j<img3d->NBY;j++)
      for (k=0;k<img3d->NBX;k++)
        img3d->image[i][j][k]=0;
}



// ---------------------------------------------------------------------------------------
//   PART. 5.4                      Le Tensor Voting
// ---------------------------------------------------------------------------------------

//adapte de l'algorithme du meme nom dans numerical recipes.
//en entree, on a la matrice 'MatIni' de dimension n*n. Elle doit etre symetrique.
//en sortie, ValP est un vecteur de taille n qui contient les valeurs propres (dans l'ordre decroissant). 
//VecP est une matrice n*n qui contient les vecteurs propres en colonne.
//remarque : la version avec n comme parametre d'entree cree beaucoup de fuites de memoire sur
//ma machine (pourquoi ???). Du coup, cette version fonctionne pour n fixe a 3.
#define ROTATE(a,i,j,k,l) g=a[i][j];h=a[k][l];a[i][j]=g-s*(h+g*tau);a[k][l]=h+s*(g-h*tau);
void jacobi3(float **MatIni,float *ValP, float **VecP){
  int j,iq,ip,i; 
  float tresh,theta,tau,t,sm,s,h,g,c;
  float b[4];
  float z[4];
  float a[4][4];   //correspond a MatIni
  float d[4];    //correspond a ValP
  float v[4][4];   //correspond a VecP
  int vTri1,vTri2;
  float TempF;
  int n;


  n=3;
  for(i=0;i<n;i++) for(j=0;j<n;j++) a[i+1][j+1]=MatIni[i][j];


//algo de numerical recipes
  for (ip=1;ip<=n;ip++) {
    for (iq=1;iq<=n;iq++) v[ip][iq]=0.0;
    v[ip][ip]=1.0;
  }
	
  for (ip=1;ip<=n;ip++) {
    b[ip]=d[ip]=a[ip][ip];
    z[ip]=0.0;
  }
	
  for (i=1;i<=50;i++) {
    sm=0.0;
    for (ip=1;ip<=n-1;ip++)
      for (iq=ip+1;iq<=n;iq++)
        sm += fabs(a[ip][iq]);
	
    if (sm == 0.0) {
		//adaptation des valeurs de l'algo de numerical recipes aux valeurs de sortie
      for(i=0;i<n;i++) ValP[i]=d[i+1];
      for(i=0;i<n;i++) for(j=0;j<n;j++) MatIni[i][j]=a[i+1][j+1];
      for(i=0;i<n;i++) for(j=0;j<n;j++) VecP[i][j]=v[i+1][j+1];
		
		//tri des donnees
      for(vTri1=0;vTri1<n-1;vTri1++) for(vTri2=vTri1+1;vTri2<n;vTri2++) if (ValP[vTri1]<ValP[vTri2]){
        TempF=ValP[vTri1]; ValP[vTri1]=ValP[vTri2]; ValP[vTri2]=TempF;
        for(i=0;i<n;i++) { TempF=VecP[i][vTri1]; VecP[i][vTri1]=VecP[i][vTri2]; VecP[i][vTri2]=TempF;}
      }
		
      return;
    }
    if (i < 4) tresh=0.2*sm/(n*n);
    else tresh=0.0;
	
    for (ip=1;ip<=n-1;ip++) {
      for (iq=ip+1;iq<=n;iq++) {
        g=100.0*fabs(a[ip][iq]);
			
        if (i > 4 && (float)(fabs(d[ip])+g) == (float)fabs(d[ip])&& (float)(fabs(d[iq])+g) == (float)fabs(d[iq]))
          a[ip][iq]=0.0;
        else if (fabs(a[ip][iq]) > tresh) {
          h=d[iq]-d[ip];
          if ((float)(fabs(h)+g) == (float)fabs(h))
            t=(a[ip][iq])/h;
          else {
            theta=0.5*h/(a[ip][iq]);
            t=1.0/(fabs(theta)+sqrt(1.0+theta*theta));
            if (theta < 0.0) t = -t;
          }
          c=1.0/sqrt(1+t*t); s=t*c; tau=s/(1.0+c); h=t*a[ip][iq];
          z[ip] -= h; z[iq] += h; d[ip] -= h; d[iq] += h; a[ip][iq]=0.0;
          for (j=1;j<=ip-1;j++) { ROTATE(a,j,ip,j,iq) }
          for (j=ip+1;j<=iq-1;j++) { ROTATE(a,ip,j,j,iq) }
          for (j=iq+1;j<=n;j++) { ROTATE(a,ip,j,iq,j) }
          for (j=1;j<=n;j++) { ROTATE(v,j,ip,j,iq)}
        }
      }
    }
    for (ip=1;ip<=n;ip++) { b[ip] += z[ip]; d[ip]=b[ip]; z[ip]=0.0; }
  }
  printf("Too many iterations in the routine jacobi\n");

//adaptation des valeurs de l'algo de numerical recipes aux valeurs de sortie
  for(i=0;i<n;i++) ValP[i]=d[i+1];
  for(i=0;i<n;i++) for(j=0;j<n;j++) MatIni[i][j]=a[i+1][j+1];
  for(i=0;i<n;i++) for(j=0;j<n;j++) VecP[i][j]=v[i+1][j+1];

//tri des donnees
  for(vTri1=0;vTri1<n-1;vTri1++) for(vTri2=vTri1+1;vTri2<n;vTri2++) if (ValP[vTri1]<ValP[vTri2]){
    TempF=ValP[vTri1]; ValP[vTri1]=ValP[vTri2]; ValP[vTri2]=TempF;
    for(i=0;i<n;i++) { TempF=VecP[i][vTri1]; VecP[i][vTri1]=VecP[i][vTri2]; VecP[i][vTri2]=TempF;}
  }

}



//initialisation a zero d'un champ de tenseur
void InitTensorField(TensorField * TF, int NBX,int NBY,int NBZ){
  int i,j,k;
  
  //initialisation des tailles
  TF->NBZ=NBZ;
  TF->NBY=NBY;
  TF->NBX=NBX;
  
  //initialisation du champ de tenseurs
  for (i=0;i<3;i++) for (j=0;j<3;j++){
    CreateImage3DFloat(&(TF->Field[i][j]),NBZ,NBY,NBX);
  }
  
  //set everything to 0 (it should be already done, it's just to be sure...)
  for (i=0;i<TF->NBZ;i++) for (j=0;j<TF->NBY;j++) for (k=0;k<TF->NBX;k++){ 
    //nettoyage du champ de tenseur
    TF->Field[0][0].image[i][j][k]=0;
    TF->Field[1][0].image[i][j][k]=0;
    TF->Field[2][0].image[i][j][k]=0;
    TF->Field[0][1].image[i][j][k]=0;
    TF->Field[1][1].image[i][j][k]=0;
    TF->Field[2][1].image[i][j][k]=0;
    TF->Field[0][2].image[i][j][k]=0;
    TF->Field[1][2].image[i][j][k]=0;
    TF->Field[2][2].image[i][j][k]=0;
    
  }
  
}



//insertion des ball voting field dans le champ de tenseurs ; 
void InsertBallFields(TensorField * TF,irtkGenericImage<unsigned char> * InputImage, double sigma,int Tboites,irtkGenericImage<float> * Saliency){
  int i,j,k;
  double V_x,V_y,V_z;
  double Dist;
  int LocX,LocY,LocZ;
  double Poids;
  int CX,CY,CZ;
  
  // Tboites devient la moitie d'un cote de boite (pour coller aux boucles_for)
  Tboites=(Tboites-1)/2;
  
  // 1 ) MISE A JOUR DU CHAMP DE TENSEUR
  for(CZ=0;CZ<InputImage->GetZ();CZ++)  for(CY=0;CY<InputImage->GetY();CY++) for(CX=0;CX<InputImage->GetX();CX++) if (InputImage->Get(CX, CY, CZ, 0)>0){
  
    //cout << CX << " " << CY << " " << CZ << "\n";
    Saliency->Put(CX, CY, CZ,0, 1.);
    
    //remplissage du champ de tenseurs
    for (i=-Tboites;i<Tboites;i++) for (j=-Tboites;j<Tboites;j++) for (k=-Tboites;k<Tboites;k++){
      LocX=CX+k;
      LocY=CY+j;
      LocZ=CZ+i;
      if ((LocX>0)&&(LocX<TF->NBX)&&(LocY>0)&&(LocY<TF->NBY)&&(LocZ>0)&&(LocZ<TF->NBZ)){
        //A ) vecteur a injecter dans le tenseur norm\'e
        //vecteur norm\'e que l'on va injecter dans le tenseur (apres ponderation)
        Dist=sqrt(pow((double)k,2.0)+pow((double)j,2.0)+pow((double)i,2.0));
        V_x=((double)k)/Dist;
        V_y=((double)j)/Dist;
        V_z=((double)i)/Dist;
  
        //B ) ponderation du vecteur
        Poids=exp(-pow(Dist,2.0)/pow(sigma,2.0));
        
        //cout << Poids << "\n";
        V_x=V_x*Poids;
        V_y=V_y*Poids;
        V_z=V_z*Poids;
        
        //C ) injection du vecteur dans le tenseur
        TF->Field[0][0].image[LocZ][LocY][LocX]+=(float)(V_x*V_x); 
        TF->Field[0][1].image[LocZ][LocY][LocX]+=(float)(V_x*V_y); 
        TF->Field[0][2].image[LocZ][LocY][LocX]+=(float)(V_x*V_z);
        TF->Field[1][0].image[LocZ][LocY][LocX]+=(float)(V_x*V_y); 
        TF->Field[1][1].image[LocZ][LocY][LocX]+=(float)(V_y*V_y); 
        TF->Field[1][2].image[LocZ][LocY][LocX]+=(float)(V_y*V_z);
        TF->Field[2][0].image[LocZ][LocY][LocX]+=(float)(V_x*V_z); 
        TF->Field[2][1].image[LocZ][LocY][LocX]+=(float)(V_z*V_y); 
        TF->Field[2][2].image[LocZ][LocY][LocX]+=(float)(V_z*V_z);
      }
    }
  }
}



//a partir d'un champ de tenseur rempli, calcul des valeurs propres
void CalcFunctionnal(TensorField * TF,irtkGenericImage<float> * lambda1,irtkGenericImage<float> * lambda2,irtkGenericImage<float> * lambda3,irtkGenericImage<float> * Saliency){
  int i,j,k;
  float ** a;
  float ** q;
  float *  d;
  
  //allocation memoire pour les variables utilisees dans l'appel de la fonction de Jacobi
  a=(float**)malloc(3*sizeof(float*));
  for(i=0;i<3;i++) a[i]=(float*)malloc(3*sizeof(float));
  q=(float**)malloc(3*sizeof(float*));
  for(i=0;i<3;i++) q[i]=(float*)malloc(3*sizeof(float));
  d=(float*)malloc(3*sizeof(float));
  
  for (i=0;i<TF->NBZ;i++) for (j=0;j<TF->NBY;j++) for (k=0;k<TF->NBX;k++){
    //cout << k << " " <<  j << " " <<  i << "\n";
    if ((TF->Field[0][0].image[i][j][k]>0.0001)||(TF->Field[1][1].image[i][j][k]>0.0001)||(TF->Field[2][2].image[i][j][k]>0.0001)){
      //remplissage de la matrice dont on extrait les valeurs propres
      a[0][0]=TF->Field[0][0].image[i][j][k];
      a[1][0]=TF->Field[1][0].image[i][j][k];
      a[2][0]=TF->Field[2][0].image[i][j][k];
      a[0][1]=TF->Field[0][1].image[i][j][k];
      a[1][1]=TF->Field[1][1].image[i][j][k];
      a[2][1]=TF->Field[2][1].image[i][j][k];
      a[0][2]=TF->Field[0][2].image[i][j][k];
      a[1][2]=TF->Field[1][2].image[i][j][k];
      a[2][2]=TF->Field[2][2].image[i][j][k];
      
      //extraction des valeurs propres
      jacobi3(a,d,q);
      
/*      d[0]=log(d[0]); d[1]=log(d[1]); d[2]=log(d[2]);
      if (d[0]>50) d[0]=50;  if ((d[0]<-50)||(isnan(d[0]))) d[0]=-50;
      if (d[1]>50) d[1]=50;  if ((d[1]<-50)||(isnan(d[1]))) d[1]=-50;
      if (d[2]>50) d[2]=50;  if ((d[2]<-50)||(isnan(d[2]))) d[2]=-50;*/
      
      
      //cout << d[0] << " " <<  d[1] << " " <<  d[2] << "\n";
      //remplissage des valeurs propres
      lambda1->Put(k, j, i,0, (float)(d[0]));
      lambda2->Put(k, j, i,0, (float)(d[1]));
      lambda3->Put(k, j, i,0, (float)(d[2]));
      if ((d[1]>0.00001)&&(Saliency->Get(k, j, i,0)<0.5)) if ((d[0]/d[1]>5.)) Saliency->Put(k, j, i,0, 0.5);
    }
  }
}




//Dans le cadre de l'algorithme de tensor voting :
// -> On souhaite perdre 1/e d'energie a la distance 'dista' de l'origine (point O) dans la direction du bout de segment (point A)
// -> La taille de la fenetre qui contient le champ de tenseur doit de meme contenir toute l'info pour laquelle l'energie est > 0.01
//Cette fonction calcule alors 'c', 'sigma' et 'Tfenetre' en fonction de 'dista' et 'angl'.
void ComputeSigmaAndBoxSize(double dista,double * sigma,int * Tfenetre){
  //correction des entrees
  dista=fabs(dista);
  
  //calcul des coefficients
  *sigma=dista;
  *Tfenetre=5*dista+1;
  
}



/*
template <class VoxelType> void anisoDiffusion<VoxelType>::Run_3D_Explicit(){
  int NBX,NBY,NBZ;
  char lambda1_output_name[] = "TV_lambda1.nii";
  char lambda2_output_name[] = "TV_lambda2.nii";
  char lambda3_output_name[] = "TV_lambda3.nii";
  double CharactDist = this->dTau;
  irtkGenericImage<unsigned char> SegImage;
  irtkGenericImage<float> lambda1;
  irtkGenericImage<float> lambda2;
  irtkGenericImage<float> lambda3;
  irtkGenericImage<float> Saliency;
  TensorField TF;
  double sigma;
  int Tboites;
  int x,y,z;
  
  //1) initialisation
  
  cout << "Compute the tensor field...\n";
  
  this->Initialize();
  NBX=this->_input->GetX();
  NBY=this->_input->GetY();
  NBZ=this->_input->GetZ();
  
  SegImage = irtkGenericImage<unsigned char>(NBX,NBY,NBZ,1);
  lambda1 = irtkGenericImage<float>(NBX,NBY,NBZ,1);
  lambda2 = irtkGenericImage<float>(NBX,NBY,NBZ,1);
  lambda3 = irtkGenericImage<float>(NBX,NBY,NBZ,1);
  Saliency = irtkGenericImage<float>(NBX,NBY,NBZ,1);
  
  //cast the input image
  for (z=0;z<NBZ;z++) for (y=0;y<NBY;y++) for (x=0;x<NBX;x++){
    SegImage.Put(x,y,z,0,(unsigned char)(this->_input->Get(x, y, z, 0)+0.00001));
  }
  
  //compute 'sigma' and the size of the boxes 'Tboites' as a function of 'CharactDist'
  ComputeSigmaAndBoxSize(CharactDist,&sigma,&Tboites);
  
  
  cout << "Sigma=" << sigma << " / Box size=" << Tboites << "\n";
  
  //Tensor field initialisation
  InitTensorField(&TF,NBX,NBY,NBZ);
  
  // 2 ) Compute the tensor field and extract the eigenvalues at each point (voxel) of the field
  InsertBallFields(&TF,&SegImage,sigma,Tboites,&Saliency);
  
  CalcFunctionnal(&TF,&lambda1,&lambda2,&lambda3,&Saliency);
  
  //N) write the 3 lambda images
  lambda1.Write(lambda1_output_name);
  lambda2.Write(lambda2_output_name);
  lambda3.Write(lambda3_output_name);
  Saliency.Write("SaliencyMap.nii");
  
}
*/

///+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
///+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
///+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
///                                           END TENSOR VOTING PROJECT
///+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
///+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
///+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++



///+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
///+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
///+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
///                                       BEGIN GRADIENT VECTOR FLOW PROJECT
///+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
///+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
///+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


/*

template <class VoxelType> void anisoDiffusion<VoxelType>::Run_3D_Explicit(){
  int NBX,NBY,NBZ,x,y,z;
  char Grad_f_x_file[] = "gradientVecX.nii";
  char Grad_f_y_file[] = "gradientVecY.nii";
  char Grad_f_z_file[] = "gradientVecZ.nii";
  irtkGenericImage<float> Grad_f_x;     //input gradient of the image I on x
  irtkGenericImage<float> Grad_f_y;     //input gradient of the image I on y
  irtkGenericImage<float> Grad_f_z;     //input gradient of the image I on z
  irtkGenericImage<float> SqNormGrad_f; //square norm of the input gradient of the image I
  irtkGenericImage<float> u_cur;        //output regularization of Grad_f_x
  irtkGenericImage<float> v_cur;        //output regularization of Grad_f_y
  irtkGenericImage<float> w_cur;        //output regularization of Grad_f_z
  irtkGenericImage<float> u_next;       //regularization of Grad_f_x at the iteration after the current one
  irtkGenericImage<float> v_next;       //regularization of Grad_f_y at the iteration after the current one
  irtkGenericImage<float> w_next;       //regularization of Grad_f_z at the iteration after the current one
  double tmpdbl1,tmpdbl2,tmpdbl3;
  float tmpf1,tmpf2,tmpf3,tmpf4,tmpf5;
  float mu;
  int IterationNb,it;
  float DeltaXsq,DeltaYsq,DeltaZsq,DeltaT;
  
  
  
  //init
  this->Initialize();
  mu=1;
  IterationNb=this->ITERATIONS_NB;
  DeltaXsq=1.;
  DeltaYsq=1.;
  DeltaZsq=1.;
  DeltaT=this->at;
  
  
  //CFL respected?
  if ((double)DeltaT>sqrt((double)DeltaXsq)*sqrt((double)DeltaYsq)*sqrt((double)DeltaZsq)/(4*mu)){
    cout << DeltaT << " " << sqrt((double)DeltaXsq)*sqrt((double)DeltaYsq)*sqrt((double)DeltaZsq)/(4*mu) << "\n";
    DeltaT=(float)(sqrt((double)DeltaXsq)*sqrt((double)DeltaYsq)*sqrt((double)DeltaZsq)/(4*mu));
  }
  
  
  
  
  //read the gradients of f
  Grad_f_x.Read(Grad_f_x_file);
  Grad_f_y.Read(Grad_f_y_file);
  Grad_f_z.Read(Grad_f_z_file);
  
  //size of the images
  NBX=Grad_f_x.GetX();
  NBY=Grad_f_x.GetY();
  NBZ=Grad_f_x.GetZ();
  
  //regularization of grad f  (current iteration and next one)
  u_cur = irtkGenericImage<float>(NBX, NBY, NBZ, 1);
  v_cur = irtkGenericImage<float>(NBX, NBY, NBZ, 1);
  w_cur = irtkGenericImage<float>(NBX, NBY, NBZ, 1);
  u_next = irtkGenericImage<float>(NBX, NBY, NBZ, 1);
  v_next = irtkGenericImage<float>(NBX, NBY, NBZ, 1);
  w_next = irtkGenericImage<float>(NBX, NBY, NBZ, 1);
  
  //compute once for all the square norm of Grad_f
  SqNormGrad_f = irtkGenericImage<float>(NBX, NBY, NBZ, 1);
  for (z=0;z<NBZ;z++) for (y=0;y<NBY;y++) for (x=0;x<NBX;x++){
    tmpdbl1=pow((double)Grad_f_x.Get(x, y, z, 0),2.0);
    tmpdbl2=pow((double)Grad_f_x.Get(x, y, z, 0),2.0);
    tmpdbl3=pow((double)Grad_f_x.Get(x, y, z, 0),2.0);
    
    SqNormGrad_f.Put(x, y, z, 0, (float)(tmpdbl1+tmpdbl2+tmpdbl3));
  }
  
  //initialisation of u, v, w
  for (z=0;z<NBZ;z++) for (y=0;y<NBY;y++) for (x=0;x<NBX;x++) u_cur.Put(x, y, z, 0, (float)Grad_f_x.Get(x, y, z, 0));
  for (z=0;z<NBZ;z++) for (y=0;y<NBY;y++) for (x=0;x<NBX;x++) v_cur.Put(x, y, z, 0, (float)Grad_f_y.Get(x, y, z, 0));
  for (z=0;z<NBZ;z++) for (y=0;y<NBY;y++) for (x=0;x<NBX;x++) w_cur.Put(x, y, z, 0, (float)Grad_f_z.Get(x, y, z, 0));
  for (z=0;z<NBZ;z++) for (y=0;y<NBY;y++) for (x=0;x<NBX;x++) u_next.Put(x, y, z, 0, (float)Grad_f_x.Get(x, y, z, 0));
  for (z=0;z<NBZ;z++) for (y=0;y<NBY;y++) for (x=0;x<NBX;x++) v_next.Put(x, y, z, 0, (float)Grad_f_y.Get(x, y, z, 0));
  for (z=0;z<NBZ;z++) for (y=0;y<NBY;y++) for (x=0;x<NBX;x++) w_next.Put(x, y, z, 0, (float)Grad_f_z.Get(x, y, z, 0));
  
  
  //resolution
  for (it=0;it<IterationNb;it++){
    cout << "Iteration " << it << "\n";
    //compute the vector field of next iteration...
    //a) upate on u
    for (z=1;z<NBZ-1;z++) for (y=1;y<NBY-1;y++) for (x=1;x<NBX-1;x++){
      tmpf1=(mu/DeltaXsq)*(u_cur.Get(x+1, y, z, 0)-2*u_cur.Get(x, y, z, 0)+u_cur.Get(x-1, y, z, 0));
      tmpf2=(mu/DeltaYsq)*(u_cur.Get(x, y+1, z, 0)-2*u_cur.Get(x, y, z, 0)+u_cur.Get(x, y-1, z, 0));
      tmpf3=(mu/DeltaZsq)*(u_cur.Get(x, y, z+1, 0)-2*u_cur.Get(x, y, z, 0)+u_cur.Get(x, y, z-1, 0));
      tmpf4=u_cur.Get(x, y, z, 0)-Grad_f_x.Get(x, y, z, 0);
      tmpf5=SqNormGrad_f.Get(x, y, z, 0);
      u_next.Put(x, y, z, 0,u_cur.Get(x, y, z, 0)+DeltaT*(tmpf1+tmpf2+tmpf3-tmpf4*tmpf5));
    }
    
    //b) upate on v
    for (z=1;z<NBZ-1;z++) for (y=1;y<NBY-1;y++) for (x=1;x<NBX-1;x++){
      tmpf1=(mu/DeltaXsq)*(v_cur.Get(x+1, y, z, 0)-2*v_cur.Get(x, y, z, 0)+v_cur.Get(x-1, y, z, 0));
      tmpf2=(mu/DeltaYsq)*(v_cur.Get(x, y+1, z, 0)-2*v_cur.Get(x, y, z, 0)+v_cur.Get(x, y-1, z, 0));
      tmpf3=(mu/DeltaZsq)*(v_cur.Get(x, y, z+1, 0)-2*v_cur.Get(x, y, z, 0)+v_cur.Get(x, y, z-1, 0));
      tmpf4=v_cur.Get(x, y, z, 0)-Grad_f_y.Get(x, y, z, 0);
      tmpf5=SqNormGrad_f.Get(x, y, z, 0);
      v_next.Put(x, y, z, 0,v_cur.Get(x, y, z, 0)+DeltaT*(tmpf1+tmpf2+tmpf3-tmpf4*tmpf5));
    }
    
    //c) upate on w
    for (z=1;z<NBZ-1;z++) for (y=1;y<NBY-1;y++) for (x=1;x<NBX-1;x++){
      tmpf1=(mu/DeltaXsq)*(w_cur.Get(x+1, y, z, 0)-2*w_cur.Get(x, y, z, 0)+w_cur.Get(x-1, y, z, 0));
      tmpf2=(mu/DeltaYsq)*(w_cur.Get(x, y+1, z, 0)-2*w_cur.Get(x, y, z, 0)+w_cur.Get(x, y-1, z, 0));
      tmpf3=(mu/DeltaZsq)*(w_cur.Get(x, y, z+1, 0)-2*w_cur.Get(x, y, z, 0)+w_cur.Get(x, y, z-1, 0));
      tmpf4=w_cur.Get(x, y, z, 0)-Grad_f_z.Get(x, y, z, 0);
      tmpf5=SqNormGrad_f.Get(x, y, z, 0);
      w_next.Put(x, y, z, 0,w_cur.Get(x, y, z, 0)+DeltaT*(tmpf1+tmpf2+tmpf3-tmpf4*tmpf5));
    }
    
    //boundary conditions
    for (z=1;z<NBZ-1;z++) for (y=1;y<NBY-1;y++) u_next.Put(0, y, z, 0,u_next.Get(1, y, z, 0));
    for (z=1;z<NBZ-1;z++) for (y=1;y<NBY-1;y++) u_next.Put(NBX-1, y, z, 0,u_next.Get(NBX-2, y, z, 0));
    for (z=1;z<NBZ-1;z++) for (x=0;x<NBX;x++)   u_next.Put(x, 0, z, 0,u_next.Get(x, 1, z, 0));
    for (z=1;z<NBZ-1;z++) for (x=0;x<NBX;x++)   u_next.Put(x, NBY-1, z, 0,u_next.Get(x, NBY-2, z, 0));
    for (y=0;y<NBY;y++)   for (x=0;x<NBX;x++)   u_next.Put(x, y, 0, 0,u_next.Get(x, y, 1, 0));
    for (y=0;y<NBY;y++)   for (x=0;x<NBX;x++)   u_next.Put(x, y, NBZ-1, 0,u_next.Get(x, y, NBZ-2, 0));
    
    for (z=1;z<NBZ-1;z++) for (y=1;y<NBY-1;y++) v_next.Put(0, y, z, 0,v_next.Get(1, y, z, 0));
    for (z=1;z<NBZ-1;z++) for (y=1;y<NBY-1;y++) v_next.Put(NBX-1, y, z, 0,v_next.Get(NBX-2, y, z, 0));
    for (z=1;z<NBZ-1;z++) for (x=0;x<NBX;x++)   v_next.Put(x, 0, z, 0,v_next.Get(x, 1, z, 0));
    for (z=1;z<NBZ-1;z++) for (x=0;x<NBX;x++)   v_next.Put(x, NBY-1, z, 0,v_next.Get(x, NBY-2, z, 0));
    for (y=0;y<NBY;y++)   for (x=0;x<NBX;x++)   v_next.Put(x, y, 0, 0,v_next.Get(x, y, 1, 0));
    for (y=0;y<NBY;y++)   for (x=0;x<NBX;x++)   v_next.Put(x, y, NBZ-1, 0,v_next.Get(x, y, NBZ-2, 0));

    for (z=1;z<NBZ-1;z++) for (y=1;y<NBY-1;y++) w_next.Put(0, y, z, 0,w_next.Get(1, y, z, 0));
    for (z=1;z<NBZ-1;z++) for (y=1;y<NBY-1;y++) w_next.Put(NBX-1, y, z, 0,w_next.Get(NBX-2, y, z, 0));
    for (z=1;z<NBZ-1;z++) for (x=0;x<NBX;x++)   w_next.Put(x, 0, z, 0,w_next.Get(x, 1, z, 0));
    for (z=1;z<NBZ-1;z++) for (x=0;x<NBX;x++)   w_next.Put(x, NBY-1, z, 0,w_next.Get(x, NBY-2, z, 0));
    for (y=0;y<NBY;y++)   for (x=0;x<NBX;x++)   w_next.Put(x, y, 0, 0,w_next.Get(x, y, 1, 0));
    for (y=0;y<NBY;y++)   for (x=0;x<NBX;x++)   w_next.Put(x, y, NBZ-1, 0,w_next.Get(x, y, NBZ-2, 0));

    //next iteration becomes current iteration
    for (z=0;z<NBZ;z++) for (y=0;y<NBY;y++) for (x=0;x<NBX;x++) u_cur.Put(x, y, z, 0,u_next.Get(x, y, z, 0));
    for (z=0;z<NBZ;z++) for (y=0;y<NBY;y++) for (x=0;x<NBX;x++) v_cur.Put(x, y, z, 0,v_next.Get(x, y, z, 0));
    for (z=0;z<NBZ;z++) for (y=0;y<NBY;y++) for (x=0;x<NBX;x++) w_cur.Put(x, y, z, 0,w_next.Get(x, y, z, 0));
  }
  
  //write the result
  u_cur.Write("u.nii");
  v_cur.Write("v.nii");
  w_cur.Write("w.nii");
  
}*/





template <class VoxelType> void anisoDiffusion<VoxelType>::Run_3D_Explicit(){
  int NBX,NBY,NBZ,x,y,z;
  char Grad_f_x_file[] = "gradientVecX.nii";
  char Grad_f_y_file[] = "gradientVecY.nii";
  char Grad_f_z_file[] = "gradientVecZ.nii";
  irtkGenericImage<float> Grad_f_x;     //input gradient of the image I on x
  irtkGenericImage<float> Grad_f_y;     //input gradient of the image I on y
  irtkGenericImage<float> Grad_f_z;     //input gradient of the image I on z
  irtkGenericImage<float> SqNormGrad_f; //square norm of the input gradient of the image I
  irtkGenericImage<float> u_cur;        //output regularization of Grad_f_x
  irtkGenericImage<float> v_cur;        //output regularization of Grad_f_y
  irtkGenericImage<float> w_cur;        //output regularization of Grad_f_z
  irtkGenericImage<float> u_next;       //regularization of Grad_f_x at the iteration after the current one
  irtkGenericImage<float> v_next;       //regularization of Grad_f_y at the iteration after the current one
  irtkGenericImage<float> w_next;       //regularization of Grad_f_z at the iteration after the current one
  double tmpdbl1,tmpdbl2,tmpdbl3;
  float tmpf1,tmpf2,tmpf3,tmpf4;
  float mu;
  int IterationNb,it;
  float DeltaXsq,DeltaYsq,DeltaZsq,DeltaT;
  float A,B,C,D;
  float *Va;
  float *Vb; 
  float *Vc;
  float *Vd;
  float *Vx;
  int n;

  
  //1) init
  this->Initialize();
  
  //1.1) parameters
  mu=this->at;
  IterationNb=this->ITERATIONS_NB;
  DeltaT=this->dTau;
  DeltaXsq=1.;
  DeltaYsq=1.;
  DeltaZsq=1.;
  
  //1.2) precomputation of fixed values 
  A=(DeltaT*mu)/(3.*DeltaXsq);
  B=(DeltaT*mu)/(3.*DeltaYsq);
  C=(DeltaT*mu)/(3.*DeltaZsq);
  D=-DeltaT/3.;
  
  //1.3) read the gradients of f
  Grad_f_x.Read(Grad_f_x_file);
  Grad_f_y.Read(Grad_f_y_file);
  Grad_f_z.Read(Grad_f_z_file);
  
  //1.4) size of the images
  NBX=Grad_f_x.GetX();  //supposed the same in Grad_f_y and Grad_f_z
  NBY=Grad_f_x.GetY();  //supposed the same in Grad_f_y and Grad_f_z
  NBZ=Grad_f_x.GetZ();  //supposed the same in Grad_f_y and Grad_f_z
  
  //1.5) variables containing the regularization of grad f  (current iteration and next one)
  u_cur = irtkGenericImage<float>(NBX, NBY, NBZ, 1);
  v_cur = irtkGenericImage<float>(NBX, NBY, NBZ, 1);
  w_cur = irtkGenericImage<float>(NBX, NBY, NBZ, 1);
  u_next = irtkGenericImage<float>(NBX, NBY, NBZ, 1);
  v_next = irtkGenericImage<float>(NBX, NBY, NBZ, 1);
  w_next = irtkGenericImage<float>(NBX, NBY, NBZ, 1);
  
  //1.6) initialisation of u, v, w
  for (z=0;z<NBZ;z++) for (y=0;y<NBY;y++) for (x=0;x<NBX;x++) u_cur.Put(x, y, z, 0, (float)Grad_f_x.Get(x, y, z, 0));
  for (z=0;z<NBZ;z++) for (y=0;y<NBY;y++) for (x=0;x<NBX;x++) v_cur.Put(x, y, z, 0, (float)Grad_f_y.Get(x, y, z, 0));
  for (z=0;z<NBZ;z++) for (y=0;y<NBY;y++) for (x=0;x<NBX;x++) w_cur.Put(x, y, z, 0, (float)Grad_f_z.Get(x, y, z, 0));
  for (z=0;z<NBZ;z++) for (y=0;y<NBY;y++) for (x=0;x<NBX;x++) u_next.Put(x, y, z, 0, (float)Grad_f_x.Get(x, y, z, 0));
  for (z=0;z<NBZ;z++) for (y=0;y<NBY;y++) for (x=0;x<NBX;x++) v_next.Put(x, y, z, 0, (float)Grad_f_y.Get(x, y, z, 0));
  for (z=0;z<NBZ;z++) for (y=0;y<NBY;y++) for (x=0;x<NBX;x++) w_next.Put(x, y, z, 0, (float)Grad_f_z.Get(x, y, z, 0));
  
  //1.7) compute once for all the square norm of Grad_f
  SqNormGrad_f = irtkGenericImage<float>(NBX, NBY, NBZ, 1);
  for (z=0;z<NBZ;z++) for (y=0;y<NBY;y++) for (x=0;x<NBX;x++){
    tmpdbl1=pow((double)Grad_f_x.Get(x, y, z, 0),2.0);
    tmpdbl2=pow((double)Grad_f_y.Get(x, y, z, 0),2.0);
    tmpdbl3=pow((double)Grad_f_z.Get(x, y, z, 0),2.0);
    SqNormGrad_f.Put(x, y, z, 0, (float)(tmpdbl1+tmpdbl2+tmpdbl3));
  }
  
  
  //1.8) recommanded order of values for DeltaT and mu
  tmpdbl1=0;
  for (z=0;z<NBZ;z++) for (y=0;y<NBY;y++) for (x=0;x<NBX;x++)
    if (tmpdbl1<fabs(SqNormGrad_f.Get(x, y, z, 0)))
          tmpdbl1=fabs(SqNormGrad_f.Get(x, y, z, 0));
  
  
  
  
  cout << "Usage: AnisoDiff toto.nii toto.nii -SemiImplicit 0 -TimeDependent 1 -dTau [Delta T] -at [mu] -iterations [Iterations number]\n";
  cout << "Recommanded order of value for Delta T: " << 1./tmpdbl1 << "\n";
  cout << "Recommanded order of value for mu: " << tmpdbl1/100. << "\n";
  
  //1.9) temporary variables dedicated to the semi implicit scheme
  n=max(max(NBX,NBY),NBZ)+4; //for boundary effects
  Va=(float*)malloc(n*sizeof(float));
  Vb=(float*)malloc(n*sizeof(float));
  Vc=(float*)malloc(n*sizeof(float));
  Vd=(float*)malloc(n*sizeof(float));
  Vx=(float*)malloc(n*sizeof(float));

  
  //2) resolution
  for (it=0;it<IterationNb;it++){
    cout << "Iteration " << it << "\n";
    
    //2.1) 1st substep - x implicit / y,z explicit
    //2.1.1) compute the vector field of next iteration...
    //2.1.1.a) upate on u
    for (z=1;z<NBZ-1;z++) for (y=1;y<NBY-1;y++){
      for (x=1;x<NBX-1;x++){
        Vb[x+1]=1+2*A;
        Va[x+1]=-A;
        Vc[x+1]=-A;
        tmpf1=A*(u_cur.Get(x+1, y, z, 0)-2*u_cur.Get(x, y, z, 0)+u_cur.Get(x-1, y, z, 0));
        tmpf2=B*(u_cur.Get(x, y+1, z, 0)-2*u_cur.Get(x, y, z, 0)+u_cur.Get(x, y-1, z, 0));
        tmpf3=C*(u_cur.Get(x, y, z+1, 0)-2*u_cur.Get(x, y, z, 0)+u_cur.Get(x, y, z-1, 0));
        tmpf4=D*(u_cur.Get(x, y, z, 0)-Grad_f_x.Get(x, y, z, 0))*SqNormGrad_f.Get(x, y, z, 0);
        Vd[x+1]=u_cur.Get(x, y, z, 0)+tmpf2+tmpf3+tmpf4;
      }
      Va[1]=Va[3]; Va[0]=Va[4]; Va[NBX]=Va[NBX-2]; Va[NBX+1]=Va[NBX-3]; //to avoid boundary effects
      Vb[1]=Vb[3]; Vb[0]=Vb[4]; Vb[NBX]=Vb[NBX-2]; Vb[NBX+1]=Vb[NBX-3]; //to avoid boundary effects
      Vc[1]=Vc[3]; Vc[0]=Vc[4]; Vc[NBX]=Vc[NBX-2]; Vc[NBX+1]=Vc[NBX-3]; //to avoid boundary effects
      Vd[1]=Vd[3]; Vd[0]=Vd[4]; Vd[NBX]=Vd[NBX-2]; Vd[NBX+1]=Vd[NBX-3]; //to avoid boundary effects
      TridiagonalSolveFloat(Va,Vb,Vc,Vd,Vx,NBX+2);
      for (x = 1; x < NBX-1; x++) u_next.Put(x, y, z, 0,Vx[x+1]);
    }
    
    //2.1.1.b) upate on v
    for (z=1;z<NBZ-1;z++) for (y=1;y<NBY-1;y++){
      for (x=1;x<NBX-1;x++){
        Vb[x+1]=1+2*A;
        Va[x+1]=-A;
        Vc[x+1]=-A;
        tmpf1=A*(v_cur.Get(x+1, y, z, 0)-2*v_cur.Get(x, y, z, 0)+v_cur.Get(x-1, y, z, 0));
        tmpf2=B*(v_cur.Get(x, y+1, z, 0)-2*v_cur.Get(x, y, z, 0)+v_cur.Get(x, y-1, z, 0));
        tmpf3=C*(v_cur.Get(x, y, z+1, 0)-2*v_cur.Get(x, y, z, 0)+v_cur.Get(x, y, z-1, 0));
        tmpf4=D*(v_cur.Get(x, y, z, 0)-Grad_f_y.Get(x, y, z, 0))*SqNormGrad_f.Get(x, y, z, 0);
        Vd[x+1]=v_cur.Get(x, y, z, 0)+tmpf2+tmpf3+tmpf4;
      }
      Va[1]=Va[3]; Va[0]=Va[4]; Va[NBX]=Va[NBX-2]; Va[NBX+1]=Va[NBX-3]; //to avoid boundary effects
      Vb[1]=Vb[3]; Vb[0]=Vb[4]; Vb[NBX]=Vb[NBX-2]; Vb[NBX+1]=Vb[NBX-3]; //to avoid boundary effects
      Vc[1]=Vc[3]; Vc[0]=Vc[4]; Vc[NBX]=Vc[NBX-2]; Vc[NBX+1]=Vc[NBX-3]; //to avoid boundary effects
      Vd[1]=Vd[3]; Vd[0]=Vd[4]; Vd[NBX]=Vd[NBX-2]; Vd[NBX+1]=Vd[NBX-3]; //to avoid boundary effects
      TridiagonalSolveFloat(Va,Vb,Vc,Vd,Vx,NBX+2);
      for (x = 1; x < NBX-1; x++) v_next.Put(x, y, z, 0,Vx[x+1]);
    }
    
    //2.1.1.c) upate on w
    for (z=1;z<NBZ-1;z++) for (y=1;y<NBY-1;y++){
      for (x=1;x<NBX-1;x++){
        Vb[x+1]=1+2*A;
        Va[x+1]=-A;
        Vc[x+1]=-A;
        tmpf1=A*(w_cur.Get(x+1, y, z, 0)-2*w_cur.Get(x, y, z, 0)+w_cur.Get(x-1, y, z, 0));
        tmpf2=B*(w_cur.Get(x, y+1, z, 0)-2*w_cur.Get(x, y, z, 0)+w_cur.Get(x, y-1, z, 0));
        tmpf3=C*(w_cur.Get(x, y, z+1, 0)-2*w_cur.Get(x, y, z, 0)+w_cur.Get(x, y, z-1, 0));
        tmpf4=D*(w_cur.Get(x, y, z, 0)-Grad_f_z.Get(x, y, z, 0))*SqNormGrad_f.Get(x, y, z, 0);
        Vd[x+1]=w_cur.Get(x, y, z, 0)+tmpf2+tmpf3+tmpf4;
      }
      Va[1]=Va[3]; Va[0]=Va[4]; Va[NBX]=Va[NBX-2]; Va[NBX+1]=Va[NBX-3]; //to avoid boundary effects
      Vb[1]=Vb[3]; Vb[0]=Vb[4]; Vb[NBX]=Vb[NBX-2]; Vb[NBX+1]=Vb[NBX-3]; //to avoid boundary effects
      Vc[1]=Vc[3]; Vc[0]=Vc[4]; Vc[NBX]=Vc[NBX-2]; Vc[NBX+1]=Vc[NBX-3]; //to avoid boundary effects
      Vd[1]=Vd[3]; Vd[0]=Vd[4]; Vd[NBX]=Vd[NBX-2]; Vd[NBX+1]=Vd[NBX-3]; //to avoid boundary effects
      TridiagonalSolveFloat(Va,Vb,Vc,Vd,Vx,NBX+2);
      for (x = 1; x < NBX-1; x++) w_next.Put(x, y, z, 0,Vx[x+1]);
    }
    
    //2.1.2) boundary conditions
    for (z=1;z<NBZ-1;z++) for (y=1;y<NBY-1;y++) u_next.Put(0, y, z, 0,u_next.Get(1, y, z, 0));
    for (z=1;z<NBZ-1;z++) for (y=1;y<NBY-1;y++) u_next.Put(NBX-1, y, z, 0,u_next.Get(NBX-2, y, z, 0));
    for (z=1;z<NBZ-1;z++) for (x=0;x<NBX;x++)   u_next.Put(x, 0, z, 0,u_next.Get(x, 1, z, 0));
    for (z=1;z<NBZ-1;z++) for (x=0;x<NBX;x++)   u_next.Put(x, NBY-1, z, 0,u_next.Get(x, NBY-2, z, 0));
    for (y=0;y<NBY;y++)   for (x=0;x<NBX;x++)   u_next.Put(x, y, 0, 0,u_next.Get(x, y, 1, 0));
    for (y=0;y<NBY;y++)   for (x=0;x<NBX;x++)   u_next.Put(x, y, NBZ-1, 0,u_next.Get(x, y, NBZ-2, 0));
    
    for (z=1;z<NBZ-1;z++) for (y=1;y<NBY-1;y++) v_next.Put(0, y, z, 0,v_next.Get(1, y, z, 0));
    for (z=1;z<NBZ-1;z++) for (y=1;y<NBY-1;y++) v_next.Put(NBX-1, y, z, 0,v_next.Get(NBX-2, y, z, 0));
    for (z=1;z<NBZ-1;z++) for (x=0;x<NBX;x++)   v_next.Put(x, 0, z, 0,v_next.Get(x, 1, z, 0));
    for (z=1;z<NBZ-1;z++) for (x=0;x<NBX;x++)   v_next.Put(x, NBY-1, z, 0,v_next.Get(x, NBY-2, z, 0));
    for (y=0;y<NBY;y++)   for (x=0;x<NBX;x++)   v_next.Put(x, y, 0, 0,v_next.Get(x, y, 1, 0));
    for (y=0;y<NBY;y++)   for (x=0;x<NBX;x++)   v_next.Put(x, y, NBZ-1, 0,v_next.Get(x, y, NBZ-2, 0));

    for (z=1;z<NBZ-1;z++) for (y=1;y<NBY-1;y++) w_next.Put(0, y, z, 0,w_next.Get(1, y, z, 0));
    for (z=1;z<NBZ-1;z++) for (y=1;y<NBY-1;y++) w_next.Put(NBX-1, y, z, 0,w_next.Get(NBX-2, y, z, 0));
    for (z=1;z<NBZ-1;z++) for (x=0;x<NBX;x++)   w_next.Put(x, 0, z, 0,w_next.Get(x, 1, z, 0));
    for (z=1;z<NBZ-1;z++) for (x=0;x<NBX;x++)   w_next.Put(x, NBY-1, z, 0,w_next.Get(x, NBY-2, z, 0));
    for (y=0;y<NBY;y++)   for (x=0;x<NBX;x++)   w_next.Put(x, y, 0, 0,w_next.Get(x, y, 1, 0));
    for (y=0;y<NBY;y++)   for (x=0;x<NBX;x++)   w_next.Put(x, y, NBZ-1, 0,w_next.Get(x, y, NBZ-2, 0));

    //2.1.3) next iteration becomes current iteration
    for (z=0;z<NBZ;z++) for (y=0;y<NBY;y++) for (x=0;x<NBX;x++) u_cur.Put(x, y, z, 0,u_next.Get(x, y, z, 0));
    for (z=0;z<NBZ;z++) for (y=0;y<NBY;y++) for (x=0;x<NBX;x++) v_cur.Put(x, y, z, 0,v_next.Get(x, y, z, 0));
    for (z=0;z<NBZ;z++) for (y=0;y<NBY;y++) for (x=0;x<NBX;x++) w_cur.Put(x, y, z, 0,w_next.Get(x, y, z, 0));
  
  
  
    //2.2) 2nd substep - y implicit / x,z explicit
    //2.2.1) compute the vector field of next iteration...
    //2.2.1.a) upate on u
    for (z=1;z<NBZ-1;z++) for (x=1;x<NBX-1;x++){
      for (y=1;y<NBY-1;y++){
        Vb[y+1]=1+2*B;
        Va[y+1]=-B;
        Vc[y+1]=-B;
        tmpf1=A*(u_cur.Get(x+1, y, z, 0)-2*u_cur.Get(x, y, z, 0)+u_cur.Get(x-1, y, z, 0));
        tmpf2=B*(u_cur.Get(x, y+1, z, 0)-2*u_cur.Get(x, y, z, 0)+u_cur.Get(x, y-1, z, 0));
        tmpf3=C*(u_cur.Get(x, y, z+1, 0)-2*u_cur.Get(x, y, z, 0)+u_cur.Get(x, y, z-1, 0));
        tmpf4=D*(u_cur.Get(x, y, z, 0)-Grad_f_x.Get(x, y, z, 0))*SqNormGrad_f.Get(x, y, z, 0);
        Vd[y+1]=u_cur.Get(x, y, z, 0)+tmpf1+tmpf3+tmpf4;
      }
      Va[1]=Va[3]; Va[0]=Va[4]; Va[NBY]=Va[NBY-2]; Va[NBY+1]=Va[NBY-3]; //to avoid boundary effects
      Vb[1]=Vb[3]; Vb[0]=Vb[4]; Vb[NBY]=Vb[NBY-2]; Vb[NBY+1]=Vb[NBY-3]; //to avoid boundary effects
      Vc[1]=Vc[3]; Vc[0]=Vc[4]; Vc[NBY]=Vc[NBY-2]; Vc[NBY+1]=Vc[NBY-3]; //to avoid boundary effects
      Vd[1]=Vd[3]; Vd[0]=Vd[4]; Vd[NBY]=Vd[NBY-2]; Vd[NBY+1]=Vd[NBY-3]; //to avoid boundary effects
      TridiagonalSolveFloat(Va,Vb,Vc,Vd,Vx,NBY+2);
      for (y=1;y<NBY-1;y++) u_next.Put(x, y, z, 0,Vx[y+1]);
    }
    
    //2.2.1.b) upate on v
    for (z=1;z<NBZ-1;z++) for (x=1;x<NBX-1;x++){
      for (y=1;y<NBY-1;y++){
        Vb[y+1]=1+2*B;
        Va[y+1]=-B;
        Vc[y+1]=-B;
        tmpf1=A*(v_cur.Get(x+1, y, z, 0)-2*v_cur.Get(x, y, z, 0)+v_cur.Get(x-1, y, z, 0));
        tmpf2=B*(v_cur.Get(x, y+1, z, 0)-2*v_cur.Get(x, y, z, 0)+v_cur.Get(x, y-1, z, 0));
        tmpf3=C*(v_cur.Get(x, y, z+1, 0)-2*v_cur.Get(x, y, z, 0)+v_cur.Get(x, y, z-1, 0));
        tmpf4=D*(v_cur.Get(x, y, z, 0)-Grad_f_y.Get(x, y, z, 0))*SqNormGrad_f.Get(x, y, z, 0);
        Vd[y+1]=v_cur.Get(x, y, z, 0)+tmpf1+tmpf3+tmpf4;
      }
      Va[1]=Va[3]; Va[0]=Va[4]; Va[NBY]=Va[NBY-2]; Va[NBY+1]=Va[NBY-3]; //to avoid boundary effects
      Vb[1]=Vb[3]; Vb[0]=Vb[4]; Vb[NBY]=Vb[NBY-2]; Vb[NBY+1]=Vb[NBY-3]; //to avoid boundary effects
      Vc[1]=Vc[3]; Vc[0]=Vc[4]; Vc[NBY]=Vc[NBY-2]; Vc[NBY+1]=Vc[NBY-3]; //to avoid boundary effects
      Vd[1]=Vd[3]; Vd[0]=Vd[4]; Vd[NBY]=Vd[NBY-2]; Vd[NBY+1]=Vd[NBY-3]; //to avoid boundary effects
      TridiagonalSolveFloat(Va,Vb,Vc,Vd,Vx,NBY+2);
      for (y=1;y<NBY-1;y++) v_next.Put(x, y, z, 0,Vx[y+1]);
    }
    
    //2.2.1.c) upate on w
    for (z=1;z<NBZ-1;z++) for (x=1;x<NBX-1;x++){
      for (y=1;y<NBY-1;y++){
        Vb[y+1]=1+2*B;
        Va[y+1]=-B;
        Vc[y+1]=-B;
        tmpf1=A*(w_cur.Get(x+1, y, z, 0)-2*w_cur.Get(x, y, z, 0)+w_cur.Get(x-1, y, z, 0));
        tmpf2=B*(w_cur.Get(x, y+1, z, 0)-2*w_cur.Get(x, y, z, 0)+w_cur.Get(x, y-1, z, 0));
        tmpf3=C*(w_cur.Get(x, y, z+1, 0)-2*w_cur.Get(x, y, z, 0)+w_cur.Get(x, y, z-1, 0));
        tmpf4=D*(w_cur.Get(x, y, z, 0)-Grad_f_z.Get(x, y, z, 0))*SqNormGrad_f.Get(x, y, z, 0);
        Vd[y+1]=w_cur.Get(x, y, z, 0)+tmpf1+tmpf3+tmpf4;
      }
      Va[1]=Va[3]; Va[0]=Va[4]; Va[NBY]=Va[NBY-2]; Va[NBY+1]=Va[NBY-3]; //to avoid boundary effects
      Vb[1]=Vb[3]; Vb[0]=Vb[4]; Vb[NBY]=Vb[NBY-2]; Vb[NBY+1]=Vb[NBY-3]; //to avoid boundary effects
      Vc[1]=Vc[3]; Vc[0]=Vc[4]; Vc[NBY]=Vc[NBY-2]; Vc[NBY+1]=Vc[NBY-3]; //to avoid boundary effects
      Vd[1]=Vd[3]; Vd[0]=Vd[4]; Vd[NBY]=Vd[NBY-2]; Vd[NBY+1]=Vd[NBY-3]; //to avoid boundary effects
      TridiagonalSolveFloat(Va,Vb,Vc,Vd,Vx,NBY+2);
      for (y=1;y<NBY-1;y++) w_next.Put(x, y, z, 0,Vx[y+1]);
    }
    
    //2.2.2) boundary conditions
    for (z=1;z<NBZ-1;z++) for (y=1;y<NBY-1;y++) u_next.Put(0, y, z, 0,u_next.Get(1, y, z, 0));
    for (z=1;z<NBZ-1;z++) for (y=1;y<NBY-1;y++) u_next.Put(NBX-1, y, z, 0,u_next.Get(NBX-2, y, z, 0));
    for (z=1;z<NBZ-1;z++) for (x=0;x<NBX;x++)   u_next.Put(x, 0, z, 0,u_next.Get(x, 1, z, 0));
    for (z=1;z<NBZ-1;z++) for (x=0;x<NBX;x++)   u_next.Put(x, NBY-1, z, 0,u_next.Get(x, NBY-2, z, 0));
    for (y=0;y<NBY;y++)   for (x=0;x<NBX;x++)   u_next.Put(x, y, 0, 0,u_next.Get(x, y, 1, 0));
    for (y=0;y<NBY;y++)   for (x=0;x<NBX;x++)   u_next.Put(x, y, NBZ-1, 0,u_next.Get(x, y, NBZ-2, 0));
    
    for (z=1;z<NBZ-1;z++) for (y=1;y<NBY-1;y++) v_next.Put(0, y, z, 0,v_next.Get(1, y, z, 0));
    for (z=1;z<NBZ-1;z++) for (y=1;y<NBY-1;y++) v_next.Put(NBX-1, y, z, 0,v_next.Get(NBX-2, y, z, 0));
    for (z=1;z<NBZ-1;z++) for (x=0;x<NBX;x++)   v_next.Put(x, 0, z, 0,v_next.Get(x, 1, z, 0));
    for (z=1;z<NBZ-1;z++) for (x=0;x<NBX;x++)   v_next.Put(x, NBY-1, z, 0,v_next.Get(x, NBY-2, z, 0));
    for (y=0;y<NBY;y++)   for (x=0;x<NBX;x++)   v_next.Put(x, y, 0, 0,v_next.Get(x, y, 1, 0));
    for (y=0;y<NBY;y++)   for (x=0;x<NBX;x++)   v_next.Put(x, y, NBZ-1, 0,v_next.Get(x, y, NBZ-2, 0));

    for (z=1;z<NBZ-1;z++) for (y=1;y<NBY-1;y++) w_next.Put(0, y, z, 0,w_next.Get(1, y, z, 0));
    for (z=1;z<NBZ-1;z++) for (y=1;y<NBY-1;y++) w_next.Put(NBX-1, y, z, 0,w_next.Get(NBX-2, y, z, 0));
    for (z=1;z<NBZ-1;z++) for (x=0;x<NBX;x++)   w_next.Put(x, 0, z, 0,w_next.Get(x, 1, z, 0));
    for (z=1;z<NBZ-1;z++) for (x=0;x<NBX;x++)   w_next.Put(x, NBY-1, z, 0,w_next.Get(x, NBY-2, z, 0));
    for (y=0;y<NBY;y++)   for (x=0;x<NBX;x++)   w_next.Put(x, y, 0, 0,w_next.Get(x, y, 1, 0));
    for (y=0;y<NBY;y++)   for (x=0;x<NBX;x++)   w_next.Put(x, y, NBZ-1, 0,w_next.Get(x, y, NBZ-2, 0));

    //2.2.3) next iteration becomes current iteration
    for (z=0;z<NBZ;z++) for (y=0;y<NBY;y++) for (x=0;x<NBX;x++) u_cur.Put(x, y, z, 0,u_next.Get(x, y, z, 0));
    for (z=0;z<NBZ;z++) for (y=0;y<NBY;y++) for (x=0;x<NBX;x++) v_cur.Put(x, y, z, 0,v_next.Get(x, y, z, 0));
    for (z=0;z<NBZ;z++) for (y=0;y<NBY;y++) for (x=0;x<NBX;x++) w_cur.Put(x, y, z, 0,w_next.Get(x, y, z, 0));
  
  
    //2.3) 1st substep - x implicit / y,z explicit
    //2.3.1) compute the vector field of next iteration...
    //2.3.1.a) upate on u
    for (y=1;y<NBY-1;y++) for (x=1;x<NBX-1;x++){
      for (z=1;z<NBZ-1;z++){
        Vb[z+1]=1+2*C;
        Va[z+1]=-C;
        Vc[z+1]=-C;
        tmpf1=A*(u_cur.Get(x+1, y, z, 0)-2*u_cur.Get(x, y, z, 0)+u_cur.Get(x-1, y, z, 0));
        tmpf2=B*(u_cur.Get(x, y+1, z, 0)-2*u_cur.Get(x, y, z, 0)+u_cur.Get(x, y-1, z, 0));
        tmpf3=C*(u_cur.Get(x, y, z+1, 0)-2*u_cur.Get(x, y, z, 0)+u_cur.Get(x, y, z-1, 0));
        tmpf4=D*(u_cur.Get(x, y, z, 0)-Grad_f_x.Get(x, y, z, 0))*SqNormGrad_f.Get(x, y, z, 0);
        Vd[z+1]=u_cur.Get(x, y, z, 0)+tmpf1+tmpf2+tmpf4;
      }
      Va[1]=Va[3]; Va[0]=Va[4]; Va[NBZ]=Va[NBZ-2]; Va[NBZ+1]=Va[NBZ-3]; //to avoid boundary effects
      Vb[1]=Vb[3]; Vb[0]=Vb[4]; Vb[NBZ]=Vb[NBZ-2]; Vb[NBZ+1]=Vb[NBZ-3]; //to avoid boundary effects
      Vc[1]=Vc[3]; Vc[0]=Vc[4]; Vc[NBZ]=Vc[NBZ-2]; Vc[NBZ+1]=Vc[NBZ-3]; //to avoid boundary effects
      Vd[1]=Vd[3]; Vd[0]=Vd[4]; Vd[NBZ]=Vd[NBZ-2]; Vd[NBZ+1]=Vd[NBZ-3]; //to avoid boundary effects
      TridiagonalSolveFloat(Va,Vb,Vc,Vd,Vx,NBZ+2);
      for (z=1;z<NBZ-1;z++) u_next.Put(x, y, z, 0,Vx[z+1]);
    }
    
    //2.3.1.b) upate on v
    for (y=1;y<NBY-1;y++) for (x=1;x<NBX-1;x++){
      for (z=1;z<NBZ-1;z++){
        Vb[z+1]=1+2*C;
        Va[z+1]=-C;
        Vc[z+1]=-C;
        tmpf1=A*(v_cur.Get(x+1, y, z, 0)-2*v_cur.Get(x, y, z, 0)+v_cur.Get(x-1, y, z, 0));
        tmpf2=B*(v_cur.Get(x, y+1, z, 0)-2*v_cur.Get(x, y, z, 0)+v_cur.Get(x, y-1, z, 0));
        tmpf3=C*(v_cur.Get(x, y, z+1, 0)-2*v_cur.Get(x, y, z, 0)+v_cur.Get(x, y, z-1, 0));
        tmpf4=D*(v_cur.Get(x, y, z, 0)-Grad_f_y.Get(x, y, z, 0))*SqNormGrad_f.Get(x, y, z, 0);
        Vd[z+1]=v_cur.Get(x, y, z, 0)+tmpf1+tmpf2+tmpf4;
      }
      Va[1]=Va[3]; Va[0]=Va[4]; Va[NBZ]=Va[NBZ-2]; Va[NBZ+1]=Va[NBZ-3]; //to avoid boundary effects
      Vb[1]=Vb[3]; Vb[0]=Vb[4]; Vb[NBZ]=Vb[NBZ-2]; Vb[NBZ+1]=Vb[NBZ-3]; //to avoid boundary effects
      Vc[1]=Vc[3]; Vc[0]=Vc[4]; Vc[NBZ]=Vc[NBZ-2]; Vc[NBZ+1]=Vc[NBZ-3]; //to avoid boundary effects
      Vd[1]=Vd[3]; Vd[0]=Vd[4]; Vd[NBZ]=Vd[NBZ-2]; Vd[NBZ+1]=Vd[NBZ-3]; //to avoid boundary effects
      TridiagonalSolveFloat(Va,Vb,Vc,Vd,Vx,NBZ+2);
      for (z=1;z<NBZ-1;z++) v_next.Put(x, y, z, 0,Vx[z+1]);
    }
    
    //2.3.1.c) upate on w
    for (y=1;y<NBY-1;y++) for (x=1;x<NBX-1;x++){
      for (z=1;z<NBZ-1;z++){
        Vb[z+1]=1+2*C;
        Va[z+1]=-C;
        Vc[z+1]=-C;
        tmpf1=A*(w_cur.Get(x+1, y, z, 0)-2*w_cur.Get(x, y, z, 0)+w_cur.Get(x-1, y, z, 0));
        tmpf2=B*(w_cur.Get(x, y+1, z, 0)-2*w_cur.Get(x, y, z, 0)+w_cur.Get(x, y-1, z, 0));
        tmpf3=C*(w_cur.Get(x, y, z+1, 0)-2*w_cur.Get(x, y, z, 0)+w_cur.Get(x, y, z-1, 0));
        tmpf4=D*(w_cur.Get(x, y, z, 0)-Grad_f_z.Get(x, y, z, 0))*SqNormGrad_f.Get(x, y, z, 0);
        Vd[z+1]=w_cur.Get(x, y, z, 0)+tmpf1+tmpf2+tmpf4;
      }
      Va[1]=Va[3]; Va[0]=Va[4]; Va[NBZ]=Va[NBZ-2]; Va[NBZ+1]=Va[NBZ-3]; //to avoid boundary effects
      Vb[1]=Vb[3]; Vb[0]=Vb[4]; Vb[NBZ]=Vb[NBZ-2]; Vb[NBZ+1]=Vb[NBZ-3]; //to avoid boundary effects
      Vc[1]=Vc[3]; Vc[0]=Vc[4]; Vc[NBZ]=Vc[NBZ-2]; Vc[NBZ+1]=Vc[NBZ-3]; //to avoid boundary effects
      Vd[1]=Vd[3]; Vd[0]=Vd[4]; Vd[NBZ]=Vd[NBZ-2]; Vd[NBZ+1]=Vd[NBZ-3]; //to avoid boundary effects
      TridiagonalSolveFloat(Va,Vb,Vc,Vd,Vx,NBZ+2);
      for (z=1;z<NBZ-1;z++) w_next.Put(x, y, z, 0,Vx[z+1]);
    }
    
    //2.3.2) boundary conditions
    for (z=1;z<NBZ-1;z++) for (y=1;y<NBY-1;y++) u_next.Put(0, y, z, 0,u_next.Get(1, y, z, 0));
    for (z=1;z<NBZ-1;z++) for (y=1;y<NBY-1;y++) u_next.Put(NBX-1, y, z, 0,u_next.Get(NBX-2, y, z, 0));
    for (z=1;z<NBZ-1;z++) for (x=0;x<NBX;x++)   u_next.Put(x, 0, z, 0,u_next.Get(x, 1, z, 0));
    for (z=1;z<NBZ-1;z++) for (x=0;x<NBX;x++)   u_next.Put(x, NBY-1, z, 0,u_next.Get(x, NBY-2, z, 0));
    for (y=0;y<NBY;y++)   for (x=0;x<NBX;x++)   u_next.Put(x, y, 0, 0,u_next.Get(x, y, 1, 0));
    for (y=0;y<NBY;y++)   for (x=0;x<NBX;x++)   u_next.Put(x, y, NBZ-1, 0,u_next.Get(x, y, NBZ-2, 0));
    
    for (z=1;z<NBZ-1;z++) for (y=1;y<NBY-1;y++) v_next.Put(0, y, z, 0,v_next.Get(1, y, z, 0));
    for (z=1;z<NBZ-1;z++) for (y=1;y<NBY-1;y++) v_next.Put(NBX-1, y, z, 0,v_next.Get(NBX-2, y, z, 0));
    for (z=1;z<NBZ-1;z++) for (x=0;x<NBX;x++)   v_next.Put(x, 0, z, 0,v_next.Get(x, 1, z, 0));
    for (z=1;z<NBZ-1;z++) for (x=0;x<NBX;x++)   v_next.Put(x, NBY-1, z, 0,v_next.Get(x, NBY-2, z, 0));
    for (y=0;y<NBY;y++)   for (x=0;x<NBX;x++)   v_next.Put(x, y, 0, 0,v_next.Get(x, y, 1, 0));
    for (y=0;y<NBY;y++)   for (x=0;x<NBX;x++)   v_next.Put(x, y, NBZ-1, 0,v_next.Get(x, y, NBZ-2, 0));

    for (z=1;z<NBZ-1;z++) for (y=1;y<NBY-1;y++) w_next.Put(0, y, z, 0,w_next.Get(1, y, z, 0));
    for (z=1;z<NBZ-1;z++) for (y=1;y<NBY-1;y++) w_next.Put(NBX-1, y, z, 0,w_next.Get(NBX-2, y, z, 0));
    for (z=1;z<NBZ-1;z++) for (x=0;x<NBX;x++)   w_next.Put(x, 0, z, 0,w_next.Get(x, 1, z, 0));
    for (z=1;z<NBZ-1;z++) for (x=0;x<NBX;x++)   w_next.Put(x, NBY-1, z, 0,w_next.Get(x, NBY-2, z, 0));
    for (y=0;y<NBY;y++)   for (x=0;x<NBX;x++)   w_next.Put(x, y, 0, 0,w_next.Get(x, y, 1, 0));
    for (y=0;y<NBY;y++)   for (x=0;x<NBX;x++)   w_next.Put(x, y, NBZ-1, 0,w_next.Get(x, y, NBZ-2, 0));

    //2.3.3) next iteration becomes current iteration
    for (z=0;z<NBZ;z++) for (y=0;y<NBY;y++) for (x=0;x<NBX;x++) u_cur.Put(x, y, z, 0,u_next.Get(x, y, z, 0));
    for (z=0;z<NBZ;z++) for (y=0;y<NBY;y++) for (x=0;x<NBX;x++) v_cur.Put(x, y, z, 0,v_next.Get(x, y, z, 0));
    for (z=0;z<NBZ;z++) for (y=0;y<NBY;y++) for (x=0;x<NBX;x++) w_cur.Put(x, y, z, 0,w_next.Get(x, y, z, 0));
  
  }
  
  //3) write the result
  u_cur.Write("u.nii");
  v_cur.Write("v.nii");
  w_cur.Write("w.nii");
  
}








/*
  //2.2.1) diffusion - x implicit / y,z,t explicit
  //2.2.1.1) explicit part
  for (t = 1; t < NBT-1; t++) for (z = 1; z < NBZ-1; z++) for (y = 1; y < NBY-1; y++) for (x = 1; x < NBX-1; x++) {
    dIdy=(imageE[t][z][y+1][x]-imageE[t][z][y-1][x])/(2*dy);
    Dyy_div_dySq=static_cast<float>((1-exp(-3.314/pow((dIdy/ay),4)))*DivPowDySqu);
    dIdz=(imageE[t][z+1][y][x]-imageE[t][z-1][y][x])/(2*dz);
    Dzz_div_dzSq=static_cast<float>((1-exp(-3.314/pow((dIdz/az),4)))*DivPowDzSqu);
    dIdt=(imageE[t+1][z][y][x]-imageE[t-1][z][y][x])/(2*dt);
    Dtt_div_dtSq=static_cast<float>((1-exp(-3.314/pow((dIdt/at),4)))*DivPowDzSqu);
          new value of the voxel
    DivDgradI=(imageE[t][z][y+1][x]-2*imageE[t][z][y][x]+imageE[t][z][y-1][x])*Dyy_div_dySq+
        (imageE[t][z+1][y][x]-2*imageE[t][z][y][x]+imageE[t][z-1][y][x])*Dzz_div_dzSq+
        (imageE[t+1][z][y][x]-2*imageE[t][z][y][x]+imageE[t-1][z][y][x])*Dtt_div_dtSq;

    imageO[t][z][y][x]=imageE[t][z][y][x]+(dTau/4.)*DivDgradI;
  }
  
  //2.2.1.2) implicit part
  for (t = 1; t < NBT-1; t++) for (z = 1; z < NBZ-1; z++) for (y = 1; y < NBY-1; y++) {
    for (x = 1; x < NBX-1; x++){
      dIdx=(imageE[t][z][y][x+1]-imageE[t][z][y][x-1])/(2*dx);
      Dxx_div_dxSq=static_cast<float>((1-exp(-3.314/pow((dIdx/ax),4)))*DivPowDxSqu);
      Va[x+1]=(dTau/4.)*Dxx_div_dxSq;
      Vb[x+1]=-1-2*(dTau/4.)*Dxx_div_dxSq;
      Vc[x+1]=(dTau/4.)*Dxx_div_dxSq;
      Vd[x+1]=imageE[t][z][y][x];
    }
    Va[1]=Va[3]; Va[0]=Va[4]; Va[NBX]=Va[NBX-2]; Va[NBX+1]=Va[NBX-3]; //to avoid boundary effects
    Vb[1]=Vb[3]; Vb[0]=Vb[4]; Vb[NBX]=Vb[NBX-2]; Vb[NBX+1]=Vb[NBX-3]; //to avoid boundary effects
    Vc[1]=Vc[3]; Vc[0]=Vc[4]; Vc[NBX]=Vc[NBX-2]; Vc[NBX+1]=Vc[NBX-3]; //to avoid boundary effects
    Vd[1]=Vd[3]; Vd[0]=Vd[4]; Vd[NBX]=Vd[NBX-2]; Vd[NBX+1]=Vd[NBX-3]; //to avoid boundary effects
    TridiagonalSolveFloat(Va,Vb,Vc,Vd,Vx,NBX+2);
    for (x = 1; x < NBX-1; x++) imageO[t][z][y][x]=-Vx[x+1];
  }
*/


///+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
///+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
///+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
///                                        END GRADIENT VECTOR FLOW PROJECT
///+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
///+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
///+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++




template <class VoxelType> void anisoDiffusion<VoxelType>::Run()
{
	if ((this->TimeDependent==true) && (this->_input->GetT()>4)){
		if (this->SemiImplicit==true)
			Run_4D_semiImplicit();
		else
			Run_4D_Explicit();
	}
	else{
		if (this->SemiImplicit==true)
			Run_3D_semiImplicit();
		else
			Run_3D_Explicit();
	}
	
}







template class anisoDiffusion<irtkBytePixel>;
template class anisoDiffusion<irtkGreyPixel>;
template class anisoDiffusion<irtkRealPixel>;








/* Solve the problem: MX=D where D is a known vector, M a tridiagonal matrix and X the unknown vector.
Inputs are a,b,c,d,n where M(i,i)=b(i), M(i,i-1)=a(i), M(i,i+1)=c(i), D(i)=d(i), D in R^n and M in R^n*R^n.
Output is X where X in R^n.  Warning: will modify c and d! */
void TridiagonalSolveFloat(const float *a, const float *b, float *c, float *d, float *x, int n){
  int i;
  double id;

  /* Modify the coefficients. */
  c[0] /= b[0];                       /* Division by zero risk. */
  d[0] /= b[0];                       /* Division by zero would imply a singular matrix. */
  for(i = 1; i < n; i++){
    id = (b[i] - c[i-1] * a[i]);      /* Division by zero risk. */
    c[i] /= id;                       /* Last value calculated is redundant. */
    d[i] = (d[i] - d[i-1] * a[i])/id;
  }
  
  /* Now back substitute. */
  x[n - 1] = d[n - 1];
  for(i = n - 2; i >= 0; i--)
    x[i] = d[i] - c[i] * x[i + 1];
}
