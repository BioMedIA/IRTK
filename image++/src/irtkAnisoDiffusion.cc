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
  TimeDependent=False;
  SemiImplicit=True;
}

template <class VoxelType> anisoDiffusion<VoxelType>::~anisoDiffusion(void)
{}

template <class VoxelType> Bool anisoDiffusion<VoxelType>::RequiresBuffering(void)
{
  return True;
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
					Vd[x+1]=imageE[z][y][x];
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


template <class VoxelType> void anisoDiffusion<VoxelType>::Run()
{
	if ((this->TimeDependent==True) && (this->_input->GetT()>4)){
		if (this->SemiImplicit==True)
			Run_4D_semiImplicit();
		else
			Run_4D_Explicit();
	}
	else{
		if (this->SemiImplicit==True)
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
	c[0] /= b[0];				/* Division by zero risk. */
	d[0] /= b[0];				/* Division by zero would imply a singular matrix. */
	for(i = 1; i < n; i++){
		id = (b[i] - c[i-1] * a[i]);	/* Division by zero risk. */
		c[i] /= id;				/* Last value calculated is redundant. */
		d[i] = (d[i] - d[i-1] * a[i])/id;
	}
 
	/* Now back substitute. */
	x[n - 1] = d[n - 1];
	for(i = n - 2; i >= 0; i--)
		x[i] = d[i] - c[i] * x[i + 1];
}
