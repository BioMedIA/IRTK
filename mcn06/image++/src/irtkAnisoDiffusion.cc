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

typedef struct
{
  double X;
  double Y;     // coordinates of a vertex
  double Z;
} Vertex;

typedef struct
{
  int V1;   //id of the first vertex
  int V2;   //id of the second vertex
  Vertex * SubiVertex;
  int NbSubdiv;
} Edge;

typedef struct
{
  int E1;    //id of the first edge
  int E2;    //id of the second edge
  int E3;    //id of the third edge
} Element;

typedef struct
{
  int NbVertexes;
  Vertex * Vertexes;
  int NbEdges;
  Edge * Edges;
  int NbElements;
  Element * Elements;
} Mesh;





//create a simple cubic mesh 
Mesh CreateBasicCubicMesh(){
  Mesh LocMesh;
  int i;
	
	
	//allocations...
  LocMesh.NbVertexes=8;
  LocMesh.Vertexes=(Vertex *)malloc(LocMesh.NbVertexes*sizeof(Vertex));
	
  LocMesh.NbEdges=18;
  LocMesh.Edges=(Edge *)malloc(LocMesh.NbEdges*sizeof(Edge));
        
  LocMesh.NbElements=12;
  LocMesh.Elements=(Element *)malloc(LocMesh.NbElements*sizeof(Element));
	
	//fill the coordinates of each point
  LocMesh.Vertexes[0].X=0; LocMesh.Vertexes[0].Y=0; LocMesh.Vertexes[0].Z=0;
  LocMesh.Vertexes[1].X=1; LocMesh.Vertexes[1].Y=0; LocMesh.Vertexes[1].Z=0;
  LocMesh.Vertexes[2].X=1; LocMesh.Vertexes[2].Y=1; LocMesh.Vertexes[2].Z=0;
  LocMesh.Vertexes[3].X=0; LocMesh.Vertexes[3].Y=1; LocMesh.Vertexes[3].Z=0;
  LocMesh.Vertexes[4].X=0; LocMesh.Vertexes[4].Y=0; LocMesh.Vertexes[4].Z=1;
  LocMesh.Vertexes[5].X=1; LocMesh.Vertexes[5].Y=0; LocMesh.Vertexes[5].Z=1;
  LocMesh.Vertexes[6].X=1; LocMesh.Vertexes[6].Y=1; LocMesh.Vertexes[6].Z=1;
  LocMesh.Vertexes[7].X=0; LocMesh.Vertexes[7].Y=1; LocMesh.Vertexes[7].Z=1;
	
	//fill each edge
  LocMesh.Edges[0].V1=0;  LocMesh.Edges[0].V2=1;
  LocMesh.Edges[1].V1=1;  LocMesh.Edges[1].V2=2;
  LocMesh.Edges[2].V1=2;  LocMesh.Edges[2].V2=3;
  LocMesh.Edges[3].V1=3;  LocMesh.Edges[3].V2=0;
  LocMesh.Edges[4].V1=4;  LocMesh.Edges[4].V2=5;
  LocMesh.Edges[5].V1=5;  LocMesh.Edges[5].V2=6;
  LocMesh.Edges[6].V1=6;  LocMesh.Edges[6].V2=7;
  LocMesh.Edges[7].V1=7;  LocMesh.Edges[7].V2=4;
  LocMesh.Edges[8].V1=0;  LocMesh.Edges[8].V2=4;
  LocMesh.Edges[9].V1=1;  LocMesh.Edges[9].V2=5;
  LocMesh.Edges[10].V1=2;  LocMesh.Edges[10].V2=6;
  LocMesh.Edges[11].V1=3;  LocMesh.Edges[11].V2=7;
  LocMesh.Edges[12].V1=1;  LocMesh.Edges[12].V2=4;
  LocMesh.Edges[13].V1=1;  LocMesh.Edges[13].V2=3;
  LocMesh.Edges[14].V1=1;  LocMesh.Edges[14].V2=6;
  LocMesh.Edges[15].V1=7;  LocMesh.Edges[15].V2=0;
  LocMesh.Edges[16].V1=7;  LocMesh.Edges[16].V2=2;
  LocMesh.Edges[17].V1=7;  LocMesh.Edges[17].V2=5;
        
  for (i=0;i<LocMesh.NbEdges;i++) LocMesh.Edges[i].NbSubdiv=1;
	
	//fill each Element
  LocMesh.Elements[0].E1=0;  LocMesh.Elements[0].E2=8;  LocMesh.Elements[0].E3=12;
  LocMesh.Elements[1].E1=9;  LocMesh.Elements[1].E2=4;  LocMesh.Elements[1].E3=12;
  LocMesh.Elements[2].E1=0;  LocMesh.Elements[2].E2=3;  LocMesh.Elements[2].E3=13;
  LocMesh.Elements[3].E1=1;  LocMesh.Elements[3].E2=2;  LocMesh.Elements[3].E3=13;
  LocMesh.Elements[4].E1=1;  LocMesh.Elements[4].E2=10;  LocMesh.Elements[4].E3=14;
  LocMesh.Elements[5].E1=9;  LocMesh.Elements[5].E2=5;  LocMesh.Elements[5].E3=14;
  LocMesh.Elements[6].E1=3;  LocMesh.Elements[6].E2=11;  LocMesh.Elements[6].E3=15;
  LocMesh.Elements[7].E1=8;  LocMesh.Elements[7].E2=7;  LocMesh.Elements[7].E3=15;
  LocMesh.Elements[8].E1=2;  LocMesh.Elements[8].E2=11;  LocMesh.Elements[8].E3=16;
  LocMesh.Elements[9].E1=10;  LocMesh.Elements[9].E2=6;  LocMesh.Elements[9].E3=16;
  LocMesh.Elements[10].E1=4;  LocMesh.Elements[10].E2=7;  LocMesh.Elements[10].E3=17;
  LocMesh.Elements[11].E1=5;  LocMesh.Elements[11].E2=6;  LocMesh.Elements[11].E3=17;
	
  return LocMesh;
}



//smooth the coordinates of the mesh LocMesh
/*void SmoothMesh(Mesh * LocMesh){
  int IdVertexLoc,IdEdgeLoc;
  double FiltX,FiltY,FiltZ;
  double NbPts;
  double * NewCoordX;
  double * NewCoordY;
  double * NewCoordZ;

  if (LocMesh->Edges[0].NbSubdiv==1){
    //intialisation  
    NewCoordX=(double*)malloc(LocMesh->NbVertexes*sizeof(double));
    NewCoordY=(double*)malloc(LocMesh->NbVertexes*sizeof(double));
    NewCoordZ=(double*)malloc(LocMesh->NbVertexes*sizeof(double));
    
    //filtering
    for (IdVertexLoc=0;IdVertexLoc<LocMesh->NbVertexes;IdVertexLoc++){
      FiltX=LocMesh->Vertexes[IdVertexLoc].X;
      FiltY=LocMesh->Vertexes[IdVertexLoc].Y;
      FiltZ=LocMesh->Vertexes[IdVertexLoc].Z;
      NbPts=1;
      for (IdEdgeLoc=0;IdEdgeLoc<LocMesh->NbEdges;IdEdgeLoc++){
        if (LocMesh->Edges[IdEdgeLoc].V1=IdVertexLoc){
          FiltX+=LocMesh->Vertexes[LocMesh->Edges[IdEdgeLoc].V2].X;
          FiltY+=LocMesh->Vertexes[LocMesh->Edges[IdEdgeLoc].V2].Y;
          FiltZ+=LocMesh->Vertexes[LocMesh->Edges[IdEdgeLoc].V2].Z;
          NbPts++;
        }
        if (LocMesh->Edges[IdEdgeLoc].V2=IdVertexLoc){
          FiltX+=LocMesh->Vertexes[LocMesh->Edges[IdEdgeLoc].V1].X;
          FiltY+=LocMesh->Vertexes[LocMesh->Edges[IdEdgeLoc].V1].Y;
          FiltZ+=LocMesh->Vertexes[LocMesh->Edges[IdEdgeLoc].V1].Z;
          NbPts++;
        }
      }
      
      NewCoordX[IdVertexLoc]=FiltX/NbPts;
      NewCoordY[IdVertexLoc]=FiltY/NbPts;
      NewCoordZ[IdVertexLoc]=FiltZ/NbPts;
    }
    
    //finalize
    for (IdVertexLoc=0;IdVertexLoc<LocMesh->NbVertexes;IdVertexLoc++){
      LocMesh->Vertexes[IdVertexLoc].X=NewCoordX[IdVertexLoc];
      LocMesh->Vertexes[IdVertexLoc].Y=NewCoordX[IdVertexLoc];
      LocMesh->Vertexes[IdVertexLoc].Z=NewCoordX[IdVertexLoc];
    }
    
    //deallocation
    free(NewCoordX);
    free(NewCoordY);
    free(NewCoordZ);
  }
}*/

//subdivide into 2 edges the edge 'IdEdge' in the mesh 'LocMesh'
void CutEdge(Mesh * LocMesh,int IdEdge){
  int i;
  int IdNgbhVortex1;
  int IdNgbhVortex2;
  int IdNgbhVortex3;
  int IdNgbhVortex4;
  int IdNewVortex;
  int IdSubdiviedEdge;
  int IdNewEdge1;
  int IdNewEdge2;
  int IdNewEdge3;
  int IdNgbhEdge1;
  int IdNgbhEdge2;
  int IdNgbhEdge3;
  int IdNgbhEdge4;
  int IdSubdiviedElement1;
  int IdSubdiviedElement2;
  int IdNewElement2;
  int IdNewElement1;
  
  //reallocations in LocMesh
  LocMesh->NbVertexes++;
  LocMesh->NbEdges=LocMesh->NbEdges+3;
  LocMesh->NbElements=LocMesh->NbElements+2;
  
  LocMesh->Vertexes=(Vertex *)realloc(LocMesh->Vertexes,LocMesh->NbVertexes*sizeof(Vertex));
  LocMesh->Edges=(Edge *)realloc(LocMesh->Edges,LocMesh->NbEdges*sizeof(Edge));
  LocMesh->Elements=(Element *)realloc(LocMesh->Elements,LocMesh->NbElements*sizeof(Element));
  
  //identifiers to take into account
  IdNgbhVortex1=LocMesh->Edges[IdEdge].V1;
  IdNgbhVortex2=LocMesh->Edges[IdEdge].V2;
  IdNewVortex=LocMesh->NbVertexes-1;
  
  IdSubdiviedEdge=IdEdge;
  IdNewEdge1=LocMesh->NbEdges-3;
  IdNewEdge2=LocMesh->NbEdges-2;
  IdNewEdge3=LocMesh->NbEdges-1;
  
  IdSubdiviedElement1=-1;
  IdSubdiviedElement2=-1;
  for (i=0;i<LocMesh->NbElements-2;i++) if ((LocMesh->Elements[i].E1==IdEdge)||(LocMesh->Elements[i].E2==IdEdge)||(LocMesh->Elements[i].E3==IdEdge)){
    if (IdSubdiviedElement1==-1){
      IdSubdiviedElement1=i;
      IdNewElement1=LocMesh->NbElements-2;
      if ((LocMesh->Elements[i].E1!=IdEdge)&&((LocMesh->Edges[LocMesh->Elements[i].E1].V1==IdNgbhVortex1)||(LocMesh->Edges[LocMesh->Elements[i].E1].V2==IdNgbhVortex1))){
        IdNgbhEdge1=LocMesh->Elements[i].E1;
        if (LocMesh->Edges[LocMesh->Elements[i].E1].V1==IdNgbhVortex1) IdNgbhVortex3=LocMesh->Edges[LocMesh->Elements[i].E1].V2;
        else IdNgbhVortex3=LocMesh->Edges[LocMesh->Elements[i].E1].V1;
      }
      if ((LocMesh->Elements[i].E2!=IdEdge)&&((LocMesh->Edges[LocMesh->Elements[i].E2].V1==IdNgbhVortex1)||(LocMesh->Edges[LocMesh->Elements[i].E2].V2==IdNgbhVortex1))){
        IdNgbhEdge1=LocMesh->Elements[i].E2;
        if (LocMesh->Edges[LocMesh->Elements[i].E2].V1==IdNgbhVortex1) IdNgbhVortex3=LocMesh->Edges[LocMesh->Elements[i].E2].V2;
        else IdNgbhVortex3=LocMesh->Edges[LocMesh->Elements[i].E2].V1;
      }
      if ((LocMesh->Elements[i].E3!=IdEdge)&&((LocMesh->Edges[LocMesh->Elements[i].E3].V1==IdNgbhVortex1)||(LocMesh->Edges[LocMesh->Elements[i].E3].V2==IdNgbhVortex1))){
        IdNgbhEdge1=LocMesh->Elements[i].E3;
        if (LocMesh->Edges[LocMesh->Elements[i].E3].V1==IdNgbhVortex1) IdNgbhVortex3=LocMesh->Edges[LocMesh->Elements[i].E3].V2;
        else IdNgbhVortex3=LocMesh->Edges[LocMesh->Elements[i].E3].V1;
      }
      if ((LocMesh->Elements[i].E1!=IdEdge)&&((LocMesh->Edges[LocMesh->Elements[i].E1].V1==IdNgbhVortex2)||(LocMesh->Edges[LocMesh->Elements[i].E1].V2==IdNgbhVortex2))){
        IdNgbhEdge3=LocMesh->Elements[i].E1;
      }
      if ((LocMesh->Elements[i].E2!=IdEdge)&&((LocMesh->Edges[LocMesh->Elements[i].E2].V1==IdNgbhVortex2)||(LocMesh->Edges[LocMesh->Elements[i].E2].V2==IdNgbhVortex2))){
        IdNgbhEdge3=LocMesh->Elements[i].E2;
      }
      if ((LocMesh->Elements[i].E3!=IdEdge)&&((LocMesh->Edges[LocMesh->Elements[i].E3].V1==IdNgbhVortex2)||(LocMesh->Edges[LocMesh->Elements[i].E3].V2==IdNgbhVortex2))){
        IdNgbhEdge3=LocMesh->Elements[i].E3;
      }
    }
    else{
      IdSubdiviedElement2=i;
      IdNewElement2=LocMesh->NbElements-1;
      if ((LocMesh->Elements[i].E1!=IdEdge)&&((LocMesh->Edges[LocMesh->Elements[i].E1].V1==IdNgbhVortex1)||(LocMesh->Edges[LocMesh->Elements[i].E1].V2==IdNgbhVortex1))){
        IdNgbhEdge2=LocMesh->Elements[i].E1;
        if (LocMesh->Edges[LocMesh->Elements[i].E1].V1==IdNgbhVortex1) IdNgbhVortex4=LocMesh->Edges[LocMesh->Elements[i].E1].V2;
        else IdNgbhVortex4=LocMesh->Edges[LocMesh->Elements[i].E1].V1;
      }
      if ((LocMesh->Elements[i].E2!=IdEdge)&&((LocMesh->Edges[LocMesh->Elements[i].E2].V1==IdNgbhVortex1)||(LocMesh->Edges[LocMesh->Elements[i].E2].V2==IdNgbhVortex1))){
        IdNgbhEdge2=LocMesh->Elements[i].E2;
        if (LocMesh->Edges[LocMesh->Elements[i].E2].V1==IdNgbhVortex1) IdNgbhVortex4=LocMesh->Edges[LocMesh->Elements[i].E2].V2;
        else IdNgbhVortex4=LocMesh->Edges[LocMesh->Elements[i].E2].V1;
      }
      if ((LocMesh->Elements[i].E3!=IdEdge)&&((LocMesh->Edges[LocMesh->Elements[i].E3].V1==IdNgbhVortex1)||(LocMesh->Edges[LocMesh->Elements[i].E3].V2==IdNgbhVortex1))){
        IdNgbhEdge2=LocMesh->Elements[i].E3;
        if (LocMesh->Edges[LocMesh->Elements[i].E3].V1==IdNgbhVortex1) IdNgbhVortex4=LocMesh->Edges[LocMesh->Elements[i].E3].V2;
        else IdNgbhVortex4=LocMesh->Edges[LocMesh->Elements[i].E3].V1;
      }
      if ((LocMesh->Elements[i].E1!=IdEdge)&&((LocMesh->Edges[LocMesh->Elements[i].E1].V1==IdNgbhVortex2)||(LocMesh->Edges[LocMesh->Elements[i].E1].V2==IdNgbhVortex2))){
        IdNgbhEdge4=LocMesh->Elements[i].E1;
      }
      if ((LocMesh->Elements[i].E2!=IdEdge)&&((LocMesh->Edges[LocMesh->Elements[i].E2].V1==IdNgbhVortex2)||(LocMesh->Edges[LocMesh->Elements[i].E2].V2==IdNgbhVortex2))){
        IdNgbhEdge4=LocMesh->Elements[i].E2;
      }
      if ((LocMesh->Elements[i].E3!=IdEdge)&&((LocMesh->Edges[LocMesh->Elements[i].E3].V1==IdNgbhVortex2)||(LocMesh->Edges[LocMesh->Elements[i].E3].V2==IdNgbhVortex2))){
        IdNgbhEdge4=LocMesh->Elements[i].E3;
      }
    }
  }
  IdNewElement1=LocMesh->NbElements-2;
  IdNewElement2=LocMesh->NbElements-1;
  
  //Coordinates of the new vertex
  LocMesh->Vertexes[IdNewVortex].X=(LocMesh->Vertexes[IdNgbhVortex1].X+LocMesh->Vertexes[IdNgbhVortex2].X)/2.;
  LocMesh->Vertexes[IdNewVortex].Y=(LocMesh->Vertexes[IdNgbhVortex1].Y+LocMesh->Vertexes[IdNgbhVortex2].Y)/2.;
  LocMesh->Vertexes[IdNewVortex].Z=(LocMesh->Vertexes[IdNgbhVortex1].Z+LocMesh->Vertexes[IdNgbhVortex2].Z)/2.;
  
  //update the vertex identifiers of the edge ends
  LocMesh->Edges[IdSubdiviedEdge].V1=IdNgbhVortex1;
  LocMesh->Edges[IdSubdiviedEdge].V2=IdNewVortex;
  
  LocMesh->Edges[IdNewEdge1].V1=IdNewVortex;
  LocMesh->Edges[IdNewEdge1].V2=IdNgbhVortex2;
  LocMesh->Edges[IdNewEdge1].NbSubdiv=1;
  
  LocMesh->Edges[IdNewEdge2].V1=IdNewVortex;
  LocMesh->Edges[IdNewEdge2].V2=IdNgbhVortex4;
  LocMesh->Edges[IdNewEdge2].NbSubdiv=1;

  LocMesh->Edges[IdNewEdge3].V1=IdNewVortex;
  LocMesh->Edges[IdNewEdge3].V2=IdNgbhVortex3;
  LocMesh->Edges[IdNewEdge3].NbSubdiv=1;


  //update the edges identifiers of the elements
  LocMesh->Elements[IdSubdiviedElement1].E1=IdSubdiviedEdge;
  LocMesh->Elements[IdSubdiviedElement1].E2=IdNewEdge3;
  LocMesh->Elements[IdSubdiviedElement1].E3=IdNgbhEdge1;
  
  LocMesh->Elements[IdSubdiviedElement2].E1=IdNgbhEdge2;
  LocMesh->Elements[IdSubdiviedElement2].E2=IdNewEdge2;
  LocMesh->Elements[IdSubdiviedElement2].E3=IdSubdiviedEdge;
  
  LocMesh->Elements[IdNewElement1].E1=IdNewEdge1;
  LocMesh->Elements[IdNewElement1].E2=IdNgbhEdge3;
  LocMesh->Elements[IdNewElement1].E3=IdNewEdge3;
  
  LocMesh->Elements[IdNewElement2].E1=IdNgbhEdge4;
  LocMesh->Elements[IdNewElement2].E2=IdNewEdge1;
  LocMesh->Elements[IdNewElement2].E3=IdNewEdge2;
}



//point projection on the surface of the segmented shape
void ProjectPointMeshSurf(double *PtX,double *PtY,double *PtZ,int *** ImSeg,int NBX,int NBY,int  NBZ,double MeanX,double MeanY,double MeanZ){
  double DirecX,DirecY,DirecZ,length;

  //project the coordinates
  DirecX=*PtX-MeanX;
  DirecY=*PtY-MeanY;
  DirecZ=*PtZ-MeanZ;
  length=sqrt(pow(DirecX,2.)+pow(DirecY,2.)+pow(DirecZ,2.));
  DirecX=DirecX/(length*2);
  DirecY=DirecY/(length*2);
  DirecZ=DirecZ/(length*2);
  if (ImSeg[(int)(*PtZ+0.5)][(int)(*PtY+0.5)][(int)(*PtX+0.5)]==1){//go outside of the surface
    while(ImSeg[(int)(*PtZ+0.5)][(int)(*PtY+0.5)][(int)(*PtX+0.5)]==1){
      *PtZ+=DirecZ;
      *PtY+=DirecY;
      *PtX+=DirecX;
    }
  }
  else{//go inside of the surface
    while(ImSeg[(int)(*PtZ+0.5)][(int)(*PtY+0.5)][(int)(*PtX+0.5)]==0){
      *PtZ-=DirecZ;
      *PtY-=DirecY;
      *PtX-=DirecX;
    }
    *PtZ+=DirecZ;
    *PtY+=DirecY;
    *PtX+=DirecX;
  }
}



//order the edges of all element so that this order is clockwise if we see it from outside the shape
void MakeClockwiseOrder(Mesh * LocMesh){
  int i,j,jswap;
  int Ed1,Ed2,Ed3;  //Order: Vertex 1 -> Edge1 -> Vertex 2 -> Edge 2 -> Vertex 3 -> Edge 3 -> Vertex 1 ...
  int temp;
  double x1,x2,x3,y1,y2,y3,z1,z2,z3; //coordinates
  double u1,u2,u3,v1,v2,v3;  //vectors
  double vp1,vp2,vp3;  //cross product
  double tempDbl;
  
  //test all elements of the mesh
  for (i=0;i<LocMesh->NbElements;i++){
    //order the vertexes and edges
    Ed1=LocMesh->Elements[i].E1;
    Ed2=LocMesh->Elements[i].E2;
    Ed3=LocMesh->Elements[i].E3;
    
    if ((LocMesh->Edges[Ed1].V2!=LocMesh->Edges[Ed2].V1)&&(LocMesh->Edges[Ed1].V2!=LocMesh->Edges[Ed2].V2)){ //swap edges
      temp=Ed2; Ed2=Ed3; Ed3=temp;
    }
    
    if (LocMesh->Edges[Ed1].V2!=LocMesh->Edges[Ed2].V1){ //swap vertexes
      temp=LocMesh->Edges[Ed2].V1; LocMesh->Edges[Ed2].V1=LocMesh->Edges[Ed2].V2; LocMesh->Edges[Ed2].V2=temp;
      if (LocMesh->Edges[Ed2].NbSubdiv>1)
        for (j=0;j<=LocMesh->Edges[Ed2].NbSubdiv/2;j++){
        jswap=LocMesh->Edges[Ed2].NbSubdiv-j;
        tempDbl=LocMesh->Edges[Ed2].SubiVertex[j].X; LocMesh->Edges[Ed2].SubiVertex[j].X=LocMesh->Edges[Ed2].SubiVertex[jswap].X;  LocMesh->Edges[Ed2].SubiVertex[jswap].X=tempDbl;
        tempDbl=LocMesh->Edges[Ed2].SubiVertex[j].Y; LocMesh->Edges[Ed2].SubiVertex[j].Y=LocMesh->Edges[Ed2].SubiVertex[jswap].Y;  LocMesh->Edges[Ed2].SubiVertex[jswap].Y=tempDbl;
        tempDbl=LocMesh->Edges[Ed2].SubiVertex[j].Z; LocMesh->Edges[Ed2].SubiVertex[j].Z=LocMesh->Edges[Ed2].SubiVertex[jswap].Z;  LocMesh->Edges[Ed2].SubiVertex[jswap].Z=tempDbl;
        }
    }
    
    if (LocMesh->Edges[Ed2].V2!=LocMesh->Edges[Ed3].V1){ //swap vertexes
      temp=LocMesh->Edges[Ed3].V1; LocMesh->Edges[Ed3].V1=LocMesh->Edges[Ed3].V2; LocMesh->Edges[Ed3].V2=temp;
      if (LocMesh->Edges[Ed3].NbSubdiv>1)
        for (j=0;j<=LocMesh->Edges[Ed3].NbSubdiv/2;j++){
        jswap=LocMesh->Edges[Ed3].NbSubdiv-j;
        tempDbl=LocMesh->Edges[Ed3].SubiVertex[j].X; LocMesh->Edges[Ed3].SubiVertex[j].X=LocMesh->Edges[Ed3].SubiVertex[jswap].X;  LocMesh->Edges[Ed3].SubiVertex[jswap].X=tempDbl;
        tempDbl=LocMesh->Edges[Ed3].SubiVertex[j].Y; LocMesh->Edges[Ed3].SubiVertex[j].Y=LocMesh->Edges[Ed3].SubiVertex[jswap].Y;  LocMesh->Edges[Ed3].SubiVertex[jswap].Y=tempDbl;
        tempDbl=LocMesh->Edges[Ed3].SubiVertex[j].Z; LocMesh->Edges[Ed3].SubiVertex[j].Z=LocMesh->Edges[Ed3].SubiVertex[jswap].Z;  LocMesh->Edges[Ed3].SubiVertex[jswap].Z=tempDbl;
        }
    }
    
    if (LocMesh->Edges[Ed3].V2!=LocMesh->Edges[Ed1].V1){
      printf("This is not good!!!\n");
    }
    
    //invert the order if necessary
    x1=LocMesh->Vertexes[LocMesh->Edges[Ed1].V1].X;  y1=LocMesh->Vertexes[LocMesh->Edges[Ed1].V1].Y;  z1=LocMesh->Vertexes[LocMesh->Edges[Ed1].V1].Z;
    x2=LocMesh->Vertexes[LocMesh->Edges[Ed1].V2].X;  y2=LocMesh->Vertexes[LocMesh->Edges[Ed1].V2].Y;  z2=LocMesh->Vertexes[LocMesh->Edges[Ed1].V2].Z;
    x3=LocMesh->Vertexes[LocMesh->Edges[Ed3].V1].X;  y3=LocMesh->Vertexes[LocMesh->Edges[Ed3].V1].Y;  z3=LocMesh->Vertexes[LocMesh->Edges[Ed3].V1].Z;
    
    u1=x2-x1; u2=y2-y1; u3=z2-z1;
    v1=x3-x1; v2=y3-y1; v3=z3-z1;
    
    vp1=u2*v3-u3*v2;
    vp2=u3*v1-u1*v3;
    vp3=u1*v2-u2*v1;
    
    //printf("%3.2lf %3.2lf %3.2lf | %3.2lf %3.2lf %3.2lf | %3.2lf %3.2lf %3.2lf || %3.2lf %3.2lf %3.2lf || %3.2lf %3.2lf %3.2lf | %lf\n",x1,y1,z1,x2,y2,z2,x3,y3,z3,u1,u2,u3,v1,v2,v3,x1*vp1+y1*vp2+z1*vp3);
    //printf("%3.2lf %3.2lf %3.2lf | %3.2lf %3.2lf %3.2lf | %3.2lf\n",x1,y1,z1,vp1,vp2,vp3,x1*vp1+y1*vp2+z1*vp3);
    
    if (x1*vp1+y1*vp2+z1*vp3<0){//we consider that the origin is within the spheric shape so the vector (x1,y1,z1) points out of the shape
      //swap edges
      temp=Ed3; Ed3=Ed2; Ed2=temp;
      //swap vertexes
      temp=LocMesh->Edges[Ed1].V1; LocMesh->Edges[Ed1].V1=LocMesh->Edges[Ed1].V2; LocMesh->Edges[Ed1].V2=temp;
      if (LocMesh->Edges[Ed1].NbSubdiv>1)
        for (j=0;j<=LocMesh->Edges[Ed1].NbSubdiv/2;j++){
        jswap=LocMesh->Edges[Ed1].NbSubdiv-j;
        tempDbl=LocMesh->Edges[Ed1].SubiVertex[j].X; LocMesh->Edges[Ed1].SubiVertex[j].X=LocMesh->Edges[Ed1].SubiVertex[jswap].X;  LocMesh->Edges[Ed1].SubiVertex[jswap].X=tempDbl;
        tempDbl=LocMesh->Edges[Ed1].SubiVertex[j].Y; LocMesh->Edges[Ed1].SubiVertex[j].Y=LocMesh->Edges[Ed1].SubiVertex[jswap].Y;  LocMesh->Edges[Ed1].SubiVertex[jswap].Y=tempDbl;
        tempDbl=LocMesh->Edges[Ed1].SubiVertex[j].Z; LocMesh->Edges[Ed1].SubiVertex[j].Z=LocMesh->Edges[Ed1].SubiVertex[jswap].Z;  LocMesh->Edges[Ed1].SubiVertex[jswap].Z=tempDbl;
        }
    }
    //new order of the elements
    LocMesh->Elements[i].E1=Ed1;
    LocMesh->Elements[i].E2=Ed2;
    LocMesh->Elements[i].E3=Ed3;
  }

}



//Project de coordinates of a mesh on a surface in  and refine the mesh until 
//the longest edge has a length smaller than 'MaxEdgeLength'.
//If NbSubdiv>1, each edge is curved into 'NbSubdiv' parts in the end
void ProjectMeshSurface(Mesh * LocMesh, double MaxEdgeLength,int NbSubdiv, int***ImSeg,int NBX,int NBY,int NBZ){
  int i,j,k,x,y,z;
  int * order;
  double SqLengthI,SqLengthJ;
  int OrigNbEdges,OrigNbVertexes;
  int OK;
  double dj,dNbSubdiv;
  double EdgeV1_X,EdgeV1_Y,EdgeV1_Z,EdgeV2_X,EdgeV2_Y,EdgeV2_Z;

  //PART 1: initialize
  order=(int*)malloc(LocMesh->NbEdges*sizeof(int));
  for (i=0;i<LocMesh->NbEdges;i++) order[i]=i;
  OK=1;
  
  
  
  //find the center of the shape
  double MeanX,MeanY,MeanZ;
  int NbPts;
  
  MeanX=0;   MeanY=0;   MeanZ=0;  NbPts=0;
  for (z = 0; z < NBZ; z++)  for (y = 0; y < NBY; y++) for (x = 0; x < NBX; x++) if (ImSeg[z][y][x]==1){ 
    NbPts++;
    MeanX+=(double)x;
    MeanY+=(double)y;
    MeanZ+=(double)z;
  }
  MeanX=MeanX/((double)NbPts);
  MeanY=MeanY/((double)NbPts);
  MeanZ=MeanZ/((double)NbPts);
  
  //project the center of the mesh on the center find the center of the shape
  for (i=0;i<LocMesh->NbVertexes;i++){
    LocMesh->Vertexes[i].X=LocMesh->Vertexes[i].X+MeanX-0.5;
    LocMesh->Vertexes[i].Y=LocMesh->Vertexes[i].Y+MeanY-0.5;
    LocMesh->Vertexes[i].Z=LocMesh->Vertexes[i].Z+MeanZ-0.5;
  }
  
  cout << MeanX << " " << MeanY << " " << MeanZ << "\n";
  
  //PART 2: mesh projection on the surface of the segmented shape
  for (i=0;i<LocMesh->NbVertexes;i++){
    ProjectPointMeshSurf(&(LocMesh->Vertexes[i].X),&(LocMesh->Vertexes[i].Y),&(LocMesh->Vertexes[i].Z),ImSeg,NBX,NBY,NBZ,MeanX,MeanY,MeanZ);
  }

  //Big loop to subdivide the mesh until each edge has a length smaller than 'MaxEdgeLength'
  while (OK==1){
    OK=0;
    
    order=(int*)realloc(order,LocMesh->NbEdges*sizeof(int));
    for (i=0;i<LocMesh->NbEdges;i++) order[i]=i;

    //order the edges as a function of their length
    for (i=0;i<LocMesh->NbEdges-1;i++) for (j=i+1;j<LocMesh->NbEdges;j++){
      SqLengthI=pow(LocMesh->Vertexes[LocMesh->Edges[order[i]].V1].X-LocMesh->Vertexes[LocMesh->Edges[order[i]].V2].X,2.0);
      SqLengthI+=pow(LocMesh->Vertexes[LocMesh->Edges[order[i]].V1].Y-LocMesh->Vertexes[LocMesh->Edges[order[i]].V2].Y,2.0);
      SqLengthI+=pow(LocMesh->Vertexes[LocMesh->Edges[order[i]].V1].Z-LocMesh->Vertexes[LocMesh->Edges[order[i]].V2].Z,2.0);
      SqLengthJ=pow(LocMesh->Vertexes[LocMesh->Edges[order[j]].V1].X-LocMesh->Vertexes[LocMesh->Edges[order[j]].V2].X,2.0);
      SqLengthJ+=pow(LocMesh->Vertexes[LocMesh->Edges[order[j]].V1].Y-LocMesh->Vertexes[LocMesh->Edges[order[j]].V2].Y,2.0);
      SqLengthJ+=pow(LocMesh->Vertexes[LocMesh->Edges[order[j]].V1].Z-LocMesh->Vertexes[LocMesh->Edges[order[j]].V2].Z,2.0);
      
      if (SqLengthI<SqLengthJ){ k=order[i]; order[i]=order[j]; order[j]=k;}
    }
    
    //subdivide the edges if their length is larger than 'MaxEdgeLength'
    OrigNbVertexes=LocMesh->NbVertexes;
    OrigNbEdges=LocMesh->NbEdges;
    for (i=0;i<OrigNbEdges;i++){
      SqLengthI=pow(LocMesh->Vertexes[LocMesh->Edges[order[i]].V1].X-LocMesh->Vertexes[LocMesh->Edges[order[i]].V2].X,2.0);
      SqLengthI+=pow(LocMesh->Vertexes[LocMesh->Edges[order[i]].V1].Y-LocMesh->Vertexes[LocMesh->Edges[order[i]].V2].Y,2.0);
      SqLengthI+=pow(LocMesh->Vertexes[LocMesh->Edges[order[i]].V1].Z-LocMesh->Vertexes[LocMesh->Edges[order[i]].V2].Z,2.0);
      
      if (sqrt(SqLengthI)>MaxEdgeLength){
        //subdvision
        OK=1;
        CutEdge(LocMesh,order[i]);    //rem: 'LocMesh' is changed here. In particular LocMesh->NbEdges is larger.
      }
    }
    
    //mesh projection on the surface of the segmented shape
    for (i=OrigNbVertexes;i<LocMesh->NbVertexes;i++)
      ProjectPointMeshSurf(&(LocMesh->Vertexes[i].X),&(LocMesh->Vertexes[i].Y),&(LocMesh->Vertexes[i].Z),ImSeg,NBX,NBY,NBZ,MeanX,MeanY,MeanZ);
    
    cout << "Nb Edges: " << LocMesh->NbEdges <<  " / Nb Elements: " << LocMesh->NbElements <<  "\n";
  }
  
  
  //PART 3: subdivision of the edges and projection on the segmented shape boundaries
  if (NbSubdiv>1){
    for (i=0;i<LocMesh->NbEdges;i++){
      //init
      LocMesh->Edges[i].NbSubdiv=NbSubdiv;
      LocMesh->Edges[i].SubiVertex=(Vertex *)malloc((LocMesh->Edges[i].NbSubdiv+1)*sizeof(Vertex));
      
      //coordinates of the subdivised edge
      EdgeV1_X=LocMesh->Vertexes[LocMesh->Edges[i].V1].X;
      EdgeV1_Y=LocMesh->Vertexes[LocMesh->Edges[i].V1].Y;
      EdgeV1_Z=LocMesh->Vertexes[LocMesh->Edges[i].V1].Z;
      EdgeV2_X=LocMesh->Vertexes[LocMesh->Edges[i].V2].X;
      EdgeV2_Y=LocMesh->Vertexes[LocMesh->Edges[i].V2].Y;
      EdgeV2_Z=LocMesh->Vertexes[LocMesh->Edges[i].V2].Z;
      dNbSubdiv=(double)LocMesh->Edges[i].NbSubdiv;
      
      for (j=0;j<LocMesh->Edges[i].NbSubdiv+1;j++){
        dj=(double)j;
        
        //interpolation
        LocMesh->Edges[i].SubiVertex[j].X=EdgeV1_X*(dNbSubdiv-dj)/dNbSubdiv+EdgeV2_X*dj/dNbSubdiv;
        LocMesh->Edges[i].SubiVertex[j].Y=EdgeV1_Y*(dNbSubdiv-dj)/dNbSubdiv+EdgeV2_Y*dj/dNbSubdiv;
        LocMesh->Edges[i].SubiVertex[j].Z=EdgeV1_Z*(dNbSubdiv-dj)/dNbSubdiv+EdgeV2_Z*dj/dNbSubdiv;
        
        //projection on the segmented shape boundaries
        ProjectPointMeshSurf(&(LocMesh->Edges[i].SubiVertex[j].X),&(LocMesh->Edges[i].SubiVertex[j].Y),&(LocMesh->Edges[i].SubiVertex[j].Z),ImSeg,NBX,NBY,NBZ,MeanX,MeanY,MeanZ);
      }
    }
  }
}






//write the geometry of mesh contained in a Mesh structure in a XML-Nektar file
void WriteMesh(Mesh * LocMesh,char FileName[256]){
  FILE * XmlMeshGeomFile;
  int i,j,NbCurvedEdge;
  
  //to have a proper order of the edges within the elements
  MakeClockwiseOrder(LocMesh);
  
  //open file
  XmlMeshGeomFile = fopen(FileName,"w");
  
  //write header
  fprintf(XmlMeshGeomFile,"<?xml version=\"1.0\" encoding=\"utf-8\"?>\n");
  fprintf(XmlMeshGeomFile,"\n");
  fprintf(XmlMeshGeomFile,"<NEKTAR>\n");
  fprintf(XmlMeshGeomFile,"<!-- Embed a 2-dimensional object in a 2-dimensional space -->\n");
  fprintf(XmlMeshGeomFile,"<!-- DIM <= SPACE -->\n");
  fprintf(XmlMeshGeomFile,"<!-- This provides a method of optimizing code for a 1-D curve embedded in 3-space. -->\n");
  fprintf(XmlMeshGeomFile,"<GEOMETRY DIM=\"2\" SPACE=\"3\">\n");
  fprintf(XmlMeshGeomFile,"\n");
  
  //write vertex
  fprintf(XmlMeshGeomFile,"  <VERTEX>\n");
  fprintf(XmlMeshGeomFile,"    <!-- Always must have four values per entry. -->\n");
  for (i=0;i<LocMesh->NbVertexes;i++)
    fprintf(XmlMeshGeomFile,"    <V ID=\"%d\">  %f   %f   %f  </V>\n",i,LocMesh->Vertexes[i].X,LocMesh->Vertexes[i].Y,LocMesh->Vertexes[i].Z);
  fprintf(XmlMeshGeomFile,"  </VERTEX>\n");
  fprintf(XmlMeshGeomFile,"  \n");
  
  //write edges
  fprintf(XmlMeshGeomFile,"  <EDGE>\n");
  fprintf(XmlMeshGeomFile,"    <!--Edges are vertex pairs -->\n");
  for (i=0;i<LocMesh->NbEdges;i++)
    fprintf(XmlMeshGeomFile,"    <E ID=\"%d\">    %d  %d   </E>\n",i,LocMesh->Edges[i].V1,LocMesh->Edges[i].V2);
  fprintf(XmlMeshGeomFile,"  </EDGE>\n");
  fprintf(XmlMeshGeomFile,"  \n");
  
  //write Element
  fprintf(XmlMeshGeomFile,"  <ELEMENT>\n");
  for (i=0;i<LocMesh->NbElements;i++)
    fprintf(XmlMeshGeomFile,"    <T ID=\"%d\">    %d     %d     %d </T>\n",i,LocMesh->Elements[i].E1,LocMesh->Elements[i].E2,LocMesh->Elements[i].E3);
  fprintf(XmlMeshGeomFile,"  </ELEMENT>\n");
  fprintf(XmlMeshGeomFile,"  \n");
  
  
  
  //write curved elements
  fprintf(XmlMeshGeomFile,"  <CURVED>\n");
  NbCurvedEdge=0;
  for (i=0;i<LocMesh->NbEdges;i++) if (LocMesh->Edges[i].NbSubdiv>1){
    fprintf(XmlMeshGeomFile,"    <E ID=\"%d\" EDGEID=\"%d\" TYPE=\"PolyEvenlySpaced\" NUMPOINTS=\"%d\">\n",NbCurvedEdge,i,LocMesh->Edges[i].NbSubdiv+1);
    for (j=0;j<LocMesh->Edges[i].NbSubdiv+1;j++){
      fprintf(XmlMeshGeomFile,"    %lf %lf %lf\n",LocMesh->Edges[i].SubiVertex[j].X,LocMesh->Edges[i].SubiVertex[j].Y,LocMesh->Edges[i].SubiVertex[j].Z);
    }
    fprintf(XmlMeshGeomFile,"   </E>\n\n");
    NbCurvedEdge++;
  }
  fprintf(XmlMeshGeomFile,"   </CURVED>\n");
  
  
  //write the end (could be more evolved)
  fprintf(XmlMeshGeomFile,"<!-- V - vertex, E - edge, F - face, L - element -->\n");
  fprintf(XmlMeshGeomFile,"  <COMPOSITE>\n");
  fprintf(XmlMeshGeomFile,"    <C ID=\"0\"> T[0-%d]\n",LocMesh->NbElements-1);
  fprintf(XmlMeshGeomFile,"    </C>\n");
  fprintf(XmlMeshGeomFile,"    <C ID=\"1\"> E[2,3,4,5,7,8]\n");
  fprintf(XmlMeshGeomFile,"    </C>\n");
  fprintf(XmlMeshGeomFile,"  </COMPOSITE>\n");
  fprintf(XmlMeshGeomFile,"  \n");
  fprintf(XmlMeshGeomFile,"  <DOMAIN> C[0] </DOMAIN>\n");
  fprintf(XmlMeshGeomFile,"  \n");
  
  //EOF
  fprintf(XmlMeshGeomFile,"</GEOMETRY>\n");
  
  //close file
  fclose(XmlMeshGeomFile);
  
  cout << "Mesh written in " << FileName << "\n";
  
}





// ---------------------------------------------------------------------------------------
//                                       main function
// ---------------------------------------------------------------------------------------


template <class VoxelType> void anisoDiffusion<VoxelType>::Run_3D_Explicit(){
  int i, j, x, y, z;
  int*** imageE;
  int NBX,NBY,NBZ,NBT;
  char output_name[] = "OutputMesh.xml";
  double Prec=11;
  int NbSubdiv=2;
  Mesh LocMesh;
  int GreyLevelToExtract;

  //1) initialisation
  this->Initialize();
  GreyLevelToExtract=255;
  NBX=this->_input->GetX();
  NBY=this->_input->GetY();
  NBZ=this->_input->GetZ();
  NBT=this->_input->GetT();
  
  //binarization of the input image
  imageE= (int***) malloc (NBZ*sizeof(int**));
  for (i=0;i<NBZ;i++) imageE[i]= (int**) malloc (NBY*sizeof(int*));
  for (i=0;i<NBZ;i++) for (j=0;j<NBY;j++) imageE[i][j]= (int*) malloc (NBX*sizeof(int));
	
  for (z = 0; z < NBZ; z++)  for (y = 0; y < NBY; y++) for (x = 0; x < NBX; x++){
    imageE[z][y][x]=static_cast<int>(this->_input->Get(x, y, z, 0)+0.001);
    if (imageE[z][y][x]==GreyLevelToExtract) imageE[z][y][x]=1;
    else imageE[z][y][x]=0;
    this->_output->Put(x, y, z,0, static_cast<VoxelType>(imageE[z][y][x]));
  }
  
  

  //2) run functions
  LocMesh=CreateBasicCubicMesh();
  ProjectMeshSurface(&LocMesh, Prec,NbSubdiv, imageE,NBX,NBY,NBZ);
  WriteMesh(&LocMesh,output_name);
  
  this->Finalize();
}






///+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
///+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
///+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


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
