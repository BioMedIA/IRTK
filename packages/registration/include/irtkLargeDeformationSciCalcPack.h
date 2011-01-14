
#ifndef _IRTKLARGEDEFORMATIONSCICALCPACK_H
#define _IRTKLARGEDEFORMATIONSCICALCPACK_H

#include <irtkImage.h>

///++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
///                           1:   FUNCTIONS FOR THE CLASS "ScalarField"
///++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


class ScalarField{
	/// ******************************************************************************
private:
	//size of the fields
	int NXtY;      //NX*NY
	int NXtYtZ;    //NX*NY*NZ
	
	//scalar field
	float * ScalField;
	
	/// ******************************************************************************
public:
	//size of the fields
	int NX;        //image size on X axis
	int NY;        //image size on Y axis
	int NZ;        //image size on Z axis
	int NT;        //image size on time
	
	/// Constructor and destructor
	ScalarField();
	~ScalarField();
	
	/// functions associated to the field
	//put a value in the scalar field
	virtual void P(float value,int x,int y,int z=0,int t=0);
	
	//Add a value to the scalar field
	virtual void Add(float value,int x,int y, int z=0,int t=0);
	
	// put a the same value at every points of the scalar field
	virtual void PutToAllVoxels(float cste,int t=0);
	
	//get a value from the scalar field
	virtual float G(int x,int y,int z=0,int t=0);
	
	//get a value from the scalar field by linear interpolation
	virtual float G(double x,double y,double z=0.,int t=0);
	virtual float G(float x, float y, float z=0., int t=0);
	
	// get the maximum absolute values out of the scalar field
	virtual float GetMaxAbsVal(int t=0);
	
	//read a scalar field (in a nifti image)
	virtual void Read(char *);
	
	//read a scalar field and perform linear interpolation to give it a specific size
	virtual void Read_and_Interpolate(char *,int,int,int);
	
	//create a void scalar field. All the values are initialize to 'cste' which is null by default
	virtual void CreateVoidField(int NBX,int NBY,int NBZ=1,int NBT=1,float cste=0.0);
	
	//write a scalar field in a nifti image
	virtual void Write(char *,char *);
	
	//write a scalar field in a vectorized ascii file
	virtual void WriteInAscii(char *);
	
	//return the number of voxels in a ScalarField
	virtual int GetNbVoxels(void);
};



///++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
///                           2:   FUNCTIONS FOR THE CLASS "VectorField"
///++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


class VectorField{
	/// ******************************************************************************
private:
	//size of the fields
	int NXtY;      //NX*NY
	int NXtYtZ;    //NX*NY*NZ
	int NXtYtZtT;  //NX*NY*NZ*NT
	
	//scalar field
	float * VecField;
	
	/// ******************************************************************************
public:
	//size of the fields
	int NX;        //image size on X axis
	int NY;        //image size on Y axis
	int NZ;        //image size on Z axis
	int NT;        //image size on time
	
	/// Constructor
	VectorField();
	
	/// Destructor
	~VectorField();
	
	//put a value in the scalar field
	virtual void P(float value,int IdDirec,int x,int y,int z=0,int t=0);
	//Add a value to the scalar field
	virtual void Add(float value,int IdDirec,int x,int y, int z=0,int t=0);
	
	//put the same value at all entries of the vector field
	virtual void PutToAllVoxels(float cste,int t=0);
	
	//get a value from the scalar field
	virtual float G(int IdDirec,int x,int y,int z=0,int t=0);
	
	//get a value from the scalar field by linear interpolation
	virtual float G(int IdDirec,double x,double y,double z=0.,int t=0);
	virtual float G(int IdDirec,float x,float y,float z=0.,int t=0);
	
	// get the maximum of the absolute values of the vector field
	virtual float GetMaxAbsVal(int t=0);
	
	//read a vector field (in 3 nifti images -> X, Y, Z)
	virtual void Read(char *,char *,char *);
	
	//read a vector field and perform linear interpolation to give it a specific size
	//* If rescaleVF!=0, the values of the vector field are rescaled proportionally to
	//  the re-sizing (usefull for velocity fields)
	virtual void Read_and_Interpolate(char *,char *,char *,int,int,int,int rescaleVF=0);
	
	//create a void scalar field
	virtual void CreateVoidField(int NBX,int NBY,int NBZ=1,int NBT=1);
	
	//write a vector field (from 3 nifti images -> X, Y, Z)
	virtual void Write(char *,char *,char *);
};

///++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
///           3: CLASS TO PERFORM CONVOLUTION AND DECONVOLUTION USING FFT
///++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

class FFTconvolver3D{
	/// ******************************************************************************
private:
    //size of the inputs (that will be copied in images having sizes = 2^{...})
    int NX;        //image size on X axis
    int NY;        //image size on Y axis
    int NZ;        //image size on Z axis
    
    //size of the fields transformed by fft
    int NXfft;
    int NYfft;
    int NZfft;
    
    //fields transformed by fft
    ScalarField RealSignalForFFT;
    ScalarField ImagSignalForFFT;
    ScalarField RealFilterForFFT;
    ScalarField ImagFilterForFFT;
    
    //temporary scalar field
    ScalarField ImageTemp;
    
    //design a kernel that is the sum of up to 4 Gaussians   (DEPRECATED -> USING MakeSumOf7AnisotropicGaussianFilters INSTEAD)
    void MakeSumOf4AnisotropicGaussianFilters(float weight1,float sigmaX1,float sigmaY1,float sigmaZ1,float weight2,float sigmaX2,float sigmaY2,float sigmaZ2,float weight3,float sigmaX3,float sigmaY3,float sigmaZ3,float weight4,float sigmaX4,float sigmaY4,float sigmaZ4);
    
    //design a kernel that is the sum of up to 7 Gaussians
    void MakeSumOf7AnisotropicGaussianFilters(float weight1,float sigmaX1,float sigmaY1,float sigmaZ1,float weight2,float sigmaX2,float sigmaY2,float sigmaZ2,float weight3,float sigmaX3,float sigmaY3,float sigmaZ3,float weight4,float sigmaX4,float sigmaY4,float sigmaZ4,float weight5,float sigmaX5,float sigmaY5,float sigmaZ5,float weight6,float sigmaX6,float sigmaY6,float sigmaZ6,float weight7,float sigmaX7,float sigmaY7,float sigmaZ7);
    
    //Fast Fourier Transform
    void DirectFFT(ScalarField * RealSignal,ScalarField * ImaginarySignal);
    
    //Inverse Fast Fourier Transform
    void InverseFFT(ScalarField * RealSignal,ScalarField * ImaginarySignal);
    
    //Fast Fourier Transform of numerical recipies (slighly modified)
    void four1NR(float data[], unsigned long nn, int isign);
    
	/// ******************************************************************************
public:
    //Constructor
    FFTconvolver3D();
	
    //Destructor
    ~FFTconvolver3D();
    
    //Initiate the complex fields for the FFT and the smoothing kernels being the sum of up to 
    //4 Gaussians (set some weights to 0 if less Gaussians are required)
    //* NX, NY, NZ: is the size of the input image
    //* w1,sX1,sY1,sZ1,: weight of the 1st Gaussian kernel and std. dev. in direction X, Y, Z
    //* w2,sX2,sY2,sZ2,: weight of the 2nd Gaussian kernel and std. dev. in direction X, Y, Z
    //* w3,sX3,sY3,sZ3,: weight of the 3rd Gaussian kernel and std. dev. in direction X, Y, Z
    //* w4,sX4,sY4,sZ4,: weight of the 4th Gaussian kernel and std. dev. in direction X, Y, Z
    //* w5,sX5,sY5,sZ5,: weight of the 5th Gaussian kernel and std. dev. in direction X, Y, Z
    //* w6,sX6,sY6,sZ6,: weight of the 6th Gaussian kernel and std. dev. in direction X, Y, Z
    //* w7,sX7,sY7,sZ7,: weight of the 7th Gaussian kernel and std. dev. in direction X, Y, Z
    virtual void InitiateConvolver(int NBX,int NBY, int NBZ, float w1=1.,float sX1=1.,float sY1=1.,float sZ1=1., float w2=0.,float sX2=1.,float sY2=1.,float sZ2=1., float w3=0.,float sX3=1.,float sY3=1.,float sZ3=1., float w4=0.,float sX4=1.,float sY4=1.,float sZ4=1., float w5=0.,float sX5=1.,float sY5=1.,float sZ5=1., float w6=0.,float sX6=1.,float sY6=1.,float sZ6=1., float w7=0.,float sX7=1.,float sY7=1.,float sZ7=1.);
    
    //change the kernel of the convolver (same notations as the constructor)
    virtual void ChangeKernel(float w1=1.,float sX1=1.,float sY1=1.,float sZ1=1., float w2=0.,float sX2=1.,float sY2=1.,float sZ2=1., float w3=0.,float sX3=1.,float sY3=1.,float sZ3=1., float w4=0.,float sX4=1.,float sY4=1.,float sZ4=1.);
	
    //convolution of a 3D scalar field using the predifined kernel
    virtual void Convolution(ScalarField *);
	
    //convolution of the real scalar field defined inside of the class using the predifined kernel
    virtual void Convolution();
	
    //deconvolution of a 3D scalar field using the predifined kernel
    // !!! NOT VALIDATED !!!
    virtual void Deconvolution(ScalarField *);
    
    //put a value in the real part of the field that is transformed by the class
    virtual void P(float value,int x,int y, int z=0);
	
    //put a value in the real part of the field that is transformed by the class
    virtual float G(int x,int y, int z=0);
};


///++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
///           4: LOW LEVEL FUNCTIONS MAKING USE OF THE CLASSES ScalarField AND VectorField 
///++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

//Compute the gradient of the scalar field "SField" and put the result in "Gradient"
//* 'DeltaX' is the spatial step between two voxels.
//* If SpecificTimeFrame<0, the calculations are done in all time frames. They are only done
//  in the time frame 'SpecificTimeFrame' otherwise.
void Cpt_Grad_ScalarField(ScalarField * SField,VectorField * Gradient,int SpecificTimeFrame=-1,float DeltaX=1);

//Compute (d VField(X) / d x) + (d VField(Y) / d y) + (d VField(Z) / d z) and put the result in 'GradScalVF'
//where 'VField' is a vector field and 'GradScalVF' a scalar field
//* 'DeltaX' is the spatial step between two voxels.
//* If SpecificTimeFrame<0, the calculations are done in all time frames. They are only done
//  in the time frame 'SpecificTimeFrame' otherwise.
void Cpt_Grad_Scal_VectorField(VectorField * VField,ScalarField * GradScalVF,int SpecificTimeFrame=-1,float DeltaX=1);

//Compute the determinant of the Jacobian of the vector field 'VField' and put the result in the scalar field 'DetJ'
//* 'DeltaX' is the spatial step between two voxels.
//* If SpecificTimeFrame<0, the calculations are done in all time frames. They are only done
//  in the time frame 'SpecificTimeFrame' otherwise.
void Cpt_JacobianDeterminant(VectorField * VField,ScalarField * DetJ,int SpecificTimeFrame=-1,float DeltaX=1);


//Compute the 3D+t mapping 'Map' from the time step 'refTimeStep' by following the velocity field 'VeloField'.
//'MappingAtRefTimeStep' is the mapping at refTimeStep (possibly not the identity).
//* An iterative leap-frog like technique is performed to compute the backward mapping. The more
//  iterations (='ConvergenceSteps') the more accurate the mapping. For most of the deformations
//  1 iteration is far enough but more iterations are suitable if the deformations have large 
//  Jacobians.
//* 'DeltaX' is the spatial step between two voxels.
void CptMappingFromVeloField(int refTimeStep,VectorField * MappingAtRefTimeStep,VectorField * VeloField,VectorField * Map,int ConvergenceSteps=1,float DeltaX=1);


//We consider here that 'PartialVeloField' contributes to 'VeloField'   (VeloField= [A velocity field] + PartialVeloField).
//This function then computes 'PartialMapping' which is the partial mapping of 'MappingAtRefTimeStep' from the time 
//subdivision 'refTimeStep' due to the contribution of 'PartialVeloField'. Note, that an Identity mapping 'MappingId' is //also defined in the inputs (to avoid defining it each time the function is used)
void CptPartialMappingFromVeloFields(int refTimeStep,VectorField * MappingAtRefTimeStep,VectorField * MappingId,VectorField * VeloField,VectorField * PartialVeloField,VectorField * PartialMapping,int ConvergenceSteps=1,float DeltaX=1);


//We consider here that 'PartialVeloField' contributes to 'VeloField'   (VeloField= [A velocity field] + PartialVeloField).
//This function then computes 'PartialLocMap' which is the partial mapping ONLY AT 'TargetSubdiv' FROM 'SourceSubdiv' due to the contribution of PartialVeloField.
//-> PartialLocMap therefore represents where are the coordinates of the points of time subdivision 'SourceSubdiv' when transported on time subdivision 'TargetSubdiv'
void ComputeLagrangianPartialMapping(int SourceSubdiv,int TargetSubdiv,VectorField * VeloField,VectorField * PartialVeloField,VectorField * PartialLocBmap,float DeltaX=1);


//Compute the projection of a 3D image 'ImagToPropag' using Mapping 'Map'.
//The image is projected at the time step 'TimeStepProj' of 'Map' and stored in 'ImageTimeT'.
//
//Importantly, the Mapping 'Map' should be an identity transformation at the time step 't' where 'ImagToPropag' is.
//It should also represent a forward mapping after 't' and a backward mapping before 't'.
void Project3Dimage(ScalarField * ImagToPropag,VectorField * Map,ScalarField * ImageTimeT,int TimeStepProj);


//By following the flow defined by the velocity field 'VeloField4Flow' measure the contribution of
//'VeloField4Measure' in the length of the flow from each point of the field. The length of flow
//is projected AT T=0 and returned in the 3D scalar field 'LengthOfFlow'
// * 'VeloField4Measure' is assumed to be part of a linear decomposition of 'VeloField4Flow'.
// * If 'VeloField4Measure'=='VeloField4Flow' then the length of the flow defined by 'VeloField4Flow'
//   is computed.
void CptLengthOfFlow(VectorField * VeloField4Flow,VectorField * VeloField4Measure,ScalarField * LengthOfFlow,int ConvergenceSteps=3,float DeltaX=1);

//By following the flow defined by the velocity field 'VeloField4Flow' measure the contribution of
//'VeloField4Measure' AT THE CURRENT TIME in the length of the flow from each point of the field. The length of flow
//is returned in the 3D+t scalar field 'LengthOfFlow'
// * 'VeloField4Measure' is assumed to be part of a linear decomposition of 'VeloField4Flow'.
// * If 'VeloField4Measure'=='VeloField4Flow' then the length of the flow defined by 'VeloField4Flow'
//   is computed.
void CptEvoLengthOfFlow(VectorField * VeloField4Flow,VectorField * VeloField4Measure,ScalarField * LengthOfFlow,int ConvergenceSteps=3,float DeltaX=1);


//compute the L_2 norm of the difference between two scalar fields
float CalcSqrtSumOfSquaredDif(ScalarField * I1,ScalarField * I2);

// Computes the transport of the initial momentum by the diffeo and stores it in Momentum
void TransportMomentum(ScalarField *InitialMomentum, VectorField *InvDiffeo, ScalarField *Momentum,float DeltaX,int t=0);

// Computes cste * transport momentum from the initial momentum by the diffeo and add it in Image.
void AddTransportMomentum(ScalarField *InitialMomentum,VectorField *TempInvDiffeo, ScalarField *Momentum,float DeltaX,float cste=1.0, int t=0);

// Computes the transport image from the initial image by the diffeo and stores it in Image.
void TransportImage(ScalarField *InitialImage, VectorField *TempInvDiffeo, ScalarField *Image,int t=0);

// Computes cste * transport image from the initial image by the diffeo and add it in Image.
void AddTransportImage(ScalarField *InitialImage, VectorField *TempInvDiffeo, ScalarField *Image,float cste=1.0,int t=0);

// Copies the values of a VectorField1(t=0) in VectorField2(t)
void DeepCopy(VectorField *VectorField1,VectorField *VectorField2,int t);

// Copies the values of a ScalarField1(t=0) in ScalarField2(t)
void DeepCopy(ScalarField *ScalarField1,ScalarField *ScalarField2,int t);

// Compute the L^2 scalar product and store it in ScalarField0
void ScalarProduct(VectorField *VectorField1, VectorField *VectorField2, ScalarField *ScalarField0, int t=0,float cste = 1.0);
void ScalarProduct(ScalarField *ScalarField1, ScalarField *ScalarField2, ScalarField *ScalarField0, int t=0, float cste = 1.0);

// Compute the L^2 scalar product between two vectorfields at time t and add it to ScalarField0
void AddScalarProduct(VectorField *VectorField1, VectorField *VectorField2, ScalarField *ScalarField0, int t=0);
void AddScalarProduct(ScalarField *ScalarField1, ScalarField *ScalarField2, ScalarField *ScalarField0, int t=0);

// Add  cste * ScalarField1 at time t1 to ScalarField2 at time t2
void AddScalarField(ScalarField *ScalarField1, ScalarField *ScalarField2,float cste,int t1 = 0,int t2=0);
// Add  cste * VectorField1 at time t1 to VectorField2 at time t2
void AddVectorField(VectorField *VectorField1, VectorField *VectorField2,float cste,int t1 = 0,int t2=0);
// Multiply a vector field by the cste
void MultiplyVectorField(VectorField *VectorField1, float cste,int t=0);
// 
void SumVectorField(VectorField *VectorField1, VectorField *VectorField2, VectorField *Output, int t1=0,int t2=0,int t3=0, float cste1 = 1.0,float cste2 =1.0);
// Compute the product element by element of ScalarField and VectorField and store it in VectorField2
void Product(ScalarField *ScalarField, VectorField *VectorField1, VectorField *VectorField2);

// Compute the dot product
float DotProduct(ScalarField *ScalarField1, ScalarField *ScalarField2,int t1=0,int t2=0);
#endif
