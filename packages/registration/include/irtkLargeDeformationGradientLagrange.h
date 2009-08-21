/*=========================================================================

  Date      : $Date: 29.06.2009$
  Changes   : $Authors: Laurent Risser, Francois-Xavier Vialard$

=========================================================================*/

#ifndef _IRTKLARGEDEFORMATIONGRADIENTLAGRANGE_H
#define _IRTKLARGEDEFORMATIONGRADIENTLAGRANGE_H

#include <irtkImageToImage.h>

#include <irtkImage.h>
#include <irtkGaussianBlurring.h>
#include <irtkConvolution.h>
#include <irtkScalarFunctionToImage.h>
#include <irtkGaussianBlurringWithPadding.h>

/**
 * Class for Large Deformation registration using Beg 05's technique
 */

template <class VoxelType> class LargeDefGradLagrange : public irtkImageToImage<VoxelType>
{
  private:
  float* JacobianMatrix;
  
  protected:
  /// Returns whether the filter requires buffering
  virtual Bool RequiresBuffering();

  /// Returns the name of the class
  virtual const char *NameOfClass();
  
  /// Function to launch the default calculations
  virtual void Run_Default();
  
  /// Subfunctions to perform the registration  (level 1)
  virtual void AllocateAllVariables();
  virtual void InitiateGradientDescent(int);
  virtual void ComputeDirectMapping();
  virtual void ComputeInverseMapping();
  virtual void ComputeJ0(int);
  virtual void ComputeJ1(int);
  virtual void ComputeGradientJ0();
  virtual void ComputeJacobianDeterminant(int);
  virtual void ComputeEnergyGradient(int);
  virtual void UpdateVelocityField();
  virtual void SpeedReparametrization();
  virtual void SaveResultGradientDescent(int);
  virtual void SaveVelocityFields(char *);
  virtual void LoadVelocityFields(char *);

  
  /// Subfunctions to perform the registration  (level 2)
  
  virtual void GetVectorFromVelocityField(int, float, float, float,float[3]);
  virtual void GetCoordFromDirectMapping(int, float, float, float,float[3]);
  virtual void GetCoordFromInverseMapping(int, float, float, float,float[3]);
  virtual float GetGreyLevelFromTemplate(float , float , float);
  virtual float GetGreyLevelFromTarget(float , float , float);
  virtual float GetGreyLevelFromFloatGenericImage(irtkGenericImage<float>*,float,float,float,float);

  
  /// Protected parameters
  //scalar and vector fields
  float* ImTemplate;
  float* ImTarget;
  float** VelocityField;
  float** DirectMapping;
  float** InverseMapping;
  float* J0;
  float* J1;
  float* GradJ0;
  float* DetJacobians;
  float** GradE;
  irtkGenericImage<float> Image3DTemp;

  //size of the fields
  int NX;
  int NY;
  int NZ;
  int NT;
  
  //temporary variables
  int NXtY;      //NX*NY
  int NXtYtZ;    //NX*NY*NZ
  int NXt3;      //NX*3
  int NXtYt3;    //NX*NY*3
  int NXtYtZt3;  //NX*NY*NZ*3
  
  
  ///Inline functions
  //returns where is the point (x,y,z) in a linearized 3D scalar field
  inline int ptSF(int x,int y, int z){
    return this->NXtY*z+this->NX*y+x;
  } 
  
  //returns where is the value of the direction IdDirec (0=x, 1=y, 2=z) of the point (x,y,z) in a linearized 3D vector field
  inline int ptVF(int x,int y, int z, int IdDirec){
    return this->NXtYt3*z+this->NXt3*y+x*3+IdDirec;
  }
  
  //returns...
  inline float CptDetJac3d(float Jac[9])
  {
    return Jac[0]*Jac[4]*Jac[8] + Jac[3]*Jac[7]*Jac[2] +Jac[6]*Jac[1]*Jac[5] -Jac[2]*Jac[4]*Jac[6] - Jac[5]*Jac[7]*Jac[0]- Jac[8]*Jac[1]*Jac[3];
  }
  
  
  
  public:

  /// Constructor
  LargeDefGradLagrange();
  
  /// Destructor
  ~LargeDefGradLagrange();

  /// Run  Large Deformation registration
  virtual void Run();
  
  
  /// public Parameters
  int iteration_nb;      //number of iterations 
  float epsilon;         //Threshold on the energy gradient convergence
  int NbTimeSubdiv;      //Number of subdivision of virtual time steps
  float MaxVelocityUpdate; //Maximum allowed velocity update at each subdivision and iteration (in voxels)
  irtkGreyImage target_image;  //explicit name (hopefully)
  float DeltaTimeSubdiv; //time step between two subdivision
  float sigma;           // width of the Gaussian kernel
  float alpha; // Weight between the norm and the penalty term
  int reparametrization; // reparameterisation frequency (in iterations)
  float DeltaVox;   //considered distance between two voxels (in order to evaluate the gradient on the grey levels)
  float UNDETERMINED_VALUE;   //value given to the voxels for which the value cannot be determined
  char PrefixInputVF[256];   //Prefix of the files containing an initial velocity field
  char PrefixOutputVF[256];   //Prefix of the files containing the final velocity field

};


#endif
