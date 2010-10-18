/*=========================================================================

  Date      : $Date: 03.03.2009$
  Changes   : $Author: Laurent Risser $

=========================================================================*/

#ifndef _ANISODIFFUSION_H
#define _ANISODIFFUSION_H

#include <irtkImageToImage.h>

/**
 * Class for anisotopic diffusion filtering
 */

template <class VoxelType> class anisoDiffusion : public irtkImageToImage<VoxelType>
{
protected:
  /// Returns whether the filter requires buffering
  virtual bool RequiresBuffering();

  /// Returns the name of the class
  virtual const char *NameOfClass();

public:

  /// Constructor
  anisoDiffusion();
  
  /// Destructor
  ~anisoDiffusion();

  /// Run anisotropic diffusion filtering
  virtual void Run();
  
  virtual void Run_4D_semiImplicit();
  virtual void Run_3D_semiImplicit();
  virtual void Run_4D_Explicit();
  virtual void Run_3D_Explicit();
  
  /// Parameters
  float ax; //caracteristic gray level variations considered significant between to neighbors voxels along direction x
  float ay; //caracteristic gray level variations considered significant between to neighbors voxels along direction y
  float az; //caracteristic gray level variations considered significant between to neighbors voxels along direction z
  float at; //caracteristic gray level variations considered significant between to neighbors voxels along time
  float dx; //distance between to neighbors voxels along direction x
  float dy; //distance between to neighbors voxels along direction y
  float dz; //distance between to neighbors voxels along direction z
  float dt; //distance between to neighbors voxels along time
  float dTau;         //virtual time step (for the diffusion pde)
  int ITERATIONS_NB; //number of virtual time iterations 
  bool TimeDependent;  //1 -> anisotrop filtering along the time / 0 -> otherwise
  bool SemiImplicit; //1 -> The scheme will be semi implicit (ADI) / 0 -> explicit scheme
};


void TridiagonalSolveFloat(const float *, const float *, float *, float *, float *, int);

#endif
