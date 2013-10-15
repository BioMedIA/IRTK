/*=========================================================================

  Library   : Image Registration Toolkit (IRTK)
  Module    : $Id$
  Copyright : Imperial College, Department of Computing
              Visual Information Processing (VIP), 2011 onwards
  Date      : $Date$
  Version   : $Revision$
  Changes   : $Author$

=========================================================================*/

#ifndef _irtkReconstructionb0_H

#define _irtkReconstructionb0_H

#include <irtkImage.h>
#include <irtkTransformation.h>
#include <irtkGaussianBlurring.h>
#include <irtkReconstruction.h>

#include <vector>
using namespace std;


/*

Reconstruction of volume from 2D slices

*/

class irtkReconstructionb0 : public irtkReconstruction
{

protected:
   //indicates into which group a stack belongs
  vector<int> _groups;
  vector<bool> _swap;
  vector<int> _stack_group;
  vector<irtkRealImage> _simulated;
  vector<irtkAffineTransformation> _shim;
  irtkMultiLevelFreeFormTransformation _fieldMap;
  double _fieldMapSpacing;
  double _smoothnessPenalty;

  
public:
  void StackRegistrations(vector<irtkRealImage>& stacks, vector<irtkRigidTransformation>& stack_transformations);
  void SetT2Template(irtkRealImage T2);
  irtkRealImage AdjustOrientation(irtkRealImage &image, bool swap);
  irtkAffineTransformation AdjustOrientationTransformation(irtkRealImage &image, bool swap);
  void ShimDistortion(irtkRealImage &acquired, irtkRealImage &simulated, irtkAffineTransformation &shim, bool swap);
  void Shim(vector<irtkRealImage> &stacks, int iter = 0);
  void FieldMapDistortion(irtkRealImage &acquired, irtkRealImage &simulated, irtkMultiLevelFreeFormTransformation &dist, bool swap);
  void FieldMap(vector<irtkRealImage> &stacks, int iter = 0);
  irtkRealImage Create4DImage(vector<irtkRealImage> &stacks);
  irtkRealImage AlignT2Template(irtkRealImage T2, double sigma=0);
  void CreateSimulated(vector<irtkRealImage> &stacks);
  void WriteSimulated();
  void SaveDistortionTransformations();
  void CorrectStacks(vector<irtkRealImage> &stacks);


  inline void SetGroups(vector<int>& stack_group, vector<int>& groups, vector<bool>& swap); 
  inline void SetReconstructed(irtkRealImage image);
  inline void SetFieldMapSpacing(double spacing);
  inline void SetSmoothnessPenalty(double penalty);

};

inline void irtkReconstructionb0::SetReconstructed(irtkRealImage image)
{
  _reconstructed = image;    
}

inline void irtkReconstructionb0::SetGroups(vector<int>& stack_group,vector<int>& groups, vector<bool>& swap)
{
  if(swap.size()!=groups.size())
  {
    cerr<<"Number of phase encoding directions does not match number of groups!"<<endl;
    exit(1);
  }
  
  _stack_group = stack_group;  
  _groups = groups;
  _swap = swap;
}

inline void irtkReconstructionb0::SetFieldMapSpacing(double spacing)
{
  _fieldMapSpacing = spacing;
}

inline void irtkReconstructionb0::SetSmoothnessPenalty(double penalty)
{
  _smoothnessPenalty = penalty;
}

#endif
