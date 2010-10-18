/*=========================================================================

  Library   : Image Registration Toolkit (IRTK)
  Module    : $Id$
  Copyright : Imperial College, Department of Computing
              Visual Information Processing (VIP), 2008 onwards
  Date      : $Date$
  Version   : $Revision$
  Changes   : $Author$

=========================================================================*/

#include <irtkImage.h>

#include <irtkGradientImage.h>

template <class VoxelType> bool irtkGradientImageZ<VoxelType>::RequiresBuffering(void)
{
  return true;
}

template <class VoxelType> const char *irtkGradientImageZ<VoxelType>::NameOfClass()
{
  return "irtkGradientImageZ";
}

template <class VoxelType> double irtkGradientImageZ<VoxelType>::Run(int x, int y, int z, int t)
{

  double previous = this->_input->Get(x, y, z-1, t);
  double next     = this->_input->Get(x, y, z+1, t);

  double gradient = previous - next;

  return gradient;
}


template <class VoxelType> void irtkGradientImageZ<VoxelType>::Run()
{
  int x, y, z, t;

  // Do the initial set up
  this->Initialize();

  // Check image dimensions....
  if (this->_input->GetZ() < 2) {
    cerr<<" irtkGradientImageZ: Dimensions of input image are wrong"<<endl;
    exit(1);
  }

  for ( t = 0; t < this->_input->GetT(); ++t) {
    for ( x = 0; x < this->_input->GetX(); ++x) {
      for ( y = 0; y < this->_input->GetY(); ++y) {
        this->_output->Put(x, y, 0, t, 0);
        for ( z = 1; z < this->_input->GetZ()-1; ++z) {
          this->_output->PutAsDouble(x, y, z, t, this->Run(x, y, z, t));
        }
        this->_output->Put(x, y, z, t, 0);
      }
    }
  }

  // Do the final cleaning up
  this->Finalize();
}


template class  irtkGradientImageZ<irtkBytePixel>;
template class  irtkGradientImageZ<irtkGreyPixel>;
template class  irtkGradientImageZ<irtkRealPixel>;
