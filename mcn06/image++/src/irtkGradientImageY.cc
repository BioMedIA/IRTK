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

template <class VoxelType> Bool irtkGradientImageY<VoxelType>::RequiresBuffering(void)
{
  return True;
}

template <class VoxelType> const char *irtkGradientImageY<VoxelType>::NameOfClass()
{
  return "irtkGradientImageY";
}

template <class VoxelType> double irtkGradientImageY<VoxelType>::Run(int x, int y, int z, int t)
{

  double previous = this->_input->Get(x, y-1, z, t);
  double next     = this->_input->Get(x, y+1, z, t);

  double gradient = previous - next;

  return gradient;
}


template <class VoxelType> void irtkGradientImageY<VoxelType>::Run()
{
  int x, y, z, t;

  // Do the initial set up
  this->Initialize();

  // Check image dimensions....
  if (this->_input->GetY() < 2) {
    cerr<<" irtkGradientImageY: Dimensions of input image are wrong"<<endl;
    exit(1);
  }

  for ( t = 0; t < this->_input->GetT(); ++t) {
    for ( z = 0; z < this->_input->GetZ(); ++z) {
      for ( x = 0; x < this->_input->GetX(); ++x) {
        this->_output->Put(x, 0, z, t, 0);
        for ( y = 1; y < this->_input->GetY()-1; ++y) {
          this->_output->PutAsDouble(x, y, z, t, this->Run(x, y, z, t));
        }
        this->_output->Put(x, y, z, t, 0);
      }
    }
  }

  // Do the final cleaning up
  this->Finalize();

}


template class  irtkGradientImageY<irtkBytePixel>;
template class  irtkGradientImageY<irtkGreyPixel>;
template class  irtkGradientImageY<irtkRealPixel>;
