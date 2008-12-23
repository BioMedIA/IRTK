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

#include <irtkConvolution.h>


template <class VoxelType> irtkConvolution_2D<VoxelType>::irtkConvolution_2D(Bool Normalization) :
    irtkConvolution<VoxelType>(Normalization)
{}

template <class VoxelType> Bool irtkConvolution_2D<VoxelType>::RequiresBuffering(void)
{
  return True;
}

template <class VoxelType> const char *irtkConvolution_2D<VoxelType>::NameOfClass()
{
  return "irtkConvolution_2D";
}

template <class VoxelType> void irtkConvolution_2D<VoxelType>::SetInput2(irtkGenericImage<irtkRealPixel> *image)
{
  if (image != NULL) {
    _input2 = image;
  } else {
    cerr << this->NameOfClass() << "::SetInput2: Input is not an image\n";
    exit(1);
  }
}


template <class VoxelType> double irtkConvolution_2D<VoxelType>::Run(int x, int y, int z, int t)
{
  float *ptr2;
  double val, sum;
  int x1, x2, y1, y2;

  // Initialize
  val = 0;
  sum = 0;
  x1 = x - this->_input2->GetX()/2;
  x2 = x + this->_input2->GetX()/2;
  y1 = y - this->_input2->GetY()/2;
  y2 = y + this->_input2->GetY()/2;

  // Check if we use normalization
  if (this->_Normalization == True) {
    // Check whether boundary checking is necessary
    if ((x1 > 0) && (x2 < this->_input->GetX()) &&
        (y1 > 0) && (y2 < this->_input->GetY())) {
      // If no, do fast convolution
      ptr2 = this->_input2->GetPointerToVoxels();
      for (y = y1; y <= y2; y++) {
        for (x = x1; x <= x2; x++) {
          val += *ptr2 * this->_input->Get(x, y, z, t);
          sum += *ptr2;
          ptr2++;
        }
      }
    } else {
      // If yes, do slow convolution which handles boundaries
      ptr2 = this->_input2->GetPointerToVoxels();
      for (y = y1; y <= y2; y++) {
        for (x = x1; x <= x2; x++) {
          if ((x >= 0) && (x < this->_input->GetX()) &&
              (y >= 0) && (y < this->_input->GetY())) {
            val += *ptr2 * this->_input->Get(x, y, z, t);
            sum += *ptr2;
          }
          ptr2++;
        }
      }
    }
    //  Normalize filter value by sum of filter elements
    if (sum > 0) {
      return val / sum;
    } else {
      return 0;
    }
  } else {
    // Check whether boundary checking is necessary
    if ((x1 > 0) && (x2 < this->_input->GetX()) &&
        (y1 > 0) && (y2 < this->_input->GetY())) {
      // If no, do fast convolution
      ptr2 = this->_input2->GetPointerToVoxels();
      for (y = y1; y <= y2; y++) {
        for (x = x1; x <= x2; x++) {
          val += *ptr2 * this->_input->Get(x, y, z, t);
          sum += *ptr2;
          ptr2++;
        }
      }
    } else {
      // If yes, do slow convolution which handles boundaries
      ptr2 = this->_input2->GetPointerToVoxels();
      for (y = y1; y <= y2; y++) {
        for (x = x1; x <= x2; x++) {
          if ((x >= 0) && (x < this->_input->GetX()) &&
              (y >= 0) && (y < this->_input->GetY())) {
            val += *ptr2 * this->_input->Get(x, y, z, t);
            sum += *ptr2;
          }
          ptr2++;
        }
      }
    }
    return val;
  }
}

template <class VoxelType> void irtkConvolution_2D<VoxelType>::Initialize()
{
  // Check kernel
  if (this->_input2 == NULL) {
    cerr << this->NameOfClass() << "::Run: Filter has no second input" << endl;
    exit(1);
  }

  // Check kernel size
  if (this->_input2->GetZ() != 1) {
    cerr << this->NameOfClass();
    cerr << "::Run: Filter dimensions should be 1 in Z" << endl;
    exit(1);
  }

  // Do the initial set up
  this->irtkImageToImage<VoxelType>::Initialize();
}

template class irtkConvolution_2D<irtkBytePixel>;
template class irtkConvolution_2D<irtkGreyPixel>;
template class irtkConvolution_2D<irtkRealPixel>;
