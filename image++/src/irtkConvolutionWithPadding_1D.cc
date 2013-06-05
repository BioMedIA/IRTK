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

template <class VoxelType> irtkConvolutionWithPadding_1D<VoxelType>::irtkConvolutionWithPadding_1D(VoxelType padding,
    bool Normalization) : irtkConvolution_1D<VoxelType>(Normalization)
{
  _padding = padding;
}

template <class VoxelType> const char *irtkConvolutionWithPadding_1D<VoxelType>::NameOfClass()
{
  return "irtkConvolutionWithPadding_1D";
}

template <class VoxelType> void irtkConvolutionWithPadding_1D<VoxelType>::PutPaddingValue(VoxelType padding)
{
  _padding = padding;
}

template <class VoxelType> VoxelType irtkConvolutionWithPadding_1D<VoxelType>::GetPaddingValue(void)
{
  return _padding;
}

template <class VoxelType> double irtkConvolutionWithPadding_1D<VoxelType>::Run(int x, int y, int z, int t)
{
  int x1, x2;
  double *ptr2;
  VoxelType *ptr;
  double val, sum;

  if (this->_input->Get(x, y, z, t) <= this->_padding) return this->_padding;

  // Initialize
  val = 0;
  sum = 0;
  x1 = x - this->_input2->GetX()/2;
  x2 = x + this->_input2->GetX()/2;

  // Check if we use normalization
  if (this->_Normalization == true) {
    // Check whether boundary checking is necessary
    if ((x1 > 0) && (x2 < this->_input->GetX())) {

      // If no, do fast Convolution
      ptr  = this->_input->GetPointerToVoxels(x1, y, z, t);
      ptr2 = this->_input2->GetPointerToVoxels();
      for (x = x1; x <= x2; x++) {
        if (*ptr > this->_padding) {
          val += *ptr2 * *ptr;
          sum += *ptr2;
        }
        ptr++;
        ptr2++;
      }
    } else {

      // If yes, do slow convolution which handles boundaries
      ptr2 = this->_input2->GetPointerToVoxels();
      for (x = x1; x <= x2; x++) {
        if ((x >= 0) && (x < this->_input->GetX())) {
          if (this->_input->Get(x, y, z, t) > this->_padding) {
            val += *ptr2 * this->_input->Get(x, y, z, t);
            sum += *ptr2;
          }
        }
        ptr2++;
      }
    }

    if (sum > 0) {
      return val / sum;
    } else {
      return 0;
    }
  } else {
    // Check whether boundary checking is necessary
    if ((x1 > 0) && (x2 < this->_input->GetX())) {

      // If no, do fast Convolution
      ptr  = this->_input->GetPointerToVoxels(x1, y, z, t);
      ptr2 = this->_input2->GetPointerToVoxels();
      for (x = x1; x <= x2; x++) {
        if (*ptr > this->_padding) {
          val += *ptr2 * *ptr;
          sum += *ptr2;
        }
        ptr++;
        ptr2++;
      }
    } else {

      // If yes, do slow convolution which handles boundaries
      ptr2 = this->_input2->GetPointerToVoxels();
      for (x = x1; x <= x2; x++) {
        if ((x >= 0) && (x < this->_input->GetX())) {
          if (this->_input->Get(x, y, z, t) > this->_padding) {
            val += *ptr2 * this->_input->Get(x, y, z, t);
            sum += *ptr2;
          }
        }
        ptr2++;
      }
    }

    return val;
  }
}

template class irtkConvolutionWithPadding_1D<unsigned char>;
template class irtkConvolutionWithPadding_1D<short>;
template class irtkConvolutionWithPadding_1D<unsigned short>;
template class irtkConvolutionWithPadding_1D<float>;
template class irtkConvolutionWithPadding_1D<double>;
