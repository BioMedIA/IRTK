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

#include <irtkImageFunction.h>

template <class VoxelType> irtkImageFunction<VoxelType>::irtkImageFunction()
{
  // Set input
  _input  = NULL;

  // Default parameters
  _DebugFlag    = False;
  _DefaultValue = 0;
}

template <class VoxelType> irtkImageFunction<VoxelType>::~irtkImageFunction()
{
  // Set input
  _input  = NULL;
}

template <class VoxelType> void irtkImageFunction<VoxelType>::SetInput(irtkGenericImage<VoxelType> *image)
{
  if (image != NULL) {
    _input = image;
  } else {
    cerr << "irtkImageFunction::SetInput: Input is not an image\n";
    exit(1);
  }
}

template <class VoxelType> void irtkImageFunction<VoxelType>::Debug(const char *message)
{
  if (_DebugFlag == True) cout << message << endl;
}

template <class VoxelType> void irtkImageFunction<VoxelType>::Initialize()
{
  // Print debugging information
  this->Debug("irtkImageFunction::Initialize");

  // Check inputs and outputs
  if (_input == NULL) {
    cerr << this->NameOfClass() << "::Run: Filter has no input" << endl;
    exit(1);
  }
}

template class irtkImageFunction<irtkBytePixel>;
template class irtkImageFunction<irtkGreyPixel>;
template class irtkImageFunction<irtkRealPixel>;
