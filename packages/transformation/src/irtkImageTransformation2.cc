/*=========================================================================

  Library   : Image Registration Toolkit (IRTK)
  Module    : $Id$
  Copyright : Imperial College, Department of Computing
              Visual Information Processing (VIP), 2008 onwards
  Date      : $Date$
  Version   : $Revision$
  Changes   : $Author$

=========================================================================*/

#include <irtkTransformation.h>

template <class VoxelType> irtkImageTransformation2<VoxelType>::irtkImageTransformation2()
{
  // Set input and output
  _input  = NULL;
  _output = NULL;

  // Set transformation
  _transformation = NULL;
  _transformation2 = NULL;

  // Set interpolator
  _interpolator = NULL;

  // Set padding value
  _TargetPaddingValue = std::numeric_limits<VoxelType>::min();

  // Set padding value
  _SourcePaddingValue = 0;

  // Set invert mode
  _Invert = False;
}

template <class VoxelType> irtkImageTransformation2<VoxelType>::~irtkImageTransformation2()
{
  // Set input and output
  _input  = NULL;
  _output = NULL;

  // Set transformation
  _transformation = NULL;
  _transformation2 = NULL;

  // Set interpolator
  _interpolator = NULL;

  // Set padding value
  _TargetPaddingValue = std::numeric_limits<VoxelType>::min();

  // Set padding value
  _SourcePaddingValue = 0;

  // Set invert mode
  _Invert = False;
}

template <class VoxelType> void irtkImageTransformation2<VoxelType>::SetTransformation(irtkTransformation *transformation, irtkTransformation *transformation2)
{
  if ((transformation != NULL) || (transformation2 != NULL)) {
    _transformation = transformation;
    _transformation2 = transformation2;
  } else {
    cerr << "irtkImageTransformation::SetInput: Transformation is NULL\n";
    exit(1);
  }
}

template <class VoxelType> void irtkImageTransformation2<VoxelType>::SetInput(irtkGenericImage<VoxelType> *image)
{
  if (image != NULL) {
    _input = image;
  } else {
    cerr << "irtkImageTransformation::SetInput: Input is NULL\n";
    exit(1);
  }
}

template <class VoxelType> void irtkImageTransformation2<VoxelType>::SetOutput(irtkGenericImage<VoxelType> *image)
{
  if (image != NULL) {
    _output = image;
  } else {
    cerr << "irtkImageTransformation2::SetOutput: Output is NULL\n";
    exit(1);
  }
}

template <class VoxelType> void irtkImageTransformation2<VoxelType>::Run()
{
  int i, j, k, l, t;
  double x, y, z;
  VoxelType *ptr2output;

  // Check inputs and outputs
  if (_input == NULL) {
    cerr << "irtkImageTransformation2::Run: Filter has no input" << endl;
    exit(1);
  }

  if (_output == NULL) {
    cerr << "irtkImageTransformation2::Run: Filter has no output" << endl;
    exit(1);
  }

  if (_transformation == NULL) {
    cerr << "irtkImageTransformation2::Run: Filter has no transformation" << endl;
    exit(1);
  }

  if (_transformation2 == NULL) {
    cerr << "irtkImageTransformation2::Run: Filter has no transformation2" << endl;
    exit(1);
  }

  if (_interpolator == NULL) {
    cerr << "irtkImageTransformation2::Run: Filter has no interpolator" << endl;
    exit(1);
  }

  if (_input->IsEmpty() == True) {
    cerr << "irtkImageTransformation2::Run: Input is empty" << endl;
    exit(1);
  }

  if (_input == _output) {
    cerr << "irtkImageTransformation2::Run: Input equals output" << endl;
    exit(1);
  }

  // Setup interpolation
  _interpolator->SetInput(_input);
  _interpolator->Initialize();

  // Pointer to voxels in output image
  ptr2output = _output->GetPointerToVoxels();

  // Calculate transformation
  for (l = 0; l < _output->GetT(); l++) {
    t = round(this->_input->TimeToImage(this->_output->ImageToTime(l)));

    if ((t >= 0) || (t < this->_input->GetT())) {

      // Calculate time
      double time = this->_output->ImageToTime(l);

      // Calculate transformation
      for (k = 0; k < _output->GetZ(); k++) {
        for (j = 0; j < _output->GetY(); j++) {
          for (i = 0; i < _output->GetX(); i++) {
            if (*ptr2output > _TargetPaddingValue) {
              x = i;
              y = j;
              z = k;
              // Transform point into world coordinates
              _output->ImageToWorld(x, y, z);
              if (_Invert == False) {
                // Transform point
                _transformation->Transform(x, y, z, time);
                _transformation2->Transform(x, y, z, time);
              } else {
                // Transform point
                _transformation2->Inverse(x, y, z, time);
                _transformation->Inverse(x, y, z, time);
              }
              // Transform point into image coordinates
              _input->WorldToImage(x, y, z);
              // Check whether transformed point is in FOV of input
              if ((x > -0.5) && (x < _input->GetX()-0.5) &&
                  (y > -0.5) && (y < _input->GetY()-0.5) &&
                  (z > -0.5) && (z < _input->GetZ()-0.5)) {
                *ptr2output = round(_interpolator->Evaluate(x, y, z, t));
              } else {
                // Fill with padding value
                *ptr2output = _SourcePaddingValue;
              }
            } else {
              // Fill with padding value
              *ptr2output = _SourcePaddingValue;
            }
            ptr2output++;
          }
        }
      }
    } else {
      for (k = 0; k < this->_output->GetZ(); k++) {
        for (j = 0; j < this->_output->GetY(); j++) {
          for (i = 0; i < this->_output->GetX(); i++) {
            *ptr2output = this->_SourcePaddingValue;
            ptr2output++;
          }
        }
      }
    }
  }
}

template class irtkImageTransformation2<irtkRealPixel>;
template class irtkImageTransformation2<irtkGreyPixel>;
template class irtkImageTransformation2<irtkBytePixel>;
