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

#include <irtkFileToImage.h>
#include <irtkNIFTI.h>

irtkBaseImage::irtkBaseImage()
{
  _matI2W = irtkMatrix(4, 4);
  _matW2I = irtkMatrix(4, 4);

  // Update transformation matrix
  this->UpdateMatrix();
}

irtkBaseImage::irtkBaseImage(const irtkBaseImage &image) : irtkObject(image)
{
  _attr   = image._attr;
  _matI2W = image._matI2W;
  _matW2I = image._matW2I;
}

irtkBaseImage::~irtkBaseImage()
{
	// Update image attributes
	_attr = irtkImageAttributes();

  // Update transformation matrix
  this->UpdateMatrix();
}

irtkBaseImage *irtkBaseImage::New(const char *filename)
{
  irtkBaseImage *image;

  irtkFileToImage *reader = irtkFileToImage::New(filename);
  image = reader->GetOutput();
  delete reader;

  return image;
}

irtkBaseImage *irtkBaseImage::New(const irtkBaseImage *image)
{
  // Convert image
  if (dynamic_cast<const irtkGenericImage<char> *>(image) != NULL) {
    return new irtkGenericImage<char>(*(dynamic_cast<const irtkGenericImage<char> *>(image)));
  } else if (dynamic_cast<const irtkGenericImage<unsigned char> *>(image) != NULL) {
    return new irtkGenericImage<unsigned char>(*(dynamic_cast<const irtkGenericImage<unsigned char> *>(image)));
  } else if (dynamic_cast<const irtkGenericImage<short> *>(image) != NULL) {
    return new irtkGenericImage<short>(*(dynamic_cast<const irtkGenericImage<short> *>(image)));
  } else if (dynamic_cast<const irtkGenericImage<unsigned short> *>(image) != NULL) {
    return new irtkGenericImage<unsigned short>(*(dynamic_cast<const irtkGenericImage<unsigned short> *>(image)));
  } else if (dynamic_cast<const irtkGenericImage<float> *>(image) != NULL) {
    return new irtkGenericImage<float>(*(dynamic_cast<const irtkGenericImage<float> *>(image)));
  } else if (dynamic_cast<const irtkGenericImage<double> *>(image) != NULL) {
    return new irtkGenericImage<double>(*(dynamic_cast<const irtkGenericImage<double> *>(image)));
  } else {
  	cerr << "irtkBaseImage::New: Cannot allocate image of unknown type" << endl;
  	exit(1);
  }
}

void irtkBaseImage::UpdateMatrix()
{
  // Note that the update of zaxis is taken place in the
  // Initialize routines, in case a specific zaxis is passed.

  // Update image to world coordinate system matrix
  _matI2W.Ident();
  _matI2W(0, 0) = _attr._xaxis[0];
  _matI2W(1, 0) = _attr._xaxis[1];
  _matI2W(2, 0) = _attr._xaxis[2];
  _matI2W(0, 1) = _attr._yaxis[0];
  _matI2W(1, 1) = _attr._yaxis[1];
  _matI2W(2, 1) = _attr._yaxis[2];
  _matI2W(0, 2) = _attr._zaxis[0];
  _matI2W(1, 2) = _attr._zaxis[1];
  _matI2W(2, 2) = _attr._zaxis[2];

  irtkMatrix tmp1(4, 4);
  tmp1.Ident();
  tmp1(0, 3) = - (_attr._x - 1) / 2.0;
  tmp1(1, 3) = - (_attr._y - 1) / 2.0;
  tmp1(2, 3) = - (_attr._z - 1) / 2.0;

  irtkMatrix tmp2(4, 4);
  tmp2.Ident();
  tmp2(0, 0) = _attr._dx;
  tmp2(1, 1) = _attr._dy;
  tmp2(2, 2) = _attr._dz;

  irtkMatrix tmp3(4, 4);
  tmp3.Ident();
  tmp3(0, 3) = _attr._xorigin;
  tmp3(1, 3) = _attr._yorigin;
  tmp3(2, 3) = _attr._zorigin;

  _matI2W = tmp3 * (_matI2W * (tmp2 * tmp1));

  // Update world to image coordinate system matrix
  _matW2I.Ident();
  _matW2I(0, 0) = _attr._xaxis[0];
  _matW2I(0, 1) = _attr._xaxis[1];
  _matW2I(0, 2) = _attr._xaxis[2];
  _matW2I(1, 0) = _attr._yaxis[0];
  _matW2I(1, 1) = _attr._yaxis[1];
  _matW2I(1, 2) = _attr._yaxis[2];
  _matW2I(2, 0) = _attr._zaxis[0];
  _matW2I(2, 1) = _attr._zaxis[1];
  _matW2I(2, 2) = _attr._zaxis[2];

  tmp1.Ident();
  tmp1(0, 3) = (_attr._x - 1) / 2.0;
  tmp1(1, 3) = (_attr._y - 1) / 2.0;
  tmp1(2, 3) = (_attr._z - 1) / 2.0;

  tmp2.Ident();
  tmp2(0, 0) = 1.0 / _attr._dx;
  tmp2(1, 1) = 1.0 / _attr._dy;
  tmp2(2, 2) = 1.0 / _attr._dz;

  tmp3.Ident();
  tmp3(0, 3) = - _attr._xorigin;
  tmp3(1, 3) = - _attr._yorigin;
  tmp3(2, 3) = - _attr._zorigin;

  _matW2I = tmp1 * (tmp2 * (_matW2I * tmp3));
}

void irtkBaseImage::Update(const irtkImageAttributes &attr)
{
  _attr = attr;
 
  // Update transformation matrix
  this->UpdateMatrix();
}

void irtkBaseImage::Orientation(int &i, int &j, int &k) const
{

#ifdef HAS_NIFTI

  // Work out orientation of axis
  mat44 mat_44;
  for (i = 0; i < 4; ++i) {
    for (j = 0; j < 4; ++j) {
      mat_44.m[i][j] = this->GetImageToWorldMatrix()(i, j);
    }
  }
  nifti_mat44_to_orientation(mat_44, &i, &j, &k);

#else

  cerr << "irtkBaseImage::Orientation: Requires NIFTI support. Please recompile with NIFTI enabled" << endl;
  exit(1);

#endif
}

void irtkBaseImage::GetMinMaxAsDouble(double *min, double *max) const
{
  int x, y, z, t;

  *min = 0;
  *max = 0;

  if (this->GetNumberOfVoxels() > 0) {
    *min = this->GetAsDouble(0, 0, 0, 0);
    *max = this->GetAsDouble(0, 0, 0, 0);
  }

  for (t = 0; t < this->GetT(); t++) {
    for (z = 0; z < this->GetZ(); z++) {
      for (y = 0; y < this->GetY(); y++) {
        for (x = 0; x < this->GetX(); x++) {
          if (this->GetAsDouble(x, y, z, t) < *min) *min = this->GetAsDouble(x, y, z, t);
          if (this->GetAsDouble(x, y, z, t) > *max) *max = this->GetAsDouble(x, y, z, t);
        }
      }
    }
  }
}

void irtkBaseImage::PutMinMaxAsDouble(double min, double max)
{
  int x, y, z, t;
  double min_val, max_val;

  // Get lower and upper bound
  this->GetMinMaxAsDouble(&min_val, &max_val);

  for (t = 0; t < this->GetT(); t++) {
    for (z = 0; z < this->GetZ(); z++) {
      for (y = 0; y < this->GetY(); y++) {
        for (x = 0; x < this->GetX(); x++) {
          this->PutAsDouble(x, y, z, t, ((this->GetAsDouble(x, y, z, t) - min_val) / (max_val - min_val)) * (max - min) + min);
        }
      }
    }
  }
}

void irtkBaseImage::Print()
{
  // Print image dimensions
  cout << "Image size is " << _attr._x << " " << _attr._y << " " << _attr._z << " " << _attr._t << endl;
  // Print voxel dimensions
  cout << "Voxel size is " << _attr._dx << " " << _attr._dy << " "
  << _attr._dz << " " << _attr._dt << endl;
  // Print origin
  cout << "Image origin is " << _attr._xorigin << " " << _attr._yorigin << " " << _attr._zorigin << " " << _attr._torigin << endl;
  // Print x-axis
  cout << "X-axis is " << _attr._xaxis[0] << " " << _attr._xaxis[1] << " " << _attr._xaxis[2] << endl;
  // Print x-axis
  cout << "Y-axis is " << _attr._yaxis[0] << " " << _attr._yaxis[1] << " " << _attr._yaxis[2] << endl;
  // Print x-axis
  cout << "Z-axis is " << _attr._zaxis[0] << " " << _attr._zaxis[1] << " " << _attr._zaxis[2] << endl;
}
