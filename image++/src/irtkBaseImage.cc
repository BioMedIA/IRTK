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
      stringstream msg;
      msg << "irtkBaseImage::New: Cannot allocate image of unknown type" << endl;
      cerr << msg.str();
      throw irtkException( msg.str(),
                           __FILE__,
                           __LINE__ );
  }
}

irtkMatrix irtkBaseImage::GetImageToWorldMatrix(const irtkImageAttributes & attr)
{
  irtkMatrix translate1(4, 4);
  translate1.Ident();
  translate1(0, 3) = - (attr._x - 1) / 2.0;
  translate1(1, 3) = - (attr._y - 1) / 2.0;
  translate1(2, 3) = - (attr._z - 1) / 2.0;
  
  irtkMatrix scale(4, 4);
  scale(0, 0) = attr._dx;
  scale(1, 1) = attr._dy;
  scale(2, 2) = attr._dz;
  scale(3, 3) = 1.0;
  
  irtkMatrix rot(4, 4);
  rot(0, 0) = attr._xaxis[0];
  rot(1, 0) = attr._xaxis[1];
  rot(2, 0) = attr._xaxis[2];
  rot(0, 1) = attr._yaxis[0];
  rot(1, 1) = attr._yaxis[1];
  rot(2, 1) = attr._yaxis[2];
  rot(0, 2) = attr._zaxis[0];
  rot(1, 2) = attr._zaxis[1];
  rot(2, 2) = attr._zaxis[2];
  rot(3, 3) = 1;

  irtkMatrix translate2(4, 4);
  translate2.Ident();
  translate2(0, 3) = attr._xorigin;
  translate2(1, 3) = attr._yorigin;
  translate2(2, 3) = attr._zorigin;
  
  return translate2 * (rot * (scale * translate1));
}

irtkMatrix irtkBaseImage::GetWorldToImageMatrix(const irtkImageAttributes & attr)
{
  irtkMatrix translate1(4, 4);
  translate1.Ident();
  translate1(0, 3) = - attr._xorigin;
  translate1(1, 3) = - attr._yorigin;
  translate1(2, 3) = - attr._zorigin;
  
  irtkMatrix rot(4, 4);
  rot(0, 0) = attr._xaxis[0];
  rot(0, 1) = attr._xaxis[1];
  rot(0, 2) = attr._xaxis[2];
  rot(1, 0) = attr._yaxis[0];
  rot(1, 1) = attr._yaxis[1];
  rot(1, 2) = attr._yaxis[2];
  rot(2, 0) = attr._zaxis[0];
  rot(2, 1) = attr._zaxis[1];
  rot(2, 2) = attr._zaxis[2];
  rot(3, 3) = 1.0;
  
  irtkMatrix scale(4, 4);
  scale(0, 0) = 1.0 / attr._dx;
  scale(1, 1) = 1.0 / attr._dy;
  scale(2, 2) = 1.0 / attr._dz;
  scale(3, 3) = 1.0;
  
  irtkMatrix translate2(4, 4);
  translate2.Ident();
  translate2(0, 3) = (attr._x - 1) / 2.0;
  translate2(1, 3) = (attr._y - 1) / 2.0;
  translate2(2, 3) = (attr._z - 1) / 2.0;
  
  return translate2 * (scale * (rot * translate1));
}

void irtkBaseImage::UpdateMatrix()
{
  _matI2W = GetImageToWorldMatrix(_attr);
  _matW2I = GetWorldToImageMatrix(_attr);
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
  stringstream msg;
  msg << "irtkBaseImage::Orientation: Requires NIFTI support. Please recompile with NIFTI enabled" << endl;
  cerr << msg.str();
  throw irtkException( msg.str(),
                           __FILE__,
                           __LINE__ );
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
  char buffer[256];
  // Print image dimensions
  cout << "Image size is " << _attr._x << " " << _attr._y << " " << _attr._z << " " << _attr._t << endl;
  // Print voxel dimensions
  cout << "Voxel size is " << _attr._dx << " " << _attr._dy << " "
  << _attr._dz << " " << _attr._dt << endl;
  // Print origin
  cout << "Image origin is " << _attr._xorigin << " " << _attr._yorigin << " " << _attr._zorigin << " " << _attr._torigin << endl;
  // Print x-axis
  sprintf(buffer,"X-axis is % 0.2f % 0.2f % 0.2f\n", _attr._xaxis[0], _attr._xaxis[1], _attr._xaxis[2]);
  cout<<buffer;
  // Print y-axis
  sprintf(buffer,"Y-axis is % 0.2f % 0.2f % 0.2f\n", _attr._yaxis[0], _attr._yaxis[1], _attr._yaxis[2]);
  cout<<buffer;
  // Print z-axis
  sprintf(buffer,"Z-axis is % 0.2f % 0.2f % 0.2f\n", _attr._zaxis[0], _attr._zaxis[1], _attr._zaxis[2]);
  cout<<buffer;
}
