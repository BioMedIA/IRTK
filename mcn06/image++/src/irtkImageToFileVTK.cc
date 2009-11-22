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

#include <irtkImageToFile.h>

#include <irtkVTK.h>

const char *irtkImageToFileVTK::NameOfClass()
{
  return "irtkImageToFileVTK";
}

void irtkImageToFileVTK::Initialize()
{
  irtkPoint p;
  int x, y, z;
  char header[2048];
  double xsize, ysize, zsize;

  // Initialize base class
  this->irtkImageToFile::Initialize();

  // Get image information from input
  x = _input->GetX();
  y = _input->GetY();
  z = _input->GetZ();

  double xaxis[3] = { 1, 0, 0 };
  double yaxis[3] = { 0, 1, 0 };
  double zaxis[3] = { 0, 0, 1 };

  _input->PutOrientation(xaxis, yaxis, zaxis);

  _input->ImageToWorld(p);
  _input->GetPixelSize(&xsize, &ysize, &zsize);

  // Construct header
  switch (this->_input->GetScalarType()) {
    case IRTK_VOXEL_CHAR: {
        sprintf(header, "%s\nCreated by %s\nBINARY\nDATASET STRUCTURED_POINTS\nDIMENSIONS %d %d %d\nORIGIN %f %f %f\nSPACING %f %f %f\nPOINT_DATA %d\nSCALARS scalars char 1\nLOOKUP_TABLE default\n", VTK_MAGIC2, this->NameOfClass(), x, y, z, p._x, p._y, p._z, xsize, ysize, zsize, x*y*z);
        break;
      }
    case IRTK_VOXEL_UNSIGNED_CHAR: {
        sprintf(header, "%s\nCreated by %s\nBINARY\nDATASET STRUCTURED_POINTS\nDIMENSIONS %d %d %d\nORIGIN %f %f %f\nSPACING %f %f %f\nPOINT_DATA %d\nSCALARS scalars unsigned char 1\nLOOKUP_TABLE default\n", VTK_MAGIC2, this->NameOfClass(), x, y, z, p._x, p._y, p._z, xsize, ysize, zsize, x*y*z);
        break;
      }
    case IRTK_VOXEL_SHORT: {
        sprintf(header, "%s\nCreated by %s\nBINARY\nDATASET STRUCTURED_POINTS\nDIMENSIONS %d %d %d\nORIGIN %f %f %f\nSPACING %f %f %f\nPOINT_DATA %d\nSCALARS scalars short 1\nLOOKUP_TABLE default\n", VTK_MAGIC2, this->NameOfClass(), x, y, z, p._x, p._y, p._z, xsize, ysize, zsize, x*y*z);
        break;
      }
    case IRTK_VOXEL_UNSIGNED_SHORT: {
        sprintf(header, "%s\nCreated by %s\nBINARY\nDATASET STRUCTURED_POINTS\nDIMENSIONS %d %d %d\nORIGIN %f %f %f\nSPACING %f %f %f\nPOINT_DATA %d\nSCALARS scalars unsigned short 1\nLOOKUP_TABLE default\n", VTK_MAGIC2, this->NameOfClass(), x, y, z, p._x, p._y, p._z, xsize, ysize, zsize, x*y*z);
        break;
      }
    case IRTK_VOXEL_FLOAT: {
        sprintf(header, "%s\nCreated by %s\nBINARY\nDATASET STRUCTURED_POINTS\nDIMENSIONS %d %d %d\nORIGIN %f %f %f\nSPACING %f %f %f\nPOINT_DATA %d\nSCALARS scalars float 1\nLOOKUP_TABLE default\n", VTK_MAGIC2, this->NameOfClass(), x, y, z, p._x, p._y, p._z, xsize, ysize, zsize, x*y*z);
        break;
      }
    case IRTK_VOXEL_DOUBLE: {
        sprintf(header, "%s\nCreated by %s\nBINARY\nDATASET STRUCTURED_POINTS\nDIMENSIONS %d %d %d\nORIGIN %f %f %f\nSPACING %f %f %f\nPOINT_DATA %d\nSCALARS scalars double 1\nLOOKUP_TABLE default\n", VTK_MAGIC2, this->NameOfClass(), x, y, z, p._x, p._y, p._z, xsize, ysize, zsize, x*y*z);
        break;
      }
    default:
      cerr << "irtkImageToFile::Run(): Unknown voxel type" << endl;
      exit(1);
  }

  // Write magic number
  this->Write(header, 0, strlen(header));

  // Calculate data address
  _start= strlen(header);
}
