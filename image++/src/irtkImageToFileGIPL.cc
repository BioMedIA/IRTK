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

#include <irtkGIPL.h>

const char *irtkImageToFileGIPL::NameOfClass()
{
  return "irtkImageToFileGIPL";
}

void irtkImageToFileGIPL::Initialize()
{
  int i;
  double min, max;
  double xsize, ysize, zsize, xaxis[3], yaxis[3], zaxis[3], pos[3];

  // Initialize base class
  this->irtkImageToFile::Initialize();

  // Write image dimensions
  this->WriteAsShort(_input->GetX(), 0);
  this->WriteAsShort(_input->GetY(), 2);
  this->WriteAsShort(_input->GetZ(), 4);
  this->WriteAsShort(1, 6);

  // Write type
  switch (this->_input->GetScalarType()) {
    case IRTK_VOXEL_CHAR: {
        this->WriteAsShort(GIPL_CHAR, 8);
      }
    case IRTK_VOXEL_UNSIGNED_CHAR: {
        this->WriteAsShort(GIPL_U_CHAR, 8);
      }
    case IRTK_VOXEL_SHORT: {
        this->WriteAsShort(GIPL_SHORT, 8);
      }
    case IRTK_VOXEL_UNSIGNED_SHORT: {
        this->WriteAsShort(GIPL_U_SHORT, 8);
      }
    case IRTK_VOXEL_FLOAT: {
        this->WriteAsShort(GIPL_FLOAT, 8);
      }
    case IRTK_VOXEL_DOUBLE: {
        this->WriteAsShort(GIPL_DOUBLE, 8);
      }
    default:
      cerr << "irtkImageToFileGIPL::Initialize: Not supported for this image type" << endl;
      exit(1);
  }
  // Write voxel dimensions
  _input->GetPixelSize(&xsize, &ysize, &zsize);
  this->WriteAsFloat(static_cast<float>(xsize), 10);
  this->WriteAsFloat(static_cast<float>(ysize), 14);
  this->WriteAsFloat(static_cast<float>(zsize), 18);
  this->WriteAsFloat(float(0), 22);

  // Write patient description
  for (i = 0; i < 80; i++) {
    this->WriteAsChar(char(0), 26+i*sizeof(char));
  }

  // Write image orientation
  _input->GetOrientation(xaxis, yaxis, zaxis);
  this->WriteAsFloat(static_cast<float>(xaxis[0]), 106);
  this->WriteAsFloat(static_cast<float>(xaxis[1]), 110);
  this->WriteAsFloat(static_cast<float>(xaxis[2]), 114);
  this->WriteAsFloat(static_cast<float>(yaxis[0]), 122);
  this->WriteAsFloat(static_cast<float>(yaxis[1]), 126);
  this->WriteAsFloat(static_cast<float>(yaxis[2]), 130);

  // Write identifier
  this->WriteAsInt(0, 154);

  // Write spare bytes
  for (i = 0; i < 28; i++) {
    this->WriteAsChar(char(0), 158+i*sizeof(char));
  }

  // Write flag and orientation
  this->WriteAsChar(char(0), 186);
  this->WriteAsChar(char(0), 187);

  // Write min and max values
  _input->GetMinMaxAsDouble(&min, &max);
  this->WriteAsDouble(min, 188);
  this->WriteAsDouble(max, 196);

  // Write origin
  _input->GetOrigin(pos[0], pos[1], pos[2]);
  this->WriteAsDouble(pos[0], 204);
  this->WriteAsDouble(pos[1], 212);
  this->WriteAsDouble(pos[2], 220);
  this->WriteAsDouble(double(0), 228);

  // Write magic number
  this->WriteAsFloat(float(0), 236);
  this->WriteAsFloat(float(0), 240);
  this->WriteAsUInt(GIPL_MAGIC_EXT, 244);
  this->WriteAsFloat(float(0), 248);
  this->WriteAsUInt(GIPL_MAGIC1, 252);

  // Calculate data address
  _start = GIPL_HEADERSIZE;
}


