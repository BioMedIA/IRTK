/*=========================================================================

  Library   : Image Registration Toolkit (IRTK)
  Module    : $Id$
  Copyright : Imperial College, Department of Computing
              Visual Information Processing (VIP), 2008 onwards
  Date      : $Date$
  Version   : $Revision$
  Changes   : $Author$

=========================================================================*/

#include <irtkCommon.h>

irtkCofstream::irtkCofstream()
{
#ifndef WORDS_BIGENDIAN
  _swapped = true;
#else
  _swapped = false;
#endif

#ifdef HAS_ZLIB
  _compressedFile = NULL;
#endif
  _uncompressedFile = NULL;
}

irtkCofstream::~irtkCofstream()
{
  this->Close();
}

void irtkCofstream::Write(char *data, long offset, long length)
{
  if (_compressed == false) {
    if (offset != -1) fseek(_uncompressedFile, offset, SEEK_SET);
    fwrite(data, length, 1, _uncompressedFile);
  } else {
#ifdef HAS_ZLIB
    if (offset != -1) {
      if (gztell(_compressedFile) > offset) {
          stringstream msg;
          msg << "Warning, writing compressed files only supports forward seek" << gztell(_compressedFile) << " " << offset << endl;
          cerr << msg.str();
          throw irtkException( msg.str(),
                               __FILE__,
                               __LINE__ );
      }
      gzseek(_compressedFile, offset, SEEK_SET);
    }
    gzwrite(_compressedFile, data, length);
#endif
  }
}

void irtkCofstream::WriteAsChar(char data, long offset)
{
  this->Write((char *)&data, offset, sizeof(char));
}

void irtkCofstream::WriteAsChar(char *data, long length, long offset)
{
  this->Write((char *)data, offset, length*sizeof(char));
}

void irtkCofstream::WriteAsUChar(unsigned char data, long offset)
{
  this->Write((char *)&data, offset, sizeof(unsigned char));
}

void irtkCofstream::WriteAsUChar(unsigned char *data, long length, long offset)
{
  this->Write((char *)data, offset, length*sizeof(unsigned char));
}

void irtkCofstream::WriteAsShort(short data, long offset)
{
  // Swap data
  if (_swapped == true) swap16((char *)&data, (char *)&data, 1);

  this->Write((char *)&data, offset, sizeof(short));
}

void irtkCofstream::WriteAsShort(short *data, long length, long offset)
{
  // Swap data
  if (_swapped == true) swap16((char *)data, (char *)data, length);

  this->Write((char *)data, offset, length*sizeof(short));

  // Swap data (BACK)
  if (_swapped == true) swap16((char *)data, (char *)data, length);
}

void irtkCofstream::WriteAsUShort(unsigned short data, long offset)
{
  // Swap data
  if (_swapped == true) swap16((char *)&data, (char *)&data, 1);

  this->Write((char *)&data, offset, sizeof(unsigned short));
}

void irtkCofstream::WriteAsUShort(unsigned short *data, long length, long offset)
{
  // Swap data
  if (_swapped == true) swap16((char *)data, (char *)data, length);

  this->Write((char *)data, offset, length*sizeof(unsigned short));

  // Swap data (BACK)
  if (_swapped == true) swap16((char *)data, (char *)data, length);
}

void irtkCofstream::WriteAsInt(int data, long offset)
{
  // Swap data
  if (_swapped == true) swap32((char *)&data, (char *)&data, 1);

  this->Write((char *)&data, offset, sizeof(int));
}

void irtkCofstream::WriteAsInt(int *data, long length, long offset)
{
  // Swap data
  if (_swapped == true) swap32((char *)data, (char *)data, length);

  this->Write((char *)data, offset, length*sizeof(int));

  // Swap data (BACK)
  if (_swapped == true) swap32((char *)data, (char *)data, length);
}

void irtkCofstream::WriteAsUInt(unsigned int data, long offset)
{
  // Swap data
  if (_swapped == true) swap32((char *)&data, (char *)&data, 1);

  this->Write((char *)&data, offset, sizeof(unsigned int));
}

void irtkCofstream::WriteAsUInt(unsigned int *data, long length, long offset)
{
  // Swap data
  if (_swapped == true) swap32((char *)data, (char *)data, length);

  this->Write((char *)data, offset, length*sizeof(unsigned int));

  // Swap data (BACK)
  if (_swapped == true) swap32((char *)data, (char *)data, length);
}

void irtkCofstream::WriteAsFloat(float data, long offset)
{
  // Swap data
  if (_swapped == true) swap32((char *)&data, (char *)&data, 1);

  this->Write((char *)&data, offset, sizeof(float));
}

void irtkCofstream::WriteAsFloat(float *data, long length, long offset)
{
  // Swap data
  if (_swapped == true) swap32((char *)data, (char *)data, length);

  this->Write((char *)data, offset, length*sizeof(float));

  // Swap data (BACK)
  if (_swapped == true) swap32((char *)data, (char *)data, length);
}

void irtkCofstream::WriteAsDouble(double *data, long length, long offset)
{
  // Swap data
  if (_swapped == true) swap64((char *)data, (char *)data, length);

  this->Write((char *)data, offset, length*sizeof(double));

  // Swap data (BACK)
  if (_swapped == true) swap64((char *)data, (char *)data, length);
}

void irtkCofstream::WriteAsDouble(double data, long offset)
{
  // Swap data
  if (_swapped == true) swap64((char *)&data, (char *)&data, 1);

  this->Write((char *)&data, offset, sizeof(double));
}

void irtkCofstream::WriteAsString(char *data, long offset)
{
  if (_compressed == false) {
    if (offset!= -1) fseek(_uncompressedFile, offset, SEEK_SET);
    fputs(data, _uncompressedFile);
  } else {
#ifdef HAS_ZLIB
    if (offset != -1) gzseek(_compressedFile, offset, SEEK_SET);
    gzputs(_compressedFile, data);
#endif
  }
}
