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

irtkCifstream::irtkCifstream()
{
  _file = NULL;
#ifndef WORDS_BIGENDIAN
  _swapped = true;
#else
  _swapped = false;
#endif
#ifdef ENABLE_UNIX_COMPRESS
  _pos = 0;
#endif
}

irtkCifstream::~irtkCifstream()
{
  this->Close();
}

void irtkCifstream::Read(char *mem, int start, int num)
{
  // Read data uncompressed
#ifdef ENABLE_UNIX_COMPRESS
  if (start == -1) {
    ReadCompressed(_file, mem, _pos, num);
    _pos = _pos + num;
  } else {
    ReadCompressed(_file, mem, start, num);
    _pos = start + num;
  }
#else
#ifdef HAS_ZLIB
  if (start != -1) gzseek(_file, start, SEEK_SET);
  gzread(_file, mem, num);
#else
  if (start != -1) fseek(_file, start, SEEK_SET);
  fread(mem, num, 1, _file);
#endif
#endif
}

void irtkCifstream::ReadAsChar(char *data, int length, int offset)
{
  // Read data (possibly compressed)
  this->Read((char *)data, offset, length * sizeof(char));
}

void irtkCifstream::ReadAsUChar(unsigned char *data, int length, int offset)
{
  // Read data (possibly compressed)
  this->Read((char *)data, offset, length * sizeof(unsigned char));
}

void irtkCifstream::ReadAsShort(short *data, int length, int offset)
{
  // Read data (possibly compressed)
  this->Read((char *)data, offset, length * sizeof(short));

  // Swap data
  if (_swapped == true) swap16((char *)data, (char *)data, length);
}

void irtkCifstream::ReadAsUShort(unsigned short *data, int length, int offset)
{
  // Read data (possibly compressed)
  this->Read((char *)data, offset,
             length * sizeof(unsigned short));

  // Swap data
  if (_swapped == true) swap16((char *)data, (char *)data, length);
}

void irtkCifstream::ReadAsInt(int *data, int length, int offset)
{
  // Read data (possibly compressed)
  this->Read((char *)data, offset, length * sizeof(int));

  // Swap data
  if (_swapped == true) swap32((char *)data, (char *)data, length);
}

void irtkCifstream::ReadAsUInt(unsigned int *data, int length, int offset)
{
  // Read data (possibly compressed)
  this->Read((char *)data, offset,
             length * sizeof(unsigned int));

  // Swap data
  if (_swapped == true) swap32((char *)data, (char *)data, length);
}

void irtkCifstream::ReadAsFloat(float *data, int length, int offset)
{
  // Read data (possibly compressed)
  this->Read((char *)data, offset, length * sizeof(float));

  // Swap data
  if (_swapped == true) swap32((char *)data, (char *)data, length);
}

void irtkCifstream::ReadAsDouble(double *data, int length, int offset)
{
  // Read data (possibly compressed)
  this->Read((char *)data, offset, length * sizeof(double));

  // Swap data
  if (_swapped == true) swap64((char *)data, (char *)data, length);
}

void irtkCifstream::ReadAsString(char *data, int length, int offset)
{
  // Read string
#ifdef HAS_ZLIB
  if (offset != -1) gzseek(_file, offset, SEEK_SET);
  gzgets(_file, data, length);
#else
  if (offset!= -1) fseek(_file, offset, SEEK_SET);
  fgets(data, length, _file);
#endif

  // Check for UNIX end-of-line char
  if ((strlen(data) > 0) && (data[strlen(data)-1] = '\n')) {
    data[strlen(data)-1] = '\0';
    return;
  }
  // Check for MAC end-of-line char
  if ((strlen(data) > 0) && (data[strlen(data)-1] = '\r')) {
    data[strlen(data)-1] = '\0';
    return;
  }
  // Check for Windows end-of-line char
  if ((strlen(data) > 1) && (data[strlen(data)-2] = '\r') && (data[strlen(data)-1] = '\n')) {
    data[strlen(data)-2] = '\0';
    return;
  }
}

