/*=========================================================================

  Library   : Image Registration Toolkit (IRTK)
  Module    : $Id$
  Copyright : Imperial College, Department of Computing
              Visual Information Processing (VIP), 2008 onwards
  Date      : $Date$
  Version   : $Revision$
  Changes   : $Author$

=========================================================================*/

#ifndef _IRTKCIFSTREAM_H

#define _IRTKCIFSTREAM_H

#include "irtkException.h"

/**
 * Class for reading compressed file streams.
 *
 * This class defines and implements functions for reading compressed file
 * streams. The file streams can be either uncompressed or compressed.
 */

class irtkCifstream : public irtkObject
{

  /// File pointer to potentially compressed file
#ifdef HAS_ZLIB
  gzFile _file;
#else
  FILE *_file;
#endif

#ifdef ENABLE_UNIX_COMPRESS
  /// Current position in the file
  long _pos;
#endif

protected:

  /// Flag whether file is swapped
  int _swapped;

public:

  /// Constructor
  irtkCifstream();

  /// Destructor
  ~irtkCifstream();

  /// Read n data as array (possibly compressed) from offset
  void Read(char *data, long length, long offset);

  /// Read n data as array of char (possibly compressed) from offset
  void ReadAsChar  (char *data, long length, long offset = -1);

  /// Read n data as array of unsigned char (possibly compressed) from offset
  void ReadAsUChar (unsigned char  *data, long length, long offset = -1);

  /// Read n data as array of short (possibly compressed) from offset
  void ReadAsShort (short *data, long length, long offset = -1);

  /// Read n data as array of unsigned short (possibly compressed) from offset
  void ReadAsUShort(unsigned short *data, long length, long offset = -1);

  /// Read n data as array of short (possibly compressed) from offset
  void ReadAsInt   (int *data, long length, long offset = -1);

  /// Read n data as array of unsigned short (possibly compressed) from offset
  void ReadAsUInt  (unsigned int *data, long length, long offset = -1);

  /// Read n data as array of short (possibly compressed) from offset
  void ReadAsFloat (float *data, long length, long offset = -1);

  /// Read n data as array of unsigned short (possibly compressed) from offset
  void ReadAsDouble(double *data, long length, long offset = -1);

  /// Read n data as string (possibly compressed) from current position
  void ReadAsString(char *data, long length, long offset = -1);

  /// Open file
  void Open(const char *);

  /// Close file
  void Close();

  /// Go to position in file
  void Seek(long);

  /// Current position in file
  long Tell();

  /// Returns whether file is swapped
  int  IsSwapped();

  /// Sets whether file is swapped
  void IsSwapped(int);

};

inline void irtkCifstream::Open(const char *filename)
{
#ifdef HAS_ZLIB
  _file = gzopen(filename, "rb");
#else
  _file = fopen(filename, "rb");
#endif

  // Check whether file was opened successful
  if (_file == NULL) {
      stringstream msg;
      msg << "cifstream::Open: Can't open file " << filename << endl;
      cerr << msg.str();
      throw irtkException( msg.str(),
                           __FILE__,
                           __LINE__ );
  }

#ifdef ENABLE_UNIX_COMPRESS
  _pos = 0;
#endif
}

inline void irtkCifstream::Close()
{
  if (_file != NULL) {
#ifdef HAS_ZLIB
    gzclose(_file);
#else
    fclose(_file);
#endif
    _file = NULL;
  }
#ifdef ENABLE_UNIX_COMPRESS
  _pos = 0;
#endif
}

inline int irtkCifstream::IsSwapped()
{
  return _swapped;
}

inline void irtkCifstream::IsSwapped(int swapped)
{
  _swapped = swapped;
}

inline long irtkCifstream::Tell()
{
#ifdef HAS_ZLIB
  return gztell(_file);
#else
  return ftell(_file);
#endif
}

inline void irtkCifstream::Seek(long offset)
{
#ifdef HAS_ZLIB
  gzseek(_file, offset, SEEK_SET);
#else
  fseek(_file, offset, SEEK_SET);
#endif
#ifdef ENABLE_UNIX_COMPRESS
  _pos = offset;
#endif
}

#endif
