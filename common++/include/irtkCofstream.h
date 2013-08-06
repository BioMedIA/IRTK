/*=========================================================================

  Library   : Image Registration Toolkit (IRTK)
  Module    : $Id$
  Copyright : Imperial College, Department of Computing
              Visual Information Processing (VIP), 2008 onwards
  Date      : $Date$
  Version   : $Revision$
  Changes   : $Author$

=========================================================================*/

#ifndef _IRTKCOFSTREAM_H

#define _IRTKCOFSTREAM_H

/**
 * Class for writing compressed file streams.
 *
 * This class defines and implements functions for writing compressed file
 * streams. At the moment only the writing of uncompressed file streams is
 * supported.
 */

#include "irtkException.h"

class irtkCofstream : public irtkObject
{

  /// File pointer to uncompressed file
  FILE *_uncompressedFile;

#ifdef HAS_ZLIB
  /// File pointer to compressed file
  gzFile _compressedFile;
#endif

protected:

  /// Flag whether file is compressed
  int _compressed;

  /// Flag whether file is swapped
  int _swapped;

public:

  /// Constructor
  irtkCofstream();

  /// Destructor
  ~irtkCofstream();

  /// Write n data as array from offset
  void Write(char *data, long offset, long length);

  /// Write data as char (possibly compressed) from offset
  void WriteAsChar  (char data, long offset = -1);
  /// Write data as char (possibly compressed) from offset
  void WriteAsChar  (char *data, long length, long offset = -1);

  /// Write data as unsigned char (possibly compressed) from offset
  void WriteAsUChar (unsigned char data, long offset = -1);
  /// Write data as unsigned char (possibly compressed) from offset
  void WriteAsUChar (unsigned char *data, long length, long offset = -1);

  /// Write data as short (possibly compressed) from offset
  void WriteAsShort (short data, long offset = -1);
  /// Write data as short (possibly compressed) from offset
  void WriteAsShort (short *data, long length, long offset = -1);

  /// Write data as unsigned short (possibly compressed) from offset
  void WriteAsUShort(unsigned short data, long offset = -1);
  /// Write data as unsigned short (possibly compressed) from offset
  void WriteAsUShort(unsigned short *data, long length, long offset = -1);

  /// Write data as int (possibly compressed) from offset
  void WriteAsInt   (int data, long offset = -1);
  /// Write data as int (possibly compressed) from offset
  void WriteAsInt   (int *data, long length, long offset = -1);

  /// Write data as unsigned int (possibly compressed) from offset
  void WriteAsUInt  (unsigned int data, long offset = -1);
  /// Write data as unsigned int (possibly compressed) from offset
  void WriteAsUInt  (unsigned int *data, long length, long offset = -1);

  /// Write data as float (possibly compressed) from offset
  void WriteAsFloat (float data, long offset = -1);
  /// Write data as float (possibly compressed) from offset
  void WriteAsFloat (float *data, long length, long offset = -1);

  /// Write data as double (possibly compressed) from offset
  void WriteAsDouble(double data, long offset = -1);
  /// Write data as double (possibly compressed) from offset
  void WriteAsDouble(double *data, long length, long offset = -1);

  /// Write data as string (possibly compressed)
  void WriteAsString(char *data, long offset = -1);

  /// Open file
  void Open(const char *);

  /// Close file
  void Close();

  /// Returns whether file is compressed
  int  IsCompressed();

  /// Sets whether file is compressed
  void IsCompressed(int);

  /// Returns whether file is swapped
  int  IsSwapped();

  /// Sets whether file is swapped
  void IsSwapped(int);

};

inline void irtkCofstream::Open(const char *filename)
{
  if (strstr(basename2(filename), ".gz") == NULL) {
    _compressed = false;
    _uncompressedFile = fopen(filename, "wb");

    // Check whether file was opened successful
    if (_uncompressedFile == NULL) {
        stringstream msg;
        msg << "cofstream::Open: Can't open file " << filename << endl;
        cerr << msg.str();
        throw irtkException( msg.str(),
                           __FILE__,
                           __LINE__ );
    }
  } else {
#ifdef HAS_ZLIB
    _compressed = true;
    _compressedFile = gzopen(filename, "wb");

    // Check whether file was opened successful
    if (_compressedFile == NULL) {
        stringstream msg;
        msg << "cofstream::Open: Can't open file " << filename << endl;
        cerr << msg.str();
        throw irtkException( msg.str(),
                             __FILE__,
                             __LINE__ );
    }
#else
    stringstream msg;
    msg << "cofstream::Open: Can't write compressed file without zlib" << endl;
    cerr << msg.str();
    throw irtkException( msg.str(),
                         __FILE__,
                         __LINE__ );
#endif
  }
}

inline void irtkCofstream::Close()
{
#ifdef HAS_ZLIB
  if (_compressedFile != NULL) {
    gzclose(_compressedFile);
    _compressedFile = NULL;
  }
#endif
  if (_uncompressedFile != NULL) {
    fclose(_uncompressedFile);
    _uncompressedFile = NULL;
  }
}

inline int irtkCofstream::IsCompressed()
{
  return _compressed;
}

inline void irtkCofstream::IsCompressed(int compressed)
{
  _compressed = compressed;
}

inline int irtkCofstream::IsSwapped()
{
  return _swapped;
}

inline void irtkCofstream::IsSwapped(int swapped)
{
  _swapped = swapped;
}

#endif


