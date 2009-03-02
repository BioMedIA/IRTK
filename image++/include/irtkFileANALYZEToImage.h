/*=========================================================================

  Library   : Image Registration Toolkit (IRTK)
  Module    : $Id$
  Copyright : Imperial College, Department of Computing
              Visual Information Processing (VIP), 2008 onwards
  Date      : $Date$
  Version   : $Revision$
  Changes   : $Author$

=========================================================================*/

#ifndef _IRTKFILEANALYZETOIMAGE_H

#define _IRTKFILEANALYZETOIMAGE_H

/**
 * Class for reading images in ANALYZE file format.
 *
 * This is a class which reads images in ANALYZE file format and converts them
 * into images. The ANALYZE file format is a file format for 3D images. It has
 * been developed by the Mayo Clinic, Rochester, MN. The header information
 * is stored in a file with extension ".hdr" while the image data is stored in
 * a file with with extension ".img".
 */

class irtkFileANALYZEToImage : public irtkFileToImage
{

  /// Filename of header
  char *_headername;

protected:

  /// Read header of ANALYZE file
  virtual void ReadHeader();

public:

  /// Contructor
  irtkFileANALYZEToImage();

  /// Destructor
  virtual ~irtkFileANALYZEToImage();

  /// Set input
  virtual void SetInput (const char *);

  /// Returns name of class
  virtual const char *NameOfClass();

  /// Returns whether file has correct header
  static int CheckHeader(const char *);
};

#endif
