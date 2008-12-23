/*=========================================================================

  Library   : Image Registration Toolkit (IRTK)
  Module    : $Id$
  Copyright : Imperial College, Department of Computing
              Visual Information Processing (VIP), 2008 onwards
  Date      : $Date$
  Version   : $Revision$
  Changes   : $Author$

=========================================================================*/

#ifndef _IRTKIMAGETOFILEPNG_H

#define _IRTKIMAGETOFILEPNG_H

#ifdef HAS_PNG

/**
 * Class for image to PNG file filter.
 *
 * This is a class which takes an image as input and produces an image file
 * in PNG file format. Note that PNG file formats support only 2D images!!!
 */

template <class VoxelType> class irtkImageToFilePNG : public irtkImageToFile<VoxelType>
{

protected:

  /// Initialize filter
  virtual void Initialize();

public:

  /// Write entire image
  virtual void Run();

  /// Return name of class
  virtual const char *NameOfClass();

};

#endif

#endif
