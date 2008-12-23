/*=========================================================================

  Library   : Image Registration Toolkit (IRTK)
  Module    : $Id$
  Copyright : Imperial College, Department of Computing
              Visual Information Processing (VIP), 2008 onwards
  Date      : $Date$
  Version   : $Revision$
  Changes   : $Author$

=========================================================================*/

#ifndef _IRTKFILEVTKTOIMAGE_H

#define _IRTKFILEVTKTOIMAGE_H

/**
 * Class for reading images in VTK file format.
 *
 * This is a class which reads images in VTK file format and converts them
 * into images. The VTK file format is a file format for 3D images. It is
 * used by the Visualization Toolkit to store structered points. Supported
 * voxel types are char, unsigned char, short and unsigned short.
 */

template <class VoxelType> class irtkFileVTKToImage : public irtkFileToImage<VoxelType>
{

protected:

  /// Read header of VTK file
  virtual void ReadHeader();

public:

  /// Returns name of class
  virtual const char *NameOfClass();

  /// Returns whether file has correct header
  static int CheckHeader(const char *);
};

#endif
