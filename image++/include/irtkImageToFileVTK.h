/*=========================================================================

  Library   : Image Registration Toolkit (IRTK)
  Module    : $Id$
  Copyright : Imperial College, Department of Computing
              Visual Information Processing (VIP), 2008 onwards
  Date      : $Date$
  Version   : $Revision$
  Changes   : $Author$

=========================================================================*/

#ifndef _IRTKIMAGETOFILEVTK_H

#define _IRTKIMAGETOFILEVTK_H

/**
 * Class for image to VTK file filter.
 *
 * This is a class which takes an image as input and produces an image file
 * in VTK file format.
 */

template <class VoxelType> class irtkImageToFileVTK : public irtkImageToFile<VoxelType>
{

protected:

  /// Initialize filter
  virtual void Initialize();

public:

  /// Return name of class
  virtual const char *NameOfClass();

};

#endif
