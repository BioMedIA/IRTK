/*=========================================================================

  Library   : Image Registration Toolkit (IRTK)
  Module    : $Id$
  Copyright : Imperial College, Department of Computing
              Visual Information Processing (VIP), 2008 onwards
  Date      : $Date$
  Version   : $Revision$
  Changes   : $Author$

=========================================================================*/

#ifndef _IRTKIMAGETOFILEGIPL_H

#define _IRTKIMAGETOFILEGIPL_H

/**
 * Class for image to GIPL file filter.
 *
 * This is a class which takes an image as input and produces an image file
 * in GIPL file format.
 */

template <class VoxelType> class irtkImageToFileGIPL : public irtkImageToFile<VoxelType>
{

protected:

  /// Initialize filter
  virtual void Initialize();

public:

  /// Return name of class
  virtual const char *NameOfClass();

};

#endif
