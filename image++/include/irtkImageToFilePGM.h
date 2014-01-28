/*=========================================================================

  Library   : Image Registration Toolkit (IRTK)
  Module    : $Id$
  Copyright : Imperial College, Department of Computing
              Visual Information Processing (VIP), 2008 onwards
  Date      : $Date$
  Version   : $Revision$
  Changes   : $Author$

Copyright (c) 1999-2014 and onwards, Imperial College London
All rights reserved.
See LICENSE for details

=========================================================================*/

#ifndef _IRTKIMAGETOFILEPGM_H

#define _IRTKIMAGETOFILEPGM_H

/**
 * Class for image to PGM file filter.
 *
 * This is a class which takes an image as input and produces an image file
 * in PGM file format. Note that PGM file formats support only 2D images!!!
 */

class irtkImageToFilePGM : public irtkImageToFile
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
