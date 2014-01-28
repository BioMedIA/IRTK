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

#ifndef _IRTKIMAGETOFILENIFTI_H

#define _IRTKIMAGETOFILENIFTI_H

#ifdef HAS_NIFTI

#include <irtkNIFTI.h>

/**
 * Class for image to NIFTI file writer.
 *
 * This is a class which takes an image as input and produces an image file
 * in NIFTI file format.
 */

class irtkImageToFileNIFTI : public irtkImageToFile
{

  /// Filename of header
  char *_headername;

  /// The nifti header
  irtkNIFTIHeader _hdr;

protected:

  /// Initialize filter
  virtual void Initialize();

  /// Finalize filter (empty, overrides parent method).
  virtual void Finalize();

public:

  /// Constructor
  irtkImageToFileNIFTI();

  /// Destructor
  virtual ~irtkImageToFileNIFTI();

  /// Set input
  virtual void SetOutput(const char *);

  /// Return name of class
  virtual const char *NameOfClass();

  /// Write entire image
  virtual void Run();
};

#endif

#endif
