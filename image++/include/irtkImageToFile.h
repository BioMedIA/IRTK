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

#ifndef _IRTKIMAGETOFILE_H

#define _IRTKIMAGETOFILE_H

#include <irtkFileToImage.h>

/**
 * Abstract base class for any general image to file filter.
 *
 * This is the abstract base class which defines a common interface for all
 * filters which take an image as input and produce an image file as output.
 * Each derived class has to implement all of the abstract member functions.
 */

class irtkImageToFile : protected irtkCofstream
{

protected:

  /// Pointer to image for input
  irtkImage *_input;

  // Pointer to filename for output
  const char *_output;

  /** Start address of the data in the image file. Should be
   *  initialized by calling Initialize().
   */
  int _start;

  /// Flag whether to reflect X axis
  int _reflectX;

  /// Flag whether to reflect Y axis
  int _reflectY;

  /// Flag whether to reflect Z axis
  int _reflectZ;

  /// Initialize filter
  virtual void Initialize();

  /// Finalize filter
  virtual void Finalize();

public:

  /// Constructor
  irtkImageToFile();

  /// Destructor
  virtual ~irtkImageToFile();

  /** Static constructor. This constructor allocates a derived class which
   *  can be used to read the image file. This involves checking the file
   *  format of the file and dynamically allocating the corresonding derived
   *  class.
   */
  static irtkImageToFile *New(const char *);

  /// Set input
  virtual void SetInput (irtkImage *);

  /// Set output
  virtual void SetOutput(const char *);

  /// Write entire image
  virtual void Run();

  /// Returns the name of the class
  virtual const char *NameOfClass() = 0;
};

#include <irtkImageToFilePGM.h>
#include <irtkImageToFilePNG.h>
#include <irtkImageToFileVTK.h>
#include <irtkImageToFileGIPL.h>
#include <irtkImageToFileANALYZE.h>
#ifdef HAS_NIFTI
#include <irtkImageToFileNIFTI.h>
#endif

#endif
