/*=========================================================================

  Library   : Image Registration Toolkit (IRTK)
  Module    : $Id$
  Copyright : Imperial College, Department of Computing
              Visual Information Processing (VIP), 2008 onwards
  Date      : $Date$
  Version   : $Revision$
  Changes   : $Author$

=========================================================================*/

#ifndef _IRTKFILETOIMAGE_H

#define _IRTKFILETOIMAGE_H

/**
 * Abstract base class for any general file to image filter.
 *
 * This is the abstract base class which defines a common interface for all
 * filters which take a filename (referrencing an image file) as input and
 * produce an image as output. For each file format a derived class should
 * be created. Each derived class has to implement abstract member functions
 * ReadHeader() and CheckHeader().
 */

class irtkFileToImage : protected irtkCifstream
{

  /// File name of image file
  char *_imagename;

protected:

	/// Image attributes
	irtkImageAttributes _attr;
	
  /// No. of bytes per voxel. Should be initialized by calling ReadHeader().
  int _bytes;

  /// Type of voxels. Should be initialized by calling ReadHeader().
  int _type;

  /// Intensity scaling parameter - slope (default: 1)
  float _slope;

  /// Intensity scaling parameter -  intercept (default: 0)
  float _intercept;

  /// Start of image data
  int _start;

  /// Flag whether to reflect X axis
  int _reflectX;

  /// Flag whether to reflect Y axis
  int _reflectY;

  /// Flag whether to reflect Z axis
  int _reflectZ;

  /// Debug flag
  int _debug;

  /** Read header. This is an abstract function. Each derived class has to
   *  implement this function in order to initialize image dimensions, voxel
   *  dimensions, voxel type and a lookup table which the address for each
   *  slice.
   */
  virtual void ReadHeader() = 0;

public:

  /// Contructor
  irtkFileToImage();

  /// Destructor
  virtual ~irtkFileToImage();

  /** Static constructor. This constructor allocates a derived class which
   *  is then  used to read the image file. This involves checking the file
   *  format of the file and dynamically allocating the corresonding derived
   *  class.
   */
  static irtkFileToImage *New(const char *);

  /// Set input
  virtual void SetInput (const char *);

  /// Get output
  virtual irtkImage *GetOutput();

  /// Get debug flag
  virtual int  GetDebugFlag();

  /// Get slope
  virtual double  GetSlope();

  /// Get intercept
  virtual double  GetIntercept();

  /// Put debug flag
  virtual void PutDebugFlag(int);

  /// Print image file information
  virtual void Print();

  /// Print  debugging info
  virtual void Debug(char *);

  // Returns the name of the class
  virtual const char *NameOfClass() = 0;
};

#include <irtkFilePGMToImage.h>
#include <irtkFileVTKToImage.h>
#include <irtkFileGIPLToImage.h>
#ifdef HAS_NIFTI
#include <irtkFileNIFTIToImage.h>
#endif
#include <irtkFileANALYZEToImage.h>

#endif
