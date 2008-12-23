/*=========================================================================

  Library   : Image Registration Toolkit (IRTK)
  Module    : $Id$
  Copyright : Imperial College, Department of Computing
              Visual Information Processing (VIP), 2008 onwards
  Date      : $Date$
  Version   : $Revision$
  Changes   : $Author$

=========================================================================*/

#ifndef _IRTKFILENIFTITOIMAGE_H

#define _IRTKFILENIFTITOIMAGE_H

#ifdef HAS_NIFTI

/**
 * Class for reading images in NIFTI file format.
 *
 * This is a class which reads images in NIFTI file format and converts them
 * into images. The NIFTI file format is a file format for 3D and 4D images.
 * More information about the format can be found at http://nifti.nimh.nih.gov/nifti-1/
 */

template <class VoxelType> class irtkFileNIFTIToImage : public irtkFileToImage<VoxelType>
{

  /// Filename of header
  char *_headername;

protected:

  /// Read header of NIFTI file
  virtual void ReadHeader();

public:

  /// Contructor
  irtkFileNIFTIToImage();

  /// Destructor
  virtual ~irtkFileNIFTIToImage();

  /// Set input
  virtual void SetInput (const char *);

  /// Returns name of class
  virtual const char *NameOfClass();

  /// Returns whether file has correct header
  static int CheckHeader(const char *);

  /// Print image file information
  virtual void Print();


};

#endif

#endif
