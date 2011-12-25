/*=========================================================================

  Library   : Image Registration Toolkit (IRTK)
  Module    : $Id$
  Copyright : Imperial College, Department of Computing
              Visual Information Processing (VIP), 2008 onwards
  Date      : $Date$
  Version   : $Revision$
  Changes   : $Author$

=========================================================================*/

#ifndef _IRTKIMAGETOOPENCV_H

#define _IRTKIMAGETOOPENCV_H

#ifdef HAS_OPENCV
#include <cv.h>

/**
 * Abstract base class for any general image to image filter.
 *
 * This is the abstract base class which defines a common interface for all
 * filters which take an image as input and produce an image as output. Each
 * derived class has to implement all abstract member functions.
 */

template <class VoxelType> class irtkImageToOpenCv : public irtkObject
{

private:

  /// Buffer
  irtkGenericImage<VoxelType> *_tmp;

protected:


  /// Input image for filter
  irtkGenericImage<VoxelType> *_input;

  /// Output image for filter
  IplImage *_output;

  /// min and max value
  VoxelType min,max;

  /** Initialize the filter. This function must be called by any derived
   *  filter class to perform some initialize tasks. */
  virtual void Initialize();

  /** Finalize the filter. This function must be called by any derived
   *  filter class to perform some initialize tasks. */
  virtual void Finalize();

public:

  /// Constructor
  irtkImageToOpenCv();

  /// Deconstuctor
  virtual ~irtkImageToOpenCv();

  /// Set input image for filter
  virtual void SetInput (irtkGenericImage<VoxelType> *);

  /// Set output image for filter
  virtual void SetOutput(IplImage *);

   /// Set input image for filter
  virtual irtkGenericImage<VoxelType> * GetInput ();

  /// Set output image for filter
  virtual IplImage* GetOutput();

  /// Generate output based on input
  virtual void GenOutput ();

  /// Generate input based on output
  virtual void GenInput ();

  /// Run filter on entire image from irtk image to OpenCV image
  virtual void   Run(int = 0);

  /// Run filter on entire image from OpenCV image to irtk image
  virtual void   Invert(int = 0);

  /// Returns the name of the class
  virtual const char *NameOfClass();
};

template <class VoxelType> const char *irtkImageToOpenCv<VoxelType>::NameOfClass()
{
  return "irtkImageToOpenCv";
}

#endif

#endif
