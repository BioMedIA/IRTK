/*=========================================================================

  Library   : Image Registration Toolkit (IRTK)
  Module    : $Id$
  Copyright : Imperial College, Department of Computing
              Visual Information Processing (VIP), 2008 onwards
  Date      : $Date$
  Version   : $Revision$
  Changes   : $Author$

=========================================================================*/

#ifndef _IRTKIMAGETOIMAGE_H

#define _IRTKIMAGETOIMAGE_H

#ifdef HAS_TBB

template <class VoxelType> class irtkMultiThreadedImageToImage;

#endif

/**
 * Abstract base class for any general image to image filter.
 *
 * This is the abstract base class which defines a common interface for all
 * filters which take an image as input and produce an image as output. Each
 * derived class has to implement all abstract member functions.
 */

template <class VoxelType> class irtkImageToImage : public irtkObject
{

#ifdef HAS_TBB

  friend class irtkMultiThreadedImageToImage<VoxelType>;

#endif

private:

  /// Debugging flag
  bool _DebugFlag;

protected:

  /// Buffer
  irtkGenericImage<VoxelType> *_tmp;

  /// Input image for filter
  irtkGenericImage<VoxelType> *_input;

  /// Output image for filter
  irtkGenericImage<VoxelType> *_output;

  /** Initialize the filter. This function must be called by any derived
   *  filter class to perform some initialize tasks. */
  virtual void Initialize();

  /** Finalize the filter. This function must be called by any derived
   *  filter class to perform some initialize tasks. */
  virtual void Finalize();

public:

  /// Constructor
  irtkImageToImage();

  /// Deconstuctor
  virtual ~irtkImageToImage();

  /// Set input image for filter
  virtual void SetInput (irtkGenericImage<VoxelType> *);

  /// Set output image for filter
  virtual void SetOutput(irtkGenericImage<VoxelType> *);

  /// Run filter on entire image
  virtual void   Run();

  /// Run filter on single voxel
  virtual double Run(int, int, int, int = 0);


  /** Returns whether the filter requires buffering. Any derived class must
   *  implement this member function to indicate whether the filter should
   *  buffer the input in case that input and output are equal. For example,
   *  filters which only require the voxel value to calculate their output
   *  should return false, otherwise true.
   */
  virtual bool RequiresBuffering() = 0;

  /// Returns the name of the class
  virtual const char *NameOfClass() = 0;

  /// Set debugging flag
  SetMacro(DebugFlag, bool);

  /// Get debugging flag
  GetMacro(DebugFlag, bool);

  /// Print debugging messages if debugging is enabled
  virtual void Debug(const char *);
};

#endif
