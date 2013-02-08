/*=========================================================================

  Library   : Image Registration Toolkit (IRTK)
  Copyright : Imperial College, Department of Computing
              Visual Information Processing (VIP), 2008 onwards
  Date      : $Date: 2013-01-14 $
  Changes   : $Author: cl6311 $

=========================================================================*/

#ifndef _IRTKMEDIAN_H

#define _IRTKMEDIAN_H

#include <irtkImageToImage.h>

/**
 * Class for median filtering an image
 *
 */

template <class VoxelType> class irtkMedianFilter : public irtkImageToImage<VoxelType>
{

protected:

  /// Returns whether the filter requires buffering
  virtual bool RequiresBuffering();

  /// Returns the name of the class
  virtual const char *NameOfClass();

  /// Initialize the filter
  virtual void Initialize();

  /// What connectivity to assume when running the filter.
  int _kernelRadius;

  irtkRealImage* _mask;

public:

  /// Constructor
  irtkMedianFilter();

  /// Destructor
  ~irtkMedianFilter();

  /// Run dilation
  virtual void Run();

  /// Set input image for filter
  void SetMask (irtkRealImage*);

  SetMacro(kernelRadius, int);

  GetMacro(kernelRadius, int);
};

#endif
