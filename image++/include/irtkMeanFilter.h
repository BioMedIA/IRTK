/*=========================================================================

  Library   : Image Registration Toolkit (IRTK)
  Copyright : Imperial College, Department of Computing
              Visual Information Processing (VIP), 2008 onwards
  Date      : $Date$
  Changes   : $Author$

=========================================================================*/

#ifndef _IRTKMEANFILTER_H

#define _IRTKMEANFILTER_H

#include <irtkImageToImage.h>

/**
 * Class for median filtering an image
 *
 */

template <class VoxelType> class irtkMeanFilter : public irtkImageToImage<VoxelType>
{

protected:

  /// Returns whether the filter requires buffering
  virtual bool RequiresBuffering();

  /// Returns the name of the class
  virtual const char *NameOfClass();

  /// Initialize the filter
  virtual void Initialize();

  /// What connectivity to assume when running the filter.
  int _radius_x;
  int _radius_y;
  int _radius_z;

public:

  /// Constructor
  irtkMeanFilter();

  /// Destructor
  ~irtkMeanFilter();

  /// Run dilation
  virtual void Run();

  virtual void SetkernelRadius(double radius);
};

#endif
