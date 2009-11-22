/*=========================================================================

  Library   : Image Registration Toolkit (IRTK)
  Module    : $Id$
  Copyright : Imperial College, Department of Computing
              Visual Information Processing (VIP), 2008 onwards
  Date      : $Date$
  Version   : $Revision$
  Changes   : $Author$

=========================================================================*/

#ifndef _IRTKCONVOLUTIONWITHGAUSSIANDERIVATIVE_H

#define _IRTKCONVOLUTIONWITHGAUSSIANDERIVATIVE_H

#include <irtkImageToImage.h>

/** 
 * Class for convolution with 1st order Gaussian derivative 
 * 
 * This class defines and implements the 1st order gaussian derivative filtering of images. 
 */

template <class VoxelType> class irtkConvolutionWithGaussianDerivative : public irtkImageToImage<VoxelType> {

protected:

  /// Sigma (standard deviation of Gaussian kernel)
  double _Sigma;

  /// Returns the name of the class
  const char *NameOfClass();

  /// Returns whether the class requires buffer (true)
  virtual Bool RequiresBuffering();

public:

  /// Constructor
  irtkConvolutionWithGaussianDerivative(double);

  /// Destructor
  ~irtkConvolutionWithGaussianDerivative();

  /// Compute derivatives
  void Ix();

  /// Compute derivatives
  void Iy();

  /// Compute derivatives
  void Iz();
  
  /// Set sigma
  SetMacro(Sigma, double);
 
  /// Get sigma
  GetMacro(Sigma, double);

};


#endif
