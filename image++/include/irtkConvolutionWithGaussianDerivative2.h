/*=========================================================================

  Library   : Image Registration Toolkit (IRTK)
  Module    : $Id$
  Copyright : Imperial College, Department of Computing
              Visual Information Processing (VIP), 2008 onwards
  Date      : $Date$
  Version   : $Revision$
  Changes   : $Author$

=========================================================================*/

#ifndef _IRTKCONVOLUTIONWITHGAUSSIANDERIVATIVE2_H

#define _IRTKCONVOLUTIONWITHGAUSSIANDERIVATIVE2_H

#include <irtkImageToImage.h>

/** 
 * Class for convolution with 2nd order Gaussian derivative 
 * 
 * This class defines and implements the 2nd order gaussian derivative filtering of images. 
 */

template <class VoxelType> class irtkConvolutionWithGaussianDerivative2 : public irtkImageToImage<VoxelType> {

protected:

  /// Sigma (standard deviation of Gaussian kernel)
  double _Sigma;

  /// Returns the name of the class
  const char *NameOfClass();

  /// Returns whether the class requires buffer (true)
  virtual bool RequiresBuffering();

public:

  /// Constructor
  irtkConvolutionWithGaussianDerivative2(double);

  /// Destructor
  ~irtkConvolutionWithGaussianDerivative2();

  /// Compute derivatives
  void Ixx();

  /// Compute derivatives
  void Iyy();

  /// Compute derivatives
  void Izz();

  /// Compute derivatives
  void Ixy();

  /// Compute derivatives
  void Ixz();

  /// Compute derivatives
  void Iyz();
  
  /// Set sigma
  SetMacro(Sigma, double);
 
  /// Get sigma
  GetMacro(Sigma, double);

};


#endif
