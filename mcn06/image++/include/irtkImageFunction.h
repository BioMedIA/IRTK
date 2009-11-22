/*=========================================================================

  Library   : Image Registration Toolkit (IRTK)
  Module    : $Id$
  Copyright : Imperial College, Department of Computing
              Visual Information Processing (VIP), 2008 onwards
  Date      : $Date$
  Version   : $Revision$
  Changes   : $Author$

=========================================================================*/

#ifndef _IRTKIMAGEFUNCTION_H

#define _IRTKIMAGEFUNCTION_H

/**
 * Abstract base class for any general image interpolation function filter.
 *
 * This is the abstract base class which defines a common interface for all
 * filters which take an image as input and sample that image at arbitrary
 * location. Each derived class has to implement all abstract member functions.
 */

class irtkImageFunction : public irtkObject
{

private:

  /// Debugging flag
  Bool _DebugFlag;

protected:

  /// Input image for filter
  irtkImage *_input;

  /// Default value to return
  double _DefaultValue;

public:

  /// Constructor
  irtkImageFunction();

  /// Deconstuctor
  virtual ~irtkImageFunction();

  /// Set input image for filter
  virtual void SetInput (irtkImage *);

  /** Initialize the filter. This function must be called by any derived
   *  filter class to perform some initialize tasks. */
  virtual void Initialize();

  /// Evaluate the filter at an arbitrary image location (in pixels)
  virtual double Evaluate(double, double, double, double = 0) = 0;

  /// Returns the name of the class
  virtual const char *NameOfClass() = 0;

  /// Set debugging flag
  SetMacro(DebugFlag, Bool);

  /// Get debugging flag
  GetMacro(DebugFlag, Bool);

  /// Print debugging messages if debugging is enabled
  virtual void Debug(const char *);
};

#include <irtkInterpolateImageFunction.h>

#endif
