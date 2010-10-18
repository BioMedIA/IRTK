/*=========================================================================

  Library   : Image Registration Toolkit (IRTK)
  Module    : $Id$
  Copyright : Imperial College, Department of Computing
              Visual Information Processing (VIP), 2008 onwards
  Date      : $Date$
  Version   : $Revision$
  Changes   : $Author$

=========================================================================*/

#ifndef _IRTKIMAGETRANSFORMATION2_H

#define _IRTKIMAGETRANSFORMATION2_H

#include <irtkImage.h>

#include <irtkImageFunction.h>

#include <irtkTransformation.h>

/**
 * Filter for image transformations.
 *
 * This class implements an image filter which takes an input image and a
 * transformation and computes the corresponding output image. The output
 * image is calculated by looping over the voxel locations and calculating
 * the corresponding voxel locations in the input image. The intensity of
 * the voxels of the output image is by interpolation from the input image.
 * Note, that the transformation is assumed to map the input image into the
 * ouput image and is therefore inverted during the execution of the filter.
 * All calculations are carried out using world coordinates rather than image
 * coordinates.
 *
 */

class irtkImageTransformation2
{

protected:

  /// Input for the image to image filter
  irtkImage *_input;

  /// Output for the image to image filter
  irtkImage *_output;

  /// First transformation
  irtkTransformation *_transformation;

  /// Second transformation
  irtkTransformation *_transformation2;

  /// Interpolation
  irtkImageFunction *_interpolator;

  /// Padding value in target (voxels in the target image with this
  /// value will be ignored)
  double _TargetPaddingValue;

  /// Padding value in source (voxels outside the source image will
  /// be set to this value)
  double _SourcePaddingValue;

  /// Scale factor for intensities in transformed image
  double _ScaleFactor;
  
  /// Offset for intensities in transformed image
  double _Offset;

  /// Flag whether to invert transformation
  int _Invert;

public:

  /** Constructor. This constructs an transformation filter with a given
   *  interpolation mode and padding value. By default the interpolation
   *  mode is set to trilinear and invert is off
   */
  irtkImageTransformation2();

  /// Destructor
  virtual ~irtkImageTransformation2();

  /// Sets input image
  virtual void SetInput (irtkImage *);

  /// Sets output image
  virtual void SetOutput(irtkImage *);

  /// Sets transformation
  virtual void SetTransformation(irtkTransformation *, irtkTransformation *);

  /// Gets the target padding value
  virtual double GetTargetPaddingValue(void);

  /// Puts the target padding value
  virtual void PutTargetPaddingValue(double);

  /// Gets the source padding value
  virtual double GetSourcePaddingValue(void);

  /// Puts the source padding value
  virtual void PutSourcePaddingValue(double);

  /// Gets the interpolator
  virtual irtkImageFunction *GetInterpolator(void);

  /// Sets the interpolator
  virtual void PutInterpolator(irtkImageFunction *);

  /// Invert on
  virtual void InvertOn(void);

  /// Sets the interpolator
  virtual void InvertOff(void);

  /// Runs the filter
  virtual void Run();

};

inline double irtkImageTransformation2::GetTargetPaddingValue()
{
  return _TargetPaddingValue;
}

inline void irtkImageTransformation2::PutTargetPaddingValue(double PaddingValue)
{
  _TargetPaddingValue = PaddingValue;
}

inline double irtkImageTransformation2::GetSourcePaddingValue()
{
  return _SourcePaddingValue;
}

inline void irtkImageTransformation2::PutSourcePaddingValue(double PaddingValue)
{
  _SourcePaddingValue = PaddingValue;
}

inline irtkImageFunction *irtkImageTransformation2::GetInterpolator()
{
  return _interpolator;
}

inline void irtkImageTransformation2::PutInterpolator(irtkImageFunction *interpolator)
{
  _interpolator = interpolator;
}

inline void irtkImageTransformation2::InvertOn()
{
  _Invert = true;
}

inline void irtkImageTransformation2::InvertOff()
{
  _Invert = false;
}

#endif
