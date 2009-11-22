/*=========================================================================

  Library   : Image Registration Toolkit (IRTK)
  Module    : $Id$
  Copyright : Imperial College, Department of Computing
              Visual Information Processing (VIP), 2008 onwards
  Date      : $Date$
  Version   : $Revision$
  Changes   : $Author$

=========================================================================*/

#ifndef _IRTKIMAGETRANSFORMATION_H

#define _IRTKIMAGETRANSFORMATION_H

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

class irtkImageTransformation
{

public:

  /// Input for the image to image filter
  irtkImage *_input;

  /// Output for the image to image filter
  irtkImage *_output;

  /// Transformation
  irtkTransformation *_transformation;

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
  irtkImageTransformation();

  /** Static constructor. This constructs an transformation filter for a
   *  given transformation. If there is a fast transformation filter, it
   *  will use this fast transformation filter, otherwise the standard
   *  transformation filter is used.
   */
  static irtkImageTransformation *New(irtkTransformation *);

  /// Destructor
  virtual ~irtkImageTransformation();

  /// Sets input image
  virtual void SetInput (irtkImage *);

  /// Sets input image and transformation
  virtual void SetInput (irtkImage *, irtkTransformation *);

  /// Sets output image
  virtual void SetOutput(irtkImage *);

  /// Sets transformation
  virtual void SetTransformation(irtkTransformation *);

  /// Puts the target padding value
  virtual void PutTargetPaddingValue(double);

  /// Gets the target padding value
  virtual double GetTargetPaddingValue();

  /// Puts the source padding value
  virtual void PutSourcePaddingValue(double);

  /// Gets the source padding value
  virtual double GetSourcePaddingValue();

  /// Gets the interpolator
  virtual irtkImageFunction *GetInterpolator(void);

  /// Sets the interpolator
  virtual void PutInterpolator(irtkImageFunction *);

  /// Sets scale and otfset
  virtual void PutScaleFactorAndOffset(double = 1, double = 0);
  
  /// Invert on
  virtual void InvertOn(void);

  /// Sets the interpolator
  virtual void InvertOff(void);

  /// Runs the filter
  virtual void Run();

};

inline void irtkImageTransformation::PutTargetPaddingValue(double PaddingValue)
{
  _TargetPaddingValue = PaddingValue;
}

inline double irtkImageTransformation::GetTargetPaddingValue()
{
  return _TargetPaddingValue;
}

inline void irtkImageTransformation::PutSourcePaddingValue(double PaddingValue)
{
  _SourcePaddingValue = PaddingValue;
}

inline double irtkImageTransformation::GetSourcePaddingValue()
{
  return _SourcePaddingValue;
}

inline irtkImageFunction *irtkImageTransformation::GetInterpolator()
{
  return _interpolator;
}

inline void irtkImageTransformation::PutInterpolator(irtkImageFunction *interpolator)
{
  _interpolator = interpolator;
}

inline void irtkImageTransformation::PutScaleFactorAndOffset(double ScaleFactor, double Offset)
{
	_ScaleFactor = ScaleFactor;
	_Offset      = Offset;
}

inline void irtkImageTransformation::InvertOn()
{
  _Invert = True;
}

inline void irtkImageTransformation::InvertOff()
{
  _Invert = False;
}

#endif
