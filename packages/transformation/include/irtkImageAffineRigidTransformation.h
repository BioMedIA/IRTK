/*=========================================================================

  Library   : Image Registration Toolkit (IRTK)
  Module    : $Id: irtkImageAffineRigidTransformation.h 8 2009-03-02 16:12:58Z dr $
  Copyright : Imperial College, Department of Computing
              Visual Information Processing (VIP), 2008 onwards
  Date      : $Date: 2009-03-02 16:12:58 +0000 (Mon, 02 Mar 2009) $
  Version   : $Revision: 8 $
  Changes   : $Author: dr $

=========================================================================*/

#ifndef _IRTKAFFINERIGIDIMAGETRANSFORMATION_H

#define _IRTKAFFINERIGIDIMAGETRANSFORMATION_H

#include <irtkImage.h>

#include <irtkTransformation.h>

class irtkImageAffineRigidTransformation  : public irtkImageTransformation
{
protected:

  /// rigid and affine transformation matrix
  irtkMatrix _rigid;
  irtkMatrix _affine;

public:

  /** Constructor. This constructs an transformation filter with a given
   *  interpolation mode and padding value. By default the interpolation
   *  mode is set to trilinear.
   */
  irtkImageAffineRigidTransformation();

  /// Destructor
  virtual ~irtkImageAffineRigidTransformation();

  /// Sets transformation
  virtual void SetTransformation(irtkTransformation *);

  /// Runs the filter
  virtual void Run();

  /// Splits the transformation matrix into rigid and affine parameters
  virtual void SplitTransformationMatrix();

  /// Merge the transformation matrizes back together
  virtual void MergeTranformationMatrix();

  /// Transforms the image origin and axis, not using the affine parameters
  virtual void TransformImageCoordSystem();
};

#endif
