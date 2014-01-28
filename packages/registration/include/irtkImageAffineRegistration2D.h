/*=========================================================================

  Library   : Image Registration Toolkit (IRTK)
  Module    : $Id$
  Copyright : Imperial College, Department of Computing
              Visual Information Processing (VIP), 2008 onwards
  Date      : $Date$
  Version   : $Revision$
  Changes   : $Author$

Copyright (c) IXICO LIMITED
All rights reserved.
See COPYRIGHT for details

=========================================================================*/

#ifndef _IRTKIMAGEAFFINEREGISTRATION2D_H

#define _IRTKIMAGEAFFINEREGISTRATION2D_H

/**
 * Filter for affine registration based on voxel similarity measures.
 *
 * This class implements a registration filter for the affine registration
 * of images. The basic algorithm is described in Studholme, Medical Image
 * Analysis, Vol. 1, No. 2, 1996.
 *
 */

class irtkImageAffineRegistration2D : public irtkImageRigidRegistration2D
{

public:

  /** Sets the output for the registration filter. The output must be a affine
   *  transformation. The current parameters of the affine transformation are
   *  used as initial guess for the affine registration. After execution of the
   *  filter the parameters of the affine transformation are updated with the
   *  optimal transformation parameters.
   */
  virtual void SetOutput(irtkTransformation *);

  /// Returns the name of the class
  virtual const char *NameOfClass();

};

inline const char *irtkImageAffineRegistration2D::NameOfClass()
{
  return "irtkImageAffineRegistration2D";
}

#endif
