/*=========================================================================

  Library   : Image Registration Toolkit (IRTK)
  Module    : $Id$
  Copyright : Imperial College, Department of Computing
              Visual Information Processing (VIP), 2009 onwards
  Date      : $Date$
  Version   : $Revision$
  Changes   : $Author$

=========================================================================*/

#ifndef _IRTKIMAGEAFFINEREGISTRATION2_H

#define _IRTKIMAGEAFFINEREGISTRATION2_H

/**
 * Filter for affine registration based on voxel similarity measures.
 *
 * This class implements a registration filter for the rigid registration of
 * two images. The basic algorithm is described in Studholme, Medical Image
 * Analysis, Vol. 1, No. 2, 1996.
 *
 */

class irtkImageAffineRegistration2 : public irtkImageRigidRegistration2
{

public:

  /** Sets the output for the registration filter. The output must be a rigid
   *  transformation. The current parameters of the rigid transformation are
   *  used as initial guess for the rigid registration. After execution of the
   *  filter the parameters of the rigid transformation are updated with the
   *  optimal transformation parameters.
   */
  virtual void SetOutput(irtkTransformation *);

  /// Returns the name of the class
  virtual const char *NameOfClass();

};

inline void irtkImageAffineRegistration2::SetOutput(irtkTransformation *transformation)
{
  if (strcmp(transformation->NameOfClass(), "irtkAffineTransformation") != 0) {
    cerr << "irtkImageAffineRegistration2::SetOutput: Transformation must be affine"
         << endl;
    exit(0);
  }
  _transformation = transformation;
}

inline const char *irtkImageAffineRegistration2::NameOfClass()
{
  return "irtkImageAffineRegistration2";
}

#endif
