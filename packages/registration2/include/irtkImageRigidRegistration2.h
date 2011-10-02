/*=========================================================================

  Library   : Image Registration Toolkit (IRTK)
  Module    : $Id$
  Copyright : Imperial College, Department of Computing
              Visual Information Processing (VIP), 2009 onwards
  Date      : $Date$
  Version   : $Revision$
  Changes   : $Author$

=========================================================================*/

#ifndef _IRTKIMAGERIGIDREGISTRATION2_H

#define _IRTKIMAGERIGIDREGISTRATION2_H

/**
 * Filter for rigid registration based on voxel similarity measures.
 *
 * This class implements a registration filter for the rigid registration of
 * two images. The basic algorithm is described in Studholme, Medical Image
 * Analysis, Vol. 1, No. 2, 1996.
 *
 */

class irtkImageRigidRegistration2 : public irtkImageRegistration2
{

protected:

  /// Update state of the registration based on current transformation estimate (source image)
  virtual void UpdateSource();

  /// Update state of the registration based on current transformation estimate (source image and source image gradient)
  virtual void UpdateSourceAndGradient();

  /// Evaluate the gradient of the similarity measure for the current transformation.
  virtual double EvaluateGradient(double *);

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

  /// Print information about the progress of the registration
  virtual void Print();

  /// Guess parameters
  virtual void GuessParameter();
};

inline void irtkImageRigidRegistration2::SetOutput(irtkTransformation *transformation)
{
  if (strcmp(transformation->NameOfClass(), "irtkRigidTransformation") != 0) {
    cerr << "irtkImageRigidRegistration2::SetOutput: Transformation must be rigid"
         << endl;
    exit(0);
  }
  _transformation = transformation;
}

inline const char *irtkImageRigidRegistration2::NameOfClass()
{
  return "irtkImageRigidRegistration2";
}

inline void irtkImageRigidRegistration2::Print()
{
  _transformation->Print();
}

#endif
