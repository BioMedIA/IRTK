/*=========================================================================

  Library   : Image Registration Toolkit (IRTK)
  Module    : $Id$
  Copyright : Imperial College, Department of Computing
              Visual Information Processing (VIP), 2008 onwards
  Date      : $Date$
  Version   : $Revision$
  Changes   : $Author$

=========================================================================*/

#ifndef _IRTKIMAGEFREEFORMREGISTRATION2D_H

#define _IRTKIMAGEFREEFORMREGISTRATION2D_H

/**
 * Filter for non-rigid registration based on voxel similarity measures.
 *
 * This class implements a registration filter for the non-rigid registration
 * of two images. The basic algorithm is described in Rueckert et al., IEEE
 * Transactions on Medical Imaging, 1999.
 *
 */

class irtkImageFreeFormRegistration2D : public irtkImageFreeFormRegistration
{

protected:

  /** Evaluates the registration. This function evaluates the registration by
   *  looping over the target image and interpolating the transformed source
   *  image while filling the joint histogram. This function returns the value
   *  of the similarity measure using Similarity().
   */
  virtual double Evaluate();

  /** Evaluates the registration. This function evaluates the registration by
   *  looping over the target image and interpolating the transformed source
   *  image while filling the joint histogram. This function returns the value
   *  of the similarity measure using Similarity(). This function uses the
   *  cached result of any previous call to Evaluate() and recalculates the
   *  similarity measure in the specified region of interest.
   */
  virtual double EvaluateDerivative(int, double);

public:

  /// Returns the name of the class
  virtual const char *NameOfClass();

  /// Guess parameters
  virtual void GuessParameter();
};

inline const char *irtkImageFreeFormRegistration2D::NameOfClass()
{
  return "irtkFreeFormRegistration2D";
}

#endif
