/*=========================================================================

  Library   : Image Registration Toolkit (IRTK)
  Module    : $Id$
  Copyright : Imperial College, Department of Computing
              Visual Information Processing (VIP), 2008 onwards
  Date      : $Date$
  Version   : $Revision$
  Changes   : $Author$

=========================================================================*/

#ifndef _IRTKIMAGEEIGENFREEFORMREGISTRATION2D_H

#define _IRTKIMAGEEIGENFREEFORMREGISTRATION2D_H

#include <irtkRegistration.h>

/**
 * Filter for adaptive and hierarchical non-rigid registration based
 * on voxel similarity measures and modal reparameterization.
 *
 * This class implements a registration filter for the adaptive hierarchical
 * non-rigid registration of two images by manipulating the shape parameter
 * vector of the eigen system rather than the control points.
 * This algorithm is yet to be published.
 *
 */

class irtkImageEigenFreeFormRegistration2D : public irtkImageFreeFormRegistration2D
{

protected:

  /// Number of modes
  int _NumberOfModes;

  /// Short cut to local transformation (includes eigen vectors and values)
  irtkEigenFreeFormTransformation *_effd;

  /// Initial set up for the registration at a multiresolution level
  virtual void Initialize(int);

  /// Final set up for the registration at a multiresolution level
  virtual void Finalize(int);

public:

  //
  // Constructor and destructor
  //

  /// Constructor
  irtkImageEigenFreeFormRegistration2D();

  /// Destructor
  virtual ~irtkImageEigenFreeFormRegistration2D();

  // Access parameters
  virtual SetMacro(NumberOfModes, int);
  virtual GetMacro(NumberOfModes, int);

  //
  // Inherited methods to be implemented
  //

  /// Returns the name of the class
  virtual const char *NameOfClass();

  /** Evaluates the gradient of the similarity metric. This function
   *  evaluates the gradient of the similarity metric of the registration
   *  by looping over the target image and interpolating the transformed
   *  source image while filling the joint histogram. The partial derivatives
   *  are approximated using a finite difference scheme. The step size for the
   *  finite difference scheme is passed as a parameter to the
   *  function. The function returns the norm of the gradient vector as
   *  well as the gradient vector containing the partial derivatives.
   */
  virtual double EvaluateGradient(float, float *);

};

//
// Inherited methods to be implemented
//

inline const char *irtkImageEigenFreeFormRegistration2D::NameOfClass()
{
  return "irtkImageEigenFreeFormRegistration2D";
}

#endif
