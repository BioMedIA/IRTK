/*=========================================================================

  Library   : Image Registration Toolkit (IRTK)
  Module    : $Id$
  Copyright : Imperial College, Department of Computing
              Visual Information Processing (VIP), 2008 onwards
  Date      : $Date$
  Version   : $Revision$
  Changes   : $Author$

=========================================================================*/

#ifndef _IRTKFLUIDFREEFORMTRANSFORMATION_H

#define _IRTKFLUIDFREEFORMTRANSFORMATION_H

#include <irtkGeometry.h>

/**
 * Class for fluid transformations
 *
 * This class implements 3D fluid transformations
 *
 */

class irtkFluidFreeFormTransformation : public irtkMultiLevelFreeFormTransformation
{

public:

  /// Constructor (default)
  irtkFluidFreeFormTransformation();

  /// Constructor (copy)
  irtkFluidFreeFormTransformation(const irtkRigidTransformation &);

  /// Constructor (copy)
  irtkFluidFreeFormTransformation(const irtkAffineTransformation &);

  /// Constructor (copy)
  irtkFluidFreeFormTransformation(const irtkFluidFreeFormTransformation &);

  /// Combine local transformations
  virtual void CombineLocalTransformation();

  /// Transforms a point
  virtual void Transform(double &, double &, double &, double = 0);

  /// Calculates displacement using global and local transformation components
  virtual void Displacement(double& x, double& y, double& z, double = 0);

  /// Transforms a single point using the global transformation component only
  virtual void GlobalTransform(double &, double &, double &, double = 0);

  /// Transforms a single point using the local transformation component only
  virtual void LocalTransform (double &, double &, double &, double = 0);

  /// Calculates displacement using the global transformation component only
  virtual void GlobalDisplacement(double &, double &, double &, double = 0);

  /// Calculates displacement using the local transformation component only
  virtual void LocalDisplacement(double &, double &, double &, double = 0);

  /// Transforms a point
  virtual void Transform(int, double &, double &, double &, double = 0);

  /// Transforms a single point using the local transformation component only
  virtual void LocalTransform (int, double &, double &, double &, double = 0);

  /// Calculates displacement using the local transformation component only
  virtual void LocalDisplacement(int, double &, double &, double &, double = 0);

  /** Convert the global transformation from a matrix representation to a
      FFD and incorporate it with any existing local displacement. **/
  virtual void MergeGlobalIntoLocalDisplacement();

  /// Calculate the Jacobian of the transformation
  virtual void Jacobian(irtkMatrix &, double, double, double, double = 0);

  /// Calculate the Jacobian of the local transformation
  virtual void LocalJacobian(irtkMatrix &, double, double, double, double = 0);

  /// Inverts the transformation
  virtual double Inverse(double &, double &, double &, double = 0, double = 0.01);

  /// Returns a string with the name of the instantiated class
  virtual const char *NameOfClass();

  /// Reads a transformation from a file
  virtual irtkCifstream& Read(irtkCifstream&);

  /// Writes a transformation to a file
  virtual irtkCofstream& Write(irtkCofstream&);
};

inline const char *irtkFluidFreeFormTransformation::NameOfClass()
{
  return "irtkFluidFreeFormTransformation";

}

#endif
