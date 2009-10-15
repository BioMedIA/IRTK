/*=========================================================================

  Library   : Image Registration Toolkit (IRTK)
  Module    : $Id$
  Copyright : Imperial College, Department of Computing
              Visual Information Processing (VIP), 2008 onwards
  Date      : $Date$
  Version   : $Revision$
  Changes   : $Author$

=========================================================================*/

#ifndef _IRTKFREEFORMTRANSFORMATION_H

#define _IRTKFREEFORMTRANSFORMATION_H

#include <irtkGeometry.h>

/**
 * Class for free form transformations
 *
 * This class implements a free form transformation
 *
 */

class irtkFreeFormTransformation : public irtkTransformation
{

public:

  /// Destructor
  virtual ~irtkFreeFormTransformation();

  /// Returns the of control points in x
  virtual int GetX() const = 0;

  /// Returns the of control points in y
  virtual int GetY() const = 0;

  /// Returns the of control points in z
  virtual int GetZ() const = 0;

  /// Returns the number of parameters of the transformation
  virtual int NumberOfDOFs() const = 0;

  /// Puts a control point value
  virtual void   Put(int, double) = 0;

  /// Gets a control point value
  virtual double Get(int) const = 0;

  /// Puts a control point status
  virtual void   PutStatus(int, _Status) = 0;

  /// Gets a control point status
  virtual _Status GetStatus(int) = 0;

  /// Subdivide FFD
  virtual void Subdivide() = 0;

  /// Transforms world coordinates (in mm) to lattice coordinates
  virtual void WorldToLattice(double &, double &, double &) const = 0;

  /// Transforms lattice coordinates to world coordinates (in mm)
  virtual void LatticeToWorld(double &, double &, double &) const = 0;

  /// Returns a string with the name of the instantiated class
  virtual const char *NameOfClass();

};

inline const char *irtkFreeFormTransformation::NameOfClass()
{
  return "irtkFreeFormTransformation";
}

#endif
