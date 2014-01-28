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

#ifndef _IRTKLINEARFREEFORMTRANSFORMATION_H

#define _IRTKLINEARFREEFORMTRANSFORMATION_H


/**
 * Class for free form transformations
 *
 * This class implements 3D free form transformation
 */

class irtkLinearFreeFormTransformation : public irtkFreeFormTransformation3D
{
protected:

    /// Calculate the bending energy of the transformation at control points (3D)
    virtual double Bending3D(int i, int j, int k);

    /// Calculate the gradient of the bending energy with respect to the parameters
    virtual void BendingGradient3D(double *gradient);

public:

  /// Constructor
  irtkLinearFreeFormTransformation();

  /// Constructor
  irtkLinearFreeFormTransformation(irtkBaseImage &, double, double, double);

  /// Constructor
  irtkLinearFreeFormTransformation(irtkGenericImage<double> &);

  /// Constructor
  irtkLinearFreeFormTransformation(double x1, double y1, double z1,
                                   double x2, double y2, double z2,
                                   double dx, double dy, double dz,
                                   double* xaxis, double* yaxis, double* zaxis);

  /// Copy Constructor
  irtkLinearFreeFormTransformation(const class irtkLinearFreeFormTransformation &);

  /// Copy Constructor
  irtkLinearFreeFormTransformation(const class irtkBSplineFreeFormTransformation3D &);

  /// Destructor
  virtual ~irtkLinearFreeFormTransformation();

  /** Approximate displacements: This function takes a set of points and a
      set of displacements and find a FFD which approximates these
      displacements. After approximatation the displacements replaced by
      the residual displacement errors at the points */
  virtual double Approximate(const double *, const double *, const double *, double *, double *, double *, int);

  /** Interpolates displacements: This function takes a set of displacements
      defined at the control points and finds a FFD which interpolates these
      displacements.
      \param dxs The x-displacements at each control point.
      \param dys The y-displacements at each control point.
      \param dzs The z-displacements at each control point. */

      /** Approximate displacements: This function takes a set of points and a
      set of displacements and find a !new! FFD which approximates these
      displacements. After approximation the displacements replaced by
      the residual displacement errors at the points */
  virtual void ApproximateGradient(const double *, const double *, const double *, double *, double *, double *, int, double *);

  virtual void Interpolate(const double* dxs, const double* dys, const double* dzs);

  /// Subdivide FFD
  virtual void Subdivide();

  /// Calculates the FFD (for a point in FFD coordinates) with checks
  virtual void FFD1(double &, double &, double &) const;

  /// Calculates the FFD (for a point in FFD coordinates) without checks
  virtual void FFD2(double &, double &, double &) const;

  /// Transforms a point
  virtual void Transform(double &, double &, double &, double = 0);

  /// Transforms a point
  virtual void Transform2(double &, double &, double &, double = 0);

  /// Transforms a point using the global transformation component only
  virtual void GlobalTransform(double &, double &, double &, double = 0);

  /// Transforms a point using the local transformation component only
  virtual void LocalTransform (double &, double &, double &, double = 0);

  /// Calculates displacement
  virtual void Displacement(double &, double &, double &, double = 0);

  /// Calculates displacement using the global transformation component only
  virtual void GlobalDisplacement(double &, double &, double &, double = 0);

  /// Calculates displacement using the local transformation component only
  virtual void LocalDisplacement(double &, double &, double &, double = 0);

  /// Calculate the Jacobian of the transformation
  virtual void Jacobian(irtkMatrix &, double, double, double, double = 0);

  /// Calculate the Jacobian of the local transformation
  virtual void LocalJacobian(irtkMatrix &, double, double, double, double = 0);

  /// Calculate the Jacobian of the global transformation
  virtual void GlobalJacobian(irtkMatrix &, double, double, double, double = 0);

  /// Calculate the bending energy of the transformation
  virtual double Bending(double x, double y, double z, double t = 0);

  /// Calculate total bending energy
  virtual double Bending();

  /// Calculate the gradient of the bending energy with respect to the parameters
  virtual void BendingGradient(double *gradient);

  /// Inverts the transformation
  virtual double Inverse(double &, double &, double &, double = 0, double = 0.0001);

  /** Returns the bounding box for a control point (in mm). The last
   *  parameter specifies what fraction of the bounding box to return. The
   *  default is 1 which equals 100% of the bounding box.
   */
  virtual void BoundingBoxCP(int, irtkPoint &, irtkPoint &, double = 1) const;

  /** Returns the bounding box for a control point (in pixels). The last
   *  parameter specifies what fraction of the bounding box to return. The
   *  default is 1 which equals 100% of the bounding box.
   */
  virtual void BoundingBoxImage(irtkGreyImage *, int, int &, int &, int &,
                                int &, int &, int &, double = 1) const;
  virtual void BoundingBoxImage(irtkRealImage *, int, int &, int &, int &,
                                  int &, int &, int &, double = 1) const;

  /** Compose this transformation (T1) with second transformation (T2). The
   *  result is defined as T = T1 o T2
   */
  virtual void Compose(irtkTransformation *);

  /// Prints the parameters of the transformation
  virtual void Print();

  /// Check file header
  static int CheckHeader(char *);

  /// Returns a string with the name of the instantiated class
  virtual const char *NameOfClass();

  /// Reads a transformation from a file
  virtual irtkCifstream& Read(irtkCifstream&);

  /// Writes a transformation to a file
  virtual irtkCofstream& Write(irtkCofstream&);

  /// Imports a transformation from a file
  virtual istream& Import(istream&);

  /// Exports a transformation to a file
  virtual ostream& Export(ostream&);

};

inline void irtkLinearFreeFormTransformation::Transform(double &x, double &y, double &z, double)
{
  double u, v, w;

  u = x;
  v = y;
  w = z;

  // Convert world coordinates in to FFD coordinates
  this->WorldToLattice(u, v, w);

  // Calculate FFD
  this->FFD1(u, v, w);

  // Add FFD to world coordinates
  x += u;
  y += v;
  z += w;

}

inline void irtkLinearFreeFormTransformation::Transform2(double &x, double &y, double &z, double)
{
  double u, v, w;

  u = x;
  v = y;
  w = z;

  // Convert world coordinates in to FFD coordinates
  this->WorldToLattice(u, v, w);

  // Calculate FFD
  this->FFD2(u, v, w);

  // Add FFD to world coordinates
  x += u;
  y += v;
  z += w;
}

inline void irtkLinearFreeFormTransformation::GlobalTransform(double &, double &, double &, double)
{}

inline void irtkLinearFreeFormTransformation::LocalTransform(double &x, double &y, double &z, double t)
{
  double u, v, w;

  u = x;
  v = y;
  w = z;

  // calculate displacement
  this->LocalDisplacement(u, v, w, t);

  // Add displacement
  x += u;
  y += v;
  z += w;
}

inline void irtkLinearFreeFormTransformation::GlobalDisplacement(double &x, double &y, double &z, double)
{
  x = 0;
  y = 0;
  z = 0;
}

inline void irtkLinearFreeFormTransformation::LocalDisplacement(double &x, double &y, double &z, double)
{
  // Convert world coordinates in to FFD coordinates
  this->WorldToLattice(x, y, z);

  // Calculate FFD
  this->FFD1(x, y, z);
}

inline void irtkLinearFreeFormTransformation::Displacement(double &x, double &y, double &z, double)
{
	this->LocalDisplacement(x, y, z);
}

inline const char *irtkLinearFreeFormTransformation::NameOfClass()
{
  return "irtkLinearFreeFormTransformation";
}

#endif
