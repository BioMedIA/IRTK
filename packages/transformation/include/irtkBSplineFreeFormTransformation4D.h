/*=========================================================================

  Library   : Image Registration Toolkit (IRTK)
  Module    : $Id$
  Copyright : Imperial College, Department of Computing
              Visual Information Processing (VIP), 2008 onwards
  Date      : $Date$
  Version   : $Revision$
  Changes   : $Author$

=========================================================================*/

#ifndef _IRTKBSPLINEFREEFORMTRANSFORMATION4D_H

#define _IRTKBSPLINEFREEFORMTRANSFORMATION4D_H


/**
 * Class for free form transformations based on tensor product B-splines.
 *
 * This class implements 4D free form transformation using B-splines.
 *
 * For more details about the implementation see Lee, Wolberg and Shin, IEEE
 * Transactions on Visualization and Computer Graphics, Vol. 3, No. 3, 1997.
 */

class irtkBSplineFreeFormTransformation4D : public irtkFreeFormTransformation4D
{

protected:

  /// B-spline basis functions
  irtkBSplineFunction _bspline;

  /// Subdivide FFD in 2D
  virtual void Subdivide2D();

  /// Subdivide FFD in 3D
  virtual void Subdivide3D();

public:

  /// Constructor
  irtkBSplineFreeFormTransformation4D();

  /// Construct FFD for given image domain
  irtkBSplineFreeFormTransformation4D(irtkBaseImage &, double, double, double, double);

  /** Construct FFD for given image domain: This constructor allows one to
   * alter the attributes of a given image before initializing the FFD, e.g.,
   * in order to add a time domain to a 3D image
   */
  irtkBSplineFreeFormTransformation4D(irtkImageAttributes &, double, double, double, double);

  /// Constructor
  irtkBSplineFreeFormTransformation4D(double x1, double y1, double z1, double t1,
                                      double x2, double y2, double z2, double t2,
                                      double dx, double dy, double dz, double dt,
                                      double* xaxis, double* yaxis, double* zaxis);

  /// Copy Constructor
  irtkBSplineFreeFormTransformation4D(const class irtkBSplineFreeFormTransformation4D &);

  /// Destructor
  virtual ~irtkBSplineFreeFormTransformation4D();

  /** Approximate displacements: This function takes a set of points and a
      set of displacements and find a FFD which approximates these
      displacements. After approximation the displacements replaced by
      the residual displacement errors at the points */
  virtual double Approximate(double *, double *, double *, double *,
                             double *, double *, double *, int);

  /** Interpolates displacements: This function takes a set of displacements
      defined at the control points and finds a FFD which interpolates these
      displacements.
      \param dxs The x-displacements at each control point.
      \param dys The y-displacements at each control point.
      \param dzs The z-displacements at each control point. */
  virtual void Interpolate(double* dxs, double* dys, double* dzs);

  /// Subdivide FFD
  virtual void Subdivide();

  /// Calculates the FFD (for a point in FFD coordinates) with checks
  void FFD1(double &, double &, double &, double) const;

  /// Calculates the FFD (for a point in FFD coordinates) without checks
  virtual void FFD2(double &, double &, double &, double) const;

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

  /// Calculate the Jacobian of the transformation with respect to the transformation parameters
  virtual void JacobianDOFs(double [3], int, double, double, double, double = 0);

  /// Calculate the bending energy of the transformation
  virtual double Bending(double, double, double, double);

  /** Returns the bounding box for a control point. The last parameter
   *  specifies what fraction of the bounding box to return. The default
   *  is 1 which equals 100% of the bounding box. Note that the bounding
   *  box is only computed in 3D (time is ignored).
   */
  virtual void BoundingBoxCP(int, irtkPoint &, irtkPoint &, double = 1) const;
  /// Bounding box in 4D
  virtual void BoundingBoxCP(int, double &, double &, double &, double &, double &,
                           double &, double &, double &, double = 1) const;

  /** Returns the bounding box for a control point (in pixels). The last
   *  parameter specifies what fraction of the bounding box to return. The
   *  default is 1 which equals 100% of the bounding box. Note that the
   *  bounding box is only computed in 3D (time is ignored).
   */
  virtual void BoundingBoxImage(irtkGreyImage *, int, int &, int &, int &,
  		                          int &, int &, int &, double = 1) const;

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
};

inline void irtkBSplineFreeFormTransformation4D::FFD1(double &x, double &y, double &z, double t) const
{
  if ((x <   -2) || (y <   -2) || (z <   -2) || (t <   -2) ||
      (x > _x+1) || (y > _y+1) || (z > _z+1) || (t > _t+1)) {
    x = 0;
    y = 0;
    z = 0;
  } else {
    this->FFD2(x, y, z, t);
  }
}

inline void irtkBSplineFreeFormTransformation4D::Transform(double &x, double &y, double &z, double t)
{
  double u, v, w;

  u = x;
  v = y;
  w = z;

  // Convert world coordinates in to FFD coordinates
  this->WorldToLattice(u, v, w);

  // Calculate FFD
  FFD1(u, v, w, this->TimeToLattice(t));

  // Add FFD to world coordinates
  x += u;
  y += v;
  z += w;
}

inline void irtkBSplineFreeFormTransformation4D::Transform2(double &x, double &y, double &z, double t)
{
  double u, v, w;

  u = x;
  v = y;
  w = z;

  // Convert world coordinates in to FFD coordinates
  this->WorldToLattice(u, v, w);

  // Calculate FFD
  this->FFD2(u, v, w, this->TimeToLattice(t));

  // Add FFD to world coordinates
  x += u;
  y += v;
  z += w;
}

inline void irtkBSplineFreeFormTransformation4D::GlobalTransform(double &, double &, double &, double)
{}

inline void irtkBSplineFreeFormTransformation4D::LocalTransform(double &x, double &y, double &z, double t)
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

inline void irtkBSplineFreeFormTransformation4D::GlobalDisplacement(double &x, double &y, double &z, double)
{
  x = 0;
  y = 0;
  z = 0;
}

inline void irtkBSplineFreeFormTransformation4D::LocalDisplacement(double &x, double &y, double &z, double t)
{
  // Convert world coordinates in to FFD coordinates
  this->WorldToLattice(x, y, z);

  // Calculate FFD
  FFD1(x, y, z, this->TimeToLattice(t));
}

inline void irtkBSplineFreeFormTransformation4D::Displacement(double &x, double &y, double &z, double t)
{
	this->LocalDisplacement(x, y, z, t);
}

inline const char *irtkBSplineFreeFormTransformation4D::NameOfClass()
{
  return "irtkBSplineFreeFormTransformation4D";
}

#endif
