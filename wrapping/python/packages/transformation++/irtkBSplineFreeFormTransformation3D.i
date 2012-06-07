//%feature("notabstract") irtkBSplineFreeFormTransformation3D;
extern class irtkBSplineFreeFormTransformation3D : public irtkFreeFormTransformation3D {
public:

  /// Returns the value of the B-spline basis function
  static double B (double);

  /// Returns the value of the i-th B-spline basis function
  static double B (int, double);

  /// Returns the 1st derivative value of the i-th B-spline basis function
  static double B_I(int, double);

  /// Returns the 2nd derivative value of the i-th B-spline basis function
  static double B_II(int, double);

  /// Constructor
  irtkBSplineFreeFormTransformation3D();

  /// Constructor
  irtkBSplineFreeFormTransformation3D(irtkBaseImage &, double, double, double);

  /// Constructor
  irtkBSplineFreeFormTransformation3D(double x1, double y1, double z1,
                                      double x2, double y2, double z2,
                                      double dx, double dy, double dz,
                                      double* xaxis, double* yaxis, double* zaxis);

  /// Copy Constructor
  irtkBSplineFreeFormTransformation3D(const class irtkBSplineFreeFormTransformation3D &);

  /// Destructor
  virtual ~irtkBSplineFreeFormTransformation3D();

  /** Approximate displacements: This function takes a set of points and a
      set of displacements and find a FFD which approximates these
      displacements. After approximation the displacements replaced by
      the residual displacement errors at the points */
  virtual double Approximate(const double *, const double *, const double *, double *, double *, double *, int);

  /** Approximate displacements: This function takes a set of points and a
      set of displacements and find a !new! FFD which approximates these
      displacements. After approximation the displacements replaced by
      the residual displacement errors at the points */
  virtual void ApproximateAsNew(const double *, const double *, const double *, double *, double *, double *, int);

  /** Interpolates displacements: This function takes a set of displacements
      defined at the control points and finds a FFD which interpolates these
      displacements.
      \param dxs The x-displacements at each control point.
      \param dys The y-displacements at each control point.
      \param dzs The z-displacements at each control point. */
  virtual void Interpolate(const double* DOUBLEARRAYINPUT, const double* DOUBLEARRAYINPUT, const double* DOUBLEARRAYINPUT);

  /// Subdivide FFD
  virtual void Subdivide();

  /// Transforms a point
  virtual void Transform(double &INOUT, double &INOUT, double &INOUT, double = 0);

  /// Transforms a point using the global transformation component only
  virtual void GlobalTransform(double &INOUT, double &INOUT, double &INOUT, double = 0);

  /// Transforms a point using the local transformation component only
  virtual void LocalTransform (double &INOUT, double &INOUT, double &INOUT, double = 0);

  /// Calculates displacement
  virtual void Displacement(double &INOUT, double &INOUT, double &INOUT, double = 0);

  /// Calculates displacement using the global transformation component only
  virtual void GlobalDisplacement(double &INOUT, double &INOUT, double &INOUT, double = 0);

  /// Calculates displacement using the local transformation component only
  virtual void LocalDisplacement(double &INOUT, double &INOUT, double &INOUT, double = 0);

  /// Calculate the Jacobian of the transformation
  virtual void Jacobian(irtkMatrix &, double, double, double, double = 0);

  /// Calculate the Jacobian of the local transformation
  virtual void LocalJacobian(irtkMatrix &, double, double, double, double = 0);

  /// Calculate the Jacobian of the local transformation
  virtual void JacobianDetDerivative(irtkMatrix *, int, int, int);

  /// Calculate the Jacobian of the global transformation
  virtual void GlobalJacobian(irtkMatrix &, double, double, double, double = 0);

  /// Calculate the Jacobian of the transformation with respect to the transformation parameters
  virtual void JacobianDOFs(double [3], int, double, double, double, double = 0);

  /// Calculate total bending energy
  virtual double Bending();

  /// Calculate bending energy
  virtual double Bending(double, double, double, double = 0);

  /// Calculate the gradient of the bending energy with respect to the parameters
  virtual void BendingGradient(double *gradient);

  /** Returns the bounding box for a control point (in mm). The last
   *  parameter specifies what fraction of the bounding box to return. The
   *  default is 1 which equals 100% of the bounding box.
   */
  virtual void BoundingBox(int, irtkPoint &, irtkPoint &, double = 1) const;

  /** Returns the bounding box for a control point (in mm). The last
   *  parameter specifies what fraction of the bounding box to return. The
   *  default is 1 which equals 100% of the bounding box.
   */
  virtual void BoundingBox(int, double &INOUT, double &INOUT, double &INOUT,
			   double &INOUT, double &INOUT, double &INOUT, double = 1) const;

  /** Returns the bounding box for a control point (in pixels). The last
   *  parameter specifies what fraction of the bounding box to return. The
   *  default is 1 which equals 100% of the bounding box.
   */
  virtual void BoundingBox(irtkGreyImage *, int, int &INOUT, int &INOUT, int &INOUT,
			   int &INOUT, int &INOUT, int &INOUT, double = 1) const;

  /** Returns the bounding box for a control point (in pixels). The last
   *  parameter specifies what fraction of the bounding box to return. The
   *  default is 1 which equals 100% of the bounding box.
   */
  virtual void MultiBoundingBox(irtkGreyImage *, int, int &, int &, int &,
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

  /// Imports a transformation from a file
  virtual istream& Import(istream&);

  /// Exports a transformation to a file
  virtual ostream& Export(ostream&);
};
