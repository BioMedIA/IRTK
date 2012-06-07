extern class irtkBSplineFreeFormTransformation : public irtkFreeFormTransformation {
public:

  /// Returns the value of the i-th B-spline basis function
  static double B (int, double);

  /// Returns the 1st derivative value of the i-th B-spline basis function
  static double B_I(int, double);

  /// Returns the 2nd derivative value of the i-th B-spline basis function
  static double B_II(int, double);

  /// Constructor
  irtkBSplineFreeFormTransformation(double, double, double, double, double, double, double, double, double, double*, double*);

  /// Copy Constructor
  irtkBSplineFreeFormTransformation(const class irtkBSplineFreeFormTransformation &);

  /// Destructor
  virtual ~irtkBSplineFreeFormTransformation();

  /** Approximate displacements: This function takes a set of points and a
      set of displacements and find a FFD which approximates these
      displacements. After approximatation the displacements replaced by
      the residual displacement errors at the points */
  virtual double Approximate(double *, double *, double *, 
			     double *, double *, double *, int);

  /** Interpolates displacements: This function takes a set of displacements
      defined at the control points and finds a FFD which interpolates these
      displacements.
      \param dxs The x-displacements at each control point.
      \param dys The y-displacements at each control point.
      \param dzs The z-displacements at each control point. */
  virtual void Interpolate(double* DOUBLEARRAYINPUT, double* DOUBLEARRAYINPUT, double* DOUBLEARRAYINPUT);

  /// Subdivide FFD
  virtual void Subdivide();

  /// Calculates the FFD (for a point in FFD coordinates) with checks
  virtual void FFD1(double &INOUT, double &INOUT, double &INOUT) const;

  /// Calculates the FFD (for a point in FFD coordinates) without checks
  virtual void FFD2(double &INOUT, double &INOUT, double &INOUT) const;

  /// Transforms a point
  virtual void Transform(double &INOUT, double &INOUT, double &INOUT);

  /// Transforms a point
  virtual void Transform2(double &INOUT, double &INOUT, double &INOUT);

  /// Transforms a point using the global transformation component only
  virtual void GlobalTransform(double &INOUT, double &INOUT, double &INOUT);

  /// Transforms a point using the local transformation component only
  virtual void LocalTransform (double &INOUT, double &INOUT, double &INOUT);

  /// Calculates displacement using the global transformation component only
  virtual void GlobalDisplacement(double &INOUT, double &INOUT, double &INOUT);

  /// Calculates displacement using the local transformation component only
  virtual void LocalDisplacement(double &INOUT, double &INOUT, double &INOUT);

  /// Calculate the Jacobian of the transformation
  virtual void Jacobian(double, double, double, irtkMatrix &);

  /// Calculate the bending energy of the transformation
  virtual double Bending(double x, double y, double z);

  /// Prints the parameters of the transformation
  virtual void Print();

  /// Check file header
  static int CheckHeader(char *);

  /// Returns a string with the name of the instantiated class
  virtual char *NameOfClass();
};
