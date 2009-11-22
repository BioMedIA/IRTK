class itkLinearFreeFormTransformation : public itkFreeFormTransformation {

public:

  /// Constructor
  itkLinearFreeFormTransformation(double = 0, double = 0, double = 0, double = 1, double = 1, double = 1, double = 1, double = 1, double = 1, double * = NULL, double * = NULL);

  /// Copy Constructor
  itkLinearFreeFormTransformation(const class itkLinearFreeFormTransformation &);

  /// Copy Constructor
  itkLinearFreeFormTransformation(const class itkBSplineFreeFormTransformation &);

  /// Destructor
  virtual ~itkLinearFreeFormTransformation();

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
  virtual void Interpolate(double* dxs, double* dys, double* dzs);

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
  virtual void Jacobian(double, double, double, itkMatrix &);

  /// Calculate the determinant of the Jacobian of the transformation
  virtual double Jacobian(double, double, double);

  /// Calculate the bending energy of the transformation
  virtual double Bending(double x, double y, double z);

  /// Reads FFD from file
  virtual void Read (char *);

  /// Writes FFD to file
  virtual void Write(char *);

  /// Prints the parameters of the transformation
  virtual void Print();

  /// Check file header
  static int CheckHeader(char *);

  /// Returns a string with the name of the instantiated class
  virtual char *NameOfClass();
};
