extern class itkMultiLevelFreeFormTransformation : public itkAffineTransformation {

public:

  /// Constructor (default)
  itkMultiLevelFreeFormTransformation();

  /// Constructor (copy)
  itkMultiLevelFreeFormTransformation(const itkAffineTransformation &);

  /// Constructor (copy)
  itkMultiLevelFreeFormTransformation(const itkMultiLevelFreeFormTransformation &);

  /// Destructor
  virtual ~itkMultiLevelFreeFormTransformation();

  /// Returns the number of levels 
  virtual int NumberOfLevels() const;

  /// Gets local transformation
  virtual itkFreeFormTransformation *GetLocalTransformation(int);

  /// Puts local transformation
  virtual void PutLocalTransformation(itkFreeFormTransformation *, int);

  /// Push local transformation on stack
  virtual void PushLocalTransformation(itkFreeFormTransformation *);

  /// Pop local transformation from stack
  virtual itkFreeFormTransformation *PopLocalTransformation();

  /// Transforms a point
  virtual void Transform(double &INOUT, double &INOUT, double &INOUT);

  /// Calculates displacement using global and local transformation components
  virtual void Displacement(double& INOUT, double& INOUT, double& INOUT);

  /// Transforms a single point using the global transformation component only
  virtual void GlobalTransform(double &INOUT, double &INOUT, double &INOUT);

  /// Transforms a single point using the local transformation component only
  virtual void LocalTransform (double &INOUT, double &INOUT, double &INOUT);

  /// Calculates displacement using the global transformation component only
  virtual void GlobalDisplacement(double &INOUT, double &INOUT, double &INOUT);

  /// Calculates displacement using the local transformation component only
  virtual void LocalDisplacement(double &INOUT, double &INOUT, double &INOUT);

  /// Transforms a point
  virtual void Transform(int, double &INOUT, double &INOUT, double &INOUT);

  /// Transforms a single point using the local transformation component only
  virtual void LocalTransform (int, double &INOUT, double &INOUT, double &INOUT);

  /// Calculates displacement using the local transformation component only
  virtual void LocalDisplacement(int, double &INOUT, double &INOUT, double &INOUT);

  /// Calculate the Jacobian of the transformation
  virtual void Jacobian(double, double, double, itkMatrix &);

  /// Calculate the Jacobian of the local transformation
  virtual void LocalJacobian(double, double, double, itkMatrix &);

  /** Approximate displacements: This function takes a set of points and a
      set of displacements and find a FreeFD which approximates these
      displacements. After approximatation the displacements replaced by
      the residual displacement errors at the points */
  virtual double Approximate(double *, double *, double *, 
			     double *, double *, double *, int);

  /// Inverts the transformation
  virtual void Invert();

  /// Checks whether transformation is an identity mapping 
  virtual Bool IsIdentity();

  /// Prints the parameters of the transformation
  virtual void Print();

  /// Check file header
  static int CheckHeader(char *);

  /// Returns a string with the name of the instantiated class
  virtual char *NameOfClass();
};
