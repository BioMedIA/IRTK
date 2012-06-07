extern class irtkAffineTransformation : public irtkRigidTransformation {
public:

  /// Constructor (default)
  irtkAffineTransformation();

  /// Constructor (copy)
  irtkAffineTransformation(const irtkRigidTransformation &);

  /// Constructor (copy)
  irtkAffineTransformation(const irtkAffineTransformation &);

  /// Destructor
  virtual ~irtkAffineTransformation();

  /// Reset transformation
  virtual void Reset();

  /// Puts scaling factor along the x-axis
  virtual void   PutScaleX(double);

  /// Gets scaling factor along the x-axis
  virtual double GetScaleX();

  /// Puts scaling factor along the y-axis
  virtual void   PutScaleY(double);

  /// Gets scaling factor along the y-axis
  virtual double GetScaleY();

  /// Puts scaling factor along the z-axis
  virtual void   PutScaleZ(double);

  /// Gets scaling factor along the z-axis
  virtual double GetScaleZ();

  /// Puts y-dependent skewing angle in the x direction (in degrees)
  virtual void   PutShearXY(double);

  /// Gets y-dependent skewing angle in the x direction (in degrees)
  virtual double GetShearXY();

  /// Puts z-dependent skewing angle in the y direction (in degrees)
  virtual void   PutShearYZ(double);

  /// Gets z-dependent skewing angle in the y direction (in degrees)
  virtual double GetShearYZ();

  /// Puts z-dependent skewing angle in the x direction (in degrees)
  virtual void   PutShearXZ(double);

  /// Gets z-dependent skewing angle in the x direction (in degrees)
  virtual double GetShearXZ();

  /// Returns the number of parameters of the transformation
  virtual int NumberOfDOFs() const;

  /// Puts a transformation parameter (transformation matrix is updated)
  virtual void   Put(int, double);

  /// Gets a transformation parameter
  virtual double Get(int) const;

  /// Calculate the Jacobian of the transformation with respect to the transformation parameters
  virtual void JacobianDOFs(double [3], int, double, double, double, double = 0);

  /// Checks whether transformation is an identity mapping
  virtual bool IsIdentity();

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

  /// Puts the transformation matrix (transformation parameters are updated)
  virtual void PutMatrix(const irtkMatrix &);

  /// Update transformation matrix
  virtual void UpdateMatrix();

  /// Updates transformation parameters based on current matrix.
  virtual void UpdateParameter();

};
