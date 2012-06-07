extern class irtkHomogeneousTransformation : public irtkTransformation {
public:

  /// Constructor (default)
  irtkHomogeneousTransformation();

  /// Constructor (from matrix)
  irtkHomogeneousTransformation(const irtkMatrix &);

  /// Constructor (copy)
  irtkHomogeneousTransformation(const irtkHomogeneousTransformation &);

  /// Destructor
  virtual ~irtkHomogeneousTransformation();

  /// Returns the number of parameters of the transformation
  virtual int    NumberOfDOFs() const;

  /// Gets a transformation parameter
  virtual double Get(int) const;

  /// Puts a transformation paramater
  virtual void   Put(int, double);

  /// Gets the transformation matrix
  virtual irtkMatrix GetMatrix() const;

  /// Puts the transformation matrix (Argument is not checked)
  virtual void      PutMatrix(const irtkMatrix &);

  /// Transforms a single point
  virtual void Transform(double& INOUT, double& INOUT, double& INOUT, double& INOUT = 0);

  /// Transforms a point using the global transformation component only
  virtual void GlobalTransform(double& INOUT, double& INOUT, double& INOUT, double = 0);

  /// Transforms a point using the local transformation component only
  virtual void LocalTransform (double& INOUT, double& INOUT, double& INOUT, double = 0);

  /// Calculates displacement
  virtual void Displacement(double& INOUT, double& INOUT, double& INOUT, double = 0);

  /// Calculates displacement using the global transformation component only
  virtual void GlobalDisplacement(double& INOUT, double& INOUT, double& INOUT, double = 0);

  /// Calculates displacement using the local transformation component only
  virtual void LocalDisplacement(double& INOUT, double& INOUT, double& INOUT, double = 0);

  /// Inverts the transformation
  virtual void Invert();

  /// Inverse transformation
  virtual double Inverse(double& INOUT, double& INOUT, double& INOUT, double = 0, double = 0.01);

  /// Calculate the Jacobian of the transformation
  virtual void Jacobian(irtkMatrix &, double, double, double, double = 0);

  /// Calculate the Jacobian of the local transformation
  virtual void LocalJacobian(irtkMatrix &, double, double, double, double = 0);

  /// Calculate the Jacobian of the global transformation
  virtual void GlobalJacobian(irtkMatrix &, double, double, double, double = 0);

  /// Checks whether transformation is an identity mapping
  virtual bool IsIdentity();

  /// Prints the parameters of the transformation
  virtual void Print();

  /// Returns a string with the name of the instantiated class
  virtual const char *NameOfClass();

  /// Reads a transformation from a file
  virtual irtkCifstream& Read(irtkCifstream&);

  /// Writes a transformation to a file
  virtual irtkCofstream& Write(irtkCofstream&);

  /// Imports a transformation from a file
  virtual void Import(char *);

  /// Exports a transformation to a file
  virtual void Export(char *);

  /// Imports a transformation from a file
  virtual istream& Import(istream&);

  /// Exports a transformation to a file
  virtual ostream& Export(ostream&);

};
