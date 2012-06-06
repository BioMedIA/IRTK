extern class itkHomogeneousTransformation : public itkTransformation {
public:

  /// Constructor (default)
  itkHomogeneousTransformation();

  /// Constructor (from matrix)
  itkHomogeneousTransformation(const itkMatrix &);
  
  /// Constructor (copy)
  itkHomogeneousTransformation(const itkHomogeneousTransformation &);
  
  /// Destructor 
  virtual ~itkHomogeneousTransformation();

  /// Returns the number of parameters of the transformation
  virtual int    NumberOfDOFs() const;

  /// Gets a transformation parameter
  virtual double Get(int) const;

  /// Puts a transformation paramater
  virtual void   Put(int, double);

  /// Gets the transformation matrix
  virtual itkMatrix GetMatrix() const;

  /// Puts the transformation matrix (Argument is not checked)
  virtual void      PutMatrix(const itkMatrix &);

  /// Transforms a single point
  virtual void Transform(double &INOUT, double &INOUT, double &INOUT);

  /// Inverts the transformation
  virtual void Invert();

  /// Calculate the Jacobian of the transformation
  virtual void Jacobian(double, double, double, itkMatrix &);

  /// Checks whether transformation is an identity mapping
  virtual Bool IsIdentity();

  /// Prints the parameters of the transformation
  virtual void Print();

  /// Returns a string with the name of the instantiated class
  virtual char *NameOfClass();
};
