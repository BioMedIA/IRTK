extern class itkRigidTransformation : public itkHomogeneousTransformation {
public:

  /// Constructor (default)
  itkRigidTransformation();

  /// Constructor (copy)
  itkRigidTransformation(const itkRigidTransformation &);
  
  /// Destructor 
  virtual ~itkRigidTransformation();

  /// Puts translation along the x-axis (transformation matrix is updated)
  void   PutTranslationX(double);

  /// Gets translation along the x-axis
  virtual double GetTranslationX();

  /// Puts translation along the y-axis (transformation matrix is updated)
  virtual void   PutTranslationY(double);

  /// Gets translation along the y-axis
  virtual double GetTranslationY();

  /// Puts translation along the z-axis (transformation matrix is updated)
  virtual void   PutTranslationZ(double);

  /// Gets translation along the z-axis
  virtual double GetTranslationZ();

  /// Puts rotation angle around the x-axis (transformation matrix is updated)
  virtual void   PutRotationX(double);

  /// Gets rotation angle around the x-axis
  virtual double GetRotationX();

  /// Puts rotation angle around the y-axis (transformation matrix is updated)
  virtual void   PutRotationY(double);

  /// Gets rotation angle around the y-axis
  virtual double GetRotationY();

  /// Puts rotation angle around the z-axis (transformation matrix is updated)
  virtual void   PutRotationZ(double);

  /// Gets rotation angle around the z-axis
  virtual double GetRotationZ();

  /// Returns the number of parameters of the transformation
  virtual int NumberOfDOFs() const;

  /// Puts a transformation parameter (transformation matrix is updated)
  virtual void   Put(int, double);

  /// Gets a transformation parameter
  virtual double Get(int) const;

  /// Puts the transformation matrix (transformation parameters are updated)
  virtual void   PutMatrix(const itkMatrix &);

  /// Inverts the transformation
  virtual void Invert();

  /// Transforms a point by the rotation part of the rigid transformation.
  virtual void Rotate(double& INOUT, double& INOUT, double& INOUT);

  /// Checks whether transformation is an identity mapping
  virtual Bool IsIdentity();

  /// Prints the parameters of the transformation
  virtual void Print();

  /// Check file header
  static int CheckHeader(char *);

  /// Returns a string with the name of the instantiated class
  virtual char *NameOfClass();
};
