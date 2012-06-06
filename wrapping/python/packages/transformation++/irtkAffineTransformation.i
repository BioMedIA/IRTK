extern class itkAffineTransformation : public itkRigidTransformation {
public:

  /// Constructor (default)
  itkAffineTransformation();

  /// Constructor (copy)
  itkAffineTransformation(const itkAffineTransformation &);

  /// Destructor 
  virtual ~itkAffineTransformation();

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

  /// Puts the transformation matrix (transformation parameters are updated)
  virtual void   PutMatrix(const itkMatrix &);

};
