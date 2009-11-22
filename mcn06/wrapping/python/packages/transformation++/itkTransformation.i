%{
  enum PyStatus {PyActive, PyPassive};
%}

// Enumeration to represent active and passive DOFs.
enum PyStatus {PyActive, PyPassive};

extern class itkTransformation {
public:
  /** Static constructor. This functions returns a pointer to a concrete
   *  transformation by reading the transformation parameters from a file
   *  and creating the approriate transformation
   */
  static itkTransformation *New(char *);

  /// Virtual destructor
  virtual ~itkTransformation();

  /// Returns the number of parameters of the transformation (abstract)
  virtual int    NumberOfDOFs() const = 0;

  /// Gets a transformation parameter (abstract)
  virtual double Get(int) const = 0;

  /// Puts a transformation paramater (abstract)
  virtual void   Put(int, double) = 0;

%extend
{
  /// Puts a control point status
  virtual void   PutStatus(int i, PyStatus status)
  {
    if (status == PyActive)
      self->PutStatus(i, _Active);
    else
      self->PutStatus(i, _Passive);
  }

  /// Gets a control point status
  virtual PyStatus GetStatus(int i) const
  {
    _Status status = self->GetStatus( i);
    if (status == _Active)
      return PyActive;
    else
      return PyPassive;
  }
}

  /// Transforms a single point
  virtual void Transform(itkPoint &);

  /// Transforms a set of points
  virtual void Transform(itkPointSet &);
  
  /// Transforms a point.
  virtual void Transform(double& INOUT, double& INOUT, double& INOUT);

  /// Calculate the Jacobian of the transformation
  virtual void Jacobian(double, double, double, itkMatrix &) = 0;

  /// Calculate the determinant of the Jacobian of the transformation
  virtual double Jacobian(double, double, double);

  /// Checks whether transformation is an identity mapping (abstract)
  virtual Bool IsIdentity() = 0;

  /// I/O
  virtual void Draw ();

  /// Prints the parameters of the transformation (abstract)
  virtual void Print() = 0;

  /// Returns a string with the name of the instantiated class (abstract)
  virtual char *NameOfClass() = 0;
};
