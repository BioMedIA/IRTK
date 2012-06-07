%{
  enum PyStatus {PyActive, PyPassive};
%}

// Enumeration to represent active and passive DOFs.
enum PyStatus {PyActive, PyPassive};

// Avoid these names disappearing down a chain of sub-classes.
%rename (Transform_point)  irtkTransformation::Transform(irtkPoint& INOUT);
%rename (Transform_point_set) irtkTransformation::Transform(irtkPointSet &);

extern class irtkTransformation {
public:

  /** Static constructor. This functions returns a pointer to a concrete
   *  transformation by reading the transformation parameters from a file
   *  and creating the approriate transformation
   */
  static irtkTransformation *New(char *);

  /** Static constructor. This functions returns a pointer to a concrete
   *  transformation by copying the transformation passed to it
   */
  static irtkTransformation *New(irtkTransformation *);

  /// Virtual destructor
  virtual ~irtkTransformation();

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
  virtual void Transform(irtkPoint& INOUT);
  
  /// Transforms a set of points
  virtual void Transform(irtkPointSet &);

  /// Transforms a single point in 4D
  virtual void Transform(double& INOUT, double& INOUT, double& INOUT, double& INOUT = 0);

  /// Transforms a point using the global transformation component only
  virtual void GlobalTransform(double& INOUT, double& INOUT, double& INOUT, double = 0) = 0;

  /// Transforms a point using the local transformation component only
  virtual void LocalTransform (double& INOUT, double& INOUT, double& INOUT, double = 0) = 0;

  /// Calculate displacement
  virtual void Displacement(double& INOUT, double& INOUT, double& INOUT, double = 0) = 0;

  /// Calculates displacement using the global transformation component only
  virtual void GlobalDisplacement(double& INOUT, double& INOUT, double& INOUT, double = 0);

  /// Calculates displacement using the local transformation component only
  virtual void LocalDisplacement(double& INOUT, double& INOUT, double& INOUT, double = 0);

  /// Inverts the transformation (abstract)
  virtual double Inverse(double& INOUT, double& INOUT, double& INOUT, double = 0, double = 0.01) = 0;

  /// Calculate the Jacobian of the transformation with respect to the transformation parameters
  virtual void JacobianDOFs(double [3], int, double, double, double, double = 0);

  /// Calculate the Jacobian of the transformation with respect to world coordinates
  virtual void Jacobian(irtkMatrix &, double, double, double, double = 0) = 0;

  /// Calculate the Jacobian of the local transformation with respect to world coordinates
  virtual void LocalJacobian(irtkMatrix &, double, double, double, double = 0) = 0;

  /// Calculate the Jacobian of the global transformation with respect to world coordinates
  virtual void GlobalJacobian(irtkMatrix &, double, double, double, double = 0) = 0;

  /// Calculate the determinant of the Jacobian of the transformation with respect to world coordinates
  virtual double Jacobian(double, double, double, double = 0);

  /// Calculate the determinant of the Jacobian of the local transformation with respect to world coordinates
  virtual double LocalJacobian(double, double, double, double = 0);

  /// Calculate the determinant of the Jacobian of the global transformation with respect to world coordinates
  virtual double GlobalJacobian(double, double, double, double = 0);

  /// Calculate displacement vectors for image
  virtual void Displacement(irtkGenericImage<double> &);

  /// Checks whether transformation is an identity mapping (abstract)
  virtual bool IsIdentity() = 0;

  /// Reads a transformation from a file
  virtual void Read (char *);

  /// Writes a transformation to a file
  virtual void Write(char *);

  /// Imports a transformation from a file
  virtual void Import(char *);

  /// Exports a transformation to a file
  virtual void Export(char *);

  /// Reads a transformation from a file (abstract)
  virtual irtkCifstream& Read(irtkCifstream&) = 0;

  /// Writes a transformation to a file (abstract)
  virtual irtkCofstream& Write( irtkCofstream&) = 0;

  /// Imports a transformation from a file
  virtual istream& Import(istream&);

  /// Exports a transformation to a file
  virtual ostream& Export(ostream&);

  /// I/O
  virtual void Draw ();

  /// Prints the parameters of the transformation (abstract)
  virtual void Print() = 0;

  /// Returns a string with the name of the instantiated class (abstract)
  virtual const char *NameOfClass() = 0;

};



