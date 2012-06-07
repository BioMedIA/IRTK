#define MAX_TRANS 200

extern class irtkMultiLevelFreeFormTransformation : public irtkAffineTransformation {

public:

  /// Local transformations
  irtkFreeFormTransformation *_localTransformation[MAX_TRANS+1];

  /// Number of local transformations
  int _NumberOfLevels;

  /// Constructor (default)
  irtkMultiLevelFreeFormTransformation();

  /// Constructor (copy)
  irtkMultiLevelFreeFormTransformation(const irtkRigidTransformation &);

  /// Constructor (copy)
  irtkMultiLevelFreeFormTransformation(const irtkAffineTransformation &);

  /// Constructor (copy)
  irtkMultiLevelFreeFormTransformation(const irtkMultiLevelFreeFormTransformation &);

  /// Destructor
  virtual ~irtkMultiLevelFreeFormTransformation();

  /// Returns the number of levels
  virtual int NumberOfLevels() const;

  /// Gets local transformation
  virtual irtkFreeFormTransformation *GetLocalTransformation(int);

  /// Puts local transformation
  virtual void PutLocalTransformation(irtkFreeFormTransformation *, int);

  /// Push local transformation on stack (append transformation)
  virtual void PushLocalTransformation(irtkFreeFormTransformation *);

  /// Insert local transformation
  virtual void InsertLocalTransformation(irtkFreeFormTransformation *, int = 0);

  /// Combine local transformation on stack
  virtual void CombineLocalTransformation();

  /// Pop local transformation from stack (remove last transformation)
  virtual irtkFreeFormTransformation *PopLocalTransformation();

  /// Remove local transformation and return the pointer (need to be deleted if not used)
  virtual irtkFreeFormTransformation *RemoveLocalTransformation(int = 0);

  /// Transforms a point
  virtual void Transform(double &INOUT, double &INOUT, double &INOUT, double = 0);

  /// Calculates displacement using global and local transformation components
  virtual void Displacement(double &INOUT, double &INOUT, double &INOUT, double = 0);

  /// Transforms a single point using the global transformation component only
  virtual void GlobalTransform(double &INOUT, double &INOUT, double &INOUT, double = 0);

  /// Transforms a single point using the local transformation component only
  virtual void LocalTransform (double &INOUT, double &INOUT, double &INOUT, double = 0);

  /// Calculates displacement using the global transformation component only
  virtual void GlobalDisplacement(double &INOUT, double &INOUT, double &INOUT, double = 0);

  /// Calculates displacement using the local transformation component only
  virtual void LocalDisplacement(double &INOUT, double &INOUT, double &INOUT, double = 0);

  /** Convert the global transformation from a matrix representation to a
      FFD and incorporate it with any existing local displacement. **/
  virtual void MergeGlobalIntoLocalDisplacement();

  // Helper function for the above.
  virtual void InterpolateGlobalDisplacement(irtkBSplineFreeFormTransformation3D *f);

  /// Transforms a point
  virtual void Transform(int, double &INOUT, double &INOUT, double &INOUT, double = 0);

  /// Transforms a single point using the local transformation component only
  virtual void LocalTransform (int, double &INOUT, double &INOUT, double &INOUT, double = 0);

  /// Calculates displacement using the local transformation component only
  virtual void LocalDisplacement(int, double &INOUT, double &INOUT, double &INOUT, double = 0);

  /// Calculate the Jacobian of the transformation
  virtual void Jacobian(irtkMatrix &, double, double, double, double = 0);

  /// Calculate the Jacobian of the local transformation
  virtual void LocalJacobian(irtkMatrix &, double, double, double, double = 0);

  /// Calculate the bending of the local transformation.
  virtual double Bending(double, double, double);

  /** Approximate displacements: This function takes a set of points and a
      set of displacements and find a FreeFD which approximates these
      displacements. After approximatation the displacements replaced by
      the residual displacement errors at the points */
  virtual double Approximate(double *, double *, double *,
                             double *, double *, double *, int);

  /// Inverts the transformation
  virtual double Inverse(double &INOUT, double &INOUT, double &INOUT, double = 0, double = 0.01);

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

%extend
{
  const char* __str__()
    {
      std::stringstream buffer;
      buffer << " MAX_TRANS " << MAX_TRANS << std::endl;
      return buffer.str().c_str();
    }
}

};
