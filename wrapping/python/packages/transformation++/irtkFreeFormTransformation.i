extern class irtkFreeFormTransformation : public irtkTransformation {
public:

  /// Destructor
  virtual ~irtkFreeFormTransformation();

  /// Returns the of control points in x
  virtual int GetX() const = 0;

  /// Returns the of control points in y
  virtual int GetY() const = 0;

  /// Returns the of control points in z
  virtual int GetZ() const = 0;

  /// Returns the number of parameters of the transformation
  virtual int NumberOfDOFs() const = 0;

  /// Puts a control point value
  virtual void   Put(int, double) = 0;

  /// Gets a control point value
  virtual double Get(int) const = 0;

%extend
{
  %feature("docstring", "PutStatus(i, Status status) -> None

Sets the control point status for the given DOF.") irtkFreeFormTransformation::PutStatus;

  virtual void PutStatus(int i, PyStatus status)
    {
      if (status == PyActive)
        self->PutStatus(i, _Active);
      else
        self->PutStatus(i, _Passive);
    }

  %apply int& OUTPUT {Status& p0OUT};
  %apply int& OUTPUT {Status& p1OUT};
  %apply int& OUTPUT {Status& p2OUT};

  %feature("docstring", "GetStatus(i) -> status

Returns the control point status for the given DOF.") irtkFreeFormTransformation::GetStatus;
  virtual PyStatus GetStatus(int i)
    {
      if (self->GetStatus(i) == _Active)
        return PyActive;
      else
        return PyPassive;
    }
}

  /// Subdivide FFD
  virtual void Subdivide() = 0;

  /// Transforms world coordinates (in mm) to lattice coordinates
  virtual void WorldToLattice(double &INOUT, double &INOUT, double &INOUT) const = 0;

  /// Transforms lattice coordinates to world coordinates (in mm)
  virtual void LatticeToWorld(double &INOUT, double &INOUT, double &INOUT) const = 0;

  /// Returns a string with the name of the instantiated class
  virtual const char *NameOfClass();

  /// Calculate the bending of the transformation.
  virtual double Bending(double, double, double, double = 0.0) = 0;

};
