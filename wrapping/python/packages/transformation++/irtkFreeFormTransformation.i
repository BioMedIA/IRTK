extern class itkFreeFormTransformation : public itkTransformation {
public:

  /** Approximate displacements: This function takes a set of points and a
      set of displacements and find a FFD which approximates these
      displacements. After approximatation the displacements replaced by
      the residual displacement errors at the points */
  virtual double Approximate(double *, double *, double *, 
			     double *, double *, double *, int) = 0;

  /** Interpolates displacements: This function takes a set of displacements
      defined at the control points and finds a FFD which interpolates these
      displacements.
      \param dxs The x-displacements at each control point.
      \param dys The y-displacements at each control point.
      \param dzs The z-displacements at each control point. */
  virtual void Interpolate(double* dxs, double* dys, double* dzs) = 0;

  /// Subdivide FFD
  virtual void Subdivide() = 0;

  /// Returns the of control points in x
  virtual int GetX() const;

  /// Returns the of control points in y
  virtual int GetY() const;

  /// Returns the of control points in z
  virtual int GetZ() const;

  /// Returns the number of parameters of the transformation
  virtual int NumberOfDOFs() const;

  /// Get the control point spacing (in mm)
  virtual void GetSpacing(double &OUTPUT, double &OUTPUT, double &OUTPUT) const;

  /// Put orientation of free-form deformation
  virtual void  PutOrientation(double *DOUBLE3INPUT, double *DOUBLE3INPUT);

%extend
{
  %apply double* OUTPUT {double* p0OUT};
  %apply double* OUTPUT {double* p1OUT};
  %apply double* OUTPUT {double* p2OUT};
  %feature("docstring", "GetXAxis() -> (x0, x1, x2)

Returns a 3-tuple representing the x-axis of the free-form transformation.");
  void GetXAxis(double* p0OUT, double* p1OUT, double* p2OUT) const
  {
    double xaxis[3], yaxis[3];

    self->GetOrientation(xaxis, yaxis);

    *p0OUT = xaxis[0];
    *p1OUT = xaxis[1];
    *p2OUT = xaxis[2];
  } 

  %feature("docstring", "GetYAxis() -> (y0, y1, y2)

Returns a 3-tuple representing the y-axis of the free-form transformation.");
  void GetYAxis(double* p0OUT, double* p1OUT, double* p2OUT) const
  {
    double xaxis[3], yaxis[3];

    self->GetOrientation(xaxis, yaxis);

    *p0OUT = yaxis[0];
    *p1OUT = yaxis[1];
    *p2OUT = yaxis[2];
  } 

  %feature("docstring", "GetZAxis() -> (z0, z1, z2)

Returns a 3-tuple representing the z-axis of the free-form transformation.");
  void GetZAxis(double* p0OUT, double* p1OUT, double* p2OUT) const
  {
    double xaxis[3], yaxis[3], zaxis[3];

    self->GetOrientation(xaxis, yaxis, zaxis);

    *p0OUT = zaxis[0];
    *p1OUT = zaxis[1];
    *p2OUT = zaxis[2];
  }
}

  /// Puts a control point value
  virtual void   Put(int, double);

  /// Gets a control point value
  virtual void   Put(int, int, int, double, double, double);

  /// Gets a control point value
  virtual double Get(int) const;

  /// Gets a control point value
  virtual void   Get(int, int, int, double &OUTPUT, double &OUTPUT, double &OUTPUT) const;

%extend
{
  %feature("docstring", "PutStatus(i, j, k, status) -> None
PutStatus(i, j, k, statusx, statusy, statusz) -> None

Sets the control point status at the given lattice coordinates.

PutStatus(i, Status status) -> None

Sets the control point status for the given DOF.") itkFreeFormTransformation::PutStatus;
  virtual void PutStatus(int i, int j, int k, PyStatus status)
    {
      if (status == PyActive)
        self->PutStatus(i, j, k, _Active);
      else
        self->PutStatus(i, j, k, _Passive);
    }

  virtual void PutStatus(int i, int j, int k, PyStatus statusx, PyStatus statusy, PyStatus statusz)
    {
      _Status statusx2 = _Passive, statusy2 = _Passive, statusz2 = _Passive;

      if (statusx == PyActive)
        statusx2 = _Active;

      if (statusy == PyActive)
        statusy2 = _Active;

      if (statusz == PyActive)
        statusz2 = _Active;

      self->PutStatus(i, j, k, statusx2, statusy2, statusz2);
    }

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

  %feature("docstring", "GetStatus(i, j, k) -> (statusx, statusy, statusz)

Returns the control point status for the given lattice coordinates.

GetStatus(i) -> status

Returns the control point status for the given DOF.") itkFreeFormTransformation::GetStatus;
  virtual void   GetStatus(int i, int j, int k, PyStatus &p0OUT, PyStatus &p1OUT, PyStatus &p2OUT)
    {
      p0OUT = PyPassive;
      p1OUT = PyPassive;
      p2OUT = PyPassive;

      _Status statusx, statusy, statusz;
      self->GetStatus(i, j, k, statusx, statusy, statusz);

      if (statusx == _Active)
        p0OUT = PyActive;

      if (statusy == _Active)
        p1OUT = PyActive;

      if (statusz == _Active)
        p2OUT = PyActive;
    }

  virtual PyStatus GetStatus(int i)
    {
      if (self->GetStatus(i) == _Active)
        return PyActive;
      else
        return PyPassive;
    }
}

  /// Transforms world coordinates (in mm) to FFD coordinates
  virtual void WorldToLattice(double &INOUT, double &INOUT, double &INOUT) const;

  /// Transforms world coordinates (in mm) to FFD coordinates
  virtual void WorldToLattice(itkPoint &) const;

  /// Transforms FFD coordinates to world coordinates (in mm)
  virtual void LatticeToWorld(double &INOUT, double &INOUT, double &INOUT) const;

  /// Transforms FFD coordinates to world coordinates (in mm)
  virtual void LatticeToWorld(itkPoint &) const;

  /// Transforms index of control points to FFD coordinates
  virtual void IndexToLattice(int index, int& INOUT, int& INOUT, int& INOUT) const;
  
  /// Transforms  FFD coordinates to index of control point
  virtual int  LatticeToIndex(int i, int j, int k) const; 

  /// Returns the control point location (in mm)
  virtual void  ControlPointLocation(int, double &OUTPUT, double &OUTPUT, double &OUTPUT) const;

  /// Put the bounding box for FFD (in mm)
  virtual void PutBoundingBox(itkPoint, itkPoint);

  /// Put the bounding box for FFD (in mm)
  virtual void PutBoundingBox(double, double, double, 
			      double, double, double);

  /// Returns the bounding box for FFD (in mm)
  virtual void BoundingBox(itkPoint &, itkPoint &) const;
  
  /// Returns the bounding box for FFD (in mm)
  virtual void BoundingBox(double &OUTPUT, double &OUTPUT, double &OUTPUT, 
 			   double &OUTPUT, double &OUTPUT, double &OUTPUT) const;

  /// Checks whether transformation is an identity mapping
  virtual Bool IsIdentity();

  /// Prints the parameters of the transformation
  virtual void Print();

  /// Returns a string with the name of the instantiated class
  virtual char *NameOfClass();

};
