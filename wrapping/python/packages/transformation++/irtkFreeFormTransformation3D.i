extern class irtkFreeFormTransformation3D : public irtkFreeFormTransformation {
public:

  /// Returns the of control points in x
  virtual int GetX() const;

  /// Returns the of control points in y
  virtual int GetY() const;

  /// Returns the of control points in z
  virtual int GetZ() const;

  /// Returns the of control point spacing in x
  virtual double GetXSpacing() const;

  /// Returns the of control point spacing in y
  virtual double GetYSpacing() const;

  /// Returns the of control point spacing in z
  virtual double GetZSpacing() const;

  /// Returns the number of parameters of the transformation
  virtual int NumberOfDOFs() const;

  /// Get the control point spacing (in mm)
  virtual void GetSpacing(double &INOUT, double &INOUT, double &INOUT) const;

%extend
{
  /// Put orientation of free-form deformation
  virtual void PutOrientation(double x0, double x1, double x2,
                      double y0, double y1, double y2,
		      double z0, double z1, double z2)
  {
    double xaxis[3], yaxis[3], zaxis[3];
    
    xaxis[0] = x0;
    xaxis[1] = x1;
    xaxis[2] = x2;
    yaxis[0] = y0;
    yaxis[1] = y1;
    yaxis[2] = y2;
    zaxis[0] = z0;
    zaxis[1] = z1;
    zaxis[2] = z2;
    
    self->PutOrientation(xaxis, yaxis, zaxis);
  }
}

  /// Get orientation of free-form deformation
  virtual void  GetOrientation(double *OUTPUT, double *OUTPUT, double *OUTPUT) const;

  /// Puts a control point value
  virtual void   Put(int, double);

  /// Gets a control point value
  virtual void   Put(int, int, int, double, double, double);

  /// Gets a control point value
  virtual double Get(int) const;

  /// Gets a control point value
  virtual void   Get(int, int, int, double &INOUT, double &INOUT, double &INOUT) const;


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

}

  /// Calculate the bending energy of the transformation
  virtual double Bending(double x, double y, double z, double t = 0) = 0;

  /// Transforms world coordinates (in mm) to FFD coordinates
  virtual void WorldToLattice(double &INOUT, double &INOUT, double &INOUT) const;

  /// Transforms world coordinates (in mm) to FFD coordinates
  virtual void WorldToLattice(irtkPoint &) const;

  /// Transforms FFD coordinates to world coordinates (in mm)
  virtual void LatticeToWorld(double &INOUT, double &INOUT, double &INOUT) const;

  /// Transforms FFD coordinates to world coordinates (in mm)
  virtual void LatticeToWorld(irtkPoint &) const;

  /// Transforms index of control points to FFD coordinates
  virtual void IndexToLattice(int index, int& OUTPUT, int& OUTPUT, int& OUTPUT) const;

  /// Transforms  FFD coordinates to index of control point
  virtual int  LatticeToIndex(int i, int j, int k) const;

  /// Returns the control point location (in mm)
  virtual void  ControlPointLocation(int, double &INOUT, double &INOUT, double &INOUT) const;

  /// Returns the control point location (in mm)
  virtual irtkPoint ControlPointLocation(int) const;

  /// Put the bounding box for FFD (in mm)
  virtual void PutBoundingBox(irtkPoint, irtkPoint);

  /// Put the bounding box for FFD (in mm)
  virtual void PutBoundingBox(double, double, double,
                              double, double, double);

  /// Returns the bounding box for FFD (in mm)
  virtual void BoundingBox(irtkPoint &, irtkPoint &) const;


  /** Returns the bounding box for a control point (in mm). The last
   *  parameter specifies what fraction of the bounding box to return. The
   *  default is 1 which equals 100% of the bounding box.
   */
  virtual void BoundingBoxCP(int, irtkPoint &, irtkPoint &, double = 1) const = 0;


  /** Returns the bounding box for a control point (in pixels). The last
   *  parameter specifies what fraction of the bounding box to return. The
   *  default is 1 which equals 100% of the bounding box.
   */
  virtual void BoundingBoxImage(irtkGreyImage *, int, int &INOUT, int &INOUT, int &INOUT,
			   int &INOUT, int &INOUT, int &INOUT, double = 1) const = 0;

  /** Approximate displacements: This function takes a set of points and a
      set of displacements and find a FFD which approximates these
      displacements. After approximatation the displacements replaced by
      the residual displacement errors at the points */
  virtual double Approximate(const double *, const double *, const double *, double *, double *, double *, int) = 0;

  /** Interpolates displacements: This function takes a set of displacements
      defined at the control points and finds a FFD which interpolates these
      displacements.
      \param dxs The x-displacements at each control point.
      \param dys The y-displacements at each control point.
      \param dzs The z-displacements at each control point. */
  virtual void Interpolate(const double* DOUBLEARRAYINPUT, const double* DOUBLEARRAYINPUT, const double* DOUBLEARRAYINPUT) = 0;
 
  /// Inverts the transformation (abstract)
  virtual double Inverse(double &INOUT, double &INOUT, double &INOUT, double = 0, double = 0.01);

  /// Checks whether transformation is an identity mapping
  virtual bool IsIdentity();

  /// Returns a string with the name of the instantiated class
  virtual const char *NameOfClass();

};
