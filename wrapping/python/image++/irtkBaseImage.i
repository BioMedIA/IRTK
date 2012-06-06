extern class itkBaseImage {
public:

  /// Default constructor
  itkBaseImage();

  /// Constructor
  itkBaseImage(int, int, int);

  /// Constructor with given size and dimensions
  itkBaseImage(int, int, int, double, double, double);

  /// Copy constructor
  itkBaseImage(const itkBaseImage &);

  /// Destructor
  ~itkBaseImage();

  /// Update transformation matrix
  void UpdateMatrix();

  /// Initialize baseimage
  void Initialize(int, int, int);

  /// Initialize baseimage with given size and dimensions
  void Initialize(int, int, int, double, double, double);

  /// Initialize baseimage with given size, dimensions and origin
  void Initialize(const itkBaseImage &);

  /// Comparison Operator == (explicit negation replaces != operator)
  int       operator==(const itkBaseImage &);

  //
  // Access functions for image dimensions
  //

  /// Returns the number of voxels in the x-direction
  int  GetX() const;

  /// Returns the number of voxels in the y-direction
  int  GetY() const;

  /// Returns the number of voxels in the z-direction
  int  GetZ() const;

  /// Returns the total number of voxels
  int  GetNumberOfVoxels() const;

  //
  // Access functions for voxel dimensions
  //

  /// Returns the number of voxels in the x-direction
  double GetXSize() const;

  /// Returns the number of voxels in the y-direction
  double GetYSize() const;

  /// Returns the number of voxels in the z-direction
  double GetZSize() const;

  /// Voxel dimensions get access
  void  GetPixelSize(double *OUTPUT, double *OUTPUT, double *OUTPUT) const;

  /// Voxel dimensions put access
  void  PutPixelSize(double,   double,   double);

  /// Image origin get access
  itkPoint GetOrigin() const;

  /// Image origin get access
  void  GetOrigin(double &, double &, double &) const;

  /// Image origin put access
  void  PutOrigin(const itkPoint &);  

  /// Image origin put access
  void  PutOrigin(double, double, double);  

  /// Put image x- and y-axis
%extend
{
  void PutOrientation(double x0, double x1, double x2,
                      double y0, double y1, double y2)
  {
    double xaxis[3], yaxis[3];
    
    xaxis[0] = x0;
    xaxis[1] = x1;
    xaxis[2] = x2;
    yaxis[0] = y0;
    yaxis[1] = y1;
    yaxis[2] = y2;
    
    self->PutOrientation(xaxis, yaxis);
  }
}

  /// Image to world coordinate conversion with a given point
  void ImageToWorld(itkPoint &) const;

  /// World to image coordinate conversion with a given point
  void WorldToImage(itkPoint &) const;

  /// Image to world coordinate conversion with three doubles
  void ImageToWorld(double &INOUT, double &INOUT, double &INOUT) const;

  /// World to image coordinate conversion with three doubles
  void WorldToImage(double &INOUT, double &INOUT, double &INOUT) const;

  /// Return transformation matrix for image to world coordinates
  itkMatrix GetImageToWorldMatrix() const;

  /// Return transformation matrix for world to image coordinates
  itkMatrix GetWorldToImageMatrix() const;

  /// Returns true if point is within the field of view of image
  int IsInFOV(double, double, double);

  /// Print function
  void Print();

  /// Returns the name of the image class 
  virtual char *NameOfClass() = 0;

%extend
{
  %apply double* OUTPUT {double* p0OUT};
  %apply double* OUTPUT {double* p1OUT};
  %apply double* OUTPUT {double* p2OUT};
  %feature("docstring", "GetXAxis() -> (x0, x1, x2)

Returns a 3-tuple representing the x-axis of the image.");
  void GetXAxis(double* p0OUT, double* p1OUT, double* p2OUT) const
  {
    double xaxis[3], yaxis[3];

    self->GetOrientation(xaxis, yaxis);

    *p0OUT = xaxis[0];
    *p1OUT = xaxis[1];
    *p2OUT = xaxis[2];
  } 

  %feature("docstring", "GetYAxis() -> (y0, y1, y2)

Returns a 3-tuple representing the y-axis of the image.");
  void GetYAxis(double* p0OUT, double* p1OUT, double* p2OUT) const
  {
    double xaxis[3], yaxis[3];

    self->GetOrientation(xaxis, yaxis);

    *p0OUT = yaxis[0];
    *p1OUT = yaxis[1];
    *p2OUT = yaxis[2];
  } 

  %feature("docstring", "GetZAxis() -> (z0, z1, z2)

Returns a 3-tuple representing the z-axis of the image.");
  void GetZAxis(double* p0OUT, double* p1OUT, double* p2OUT) const
  {
    double xaxis[3], yaxis[3], zaxis[3];

    self->GetOrientation(xaxis, yaxis, zaxis);

    *p0OUT = zaxis[0];
    *p1OUT = zaxis[1];
    *p2OUT = zaxis[2];
  } 

  const char* __str__()
    {
      std::stringstream buffer;
      double origin[3];
      double xaxis[3], yaxis[3], zaxis[3];

      self->GetOrigin(origin[0], origin[1], origin[2]);
      self->GetOrientation(xaxis, yaxis, zaxis);

      buffer << "Image size is " << self->GetX() << " " << self->GetY() << " " << self->GetZ() << std::endl
             << "Voxel size is " << self->GetXSize() << " " << self->GetYSize() << " " << self->GetZSize() << std::endl
             << "Image origin is " << origin[0] << " " << origin[1] << " " << origin[2] << std::endl
             << "X-axis is " << xaxis[0] << " " << xaxis[1] << " " << xaxis[2] << std::endl
             << "Y-axis is " << yaxis[0] << " " << yaxis[1] << " " << yaxis[2] << std::endl
             << "Z-axis is " << zaxis[0] << " " << zaxis[1] << " " << zaxis[2];

      return buffer.str().c_str();
    }
}
};
