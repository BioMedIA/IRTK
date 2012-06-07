extern class irtkBaseImage {
public:

  /// Destructor
  virtual ~irtkBaseImage();

  //
  // Access functions for image dimensions
  //

  /// Returns the number of voxels in the x-direction
  virtual int  GetX() const;

  /// Returns the number of voxels in the y-direction
  virtual int  GetY() const;

  /// Returns the number of voxels in the z-direction
  virtual int  GetZ() const;

  /// Returns the number of voxels in the t-direction
  virtual int  GetT() const;

  /// Returns the total number of voxels
  virtual int  GetNumberOfVoxels() const;

  /// Gets the image attributes
  virtual irtkImageAttributes GetImageAttributes() const;
  
  /// Initialize image from attributes
  virtual void Initialize(const irtkImageAttributes &) = 0;

  /// Clear an image
  virtual void Clear() = 0;

  //
  // Access functions for voxel dimensions
  //

  /// Returns the number of voxels in the x-direction
  virtual double GetXSize() const;

  /// Returns the number of voxels in the y-direction
  virtual double GetYSize() const;

  /// Returns the number of voxels in the z-direction
  virtual double GetZSize() const;

  /// Returns the number of voxels in the t-direction
  virtual double GetTSize() const;

  /// Voxel dimensions get access
  // virtual void  GetPixelSize(double *OUTPUT, double *OUTPUT, double *OUTPUT) const;

  /// Voxel dimensions get access
  virtual void  GetPixelSize(double *OUTPUT, double *OUTPUT, double *OUTPUT, double *OUTPUT) const;

  /// Voxel dimensions put access
  virtual void  PutPixelSize(double, double, double);

  /// Voxel dimensions put access
  virtual void  PutPixelSize(double, double, double, double);

  /// Image origin get access
  virtual irtkPoint GetOrigin() const;

  /// Image origin get access
  virtual void  GetOrigin(double &, double &, double &) const;

  /// Image origin get access
  virtual void  GetOrigin(double &, double &, double &, double &) const;

  /// Image origin put access
  virtual void  PutOrigin(const irtkPoint &);

  /// Image origin put access
  virtual void  PutOrigin(double, double, double);

  /// Image origin put access
  virtual void  PutOrigin(double, double, double, double);

  /// Put image x- and y-axis and z-axis
%extend
{
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

  /// Get image x- and y-axis and z-axis
  // virtual void  GetOrientation(double *, double *, double *) const;

%extend
{
  %apply int& OUTPUT {int& xOUT};
  %apply int& OUTPUT {int& yOUT};
  %apply int& OUTPUT {int& zOUT};

  
  %feature("docstring", "orientation() -> (x, y, z)

  Returns a 3 3-tuples representing the axes of the image.");

  /// Get orientation of axis relative to patient
  virtual void  Orientation(int &xOUT, int &yOUT, int &zOUT) const
  {
    int x, y, z;
    self->Orientation(x, y, z);
    xOUT = x;
    yOUT = y;
    zOUT = z;
  }
}

  /// Image to world coordinate conversion with a given point
  virtual void ImageToWorld(irtkPoint &) const;

  /// World to image coordinate conversion with a given point
  virtual void WorldToImage(irtkPoint &) const;

  /// Image to world coordinate conversion with three doubles
  virtual void ImageToWorld(double &INOUT, double &INOUT, double &INOUT) const;

  /// World to image coordinate conversion with three doubles
  virtual void WorldToImage(double &INOUT, double &INOUT, double &INOUT) const;

  /// Return transformation matrix for image to world coordinates
  virtual irtkMatrix GetImageToWorldMatrix() const;

  /// Return transformation matrix for world to image coordinates
  virtual irtkMatrix GetWorldToImageMatrix() const;

  /// Image to time coordinate conversion
  virtual double ImageToTime(double) const;

  /// Time to image coordinate conversion
  virtual double TimeToImage(double) const;

  /// Returns true if point is within the field of view of image
  virtual bool IsInFOV(double, double, double);

  /// boolean operation for empty
  virtual bool IsEmpty() const;

%extend
{
  %apply double* OUTPUT {double* minOUT};
  %apply double* OUTPUT {double* maxOUT};
  
  %feature("docstring", "GetMinMaxAsDouble() -> (min, max)

  Returns a 2-tuple representing the min and max values of the image.");

  /// Minimum and maximum pixel values get accessor
  virtual void GetMinMaxAsDouble(double* minOUT, double* maxOUT) const
  {
    double minVal, maxVal;
    self->GetMinMaxAsDouble(&minVal, &maxVal);
    *minOUT = minVal;
    *maxOUT = maxVal;
  }
}

  /// Minimum and maximum pixel values put accessor
  virtual void PutMinMaxAsDouble(double, double);

  /// Function for pixel get access as double
  virtual double GetAsDouble(int, int, int, int = 0) const = 0;

  /// Function for pixel put access
  virtual void   PutAsDouble(int, int, int, double) = 0;

  /// Function for pixel put access
  virtual void   PutAsDouble(int, int, int, int, double) = 0;

  /// Returns the name of the image class
  // virtual const char *NameOfClass() = 0;

  /// Function for pixel access via pointers
  virtual void *GetScalarPointer(int = 0, int = 0, int = 0, int = 0) const = 0;

  /// Function which returns pixel scalar type
  virtual int GetScalarType() const = 0;

  /// Function which returns the minimum value the pixel can hold without overflowing
  virtual double GetScalarTypeMin() const = 0;

  /// Function which returns the minimum value the pixel can hold without overflowing
  virtual double GetScalarTypeMax() const = 0;

  /// Reflect image around x
  virtual void ReflectX() = 0;

  /// Reflect image around y
  virtual void ReflectY() = 0;

  /// Reflect image around z
  virtual void ReflectZ() = 0;

  /// Flip x and y axis
  virtual void FlipXY(int) = 0;

  /// Flip x and z axis
  virtual void FlipXZ(int) = 0;

  /// Flip y and z axis
  virtual void FlipYZ(int) = 0;

  /// Flip x and t axis
  virtual void FlipXT(int) = 0;

  /// Flip y and t axis
  virtual void FlipYT(int) = 0;

  /// Flip z and t axis
  virtual void FlipZT(int) = 0;

  /// Read file and construct image
  static irtkBaseImage *New(const char *);

  /// Read file and construct image
  static irtkBaseImage *New(const irtkBaseImage *);

  /// Write file
  virtual void Write(const char *) = 0;

  /// Print function
  virtual void Print();


%extend
{
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


