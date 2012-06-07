template <class VoxelType> extern class irtkGenericImage : public irtkBaseImage
{
public:


  /// Default constructor
  irtkGenericImage();

  /// Constructor from image file
  irtkGenericImage(char *);

  /// Constructor for given image size
  irtkGenericImage(int, int, int, int = 1);

  /// Copy constructor for image 
  irtkGenericImage(const irtkGenericImage &);

  /// Constructor for given image attributes
  irtkGenericImage(const irtkImageAttributes &);

  /// Copy constructor for image of different type
  template <class T> irtkGenericImage(const irtkGenericImage<T> &);

  /// Destructor
  ~irtkGenericImage(void);

  /// Initialize an image
  void Initialize(const irtkImageAttributes &);

  /// Clear an image
  void Clear();

  /// Read image from file
  void Read (const char *);

  /// Write image to file
  void Write(const char *);

  /// Minimum and maximum pixel values get accessor
  void GetMinMax(VoxelType *OUTPUT, VoxelType *OUTPUT) const;

  /// Average pixel values get accessor
  VoxelType GetAverage(int = 1) const;

  /// Standard Deviation of the pixels
  VoxelType GetSD(int = 1) const;

  /// Get Max Intensity position around the point
  void GetMaxPosition(irtkPoint &, int = 1, int = 0) const;

  /// Get Gravity center position of a given window
  void GravityCenter(irtkPoint &, int = 1, int = 0) const;

  /// Minimum and maximum pixel values get accessor with padding
  void GetMinMaxPad(VoxelType *OUTPUT, VoxelType *OUTPUT, VoxelType pad) const;

  /// Minimum and maximum pixel values put accessor
  void PutMinMax(VoxelType, VoxelType);

  /// Function for pixel access via pointers
  VoxelType * GetPointerToVoxels(int = 0, int = 0, int = 0, int = 0) const;

  /// Funnction to convert pixel to index
  int VoxelToIndex(int, int, int, int = 0) const;

  /// Function for pixel get access
  VoxelType   Get(int, int, int, int = 0) const;

  /// Function for pixel put access
  void   Put(int, int, int, VoxelType);

  /// Function for pixel put access
  void   Put(int, int, int, int, VoxelType);

  /// Function for pixel access from via operators
  VoxelType& operator()(int, int, int, int = 0);

  /// Function for image slice get access
  irtkGenericImage GetRegion(int z, int t) const;

  /// Function for image slice get access in certain region
  irtkGenericImage GetRegion(int x1, int y1, int z1, int x2, int y2, int z2) const;

  /// Function for image slice get access in certain region
  irtkGenericImage GetRegion(int x1, int y1, int z1, int t1, int x2, int y2, int z2, int t2) const;

  /// Function for image frame get access
  irtkGenericImage GetFrame(int t) const;

  //
  // Operators for image arithmetics
  //
  
  /// Addition operator
  irtkGenericImage  operator+ (const irtkGenericImage &);

  /// Addition operator (stores result)
  irtkGenericImage& operator+=(const irtkGenericImage &);

  /// Subtraction operator
  irtkGenericImage  operator- (const irtkGenericImage &);

  /// Subtraction operator (stores result)
  irtkGenericImage& operator-=(const irtkGenericImage &);

  /// Multiplication operator
  irtkGenericImage  operator* (const irtkGenericImage &);

  /// Multiplication operator (stores result)
  irtkGenericImage& operator*=(const irtkGenericImage &);

  /// Division operator
  irtkGenericImage  operator/ (const irtkGenericImage &);

  /// Division operator (stores result)
  irtkGenericImage& operator/=(const irtkGenericImage &);

  //
  // Operators for image and Type arithmetics
  //

  /// Addition operator for type
  irtkGenericImage  operator+ (VoxelType);
  /// Addition operator for type (stores result)
  irtkGenericImage& operator+=(VoxelType);
  /// Subtraction operator for type
  irtkGenericImage  operator- (VoxelType);
  /// Subtraction operator for type (stores result)
  irtkGenericImage& operator-=(VoxelType);
  /// Multiplication operator for type
  irtkGenericImage  operator* (VoxelType);
  /// Multiplication operator for type (stores result)
  irtkGenericImage& operator*=(VoxelType);
  /// Division operator for type
  irtkGenericImage  operator/ (VoxelType);
  /// Division operator for type (stores result)
  irtkGenericImage& operator/=(VoxelType);

  //
  // Operators for image thresholding
  //

  /// Threshold operator >  (sets all values >  given value to that value)
  irtkGenericImage  operator> (VoxelType);
  /// Threshold operator >= (sets all values >= given value to that value)
  irtkGenericImage& operator>=(VoxelType);
  /// Threshold operator <  (sets all values <  given value to that value)
  irtkGenericImage  operator< (VoxelType);
  /// Threshold operator <= (sets all values <= given value to that value)
  irtkGenericImage& operator<=(VoxelType);

  /// Comparison operators == (explicit negation yields != operator)
  bool operator==(const irtkGenericImage &);

  /// Comparison operator != (if _HAS_STL is defined, negate == operator)
  ///  bool operator!=(const irtkGenericImage &);

  //
  // Reflections and axis flipping
  //

  /// Reflect image around x
  void ReflectX();
  /// Reflect image around y
  void ReflectY();
  /// Reflect image around z
  void ReflectZ();
  /// Flip x and y axis
  void FlipXY(int);
  /// Flip x and z axis
  void FlipXZ(int);
  /// Flip y and z axis
  void FlipYZ(int);
  /// Flip x and t axis
  void FlipXT(int);
  /// Flip y and t axis
  void FlipYT(int);
  /// Flip z and t axis
  void FlipZT(int);

  //
  // Conversions from and to VTK
  //

#ifdef HAS_VTK

  /// Return the VTK scalar image type of an IRTK image
  int ImageToVTKScalarType();

  /// Conversion to VTK structured points
  //  void ImageToVTK(vtkStructuredPoints *);

  /// Conversion from VTK structured points
  //  void VTKToImage(vtkStructuredPoints *);

#endif

  /// Function for pixel get access as double
  double GetAsDouble(int, int, int, int = 0) const;

  /// Function for pixel put access
  void   PutAsDouble(int, int, int, double);

  /// Function for pixel put access
  void   PutAsDouble(int, int, int, int, double);

  /// Returns the name of the image class
  const char *NameOfClass();
  
  /// Function for pixel access via pointers
  void *GetScalarPointer(int = 0, int = 0, int = 0, int = 0) const;

  /// Function which returns pixel scalar type
  virtual int GetScalarType() const;
  
  /// Function which returns the minimum value the pixel can hold without overflowing
  virtual double GetScalarTypeMin() const;

  /// Function which returns the minimum value the pixel can hold without overflowing
  virtual double GetScalarTypeMax() const;


};
