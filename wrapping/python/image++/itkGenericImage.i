template <class VoxelType> extern class itkGenericImage : public itkBaseImage
{
public:
  
  /// Default constructor
  itkGenericImage(void);   
                      
  /// Constructor from image file
  itkGenericImage(char *);                       

  /// Constructor for given image dimensions
  itkGenericImage(int, int, int = 1);    
  
  /// Constructor for given image dimensions and voxel dimensions
  itkGenericImage(int, int, int, double, double, double);    

  /// Copy constructor for image of type unsigned char
  itkGenericImage(const itkGenericImage<itkBytePixel> &);

  /// Copy constructor for image of type short
  itkGenericImage(const itkGenericImage<itkGreyPixel> &);

  /// Copy constructor for image of type float
  itkGenericImage(const itkGenericImage<itkRealPixel> &);

  /// Destructor
  ~itkGenericImage(void);                      

  /// Initialize an image
  void Initialize(int, int, int);

  /// Initialize baseimage with given size and dimensions
  void Initialize(int, int, int, double, double, double);

  /// Initialize an image
  void Initialize(const itkBaseImage &);

  /// Read image from file
  void Read (const char *); 

  /// Write image to file
  void Write(const char *);

  /// Minimum and maximum pixel values get accessor
  void GetMinMax(VoxelType *OUTPUT, VoxelType *OUTPUT) const;

  /// Minimum and maximum pixel values get accessor within a slice
  void GetMinMax(int,  VoxelType *OUTPUT, VoxelType *OUTPUT) const;

  /// Minimum and maximum pixel values put accessor
  void PutMinMax(VoxelType, VoxelType);

  /// Minimum and maximum pixel values put accessor within a slice
  void PutMinMax(int, VoxelType,  VoxelType);

  /// Funnction to convert pixel to index
  int VoxelToIndex(int, int, int) const;
  
  /// Function for pixel get access 
  VoxelType   Get(int, int, int = 0) const;
  
  /// Function for pixel get access as double
  double GetAsDouble(int, int, int = 0) const;
  
  /// Function for pixel put access 
  void   Put(int, int, int, VoxelType);

  /// Function for pixel put access 
  void   PutAsDouble(int, int, int, double);
  
  /// Function for image slice put access
  void     PutRegion(int z, const itkGenericImage &);

  /// Function for image slice put access in certain region
  void     PutRegion(int x1, int y1, int z1, int x2, int y2, int z2,
		     const itkGenericImage &);

  /// Function for image slice get access
  itkGenericImage GetRegion(int z) const;

  /// Function for image slice get access in certain region
  itkGenericImage GetRegion(int x1, int y1, int z1, int x2, int y2, int z2) const;

  //
  // Operators for image arithmetics
  //

  /// Addition operator
  itkGenericImage  operator+ (const itkGenericImage &);
  /// Addition operator (stores result)
  itkGenericImage& operator+=(const itkGenericImage &);
  /// Subtraction operator
  itkGenericImage  operator- (const itkGenericImage &);
  /// Subtraction operator (stores result)
  itkGenericImage& operator-=(const itkGenericImage &);
  /// Multiplication operator
  itkGenericImage  operator* (const itkGenericImage &);
  /// Multiplication operator (stores result)
  itkGenericImage& operator*=(const itkGenericImage &);
  /// Division operator
  itkGenericImage  operator/ (const itkGenericImage &);
  /// Division operator (stores result)
  itkGenericImage& operator/=(const itkGenericImage &);

  //
  // Operators for image and Type arithmetics
  //

  /// Addition operator for type
  itkGenericImage  operator+ (VoxelType);
  /// Addition operator for type (stores result)
  itkGenericImage& operator+=(VoxelType);
  /// Subtraction operator for type
  itkGenericImage  operator- (VoxelType);
  /// Subtraction operator for type (stores result)
  itkGenericImage& operator-=(VoxelType);
  /// Multiplication operator for type
  itkGenericImage  operator* (VoxelType);
  /// Multiplication operator for type (stores result)
  itkGenericImage& operator*=(VoxelType);
  /// Division operator for type
  itkGenericImage  operator/ (VoxelType);
  /// Division operator for type (stores result)
  itkGenericImage& operator/=(VoxelType);

  //
  // Operators for image thresholding
  //

  /// Threshold operator >  (sets all values >  given value to that value)
  itkGenericImage  operator> (VoxelType);
  /// Threshold operator >= (sets all values >= given value to that value)
  itkGenericImage& operator>=(VoxelType);
  /// Threshold operator <  (sets all values <  given value to that value)
  itkGenericImage  operator< (VoxelType);
  /// Threshold operator <= (sets all values <= given value to that value)
  itkGenericImage& operator<=(VoxelType);

  /// Comparison operators == (explicit negation yields != operator)
  int operator==(const itkGenericImage &);
  /// Comparison operator != (if _HAS_STL is defined, negate == operator)
  ///  Bool operator!=(const itkGenericImage &);

  /// Boolean operation for empty
  int IsEmpty() const;
  
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
  void FlipXY();
  /// Flip x and z axis
  void FlipXZ();
  /// Flip y and z axis
  void FlipYZ();
  
  /// Returns the name of the image class 
  char *NameOfClass();

};
