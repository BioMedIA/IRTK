extern class itkPoint {
public:

  /// x coordinate of Point
  double _x;

  /// y coordinate of Point
  double _y;

  /// z coordinate of Point
  double _z;

  //
  // Constructors and destructor 
  //

  /// Constructor
  itkPoint();

  /// Constructor with three coordinates
  itkPoint(double, double, double);

  /// Constructor with Point
  itkPoint(const itkPoint &);

  /// Constructor with Vector
  itkPoint(const itkVector&);

  /// Default destructor
  virtual ~itkPoint(void);

  //
  // Operators for Point
  //

  /// Substraction operator for point
  itkPoint& operator-=(const itkPoint&);

  /// Addition operator for point
  itkPoint& operator+=(const itkPoint&);

  /// Multiplication operator for point
  itkPoint& operator*=(const itkPoint&);

  /// Division operator for point
  itkPoint& operator/=(const itkPoint&);

  /// Return result of point substraction
  itkPoint  operator- (const itkPoint&);

  /// Return result of point addition
  itkPoint  operator+ (const itkPoint&);

  /// Return result of point multiplication
  itkPoint  operator* (const itkPoint&);

  /// Return result of point division
  itkPoint  operator/ (const itkPoint&);

  // 
  // Operators for comparison
  //

  /// Comparison operator ==
  int    operator==(const itkPoint&);

  /// Comparison operator != (if USE_STL is defined, negate == operator)
  int    operator!=(const itkPoint&);

  /// Comparison operator <
  int    operator<(const itkPoint&);

  /// Comparison operator >
  int    operator>(const itkPoint&);

  //
  // Operators for double
  //

  /// Substraction of double
  itkPoint& operator-=(double);

  /// Addition of double
  itkPoint& operator+=(double);

  /// Multiplication with double
  itkPoint& operator*=(double);

  /// Division by double
  itkPoint& operator/=(double);

  // Return result of substraction of double
  itkPoint  operator- (double);

  // Return result of addition of double
  itkPoint  operator+ (double);

  // Return result of multiplication with double
  itkPoint  operator* (double);

  // Return result of division by double
  itkPoint  operator/ (double);

  //
  // Operators for Vector
  //

  /// Substraction operator for vectors
  itkPoint& operator-=(const itkVector&);

  /// Addition operator for vectors
  itkPoint& operator+=(const itkVector&);

  /// Multiplication operator for vectors (componentwise)
  itkPoint& operator*=(const itkVector&);

  /// Division operator for vectors (componentwise)
  itkPoint& operator/=(const itkVector&);

  // Return result of vector substraction
  itkPoint  operator- (const itkVector&);

  // Return result of vector addition
  itkPoint  operator+ (const itkVector&);

  // Return result of vector multiplication
  itkPoint  operator* (const itkVector&);

  // Return result of vector division
  itkPoint  operator/ (const itkVector&);

  //
  // Operators for Matrix
  //

  /// Point multiplication operator for matrices
  itkPoint& operator*=(const itkMatrix&);

  /// Return result from Matrix multiplication
  itkPoint  operator* (const itkMatrix&);

  //
  // Distance methods 
  //

  /// Distance from origin
  double  Distance(void) const;

  /// Distance from point
  double  Distance(const itkPoint&) const;

%extend
{
  const char* __str__()
    {
      static char buffer[1024];
      sprintf(buffer, "%f %f %f", self->_x, self->_y, self->_z);
      return buffer;
    }
}
};
