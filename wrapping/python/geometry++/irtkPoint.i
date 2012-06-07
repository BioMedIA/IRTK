extern class irtkPoint {
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
  irtkPoint();

  /// Constructor with three coordinates
  irtkPoint(double, double, double);

  /// Constructor with Point
  irtkPoint(const irtkPoint &);

  /// Constructor with Vector
  irtkPoint(const irtkVector&);

  /// Default destructor
  virtual ~irtkPoint(void);

  //
  // Operators for Point
  //

  /// Substraction operator for point
  irtkPoint& operator-=(const irtkPoint&);

  /// Addition operator for point
  irtkPoint& operator+=(const irtkPoint&);

  /// Multiplication operator for point
  irtkPoint& operator*=(const irtkPoint&);

  /// Division operator for point
  irtkPoint& operator/=(const irtkPoint&);

  /// Return result of point substraction
  irtkPoint  operator- (const irtkPoint&);

  /// Return result of point addition
  irtkPoint  operator+ (const irtkPoint&);

  /// Return result of point multiplication
  irtkPoint  operator* (const irtkPoint&);

  /// Return result of point division
  irtkPoint  operator/ (const irtkPoint&);

  // 
  // Operators for comparison
  //

  /// Comparison operator ==
  int    operator==(const irtkPoint&);

  /// Comparison operator != (if USE_STL is defined, negate == operator)
  int    operator!=(const irtkPoint&);

  /// Comparison operator <
  int    operator<(const irtkPoint&);

  /// Comparison operator >
  int    operator>(const irtkPoint&);

  //
  // Operators for double
  //

  /// Substraction of double
  irtkPoint& operator-=(double);

  /// Addition of double
  irtkPoint& operator+=(double);

  /// Multiplication with double
  irtkPoint& operator*=(double);

  /// Division by double
  irtkPoint& operator/=(double);

  // Return result of substraction of double
  irtkPoint  operator- (double);

  // Return result of addition of double
  irtkPoint  operator+ (double);

  // Return result of multiplication with double
  irtkPoint  operator* (double);

  // Return result of division by double
  irtkPoint  operator/ (double);

  //
  // Operators for Vector
  //

  /// Substraction operator for vectors
  irtkPoint& operator-=(const irtkVector&);

  /// Addition operator for vectors
  irtkPoint& operator+=(const irtkVector&);

  /// Multiplication operator for vectors (componentwise)
  irtkPoint& operator*=(const irtkVector&);

  /// Division operator for vectors (componentwise)
  irtkPoint& operator/=(const irtkVector&);

  // Return result of vector substraction
  irtkPoint  operator- (const irtkVector&);

  // Return result of vector addition
  irtkPoint  operator+ (const irtkVector&);

  // Return result of vector multiplication
  irtkPoint  operator* (const irtkVector&);

  // Return result of vector division
  irtkPoint  operator/ (const irtkVector&);

  //
  // Operators for Matrix
  //

  /// Point multiplication operator for matrices
  irtkPoint& operator*=(const irtkMatrix&);

  /// Return result from Matrix multiplication
  irtkPoint  operator* (const irtkMatrix&);

  //
  // Distance methods 
  //

  /// Distance from origin
  double  Distance(void) const;

  /// Distance from point
  double  Distance(const irtkPoint&) const;

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
