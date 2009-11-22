template <typename T> extern class itkVector3D
{
 public:
  /** The x-component. */
  double _x;

  /** The y-component. */
  double _y;

  /** The z-component. */
  double _z;

  /** Constructor. */
  itkVector3D(double x = 0, double y = 0, double z = 1);

  /** Normalizes the vector. */
  void Normalize();

  /** Operator for multiplying by a scalar. */
  itkVector3D operator*(double s);

  /** Operator for adding two vectors. */
  itkVector3D operator+(const itkVector3D& v);

  /** Operator for subtraction. */
  itkVector3D operator-(const itkVector3D& v);

  /** Operator for multiplying by a scalar. */
  itkVector3D& operator*=(double s);

  /** Operator for mulityplying by a vector. */
  itkVector3D& operator*=(const itkVector3D& v);

  /** Operator for adding a vector. */
  itkVector3D& operator+=(const itkVector3D& v);

  /** Operator for subtracting a vector. */
  itkVector3D& operator-=(const itkVector3D& v);

  /** Operator for testing equality of two vectors. */
  bool operator==(const itkVector3D& v);

  /** Operator for testing non-equality of two vector. */
  bool operator!=(const itkVector3D& v);

  /** Operator for comparing sizes of vectors. */
  bool operator<(const itkVector3D& v);

  /** Operator for comparing sizes of vectors. */
  bool operator>(const itkVector3D& v);

  /** Operator for comparing sizes of vectors. */
  bool operator<=(const itkVector3D& v);

  /** Operator for comparing sizes of vectors. */
  bool operator>=(const itkVector3D& v);

  /** Operator for dividing one vector by another. */
  itkVector3D& operator/=(const itkVector3D& v);

  /** Operator for dividing one vector by another. */
  itkVector3D operator/(const itkVector3D& v);

  /** Takes the cross-product of two vectors. */
  static itkVector3D CrossProduct(const itkVector3D& v1, const itkVector3D& v2);

  /** Takes the dot-product of two vectors. */
  static double DotProduct(const itkVector3D& v1, const itkVector3D& v2);
  
%extend
{
  const char* __str__()
  {
    static char buffer[1024];
    sprintf(buffer, "%f %f %f", static_cast<double>(self->_x), static_cast<double>(self->_y), static_cast<double>(self->_z));
    return buffer;
  }
}
};
