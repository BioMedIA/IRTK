template <typename T> extern class irtkVector3D
{
 public:
  /** The x-component. */
  double _x;

  /** The y-component. */
  double _y;

  /** The z-component. */
  double _z;

  /** Constructor. */
  irtkVector3D(double x, double y, double z);

  /** Normalizes the vector. */
  void Normalize();

  /** Operator for multiplying by a scalar. */
  irtkVector3D operator*(double s);

  /** Operator for adding two vectors. */
  irtkVector3D operator+(const irtkVector3D& v);

  /** Operator for subtraction. */
  irtkVector3D operator-(const irtkVector3D& v);

  /** Operator for multiplying by a scalar. */
  irtkVector3D& operator*=(double s);

  /** Operator for mulityplying by a vector. */
  irtkVector3D& operator*=(const irtkVector3D& v);

  /** Operator for adding a vector. */
  irtkVector3D& operator+=(const irtkVector3D& v);

  /** Operator for subtracting a vector. */
  irtkVector3D& operator-=(const irtkVector3D& v);

  /** Operator for testing equality of two vectors. */
  bool operator==(const irtkVector3D& v);

  /** Operator for testing non-equality of two vector. */
  bool operator!=(const irtkVector3D& v);

  /** Operator for comparing sizes of vectors. */
  bool operator<(const irtkVector3D& v);

  /** Operator for comparing sizes of vectors. */
  bool operator>(const irtkVector3D& v);

  /** Operator for comparing sizes of vectors. */
  bool operator<=(const irtkVector3D& v);

  /** Operator for comparing sizes of vectors. */
  bool operator>=(const irtkVector3D& v);

  /** Operator for dividing one vector by another. */
  irtkVector3D& operator/=(const irtkVector3D& v);

  /** Operator for dividing one vector by another. */
  irtkVector3D operator/(const irtkVector3D& v);

  /** Takes the cross-product of two vectors. */
  static irtkVector3D CrossProduct(const irtkVector3D& v1, const irtkVector3D& v2);

  /** Takes the dot-product of two vectors. */
  static double DotProduct(const irtkVector3D& v1, const irtkVector3D& v2);
  
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
