extern class itkQuaternion
{
public:
  double _t;
  double _u;
  double _v;
  double _w;
  
  itkQuaternion();
  itkQuaternion(double t, double u, double v, double w);
  itkQuaternion(const itkQuaternion& q);
%extend
{
  itkQuaternion operator+(const itkQuaternion& a)
  {
    return *self + a;
  }
  
  itkQuaternion operator-(const itkQuaternion& a)
  {
    return *self - a;
  }
  
  itkQuaternion operator*(const itkQuaternion& a)
  {
    return *self*a;
  }
  
  itkQuaternion operator*(double b)
  {
    return *self*b;
  }
  
  const char* __str__()
  {
    static char buffer[1024];
    
    sprintf(buffer, "t = %.12f, u = %.12f, v = %.12f, w = %.12f", self->_t, self->_u, self->_v, self->_w);
    
    return buffer;
  }
}
  itkQuaternion Conjugate() const;
  itkQuaternion Inverse() const;
  void Normalize();
  double Length() const;
  itkMatrix ComputeRotationMatrix() const;
  static itkMatrix ComputeRotationMatrix(double alpha, double u, double v, double w);
  static itkQuaternion ComputeRotationQuaternion(double alpha, double u, double v, double w);
};