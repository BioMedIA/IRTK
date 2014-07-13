extern class irtkVector {
public:

  irtkVector();
  irtkVector(int);
  irtkVector(const irtkVector &);
  ~irtkVector();
  void Initialize(int);
  //
  // Vector access functions
  //
  int Rows() const;
  void   Put(int, double);
  double Get(int) const;
  //
  // Operators for vector access
  //
  double &operator()(int);
  //
  // Vector operators for doubles
  //
  irtkVector& operator-=(double);
  irtkVector& operator+=(double);
  irtkVector& operator*=(double);
  irtkVector& operator/=(double);
  irtkVector  operator- (double);
  irtkVector  operator+ (double);
  irtkVector  operator* (double);
  irtkVector  operator/ (double);
  //
  // Vector operators for vectors
  //
  irtkVector& operator-=(const irtkVector&);
  irtkVector& operator+=(const irtkVector&);
  irtkVector& operator*=(const irtkVector&);
  irtkVector& operator/=(const irtkVector&);
  irtkVector  operator- (const irtkVector&);
  irtkVector  operator+ (const irtkVector&);
  irtkVector  operator* (const irtkVector&);
  irtkVector  operator/ (const irtkVector&);
  bool operator==(const irtkVector &);
  bool operator<(const irtkVector &);
  double ScalarProduct(const irtkVector&);
  irtkVector CrossProduct(const irtkVector&);
  double Norm(void) const;
  void   Normalize(void);
  //
  // Vector in- and output
  //
  void Print();
  void Read (char *);
  void Write(char *);

  void Vector2GSL(gsl_vector *) const;
  void GSL2Vector(gsl_vector *);

%extend
{
  const char* __str__()
    {
      std::stringstream buffer;

      for (int i = 0; i < self->Rows(); i++)
        buffer << self->Get(i) << " ";
      
      return buffer.str().c_str();
    }
}
};
