extern class itkMatrix
{
public:
  itkMatrix();
  itkMatrix(int, int);
  itkMatrix(const itkMatrix &);
  ~itkMatrix();
  void Initialize(int, int);
  int Rows() const;
  int Cols() const;
  void   Put(int, int, double);
  double Get(int, int) const;
  itkMatrix& operator-=(const double&);
  itkMatrix& operator+=(const double&);
  itkMatrix& operator*=(const double&);
  itkMatrix& operator/=(const double&);
  itkMatrix  operator- (const double&);
  itkMatrix  operator+ (const double&);
  itkMatrix  operator* (const double&);
  itkMatrix  operator/ (const double&);
  itkMatrix& operator-=(const itkMatrix&);
  itkMatrix& operator+=(const itkMatrix&);
  itkMatrix& operator*=(const itkMatrix&);
  itkMatrix  operator- (const itkMatrix&);
  itkMatrix  operator+ (const itkMatrix&);
  itkMatrix  operator* (const itkMatrix&);
  int operator==(const itkMatrix &);
  double Norm(void) const;
  double Det() const;
  void   SVD(itkMatrix &, itkVector &, itkMatrix &) const;
  void   Ident();
  bool IsIdentity() const;
  void   Invert();
  void   Transpose();
  void   Eigenvalues(itkMatrix &, itkVector &);
  void   LeastSquaresFit(const itkVector &, itkVector &);
  void Print();
  void Read (char *);
  void Write(char *);
  void Import (char *, int, int);

%extend
{
  const char* __str__()
    {
      std::stringstream buffer;

      for (int i = 0; i < self->Rows(); i++)
      {
        for (int j = 0; j < self->Cols(); j++)
          buffer << self->Get(i, j) << " ";
        
        if (i < self->Rows() - 1)
          buffer << std::endl;
      }
      
      char* buffer2 = new char[buffer.str().length() + 1];
      strcpy(buffer2, buffer.str().c_str());
      
      return buffer2;
    }
}     
};
