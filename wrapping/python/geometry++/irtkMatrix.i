extern class irtkMatrix
{
public:

  irtkMatrix();
  irtkMatrix(int, int);
  irtkMatrix(const irtkMatrix &);
  ~irtkMatrix();
  void Initialize(int, int);
  //
  // Matrix access functions
  //
  int Rows() const;
  int Cols() const;
  void   Put(int, int, double);
  double Get(int, int) const;
  //
  // Operators for matrix access
  //
  double& operator()(int, int);
  irtkMatrix  operator()(int, int, int, int);
  void    operator()(irtkMatrix &, int, int);
  //
  // Matrix operators for doubles
  //
  irtkMatrix& operator-=(const double&);
  irtkMatrix& operator+=(const double&);
  irtkMatrix& operator*=(const double&);
  irtkMatrix& operator/=(const double&);
  irtkMatrix  operator- (const double&);
  irtkMatrix  operator+ (const double&);
  irtkMatrix  operator* (const double&);
  irtkMatrix  operator/ (const double&);
  //
  // Matrix operators for matrices
  //
  irtkMatrix& operator-=(const irtkMatrix&);
  irtkMatrix& operator+=(const irtkMatrix&);
  irtkMatrix& operator*=(const irtkMatrix&);
  irtkMatrix  operator- (const irtkMatrix&);
  irtkMatrix  operator+ (const irtkMatrix&);
  irtkMatrix  operator* (const irtkMatrix&);
  irtkMatrix  operator~ (void);
  int operator==(const irtkMatrix &);
  // Matrix exponential via Pade approximation.
  // See Golub and Van Loan, Matrix Computations, Algorithm 11.3-1.
  friend irtkMatrix expm(irtkMatrix);
  // Matrix logarithm.
  friend irtkMatrix logm(irtkMatrix);
  // Matrix square root.
  friend irtkMatrix sqrtm(irtkMatrix);
  friend irtkMatrix FrechetMean(irtkMatrix *, int, int = 10);
  friend irtkMatrix FrechetMean(irtkMatrix *, double *, int, int = 10);
  //
  // Matrix functions
  //
  double Norm(void) const;
  double Trace(void) const;
  // The infinity norm is the maximum of the absolute value row sums.
  double InfinityNorm(void) const;
  double Det() const;
  void   SVD(irtkMatrix &, irtkVector &, irtkMatrix &) const;
  void   Ident();
  bool IsIdentity() const;
  void   Invert();
  void   Adjugate(double &);
  void   Transpose();
  void   Eigenvalues(irtkMatrix &, irtkVector &);
  void   LeastSquaresFit(const irtkVector &, irtkVector &);
  //
  // Matrix in- and output
  //
  void Print();
  void Read (char *);
  void Write(char *);
  void Import (char *, int, int);

  void Matrix2NR(float **) const;
  void Matrix2NR(double **) const;
  void NR2Matrix(float **);
  void NR2Matrix(double **);
  //
  // Matrix operators for vectors
  //
  irtkVector  operator* (const irtkVector&);

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
