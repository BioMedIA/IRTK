extern class itkVector {
public:

  /// Default constructor
  itkVector();

  /// Constructor for given row dimensions
  itkVector(int);

  /// Copy constructor
  itkVector(const itkVector &);

  /// Destructor
  ~itkVector();

  /// Intialize matrix with number of rows
  void Initialize(int);

  //
  // Vector access functions
  // 

  /// Returns number of rows
  int Rows() const;

  /// Puts vector value
  void   Put(int, double);

  /// Gets vector value
  double Get(int) const;

  //
  // Vector operators for doubles
  //
  
  /// Subtraction of a double
  itkVector& operator-=(double);
  
  /// Addition of a double
  itkVector& operator+=(double);
  
  /// Multiplication with a double
  itkVector& operator*=(double);

  /// Division by a double
  itkVector& operator/=(double);
 
  /// Return result of subtraction of a double
  itkVector  operator- (double);
 
  /// Return result of addition of a double
  itkVector  operator+ (double);
 
  /// Return result of multiplication with a double
  itkVector  operator* (double);
 
  /// Return result of division by a double
  itkVector  operator/ (double);

  //
  // Vector operators for vectors
  //

  /// Comparison operator ==
  int operator==(const itkVector &);

  /// Vector subtraction operator
  itkVector& operator-=(const itkVector&);

  /// Vector addition operator
  itkVector& operator+=(const itkVector&);

  /// Vector componentwise multiplication operator (no scalar nor cross product)
  itkVector& operator*=(const itkVector&);

  /// Vector componentwise division operator
  itkVector& operator/=(const itkVector&);

  /// Return result for vector subtraction
  itkVector  operator- (const itkVector&);

  /// Return result for vector addition
  itkVector  operator+ (const itkVector&);

  /// Return result for componentwise vector multiplication (no scalar nor cross product)
  itkVector  operator* (const itkVector&);

  /// Return result for componentwise vector division
  itkVector  operator/ (const itkVector&);

  /// Scalar/dot product
  double ScalarProduct(const itkVector&);

  /// Vector/cross product
  itkVector CrossProduct(const itkVector&);

  /// Returns norm of a vector
  double Norm(void) const;

  /// Vector normalization
  void   Normalize(void);

  //
  // Vector in- and output
  //

  /// Print vector
  void Print();

  /// Read vector from file
  void Read (char *);

  /// Write vector to file
  void Write(char *);

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
