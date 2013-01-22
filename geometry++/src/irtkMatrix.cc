/*=========================================================================

  Library   : Image Registration Toolkit (IRTK)
  Module    : $Id$
  Copyright : Imperial College, Department of Computing
              Visual Information Processing (VIP), 2008 onwards
  Date      : $Date$
  Version   : $Revision$
  Changes   : $Author$

=========================================================================*/

#include <irtkGeometry.h>

#ifdef USE_VXL
// VXL header file was loaded in irtkMatrix.h
#include <vnl/algo/vnl_determinant.h>
#include <vnl/algo/vnl_matrix_inverse.h>
#include <vnl/algo/vnl_real_eigensystem.h>
#include <vnl/algo/vnl_symmetric_eigensystem.h>
#else
#include <nr.h>
#include <nrutil.h>
#endif

#define NR_TOL 1.0e-5

//#define JACOBI

irtkMatrix::irtkMatrix()
{
  _rows = 0;
  _cols = 0;
  _matrix = NULL;
}

irtkMatrix::irtkMatrix(int rows, int cols)
{
  int i, j;

  _rows = rows;
  _cols = cols;
  _matrix = NULL;
  if ((_rows > 0) && (_cols > 0)) _matrix = Allocate(_matrix, _rows, _cols);
  for (j = 0; j < _cols; j++) {
    for (i = 0; i < _rows; i++) {
      _matrix[j][i] = 0;
    }
  }
}

irtkMatrix::irtkMatrix(const irtkMatrix& m) : irtkObject(m)
{
  int i, j;

  _rows = m._rows;
  _cols = m._cols;
  _matrix = NULL;
  if ((_rows > 0) && (_cols > 0)) _matrix = Allocate(_matrix, _rows, _cols);
  for (j = 0; j < _cols; j++) {
    for (i = 0; i < _rows; i++) {
      _matrix[j][i] = m._matrix[j][i];
    }
  }
}

irtkMatrix::~irtkMatrix()
{
  if (_matrix != NULL) Deallocate(_matrix);
  _matrix = NULL;
  _rows = 0;
  _cols = 0;
}

void irtkMatrix::Initialize(int rows, int cols)
{
  int i, j;

  if (_matrix != NULL) Deallocate(_matrix);
  _rows = rows;
  _cols = cols;
  _matrix = NULL;
  if ((_rows > 0) && (_cols > 0)) _matrix = Allocate(_matrix, _rows, _cols);

  for (j = 0; j < _cols; j++) {
    for (i = 0; i < _rows; i++) {
      _matrix[j][i] = 0;
    }
  }
}

void irtkMatrix::operator()(irtkMatrix &m, int rows, int cols)
{
  int i, j;

  if ((m._rows + rows > _rows) || (m._cols + cols > _cols)) {
    cerr << "irtkMatrix::operator(): Invalid range" << endl;
    exit(1);
  }
  for (j = cols; j < m._cols + cols; j++) {
    for (i = rows; i < m._rows + rows; i++) {
      _matrix[j][i] = m._matrix[j - cols][i - rows];
    }
  }
}

irtkMatrix irtkMatrix::operator()(int rows1, int cols1, int rows2, int cols2)
{
  int i, j;

  if ((rows1 < 0) || (rows2 > _rows) || (cols1 < 0) || (cols2 > _cols)) {
    cerr << "irtkMatrix::operator(): Invalid range" << endl;
    exit(1);
  }

  // Create new matrix
  irtkMatrix m(rows2-rows1, cols2-cols1);

  for (j = cols1; j < cols2; j++) {
    for (i = rows1; i < rows2; i++) {
      m._matrix[j-cols1][i-rows1] = _matrix[j][i];
    }
  }

  return m;
}

irtkMatrix& irtkMatrix::operator =(const irtkMatrix& m)
{
  int i, j;

  if (_matrix != NULL) Deallocate(_matrix);
  _rows = m._rows;
  _cols = m._cols;
  _matrix = NULL;
  if ((_rows > 0) && (_cols > 0)) _matrix = Allocate(_matrix, _rows, _cols);

  for (j = 0; j < _cols; j++) {
    for (i = 0; i < _rows; i++) {
      _matrix[j][i] = m._matrix[j][i];
    }
  }
  return *this;
}

irtkMatrix& irtkMatrix::operator-=(const irtkMatrix& m)
        {
  int i, j;

  if ((_rows != m._rows) || (_cols != m._cols)) {
    cerr << "irtkMatrix::operator-=: Size mismatch" << endl;
    exit(1);
  }
  for (j = 0; j < _cols; j++) {
    for (i = 0; i < _rows; i++) {
      _matrix[j][i] -= m._matrix[j][i];
    }
  }
  return *this;
        }

irtkMatrix& irtkMatrix::operator+=(const irtkMatrix& m)
        {
  int i, j;

  if ((_rows != m._rows) || (_cols != m._cols)) {
    cerr << "irtkMatrix::operator+=: Size mismatch" << endl;
    exit(1);
  }
  for (j = 0; j < _cols; j++) {
    for (i = 0; i < _rows; i++) {
      _matrix[j][i] += m._matrix[j][i];
    }
  }
  return *this;
        }

irtkMatrix& irtkMatrix::operator*=(const irtkMatrix& m)
        {
  *this = *this * m;
  return *this;
        }

irtkMatrix  irtkMatrix::operator- (const irtkMatrix& m)
{
  int i, j;

  if ((_rows != m._rows) || (_cols != m._cols)) {
    cerr << "irtkMatrix::operator-: Size mismatch" << endl;
    exit(1);
  }
  irtkMatrix tmp(_rows, _cols);
  for (j = 0; j < _cols; j++) {
    for (i = 0; i < _rows; i++) {
      tmp._matrix[j][i] = _matrix[j][i] - m._matrix[j][i];
    }
  }
  return tmp;
}

irtkMatrix  irtkMatrix::operator+ (const irtkMatrix& m)
{
  int i, j;

  if ((_rows != m._rows) || (_cols != m._cols)) {
    cerr << "irtkMatrix::operator+: Size mismatch" << endl;
    exit(1);
  }
  irtkMatrix tmp(_rows, _cols);
  for (j = 0; j < _cols; j++) {
    for (i = 0; i < _rows; i++) {
      tmp._matrix[j][i] = _matrix[j][i] + m._matrix[j][i];
    }
  }
  return tmp;
}

irtkMatrix  irtkMatrix::operator* (const irtkMatrix& m)
        {
  int i, j, k;

  if (_cols != m._rows) {
    cerr << "irtkMatrix::operator*: Size mismatch" << endl;
    exit(1);
  }

  irtkMatrix tmp(_rows, m._cols);
  for (i = 0; i < _rows; i++) {
    for (j = 0; j < m._cols; j++) {
      tmp._matrix[j][i] = 0.0;
      for (k = 0; k < _cols; k++) {
        tmp._matrix[j][i] += _matrix[k][i] * m._matrix[j][k];
      }
    }
  }

  return tmp;
        }

irtkMatrix expm (irtkMatrix A)
{
  // Matrix exponential via Pade approximation.
  // See Golub and Van Loan, Matrix Computations, Algorithm 11.3-1.

  double norm, c;
  int rows, cols, q, e, j, k, minusOnePower;

  rows = A._rows;
  cols = A._cols;
  irtkMatrix D(rows, cols);
  irtkMatrix N(rows, cols);
  irtkMatrix X(rows, cols);

  // First A is scaled by a power of 2 so its norm is < 1/2.
  // j is the index of the power required.
  norm = A.InfinityNorm();
  e = (int) ceil(log(norm) / log(2.0));
  j = max(0, 1 + e);

  A = A / (pow(2.0, j));

  // Number of iterations (determines accuracy).
  q = 6;

  D.Ident();
  N.Ident();
  X.Ident();
  c = 1.0;
  minusOnePower = 1;

  for (k = 1; k <= q; ++k) {
    c = c * (q - k + 1) / ((double) k * (2*q - k + 1));
    X = A * X;
    N = N +  X * c;
    minusOnePower *= -1;
    D = D + X * minusOnePower * c;
  }

  D.Invert();
  X = D*N;

  for (k = 1; k <= j; ++k)
    X = X*X;

  return X;
}

irtkMatrix logm(irtkMatrix A)
{
  int k;
  int i, n, maxIterations = 100;
  double tol = 0.00001;

  if (A._rows != A._cols) {
    cerr << "irtkMatrix logm(irtkMatrix A): argument must be square." << endl;
    exit(1);
  }

  irtkVector eValues(A._rows);
  irtkMatrix eVectors(A), tmp(A);
  A.Eigenvalues(eVectors, eValues, tmp);
  for (i = 0; i < A._rows; ++i) {
    if (eValues(i) <= 0) {
      cerr << "irtkMatrix logm(irtkMatrix A): eigenvalue <= 0." << endl;
      exit(1);
    }
  }

  irtkMatrix I(A._rows, A._cols);
  irtkMatrix Z(A._rows, A._cols);
  irtkMatrix X(A._rows, A._cols);
  irtkMatrix D(A._rows, A._cols);

  I.Ident();

  D = A - I;
  n = 0;
  k = 0;
  while (D.InfinityNorm() > 0.5 && n < maxIterations) {
    A = sqrtm(A);
    ++k;
    D = A - I;
    ++n;
  }

  A = I - A;
  Z = A;
  X = A;
  i = 1;
  n = 0;
  while (Z.InfinityNorm() > tol && n < maxIterations) {
    Z = Z * A;
    ++i;
    X = X + (Z / i);
    ++n;
  }

  X = X * pow(2.0, k);
  X = X * -1.0;
  return X;

}

irtkMatrix sqrtm(irtkMatrix A)
{
  double tol = 0.0001;
  int maxIterations = 100;
  int i = 0;

  irtkMatrix X(A);
  irtkMatrix Y(A._rows, A._cols);
  irtkMatrix D(A._rows, A._cols);
  irtkMatrix invX(A._rows, A._cols);
  irtkMatrix invY(A._rows, A._cols);

  Y.Ident();

  D = X*X - A;

  while (D.InfinityNorm() > tol && i < maxIterations) {
    invX = X; invX.Invert();
    invY = Y; invY.Invert();
    X = (X + invY) * 0.5;
    Y = (Y + invX) * 0.5;
    D = X*X - A;
    ++i;
  }

  return X;
}

irtkMatrix FrechetMean (irtkMatrix *matrices, int number, int iterations)
{

  int i, n = 0;
  double normLogDeltaMu = 0;
  double tolerance = 0.0000000000001;

  irtkMatrix mu(4, 4);
  irtkMatrix muInv(4, 4);
  irtkMatrix deltaMu(4, 4);
  irtkMatrix logDeltaMu(4, 4);

  irtkMatrix deltaM(4, 4);
  irtkMatrix sumLogs(4, 4);

  mu = matrices[0];

  do {

    muInv = mu; muInv.Invert();
    // Reset sum to zero.
    sumLogs.Initialize(4, 4);

    for (i = 0; i < number; i++) {
      deltaM = muInv * matrices[i];
      sumLogs += logm( deltaM );
    }

    sumLogs /= ((double) number);
    deltaMu = expm(sumLogs);

    mu = mu * deltaMu;

    logDeltaMu = logm(deltaMu);
    normLogDeltaMu = logDeltaMu.InfinityNorm();

    n++;

  } while (normLogDeltaMu > tolerance && n < iterations);

  if (n == iterations){
    cerr << "irtkMatrix::FrecheMean : Warning, reached maximum iterations." << endl;
  }

  return mu;
}

irtkMatrix FrechetMean (irtkMatrix *matrices, double *weights, int number, int iterations)
{

  int i, n = 0;
  double normLogDeltaMu = 0;
  double tolerance = 0.0000000000001;
  double totalWeight = 0.0f;

  irtkMatrix mu(4, 4);
  irtkMatrix muInv(4, 4);
  irtkMatrix deltaMu(4, 4);
  irtkMatrix logDeltaMu(4, 4);

  irtkMatrix deltaM(4, 4);
  irtkMatrix sumLogs(4, 4);

  for (i = 0; i < number; i++){
    totalWeight += weights[i];
  }

  if (totalWeight <= 0.0){
    cerr << "itkMatrix::FrechetMean : Weight sum must be positive." << endl;
    exit(1);
  }

  mu = matrices[0];

  do {

    muInv = mu; muInv.Invert();
    // Reset sum to zero.
    sumLogs.Initialize(4, 4);

    for (i = 0; i < number; i++) {
      deltaM = muInv * matrices[i];
      sumLogs += logm( deltaM ) * weights[i];
    }

    sumLogs /= ((double) totalWeight);
    deltaMu = expm(sumLogs);

    mu = mu * deltaMu;

    logDeltaMu = logm(deltaMu);
    normLogDeltaMu = logDeltaMu.InfinityNorm();

    n++;

  } while (normLogDeltaMu > tolerance && n < iterations);

  return mu;
}

irtkMatrix irtkMatrix::operator~ (void)
        {
  irtkMatrix m;

  m = *this;
  m.Transpose();
  return m;
        }

irtkMatrix irtkMatrix::operator! (void)
        {
  irtkMatrix m;

  m = *this;
  m.Invert();
  return m;
        }

double irtkMatrix::Det() const
{
  if (_rows != _cols) {
    cerr << "irtkMatrix::Det: Must be square" << endl;
    exit(1);
  }

#ifdef USE_VXL
  vnl_matrix<float> vnl(_rows,_cols);
  Matrix2Vnl(&vnl);
  d = vnl_determinant(vnl);
#else
  double **a, d;
  int i, *index;
  a = dmatrix(1, _rows, 1, _rows);
  Matrix2NR(a);
  index = ivector(1, _rows);
  ludcmp(a, _rows, index, &d);
  for (i = 1; i <= _rows; i++) {
    d *= a[i][i];
  }
  free_dmatrix(a, 1, _rows, 1, _rows);
  free_ivector( index, 1, _rows);
#endif

  return d;
}

void irtkMatrix::Invert(void)
{
  if (_rows != _cols) {
    cerr << "irtkMatrix::Invert: Must be square" << endl;
    exit(1);
  }

#ifdef USE_VXL
  vnl_matrix<float> input(_rows,_cols);
  Matrix2Vnl(&input);
  vnl_matrix<float> output = vnl_matrix_inverse<float>(input);
  Vnl2Matrix(&output);
#else
  double **a, **b, *v, d;
  int i, j, *index;

  // Allocate memory
  a = dmatrix(1, _rows, 1, _rows);
  b = dmatrix(1, _rows, 1, _rows);
  v = dvector(1, _rows);

  // Convert matrix to NR format
  Matrix2NR(a);

  index = ivector(1, _rows);

  ludcmp(a, _rows, index, &d);

  for (j = 1; j <= _rows; j++) {
    d *= a[j][j];
  }
  if (d == 0) {
    cerr << "irtkMatrix::Invert: Zero determinant\n";
    //exit(1);
  }
  for (j = 1; j <= _rows; j++) {
    for (i = 1; i <= _rows; i++) {
      v[i] = 0.0;
    }
    v[j] = 1.0;

    lubksb(a, _rows, index, v);

    for (i = 1; i <= _rows; i++) {
      b[i][j] = v[i];
    }
  }

  // Convert NR format back
  NR2Matrix(b);

  // Deallocate memory
  free_dmatrix(a, 1, _rows, 1, _rows);
  free_dmatrix(b, 1, _rows, 1, _rows);
  free_dvector(v, 1, _rows);
  free_ivector( index, 1, _rows);
#endif
}

void irtkMatrix::Adjugate(double &d)
{
  if (_rows != _cols) {
    cerr << "irtkMatrix::Adjugate: Must be square" << endl;
    exit(1);
  }

#ifdef USE_VXL
  vnl_matrix<float> input(_rows,_cols);
  Matrix2Vnl(&input);
  vnl_matrix<float> output = vnl_matrix_inverse<float>(input);
  Vnl2Matrix(&output);
#else
  double **a, **b, *v;
  int i, j, *index;

  // Allocate memory
  a = dmatrix(1, _rows, 1, _rows);
  b = dmatrix(1, _rows, 1, _rows);
  v = dvector(1, _rows);

  // Convert matrix to NR format
  Matrix2NR(a);

  index = ivector(1, _rows);

  ludcmp(a, _rows, index, &d);

  for (j = 1; j <= _rows; j++) {
    d *= a[j][j];
  }
  if (d == 0) {
    cerr << "irtkMatrix::Invert: Zero determinant\n";
    //exit(1);
  }
  for (j = 1; j <= _rows; j++) {
    for (i = 1; i <= _rows; i++) {
      v[i] = 0.0;
    }
    v[j] = 1.0;

    lubksb(a, _rows, index, v);

    for (i = 1; i <= _rows; i++) {
      b[i][j] = v[i]*d;
    }
  }

  // Convert NR format back
  NR2Matrix(b);

  // Deallocate memory
  free_dmatrix(a, 1, _rows, 1, _rows);
  free_dmatrix(b, 1, _rows, 1, _rows);
  free_dvector(v, 1, _rows);
  free_ivector( index, 1, _rows);
#endif
}

void irtkMatrix::Transpose(void)
{
  int i, j;

  if (_rows == _cols) {
    double tmp;
    for (i = 0; i < _rows; i++) {
      for (j = 0; j < i; j++) {
        tmp=_matrix[i][j];
        _matrix[i][j] = _matrix[j][i];
        _matrix[j][i] = tmp;
      }
    }
  } else {
    irtkMatrix tmp(_rows, _cols);
    for (j = 0; j < _cols; j++) {
      for (i = 0; i < _rows; i++) {
        tmp._matrix[j][i] = _matrix[j][i];
      }
    }
    Deallocate(_matrix);
    _matrix = Allocate(_matrix, tmp._cols, tmp._rows);
    _rows = tmp._cols;
    _cols = tmp._rows;
    for (j = 0; j < _cols; j++) {
      for (i = 0; i < _rows; i++) {
        _matrix[j][i] = tmp._matrix[i][j];
      }
    }
  }
}

void irtkMatrix::Eigenvalues(irtkMatrix &E1, irtkVector &e, irtkMatrix &E2)
{
  int i, j, sym = true, diag = true, square = true;

#ifdef USE_VXL

  vnl_matrix<float> M(_rows,_cols);
  Matrix2Vnl(&M);

  vnl_svd<float> svd(M);

  E1 = irtkMatrix(_rows, _cols);
  E1.Vnl2Matrix(&svd.U());
  E2 = irtkMatrix(_cols, _cols);
  E2.Vnl2Matrix(&svd.V());
  e = irtkVector(_cols);
  e.Vnl2Vector(&svd.W());

  // eigenvalues are still in square root
  for (i = 0; i < _rows; i++) {
    e(i) = e(i)*e(i);
  }
#else
  // check for square matrix
  if (_rows != _cols) {
    square = false;
  }

  // check for symmetric matrix
  if (square) {
    for (i = 0; i < _rows; i++) {
      for (j = i+1; j < _cols; j++) {
        if (abs(_matrix[i][j] - _matrix[j][i]) > 0.001) {
          sym = false;
        }
      }
    }
  } else {
    sym = false;
    diag = false;
  }

  // check for diagonizable matrix (if it commutes with its conjugate transpose: A*A = AA*)
  // only needed when matrix not symmetric (every symmetric matrix is diagonizable)
  if (!sym && square) {
    irtkMatrix AT = *this; AT.Transpose();
    irtkMatrix ATA = AT * *this;
    irtkMatrix AAT = *this * AT;
    for (i = 0; i < _rows; i++) {
      for (j = 0; j < _cols; j++) {
        if (abs(AAT(i, j) - ATA(i, j)) > 0.001) {
          diag = false;
        }
      }
    }
  }

  // compute eigenvalues with algorithm, depending on matrix properties (symmetry, diagonizability)
  if (!sym && !diag) { // SVD
    double **NR_u, **NR_v, *NR_w;

    NR_w = dvector(1, _cols);
    NR_u = dmatrix(1, _rows, 1, _cols);
    NR_v = dmatrix(1, _cols, 1, _cols);

    // Convert matrix to NR
    this->Matrix2NR(NR_u);

    // SVD
    svdcmp(NR_u, _rows, _cols, NR_w, NR_v);

    // Convert matrix to NR
    E1 = irtkMatrix(_rows, _cols);
    E1.NR2Matrix(NR_u);
    E2 = irtkMatrix(_cols, _cols);
    E2.NR2Matrix(NR_v);
    e = irtkVector(_cols);
    e.NR2Vector(NR_w);

    // Free memory
    free_dvector(NR_w, 1, _cols);
    free_dmatrix(NR_u, 1, _rows, 1, _cols);
    free_dmatrix(NR_v, 1, _cols, 1, _cols);

    // eigenvalues are still in square root
    for (i = 0; i < _rows; i++) {
      e(i) = e(i)*e(i);
    }

  } else {
    float *eigen_value, **eigen_vector, **m;

    // Allocate menory
    m            = ::matrix(1, _rows, 1, _rows);
    eigen_vector = ::matrix(1, _rows, 1, _rows);
    eigen_value  = ::vector(1, _rows);

    // Convert matrix to NR format
    Matrix2NR(m);

    if (sym) { // Jacobi
      int dummy;

      jacobi(m, _rows, eigen_value, eigen_vector, &dummy);
    } else { // Eigenvalue Decomposition
      float *dummy = ::vector(1, _rows);

      tred2(m, _rows, eigen_value, dummy);
      tqli(eigen_value, dummy, _rows, m);

      for (i = 1; i <= _rows; i++) {
        for (j = 1; j <= _rows; j++) {
          eigen_vector[i][j] = m[i][j];
        }
      }
      free_vector(dummy, 1, _rows);
    }
    // Convert NR format back
    E1 = irtkMatrix(_rows, _rows);
    E1.NR2Matrix(eigen_vector);
    E2 = E1;
    E2.Invert();
    e = irtkVector(_rows);
    e.NR2Vector(eigen_value);

    // Free memory
    free_matrix(m, 1, _rows, 1, _rows);
    free_matrix(eigen_vector, 1, _rows, 1, _rows);
    free_vector(eigen_value, 1, _rows);
  }

  // Sort eigenvectors only by absolute values
  //  for (i = 1; i <= _rows; i++) {
  //    eigen_value[i] = fabs(eigen_value[i]);
  //  }

  // Sort eigenvectors
  //  eigsrt(eigen_value, eigen_vector, _rows);

#endif
}

void irtkMatrix::SVD(irtkMatrix &u, irtkVector &w, irtkMatrix &v)
{
  int i;

  this->Eigenvalues(u, w, v);

  // eigenvalues are not in square root
  for (i = 0; i < _rows; i++) {
    w(i) = sqrt(w(i));
  }
}

void irtkMatrix::LeastSquaresFit(const irtkVector &y, irtkVector &x)
{
#ifdef USE_VXL
  cerr << "Not provided by VXL library." << endl;
#else
  int j;
  float wmax, thresh, *NR_b, *NR_w, *NR_x, *NR_y, **NR_u, **NR_v;

  // nmatrix should be rows
  // ma should be cols

  if (y.Rows() != _rows) {
    cerr << "irtkMatrix::LeastSquaresFit: Y has wrong dimensions" << endl;
    exit(1);
  }

  if (x.Rows() != _cols) {
    cerr << "irtkMatrix::LeastSquaresFit: X has wrong dimensions" << endl;
    exit(1);
  }

  // Allocate NR memory
  NR_b = ::vector(1, _rows);
  NR_w = ::vector(1, _cols);
  NR_x = ::vector(1, _cols);
  NR_y = ::vector(1, _rows);
  NR_v = ::matrix(1, _cols, 1, _cols);
  NR_u = ::matrix(1, _rows, 1, _cols);

  // Convert matrix to NR
  this->Matrix2NR(NR_u);

  // Convert vector to NR
  y.Vector2NR(NR_y);

  // Calculate least squares fit via SVD
  wmax = 0.0;

  svdcmp(NR_u, _rows, _cols, NR_w, NR_v);

  for (j = 1; j <= _cols; j++) if (NR_w[j] > wmax) wmax = NR_w[j];
  thresh = NR_TOL * wmax;
  for (j = 1; j <= _cols; j++) if (NR_w[j] < thresh) NR_w[j] = 0.0;

  svbksb(NR_u, NR_w, NR_v, _rows, _cols, NR_y, NR_x);

  // Convert NR to vector
  x.NR2Vector(NR_x);

  // Deallocate NR memory
  free_vector(NR_b, 1, _rows);
  free_vector(NR_w, 1, _cols);
  free_vector(NR_x, 1, _cols);
  free_vector(NR_y, 1, _rows);
  free_matrix(NR_v, 1, _cols, 1 ,_cols);
  free_matrix(NR_u, 1, _rows, 1, _cols);
#endif
}

void irtkMatrix::Ident()
{
  int i, j;

  for (j = 0; j < _cols; j++) {
    for (i = 0; i < _rows; i++) {
      if (i == j) {
        _matrix[j][i] = 1;
      } else {
        _matrix[j][i] = 0;
      }
    }
  }
}

bool irtkMatrix::IsIdentity() const
{
  if (_rows == 0 || _cols == 0 || _rows != _cols)
    return false;

  for (int j = 0; j < _cols; j++)
    for (int i = 0; i < _rows; i++)
      if (i == j && _matrix[j][i] != 1)
        return false;
      else if (i != j && _matrix[j][i] != 0)
        return false;

  return true;
}

ostream& operator<< (ostream& os, const irtkMatrix &m)
{
  int index, i, j;

  // Write header in ascii
  os << "irtkMatrix " << m._rows << " x " << m._cols << endl;

  // Allocate temporary memory
  double *data  = new double [m._rows*m._cols];

  // Convert data
  index = 0;
  for (j = 0; j < m._cols; j++) {
    for (i = 0; i < m._rows; i++) {
      data[index] = m._matrix[j][i];
      index++;
    }
  }

#ifndef WORDS_BIGENDIAN
  swap64((char *)data, (char *)data, m._rows*m._cols);
#endif

  // Write binary data
  os.write((char *)data, m._rows*m._cols*sizeof(double));

  // Free temporary memory
  delete [] data;

  return os;
}

istream& operator>> (istream& is, irtkMatrix &m)
{
  int index, i, j, cols, rows;
  char buffer[255];

  // Read header
  is >> buffer;
  if (strcmp(buffer, "irtkMatrix") != 0) {
    cerr << "irtkMatrix: Can't read file " << buffer << endl;
    exit(1);
  }

  // Read size
  is >> rows;
  is >> buffer;
  is >> cols;

  // Allocate matrix
  m = irtkMatrix(rows, cols);

  // Read header, skip comments
  is.get(buffer, 255);
  is.clear();
  is.seekg(1, ios::cur);

  // Allocate temporary memory
  double *data  = new double [m._rows*m._cols];

  // Read binary data
  is.read((char *)data, m._rows*m._cols*sizeof(double));

#ifndef WORDS_BIGENDIAN
  swap64((char *)data, (char *)data, m._rows*m._cols);
#endif

  // Convert data
  index = 0;
  for (j = 0; j < m._cols; j++) {
    for (i = 0; i < m._rows; i++) {
      m._matrix[j][i] = data[index];
      index++;
    }
  }

  // Free temporary memory
  delete []data;

  return is;
}

void irtkMatrix::Print()
{
  int i, j;

  cout << "irtkMatrix " << _rows << " x " << _cols << endl;
  cout.setf(ios::right);
  cout.setf(ios::fixed);
  cout.precision(4);
  for (i = 0; i < _rows; i++) {
    for (j = 0; j < _cols; j++) {
      cout << setw(15) << _matrix[j][i] << " ";
    }
    cout << endl;
  }
  cout.precision(6);
  cout.unsetf(ios::right);
  cout.unsetf(ios::fixed);
}

void irtkMatrix::Read(char *filename)
{
  // Open file stream
  ifstream from(filename, ios::in | ios::binary);

  // Check whether file opened ok
  if (!from) {
    cerr << "irtkMatrix::Read: Can't open file " << filename << endl;
    exit(1);
  }

  // Read matrix
  from >> *this;
}

void irtkMatrix::Write(char *filename)
{
  // Open file stream
  ofstream to(filename, ios::out | ios::binary);

  // Check whether file opened ok
  if (!to) {
    cerr << "irtkMatrix::Write: Can't open file " << filename << endl;
    exit(1);
  }

  // Write matrix
  to << *this;
}

void irtkMatrix::Import(char *filename, int rows, int cols)
{
  int i, j;

  // Open file stream
  ifstream from(filename);

  // Check whether file opened ok
  if (!from) {
    cerr << "irtkMatrix::Read: Can't open file " << filename << endl;
    exit(1);
  }

  // Initialize matrix
  this->Initialize(rows, cols);

  // Read matrix
  for (i = 0; i < _rows; i++) {
    for (j = 0; j < _cols; j++) {
      from >> _matrix[j][i];
    }
  }
}

#ifdef USE_VXL

template <class T>
void irtkMatrix::Matrix2Vnl(vnl_matrix<T> *m) const
{
  unsigned r, c;

  for (c = 0; c < (unsigned) _cols; c++) {
    for (r = 0; r < (unsigned) _rows; r++) {
      (*m)(r,c) = (T) _matrix[c][r];
    }
  }
}

template <class T>
void irtkMatrix::Vnl2Matrix(vnl_matrix<T> *m)
{
  unsigned r, c;

  for (c = 0; c < (unsigned) _cols; c++) {
    for (r = 0; r < (unsigned) _rows; r++) {
      _matrix[c][r] = (float) (*m)(r,c);
    }
  }
}

#else

void irtkMatrix::Matrix2NR(float **m) const
{
  int i, j;

  for (j = 0; j < _cols; j++) {
    for (i = 0; i < _rows; i++) {
      m[i+1][j+1] = _matrix[j][i];
    }
  }
}

void irtkMatrix::NR2Matrix(float **m)
{
  int i, j;

  for (j = 0; j < _cols; j++) {
    for (i = 0; i < _rows; i++) {
      _matrix[j][i] = m[i+1][j+1];
    }
  }
}

void irtkMatrix::Matrix2NR(double **m) const
{
  int i, j;

  for (j = 0; j < _cols; j++) {
    for (i = 0; i < _rows; i++) {
      m[i+1][j+1] = _matrix[j][i];
    }
  }
}

void irtkMatrix::NR2Matrix(double **m)
{
  int i, j;

  for (j = 0; j < _cols; j++) {
    for (i = 0; i < _rows; i++) {
      _matrix[j][i] = m[i+1][j+1];
    }
  }
}


irtkVector  irtkMatrix::operator* (const irtkVector& v)
        {
  int i, j;

  if (_cols != v.Rows()) {
    cerr << "irtkMatrix::operator*(irtkVector): Size mismatch" << endl;
    exit(1);
  }

  irtkVector result(_rows);
  double tmp;
  for (i = 0; i < _rows; i++) {
    tmp = 0.0;
    for (j = 0; j < _cols; j++) {
      tmp += _matrix[j][i] * v.Get(j);
    }
    result.Put(i,tmp);
  }

  return result;
        }

#endif
