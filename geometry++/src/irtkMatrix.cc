/*=========================================================================

  Library   : Image Registration Toolkit (IRTK)
  Module    : $Id$
  Copyright : Imperial College, Department of Computing
              Visual Information Processing (VIP), 2008 onwards
  Date      : $Date$
  Version   : $Revision$
  Changes   : $Author$

Copyright (c) 1999-2014 and onwards, Imperial College London
All rights reserved.
See LICENSE for details

=========================================================================*/

#include <irtkGeometry.h>

#ifdef USE_VXL
// VXL header file was loaded in irtkMatrix.h
#include <vnl/algo/vnl_determinant.h>
#include <vnl/algo/vnl_matrix_inverse.h>
#include <vnl/algo/vnl_real_eigensystem.h>
#include <vnl/algo/vnl_symmetric_eigensystem.h>
#else
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_eigen.h>
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
    cerr << "irtkMatrix::FrechetMean : Weight sum must be positive." << endl;
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
  double d;
  int i, s;

  gsl_matrix *m = gsl_matrix_alloc(_rows, _rows);
  Matrix2GSL(m);

  gsl_permutation *p = gsl_permutation_alloc(_rows);

  gsl_linalg_LU_decomp(m, p, &s);
  d = s;
  for (i = 0; i < _rows; i++) {
     d *= gsl_matrix_get(m, i, i);
  }

  gsl_permutation_free(p);
  gsl_matrix_free(m);
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
  double d;
  int i, j, s;

  // Allocate memory
  gsl_matrix *a = gsl_matrix_alloc(_rows, _rows);
  gsl_matrix *b = gsl_matrix_alloc(_rows, _rows);
  gsl_vector *v = gsl_vector_alloc(_rows);
  gsl_vector *x = gsl_vector_alloc(_rows);

  // Convert matrix to GSL format
  Matrix2GSL(a);

  // Do the LU decomposition
  gsl_permutation *p = gsl_permutation_alloc(_rows);
  gsl_linalg_LU_decomp(a, p, &s);

  d = s;
  for (j = 0; j < _rows; j++) {
    d *= gsl_matrix_get(a, j, j);
  }

  if (d == 0) {
    cerr << "irtkMatrix::Invert: Zero determinant\n";
    //exit(1);
  }

  for (j = 0; j < _rows; j++) {
    for (i = 0; i < _rows; i++) {
      gsl_vector_set(v, i, 0.0);
    }
    gsl_vector_set(v, j, 1.0);

    gsl_linalg_LU_solve(a, p, v, x);

    for (i = 0; i < _rows; i++) {
      gsl_matrix_set(b, i, j, gsl_vector_get(x, i));
    }
  }

  // Convert GSL format back
  GSL2Matrix(b);

  // Deallocate memory
  gsl_permutation_free(p);
  gsl_matrix_free(a);
  gsl_matrix_free(b);
  gsl_vector_free(v);
  gsl_vector_free(x);
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
  int i, j, s;

  // Allocate memory
  gsl_matrix *a = gsl_matrix_alloc(_rows, _rows);
  gsl_matrix *b = gsl_matrix_alloc(_rows, _rows);
  gsl_vector *v = gsl_vector_alloc(_rows);
  gsl_vector *x = gsl_vector_alloc(_rows);

  // Convert matrix to GSL format
  Matrix2GSL(a);

  // Dp the LU decomposition
  gsl_permutation *p = gsl_permutation_alloc(_rows);
  gsl_linalg_LU_decomp(a, p, &s);

  d = s;
  for (j = 0; j < _rows; j++) {
    d *= gsl_matrix_get(a, j, j);
  } 

  if (d == 0) {
    cerr << "irtkMatrix::Invert: Zero determinant\n";
    //exit(1);
  }

  for (j = 0; j < _rows; j++) {
    for (i = 0; i < _rows; i++) {
      gsl_vector_set(v, i, 0.0);
    }
    gsl_vector_set(v, j, 1.0);

    gsl_linalg_LU_solve(a, p, v, x);

    for (i = 0; i < _rows; i++) {
      gsl_matrix_set(b, i, j, gsl_vector_get(x, i) * d);
    }
  }

  // Convert GSL format back
  GSL2Matrix(b);

  // Deallocate memory
  gsl_permutation_free(p);
  gsl_vector_free(x);
  gsl_vector_free(v);
  gsl_matrix_free(a);
  gsl_matrix_free(b);
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

  // check for diagonizable matrix (if it commutes with its conjugate transpose: A*A = AA*) - normality test equivalent
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
    gsl_matrix *GSL_u = gsl_matrix_alloc(_rows, _cols);
    gsl_matrix *GSL_v = gsl_matrix_alloc(_cols, _cols);
    gsl_vector *GSL_w = gsl_vector_alloc(_cols);
    gsl_vector *work = gsl_vector_alloc(_cols);

    // Convert matrix to GSL
    this->Matrix2GSL(GSL_u);

    // SVD
    gsl_linalg_SV_decomp(GSL_u, GSL_v, GSL_w, work);

    // Convert matrix to GSL
    E1 = irtkMatrix(_rows, _cols);
    E1.GSL2Matrix(GSL_u);
    E2 = irtkMatrix(_cols, _cols);
    E2.GSL2Matrix(GSL_v);
    e = irtkVector(_cols);
    e.GSL2Vector(GSL_w);

    // Free memory
    gsl_vector_free(GSL_w);
    gsl_vector_free(work);
    gsl_matrix_free(GSL_u);
    gsl_matrix_free(GSL_v);

    // eigenvalues are still in square root
    for (i = 0; i < _rows; i++) {
      e(i) = e(i)*e(i);
    }
    
  } else {
    // Allocate menory
    gsl_matrix *GSL_m = gsl_matrix_alloc(_rows, _rows);
    gsl_matrix *evec = gsl_matrix_alloc(_rows, _rows);
    gsl_vector *eval = gsl_vector_alloc(_rows);

    // Convert matrix to GSL format
    Matrix2GSL(GSL_m);

    if (sym) { // Jacobi
      gsl_eigen_symmv_workspace *w = gsl_eigen_symmv_alloc(_rows);
      gsl_eigen_symmv(GSL_m, eval, evec, w);
      gsl_eigen_symmv_free(w);
    } else { // Eigenvalue Decomposition
      irtkException e("Cannot calculate eigenvalues, matrix is not symmetric.");
      throw e;
    }
    // Convert GSL format back
    E1 = irtkMatrix(_rows, _rows);
    E1.GSL2Matrix(evec);
    E2 = E1;
    E2.Invert();
    e = irtkVector(_rows);
    e.GSL2Vector(eval);

    // Free memory
    gsl_matrix_free(GSL_m);
    gsl_matrix_free(evec);
    gsl_vector_free(eval);
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
  gsl_matrix *GSL_u = gsl_matrix_alloc(_rows, _cols);

  // Convert matrix to GSL compatible format
  this->Matrix2GSL(GSL_u);

  gsl_matrix *GSL_v = gsl_matrix_alloc(_cols, _cols);
  gsl_vector *GSL_w = gsl_vector_alloc(_cols);
  gsl_vector *work = gsl_vector_alloc(_cols);

  // SVD
  gsl_linalg_SV_decomp(GSL_u, GSL_v, GSL_w, work);

  // Convert matrix to GSL
  u = irtkMatrix(_rows, _cols);
  u.GSL2Matrix(GSL_u);
  v = irtkMatrix(_cols, _cols);
  v.GSL2Matrix(GSL_v);
  w = irtkVector(_cols);
  w.GSL2Vector(GSL_w);

  // Free memory
  gsl_matrix_free(GSL_u);
  gsl_matrix_free(GSL_v);
  gsl_vector_free(GSL_w);
  gsl_vector_free(work);
}

void irtkMatrix::LeastSquaresFit(const irtkVector &y, irtkVector &x)
{
#ifdef USE_VXL
  cerr << "Not provided by VXL library." << endl;
#else
  int j;
  float wmax, thresh;

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

  // Allocate GSL memory
  gsl_vector *GSL_b = gsl_vector_alloc(_rows);
  gsl_vector *GSL_w = gsl_vector_alloc(_cols);
  gsl_vector *GSL_x = gsl_vector_alloc(_cols);
  gsl_vector *GSL_y = gsl_vector_alloc(_rows);
  gsl_matrix *GSL_v = gsl_matrix_alloc(_cols, _cols);
  gsl_matrix *GSL_u = gsl_matrix_alloc(_rows, _cols);

  gsl_vector *work = gsl_vector_alloc(_cols);

  // Convert matrix to GSL
  this->Matrix2GSL(GSL_u);

  // Convert vector to GSL
  y.Vector2GSL(GSL_y);

  // Calculate least squares fit via SVD
  wmax = 0.0;

  gsl_linalg_SV_decomp(GSL_u, GSL_v, GSL_w, work);

  for (j = 0; j < _cols; j++) if (gsl_vector_get(GSL_w, j) > wmax) wmax = gsl_vector_get(GSL_w, j);
  thresh = NR_TOL * wmax;
  for (j = 0; j < _cols; j++) if (gsl_vector_get(GSL_w, j) < thresh) gsl_vector_set(GSL_w, j, 0.0);

  gsl_linalg_SV_solve(GSL_u, GSL_v, GSL_w, GSL_y, GSL_x);

  // Convert GSL to vector
  x.GSL2Vector(GSL_x);

  // Deallocate GSL memory
  gsl_vector_free(GSL_b);
  gsl_vector_free(GSL_w);
  gsl_vector_free(GSL_x);
  gsl_vector_free(GSL_y);
  gsl_vector_free(work);
  gsl_matrix_free(GSL_v);
  gsl_matrix_free(GSL_u);
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

void irtkMatrix::Matrix2GSL(gsl_matrix *m) const
{
   int i, j;

   for (i = 0; i < _rows; i++) {
     for (j = 0; j < _cols; j++) {
	gsl_matrix_set(m, i, j, _matrix[j][i]);
     }
   }
}

void irtkMatrix::GSL2Matrix(gsl_matrix *m)
{
  int i, j;

  for (j = 0; j < _cols; j++) {
    for (i = 0; i < _rows; i++) {
      _matrix[j][i] = gsl_matrix_get(m, i, j);
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
