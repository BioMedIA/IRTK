/*=========================================================================

  Library   : Image Registration Toolkit (IRTK)
  Module    : $Id$
  Copyright : Imperial College, Department of Computing
              Visual Information Processing (VIP), 2008 onwards
  Date      : $Date$
  Version   : $Revision$
  Changes   : $Author$

=========================================================================*/

#ifndef IRTKEIGENANALYSIS_H

#define IRTKEIGENANALYSIS_H

#include <irtkImage.h>

class irtkEigenAnalysis
{
public:
  irtkEigenAnalysis (int _size);
  ~irtkEigenAnalysis ();

  // set the matrix for eigensolving
  float& Matrix (int row, int col) {
    return mat[row][col];
  }
  irtkEigenAnalysis& Matrix (float** inmat);

  // get the results of eigensolving
  float Eigenvalue (int d) {
    return diag[d];
  }
  float Eigenvector (int row, int col) {
    return mat[row][col];
  }
  const float* Eigenvalue () {
    return diag;
  }
  const float** Eigenvector () {
    return (const float**) mat;
  }

  // solve eigensystem
  void EigenStuff2 ();  // uses TriDiagonal2
  void EigenStuff3 ();  // uses TriDiagonal3
  void EigenStuff4 ();  // uses TriDiagonal4
  void EigenStuffN ();  // uses TriDiagonalN
  void EigenStuff  ();  // uses switch statement

  // solve eigensystem, use decreasing sort on eigenvalues
  void DecrSortEigenStuff2 ();
  void DecrSortEigenStuff3 ();
  void DecrSortEigenStuff4 ();
  void DecrSortEigenStuffN ();
  void DecrSortEigenStuff  ();

  // solve eigensystem, use increasing sort on eigenvalues
  void IncrSortEigenStuff2 ();
  void IncrSortEigenStuff3 ();
  void IncrSortEigenStuff4 ();
  void IncrSortEigenStuffN ();
  void IncrSortEigenStuff  ();

  // solve eigensystem, use decreasing sort on eigenvalues
  void DecrMagSortEigenStuff2 ();
  void DecrMagSortEigenStuff3 ();
  void DecrMagSortEigenStuff4 ();
  void DecrMagSortEigenStuffN ();
  void DecrMagSortEigenStuff  ();

  // solve eigensystem, use increasing sort on eigenvalues
  void IncrMagSortEigenStuff2 ();
  void IncrMagSortEigenStuff3 ();
  void IncrMagSortEigenStuff4 ();
  void IncrMagSortEigenStuffN ();
  void IncrMagSortEigenStuff  ();

  // debugging output?
  float& Tdiag (int i) {
    return diag[i];
  }
  float& Tsubdiag (int i) {
    return subd[i];
  }
  void Reduce () {
    TridiagonalN(size,mat,diag,subd);
  }

private:
  int size;
  float** mat;
  float* diag;
  float* subd;

  // Householder reduction to tridiagonal form
  void Tridiagonal2 (float** mat, float* diag, float* subd);
  void Tridiagonal3 (float** mat, float* diag, float* subd);
  void Tridiagonal4 (float** mat, float* diag, float* subd);
  void TridiagonalN (int n, float** mat, float* diag, float* subd);

  // QL algorithm with implicit shifting, applies to tridiagonal matrices
  void QLAlgorithm (int n, float* diag, float* subd, float** mat);

  // sort eigenvalues from largest to smallest
  void DecreasingSort (int n, float* eigval, float** eigvec);

  // sort eigenvalues from smallest to largest
  void IncreasingSort (int n, float* eigval, float** eigvec);

  // sort eigenvalues from largest to smallest (in magnitude)
  void DecreasingMagnitudeSort (int n, float* eigval, float** eigvec);

  // sort eigenvalues from smallest to largest (in magnitude)
  void IncreasingMagnitudeSort (int n, float* eigval, float** eigvec);

// error handling
public:
  static int verbose;
  static unsigned error;
  static void Report (ostream& ostr);
private:
  static const unsigned invalid_size;
  static const unsigned allocation_failed;
  static const unsigned ql_exceeded;
  static const char* message[3];
  static int Number (unsigned single_error);
  static void Report (unsigned single_error);
};


class irtkEigenAnalysisD
{
public:
  irtkEigenAnalysisD (int _size);
  ~irtkEigenAnalysisD ();

  // set the matrix for eigensolving
  double& Matrix (int row, int col) {
    return mat[row][col];
  }
  irtkEigenAnalysisD& Matrix (double** inmat);

  // get the results of eigensolving
  double Eigenvalue (int d) {
    return diag[d];
  }
  double Eigenvector (int row, int col) {
    return mat[row][col];
  }
  const double* Eigenvalue () {
    return diag;
  }
  const double** Eigenvector () {
    return (const double**) mat;
  }

  // solve eigensystem
  void EigenStuff2 ();  // uses TriDiagonal2
  void EigenStuff3 ();  // uses TriDiagonal3
  void EigenStuff4 ();  // uses TriDiagonal4
  void EigenStuffN ();  // uses TriDiagonalN
  void EigenStuff  ();  // uses switch statement

  // solve eigensystem, use decreasing sort on eigenvalues
  void DecrSortEigenStuff2 ();
  void DecrSortEigenStuff3 ();
  void DecrSortEigenStuff4 ();
  void DecrSortEigenStuffN ();
  void DecrSortEigenStuff  ();

  // solve eigensystem, use increasing sort on eigenvalues
  void IncrSortEigenStuff2 ();
  void IncrSortEigenStuff3 ();
  void IncrSortEigenStuff4 ();
  void IncrSortEigenStuffN ();
  void IncrSortEigenStuff  ();

  // solve eigensystem, use decreasing sort on eigenvalues
  void DecrMagSortEigenStuff2 ();
  void DecrMagSortEigenStuff3 ();
  void DecrMagSortEigenStuff4 ();
  void DecrMagSortEigenStuffN ();
  void DecrMagSortEigenStuff  ();

  // solve eigensystem, use increasing sort on eigenvalues
  void IncrMagSortEigenStuff2 ();
  void IncrMagSortEigenStuff3 ();
  void IncrMagSortEigenStuff4 ();
  void IncrMagSortEigenStuffN ();
  void IncrMagSortEigenStuff  ();

  // debugging output?
  double& Tdiag (int i) {
    return diag[i];
  }
  double& Tsubdiag (int i) {
    return subd[i];
  }
  void Reduce () {
    TridiagonalN(size,mat,diag,subd);
  }

private:
  int size;
  double** mat;
  double* diag;
  double* subd;

  // Householder reduction to tridiagonal form
  void Tridiagonal2 (double** mat, double* diag, double* subd);
  void Tridiagonal3 (double** mat, double* diag, double* subd);
  void Tridiagonal4 (double** mat, double* diag, double* subd);
  void TridiagonalN (int n, double** mat, double* diag, double* subd);

  // QL algorithm with implicit shifting, applies to tridiagonal matrices
  void QLAlgorithm (int n, double* diag, double* subd, double** mat);

  // sort eigenvalues from largest to smallest
  void DecreasingSort (int n, double* eigval, double** eigvec);

  // sort eigenvalues from smallest to largest
  void IncreasingSort (int n, double* eigval, double** eigvec);

  // sort eigenvalues from largest to smallest (in magnitude)
  void DecreasingMagnitudeSort (int n, double* eigval, double** eigvec);

  // sort eigenvalues from smallest to largest (in magnitude)
  void IncreasingMagnitudeSort (int n, double* eigval, double** eigvec);
// error handling
public:
  static int verbose;
  static unsigned error;
  static void Report (ostream& ostr);
private:
  static const unsigned invalid_size;
  static const unsigned allocation_failed;
  static const unsigned ql_exceeded;
  static const char* message[3];
  static int Number (unsigned single_error);
  static void Report (unsigned single_error);
};
#endif
