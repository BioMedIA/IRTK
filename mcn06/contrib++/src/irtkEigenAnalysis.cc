/*=========================================================================

  Library   : Image Registration Toolkit (IRTK)
  Module    : $Id$
  Copyright : Imperial College, Department of Computing
              Visual Information Processing (VIP), 2008 onwards
  Date      : $Date$
  Version   : $Revision$
  Changes   : $Author$

=========================================================================*/

#include<irtkEigenAnalysis.h>

//===========================================================================
// error handling
int irtkEigenAnalysis::verbose = 0;
unsigned irtkEigenAnalysis::error = 0;
const unsigned irtkEigenAnalysis::invalid_size      = 0x00000001;
const unsigned irtkEigenAnalysis::allocation_failed = 0x00000002;
const unsigned irtkEigenAnalysis::ql_exceeded       = 0x00000004;
const char* irtkEigenAnalysis::message[3] = {
  "invalid matrix size",
  "allocation failed",
  "QL algorithm - exceeded maximum iterations"
};
//---------------------------------------------------------------------------
irtkEigenAnalysis::
irtkEigenAnalysis (int _size)
{
  if ( (size = _size) <= 1 ) {
    Report(invalid_size);
    return;
  }
  if ( (mat = new float*[size]) == 0 ) {
    Report(allocation_failed);
    return;
  }
  for (int d = 0; d < size; d++)
    if ( (mat[d] = new float[size]) == 0 ) {
      Report(allocation_failed);
      return;
    }
  if ( (diag = new float[size]) == 0 ) {
    Report(allocation_failed);
    return;
  }
  if ( (subd = new float[size]) == 0 ) {
    Report(allocation_failed);
    return;
  }
}
//---------------------------------------------------------------------------
irtkEigenAnalysis::
~irtkEigenAnalysis ()
{
  delete[] subd;
  delete[] diag;
  for (int d = 0; d < size; d++)
    delete[] mat[d];
  delete[] mat;
}


//---------------------------------------------------------------------------
void irtkEigenAnalysis::
Tridiagonal2 (float** mat, float* diag, float* subd)
{
  // matrix is already tridiagonal

  diag[0] = mat[0][0];
  diag[1] = mat[1][1];
  subd[0] = mat[0][1];
  subd[1] = 0;
  mat[0][0] = 1;  mat[0][1] = 0;
  mat[1][0] = 0;  mat[1][1] = 1;
}
//---------------------------------------------------------------------------
void irtkEigenAnalysis::
Tridiagonal3 (float** mat, float* diag, float* subd)
{
  float a = mat[0][0], b = mat[0][1], c = mat[0][2],
                                          d = mat[1][1], e = mat[1][2],
                                                             f = mat[2][2];

  diag[0] = a;
  subd[2] = 0;
  if ( c != 0 ) {
    float ell = float(sqrt(b*b+c*c));
    b /= ell;
    c /= ell;
    float q = 2*b*e+c*(f-d);
    diag[1] = d+c*q;
    diag[2] = f-c*q;
    subd[0] = ell;
    subd[1] = e-b*q;
    mat[0][0] = 1; mat[0][1] = 0; mat[0][2] = 0;
    mat[1][0] = 0; mat[1][1] = b; mat[1][2] = c;
    mat[2][0] = 0; mat[2][1] = c; mat[2][2] = -b;
  } else {
    diag[1] = d;
    diag[2] = f;
    subd[0] = b;
    subd[1] = e;
    mat[0][0] = 1; mat[0][1] = 0; mat[0][2] = 0;
    mat[1][0] = 0; mat[1][1] = 1; mat[1][2] = 0;
    mat[2][0] = 0; mat[2][1] = 0; mat[2][2] = 1;
  }
}
//---------------------------------------------------------------------------
void irtkEigenAnalysis::
Tridiagonal4 (float** mat, float* diag, float* subd)
{
  // save matrix M
  float
  a = mat[0][0], b = mat[0][1], c = mat[0][2], d = mat[0][3],
                                    e = mat[1][1], f = mat[1][2], g = mat[1][3],
                                                                      h = mat[2][2], i = mat[2][3],
                                                                                         j = mat[3][3];

  diag[0] = a;
  subd[3] = 0;

  mat[0][0] = 1; mat[0][1] = 0; mat[0][2] = 0; mat[0][3] = 0;
  mat[1][0] = 0;
  mat[2][0] = 0;
  mat[3][0] = 0;

  if ( c != 0 || d != 0 ) {
    float q11, q12, q13;
    float q21, q22, q23;
    float q31, q32, q33;

    // build column Q1
    float len = float(sqrt(b*b+c*c+d*d));
    q11 = b/len;
    q21 = c/len;
    q31 = d/len;

    subd[0] = len;

    // compute S*Q1
    float v0 = e*q11+f*q21+g*q31;
    float v1 = f*q11+h*q21+i*q31;
    float v2 = g*q11+i*q21+j*q31;

    diag[1] = q11*v0+q21*v1+q31*v2;

    // build column Q3 = Q1x(S*Q1)
    q13 = q21*v2-q31*v1;
    q23 = q31*v0-q11*v2;
    q33 = q11*v1-q21*v0;
    len = float(sqrt(q13*q13+q23*q23+q33*q33));
    if ( len > 0 ) {
      q13 /= len;
      q23 /= len;
      q33 /= len;

      // build column Q2 = Q3xQ1
      q12 = q23*q31-q33*q21;
      q22 = q33*q11-q13*q31;
      q32 = q13*q21-q23*q11;

      v0 = q12*e+q22*f+q32*g;
      v1 = q12*f+q22*h+q32*i;
      v2 = q12*g+q22*i+q32*j;
      subd[1] = q11*v0+q21*v1+q31*v2;
      diag[2] = q12*v0+q22*v1+q32*v2;
      subd[2] = q13*v0+q23*v1+q33*v2;

      v0 = q13*e+q23*f+q33*g;
      v1 = q13*f+q23*h+q33*i;
      v2 = q13*g+q23*i+q33*j;
      diag[3] = q13*v0+q23*v1+q33*v2;
    } else { // S*Q1 parallel to Q1, choose any valid Q2 and Q3
      subd[1] = 0;

      len = q21*q21+q31*q31;
      if ( len > 0 ) {
        float tmp = q11-1;
        q12 = -q21;
        q22 = 1+tmp*q21*q21/len;
        q32 = tmp*q21*q31/len;

        q13 = -q31;
        q23 = q32;
        q33 = 1+tmp*q31*q31/len;

        v0 = q12*e+q22*f+q32*g;
        v1 = q12*f+q22*h+q32*i;
        v2 = q12*g+q22*i+q32*j;
        diag[2] = q12*v0+q22*v1+q32*v2;
        subd[2] = q13*v0+q23*v1+q33*v2;

        v0 = q13*e+q23*f+q33*g;
        v1 = q13*f+q23*h+q33*i;
        v2 = q13*g+q23*i+q33*j;
        diag[3] = q13*v0+q23*v1+q33*v2;
      } else { // Q1 = (+-1,0,0)
        q12 = 0; q22 = 1; q32 = 0;
        q13 = 0; q23 = 0; q33 = 1;

        diag[2] = h;
        diag[3] = j;
        subd[2] = i;
      }
    }

    mat[1][1] = q11; mat[1][2] = q12; mat[1][3] = q13;
    mat[2][1] = q21; mat[2][2] = q22; mat[2][3] = q23;
    mat[3][1] = q31; mat[3][2] = q32; mat[3][3] = q33;
  } else {
    diag[1] = e;
    subd[0] = b;
    mat[1][1] = 1;
    mat[2][1] = 0;
    mat[3][1] = 0;

    if ( g != 0 ) {
      float ell = float(sqrt(f*f+g*g));
      f /= ell;
      g /= ell;
      float Q = 2*f*i+g*(j-h);

      diag[2] = h+g*Q;
      diag[3] = j-g*Q;
      subd[1] = ell;
      subd[2] = i-f*Q;
      mat[1][2] = 0;  mat[1][3] = 0;
      mat[2][2] = f;  mat[2][3] = g;
      mat[3][2] = g;  mat[3][3] = -f;
    } else {
      diag[2] = h;
      diag[3] = j;
      subd[1] = f;
      subd[2] = i;
      mat[1][2] = 0;  mat[1][3] = 0;
      mat[2][2] = 1;  mat[2][3] = 0;
      mat[3][2] = 0;  mat[3][3] = 1;
    }
  }
}
//---------------------------------------------------------------------------
void irtkEigenAnalysis::
TridiagonalN (int n, float** mat, float* diag, float* subd)
{
  int i, j, k, ell;

  for (i = n-1, ell = n-2; i >= 1; i--, ell--) {
    float h = 0, scale = 0;

    if ( ell > 0 ) {
      for (k = 0; k <= ell; k++)
        scale += float(fabs(mat[i][k]));
      if ( scale == 0 )
        subd[i] = mat[i][ell];
      else {
        for (k = 0; k <= ell; k++) {
          mat[i][k] /= scale;
          h += mat[i][k]*mat[i][k];
        }
        float f = mat[i][ell];
        float g = ( f > 0 ? -float(sqrt(h)) : float(sqrt(h)) );
        subd[i] = scale*g;
        h -= f*g;
        mat[i][ell] = f-g;
        f = 0;
        for (j = 0; j <= ell; j++) {
          mat[j][i] = mat[i][j]/h;
          g = 0;
          for (k = 0; k <= j; k++)
            g += mat[j][k]*mat[i][k];
          for (k = j+1; k <= ell; k++)
            g += mat[k][j]*mat[i][k];
          subd[j] = g/h;
          f += subd[j]*mat[i][j];
        }
        float hh = f/(h+h);
        for (j = 0; j <= ell; j++) {
          f = mat[i][j];
          subd[j] = g = subd[j] - hh*f;
          for (k = 0; k <= j; k++)
            mat[j][k] -= f*subd[k]+g*mat[i][k];
        }
      }
    } else
      subd[i] = mat[i][ell];

    diag[i] = h;
  }

  diag[0] = subd[0] = 0;
  for (i = 0, ell = -1; i <= n-1; i++, ell++) {
    if ( diag[i] ) {
      for (j = 0; j <= ell; j++) {
        float sum = 0;
        for (k = 0; k <= ell; k++)
          sum += mat[i][k]*mat[k][j];
        for (k = 0; k <= ell; k++)
          mat[k][j] -= sum*mat[k][i];
      }
    }
    diag[i] = mat[i][i];
    mat[i][i] = 1;
    for (j = 0; j <= ell; j++)
      mat[j][i] = mat[i][j] = 0;
  }

  // re-ordering if irtkEigenAnalysis::QLAlgorithm is used subsequently
  for (i = 1, ell = 0; i < n; i++, ell++)
    subd[ell] = subd[i];
  subd[n-1] = 0;
}
//---------------------------------------------------------------------------
void irtkEigenAnalysis::
QLAlgorithm (int n, float* diag, float* subd, float** mat)
{
  const int eigen_maxiter = 30;

  for (int ell = 0; ell < n; ell++) {
    int iter;
    for (iter = 0; iter < eigen_maxiter; iter++) {
      int m;
      for (m = ell; m <= n-2; m++) {
        float dd = float(fabs(diag[m])+fabs(diag[m+1]));
        if ( (float)(fabs(subd[m])+dd) == dd )
          break;
      }
      if ( m == ell )
        break;

      float g = (diag[ell+1]-diag[ell])/(2*subd[ell]);
      float r = float(sqrt(g*g+1));
      if ( g < 0 )
        g = diag[m]-diag[ell]+subd[ell]/(g-r);
      else
        g = diag[m]-diag[ell]+subd[ell]/(g+r);
      float s = 1, c = 1, p = 0;
      for (int i = m-1; i >= ell; i--) {
        float f = s*subd[i], b = c*subd[i];
        if ( fabs(f) >= fabs(g) ) {
          c = g/f;
          r = float(sqrt(c*c+1));
          subd[i+1] = f*r;
          c *= (s = 1/r);
        } else {
          s = f/g;
          r = float(sqrt(s*s+1));
          subd[i+1] = g*r;
          s *= (c = 1/r);
        }
        g = diag[i+1]-p;
        r = (diag[i]-g)*s+2*b*c;
        p = s*r;
        diag[i+1] = g+p;
        g = c*r-b;

        for (int k = 0; k < n; k++) {
          f = mat[k][i+1];
          mat[k][i+1] = s*mat[k][i]+c*f;
          mat[k][i] = c*mat[k][i]-s*f;
        }
      }
      diag[ell] -= p;
      subd[ell] = g;
      subd[m] = 0;
    }
    if ( iter == eigen_maxiter ) {
      Report(ql_exceeded);
      return;
    }
  }
}
//---------------------------------------------------------------------------
void irtkEigenAnalysis::
DecreasingSort (int n, float* eigval, float** eigvec)
{
  // sort eigenvalues in decreasing order, e[0] >= ... >= e[n-1]
  for (int i = 0, k; i <= n-2; i++) {
    // locate maximum eigenvalue
    float max = eigval[k=i];
    int j;
    for (j = i+1; j < n; j++)
      if ( eigval[j] > max )
        max = eigval[k=j];

    if ( k != i ) {
      // swap eigenvalues
      eigval[k] = eigval[i];
      eigval[i] = max;

      // swap eigenvectors
      for (j = 0; j < n; j++) {
        float tmp = eigvec[j][i];
        eigvec[j][i] = eigvec[j][k];
        eigvec[j][k] = tmp;
      }
    }
  }
}
//---------------------------------------------------------------------------
void irtkEigenAnalysis::
IncreasingSort (int n, float* eigval, float** eigvec)
{
  // sort eigenvalues in increasing order, e[0] <= ... <= e[n-1]
  for (int i = 0, k; i <= n-2; i++) {
    // locate minimum eigenvalue
    float min = eigval[k=i];
    int j;
    for (j = i+1; j < n; j++)
      if ( eigval[j] < min )
        min = eigval[k=j];

    if ( k != i ) {
      // swap eigenvalues
      eigval[k] = eigval[i];
      eigval[i] = min;

      // swap eigenvectors
      for (j = 0; j < n; j++) {
        float tmp = eigvec[j][i];
        eigvec[j][i] = eigvec[j][k];
        eigvec[j][k] = tmp;
      }
    }
  }
}
//---------------------------------------------------------------------------
void irtkEigenAnalysis::
IncreasingMagnitudeSort (int n, float* eigval, float** eigvec)
{
  // sort eigenvalues in increasing magnitude order, |e[0]| <= ... <= |e[n-1]|
  for (int i = 0, k; i <= n-2; i++) {
    // locate minimum eigenvalue
    float magmin = eigval[k=i];
    int j;
    for (j = i+1; j < n; j++)
      if ( fabs(eigval[j]) < fabs(magmin) )
        magmin = eigval[k=j];

    if ( k != i ) {
      // swap eigenvalues
      eigval[k] = eigval[i];
      eigval[i] = magmin;

      // swap eigenvectors
      for (j = 0; j < n; j++) {
        float tmp = eigvec[j][i];
        eigvec[j][i] = eigvec[j][k];
        eigvec[j][k] = tmp;
      }
    }
  }
}
//---------------------------------------------------------------------------
void irtkEigenAnalysis::
DecreasingMagnitudeSort (int n, float* eigval, float** eigvec)
{
  // sort eigenvalues in decreasing magnitude order, |e[0]| >= ... >= |e[n-1]|
  for (int i = 0, k; i <= n-2; i++) {
    // locate maximum eigenvalue
    float magmax = eigval[k=i];
    int j;
    for (j = i+1; j < n; j++)
      if ( fabs(eigval[j]) > fabs(magmax) )
        magmax = eigval[k=j];

    if ( k != i ) {
      // swap eigenvalues
      eigval[k] = eigval[i];
      eigval[i] = magmax;

      // swap eigenvectors
      for (j = 0; j < n; j++) {
        float tmp = eigvec[j][i];
        eigvec[j][i] = eigvec[j][k];
        eigvec[j][k] = tmp;
      }
    }
  }
}
//---------------------------------------------------------------------------
irtkEigenAnalysis& irtkEigenAnalysis::
Matrix (float** inmat)
{
  for (int row = 0; row < size; row++)
    for (int col = 0; col < size; col++)
      mat[row][col] = inmat[row][col];
  return *this;
}
//---------------------------------------------------------------------------
void irtkEigenAnalysis::
EigenStuff2 ()
{
  Tridiagonal2(mat,diag,subd);
  QLAlgorithm(size,diag,subd,mat);
}
//---------------------------------------------------------------------------
void irtkEigenAnalysis::
EigenStuff3 ()
{
  Tridiagonal3(mat,diag,subd);
  QLAlgorithm(size,diag,subd,mat);
}
//---------------------------------------------------------------------------
void irtkEigenAnalysis::
EigenStuff4 ()
{
  Tridiagonal4(mat,diag,subd);
  QLAlgorithm(size,diag,subd,mat);
}
//---------------------------------------------------------------------------
void irtkEigenAnalysis::
EigenStuffN ()
{
  TridiagonalN(size,mat,diag,subd);
  QLAlgorithm(size,diag,subd,mat);
}
//---------------------------------------------------------------------------
void irtkEigenAnalysis::
EigenStuff ()
{
  switch ( size ) {
  case 2 : Tridiagonal2(mat,diag,subd);       break;
  case 3 : Tridiagonal3(mat,diag,subd);       break;
  case 4 : Tridiagonal4(mat,diag,subd);       break;
  default: TridiagonalN(size,mat,diag,subd);  break;
  }
  QLAlgorithm(size,diag,subd,mat);
}
//---------------------------------------------------------------------------
void irtkEigenAnalysis::
DecrSortEigenStuff2 ()
{
  Tridiagonal2(mat,diag,subd);
  QLAlgorithm(size,diag,subd,mat);
  DecreasingSort(size,diag,mat);
}
//---------------------------------------------------------------------------
void irtkEigenAnalysis::
DecrSortEigenStuff3 ()
{
  Tridiagonal3(mat,diag,subd);
  QLAlgorithm(size,diag,subd,mat);
  DecreasingSort(size,diag,mat);
}
//---------------------------------------------------------------------------
void irtkEigenAnalysis::
DecrSortEigenStuff4 ()
{
  Tridiagonal4(mat,diag,subd);
  QLAlgorithm(size,diag,subd,mat);
  DecreasingSort(size,diag,mat);
}
//---------------------------------------------------------------------------
void irtkEigenAnalysis::
DecrSortEigenStuffN ()
{
  TridiagonalN(size,mat,diag,subd);
  QLAlgorithm(size,diag,subd,mat);
  DecreasingSort(size,diag,mat);
}
//---------------------------------------------------------------------------
void irtkEigenAnalysis::
DecrSortEigenStuff ()
{
  switch ( size ) {
  case 2 : Tridiagonal2(mat,diag,subd);       break;
  case 3 : Tridiagonal3(mat,diag,subd);       break;
  case 4 : Tridiagonal4(mat,diag,subd);       break;
  default: TridiagonalN(size,mat,diag,subd);  break;
  }
  QLAlgorithm(size,diag,subd,mat);
  DecreasingSort(size,diag,mat);
}
//---------------------------------------------------------------------------
void irtkEigenAnalysis::
IncrSortEigenStuff2 ()
{
  Tridiagonal2(mat,diag,subd);
  QLAlgorithm(size,diag,subd,mat);
  IncreasingSort(size,diag,mat);
}
//---------------------------------------------------------------------------
void irtkEigenAnalysis::
IncrSortEigenStuff3 ()
{
  Tridiagonal3(mat,diag,subd);
  QLAlgorithm(size,diag,subd,mat);
  IncreasingSort(size,diag,mat);
}
//---------------------------------------------------------------------------
void irtkEigenAnalysis::
IncrSortEigenStuff4 ()
{
  Tridiagonal4(mat,diag,subd);
  QLAlgorithm(size,diag,subd,mat);
  IncreasingSort(size,diag,mat);
}
//---------------------------------------------------------------------------
void irtkEigenAnalysis::
IncrSortEigenStuffN ()
{
  TridiagonalN(size,mat,diag,subd);
  QLAlgorithm(size,diag,subd,mat);
  IncreasingSort(size,diag,mat);
}
//---------------------------------------------------------------------------
void irtkEigenAnalysis::
IncrSortEigenStuff ()
{
  switch ( size ) {
  case 2 : Tridiagonal2(mat,diag,subd);       break;
  case 3 : Tridiagonal3(mat,diag,subd);       break;
  case 4 : Tridiagonal4(mat,diag,subd);       break;
  default: TridiagonalN(size,mat,diag,subd);  break;
  }
  QLAlgorithm(size,diag,subd,mat);
  IncreasingSort(size,diag,mat);
}
//---------------------------------------------------------------------------
void irtkEigenAnalysis::
DecrMagSortEigenStuff2 ()
{
  Tridiagonal2(mat,diag,subd);
  QLAlgorithm(size,diag,subd,mat);
  DecreasingMagnitudeSort(size,diag,mat);
}
//---------------------------------------------------------------------------
void irtkEigenAnalysis::
DecrMagSortEigenStuff3 ()
{
  Tridiagonal3(mat,diag,subd);
  QLAlgorithm(size,diag,subd,mat);
  DecreasingMagnitudeSort(size,diag,mat);
}
//---------------------------------------------------------------------------
void irtkEigenAnalysis::
DecrMagSortEigenStuff4 ()
{
  Tridiagonal4(mat,diag,subd);
  QLAlgorithm(size,diag,subd,mat);
  DecreasingMagnitudeSort(size,diag,mat);
}
//---------------------------------------------------------------------------
void irtkEigenAnalysis::
DecrMagSortEigenStuffN ()
{
  TridiagonalN(size,mat,diag,subd);
  QLAlgorithm(size,diag,subd,mat);
  DecreasingMagnitudeSort(size,diag,mat);
}
//---------------------------------------------------------------------------
void irtkEigenAnalysis::
DecrMagSortEigenStuff ()
{
  switch ( size ) {
  case 2 : Tridiagonal2(mat,diag,subd);       break;
  case 3 : Tridiagonal3(mat,diag,subd);       break;
  case 4 : Tridiagonal4(mat,diag,subd);       break;
  default: TridiagonalN(size,mat,diag,subd);  break;
  }
  QLAlgorithm(size,diag,subd,mat);
  DecreasingMagnitudeSort(size,diag,mat);
}
//---------------------------------------------------------------------------
void irtkEigenAnalysis::
IncrMagSortEigenStuff2 ()
{
  Tridiagonal2(mat,diag,subd);
  QLAlgorithm(size,diag,subd,mat);
  IncreasingMagnitudeSort(size,diag,mat);
}
//---------------------------------------------------------------------------
void irtkEigenAnalysis::
IncrMagSortEigenStuff3 ()
{
  Tridiagonal3(mat,diag,subd);
  QLAlgorithm(size,diag,subd,mat);
  IncreasingMagnitudeSort(size,diag,mat);
}
//---------------------------------------------------------------------------
void irtkEigenAnalysis::
IncrMagSortEigenStuff4 ()
{
  Tridiagonal4(mat,diag,subd);
  QLAlgorithm(size,diag,subd,mat);
  IncreasingMagnitudeSort(size,diag,mat);
}
//---------------------------------------------------------------------------
void irtkEigenAnalysis::
IncrMagSortEigenStuffN ()
{
  TridiagonalN(size,mat,diag,subd);
  QLAlgorithm(size,diag,subd,mat);
  IncreasingMagnitudeSort(size,diag,mat);
}
//---------------------------------------------------------------------------
void irtkEigenAnalysis::
IncrMagSortEigenStuff ()
{
  switch ( size ) {
  case 2 : Tridiagonal2(mat,diag,subd);       break;
  case 3 : Tridiagonal3(mat,diag,subd);       break;
  case 4 : Tridiagonal4(mat,diag,subd);       break;
  default: TridiagonalN(size,mat,diag,subd);  break;
  }
  QLAlgorithm(size,diag,subd,mat);
  IncreasingMagnitudeSort(size,diag,mat);
}
//---------------------------------------------------------------------------
int irtkEigenAnalysis::
Number (unsigned single_error)
{
  int result;
  for (result = -1; single_error; single_error >>= 1)
    result++;
  return result;
}
//---------------------------------------------------------------------------
void irtkEigenAnalysis::
Report (unsigned single_error)
{
  if ( irtkEigenAnalysis::verbose )
    cout << "irtkEigenAnalysis: " << message[Number(single_error)] << endl;
  else
    ofstream("eigen.err",ios::out|ios::app)
    << "irtkEigenAnalysis: " << message[Number(single_error)] << endl;

  error |= single_error;
}
//---------------------------------------------------------------------------
void irtkEigenAnalysis::
Report (ostream& ostr)
{
  for (unsigned single_error = 1; single_error; single_error <<= 1)
    if ( error & single_error )
      ostr << "irtkEigenAnalysis: " << message[Number(single_error)] << endl;

  error = 0;
}



//===========================================================================
// error handling
int irtkEigenAnalysisD::verbose = 0;
unsigned irtkEigenAnalysisD::error = 0;
const unsigned irtkEigenAnalysisD::invalid_size      = 0x00000001;
const unsigned irtkEigenAnalysisD::allocation_failed = 0x00000002;
const unsigned irtkEigenAnalysisD::ql_exceeded       = 0x00000004;
const char* irtkEigenAnalysisD::message[3] = {
  "invalid matrix size",
  "allocation failed",
  "QL algorithm - exceeded maximum iterations"
};
irtkEigenAnalysisD::
irtkEigenAnalysisD (int _size)
{
  if ( (size = _size) <= 1 ) {
    Report(invalid_size);
    return;
  }
  if ( (mat = new double*[size]) == 0 ) {
    Report(allocation_failed);
    return;
  }
  for (int d = 0; d < size; d++)
    if ( (mat[d] = new double[size]) == 0 ) {
      Report(allocation_failed);
      return;
    }
  if ( (diag = new double[size]) == 0 ) {
    Report(allocation_failed);
    return;
  }
  if ( (subd = new double[size]) == 0 ) {
    Report(allocation_failed);
    return;
  }
}
//---------------------------------------------------------------------------
irtkEigenAnalysisD::
~irtkEigenAnalysisD ()
{
  delete[] subd;
  delete[] diag;
  for (int d = 0; d < size; d++)
    delete[] mat[d];
  delete[] mat;
}

//---------------------------------------------------------------------------
irtkEigenAnalysisD& irtkEigenAnalysisD::
Matrix (double** inmat)
{
  for (int row = 0; row < size; row++)
    for (int col = 0; col < size; col++)
      mat[row][col] = inmat[row][col];
  return *this;
}


//---------------------------------------------------------------------------
void irtkEigenAnalysisD::
QLAlgorithm (int n, double* diag, double* subd, double** mat)
{
  const int eigen_maxiter = 30;

  for (int ell = 0; ell < n; ell++) {
    int iter;
    for (iter = 0; iter < eigen_maxiter; iter++) {
      int m;
      for (m = ell; m <= n-2; m++) {
        double dd = double(fabs(diag[m])+fabs(diag[m+1]));
        if ( (double)(fabs(subd[m])+dd) == dd )
          break;
      }
      if ( m == ell )
        break;

      double g = (diag[ell+1]-diag[ell])/(2*subd[ell]);
      double r = double(sqrt(g*g+1));
      if ( g < 0 )
        g = diag[m]-diag[ell]+subd[ell]/(g-r);
      else
        g = diag[m]-diag[ell]+subd[ell]/(g+r);
      double s = 1, c = 1, p = 0;
      for (int i = m-1; i >= ell; i--) {
        double f = s*subd[i], b = c*subd[i];
        if ( fabs(f) >= fabs(g) ) {
          c = g/f;
          r = double(sqrt(c*c+1));
          subd[i+1] = f*r;
          c *= (s = 1/r);
        } else {
          s = f/g;
          r = double(sqrt(s*s+1));
          subd[i+1] = g*r;
          s *= (c = 1/r);
        }
        g = diag[i+1]-p;
        r = (diag[i]-g)*s+2*b*c;
        p = s*r;
        diag[i+1] = g+p;
        g = c*r-b;

        for (int k = 0; k < n; k++) {
          f = mat[k][i+1];
          mat[k][i+1] = s*mat[k][i]+c*f;
          mat[k][i] = c*mat[k][i]-s*f;
        }
      }
      diag[ell] -= p;
      subd[ell] = g;
      subd[m] = 0;
    }
    if ( iter == eigen_maxiter ) {
      Report(ql_exceeded);
      return;
    }
  }
}
//---------------------------------------------------------------------------
void irtkEigenAnalysisD::
EigenStuff2 ()
{
  Tridiagonal2(mat,diag,subd);
  QLAlgorithm(size,diag,subd,mat);
}
//---------------------------------------------------------------------------
void irtkEigenAnalysisD::
EigenStuff3 ()
{
  Tridiagonal3(mat,diag,subd);
  QLAlgorithm(size,diag,subd,mat);
}
//---------------------------------------------------------------------------
void irtkEigenAnalysisD::
EigenStuff4 ()
{
  Tridiagonal4(mat,diag,subd);
  QLAlgorithm(size,diag,subd,mat);
}
//---------------------------------------------------------------------------
void irtkEigenAnalysisD::
EigenStuffN ()
{
  TridiagonalN(size,mat,diag,subd);
  QLAlgorithm(size,diag,subd,mat);
}
//---------------------------------------------------------------------------
void irtkEigenAnalysisD::
EigenStuff ()
{
  switch ( size ) {
  case 2 : Tridiagonal2(mat,diag,subd);       break;
  case 3 : Tridiagonal3(mat,diag,subd);       break;
  case 4 : Tridiagonal4(mat,diag,subd);       break;
  default: TridiagonalN(size,mat,diag,subd);  break;
  }
  QLAlgorithm(size,diag,subd,mat);
}
//---------------------------------------------------------------------------
void irtkEigenAnalysisD::
DecrSortEigenStuff2 ()
{
  Tridiagonal2(mat,diag,subd);
  QLAlgorithm(size,diag,subd,mat);
  DecreasingSort(size,diag,mat);
}
//---------------------------------------------------------------------------
void irtkEigenAnalysisD::
DecrSortEigenStuff3 ()
{
  Tridiagonal3(mat,diag,subd);
  QLAlgorithm(size,diag,subd,mat);
  DecreasingSort(size,diag,mat);
}
//---------------------------------------------------------------------------
void irtkEigenAnalysisD::
DecrSortEigenStuff4 ()
{
  Tridiagonal4(mat,diag,subd);
  QLAlgorithm(size,diag,subd,mat);
  DecreasingSort(size,diag,mat);
}
//---------------------------------------------------------------------------
void irtkEigenAnalysisD::
DecrSortEigenStuffN ()
{
  TridiagonalN(size,mat,diag,subd);
  QLAlgorithm(size,diag,subd,mat);
  DecreasingSort(size,diag,mat);
}
//---------------------------------------------------------------------------
void irtkEigenAnalysisD::
DecrSortEigenStuff ()
{
  switch ( size ) {
  case 2 : Tridiagonal2(mat,diag,subd);       break;
  case 3 : Tridiagonal3(mat,diag,subd);       break;
  case 4 : Tridiagonal4(mat,diag,subd);       break;
  default: TridiagonalN(size,mat,diag,subd);  break;
  }
  QLAlgorithm(size,diag,subd,mat);
  DecreasingSort(size,diag,mat);
}
//---------------------------------------------------------------------------
void irtkEigenAnalysisD::
IncrSortEigenStuff2 ()
{
  Tridiagonal2(mat,diag,subd);
  QLAlgorithm(size,diag,subd,mat);
  IncreasingSort(size,diag,mat);
}
//---------------------------------------------------------------------------
void irtkEigenAnalysisD::
IncrSortEigenStuff3 ()
{
  Tridiagonal3(mat,diag,subd);
  QLAlgorithm(size,diag,subd,mat);
  IncreasingSort(size,diag,mat);
}
//---------------------------------------------------------------------------
void irtkEigenAnalysisD::
IncrSortEigenStuff4 ()
{
  Tridiagonal4(mat,diag,subd);
  QLAlgorithm(size,diag,subd,mat);
  IncreasingSort(size,diag,mat);
}
//---------------------------------------------------------------------------
void irtkEigenAnalysisD::
IncrSortEigenStuffN ()
{
  TridiagonalN(size,mat,diag,subd);
  QLAlgorithm(size,diag,subd,mat);
  IncreasingSort(size,diag,mat);
}
//---------------------------------------------------------------------------
void irtkEigenAnalysisD::
IncrSortEigenStuff ()
{
  switch ( size ) {
  case 2 : Tridiagonal2(mat,diag,subd);       break;
  case 3 : Tridiagonal3(mat,diag,subd);       break;
  case 4 : Tridiagonal4(mat,diag,subd);       break;
  default: TridiagonalN(size,mat,diag,subd);  break;
  }
  QLAlgorithm(size,diag,subd,mat);
  IncreasingSort(size,diag,mat);
}

//---------------------------------------------------------------------------
void irtkEigenAnalysisD::
DecrMagSortEigenStuff2 ()
{
  Tridiagonal2(mat,diag,subd);
  QLAlgorithm(size,diag,subd,mat);
  DecreasingMagnitudeSort(size,diag,mat);
}
//---------------------------------------------------------------------------
void irtkEigenAnalysisD::
DecrMagSortEigenStuff3 ()
{
  Tridiagonal3(mat,diag,subd);
  QLAlgorithm(size,diag,subd,mat);
  DecreasingMagnitudeSort(size,diag,mat);
}
//---------------------------------------------------------------------------
void irtkEigenAnalysisD::
DecrMagSortEigenStuff4 ()
{
  Tridiagonal4(mat,diag,subd);
  QLAlgorithm(size,diag,subd,mat);
  DecreasingMagnitudeSort(size,diag,mat);
}
//---------------------------------------------------------------------------
void irtkEigenAnalysisD::
DecrMagSortEigenStuffN ()
{
  TridiagonalN(size,mat,diag,subd);
  QLAlgorithm(size,diag,subd,mat);
  DecreasingMagnitudeSort(size,diag,mat);
}
//---------------------------------------------------------------------------
void irtkEigenAnalysisD::
DecrMagSortEigenStuff ()
{
  switch ( size ) {
  case 2 : Tridiagonal2(mat,diag,subd);       break;
  case 3 : Tridiagonal3(mat,diag,subd);       break;
  case 4 : Tridiagonal4(mat,diag,subd);       break;
  default: TridiagonalN(size,mat,diag,subd);  break;
  }
  QLAlgorithm(size,diag,subd,mat);
  DecreasingMagnitudeSort(size,diag,mat);
}
//---------------------------------------------------------------------------
void irtkEigenAnalysisD::
IncrMagSortEigenStuff2 ()
{
  Tridiagonal2(mat,diag,subd);
  QLAlgorithm(size,diag,subd,mat);
  IncreasingMagnitudeSort(size,diag,mat);
}
//---------------------------------------------------------------------------
void irtkEigenAnalysisD::
IncrMagSortEigenStuff3 ()
{
  Tridiagonal3(mat,diag,subd);
  QLAlgorithm(size,diag,subd,mat);
  IncreasingMagnitudeSort(size,diag,mat);
}
//---------------------------------------------------------------------------
void irtkEigenAnalysisD::
IncrMagSortEigenStuff4 ()
{
  Tridiagonal4(mat,diag,subd);
  QLAlgorithm(size,diag,subd,mat);
  IncreasingMagnitudeSort(size,diag,mat);
}
//---------------------------------------------------------------------------
void irtkEigenAnalysisD::
IncrMagSortEigenStuffN ()
{
  TridiagonalN(size,mat,diag,subd);
  QLAlgorithm(size,diag,subd,mat);
  IncreasingMagnitudeSort(size,diag,mat);
}
//---------------------------------------------------------------------------
void irtkEigenAnalysisD::
IncrMagSortEigenStuff ()
{
  switch ( size ) {
  case 2 : Tridiagonal2(mat,diag,subd);       break;
  case 3 : Tridiagonal3(mat,diag,subd);       break;
  case 4 : Tridiagonal4(mat,diag,subd);       break;
  default: TridiagonalN(size,mat,diag,subd);  break;
  }
  QLAlgorithm(size,diag,subd,mat);
  IncreasingMagnitudeSort(size,diag,mat);
}


//---------------------------------------------------------------------------
void irtkEigenAnalysisD::
DecreasingSort (int n, double* eigval, double** eigvec)
{
  // sort eigenvalues in decreasing order, e[0] >= ... >= e[n-1]
  for (int i = 0, k; i <= n-2; i++) {
    // locate maximum eigenvalue
    double max = eigval[k=i];
    int j;
    for (j = i+1; j < n; j++)
      if ( eigval[j] > max )
        max = eigval[k=j];

    if ( k != i ) {
      // swap eigenvalues
      eigval[k] = eigval[i];
      eigval[i] = max;

      // swap eigenvectors
      for (j = 0; j < n; j++) {
        double tmp = eigvec[j][i];
        eigvec[j][i] = eigvec[j][k];
        eigvec[j][k] = tmp;
      }
    }
  }
}
//---------------------------------------------------------------------------
void irtkEigenAnalysisD::
IncreasingSort (int n, double* eigval, double** eigvec)
{
  // sort eigenvalues in increasing order, e[0] <= ... <= e[n-1]
  for (int i = 0, k; i <= n-2; i++) {
    // locate minimum eigenvalue
    double min = eigval[k=i];
    int j;
    for (j = i+1; j < n; j++)
      if ( eigval[j] < min )
        min = eigval[k=j];

    if ( k != i ) {
      // swap eigenvalues
      eigval[k] = eigval[i];
      eigval[i] = min;

      // swap eigenvectors
      for (j = 0; j < n; j++) {
        double tmp = eigvec[j][i];
        eigvec[j][i] = eigvec[j][k];
        eigvec[j][k] = tmp;
      }
    }
  }
}
//---------------------------------------------------------------------------
void irtkEigenAnalysisD::
IncreasingMagnitudeSort (int n, double* eigval, double** eigvec)
{
  // sort eigenvalues in increasing magnitude order, |e[0]| <= ... <= |e[n-1]|
  for (int i = 0, k; i <= n-2; i++) {
    // locate minimum eigenvalue
    double magmin = eigval[k=i];
    int j;
    for (j = i+1; j < n; j++)
      if ( fabs(eigval[j]) < fabs(magmin) )
        magmin = eigval[k=j];

    if ( k != i ) {
      // swap eigenvalues
      eigval[k] = eigval[i];
      eigval[i] = magmin;

      // swap eigenvectors
      for (j = 0; j < n; j++) {
        double tmp = eigvec[j][i];
        eigvec[j][i] = eigvec[j][k];
        eigvec[j][k] = tmp;
      }
    }
  }
}
//---------------------------------------------------------------------------
void irtkEigenAnalysisD::
DecreasingMagnitudeSort (int n, double* eigval, double** eigvec)
{
  // sort eigenvalues in decreasing magnitude order, |e[0]| >= ... >= |e[n-1]|
  for (int i = 0, k; i <= n-2; i++) {
    // locate maximum eigenvalue
    double magmax = eigval[k=i];
    int j;
    for (j = i+1; j < n; j++)
      if ( fabs(eigval[j]) > fabs(magmax) )
        magmax = eigval[k=j];

    if ( k != i ) {
      // swap eigenvalues
      eigval[k] = eigval[i];
      eigval[i] = magmax;

      // swap eigenvectors
      for (j = 0; j < n; j++) {
        double tmp = eigvec[j][i];
        eigvec[j][i] = eigvec[j][k];
        eigvec[j][k] = tmp;
      }
    }
  }
}
//---------------------------------------------------------------------------
int irtkEigenAnalysisD::
Number (unsigned single_error)
{
  int result;
  for (result = -1; single_error; single_error >>= 1)
    result++;
  return result;
}

//---------------------------------------------------------------------------
void irtkEigenAnalysisD::
Tridiagonal2 (double** mat, double* diag, double* subd)
{
  // matrix is already tridiagonal

  diag[0] = mat[0][0];
  diag[1] = mat[1][1];
  subd[0] = mat[0][1];
  subd[1] = 0;
  mat[0][0] = 1;  mat[0][1] = 0;
  mat[1][0] = 0;  mat[1][1] = 1;
}
//---------------------------------------------------------------------------
void irtkEigenAnalysisD::
Tridiagonal3 (double** mat, double* diag, double* subd)
{
  double a = mat[0][0], b = mat[0][1], c = mat[0][2],
                            d = mat[1][1], e = mat[1][2],
                                               f = mat[2][2];

  diag[0] = a;
  subd[2] = 0;
  if ( c != 0 ) {
    double ell = double(sqrt(b*b+c*c));
    b /= ell;
    c /= ell;
    double q = 2*b*e+c*(f-d);
    diag[1] = d+c*q;
    diag[2] = f-c*q;
    subd[0] = ell;
    subd[1] = e-b*q;
    mat[0][0] = 1; mat[0][1] = 0; mat[0][2] = 0;
    mat[1][0] = 0; mat[1][1] = b; mat[1][2] = c;
    mat[2][0] = 0; mat[2][1] = c; mat[2][2] = -b;
  } else {
    diag[1] = d;
    diag[2] = f;
    subd[0] = b;
    subd[1] = e;
    mat[0][0] = 1; mat[0][1] = 0; mat[0][2] = 0;
    mat[1][0] = 0; mat[1][1] = 1; mat[1][2] = 0;
    mat[2][0] = 0; mat[2][1] = 0; mat[2][2] = 1;
  }
}
//---------------------------------------------------------------------------
void irtkEigenAnalysisD::
Tridiagonal4 (double** mat, double* diag, double* subd)
{
  // save matrix M
  double
  a = mat[0][0], b = mat[0][1], c = mat[0][2], d = mat[0][3],
                                    e = mat[1][1], f = mat[1][2], g = mat[1][3],
                                                                      h = mat[2][2], i = mat[2][3],
                                                                                         j = mat[3][3];

  diag[0] = a;
  subd[3] = 0;

  mat[0][0] = 1; mat[0][1] = 0; mat[0][2] = 0; mat[0][3] = 0;
  mat[1][0] = 0;
  mat[2][0] = 0;
  mat[3][0] = 0;

  if ( c != 0 || d != 0 ) {
    double q11, q12, q13;
    double q21, q22, q23;
    double q31, q32, q33;

    // build column Q1
    double len = double(sqrt(b*b+c*c+d*d));
    q11 = b/len;
    q21 = c/len;
    q31 = d/len;

    subd[0] = len;

    // compute S*Q1
    double v0 = e*q11+f*q21+g*q31;
    double v1 = f*q11+h*q21+i*q31;
    double v2 = g*q11+i*q21+j*q31;

    diag[1] = q11*v0+q21*v1+q31*v2;

    // build column Q3 = Q1x(S*Q1)
    q13 = q21*v2-q31*v1;
    q23 = q31*v0-q11*v2;
    q33 = q11*v1-q21*v0;
    len = double(sqrt(q13*q13+q23*q23+q33*q33));
    if ( len > 0 ) {
      q13 /= len;
      q23 /= len;
      q33 /= len;

      // build column Q2 = Q3xQ1
      q12 = q23*q31-q33*q21;
      q22 = q33*q11-q13*q31;
      q32 = q13*q21-q23*q11;

      v0 = q12*e+q22*f+q32*g;
      v1 = q12*f+q22*h+q32*i;
      v2 = q12*g+q22*i+q32*j;
      subd[1] = q11*v0+q21*v1+q31*v2;
      diag[2] = q12*v0+q22*v1+q32*v2;
      subd[2] = q13*v0+q23*v1+q33*v2;

      v0 = q13*e+q23*f+q33*g;
      v1 = q13*f+q23*h+q33*i;
      v2 = q13*g+q23*i+q33*j;
      diag[3] = q13*v0+q23*v1+q33*v2;
    } else { // S*Q1 parallel to Q1, choose any valid Q2 and Q3
      subd[1] = 0;

      len = q21*q21+q31*q31;
      if ( len > 0 ) {
        double tmp = q11-1;
        q12 = -q21;
        q22 = 1+tmp*q21*q21/len;
        q32 = tmp*q21*q31/len;

        q13 = -q31;
        q23 = q32;
        q33 = 1+tmp*q31*q31/len;

        v0 = q12*e+q22*f+q32*g;
        v1 = q12*f+q22*h+q32*i;
        v2 = q12*g+q22*i+q32*j;
        diag[2] = q12*v0+q22*v1+q32*v2;
        subd[2] = q13*v0+q23*v1+q33*v2;

        v0 = q13*e+q23*f+q33*g;
        v1 = q13*f+q23*h+q33*i;
        v2 = q13*g+q23*i+q33*j;
        diag[3] = q13*v0+q23*v1+q33*v2;
      } else { // Q1 = (+-1,0,0)
        q12 = 0; q22 = 1; q32 = 0;
        q13 = 0; q23 = 0; q33 = 1;

        diag[2] = h;
        diag[3] = j;
        subd[2] = i;
      }
    }

    mat[1][1] = q11; mat[1][2] = q12; mat[1][3] = q13;
    mat[2][1] = q21; mat[2][2] = q22; mat[2][3] = q23;
    mat[3][1] = q31; mat[3][2] = q32; mat[3][3] = q33;
  } else {
    diag[1] = e;
    subd[0] = b;
    mat[1][1] = 1;
    mat[2][1] = 0;
    mat[3][1] = 0;

    if ( g != 0 ) {
      double ell = double(sqrt(f*f+g*g));
      f /= ell;
      g /= ell;
      double Q = 2*f*i+g*(j-h);

      diag[2] = h+g*Q;
      diag[3] = j-g*Q;
      subd[1] = ell;
      subd[2] = i-f*Q;
      mat[1][2] = 0;  mat[1][3] = 0;
      mat[2][2] = f;  mat[2][3] = g;
      mat[3][2] = g;  mat[3][3] = -f;
    } else {
      diag[2] = h;
      diag[3] = j;
      subd[1] = f;
      subd[2] = i;
      mat[1][2] = 0;  mat[1][3] = 0;
      mat[2][2] = 1;  mat[2][3] = 0;
      mat[3][2] = 0;  mat[3][3] = 1;
    }
  }
}
//---------------------------------------------------------------------------
void irtkEigenAnalysisD::
TridiagonalN (int n, double** mat, double* diag, double* subd)
{
  int i, j, k, ell;

  for (i = n-1, ell = n-2; i >= 1; i--, ell--) {
    double h = 0, scale = 0;

    if ( ell > 0 ) {
      for (k = 0; k <= ell; k++)
        scale += double(fabs(mat[i][k]));
      if ( scale == 0 )
        subd[i] = mat[i][ell];
      else {
        for (k = 0; k <= ell; k++) {
          mat[i][k] /= scale;
          h += mat[i][k]*mat[i][k];
        }
        double f = mat[i][ell];
        double g = ( f > 0 ? -double(sqrt(h)) : double(sqrt(h)) );
        subd[i] = scale*g;
        h -= f*g;
        mat[i][ell] = f-g;
        f = 0;
        for (j = 0; j <= ell; j++) {
          mat[j][i] = mat[i][j]/h;
          g = 0;
          for (k = 0; k <= j; k++)
            g += mat[j][k]*mat[i][k];
          for (k = j+1; k <= ell; k++)
            g += mat[k][j]*mat[i][k];
          subd[j] = g/h;
          f += subd[j]*mat[i][j];
        }
        double hh = f/(h+h);
        for (j = 0; j <= ell; j++) {
          f = mat[i][j];
          subd[j] = g = subd[j] - hh*f;
          for (k = 0; k <= j; k++)
            mat[j][k] -= f*subd[k]+g*mat[i][k];
        }
      }
    } else
      subd[i] = mat[i][ell];

    diag[i] = h;
  }

  diag[0] = subd[0] = 0;
  for (i = 0, ell = -1; i <= n-1; i++, ell++) {
    if ( diag[i] ) {
      for (j = 0; j <= ell; j++) {
        double sum = 0;
        for (k = 0; k <= ell; k++)
          sum += mat[i][k]*mat[k][j];
        for (k = 0; k <= ell; k++)
          mat[k][j] -= sum*mat[k][i];
      }
    }
    diag[i] = mat[i][i];
    mat[i][i] = 1;
    for (j = 0; j <= ell; j++)
      mat[j][i] = mat[i][j] = 0;
  }

  // re-ordering if irtkEigenAnalysisD::QLAlgorithm is used subsequently
  for (i = 1, ell = 0; i < n; i++, ell++)
    subd[ell] = subd[i];
  subd[n-1] = 0;
}

//---------------------------------------------------------------------------
void irtkEigenAnalysisD::
Report (unsigned single_error)
{
  if ( irtkEigenAnalysisD::verbose )
    cout << "irtkEigenAnalysisD: " << message[Number(single_error)] << endl;
  else
    ofstream("eigen.err",ios::out|ios::app)
    << "irtkEigenAnalysisD: " << message[Number(single_error)] << endl;

  error |= single_error;
}
//---------------------------------------------------------------------------
void irtkEigenAnalysisD::
Report (ostream& ostr)
{
  for (unsigned single_error = 1; single_error; single_error <<= 1)
    if ( error & single_error )
      ostr << "irtkEigenAnalysisD: " << message[Number(single_error)] << endl;

  error = 0;
}
//===========================================================================

#ifdef EIGEN_TEST

int main ()
{
  irtkEigenAnalysisD eig(3);

  eig.Matrix(0,0) = 2;  eig.Matrix(0,1) = 1;  eig.Matrix(0,2) = 1;
  eig.Matrix(1,0) = 1;  eig.Matrix(1,1) = 2;  eig.Matrix(1,2) = 1;
  eig.Matrix(2,0) = 1;  eig.Matrix(2,1) = 1;  eig.Matrix(2,2) = 2;

  eig.IncrSortEigenStuff3();

  cout.setf(ios::fixed);

  cout << "eigenvalues = " << endl;
  for (int row = 0; row < 3; row++)
    cout << eig.Eigenvalue(row) << ' ';
  cout << endl;

  cout << "eigenvectors = " << endl;
  for (row = 0; row < 3; row++) {
    for (int col = 0; col < 3; col++)
      cout << eig.Eigenvector(row,col) << ' ';
    cout << endl;
  }

  // eigenvalues =
  //    1.000000 1.000000 4.000000
  // eigenvectors =
  //    0.411953  0.704955 0.577350
  //    0.404533 -0.709239 0.577350
  //   -0.816485  0.004284 0.577350

  return 0;
}

#endif
