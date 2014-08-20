/*=========================================================================

  Library   : Image Registration Toolkit (IRTK)
  Module    : $Id$
  Copyright : Imperial College, Department of Computing
              Visual Information Processing (VIP), 2008 onwards
  Date      : $Date$
  Version   : $Revision$
  Changes   : $Author$

Copyright (c) IXICO LIMITED
All rights reserved.
See COPYRIGHT for details

=========================================================================*/

#include <irtkImage.h>

#include <irtkTransformation.h>

#include <gsl/gsl_vector.h>
#include <gsl/gsl_multiroots.h>

#define MAXITS 100
#define TOLF 1.0e-4
#define ALF 1.0e-4
#define TOLMIN 1.0e-6
#define TOLX 1.0e-7
#define STPMX 100.0

irtkTransformation *irtkTransformationPointer;

double x_invert, y_invert, z_invert;


void irtkTransformationEvaluate(int, float point[], float fval[])
{
  double xpoint,ypoint,zpoint;

  xpoint = point[1];
  ypoint = point[2];
  zpoint = point[3];

  irtkTransformationPointer->Transform(xpoint, ypoint, zpoint);

  fval[1] = xpoint-x_invert;
  fval[2] = ypoint-y_invert;
  fval[3] = zpoint-z_invert;
}

int evaluate_f (const gsl_vector * x, void * p, gsl_vector * f) {
  float point[3], func[3];
  int i;
  
  for (i = 0; i < 3; i++) {
     point[i] = gsl_vector_get(x, i);
  }

  irtkTransformationEvaluate(3, point-1, func-1);

  for (i = 0; i < 3; i++) {
     gsl_vector_set(f, i, func[i]);
  }

  return GSL_SUCCESS;
}

int evaluate_df (const gsl_vector * x, void * p, gsl_matrix * J) {
  irtkMatrix jacobian(3, 3);

  const double xpoint = gsl_vector_get(x, 0);
  const double ypoint = gsl_vector_get(x, 1);
  const double zpoint = gsl_vector_get(x, 2);

  irtkTransformationPointer->Jacobian(jacobian,xpoint,ypoint,zpoint);

  int i, j;
  for (i = 0; i < 3; i++) {
     for (j = 0; j < 3; j++) {
        gsl_matrix_set (J, i, j, jacobian(i, j));
     }
  }
  return GSL_SUCCESS;
}

int evaluate_fdf (const gsl_vector * x, void * params,
                  gsl_vector *f, gsl_matrix * J) {
  evaluate_f(x, params, f);
  evaluate_df(x, params, J);
  return GSL_SUCCESS;
}

void newt2(float x[], int n, int *check, void (*vecfunc)(int, float [], float [])) {
  const gsl_multiroot_fdfsolver_type * T;
  gsl_multiroot_fdfsolver * s;

  int status, i;
  size_t iter = 0;

  gsl_multiroot_function_fdf f = {&evaluate_f,
                                  &evaluate_df,
                                  &evaluate_fdf,
                                  n, NULL};

  gsl_vector *xvec = gsl_vector_alloc(n);
  for (i = 0; i < n; i++) {
     gsl_vector_set(xvec, i, x[i+1]);
  }
  
  T = gsl_multiroot_fdfsolver_newton;
  s = gsl_multiroot_fdfsolver_alloc(T, n);
  gsl_multiroot_fdfsolver_set(s, &f, xvec);

  do {
     iter++;

     status = gsl_multiroot_fdfsolver_iterate(s);

     if (status)
        break;

     status = gsl_multiroot_test_residual(s->f, TOLF);
  } while (status == GSL_CONTINUE && iter < MAXITS);

  for (i = 0; i < n; i++) {
     x[i+1] = gsl_vector_get(s->x, i);
  } 

  gsl_multiroot_fdfsolver_free(s);
  gsl_vector_free(xvec);
}
