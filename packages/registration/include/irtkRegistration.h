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

#ifndef _IRTKREGISTRATION_H

#define _IRTKREGISTRATION_H

#define MAX_NO_RESOLUTIONS 10

// Definition of available states for individual DOFs
typedef enum { Active, Passive } DOFStatus;

// Definition of available optimization m
typedef enum { DownhillDescent,
               GradientDescent,
               GradientDescentConstrained,
               SteepestGradientDescent,
               ConjugateGradientDescent,
               ClosedForm
             } irtkOptimizationMethod;

// Definition of available similarity measures
typedef enum { JE, CC, MI, NMI, SSD, CR_XY, CR_YX, LC, K, ML, NGD, NGP, NGS }
irtkSimilarityMeasure;

#include <irtkImage.h>

#include <irtkHistogram.h>

#include <irtkResampling.h>

#include <irtkImageFunction.h>

#include <irtkTransformation.h>

#include <irtkSimilarityMetric.h>

#include <irtkOptimizer.h>

#include <irtkUtil.h>

class irtkRegistration : public irtkObject
{

public:

    irtkHistory *history;

  /// Evaluate similarity metric
  virtual double Evaluate() = 0;

  /// Evaluate gradient of similarity metric
  virtual double EvaluateGradient(float, float *) = 0;

};

double combine_mysimilarity(double,double,double,double);
double combine_mysimilarity(irtkSimilarityMetric **, double *, double);

#include <irtkPointRegistration.h>
#include <irtkSurfaceRegistration.h>
#include <irtkModelRegistration.h>
#include <irtkImageRegistration.h>
#include <irtkSymmetricImageRegistration.h>

#endif
