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

// Image Registration Toolkit (IRTK) includes
#include <irtkRegistration.h>

irtkImageEigenFreeFormRegistration2D::irtkImageEigenFreeFormRegistration2D() :
    irtkImageFreeFormRegistration2D()
{
  // Print debugging information
  this->Debug("irtkImageEigenFreeFormRegistration2D::irtkImageEigenFreeFormRegistration2D");

  // Initialize fields
  _NumberOfModes = 0;
  _effd = NULL;
}

irtkImageEigenFreeFormRegistration2D::~irtkImageEigenFreeFormRegistration2D()
{
  // empty (_effd allocated externally)
}

// Initial set up for the registration at a multiresolution level
void irtkImageEigenFreeFormRegistration2D::Initialize(int level)
{
  // Print debugging information
  this->Debug("irtkImageEigenFreeFormRegistration2D::Initialize(int)");

  // Call parent class
  this->irtkImageFreeFormRegistration2D::Initialize(level);

  // Check level
  _effd = (irtkEigenFreeFormTransformation *)_affd;
  if (_effd == NULL) {
    cerr << this->NameOfClass()
         << "::Initialize(int): Invalid level\n";
    exit(1);
  }

  // Check level parameters
  if (_NumberOfModes == 0) {
    _NumberOfModes = _effd->NumberOfDOFs();
  } else if (_effd->NumberOfDOFs() != _NumberOfModes) {
    irtkVector *ShapeVector = _effd->GetShapeVector();
    ShapeVector->Initialize(_NumberOfModes);
  }
  cout << "Number of Modes: " << _NumberOfModes << endl; cout.flush();
}

// Final set up for the registration at a multiresolution level
void irtkImageEigenFreeFormRegistration2D::Finalize(int level)
{
  // Print debugging information
  this->Debug("irtkImageEigenFreeFormRegistration2D::Finalize(int)");

  // Finalize base class
  this->irtkImageRegistration::Finalize(level);
}

double irtkImageEigenFreeFormRegistration2D::EvaluateGradient(float step, float *dx)
{
  int i;
  double s1, s2, norm, parameterValue;

  for (i = 0; i < _transformation->NumberOfDOFs(); i++) {
    if (_transformation->irtkTransformation::GetStatus(i) == _Active) {
      parameterValue = _transformation->Get(i);
      _transformation->Put(i, parameterValue + step);
      s1 = this->Evaluate();
      _transformation->Put(i, parameterValue - step);
      s2 = this->Evaluate();
      _transformation->Put(i, parameterValue);
      dx[i] = s1 - s2;
    } else {
      dx[i] = 0;
    }
  }

  // Calculate norm of vector
  norm = 0;
  for (i = 0; i < _transformation->NumberOfDOFs(); i++) {
    norm += dx[i] * dx[i];
  }

  // Normalize vector
  norm = sqrt(norm);
  if (norm > 0) {
    for (i = 0; i < _transformation->NumberOfDOFs(); i++) {
      dx[i] /= norm;
    }
  } else {
    for (i = 0; i < _transformation->NumberOfDOFs(); i++) {
      dx[i] = 0;
    }
  }

  return norm;
}
