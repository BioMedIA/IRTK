/*=========================================================================

  Library   : Image Registration Toolkit (IRTK)
  Module    : $Id$
  Copyright : Imperial College, Department of Computing
              Visual Information Processing (VIP), 2008 onwards
  Date      : $Date$
  Version   : $Revision$
  Changes   : $Author$

=========================================================================*/

#include <irtkRegistration.h>

void irtkImageAffineRegistration::SetOutput(irtkTransformation *transformation)
{
  if (strcmp(transformation->NameOfClass(), "irtkAffineTransformation") != 0) {
    cerr << "irtkImageAffineRegistration::SetOutput: Transformation must be affine"
         << endl;
    exit(0);
  }
  _transformation = transformation;
}

