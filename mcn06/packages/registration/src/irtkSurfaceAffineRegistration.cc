/*=========================================================================

  Library   : Image Registration Toolkit (IRTK)
  Module    : $Id$
  Copyright : Imperial College, Department of Computing
              Visual Information Processing (VIP), 2008 onwards
  Date      : $Date$
  Version   : $Revision$
  Changes   : $Author$

=========================================================================*/

#ifdef HAS_VTK

#include <irtkSurfaceRegistration.h>

irtkSurfaceAffineRegistration::irtkSurfaceAffineRegistration (): irtkSurfaceRegistration()
{
}

const char *irtkSurfaceAffineRegistration::NameOfClass ()
{
  return "irtkSurfaceAffineRegistration";
}

void irtkSurfaceAffineRegistration::SetOutput (irtkTransformation * transformation)
{
  if (strcmp (transformation->NameOfClass (), "irtkAffineTransformation") != 0) {
    cerr << this->NameOfClass ()
         << "::SetOutput: Transformation must be affine" << endl;
    exit (0);
  }
  _transformation = transformation;
}

void irtkSurfaceAffineRegistration::Initialize ()
{
  // Initialize base class
  this->irtkSurfaceRegistration::Initialize ();

  // Create point-based registration
  _preg = new irtkPointAffineRegistration;
}

void irtkSurfaceAffineRegistration::Finalize()
{
  // Delete point registration
  delete _preg;
}

#endif
