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

irtkSurfaceFreeFormRegistration::irtkSurfaceFreeFormRegistration (): irtkSurfaceRegistration()
{
}

const char *irtkSurfaceFreeFormRegistration::NameOfClass ()
{
  return "irtkSurfaceFreeFormRegistration";
}

void irtkSurfaceFreeFormRegistration::SetOutput(irtkTransformation *transformation)
{
  if (strcmp(transformation->NameOfClass(),
             "irtkMultiLevelFreeFormTransformation") != 0) {
    cerr << "irtkSurfaceFreeFormRegistration::SetOutput: Transformation must be "
         << "irtkMultiLevelFreeFormTransformation" << endl;
    exit(0);
  }
  _transformation = transformation;
}

void irtkSurfaceFreeFormRegistration::Initialize ()
{
  // Initialize base class
  this->irtkSurfaceRegistration::Initialize ();

  // Create point-based registration
  _preg = new irtkPointFreeFormRegistration;
}

void irtkSurfaceFreeFormRegistration::Finalize ()
{
  // Delete point registration
  delete _preg;
}

#endif
