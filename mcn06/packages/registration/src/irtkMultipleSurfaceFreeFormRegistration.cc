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

#include <irtkMultipleSurfaceRegistration.h>

irtkMultipleSurfaceFreeFormRegistration::irtkMultipleSurfaceFreeFormRegistration (): irtkMultipleSurfaceRegistration()
{
}

const char *irtkMultipleSurfaceFreeFormRegistration::NameOfClass ()
{
  return "irtkMultipleSurfaceFreeFormRegistration";
}

void irtkMultipleSurfaceFreeFormRegistration::SetOutput(irtkTransformation *transformation)
{
  if (strcmp(transformation->NameOfClass(),
             "irtkMultiLevelFreeFormTransformation") != 0) {
    cerr << "irtkMultipleSurfaceFreeFormRegistration::SetOutput: Transformation must be "
         << "irtkMultiLevelFreeFormTransformation" << endl;
    exit(0);
  }
  _transformation = transformation;
}

void irtkMultipleSurfaceFreeFormRegistration::Initialize ()
{
  // Initialize base class
  this->irtkMultipleSurfaceRegistration::Initialize ();

  // Create point-based registration
  _preg = new irtkPointFreeFormRegistration;
}

void irtkMultipleSurfaceFreeFormRegistration::Finalize ()
{
  // Delete point registration
  delete _preg;
}

#endif
