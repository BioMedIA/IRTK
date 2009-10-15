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

#ifndef _IRTKSURFACERIGIDREGISTRATION_H

#define _IRTKSURFACERIGIDREGISTRATION_H

#include <irtkImage.h>
#include <irtkTransformation.h>
#include <irtkPointRegistration.h>
#include <irtkSurfaceRegistration.h>

#include <vtkCellLocator.h>
#include <vtkPointLocator.h>
#include <irtkLocator.h>
#include <vtkFeatureEdges.h>

class irtkSurfaceRigidRegistration : public irtkSurfaceRegistration
{

protected:

  /// Initial set up for the registration
  virtual void Initialize();

  /// Final set up for the registration
  virtual void Finalize();

public:

  /// Constructor
  irtkSurfaceRigidRegistration();

  /// Sets output for the registration filter
  virtual void SetOutput(irtkTransformation *);

  /// Returns the name of the class
  virtual const char *NameOfClass();

};

#endif

#endif
