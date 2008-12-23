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

#ifndef _IRTKSURFACEFREEFORMREGISTRATION_H

#define _IRTKSURFACEFREEFORMREGISTRATION_H

class irtkSurfaceFreeFormRegistration : public irtkSurfaceRegistration
{

protected:

  /// Initial set up for the registration
  virtual void Initialize();

  /// Final set up for the registration
  virtual void Finalize();

public:

  /// Constructor
  irtkSurfaceFreeFormRegistration();

  /// Sets output for the registration filter
  virtual void SetOutput(irtkTransformation *);

  /// Returns the name of the class
  virtual const char *NameOfClass();
};

#endif

#endif
