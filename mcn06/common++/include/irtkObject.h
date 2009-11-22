/*=========================================================================

  Library   : Image Registration Toolkit (IRTK)
  Module    : $Id$
  Copyright : Imperial College, Department of Computing
              Visual Information Processing (VIP), 2008 onwards
  Date      : $Date$
  Version   : $Revision$
  Changes   : $Author$

=========================================================================*/

#ifndef _IRTKOBJECT_H

#define _IRTKOBJECT_H

/**

  Object class.

  @author     (C)opyright Daniel Rueckert and Julia Schnabel 1994-1998++

*/

class irtkObject
{

public:

  /// Constructor
  irtkObject();

  /// Destructor
  virtual ~irtkObject();

  /// Set delete method needed for TCL wrapping
  void SetDeleteMethod(void (*f)(void *));

protected:

  /// Pointer to delete method needed for TCL wrapping
  void (*DeleteMethod)(void *);

};

#endif

