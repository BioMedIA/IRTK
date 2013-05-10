/*=========================================================================

  Library   : Image Registration Toolkit (IRTK)
  Module    : $Id$
  Copyright : Imperial College, Department of Computing
              Visual Information Processing (VIP), 2008 onwards
  Date      : $Date$
  Version   : $Revision$
  Changes   : $Author$

=========================================================================*/

#include <irtkBSplineFunction.h>

double irtkBSplineFunction::LookupTable   [irtkBSplineFunction::LookupTableSize][4];
double irtkBSplineFunction::LookupTable_I [irtkBSplineFunction::LookupTableSize][4];
double irtkBSplineFunction::LookupTable_II[irtkBSplineFunction::LookupTableSize][4];
bool   irtkBSplineFunction::_initialized = false;
