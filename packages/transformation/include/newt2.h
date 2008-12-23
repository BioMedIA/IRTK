/*=========================================================================

  Library   : Image Registration Toolkit (IRTK)
  Module    : $Id$
  Copyright : Imperial College, Department of Computing
              Visual Information Processing (VIP), 2008 onwards
  Date      : $Date$
  Version   : $Revision$
  Changes   : $Author$

=========================================================================*/

extern void newt2(float x[], int n, int *check, void (*vecfunc)(int, float [], float []));

extern double x_invert, y_invert, z_invert;

extern irtkTransformation *irtkTransformationPointer;

extern void irtkTransformationEvaluate(int, float point[], float fval[]);
