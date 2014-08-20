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

#ifdef HAS_VTK

#include <vtkPolyData.h>
#include <vtkCellLocator.h>
#include <vtkPointLocator.h>
#include "vtkKDTreePointLocator.h"
#include <vtkGenericCell.h>
#include <irtkImage.h>

#ifndef _IRTKLOCATOR_H

#define _IRTKLOCATOR_H

class irtkLocator : public irtkObject
{

protected:

  /// VTK cell locator
  vtkCellLocator *cell_locator;

  /// VTK point locator
  vtkPointLocator *point_locator;

  /// VTK KD-tree point locator
#if VTK_MAJOR_VERSION >= 6
  vtkKDTreePointLocator *kd_locator;
#else
  vtkKDTreePointLocator *kd_locator;
#endif
 
  /// VTK Cell
  vtkGenericCell *cell;

  /// Cell id
  int cellId;

  /// Sub id
  int subId;

  /// Distance
  double dist2;

  /// Coordinates of closest point
  double closestPoint[3];

  /// Locator type
  int loc_type;

  /// Pointer to VTK dataset
  vtkPolyData *_dataset;

  /// Temp id
  int temp_id;

public:

  /// Default contructor
  irtkLocator(void);

  /// Set number of elements per bucket
  void SetElementsPerBucket(int);

  /// Set the locator type
  void SelectLocatorType(int type);

  int GetLocatorType();

  /// Set dataset
  void SetDataSet(vtkPolyData *);

  /// Find the closest point
  int FindClosestPoint(double *xyz);

  /// Return name of class
  const char *NameOfClass();
};

inline int irtkLocator::GetLocatorType()
{
  return loc_type;

}

#endif

#endif
