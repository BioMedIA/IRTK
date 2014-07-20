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

#include <irtkLocator.h>

#include <vtkKDTreePointLocator.h>

irtkLocator::irtkLocator()
{
  cell = vtkGenericCell::New();
}

int irtkLocator::FindClosestPoint(double *xyz)
{
  if (loc_type == 0) {
    cell_locator->FindClosestPoint(xyz, &closestPoint[0], cell, (vtkIdType&) cellId, subId,
                                   dist2);
    xyz[0] = closestPoint[0];
    xyz[1] = closestPoint[1];
    xyz[2] = closestPoint[2];
    return cellId;
  } else if (loc_type == 1) {
    temp_id = point_locator->FindClosestPoint(xyz);
    _dataset->GetPoint(temp_id, xyz);
    return temp_id;
  } else if (loc_type == 2) {
    temp_id = kd_locator->FindClosestPoint(xyz);
    _dataset->GetPoint(temp_id, xyz);
    return temp_id;
  } else {
    return 0;
  }
}

void irtkLocator::SelectLocatorType(int type)
{
  if (type==0) {
    cell_locator = vtkCellLocator::New();
    cell_locator->CacheCellBoundsOn();
    cell_locator->SetNumberOfCellsPerBucket(5);
    loc_type=0;
  } else if (type==1) {
    point_locator = vtkPointLocator::New();
    point_locator->SetNumberOfPointsPerBucket(5);
    loc_type=1;
  } else if (type==2) {
#if VTK_MAJOR_VERSION >= 6
    kd_locator = vtkKDTreePointLocator::New();
#else
    kd_locator = vtkKDTreePointLocator::New();
#endif

    kd_locator->SetNumberOfPointsPerBucket(50);
    loc_type=2;
  } else {
    cerr << "Unkown locator" << endl;
    exit(1);
  }
}

void irtkLocator::SetElementsPerBucket(int elements)
{
  if (loc_type == 0) {
    cell_locator->SetNumberOfCellsPerBucket(elements);
  }
  if (loc_type == 1) {
    point_locator->SetNumberOfPointsPerBucket(elements);
  }
  if (loc_type == 2) {
    kd_locator->SetNumberOfPointsPerBucket(elements);
  }
}

void irtkLocator::SetDataSet(vtkPolyData *dataset)
{
  if (loc_type == 0) {
    cell_locator->SetDataSet(dataset);
    cell_locator->BuildLocator();
    _dataset = dataset;
  }
  if (loc_type == 1) {
    point_locator->SetDataSet(dataset);
    point_locator->BuildLocator();
    _dataset = dataset;
  }
  if (loc_type == 2) {
    kd_locator->SetDataSet(dataset);
    kd_locator->BuildLocator();
    _dataset = dataset;
  }
}

const char *irtkLocator::NameOfClass()
{
  if (loc_type == 0) {
    return cell_locator->GetClassName();
  } else if (loc_type == 1) {
    return point_locator->GetClassName();
  } else if (loc_type == 2) {
    return kd_locator->GetClassName();
  } else return "Unspecified locator";
}

#endif
