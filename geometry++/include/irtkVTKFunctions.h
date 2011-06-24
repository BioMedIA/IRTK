/*=========================================================================

  Library   : Image Registration Toolkit (IRTK)
  Module    : $Id: irtkUtil.h 250 2010-11-10 22:54:31Z ws207 $
  Copyright : Imperial College, Department of Computing
              Visual Information Processing (VIP), 2008 onwards
  Date      : $Date: 2010-11-10 22:54:31 +0000 (三, 10 十一月 2010) $
  Version   : $Revision: 250 $
  Changes   : $Author: ws207 $

=========================================================================*/

#ifndef _IRTKVTKFUNCTIONS_H

#define _IRTKVTKFUNCTIONS_H


#ifdef HAS_VTK

#include <vtkPointData.h>
#include <vtkPolyData.h>
#include <vtkPolyDataReader.h>
#include <vtkPolyDataWriter.h>
#include <vtkDoubleArray.h>
#include <vtkPointLocator.h>
#include <vtkIdList.h>
#include <vtkSmartPointer.h>
#include <vtkCell.h>
#include <vtkDelaunay2D.h>

void GetConnectedVertices(vtkSmartPointer<vtkPolyData> mesh, int seed, vtkSmartPointer<vtkIdList> connectedVertices);

#endif

#endif
