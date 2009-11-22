/*=========================================================================

  Library   : Image Registration Toolkit (IRTK)
  Module    : $Id$
  Copyright : Imperial College, Department of Computing
              Visual Information Processing (VIP), 2008 onwards
  Date      : $Date$
  Version   : $Revision$
  Changes   : $Author$

=========================================================================*/

#ifndef _IRTKCONTOUR_H

#define _IRTKCONTOUR_H

#include <irtkImage.h>
#include <vector>

enum ContourMode { NewContour, NewPoint, CloseContour };

/// Class for storing the contour in viewer coordinates
class irtkContour
{

protected:

  /// Contour parts
  vector<irtkPointSet> _pointSets;
  irtkPointSet _allPoints;
  int _updateAllPoints;

public:

  /// Constructor
  irtkContour();

  /// Destructor
  virtual ~irtkContour() {};

  virtual void Add(irtkPoint p);
  int IsEmpty();
  int Size();
  void Clear();
  int IsInside(double x, double y);
  virtual irtkPoint   &operator()(int);
  void AddNewSet(irtkPoint p);
  virtual void DeleteLastSet();
  void Print();

private:

  void AllPoints();

protected:

  void AddPointSet();

};

#endif
