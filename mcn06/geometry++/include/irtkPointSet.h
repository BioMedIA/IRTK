/*=========================================================================

  Library   : Image Registration Toolkit (IRTK)
  Module    : $Id$
  Copyright : Imperial College, Department of Computing
              Visual Information Processing (VIP), 2008 onwards
  Date      : $Date$
  Version   : $Revision$
  Changes   : $Author$

=========================================================================*/

#ifndef _IRTKPOINTSET_H

#define _IRTKPOINTSET_H

#define POINTSET_SIZE 4096

/**

  Point set class.

*/

class irtkPointSet : public irtkObject
{

private:

  /// Allocated size of irtkPointSet
  int _m;

  /// Actual size of irtkPointSet
  int _n;

  /// Pointer to Points
  irtkPoint *_data;

public:

  //
  // Constructors and destructor
  //

  /// Default constructor
  irtkPointSet();

  /// Copy constructor
  irtkPointSet(const irtkPointSet &);

  /// Destructor
  virtual ~irtkPointSet();

  //
  // Size get acccessor and clearing of irtkPointSet
  //

  /// Access function for size
  int  Size() const;

  /// Clearing of irtkPointSet
  void Clear();

  //
  // Operators for access
  //

  /// Operator for Point put access
  irtkPoint   &operator()(int);

  /// Operator for Point get access
  irtkPoint    operator()(int) const;

  /// Operator
  irtkPointSet operator()(int, int) const;

  //
  // Operators for irtkPointSet arithmetic
  //

  // irtkPointSet operation for =
  irtkPointSet &operator=(const irtkPointSet &);

  // irtkPointSet operator for Point adding
  irtkPointSet& operator+=(const irtkPoint&);

  // irtkPointSet operator for Point substraction
  irtkPointSet& operator-=(const irtkPoint&);

  // irtkPointSet operator for irtkPointSet adding
  irtkPointSet& operator+=(const irtkPointSet&);

  // irtkPointSet operator for irtkPointSet substraction
  irtkPointSet& operator-=(const irtkPointSet&);

  /// Centre of gravity
  virtual irtkPoint CenterOfGravity() const;

  /// Bounding box
  virtual void BoundingBox(irtkPoint &, irtkPoint &) const;

  //
  // Explicit methods to add or delete points
  //

  /// Adding of a Point to Pointset
  void Add(const irtkPoint&);

  /// Deleting of a Point from Pointset
  void Del(const irtkPoint&);

  /// Adding of a irtkPointSet to Pointset
  void Add(const irtkPointSet &);

  /// Deleting of a irtkPointSet from Pointset
  void Del(const irtkPointSet &);

  /// Adding of a Point to Pointset
  void Add(double *);

  /// Deleting of a Point from Pointset
  void Del(double *);

  //
  // irtkPointSet in- and output
  //

  /// Interface to output stream
  friend ostream& operator<< (ostream&, const irtkPointSet&);

  /// Interface to input stream
  friend istream& operator>> (istream&, irtkPointSet&);

  /// Read pointset from file
  void Read (char *);

  /// Write pointset to file
  void Write(char *);

  /// Read pointset from file in VTK format
  void ReadVTK (char *);

  /// Write pointset to file in VTK format
  void WriteVTK(char *);

  //
  // Misc. functions
  //

  /// Computes the standard deviation ellipsoid about the centroid of a point set.
  irtkPoint StandardDeviationEllipsoid() const;

  /// Tests if a point is inside the polygon defined by the point set
  int IsInside(double, double) const;

};

inline irtkPointSet::irtkPointSet()
{
  // Initialize
  _n = 0;
  _m = POINTSET_SIZE;
  _data = new irtkPoint[POINTSET_SIZE];
}

inline irtkPointSet::irtkPointSet(const irtkPointSet &pset)
{
  int i;

  // Initialize
  _n = 0;
  _m = 0;
  _data = NULL;

  // Allocate memory
  if (pset._m  > 0) {
    _data = new irtkPoint[pset._m];
    _m = pset._m;
    _n = pset._n;
  }

  // Copy points
  for (i = 0; i < _n; i++) {
    _data[i] = pset._data[i];
  }
}

inline irtkPointSet::~irtkPointSet()
{
  if (_m > 0) {
    delete []_data;
    _m = 0;
    _n = 0;
  }
  _m = 0;
  _n = 0;
  _data = NULL;
}

inline irtkPoint& irtkPointSet::operator()(int j)
{
#ifdef NO_BOUNDS
  return _data[j];
#else
  if ((j >= 0) && (j < _n)) {
    return _data[j];
  } else {
    cerr << "irtkPointSet::operator(int) parameter out of range\n";
    exit(1);
  }
#endif
}

inline irtkPoint  irtkPointSet::operator()(int j) const
{
#ifdef NO_BOUNDS
  return _data[j];
#else
  if ((j >= 0) && (j < _n)) {
    return _data[j];
  } else {
    cerr << "irtkPointSet::operator(int) parameter out of range\n";
    exit(1);
  }
#endif
}

inline irtkPointSet& irtkPointSet::operator+=(const irtkPoint &p)
{
  this->Add(p);
  return *this;
}

inline irtkPointSet& irtkPointSet::operator-=(const irtkPoint &p)
{
  this->Del(p);
  return *this;
}

inline irtkPointSet& irtkPointSet::operator+=(const irtkPointSet &pset)
{
  this->Add(pset);
  return *this;
}

inline irtkPointSet& irtkPointSet::operator-=(const irtkPointSet &pset)
{
  this->Del(pset);
  return *this;
}

inline irtkPointSet irtkPointSet::operator()(int j, int k) const
{
  int i;
  irtkPointSet pset;

#ifdef NO_BOUNDS
  for (i = j; i < k; j++) {
    pset += _data[i];
  }
  return pset;
#else
  if ((j >= 0) && (k < _n)) {
    for (i = j; i < k; j++) {
      pset += _data[i];
    }
    return pset;
  } else {
    cerr << "irtkPointSet::operator(int, int) parameter out of range\n";
    exit(1);
  }
  return pset;
#endif
}

inline int irtkPointSet::Size() const
{
  return _n;
}

inline void irtkPointSet::Add(double *p)
{
  this->Add(irtkPoint(p[0], p[1], p[2]));
}

inline void irtkPointSet::Del(double *p)
{
  this->Del(irtkPoint(p[0], p[1], p[2]));
}

#endif

