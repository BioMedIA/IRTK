/*=========================================================================

  Library   : Image Registration Toolkit (IRTK)
  Module    : $Id$
  Copyright : Imperial College, Department of Computing
              Visual Information Processing (VIP), 2008 onwards
  Date      : $Date$
  Version   : $Revision$
  Changes   : $Author$

=========================================================================*/

#ifndef _IRTKSEGMENTTABLE_H

#define _IRTKSEGMENTTABLE_H

#include <irtkSegment.h>

class irtkSegmentTable
{

  friend class irtkRView;
  friend class irtkViewer;

protected:

  /// Segment Table
  irtkSegment _entry[SHRT_MAX+1];

public:

/// Constructor (basic)
  irtkSegmentTable();

  /// Constructor (default)
  irtkSegmentTable(int);

  /// Destructor
  virtual ~irtkSegmentTable();

  /// Size of lookup table
  int Size();

  /// Sets all values for a segment
  void Set(int, char*, unsigned char, unsigned char, unsigned char, double, int);

  /// Sets the label value of sid (index)
  void SetLabel(int, char*);

  /// Sets the value of r,g,b as unsigned char of sid (index)
  void SetColor(int, unsigned char, unsigned char, unsigned char);

  /// Sets the transparency of sid (index)
  void SetTrans(int, double);

  /// Sets the visibility
  void SetVisibility(int, int);

  /// Get values at index id (int)
  char *Get(int, unsigned char*, unsigned char*, unsigned char*, double*, int*) const;

  /// Get segment colours at index id (int)
  void GetColor(int, unsigned char*, unsigned char*, unsigned char*);

  /// Get segment transparency value (double)
  void GetTrans(int, double*);

  /// Get segment hex value (char *)
  void GetHex(int, char*);

  /// Get segment label (char *)
  char *GetLabel(int);

  /// Is method that finds whether a segment is visible
  int GetVisibility(int);

  /// Find if entry contains valid value
  int IsValid(int);

  /// Remove segment with index
  void Clear(int);

  /// Remove all segments
  void Clear();

  /// Reads a segmentTable from a file
  void Read(char *);

  /// Writes a segmentTable to a file
  void Write(char *);
};

inline void irtkSegmentTable::GetColor(int id, unsigned char* r, unsigned char* g, unsigned char* b)
{
  _entry[id].getColor(r, g, b);
}

inline void irtkSegmentTable::GetTrans(int id, double* d)
{
  *d = _entry[id].getTrans();
}

inline void irtkSegmentTable::GetHex(int id, char* h)
{
  _entry[id].getHex(h);
}

inline char *irtkSegmentTable::GetLabel(int id)
{
  return _entry[id].getLabel();
}

inline int irtkSegmentTable::GetVisibility(int id)
{
  return _entry[id].getVisibility();
}

inline int irtkSegmentTable::IsValid(int id)
{
  if (_entry[id].getLabel() != NULL) {
    return true;
  } else {
    return false;
  }
}

inline int irtkSegmentTable::Size()
{
  return SHRT_MAX+1;
}

#endif

