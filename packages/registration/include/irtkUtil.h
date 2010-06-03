/*=========================================================================

  Library   : Image Registration Toolkit (IRTK)
  Module    : $Id$
  Copyright : Imperial College, Department of Computing
              Visual Information Processing (VIP), 2008 onwards
  Date      : $Date$
  Version   : $Revision$
  Changes   : $Author$

=========================================================================*/

#ifndef _IRTKUTIL_H

#define _IRTKUTIL_H

/// Struct to keep registration parameters and associated similarities.
typedef struct _irtkHistoryEntry {

public:

  /// Parameters
  double *_parameter;

  /// Value of similarity measure
  double _similarity;

  /// Pointer to the next entry in the list
  struct _irtkHistoryEntry *_succ;

} irtkHistoryEntry;

/// History class to avoid loops in registration optimization.
class irtkHistory
{

  /// List of entries
  irtkHistoryEntry *_entries;

  /// Length of list
  int _no;

/// Some statistics
  int _success;

  /// Some statistics
  int _failed;

  /// Are two transformations the same
  int IsEqual(const irtkHistoryEntry *, const irtkTransformation *);

public:

  /// Constructor
  irtkHistory();

  /// Destrcutor
  ~irtkHistory();

  /// Add a transformation and its corresponding similarity value to history
  void Add(const irtkTransformation *, double);

  /** Check whether a transformation is in history. If not, return false
      otherwise return true and its corresponding similarity value.
  */
  int Find(const irtkTransformation *, double &);

  /// Print some debugging information and statistics
  void Print();

  /// Clear list
  void Clear();
};

inline irtkHistory::irtkHistory()
{
  _entries = NULL;
  _no      = 0;
  _success = 0;
  _failed  = 0;
}

inline void irtkHistory::Clear()
{
  irtkHistoryEntry *entry;

  entry = _entries;
  while (entry != NULL) {
    _entries = entry->_succ;
    delete [] entry->_parameter;
    delete entry;
    entry = _entries;
  }
  _entries = NULL;
  _no      = 0;
  _success = 0;
  _failed  = 0;
}

inline irtkHistory::~irtkHistory()
{
  this->Clear();
}

inline void irtkHistory::Add(const irtkTransformation *transformation,
                             double similarity)
{
  // Allocate a new entry
  irtkHistoryEntry *entry = new irtkHistoryEntry;

  // Create a new entry
  entry->_similarity = similarity;
  entry->_parameter  = new double[transformation->NumberOfDOFs()];
  for (int i = 0; i < transformation->NumberOfDOFs(); i++) {
    entry->_parameter[i] = transformation->Get(i);
  }
  entry->_succ = _entries;

  // Insert the new entry at the beginning of the list
  _entries = entry;
  _no++;
}

inline int  irtkHistory::IsEqual(const irtkHistoryEntry *entry,
                                 const irtkTransformation *transformation)
{
  for (int i = 0; i < transformation->NumberOfDOFs(); i++) {
    if (entry->_parameter[i] != transformation->Get(i)) return False;
  }
  return True;
}

inline int  irtkHistory::Find(const irtkTransformation *transformation,
                              double &similarity)
{
  for (irtkHistoryEntry *entry = _entries; entry != NULL; entry = entry->_succ) {
    if (IsEqual(entry, transformation) == True) {
      similarity = entry->_similarity;
      _success++;
      return True;
    }
  }
  _failed++;
  return False;
}

inline void irtkHistory::Print()
{
  if (_no > 0) {
    cout << "Cache has " << _no << " entries" << endl;
    cout << "Cache statistics: " << 100 * _success / (_success + _failed)
         << " % hits (" << _success << " out of " << (_success + _failed)
         << ")" << endl;
  }
}

extern void irtkPadding(irtkGreyImage &, irtkGreyPixel);
extern void irtkPadding(irtkGreyImage &, irtkGreyPixel, irtkFreeFormTransformation3D *ffd);
extern void irtkPadding(irtkGreyImage **, irtkGreyPixel, irtkFreeFormTransformation3D *ffd, int numberOfImages);
extern int  irtkCalculateNumberOfBins(irtkGreyImage *, int, int, int);
extern int  irtkCalculateNumberOfBins(irtkGreyImage **, int, int, int, int);
extern double GuessResolution(double, double);
extern double GuessResolution(double, double, double);
extern int GuessPadding(irtkGreyImage &);
extern int read_line(istream &, char *, char *&);


#ifdef HAS_VTK

#include <vtkPolyData.h>

void MarkBoundary(vtkPolyData *polydata);

#endif

#endif
