/*=========================================================================

  Library   : Image Registration Toolkit (IRTK)
  Module    : $Id$
  Copyright : Imperial College, Department of Computing
              Visual Information Processing (VIP), 2008 onwards
  Date      : $Date$
  Version   : $Revision$
  Changes   : $Author$

=========================================================================*/

#ifndef _IRTKEIGENFREEFORMTRANSFORMATION_H

#define _IRTKEIGENFREEFORMTRANSFORMATION_H

#include <irtkGeometry.h>

#undef NORMAL

/**
 * Adaptive FFD class definition based on eigen parameterization.
 *
 * This class implements adaptive hierarchical FFD's using
 * the eigenmodes (e.g. derived from biomechanics or statistics).
 *
 */

class irtkEigenFreeFormTransformation : public irtkBSplineFreeFormTransformation3D
{

protected:

  /// Eigenvalues
  irtkVector *_EigenValues;

  /// Eigenvectors
  irtkMatrix *_EigenVectors;

  /// Shape vector
  irtkVector *_ShapeVector;

  /// Label of each element (0=unlabelled by default)
  int *_label;

public:

  //
  // Constructors and destructor
  //

  /// Default constructor
  irtkEigenFreeFormTransformation();

  /// Constructor (FFD, eigen system, shape vector)
  irtkEigenFreeFormTransformation(const irtkBSplineFreeFormTransformation3D &,
                                  const irtkVector&,
                                  const irtkMatrix&,
                                  const irtkVector&);

  /// Copy Constructor
  irtkEigenFreeFormTransformation(const irtkEigenFreeFormTransformation &);

  /// Destructor
  virtual ~irtkEigenFreeFormTransformation();

  // Access parameters
  virtual SetMacro(EigenValues,  irtkVector *);
  virtual GetMacro(EigenValues,  irtkVector *);
  virtual SetMacro(EigenVectors, irtkMatrix *);
  virtual GetMacro(EigenVectors, irtkMatrix *);
  virtual SetMacro(ShapeVector,  irtkVector *);
  virtual GetMacro(ShapeVector,  irtkVector *);

  //
  // Methods to manipulate the status/labels of the control points
  //

  /// Gets a DOF status (overloaded, always returns _Active)
  virtual _Status GetStatus(int) const;

  /// Returns the label of the indexed control point
  virtual int GetLabel(int);

  /// Sets the label of the indexed control point
  void PutLabel(int, int);

  //
  // Methods for element access
  //

  /// Returns the number of elements of the transformation
  virtual int  NumberOfElements() const;

  /// Gets an element by filling a pre-allocated array with node indices
  virtual void GetElement(int, int*);

  /// Gets an element by allocating array and performing overloaded call
  virtual int *GetElement(int);

  //
  // Inherited methods from irtkTransformation to be implemented
  //

  /// Returns the number of DOFs (rows of _ShapeVector)
  virtual int    NumberOfDOFs() const;

  /// Gets a shape parameter
  virtual double Get(int) const;

  /// Puts a shape parameter
  virtual void   Put(int, double);

  /// Prints the transformation info
  void Info();

  /// Prints the parameters of the transformation
  virtual void   Print();

  /// Check file header
  static int     CheckHeader(char *);

  /// Check keyword from file header
  static int     CheckKeyword(char *);

  /// Returns a string with the name of the instantiated class
  virtual const char  *NameOfClass();

  /// Reads a transformation from a file
  virtual irtkCifstream& Import(irtkCifstream&);

  /// Writes a transformation to a file
  virtual irtkCofstream& Export(irtkCofstream&);

  /// Imports a transformation from a file
  virtual istream& Import(istream&);

  /// Exports a transformation to a file
  virtual ostream& Export(ostream&);
};

//
// Instantiated inline methods
//

inline int irtkEigenFreeFormTransformation::NumberOfElements() const
{
  if (_z > 1) {
    return (_x - 1) * (_y - 1) * (_z - 1);
  } else {
    return (_x - 1) * (_y - 1);
  }
}

inline void irtkEigenFreeFormTransformation::GetElement(int index, int *element)
{
  int i, j, k;

  if ((index < 0) || (index >= this->NumberOfElements())) {
    cerr << "irtkEigenFreeFormTransformation::GetElement() : index out of range\n";
    exit(0);
  }
  if (element == NULL) {
    cerr << "irtkEigenFreeFormTransformation::GetElement() : invalid field passsed\n";
    exit(0);
  }

  // Find top left node of element (=N1)
  // See vtkStructuredGrid::GetCell()
  if (_z > 1) {
    // 3D (case VTK_XYZ_GRID) - still needs checking!
    i =  index % ( _x - 1);
    j = (index / ( _x - 1)) % (_y - 1);
    k =  index / ((_x - 1)  * (_y - 1));
  } else {
    // 2D (case VTK_XY_PLANE)
    i = index % ( _x - 1);
    j = index / ( _x - 1);
    k = 0;
  }

  // In-plane element node indices
  element[0] = this->LatticeToIndex(i+1, j,   k);
  element[1] = this->LatticeToIndex(i,   j,   k);
  element[2] = this->LatticeToIndex(i,   j+1, k);
  element[3] = this->LatticeToIndex(i+1, j+1, k);

  // Through-plane element node indices, if applicable
  if (_z > 1) {
    element[4] = this->LatticeToIndex(i+1, j,   k+1);
    element[5] = this->LatticeToIndex(i,   j,   k+1);
    element[6] = this->LatticeToIndex(i,   j+1, k+1);
    element[7] = this->LatticeToIndex(i+1, j+1, k+1);
  }

  return;
}

inline int *irtkEigenFreeFormTransformation::GetElement(int index)
{
  int *element = NULL;

  // Allocate 2D or 3D element
  if (_z > 1) {
    element = new int[8];
  } else {
    element = new int[4];
  }

  // Overloaded call
  this->GetElement(index, element);

  // Return element
  return element;
}

inline _Status irtkEigenFreeFormTransformation::GetStatus(int index) const
{
  return _Active;
}

inline int irtkEigenFreeFormTransformation::GetLabel(int index)
{
  // Bound checking
  if ((index < 0) || (index >= this->NumberOfElements())) {
    cerr << "irtkEigenFreeFormTransformation::GetLabel(): index out of range\n";
    exit(1);
  }

  // One label per element
  return _label[index];
}

inline void irtkEigenFreeFormTransformation::PutLabel(int index, int label)
{
  // Bound checking
  if ((index < 0) || (index >= this->NumberOfElements())) {
    cerr << "irtkEigenFreeFormTransformation::PutLabel(): index out of range\n";
    exit(1);
  }

  // One label per element
  _label[index] = label;
}

inline int irtkEigenFreeFormTransformation::NumberOfDOFs() const
{
  return _ShapeVector->Rows();
}

inline double irtkEigenFreeFormTransformation::Get(int index) const
{
  if (index >= this->NumberOfDOFs()) {
    cerr << "irtkEigenFreeFormTransformation::Get: No such dof" << endl;
    exit(1);
  }
  return _ShapeVector->Get(index);
}

inline int irtkEigenFreeFormTransformation::CheckKeyword(char *header)
{
  if (strcmp(header, "EFFD:") == 0) {
    return True;
  } else {
    return False;
  }
}

inline const char *irtkEigenFreeFormTransformation::NameOfClass()
{
  return "irtkEigenFreeFormTransformation";
}

#endif
