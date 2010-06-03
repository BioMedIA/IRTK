/*=========================================================================

  Library   : Image Registration Toolkit (IRTK)
  Module    : $Id$
  Copyright : Imperial College, Department of Computing
              Visual Information Processing (VIP), 2008 onwards
  Date      : $Date$
  Version   : $Revision$
  Changes   : $Author$

=========================================================================*/

#ifndef _IRTKMULTILEVELFREEFORMTRANSFORMATION_H

#define _IRTKMULTILEVELFREEFORMTRANSFORMATION_H

#include <irtkGeometry.h>

#define MAX_TRANS 200

/**
 * Class for free form transformations based on tensor product B-splines.
 *
 * This class implements 3D free form transformation using B-splines.
 *
 * For more details about the implementation see Lee, Wolberg and Shin, IEEE
 * Transactions on Visualization and Computer Graphics, Vol. 3, No. 3, 1997.
 */

class irtkMultiLevelFreeFormTransformation : public irtkAffineTransformation
{

public:

  /// Local transformations
  irtkFreeFormTransformation *_localTransformation[MAX_TRANS];

  /// Number of local transformations
  int _NumberOfLevels;

  /// Constructor (default)
  irtkMultiLevelFreeFormTransformation();

  /// Constructor (copy)
  irtkMultiLevelFreeFormTransformation(const irtkRigidTransformation &);

  /// Constructor (copy)
  irtkMultiLevelFreeFormTransformation(const irtkAffineTransformation &);

  /// Constructor (copy)
  irtkMultiLevelFreeFormTransformation(const irtkMultiLevelFreeFormTransformation &);

  /// Destructor
  virtual ~irtkMultiLevelFreeFormTransformation();

  /// Returns the number of levels
  virtual int NumberOfLevels() const;

  /// Gets local transformation
  virtual irtkFreeFormTransformation *GetLocalTransformation(int);

  /// Puts local transformation
  virtual void PutLocalTransformation(irtkFreeFormTransformation *, int);

  /// Push local transformation on stack
  virtual void PushLocalTransformation(irtkFreeFormTransformation *);

  /// Combine local transformation on stack
  virtual void CombineLocalTransformation();

  /// Pop local transformation from stack
  virtual irtkFreeFormTransformation *PopLocalTransformation();

  /// Transforms a point
  virtual void Transform(double &, double &, double &, double = 0);

  /// Calculates displacement using global and local transformation components
  virtual void Displacement(double& x, double& y, double& z, double = 0);

  /// Transforms a single point using the global transformation component only
  virtual void GlobalTransform(double &, double &, double &, double = 0);

  /// Transforms a single point using the local transformation component only
  virtual void LocalTransform (double &, double &, double &, double = 0);

  /// Calculates displacement using the global transformation component only
  virtual void GlobalDisplacement(double &, double &, double &, double = 0);

  /// Calculates displacement using the local transformation component only
  virtual void LocalDisplacement(double &, double &, double &, double = 0);

  /// Transforms a point
  virtual void Transform(int, double &, double &, double &, double = 0);

  /// Transforms a single point using the local transformation component only
  virtual void LocalTransform (int, double &, double &, double &, double = 0);

  /// Calculates displacement using the local transformation component only
  virtual void LocalDisplacement(int, double &, double &, double &, double = 0);

  /// Calculate the Jacobian of the transformation
  virtual void Jacobian(irtkMatrix &, double, double, double, double = 0);

  /// Calculate the Jacobian of the local transformation
  virtual void LocalJacobian(irtkMatrix &, double, double, double, double = 0);

  /** Approximate displacements: This function takes a set of points and a
      set of displacements and find a FreeFD which approximates these
      displacements. After approximatation the displacements replaced by
      the residual displacement errors at the points */
  virtual double Approximate(double *, double *, double *,
                             double *, double *, double *, int);

  /// Inverts the transformation
  virtual double Inverse(double &, double &, double &, double = 0, double = 0.01);

  /// Inverts the transformation
  virtual double Inverse(int, double &, double &, double &, double = 0, double = 0.01);

  /// Checks whether transformation is an identity mapping
  virtual Bool IsIdentity();

  /// Prints the parameters of the transformation
  virtual void Print();

  /// Check file header
  static int CheckHeader(char *);

  /// Returns a string with the name of the instantiated class
  virtual const char *NameOfClass();

  /// Reads a transformation from a file
  virtual irtkCifstream& Read(irtkCifstream&);

  /// Writes a transformation to a file
  virtual irtkCofstream& Write(irtkCofstream&);

  /// Imports a transformation from a file
  virtual istream& Import(istream&);

  /// Exports a transformation to a file
  virtual ostream& Export(ostream&);

};

inline int irtkMultiLevelFreeFormTransformation::NumberOfLevels() const
{
  return _NumberOfLevels;
}

inline irtkFreeFormTransformation *irtkMultiLevelFreeFormTransformation::GetLocalTransformation(int level)
{
  if (level < _NumberOfLevels) {
    return _localTransformation[level];
  } else {
    cerr << "irtkMultiLevelFreeFormTransformation::GetLocalTransformation: No such "
         << "transformation" << endl;
    exit(1);
  }
}

inline void irtkMultiLevelFreeFormTransformation::PutLocalTransformation(irtkFreeFormTransformation *transformation, int i)
{
  if (i < _NumberOfLevels) {
    _localTransformation[i] = transformation;
  } else {
    cerr << "irtkMultiLevelFreeFormTransformation::PutLocalTransformation: No such "
         << "transformation" << endl;
    exit(1);
  }
}

inline void irtkMultiLevelFreeFormTransformation::PushLocalTransformation(irtkFreeFormTransformation *transformation)
{
  if (_NumberOfLevels < MAX_TRANS) {
    _localTransformation[_NumberOfLevels] = transformation;
    _NumberOfLevels++;
  } else {
    cerr << "irtkMultiLevelFreeFormTransformation::PushLocalTransformation: Stack "
         << "overflow" << endl;
    exit(1);
  }
}

inline irtkFreeFormTransformation *irtkMultiLevelFreeFormTransformation::PopLocalTransformation()
{
  irtkFreeFormTransformation *localTransformation;

  if (_NumberOfLevels > 0) {
    localTransformation = _localTransformation[_NumberOfLevels-1];
    _localTransformation[_NumberOfLevels-1] = NULL;
    _NumberOfLevels--;
  } else {
    cerr << "irtkMultiLevelFreeFormTransformation:PopLocalTransformation: Stack "
         << "underflow" << endl;
    exit(1);
  }
  return localTransformation;
}

inline const char *irtkMultiLevelFreeFormTransformation::NameOfClass()
{
  return "irtkMultiLevelFreeFormTransformation";

}

#endif
