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
  irtkFreeFormTransformation *_localTransformation[MAX_TRANS+1];

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

  /// Push local transformation on stack (append transformation)
  virtual void PushLocalTransformation(irtkFreeFormTransformation *);

  /// Insert local transformation
  virtual void InsertLocalTransformation(irtkFreeFormTransformation *, int = 0);

  /// Combine local transformation on stack
  virtual void CombineLocalTransformation();

  /// Pop local transformation from stack (remove last transformation)
  virtual irtkFreeFormTransformation *PopLocalTransformation();

  /// Remove local transformation and return the pointer (need to be deleted if not used)
  virtual irtkFreeFormTransformation *RemoveLocalTransformation(int = 0);

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

  /** Convert the global transformation from a matrix representation to a
      FFD and incorporate it with any existing local displacement. **/
  virtual void MergeGlobalIntoLocalDisplacement();

  // Helper function for the above.
  virtual void InterpolateGlobalDisplacement(irtkBSplineFreeFormTransformation3D *f);

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

  /// Calculate the bending of the local transformation.
  virtual double Bending(double, double, double);

  /** Approximate displacements: This function takes a set of points and a
      set of displacements and find a FreeFD which approximates these
      displacements. After approximatation the displacements replaced by
      the residual displacement errors at the points */
  virtual double Approximate(double *, double *, double *,
                             double *, double *, double *, int);

  /** Approximate displacements: This function takes a set of points from a complete image
      and a set of displacements and find a FreeFD which approximates these
      displacements. After approximatation the displacements replaced by
      the residual displacement errors at the points */
  virtual void ApproximateAsNew(double *, double *, double *,
                             double *, double *, double *, int);

  /// Inverts the transformation
  virtual double Inverse(double &, double &, double &, double = 0, double = 0.01);

  /// Checks whether transformation is an identity mapping
  virtual bool IsIdentity();

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

inline void irtkMultiLevelFreeFormTransformation::InsertLocalTransformation(irtkFreeFormTransformation *transformation, int pos)
{
	int i;

  if (_NumberOfLevels < MAX_TRANS) {
  	for (i = pos; i < _NumberOfLevels + 1; i++) _localTransformation[i+1] = _localTransformation[i];
  	_localTransformation[pos] = transformation;
    _NumberOfLevels++;
  } else {
    cerr << "irtkMultiLevelFreeFormTransformation::Insert"
    		"LocalTransformation: Stack overflow" << endl;
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
    cerr << "irtkMultiLevelFreeFormTransformation:PopLocalTransformation: Stack underflow" << endl;
    exit(1);
  }
  return localTransformation;
}

inline irtkFreeFormTransformation *irtkMultiLevelFreeFormTransformation::RemoveLocalTransformation(int pos)
{
	int i;
  irtkFreeFormTransformation *localTransformation;

  if (_NumberOfLevels > pos) {
    localTransformation = _localTransformation[_NumberOfLevels-1];
  	for (i = _NumberOfLevels-1; i > pos; i--) _localTransformation[i-1] = _localTransformation[i];
    _localTransformation[_NumberOfLevels-1] = NULL;
    _NumberOfLevels--;
  } else {
    cerr << "irtkMultiLevelFreeFormTransformation:DeleteLocalTransformation: No such " << "transformation" << endl;
    exit(1);
  }
  return localTransformation;
}

inline const char *irtkMultiLevelFreeFormTransformation::NameOfClass()
{
  return "irtkMultiLevelFreeFormTransformation";

}

#endif
