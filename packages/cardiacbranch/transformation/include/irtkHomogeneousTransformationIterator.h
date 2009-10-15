/*=========================================================================

  Library   : Image Registration Toolkit (IRTK)
  Module    : $Id$
  Copyright : Imperial College, Department of Computing
              Visual Information Processing (VIP), 2008 onwards
  Date      : $Date$
  Version   : $Revision$
  Changes   : $Author$

=========================================================================*/

#ifndef _IRTKHOMOGENEOUSTRANSFORMATION_ITERATOR_H

#define _IRTKHOMOGENEOUSTRANSFORMATION_ITERATOR_H

/**
 * Class for iterator for homogeneous matrix transformations.
 *
 * This class implements a fast access iterator 3D for homogeneous
 * matrix transformations.
 *
 * NOTE: This class has NO copy constructor
 */

class irtkHomogeneousTransformationIterator : public irtkPoint
{

  /// Pointer to transformation
  irtkHomogeneousTransformation *_transformation;

public:

  /// Current x,y,z position in x-direction
  double _xx, _xy, _xz;

  /// Current x,y,z position in y-direction
  double _yx, _yy, _yz;

  /// Current x,y,z position in z-direction
  double _zx, _zy, _zz;

  /// Current x,y,z offsets in x-direction
  double _xdx, _xdy, _xdz;

  /// Current x,y,z offsets in y-direction
  double _ydx, _ydy, _ydz;

  /// Current x,y,z offsets in z-direction
  double _zdx, _zdy, _zdz;

  /// Constructor
  irtkHomogeneousTransformationIterator(irtkHomogeneousTransformation * = NULL);

  /** Initialize iterator. This function initializes the transformation
   *  iterator at point x, y, z with voxel offsets 1, 1, 1 */
  void Initialize(irtkBaseImage *target, irtkBaseImage *source, double x = 0, double y = 0, double z = 0);

  /// Advance iterator in x-direction
  void NextX();

  /// Advance iterator in x-direction by a certain amount
  void NextX(double);

  /// Advance iterator in y-direction
  void NextY();

  /// Advance iterator in y-direction by a certain amount
  void NextY(double);

  /// Advance iterator in z-direction
  void NextZ();

  /// Advance iterator in z-direction by a certain amount
  void NextZ(double);

  /** Sets the transformation for the iterator. */
  void SetTransformation(irtkHomogeneousTransformation* pTransformation);
};

inline irtkHomogeneousTransformationIterator::irtkHomogeneousTransformationIterator(irtkHomogeneousTransformation *transformation)
{
  _transformation = transformation;
}

inline void irtkHomogeneousTransformationIterator::Initialize(irtkBaseImage *target, irtkBaseImage *source, double x, double y, double z)
{

  if (_transformation == NULL) {
    cout << "irtkHomogeneousTransformationIterator::Initialize(): Transformation has not been set." << endl;
    exit(1);
  }

  // Transform point
  irtkMatrix matrix = source->GetWorldToImageMatrix() * _transformation->GetMatrix() * target->GetImageToWorldMatrix();

  _x = matrix(0, 0) * x + matrix(0, 1) * y + matrix(0, 2) * z + matrix(0, 3);
  _y = matrix(1, 0) * x + matrix(1, 1) * y + matrix(1, 2) * z + matrix(1, 3);
  _z = matrix(2, 0) * x + matrix(2, 1) * y + matrix(2, 2) * z + matrix(2, 3);

  // Calculate starting point
  _xx = _yx = _zx = _x;
  _xy = _yy = _zy = _y;
  _xz = _yz = _zz = _z;

  // Calculate offsets
  _zdx = matrix(0, 2);
  _zdy = matrix(1, 2);
  _zdz = matrix(2, 2);
  _ydx = matrix(0, 1);
  _ydy = matrix(1, 1);
  _ydz = matrix(2, 1);
  _xdx = matrix(0, 0);
  _xdy = matrix(1, 0);
  _xdz = matrix(2, 0);
}

inline void irtkHomogeneousTransformationIterator::SetTransformation(irtkHomogeneousTransformation* transformation)
{
  _transformation = transformation;
}

inline void irtkHomogeneousTransformationIterator::NextX()
{
  _x = _xx += _xdx;
  _y = _xy += _xdy;
  _z = _xz += _xdz;
}

inline void irtkHomogeneousTransformationIterator::NextX(double offset)
{
  _x = _xx += _xdx * offset;
  _y = _xy += _xdy * offset;
  _z = _xz += _xdz * offset;
}

inline void irtkHomogeneousTransformationIterator::NextY()
{
  _yx += _ydx;
  _yy += _ydy;
  _yz += _ydz;
  _x = _xx = _yx;
  _y = _xy = _yy;
  _z = _xz = _yz;
}

inline void irtkHomogeneousTransformationIterator::NextY(double offset)
{
  _yx += _ydx * offset;
  _yy += _ydy * offset;
  _yz += _ydz * offset;
  _x = _xx = _yx;
  _y = _xy = _yy;
  _z = _xz = _yz;
}

inline void irtkHomogeneousTransformationIterator::NextZ()
{
  _zx += _zdx;
  _zy += _zdy;
  _zz += _zdz;
  _x = _xx = _yx = _zx;
  _y = _xy = _yy = _zy;
  _z = _xz = _yz = _zz;
}

inline void irtkHomogeneousTransformationIterator::NextZ(double offset)
{
  _zx += _zdx * offset;
  _zy += _zdy * offset;
  _zz += _zdz * offset;
  _x = _xx = _yx = _zx;
  _y = _xy = _yy = _zy;
  _z = _xz = _yz = _zz;
}

#endif
