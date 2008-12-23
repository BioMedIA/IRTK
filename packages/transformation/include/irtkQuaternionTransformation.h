/*=========================================================================

  Library   : Image Registration Toolkit (IRTK)
  Module    : $Id$
  Copyright : Imperial College, Department of Computing
              Visual Information Processing (VIP), 2008 onwards
  Date      : $Date$
  Version   : $Revision$
  Changes   : $Author$

=========================================================================*/

#ifndef _IRTKQUATERNIONTRANSFORMATION_H

#define _IRTKQUATERNIONTRANSFORMATION_H

#include <irtkTransformation.h>

class irtkQuaternionTransformation : public irtkTransformation
{
private:
  /** The origin about which the rotation is performed. */
  double _origin[3];

  /** The translation. */
  double _translation[3];

  /** The scaling in the alpha, beta, and gamma directions. */
  double _scale[3];

  /** The alpha axis. */
  double _alphaAxis[3];

  /** The beta axis. */
  double _betaAxis[3];

  /** The gamma axis. */
  double _gammaAxis[3];

  /** The alpha axis after it has been rotated by the transformation. */
  double _finalAlphaAxis[3];

  /** The beta axis after it has been rotated by the transformation. */
  double _finalBetaAxis[3];

  /** The gamma axis after it has been rotated by the transformation. */
  double _finalGammaAxis[3];

  /** The alpha angle of rotation (in radians). */
  double _alpha;

  /** The beta angle of rotation (in radians). */
  double _beta;

  /** The gamma angle of rotation (in radians). */
  double _gamma;

  /** The matrix used for the transformation. */
  irtkMatrix _matrix;

  /** The number of DOFs. */
  static const int NUMBER_OF_DOFS = 9;

protected:
  /** Updates the gamma axis. */
  void UpdateGammaAxis();

  /** Checks that the alpha and beta axes are not in the same direction and
      projects the beta axis into the plane normal to the alpha axis. Then
      calculates the gamma axis which is equal to the cross product of the
      alpha and beta axes. */
  void UpdateAxes();

  /** Normalizes a vector. */
  void Normalize(double* vec);

public:
  /** Default constructor. */
  irtkQuaternionTransformation();

  /** Creates a quaternion transformation with a particular origin and set of
      axes. */
  irtkQuaternionTransformation(double* origin, double* alphaAxis, double* betaAxis);

  /** Copies a transformation. */
  irtkQuaternionTransformation(const irtkQuaternionTransformation& trans);

  /** Assignment operator. */
  irtkQuaternionTransformation& operator=(const irtkQuaternionTransformation& trans);

  /** Sets the alpha angle of rotation. */
  void SetAlpha(double angle);

  /** Sets the beta angle of rotation. */
  void SetBeta(double angle);

  /** Sets the gamma angle of rotation. */
  void SetGamma(double angle);

  /** Gets the alpha angle of rotation. */
  double GetAlpha() const;

  /** Gets the beta angle of rotation. */
  double GetBeta() const;

  /** Gets the gamma angle of rotation. */
  double GetGamma() const;

  /** Sets the axes of rotation. */
  void SetAxes(double* alphaAxis, double* betaAxis);

  /** Returns the axes of rotation. */
  void GetAxes(double* alphaAxis, double* betaAxis) const;

  /** Returns the axes of rotation. */
  void GetAxes(double* alphaAxis, double* betaAxis, double* gammaAxis) const;

  /** Returns the axes of rotation after they have been rotated by the
      transformation. */
  void GetRotatedAxes(double* alphaAxis, double* betaAxis) const;

  /** Returns the axes of rotation after they have been rotated by the
      transformation. */
  void GetRotatedAxes(double* alphaAxis, double* betaAxis, double* gammaAxis) const;

  /** Sets the origin of rotation. */
  void SetOrigin(double* origin);

  /** Gets the origin of rotation. */
  void GetOrigin(double* origin) const;

  /** Sets the translation. */
  void SetTranslation(double* translation);

  /** Gets the translation. */
  void GetTranslation(double* translation) const;

  /** Sets the scale factors. */
  void SetScale(double* scale);

  /** Returns the scale factors. */
  void GetScale(double* scale) const;

  /** Resets the origin of the rotation to the final position of the
      origin after translation. The translation parameters are set to 0
      after calling this function. */
  void ResetOrigin();

  /** Returns the matrix for the transformation. */
  const irtkMatrix& GetMatrix() const;

  /** Updates the matrix for the rotation. */
  void UpdateMatrix();

  /** Returns the number of DOFs for the transformation. */
  virtual int NumberOfDOFs() const;

  /** Returns a parameter. */
  virtual double Get(int param) const;

  /** Puts a parameter. */
  virtual void Put(int param, double val);

  /** Transforms a point. */
  virtual void Transform(double& x, double& y, double& z);

  /// Transforms a point using the global transformation component only
  virtual void GlobalTransform(double &, double &, double &);

  /// Transforms a point using the local transformation component only
  virtual void LocalTransform (double &, double &, double &);

  /// Calculates displacement using the global transformation component only
  virtual void GlobalDisplacement(double &, double &, double &);

  /// Calculates displacement using the local transformation component only
  virtual void LocalDisplacement(double &, double &, double &);

  /** Inverts the transformation. */
  virtual void Invert();

  /** Computes the Jacobian of the transformation at a point. */
  virtual void Jacobian(double x, double y, double z, irtkMatrix& jac);

  /** Computes the determinant of the Jacobian at a point. */
  virtual double Jacobian(double x, double y, double z);

  /** Returns the local Jacobian of the transformation at a point.
      \param x The x-coordinate of the point.
      \param y The y-coordinate of the point.
      \param z The z-coordinate of the point.
      \param jac The matrix to return the Jacobian in. */

  virtual void LocalJacobian(double x, double y, double z, irtkMatrix& jac);

  /** Returns the global Jacobian of the transformation at a point.
      \param x The x-coordinate of the point.
      \param y The y-coordinate of the point.
      \param z The z-coordinate of the point.
      \param jac The matrix to return the Jacobian in. */
  virtual void GlobalJacobian(double x, double y, double z, irtkMatrix& jac);

  /** Returns true if the file passed to the function is an
      irtkQuaternionTransformation. */
  static int CheckHeader(char* pFileName);

  /** Returns true if the transformation is an identity transformation. */
  virtual Bool IsIdentity();

  /** Reads the transformation. */
  virtual irtkCifstream& Read(irtkCifstream& is);

  /** Imports the transformation. */
  virtual istream& Import(istream& is);

  /** Writes the transformation. */
  virtual irtkCofstream& Write(irtkCofstream& os);

  /** Exports the transformation. */
  virtual ostream& Export(ostream& os);

  /** Prints the parameters of the transformation. */
  virtual void Print();

  /** Returns the name of the class. */
  virtual const char* NameOfClass();
};

inline void irtkQuaternionTransformation::UpdateGammaAxis()
{
  _gammaAxis[0] = _alphaAxis[1]*_betaAxis[2] - _alphaAxis[2]*_betaAxis[1];
  _gammaAxis[1] = _alphaAxis[2]*_betaAxis[0] - _alphaAxis[0]*_betaAxis[2];
  _gammaAxis[2] = _alphaAxis[0]*_betaAxis[1] - _alphaAxis[1]*_betaAxis[0];
}

inline void irtkQuaternionTransformation::UpdateAxes()
{
  Normalize(_alphaAxis);

  UpdateGammaAxis();

  if (_gammaAxis[0] == 0 && _gammaAxis[1] == 0 && _gammaAxis[2] == 0) {
    std::stringstream mesg;
    mesg << "Alpha axis and beta axis are collinear:" << std::endl
    << "_alphaAxis = (" << _alphaAxis[0] << ", " << _alphaAxis[1] << ", " << _alphaAxis[2] << ")" << std::endl
    << "_betaAxis = (" << _betaAxis[0] << ", " << _betaAxis[1] << ", " << _betaAxis[2] << ")";
  }

  double dp = _alphaAxis[0]*_betaAxis[0] + _alphaAxis[1]*_betaAxis[1] +
              _alphaAxis[2]*_betaAxis[2];

  _betaAxis[0] -= dp*_alphaAxis[0];
  _betaAxis[1] -= dp*_alphaAxis[1];
  _betaAxis[2] -= dp*_alphaAxis[2];

  Normalize(_betaAxis);
  Normalize(_gammaAxis);
}

inline void irtkQuaternionTransformation::Normalize(double* vec)
{
  double len = sqrt(vec[0]*vec[0] + vec[1]*vec[1] + vec[2]*vec[2]);

  if (len == 0)
    throw irtkException("Tried to normalize vector of length 0.", __FILE__,
                        __LINE__);

  vec[0] /= len;
  vec[1] /= len;
  vec[2] /= len;
}

inline void irtkQuaternionTransformation::SetAlpha(double angle)
{
  _alpha = angle*M_PI/180;
  UpdateMatrix();
}

inline void irtkQuaternionTransformation::SetBeta(double angle)
{
  _beta = angle*M_PI/180;
  UpdateMatrix();
}

inline void irtkQuaternionTransformation::SetGamma(double angle)
{
  _gamma = angle*M_PI/180;
  UpdateMatrix();
}

inline double irtkQuaternionTransformation::GetAlpha() const
{
  return _alpha*180/M_PI;
}

inline double irtkQuaternionTransformation::GetBeta() const
{
  return _beta*180/M_PI;
}

inline double irtkQuaternionTransformation::GetGamma() const
{
  return _gamma*180/M_PI;
}

inline void irtkQuaternionTransformation::SetAxes(double* alphaAxis, double* betaAxis)
{
  for (int i = 0; i < 3; i++) {
    _alphaAxis[i] = alphaAxis[i];
    _betaAxis[i] = betaAxis[i];
  }

  UpdateAxes();
  UpdateMatrix();
}

inline void irtkQuaternionTransformation::GetAxes(double* alphaAxis, double* betaAxis) const
{
  for (int i = 0; i < 3; i++) {
    alphaAxis[i] = _alphaAxis[i];
    betaAxis[i] = _betaAxis[i];
  }
}

inline void irtkQuaternionTransformation::GetAxes(double* alphaAxis, double* betaAxis, double* gammaAxis) const
{
  for (int i = 0; i < 3; i++) {
    alphaAxis[i] = _alphaAxis[i];
    betaAxis[i] = _betaAxis[i];
    gammaAxis[i] = _gammaAxis[i];
  }
}

inline void irtkQuaternionTransformation::GetRotatedAxes(double* alphaAxis, double* betaAxis) const
{
  for (int i = 0; i < 3; i++) {
    alphaAxis[i] = _finalAlphaAxis[i];
    betaAxis[i] = _finalBetaAxis[i];
  }
}

inline void irtkQuaternionTransformation::GetRotatedAxes(double* alphaAxis, double* betaAxis, double* gammaAxis) const
{
  for (int i = 0; i < 3; i++) {
    alphaAxis[i] = _finalAlphaAxis[i];
    betaAxis[i] = _finalBetaAxis[i];
    gammaAxis[i] = _finalGammaAxis[i];
  }
}

inline void irtkQuaternionTransformation::SetOrigin(double* origin)
{
  _origin[0] = origin[0];
  _origin[1] = origin[1];
  _origin[2] = origin[2];

  UpdateMatrix();
}

inline void irtkQuaternionTransformation::GetOrigin(double* origin) const
{
  for (int i = 0; i < 3; i++)
    origin[i] = _origin[i];
}

inline void irtkQuaternionTransformation::SetTranslation(double* translation)
{
  _translation[0] = translation[0];
  _translation[1] = translation[1];
  _translation[2] = translation[2];

  UpdateMatrix();
}

inline void irtkQuaternionTransformation::GetTranslation(double* translation) const
{
  translation[0] = _translation[0];
  translation[1] = _translation[1];
  translation[2] = _translation[2];
}

inline void irtkQuaternionTransformation::SetScale(double* scale)
{
  _scale[0] = scale[0];
  _scale[1] = scale[1];
  _scale[2] = scale[2];
}

inline void irtkQuaternionTransformation::GetScale(double* scale) const
{
  scale[0] = _scale[0];
  scale[1] = _scale[1];
  scale[2] = _scale[2];
}

inline const irtkMatrix& irtkQuaternionTransformation::GetMatrix() const
{
  return _matrix;
}

inline int irtkQuaternionTransformation::NumberOfDOFs() const
{
  return NUMBER_OF_DOFS;
}

inline void irtkQuaternionTransformation::GlobalDisplacement(double &x, double &y, double &z)
{
  cerr << "irtkQuaternionTransformation::GlobalDisplacement: No implemented yet" << endl;
}

inline void irtkQuaternionTransformation::LocalDisplacement(double &x, double &y, double &z)
{
  cerr << "irtkQuaternionTransformation::LocalDisplacement: No implemented yet" << endl;
}

inline void irtkQuaternionTransformation::GlobalTransform(double &x, double &y, double &z)
{
  cerr << "irtkQuaternionTransformation::GlobalTransform: No implemented yet" << endl;
}

inline void irtkQuaternionTransformation::LocalTransform(double &x, double &y, double &z)
{
  cerr << "irtkQuaternionTransformation::LocalTransform: No implemented yet" << endl;
}

inline double irtkQuaternionTransformation::Get(int param) const
{
  std::stringstream mesg;

  switch (param) {
  case 0:
    return _alpha*180/M_PI;
    break;

  case 1:
    return _beta*180/M_PI;
    break;

  case 2:
    return _gamma*180/M_PI;
    break;

  case 3:
    return _translation[0];
    break;

  case 4:
    return _translation[1];
    break;

  case 5:
    return _translation[2];
    break;

  case 6:
    return _scale[0];
    break;

  case 7:
    return _scale[1];
    break;

  case 8:
    return _scale[2];
    break;

  default:
    mesg << "Invalid parameter: " << param;
    throw irtkException(mesg.str(), __FILE__, __LINE__);
    break;
  }
}

inline void irtkQuaternionTransformation::Put(int param, double val)
{
  std::stringstream mesg;

  switch (param) {
  case 0:
    _alpha = val*M_PI/180;
    break;

  case 1:
    _beta = val*M_PI/180;
    break;

  case 2:
    _gamma = val*M_PI/180;
    break;

  case 3:
    _translation[0] = val;
    break;

  case 4:
    _translation[1] = val;
    break;

  case 5:
    _translation[2] = val;
    break;

  case 6:
    _scale[0] = val;
    break;

  case 7:
    _scale[1] = val;
    break;

  case 8:
    _scale[2] = val;
    break;

  default:
    mesg << "Invalid parameter: " << param;
    throw irtkException(mesg.str(), __FILE__, __LINE__);
    break;
  }

  UpdateMatrix();
}

inline const char* irtkQuaternionTransformation::NameOfClass()
{
  return "irtkQuaternionTransformation";
}

#endif // __IRTKQUATERNIONTRANSFORMATION_H
