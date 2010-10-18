/*=========================================================================

  Library   : Image Registration Toolkit (IRTK)
  Module    : $Id$
  Copyright : Imperial College, Department of Computing
              Visual Information Processing (VIP), 2008 onwards
  Date      : $Date$
  Version   : $Revision$
  Changes   : $Author$

=========================================================================*/

#include <irtkQuaternion.h>

#include <irtkQuaternionTransformation.h>

const unsigned int MAX_BUFFER_LENGTH = 1024;

irtkQuaternionTransformation::irtkQuaternionTransformation()
{
  _origin[0] = 0;
  _origin[1] = 0;
  _origin[2] = 0;

  _translation[0] = 0;
  _translation[1] = 0;
  _translation[2] = 0;

  _alphaAxis[0] = 1;
  _alphaAxis[1] = 0;
  _alphaAxis[2] = 0;

  _betaAxis[0] = 0;
  _betaAxis[1] = 1;
  _betaAxis[2] = 0;

  _scale[0] = 1;
  _scale[1] = 1;
  _scale[2] = 1;

  UpdateAxes();

  _alpha = 0;
  _beta = 0;
  _gamma = 0;

  _status = new _Status[NumberOfDOFs()];
  for (int i = 0; i < NumberOfDOFs(); i++) {
    _status[i] = _Active;
  }

  UpdateMatrix();
}

irtkQuaternionTransformation::irtkQuaternionTransformation(double* origin, double* alphaAxis, double* betaAxis)
{
  for (int i = 0; i < 3; i++) {
    _origin[i] = origin[i];
    _translation[i] = 0;
    _alphaAxis[i] = alphaAxis[i];
    _betaAxis[i] = betaAxis[i];
    _scale[i] = 1;
  }

  UpdateAxes();

  _alpha = 0;
  _beta = 0;
  _gamma = 0;

  _status = new _Status[NumberOfDOFs()];
  for (int i = 0; i < NumberOfDOFs(); i++) {
    _status[i] = _Active;
  }

  UpdateMatrix();
}

irtkQuaternionTransformation::irtkQuaternionTransformation(const irtkQuaternionTransformation& trans) : irtkTransformation(trans)
{
  for (int i = 0; i < 3; i++) {
    _origin[i] = trans._origin[i];
    _translation[i] = trans._translation[i];
    _alphaAxis[i] = trans._alphaAxis[i];
    _betaAxis[i] = trans._betaAxis[i];
  }

  UpdateAxes();

  _alpha = trans._alpha;
  _beta = trans._beta;
  _gamma = trans._gamma;

  _status = new _Status[NumberOfDOFs()];
  for (int i = 0; i < NumberOfDOFs(); i++) {
    _status[i] = trans._status[i];
  }

  UpdateMatrix();
}

irtkQuaternionTransformation& irtkQuaternionTransformation::operator=(const irtkQuaternionTransformation& trans)
{
  if (this == &trans) {
    return *this;
  }

  for (int i = 0; i < 3; i++) {
    _origin[i] = trans._origin[i];
    _translation[i] = trans._translation[i];
    _alphaAxis[i] = trans._alphaAxis[i];
    _betaAxis[i] = trans._betaAxis[i];
  }

  UpdateAxes();

  _alpha = trans._alpha;
  _beta = trans._beta;
  _gamma = trans._gamma;

  if (_status != NULL) {
    delete[] _status;
  }
  _status = new _Status[NumberOfDOFs()];
  for (int i = 0; i < NumberOfDOFs(); i++) {
    _status[i] = trans._status[i];
  }

  UpdateMatrix();

  return *this;
}

void irtkQuaternionTransformation::UpdateMatrix()
{
  // Compute matrix needed to translate origin of rotation axes to global
  // origin.
  irtkMatrix trans1(4, 4);
  trans1.Ident();
  trans1(0, 3) = -_origin[0];
  trans1(1, 3) = -_origin[1];
  trans1(2, 3) = -_origin[2];

  // Compute the scaling matrix.
  irtkMatrix scale(4, 4);
  scale(0, 0) = _scale[0];
  scale(1, 1) = _scale[1];
  scale(2, 2) = _scale[2];
  scale(3, 3) = 1;

  // Compute the three matrices needed to represent the three rotations:
  //
  // (1) Rotation by alpha radians about alpha axis.
  // (2) Rotation by beta radians about beta axis after beta axis has been
  // rotated by first rotation.
  // (3) Rotation by gamma radians about gamma axis after gamma axis has been
  // rotated by first two rotations.
  irtkMatrix rotMatrix1 = irtkQuaternion::ComputeRotationMatrix(_alpha, _alphaAxis[0], _alphaAxis[1], _alphaAxis[2]);

  double betaAxis2[3];
  betaAxis2[0] = rotMatrix1(0, 0)*_betaAxis[0] + rotMatrix1(0, 1)*_betaAxis[1] + rotMatrix1(0, 2)*_betaAxis[2];
  betaAxis2[1] = rotMatrix1(1, 0)*_betaAxis[0] + rotMatrix1(1, 1)*_betaAxis[1] + rotMatrix1(1, 2)*_betaAxis[2];
  betaAxis2[2] = rotMatrix1(2, 0)*_betaAxis[0] + rotMatrix1(2, 1)*_betaAxis[1] + rotMatrix1(2, 2)*_betaAxis[2];

  irtkMatrix rotMatrix2 = irtkQuaternion::ComputeRotationMatrix(_beta, betaAxis2[0], betaAxis2[1], betaAxis2[2]);

  irtkMatrix secondFirstRot = rotMatrix2*rotMatrix1;

  double gammaAxis2[3];
  gammaAxis2[0] = secondFirstRot(0, 0)*_gammaAxis[0] + secondFirstRot(0, 1)*_gammaAxis[1] + secondFirstRot(0, 2)*_gammaAxis[2];
  gammaAxis2[1] = secondFirstRot(1, 0)*_gammaAxis[0] + secondFirstRot(1, 1)*_gammaAxis[1] + secondFirstRot(1, 2)*_gammaAxis[2];
  gammaAxis2[2] = secondFirstRot(2, 0)*_gammaAxis[0] + secondFirstRot(2, 1)*_gammaAxis[1] + secondFirstRot(2, 2)*_gammaAxis[2];

  irtkMatrix rotMatrix3 = irtkQuaternion::ComputeRotationMatrix(_gamma, gammaAxis2[0], gammaAxis2[1], gammaAxis2[2]);

  // Compute matrix needed to translate axis back to original location with
  // translation.
  irtkMatrix trans2(4, 4);
  trans2.Ident();
  trans2(0, 3) = _origin[0] + _translation[0];
  trans2(1, 3) = _origin[1] + _translation[1];
  trans2(2, 3) = _origin[2] + _translation[2];

  // Compute the final orientation of the alpha, beta, and gamma axes after the
  // rotation.
  _matrix = rotMatrix3*rotMatrix2*rotMatrix1;
  for (int i = 0; i < 3; i++) {
    _finalAlphaAxis[i] = _matrix(i, 0)*_alphaAxis[0] + _matrix(i, 1)*_alphaAxis[1] + _matrix(i, 2)*_alphaAxis[2];
    _finalBetaAxis[i] = _matrix(i, 0)*_betaAxis[0] + _matrix(i, 1)*_betaAxis[1] + _matrix(i, 2)*_betaAxis[2];
    _finalGammaAxis[i] = _matrix(i, 0)*_gammaAxis[0] + _matrix(i, 1)*_gammaAxis[1] + _matrix(i, 2)*_gammaAxis[2];
  }

  // Compute the total matrix of transformations.
  _matrix = trans2*_matrix*scale*trans1;
}

void irtkQuaternionTransformation::ResetOrigin()
{
  _origin[0] = _origin[0] + _translation[0];
  _origin[1] = _origin[1] + _translation[1];
  _origin[2] = _origin[2] + _translation[2];
  _translation[0] = 0;
  _translation[1] = 0;
  _translation[2] = 0;
  UpdateMatrix();
}

void irtkQuaternionTransformation::Transform(double& x, double& y, double& z)
{
  double nx, ny, nz;

  nx = _matrix(0, 0)*x + _matrix(0, 1)*y + _matrix(0, 2)*z + _matrix(0, 3);
  ny = _matrix(1, 0)*x + _matrix(1, 1)*y + _matrix(1, 2)*z + _matrix(1, 3);
  nz = _matrix(2, 0)*x + _matrix(2, 1)*y + _matrix(2, 2)*z + _matrix(2, 3);

  x = nx;
  y = ny;
  z = nz;
}

void irtkQuaternionTransformation::Invert()
{
  irtkMatrix trans2i(4, 4);
  trans2i.Ident();
  trans2i(0, 3) = -_origin[0] - _translation[0];
  trans2i(1, 3) = -_origin[1] - _translation[1];
  trans2i(2, 3) = -_origin[2] - _translation[2];

  irtkMatrix rotMatrix1 = irtkQuaternion::ComputeRotationMatrix(-_gamma, _finalGammaAxis[0], _finalGammaAxis[1], _finalGammaAxis[2]);

  double betaAxis2[3];
  betaAxis2[0] = rotMatrix1(0, 0)*_finalBetaAxis[0] + rotMatrix1(0, 1)*_finalBetaAxis[1] + rotMatrix1(0, 2)*_finalBetaAxis[2];
  betaAxis2[1] = rotMatrix1(1, 0)*_finalBetaAxis[0] + rotMatrix1(1, 1)*_finalBetaAxis[1] + rotMatrix1(1, 2)*_finalBetaAxis[2];
  betaAxis2[2] = rotMatrix1(2, 0)*_finalBetaAxis[0] + rotMatrix1(2, 1)*_finalBetaAxis[1] + rotMatrix1(2, 2)*_finalBetaAxis[2];

  irtkMatrix rotMatrix2 = irtkQuaternion::ComputeRotationMatrix(-_beta, betaAxis2[0], betaAxis2[1], betaAxis2[2]);

  irtkMatrix secondFirstRot = rotMatrix2*rotMatrix1;

  double alphaAxis2[3];
  alphaAxis2[0] = secondFirstRot(0, 0)*_finalAlphaAxis[0] + secondFirstRot(0, 1)*_finalAlphaAxis[1] + secondFirstRot(0, 2)*_finalAlphaAxis[2];
  alphaAxis2[1] = secondFirstRot(1, 0)*_finalAlphaAxis[0] + secondFirstRot(1, 1)*_finalAlphaAxis[1] + secondFirstRot(1, 2)*_finalAlphaAxis[2];
  alphaAxis2[2] = secondFirstRot(2, 0)*_finalAlphaAxis[0] + secondFirstRot(2, 1)*_finalAlphaAxis[1] + secondFirstRot(2, 2)*_finalAlphaAxis[2];

  irtkMatrix rotMatrix3 = irtkQuaternion::ComputeRotationMatrix(-_alpha, alphaAxis2[0], alphaAxis2[1], alphaAxis2[2]);

  irtkMatrix trans1i(4, 4);
  trans1i.Ident();
  trans1i(0, 3) = _origin[0];
  trans1i(1, 3) = _origin[1];
  trans1i(2, 3) = _origin[2];

  for (int i = 0; i < 3; i++) {
    double swap;

    swap = _alphaAxis[i];
    _alphaAxis[i] = _finalAlphaAxis[i];
    _finalAlphaAxis[i] = swap;

    swap = _betaAxis[i];
    _betaAxis[i] = _finalBetaAxis[i];
    _finalBetaAxis[i] = swap;

    swap = _gammaAxis[i];
    _gammaAxis[i] = _finalGammaAxis[i];
    _finalGammaAxis[i] = swap;
  }

  _alpha = -_alpha;
  _beta = -_beta;
  _gamma = -_gamma;
  for (int i = 0; i < 3; i++) {
    double swap;

    swap = _translation[i];
    _translation[i] = -_translation[i];
    _origin[i] = _origin[i] + swap;
  }

  _matrix = trans1i*rotMatrix3*rotMatrix2*rotMatrix1*trans2i;
}

void irtkQuaternionTransformation::Jacobian(double, double, double , irtkMatrix& jac)
{
  jac.Initialize(3, 3);

  for (int i = 0; i < 3; i++)
    for (int j = 0; j < 3; j++)
      jac(i, j) = _matrix(i, j);
}

double irtkQuaternionTransformation::Jacobian(double, double, double)
{
  return (_matrix(0, 0)*(_matrix(1, 1)*_matrix(2, 2) - _matrix(1, 2)*_matrix(2, 1)) - _matrix(0, 1)*(_matrix(1, 0)*_matrix(2, 2) - _matrix(1, 2)*_matrix(2, 0)) + _matrix(0, 2)*(_matrix(1, 0)*_matrix(2, 1) - _matrix(1, 1)*_matrix(2, 0)));
}

void irtkQuaternionTransformation::LocalJacobian(double, double, double, irtkMatrix&)
{
  cerr << "irtkQuaternionTransformation::LocalJacobian: Not implemented yet" << endl;
  exit(1);
}

void irtkQuaternionTransformation::GlobalJacobian(double, double, double, irtkMatrix&)
{
  cerr << "irtkQuaternionTransformation::GlobalJacobian: Not implemented yet" << endl;
  exit(1);
}

bool irtkQuaternionTransformation::IsIdentity()
{
  for (int i = 0; i < 3; i++)
    if (_translation[i] != 0)
      return false;

  return (_alpha == 0 && _beta == 0 && _gamma == 0);
}

int irtkQuaternionTransformation::CheckHeader(char* pFileName)
{
  gzFile inputFile = gzopen(pFileName, "rb");

  if (inputFile == NULL) {
    std::stringstream mesg;

    mesg << "Couldn't open file " << pFileName << " for reading.";

    throw irtkException(mesg.str(), __FILE__, __LINE__);
  }

  // Read the name of the class.
  char buffer[MAX_BUFFER_LENGTH];
  gzgets(inputFile, buffer, strlen("irtkQuaternionTransformation") + 1);

  if (strcmp(buffer, "irtkQuaternionTransformation") != 0) {
    gzclose(inputFile);
    return false;
  }

  gzclose(inputFile);
  return true;
}

irtkCifstream& irtkQuaternionTransformation::Read(irtkCifstream& is)
{
  unsigned int magicNo, transType;

  // Read magic number for transformations.
  is.ReadAsUInt(&magicNo, 1);
  if (magicNo != IRTKTRANSFORMATION_MAGIC) {
    std::cerr << "itQuaternionTransformation::Read: Not a valid transformation file." << std::endl;
    exit(EXIT_FAILURE);
  }

  // Read transformation type.
  is.ReadAsUInt(&transType, 1);
  if (transType != IRTKTRANSFORMATION_QUATERNION) {
    std::cerr << "irtkQuaternionTransformation::Read: Not a valid quaternion transformation." << std::endl;
    exit(EXIT_FAILURE);
  }

  if (_status != NULL)
    delete[] _status;

  _status = new _Status[NumberOfDOFs()];

  // Read the parameters of the transformation.
  is.ReadAsDouble(_origin, 3);
  is.ReadAsDouble(_translation, 3);
  is.ReadAsDouble(_alphaAxis, 3);
  is.ReadAsDouble(_betaAxis, 3);
  is.ReadAsDouble(&_alpha, 1);
  is.ReadAsDouble(&_beta, 1);
  is.ReadAsDouble(&_gamma, 1);
  is.ReadAsInt((int*) _status, NumberOfDOFs());

  UpdateAxes();
  UpdateMatrix();

  return is;
}

istream& irtkQuaternionTransformation::Import(istream& is)
{
  char buffer[1024];
  is >> buffer;
  if (strcmp(buffer, "irtkQuaternionTransformation") != 0) {
    std::cerr << "irtkQuaternionTransformation::Import: Not a valid quaternion transformation file." << std::endl;
    exit(EXIT_FAILURE);
  }

  if (_status != NULL)
    delete[] _status;

  _status = new _Status[NumberOfDOFs()];

  // Read the parameters of the transformation.
  is >> _origin[0] >> _origin[1] >> _origin[2];
  is >> _translation[0] >> _translation[1] >> _translation[2];
  is >> _alphaAxis[0] >> _alphaAxis[1] >> _alphaAxis[2];
  is >> _betaAxis[0] >> _betaAxis[1] >> _betaAxis[2];
  is >> _alpha >> _beta >> _gamma;
  // Read the status values.
  for (int i = 0; i < NumberOfDOFs(); i++) {
    int val;
    is >> val;
    if (val == 0) {
      _status[i] = _Passive;
    } else {
      _status[i] = _Active;
    }
  }

  UpdateAxes();
  UpdateMatrix();

  return is;
}

irtkCofstream& irtkQuaternionTransformation::Write(irtkCofstream& os)
{
  unsigned int magicNumber, transType;

  // Write the magic number for transformations.
  magicNumber = IRTKTRANSFORMATION_MAGIC;
  os.WriteAsUInt(&magicNumber, 1);

  // Write the transformation type.
  transType = IRTKTRANSFORMATION_QUATERNION;
  os.WriteAsUInt(&transType, 1);

  // Write the parameters of the transformation.
  os.WriteAsDouble(_origin, 3);
  os.WriteAsDouble(_translation, 3);
  os.WriteAsDouble(_alphaAxis, 3);
  os.WriteAsDouble(_betaAxis, 3);
  os.WriteAsDouble(&_alpha, 1);
  os.WriteAsDouble(&_beta, 1);
  os.WriteAsDouble(&_gamma, 1);
  os.WriteAsInt((int*) _status, NumberOfDOFs());

  return os;
}

ostream& irtkQuaternionTransformation::Export(ostream& os)
{
  os << "irtkQuaternionTransformation" << std::endl;
  os << _origin[0] << " " << _origin[1] << " " << _origin[2] << std::endl;
  os << _alphaAxis[0] << " " << _alphaAxis[1] << " " << _alphaAxis[2] << std::endl;
  os << _betaAxis[0] << " " << _betaAxis[1] << " " << _betaAxis[2] << std::endl;
  os << _alpha << " " << _beta << " " << _gamma << std::endl;
  // Write the status values.
  for (int i = 0; i < NumberOfDOFs(); i++) {
    int val;
    if (_status[i] == _Passive) {
      val = 0;
    } else {
      val = 1;
    }
    os << val << " ";
  }
  os << std::endl;

  return os;
}

void irtkQuaternionTransformation::Print()
{
  std::cout << "_origin = (" << _origin[0] << ", " << _origin[1] << ", " << _origin[2] << ")" << std::endl
            << "_translation = (" << _translation[0] << ", " << _translation[1] << ", " << _translation[2] << ")" << std::endl
            << "_alphaAxis = (" << _alphaAxis[0] << ", " << _alphaAxis[1] << ", " << _alphaAxis[2] << ")" << std::endl
            << "_betaAxis = (" << _betaAxis[0] << ", " << _betaAxis[1] << ", " << _betaAxis[2] << ")" << std::endl
            << "_gammaAxis = (" << _gammaAxis[0] << ", " << _gammaAxis[1] << ", " << _gammaAxis[2] << ")" << std::endl
            << "_alpha = " << _alpha*180/M_PI << std::endl
            << "_beta = " << _beta*180/M_PI << std::endl
            << "_gamma = " << _gamma*180/M_PI << std::endl;
}
