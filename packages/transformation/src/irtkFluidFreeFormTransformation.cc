/*=========================================================================
 
  Library   : Image Registration Toolkit (IRTK)
  Module    : $Id$
  Copyright : Imperial College, Department of Computing
              Visual Information Processing (VIP), 2008 onwards
  Date      : $Date$
  Version   : $Revision$
  Changes   : $Author$

Copyright (c) IXICO LIMITED
All rights reserved.
See COPYRIGHT for details
 
=========================================================================*/

#include <irtkTransformation.h>

#include <irtkFluidFreeFormTransformation.h>

#include <newt2.h>

irtkFluidFreeFormTransformation::irtkFluidFreeFormTransformation() : irtkMultiLevelFreeFormTransformation()
{}

irtkFluidFreeFormTransformation::irtkFluidFreeFormTransformation(const irtkFluidFreeFormTransformation &transformation) : irtkMultiLevelFreeFormTransformation(transformation)
{
}

irtkFluidFreeFormTransformation::irtkFluidFreeFormTransformation(const irtkRigidTransformation &transformation) : irtkMultiLevelFreeFormTransformation(transformation)
{
  int i;

  // Initialize local transformation
  for (i = 0; i < MAX_TRANS; i++) {
    _localTransformation[i] = NULL;
  }

  _NumberOfLevels = 0;
}

irtkFluidFreeFormTransformation::irtkFluidFreeFormTransformation(const irtkAffineTransformation &transformation) : irtkMultiLevelFreeFormTransformation(transformation)
{
  int i;

  // Initialize local transformation
  for (i = 0; i < MAX_TRANS; i++) {
    _localTransformation[i] = NULL;
  }

  _NumberOfLevels = 0;
}

void irtkFluidFreeFormTransformation::CombineLocalTransformation()
{
  cerr << "irtkFluidFreeFormTransformation::CombineLocalTransformation: Does not make sense" << endl;
  exit(1);
}

void irtkFluidFreeFormTransformation::Transform(double &x, double &y, double &z, double t)
{
  int i;

  // Compute global transformation
  this->irtkAffineTransformation::Transform(x, y, z, t);

  // Compute local transformation
  for (i = 0; i < _NumberOfLevels; i++) {

    // Displacement
    _localTransformation[i]->Transform(x, y, z, t);
  }
}

void irtkFluidFreeFormTransformation::Displacement(double& x, double& y, double& z, double t)
{
  double u, v, w;

  // Save current position
  u = x;
  v = y;
  w = z;

  // Transform current position
  this->Transform(u, v, w, t);

  x = u - x;
  y = v - y;
  z = w - z;
}

void irtkFluidFreeFormTransformation::Transform(int n, double &x, double &y, double &z, double t)
{
  int i;

  if (n > _NumberOfLevels) {
    cerr << "irtkFluidFreeFormTransformation::Transform: No such "
    << "transformation" << endl;
    exit(1);
  }

  // Compute global transformation
  this->irtkAffineTransformation::Transform(x, y, z, t);

  // Compute local transformation
  for (i = 0; i < n; i++) {

    // Displacement
    _localTransformation[i]->Transform(x, y, z, t);
  }
}

void irtkFluidFreeFormTransformation::GlobalTransform(double &, double &, double &, double)
{
  cerr << "irtkFluidFreeFormTransformation::GlobalTransform: Does not make sense" << endl;
  exit(1);
}

void irtkFluidFreeFormTransformation::GlobalDisplacement(double &, double &, double &, double)
{
  cerr << "irtkFluidFreeFormTransformation::GlobalDisplacement: Does not make sense" << endl;
  exit(1);
}

void irtkFluidFreeFormTransformation::LocalTransform(double &, double &, double &, double)
{
  cerr << "irtkFluidFreeFormTransformation::LocalTransform: Does not make sense" << endl;
  exit(1);
}

void irtkFluidFreeFormTransformation::LocalTransform(int, double &, double &, double &, double)
{
  cerr << "irtkFluidFreeFormTransformation::LocalTransform: Does not make sense" << endl;
  exit(1);
}

void irtkFluidFreeFormTransformation::LocalDisplacement(int, double &, double &, double &, double)
{
  cerr << "irtkFluidFreeFormTransformation::LocalDisplacement: Does not make sense" << endl;
  exit(1);
}

void irtkFluidFreeFormTransformation::LocalDisplacement(double &, double &, double &, double)
{
  cerr << "irtkFluidFreeFormTransformation::LocalDisplacement: Does not make sense" << endl;
  exit(1);
}

void irtkFluidFreeFormTransformation::MergeGlobalIntoLocalDisplacement()
{
	int i;

  irtkBSplineFreeFormTransformation3D *ffd =
  		dynamic_cast<irtkBSplineFreeFormTransformation3D *> (this->GetLocalTransformation(0));

  if (ffd == NULL){
  	cerr << "irtkFluidFreeFormTransformation::MergeGlobalIntoLocalDisplacement: ";
  	cerr << "FFD should be of type irtkBSplineFreeFormTransformation3D" << endl;
  	exit(1);
  }

  // Get a copy for making a FFD interpolation of the global affine component.
  irtkBSplineFreeFormTransformation3D *ffdCopy = new irtkBSplineFreeFormTransformation3D(*ffd);

  this->InterpolateGlobalDisplacement(ffdCopy);

  // Shift all existing FFDs in the stack along by one place.
  for (i = this->NumberOfLevels(); i > 0; i--){
  	this->PutLocalTransformation(this->GetLocalTransformation(i-1), i);
  }

  // Insert the new FFD for the global transformation into the first position.
  this->PutLocalTransformation(ffdCopy, 0);

  // Reset matrix previously used for global transform to identity
  this->Reset();
}

void irtkFluidFreeFormTransformation::Jacobian(irtkMatrix &jac, double x, double y, double z, double t)
{
  int i;
  irtkMatrix tmp_jac;

  // Initialize to identity
  jac.Ident();

  // Compute global jacobian
  this->GlobalJacobian(tmp_jac, x, y, z, t);
  jac = tmp_jac * jac;

  // Transform
  this->irtkAffineTransformation::Transform(x, y, z, t);

  // Compute local jacobian
  for (i = 0; i < _NumberOfLevels; i++) {

    // Calculate jacobian
    _localTransformation[i]->Jacobian(tmp_jac, x, y, z, t);

    // Multiply jacobian
    jac = tmp_jac * jac;

    // Transform
    _localTransformation[i]->Transform(x, y, z, t);

  }
}

void irtkFluidFreeFormTransformation::LocalJacobian(irtkMatrix &jac, double x, double y, double z, double t)
{
  cerr << "irtkFluidFreeFormTransformation::LocalJacobian: Does not make sense" << endl;
  exit(1);
}

irtkCifstream& irtkFluidFreeFormTransformation::Read(irtkCifstream& from)
{
  int i, offset;
  unsigned int magic_no, trans_type;

  // Delete old local transformations
  for (i = 0; i < _NumberOfLevels; i++) {
    if (_localTransformation[i] != NULL) delete _localTransformation[i];
  }

  // Read magic no. for transformations
  from.ReadAsUInt(&magic_no, 1);
  if (magic_no != IRTKTRANSFORMATION_MAGIC) {
    cerr << "irtkFluidFreeFormTransformation::Read: Not a vaild transformation file" << endl;
    exit(1);
  }

  // Read transformation type
  from.ReadAsUInt(&trans_type, 1);
  if (trans_type != IRTKTRANSFORMATION_FLUID) {
    cerr << "irtkFluidFreeFormTransformation::Read: Not a vaild fluid transformation" << endl;
    exit(1);
  }

  // Read number of local transformations
  from.ReadAsInt(&_NumberOfLevels, 1);

  // Read global transformation
  this->irtkAffineTransformation::Read(from);

  // Read local transformations
  for (i = 0; i < _NumberOfLevels; i++) {

    // Remember current file location
    offset = from.Tell();

    // Read magic no. for transformations
    from.ReadAsUInt(&magic_no, 1);
    if (magic_no != IRTKTRANSFORMATION_MAGIC) {
      cerr << "irtkFluidFreeFormTransformation::Read: Not a valid fluid transformation found at = " << offset << endl;
      exit(1);
    }

    // Read transformation type
    from.ReadAsUInt(&trans_type, 1);
    if ((trans_type == IRTKTRANSFORMATION_LINEAR_FFD) || (trans_type == IRTKTRANSFORMATION_LINEAR_FFD_EXT1)) {
      from.Seek(offset);
      _localTransformation[i] = new irtkLinearFreeFormTransformation;
      ((irtkLinearFreeFormTransformation *)_localTransformation[i])->Read(from);
    } else {
      if ((trans_type == IRTKTRANSFORMATION_BSPLINE_FFD) || (trans_type == IRTKTRANSFORMATION_BSPLINE_FFD_EXT1)) {
        from.Seek(offset);
        _localTransformation[i] = new irtkBSplineFreeFormTransformation;
        ((irtkBSplineFreeFormTransformation *)_localTransformation[i])->Read(from);
      } else {
        cerr << "irtkFluidFreeFormTransformation::Read: No a valid FFD transformation type at = " << offset << " " << endl;
        exit(1);
      }
    }
  }
  return from;
}

irtkCofstream& irtkFluidFreeFormTransformation::Write(irtkCofstream& to)
{
  int i;
  unsigned int magic_no, trans_type;

  // Write magic no. for transformations
  magic_no = IRTKTRANSFORMATION_MAGIC;
  to.WriteAsUInt(&magic_no, 1);

  // Write transformation type
  trans_type = IRTKTRANSFORMATION_FLUID;
  to.WriteAsUInt(&trans_type, 1);

  // Write transformation type
  to.WriteAsInt(&_NumberOfLevels, 1);

  // Write global transformation
  this->irtkAffineTransformation::Write(to);

  // Write local transformations
  for (i = 0; i < _NumberOfLevels; i++) {
    _localTransformation[i]->Write(to);
  }

  return to;
}

double irtkFluidFreeFormTransformation::Inverse(double &x, double &y, double &z, double t, double tolerance)
{
  int check;
  double error;

  // Initialize global variables
  irtkTransformationPointer  = this;
  x_invert = x;
  y_invert = y;
  z_invert = z;

  // Pointer to B-spline wrapper
  void (*Newton_function)(int, float [], float []) = irtkTransformationEvaluate;

  // Calculate initial estimate using affine transformation
  this->irtkHomogeneousTransformation::Inverse(x, y, z, t);

  // Inverse
  float invert[3], f_invert[3];
  invert[0] = x;
  invert[1] = y;
  invert[2] = z;

  // Numerically approximate the inverse transformation
  newt2(invert-1, 3, &check, Newton_function);

  // Calculate error
  irtkTransformationEvaluate(3, invert-1, f_invert-1);
  error = sqrt(f_invert[0]*f_invert[0]+f_invert[1]*f_invert[1]+f_invert[2]*f_invert[2]);
  if (error > tolerance) {
    cout << "irtkFluidFreeFormTransformation::Inverse: RMS error = " << error << "\n";
  }

  // Set output to solution
  x = invert[0];
  y = invert[1];
  z = invert[2];

  return error;
}
