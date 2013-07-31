/*=========================================================================

  Library   : Image Registration Toolkit (IRTK)
  Module    : $Id: irtkImageAffineRigidTransformation.cc 235 2010-10-18 09:25:20Z dr $
  Copyright : Imperial College, Department of Computing
              Visual Information Processing (VIP), 2008 onwards
  Date      : $Date: 2010-10-18 10:25:20 +0100 (Mon, 18 Oct 2010) $
  Version   : $Revision: 235 $
  Changes   : $Author: dr $

=========================================================================*/

#include <irtkTransformation.h>

#include <irtkImageAffineRigidTransformation.h>

#include <irtkHomogeneousTransformationIterator.h>

irtkImageAffineRigidTransformation::irtkImageAffineRigidTransformation() : irtkImageTransformation()
{
  _rigid = irtkMatrix(4, 4);
  _affine = irtkMatrix(4, 4);
}

irtkImageAffineRigidTransformation::~irtkImageAffineRigidTransformation()
{}

void irtkImageAffineRigidTransformation::SetTransformation(irtkTransformation *transformation)
{
  if (strcmp(transformation->NameOfClass(), "irtkRigidTransformation") == 0) {
    this->_transformation = (irtkRigidTransformation *)transformation;
    return;
  }
  if (strcmp(transformation->NameOfClass(), "irtkAffineTransformation") == 0) {
    this->_transformation = (irtkAffineTransformation *)transformation;
    return;
  }
  cerr << "irtkImageAffineRigidTransformation::SetTransform: Not a ";
  cerr << "rigid or affine transformation: " << transformation->NameOfClass()
       << endl;
  exit(1);
}

void irtkImageAffineRigidTransformation::Run()
{
  int i, j, k, l;

  // Check inputs and outputs
  if (this->_input == NULL) {
    cerr << "irtkImageAffineRigidTransformation::Run: Filter has no input" << endl;
    exit(1);
  }

  if (this->_output == NULL) {
    cerr << "irtkImageAffineRigidTransformation::Run: Filter has no output" << endl;
    exit(1);
  }

  if (this->_transformation == NULL) {
    cerr << "irtkImageAffineRigidTransformation::Run: Filter has no transformation" << endl;
    exit(1);
  }

  if (this->_interpolator == NULL) {
    cerr << "irtkImageAffineRigidTransformation::Run: Filter has no interpolator" << endl;
    exit(1);
  }

  if (this->_input->IsEmpty() == true) {
    cerr << "irtkImageAffineRigidTransformation::Run: Input is empty" << endl;
    this->_input->Print();
    exit(1);
  }

  if (this->_input == this->_output) {
    cerr << "irtkImageAffineRigidTransformation::Run: Input equals output" << endl;
    exit(1);
  }

  // Split image transformation matrix into rigid and affine part
  this->SplitTransformationMatrix();

  // Setup interpolation
  this->_interpolator->SetInput(this->_input);
  this->_interpolator->Initialize();

  // Check for inverse transformation
  if (this->_Invert == true) {
	((irtkAffineTransformation *)this->_transformation)->Invert();
	_rigid.Invert();

    // Create iterator
    irtkHomogeneousTransformationIterator
    iterator((irtkAffineTransformation *)this->_transformation);

    // Loop over all voxels in the output (reference) volume
    for (l = 0; l < this->_output->GetT(); l++) {
      int t = round(this->_input->TimeToImage(this->_output->ImageToTime(l)));
      if ((t >= 0) && (t < this->_input->GetT())) {

    	// Initialize iterator
        iterator.Initialize(this->_output, this->_input);
        for (k = 0; k < this->_output->GetZ(); k++) {
          for (j = 0; j < this->_output->GetY(); j++) {
            for (i = 0; i < this->_output->GetX(); i++) {
              if (this->_output->GetAsDouble(i, j, k, l) > this->_TargetPaddingValue) {
                // Check whether transformed point is inside input volume
                if ((iterator._x > -0.5) && (iterator._x < this->_input->GetX()-0.5) &&
                    (iterator._y > -0.5) && (iterator._y < this->_input->GetY()-0.5) &&
                    (iterator._z > -0.5) && (iterator._z < this->_input->GetZ()-0.5)) {
              	  this->_output->PutAsDouble(i, j, k, l, _ScaleFactor * this->_interpolator->Evaluate(iterator._x, iterator._y, iterator._z, t) + _Offset);
                } else {
              	  this->_output->PutAsDouble(i, j, k, l, this->_SourcePaddingValue);
                }
              } else {
            	this->_output->PutAsDouble(i, j, k, l, this->_SourcePaddingValue);
              }
              iterator.NextX();
            }
            iterator.NextY();
          }
          iterator.NextZ();
        }
      } else {
        for (k = 0; k < this->_output->GetZ(); k++) {
          for (j = 0; j < this->_output->GetY(); j++) {
            for (i = 0; i < this->_output->GetX(); i++) {
          	  this->_output->PutAsDouble(i, j, k, l, this->_SourcePaddingValue);
            }
          }
        }
      }
    }

    // Transforms the image origin and axis, not using the affine parameters
    this->TransformImageCoordSystem();

    // Reverse invert matrix
	((irtkAffineTransformation *)this->_transformation)->Invert();
	_rigid.Invert();

  } else {
	// Transforms the image origin and axis, not using the affine parameters
	this->TransformImageCoordSystem();

    // Create iterator
    irtkHomogeneousTransformationIterator
    iterator((irtkAffineTransformation *)this->_transformation);

    // Loop over all voxels in the output (reference) volume
    for (l = 0; l < this->_output->GetT(); l++) {
      int t = round(this->_input->TimeToImage(this->_output->ImageToTime(l)));
      if ((t >= 0) && (t < this->_input->GetT())) {

    	// Initialize iterator
        iterator.Initialize(this->_output, this->_input);
        for (k = 0; k < this->_output->GetZ(); k++) {
          for (j = 0; j < this->_output->GetY(); j++) {
            for (i = 0; i < this->_output->GetX(); i++) {
              if (this->_output->GetAsDouble(i, j, k, l) > this->_TargetPaddingValue) {
                // Check whether transformed point is inside input volume
                if ((iterator._x > -0.5) && (iterator._x < this->_input->GetX()-0.5) &&
                    (iterator._y > -0.5) && (iterator._y < this->_input->GetY()-0.5) &&
                    (iterator._z > -0.5) && (iterator._z < this->_input->GetZ()-0.5)) {
              	  this->_output->PutAsDouble(i, j, k, l, _ScaleFactor * this->_interpolator->Evaluate(iterator._x, iterator._y, iterator._z, t) + _Offset);
                } else {
              	  this->_output->PutAsDouble(i, j, k, l, this->_SourcePaddingValue);
                }
              } else {
            	this->_output->PutAsDouble(i, j, k, l, this->_SourcePaddingValue);
              }
              iterator.NextX();
            }
            iterator.NextY();
          }
          iterator.NextZ();
        }
      } else {
        for (k = 0; k < this->_output->GetZ(); k++) {
          for (j = 0; j < this->_output->GetY(); j++) {
            for (i = 0; i < this->_output->GetX(); i++) {
          	  this->_output->PutAsDouble(i, j, k, l, this->_SourcePaddingValue);
            }
          }
        }
      }
    }
  }

  // Merge transformation matrizes back together
  MergeTranformationMatrix();

}


void irtkImageAffineRigidTransformation::SplitTransformationMatrix()
{
  double rx, ry, rz, cosrx, cosry, cosrz, sinrx, sinry, sinrz, tx, ty, tz;
  double sxy, sxz, syz, tansxy, tansxz, tansyz, sx, sy, sz;

  if (strcmp(_transformation->NameOfClass(), "irtkRigidTransformation") == 0) {
	_rigid = ((irtkRigidTransformation *)this->_transformation)->GetMatrix();
	_affine.Ident();
  } else {
	// build rigid matrix
    _rigid.Ident();

    // Get parameters
    rx = ((irtkAffineTransformation *)this->_transformation)->GetRotationX();
    ry = ((irtkAffineTransformation *)this->_transformation)->GetRotationY();
    rz = ((irtkAffineTransformation *)this->_transformation)->GetRotationZ();
    tx = ((irtkAffineTransformation *)this->_transformation)->GetTranslationX();
    ty = ((irtkAffineTransformation *)this->_transformation)->GetTranslationY();
    tz = ((irtkAffineTransformation *)this->_transformation)->GetTranslationZ();

    // Update sines and cosines
    cosrx = cos(rx*(M_PI/180.0));
    cosry = cos(ry*(M_PI/180.0));
    cosrz = cos(rz*(M_PI/180.0));
    sinrx = sin(rx*(M_PI/180.0));
    sinry = sin(ry*(M_PI/180.0));
    sinrz = sin(rz*(M_PI/180.0));

    // Add other transformation parameters to transformation matrix
    _rigid(0,0) = cosry*cosrz;
    _rigid(0,1) = cosry*sinrz;
    _rigid(0,2) = -sinry;
    _rigid(0,3) = tx;
    _rigid(1,0) = (sinrx*sinry*cosrz-cosrx*sinrz);
    _rigid(1,1) = (sinrx*sinry*sinrz+cosrx*cosrz);
    _rigid(1,2) = sinrx*cosry;
    _rigid(1,3) = ty;
    _rigid(2,0) = (cosrx*sinry*cosrz+sinrx*sinrz);
    _rigid(2,1) = (cosrx*sinry*sinrz-sinrx*cosrz);
    _rigid(2,2) = cosrx*cosry;
    _rigid(2,3) = tz;
    _rigid(3,3) = 1.0;


    // build affine matrix
    _affine.Ident();

    // Get parameters
    sxy = ((irtkAffineTransformation *)this->_transformation)->GetShearXY();
    sxz = ((irtkAffineTransformation *)this->_transformation)->GetShearXZ();
    syz = ((irtkAffineTransformation *)this->_transformation)->GetShearYZ();
    sx = ((irtkAffineTransformation *)this->_transformation)->GetScaleX();
    sy = ((irtkAffineTransformation *)this->_transformation)->GetScaleY();
    sz = ((irtkAffineTransformation *)this->_transformation)->GetScaleZ();

    // Update affine transformation: Add shearing
    irtkMatrix skew(4, 4);
    skew.Ident();
    tansxy = tan(sxy*(M_PI/180.0));
    tansxz = tan(sxz*(M_PI/180.0));
    tansyz = tan(syz*(M_PI/180.0));
    skew(0, 1) = tansxy;
    skew(0, 2) = tansxz;
    skew(1, 2) = tansyz;
    _affine *= skew;

    // Update affine transformation: Add scaling
    irtkMatrix scale(4, 4);
    scale.Ident();
    scale(0, 0) = sx / 100.0;
    scale(1, 1) = sy / 100.0;
    scale(2, 2) = sz / 100.0;
    _affine *= scale;
  }

  ((irtkAffineTransformation *)this->_transformation)->PutMatrix(_affine);
}

void irtkImageAffineRigidTransformation::MergeTranformationMatrix()
{
  if (strcmp(_transformation->NameOfClass(), "irtkRigidTransformation") == 0) {
	((irtkRigidTransformation *)this->_transformation)->PutMatrix(_rigid);
  } else {
	_affine = _affine * _rigid;
    ((irtkAffineTransformation *)this->_transformation)->PutMatrix(_affine);
  }
}

void irtkImageAffineRigidTransformation::TransformImageCoordSystem()
{
  double xorigin, yorigin, zorigin, * xaxis, * yaxis, * zaxis;
  irtkMatrix axis(3, 3), rot(3, 3);
  int i, j;

  xaxis = new double [3];
  yaxis = new double [3];
  zaxis = new double [3];
  this->_output->GetOrientation(xaxis, yaxis, zaxis);
  this->_output->GetOrigin(xorigin, yorigin, zorigin);

  xorigin += _rigid(0,3);
  yorigin += _rigid(1,3);
  zorigin += _rigid(2,3);

  axis.Ident();
  for (i = 0; i < 3; i++) {
	axis(i, 0) = xaxis[i];
	axis(i, 1) = yaxis[i];
	axis(i, 2) = zaxis[i];

	for (j = 0; j < 3; j++) {
	  rot(i, j) = _rigid(i, j);
	}
  }
//  cout<<"axis_before: "<<endl;
//  cout<<axis(0, 0)<<" , "<<axis(0, 1)<<" , "<<axis(0, 2)<<endl;
//  cout<<axis(1, 0)<<" , "<<axis(1, 1)<<" , "<<axis(1, 2)<<endl;
//  cout<<axis(2, 0)<<" , "<<axis(2, 1)<<" , "<<axis(2, 2)<<endl;
//  cout<<"rot: "<<endl;
//  cout<<rot(0, 0)<<" , "<<rot(0, 1)<<" , "<<rot(0, 2)<<endl;
//  cout<<rot(1, 0)<<" , "<<rot(1, 1)<<" , "<<rot(1, 2)<<endl;
//  cout<<rot(2, 0)<<" , "<<rot(2, 1)<<" , "<<rot(2, 2)<<endl;
//  axis = rot*axis;
//  cout<<"axis_after: "<<endl;
//  cout<<axis(0, 0)<<" , "<<axis(0, 1)<<" , "<<axis(0, 2)<<endl;
//  cout<<axis(1, 0)<<" , "<<axis(1, 1)<<" , "<<axis(1, 2)<<endl;
//  cout<<axis(2, 0)<<" , "<<axis(2, 1)<<" , "<<axis(2, 2)<<endl;
  for (i = 0; i < 3; i++) {
	xaxis[i] = axis(i, 0);
	yaxis[i] = axis(i, 1);
	zaxis[i] = axis(i, 2);
  }

  this->_output->PutOrigin(xorigin, yorigin, zorigin);
  this->_output->PutOrientation(xaxis, yaxis, zaxis);
}












