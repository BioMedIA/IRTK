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

irtkMatrix irtkQuaternion::ComputeRotationMatrix() const
{
  irtkMatrix rot(4, 4);
  rot(0, 0) = 1 - 2*_w*_w - 2*_v*_v;
  rot(0, 1) = 2*_u*_v - 2*_t*_w;
  rot(0, 2) = 2*_u*_w + 2*_t*_v;
  rot(0, 3) = 0;
  rot(1, 0) = 2*_u*_v + 2*_t*_w;
  rot(1, 1) = 1 - 2*_u*_u - 2*_w*_w;
  rot(1, 2) = 2*_v*_w - 2*_t*_u;
  rot(1, 3) = 0;
  rot(2, 0) = 2*_u*_w - 2*_t*_v;
  rot(2, 1) = 2*_v*_w + 2*_t*_u;
  rot(2, 2) = 1 - 2*_u*_u - 2*_v*_v;
  rot(2, 3) = 0;
  rot(3, 0) = 0;
  rot(3, 1) = 0;
  rot(3, 2) = 0;
  rot(3, 3) = 1;

  return rot;
}

irtkMatrix irtkQuaternion::ComputeTransformationMatrix(double alpha,
    double u, double v, double w, double ox, double oy, double oz)
{
  irtkMatrix translateToOrigin(4, 4);
  translateToOrigin.Ident();
  translateToOrigin(0, 3) = -ox;
  translateToOrigin(1, 3) = -oy;
  translateToOrigin(2, 3) = -oz;

  irtkMatrix rot = irtkQuaternion::ComputeRotationMatrix(alpha, u, v, w);

  irtkMatrix translateToPoint(4, 4);
  translateToPoint.Ident();
  translateToPoint(0, 3) = ox;
  translateToPoint(1, 3) = oy;
  translateToPoint(2, 3) = oz;

  return translateToPoint*rot*translateToOrigin;
}

irtkMatrix irtkQuaternion::ComputeRotationMatrix(double alpha, double u, double v, double w)
{
  irtkQuaternion q = irtkQuaternion::ComputeRotationQuaternion(alpha, u, v, w);

  irtkMatrix rot(4, 4);
  rot(0, 0) = 1 - 2*q._w*q._w - 2*q._v*q._v;
  rot(0, 1) = 2*q._u*q._v - 2*q._t*q._w;
  rot(0, 2) = 2*q._u*q._w + 2*q._t*q._v;
  rot(0, 3) = 0;
  rot(1, 0) = 2*q._u*q._v + 2*q._t*q._w;
  rot(1, 1) = 1 - 2*q._u*q._u - 2*q._w*q._w;
  rot(1, 2) = 2*q._v*q._w - 2*q._t*q._u;
  rot(1, 3) = 0;
  rot(2, 0) = 2*q._u*q._w - 2*q._t*q._v;
  rot(2, 1) = 2*q._v*q._w + 2*q._t*q._u;
  rot(2, 2) = 1 - 2*q._u*q._u - 2*q._v*q._v;
  rot(2, 3) = 0;
  rot(3, 0) = 0;
  rot(3, 1) = 0;
  rot(3, 2) = 0;
  rot(3, 3) = 1;

  return rot;
}

irtkQuaternion irtkQuaternion::ComputeRotationQuaternion(double alpha, double u, double v, double w)
{
  double length = sqrt(u*u + v*v + w*w);

  double t2 = cos(alpha/2);
  double u2 = (u*sin(alpha/2))/length;
  double v2 = (v*sin(alpha/2))/length;
  double w2 = (w*sin(alpha/2))/length;

  return irtkQuaternion(t2, u2, v2, w2);
}
