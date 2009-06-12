/*=========================================================================

  Library   : Image Registration Toolkit (IRTK)
  Module    : $Id$
  Copyright : Imperial College, Department of Computing
              Visual Information Processing (VIP), 2008 onwards
  Date      : $Date$
  Version   : $Revision$
  Changes   : $Author$

=========================================================================*/

#include <irtkImage.h>

#include <irtkImageFunction.h>

irtkBSplineInterpolateImageFunction2D::irtkBSplineInterpolateImageFunction2D(int SplineDegree)
{
  if ((SplineDegree < 2) || (SplineDegree > 5)) {
    cerr << "irtkBSplineInterpolateImageFunction2D: Unsupported spline degree\n";
    exit(1);
  }
  _SplineDegree = SplineDegree;
}

irtkBSplineInterpolateImageFunction2D::~irtkBSplineInterpolateImageFunction2D(void)
{}

const char *irtkBSplineInterpolateImageFunction2D::NameOfClass()
{
  return "irtkBSplineInterpolateImageFunction2D";
}

void irtkBSplineInterpolateImageFunction2D::Initialize()
{
  /// Initialize baseclass
  this->irtkImageFunction::Initialize();

  // Check if image is 2D
  if (this->_input->GetZ() != 1) {
    cerr << "irtkBSplineInterpolateImageFunction2D::Initialize(): ";
    cerr << "Input image is not 2D" << endl;
    exit(1);
  }

  // Compute size of image
  this->_x = this->_input->GetX();
  this->_y = this->_input->GetY();
  this->_t = this->_input->GetT();

  // Compute size of image divided by 2
  this->_xhalf = 2 * this->_x - 2;
  this->_yhalf = 2 * this->_y - 2;

  // Compute domain on which the B-spline is defined
  this->_x1 = round(_SplineDegree/2.0);
  this->_y1 = round(_SplineDegree/2.0);
  this->_x2 = this->_input->GetX() - round(_SplineDegree/2.0 + 1);
  this->_y2 = this->_input->GetY() - round(_SplineDegree/2.0 + 1);

  // Compute min and max values
  this->_input->GetMinMaxAsDouble(&this->_min, &this->_max);

  // Allocate coefficient image
  this->_coeff = irtkRealImage(this->_x, this->_y, 1);

  // Compute B-Spline interpolation coefficients
  this->ComputeCoefficients();
}

double irtkBSplineInterpolateImageFunction2D::InitialAntiCausalCoefficient(double c[], int DataLength, double z)
{
  /* this initialization corresponds to mirror boundaries */
  return((z / (z * z - 1.0)) * (z * c[DataLength - 2] + c[DataLength - 1]));
}

double irtkBSplineInterpolateImageFunction2D::InitialCausalCoefficient(double c[], int DataLength, double z, double Tolerance)
{
  double Sum, zn, z2n, iz;
  int n, Horizon;

  /* this initialization corresponds to mirror boundaries */
  Horizon = DataLength;
  if (Tolerance > 0.0) {
    Horizon = (int)ceil(log(Tolerance) / log(fabs(z)));
  }
  if (Horizon < DataLength) {
    /* accelerated loop */
    zn = z;
    Sum = c[0];
    for (n = 1; n < Horizon; n++) {
      Sum += zn * c[n];
      zn *= z;
    }
    return(Sum);
  } else {
    /* full loop */
    zn = z;
    iz = 1.0 / z;
    z2n = pow(z, (double)(DataLength - 1));
    Sum = c[0] + z2n * c[DataLength - 1];
    z2n *= z2n * iz;
    for (n = 1; n <= DataLength - 2; n++) {
      Sum += (zn + z2n) * c[n];
      zn *= z;
      z2n *= iz;
    }
    return(Sum / (1.0 - zn * zn));
  }
}

void irtkBSplineInterpolateImageFunction2D::ConvertToInterpolationCoefficients(double *c, int DataLength, double *z, int NbPoles, double Tolerance)
{
  double Lambda = 1.0;
  int n, k;

  /* special case required by mirror boundaries */
  if (DataLength == 1) {
    return;
  }

  /* compute the overall gain */
  for (k = 0; k < NbPoles; k++) {
    Lambda = Lambda * (1.0 - z[k]) * (1.0 - 1.0 / z[k]);
  }

  /* apply the gain */
  for (n = 0; n < DataLength; n++) {
    c[n] *= Lambda;
  }

  /* loop over all poles */
  for (k = 0; k < NbPoles; k++) {
    /* causal initialization */
    c[0] = InitialCausalCoefficient(c, DataLength, z[k], Tolerance);
    /* causal recursion */
    for (n = 1; n < DataLength; n++) {
      c[n] += z[k] * c[n - 1];
    }
    /* anticausal initialization */
    c[DataLength - 1] = InitialAntiCausalCoefficient(c, DataLength, z[k]);
    /* anticausal recursion */
    for (n = DataLength - 2; 0 <= n; n--) {
      c[n] = z[k] * (c[n + 1] - c[n]);
    }
  }
}

void irtkBSplineInterpolateImageFunction2D::ComputeCoefficients()
{
  double *data;
  double Pole[2];
  int NbPoles;
  int x, y, z, t;

  /* recover the poles from a lookup table */
  switch (_SplineDegree) {
  case 2:
    NbPoles = 1;
    Pole[0] = sqrt(8.0) - 3.0;
    break;
  case 3:
    NbPoles = 1;
    Pole[0] = sqrt(3.0) - 2.0;
    break;
  case 4:
    NbPoles = 2;
    Pole[0] = sqrt(664.0 - sqrt(438976.0)) + sqrt(304.0) - 19.0;
    Pole[1] = sqrt(664.0 + sqrt(438976.0)) - sqrt(304.0) - 19.0;
    break;
  case 5:
    NbPoles = 2;
    Pole[0] = sqrt(135.0 / 2.0 - sqrt(17745.0 / 4.0)) + sqrt(105.0 / 4.0)
              - 13.0 / 2.0;
    Pole[1] = sqrt(135.0 / 2.0 + sqrt(17745.0 / 4.0)) - sqrt(105.0 / 4.0)
              - 13.0 / 2.0;
    break;
  default:
    cerr << "Invalid spline degree" << endl;
    exit(1);
  }

  /* convert the image samples into interpolation coefficients */

  // Set z to fixed value.
  z = 0;
  
  /* in-place separable process, along x */
  data = new double[_x];
  for (t = 0; t < this->_t; t++) {
      for (y = 0; y < this->_y; y++) {
        for (x = 0; x < this->_x; x++) {
          data[x] = this->_input->GetAsDouble(x, y, z, t);
        }
        ConvertToInterpolationCoefficients(data, this->_x, Pole, NbPoles, DBL_EPSILON);
        for (x = 0; x < this->_x; x++) {
          this->_coeff(x, y, z, t) = data[x];
        }
      }
  }
  delete []data;

  /* in-place separable process, along y */
  data = new double[_y];
  for (t = 0; t < this->_t; t++) {
      for (x = 0; x < this->_x; x++) {
        for (y = 0; y < this->_y; y++) {
          data[y] = this->_coeff(x, y, z, t);
        }
        ConvertToInterpolationCoefficients(data, this->_y, Pole, NbPoles, DBL_EPSILON);
        for (y = 0; y < this->_y; y++) {
          this->_coeff(x, y, z, t) = data[y];
        }
      }
  }

  delete []data;

}

double irtkBSplineInterpolateImageFunction2D::Evaluate(double x, double y, double z, double time)
{
  int i, j, k, l, m;
  int xIndex[6], yIndex[6];
  double xWeight[6], yWeight[6];
  double value, w, w2, w4, t, t0, t1;

  /* compute the interpolation indexes */
  if (_SplineDegree & 1) {
    i = (int)floor(x) - this->_SplineDegree / 2;
    j = (int)floor(y) - this->_SplineDegree / 2;
    for (m = 0; m <= this->_SplineDegree; m++) {
      xIndex[m] = i++;
      yIndex[m] = j++;
    }
  } else {
    i = (int)floor(x + 0.5) - this->_SplineDegree / 2;
    j = (int)floor(y + 0.5) - this->_SplineDegree / 2;
    for (m = 0; m <= this->_SplineDegree; m++) {
      xIndex[m] = i++;
      yIndex[m] = j++;
    }
  }

  /* compute the interpolation weights */
  switch (_SplineDegree) {
  case 2:
    /* x */
    w = x - (double)xIndex[1];
    xWeight[1] = 3.0 / 4.0 - w * w;
    xWeight[2] = (1.0 / 2.0) * (w - xWeight[1] + 1.0);
    xWeight[0] = 1.0 - xWeight[1] - xWeight[2];
    /* y */
    w = y - (double)yIndex[1];
    yWeight[1] = 3.0 / 4.0 - w * w;
    yWeight[2] = (1.0 / 2.0) * (w - yWeight[1] + 1.0);
    yWeight[0] = 1.0 - yWeight[1] - yWeight[2];
    break;
  case 3:
    /* x */
    w = x - (double)xIndex[1];
    xWeight[3] = (1.0 / 6.0) * w * w * w;
    xWeight[0] = (1.0 / 6.0) + (1.0 / 2.0) * w * (w - 1.0) - xWeight[3];
    xWeight[2] = w + xWeight[0] - 2.0 * xWeight[3];
    xWeight[1] = 1.0 - xWeight[0] - xWeight[2] - xWeight[3];
    /* y */
    w = y - (double)yIndex[1];
    yWeight[3] = (1.0 / 6.0) * w * w * w;
    yWeight[0] = (1.0 / 6.0) + (1.0 / 2.0) * w * (w - 1.0) - yWeight[3];
    yWeight[2] = w + yWeight[0] - 2.0 * yWeight[3];
    yWeight[1] = 1.0 - yWeight[0] - yWeight[2] - yWeight[3];
    break;
  case 4:
    /* x */
    w = x - (double)xIndex[2];
    w2 = w * w;
    t = (1.0 / 6.0) * w2;
    xWeight[0] = 1.0 / 2.0 - w;
    xWeight[0] *= xWeight[0];
    xWeight[0] *= (1.0 / 24.0) * xWeight[0];
    t0 = w * (t - 11.0 / 24.0);
    t1 = 19.0 / 96.0 + w2 * (1.0 / 4.0 - t);
    xWeight[1] = t1 + t0;
    xWeight[3] = t1 - t0;
    xWeight[4] = xWeight[0] + t0 + (1.0 / 2.0) * w;
    xWeight[2] = 1.0 - xWeight[0] - xWeight[1] - xWeight[3] - xWeight[4];
    /* y */
    w = y - (double)yIndex[2];
    w2 = w * w;
    t = (1.0 / 6.0) * w2;
    yWeight[0] = 1.0 / 2.0 - w;
    yWeight[0] *= yWeight[0];
    yWeight[0] *= (1.0 / 24.0) * yWeight[0];
    t0 = w * (t - 11.0 / 24.0);
    t1 = 19.0 / 96.0 + w2 * (1.0 / 4.0 - t);
    yWeight[1] = t1 + t0;
    yWeight[3] = t1 - t0;
    yWeight[4] = yWeight[0] + t0 + (1.0 / 2.0) * w;
    yWeight[2] = 1.0 - yWeight[0] - yWeight[1] - yWeight[3] - yWeight[4];
    break;
  case 5:
    /* x */
    w = x - (double)xIndex[2];
    w2 = w * w;
    xWeight[5] = (1.0 / 120.0) * w * w2 * w2;
    w2 -= w;
    w4 = w2 * w2;
    w -= 1.0 / 2.0;
    t = w2 * (w2 - 3.0);
    xWeight[0] = (1.0 / 24.0) * (1.0 / 5.0 + w2 + w4) - xWeight[5];
    t0 = (1.0 / 24.0) * (w2 * (w2 - 5.0) + 46.0 / 5.0);
    t1 = (-1.0 / 12.0) * w * (t + 4.0);
    xWeight[2] = t0 + t1;
    xWeight[3] = t0 - t1;
    t0 = (1.0 / 16.0) * (9.0 / 5.0 - t);
    t1 = (1.0 / 24.0) * w * (w4 - w2 - 5.0);
    xWeight[1] = t0 + t1;
    xWeight[4] = t0 - t1;
    /* y */
    w = y - (double)yIndex[2];
    w2 = w * w;
    yWeight[5] = (1.0 / 120.0) * w * w2 * w2;
    w2 -= w;
    w4 = w2 * w2;
    w -= 1.0 / 2.0;
    t = w2 * (w2 - 3.0);
    yWeight[0] = (1.0 / 24.0) * (1.0 / 5.0 + w2 + w4) - yWeight[5];
    t0 = (1.0 / 24.0) * (w2 * (w2 - 5.0) + 46.0 / 5.0);
    t1 = (-1.0 / 12.0) * w * (t + 4.0);
    yWeight[2] = t0 + t1;
    yWeight[3] = t0 - t1;
    t0 = (1.0 / 16.0) * (9.0 / 5.0 - t);
    t1 = (1.0 / 24.0) * w * (w4 - w2 - 5.0);
    yWeight[1] = t0 + t1;
    yWeight[4] = t0 - t1;
    break;
  default:
    printf("Invalid spline degree\n");
    return(0.0);
  }

  /* apply the mirror boundary conditions */
  for (m = 0; m <= this->_SplineDegree; m++) {
    xIndex[m] = (_x == 1) ? (0) : ((xIndex[m] < 0) ? (-xIndex[m] - this->_xhalf * ((-xIndex[m]) / this->_xhalf)) : (xIndex[m] - this->_xhalf * (xIndex[m] / this->_xhalf)));
    if (_x <= xIndex[m]) {
      xIndex[m] = this->_xhalf - xIndex[m];
    }
    yIndex[m] = (_y == 1) ? (0) : ((yIndex[m] < 0) ? (-yIndex[m] - this->_yhalf * ((-yIndex[m]) / this->_yhalf)) : (yIndex[m] - this->_yhalf * (yIndex[m] / this->_yhalf)));
    if (_y <= yIndex[m]) {
      yIndex[m] = this->_yhalf - yIndex[m];
    }
  }

  /* perform interpolation */
  k = round(z);
  l = round(time);
  value = 0.0;
  for (j = 0; j <= this->_SplineDegree; j++) {
    for (i = 0; i <= this->_SplineDegree; i++) {
      value += xWeight[i] * yWeight[j] * this->_coeff(xIndex[i], yIndex[j], k, l);
    }
  }

  if (this->_clamped) {
    if (value > this->_max) return this->_max;
    if (value < this->_min) return this->_min;
  }

  return(value);
}

double irtkBSplineInterpolateImageFunction2D::EvaluateInside(double x, double y, double z, double time)
{
  int i, j, k, l, m;
  int xIndex[6], yIndex[6];
  double xWeight[6], yWeight[6];
  double value, w, w2, w4, t, t0, t1;

  /* compute the interpolation indexes */
  if (_SplineDegree & 1) {
    i = (int)floor(x) - this->_SplineDegree / 2;
    j = (int)floor(y) - this->_SplineDegree / 2;
    for (m = 0; m <= this->_SplineDegree; m++) {
      xIndex[m] = i++;
      yIndex[m] = j++;
    }
  } else {
    i = (int)floor(x + 0.5) - this->_SplineDegree / 2;
    j = (int)floor(y + 0.5) - this->_SplineDegree / 2;
    for (m = 0; m <= this->_SplineDegree; m++) {
      xIndex[m] = i++;
      yIndex[m] = j++;
    }
  }

  /* compute the interpolation weights */
  switch (_SplineDegree) {
  case 2:
    /* x */
    w = x - (double)xIndex[1];
    xWeight[1] = 3.0 / 4.0 - w * w;
    xWeight[2] = (1.0 / 2.0) * (w - xWeight[1] + 1.0);
    xWeight[0] = 1.0 - xWeight[1] - xWeight[2];
    /* y */
    w = y - (double)yIndex[1];
    yWeight[1] = 3.0 / 4.0 - w * w;
    yWeight[2] = (1.0 / 2.0) * (w - yWeight[1] + 1.0);
    yWeight[0] = 1.0 - yWeight[1] - yWeight[2];
    break;
  case 3:
    /* x */
    w = x - (double)xIndex[1];
    xWeight[3] = (1.0 / 6.0) * w * w * w;
    xWeight[0] = (1.0 / 6.0) + (1.0 / 2.0) * w * (w - 1.0) - xWeight[3];
    xWeight[2] = w + xWeight[0] - 2.0 * xWeight[3];
    xWeight[1] = 1.0 - xWeight[0] - xWeight[2] - xWeight[3];
    /* y */
    w = y - (double)yIndex[1];
    yWeight[3] = (1.0 / 6.0) * w * w * w;
    yWeight[0] = (1.0 / 6.0) + (1.0 / 2.0) * w * (w - 1.0) - yWeight[3];
    yWeight[2] = w + yWeight[0] - 2.0 * yWeight[3];
    yWeight[1] = 1.0 - yWeight[0] - yWeight[2] - yWeight[3];
    break;
  case 4:
    /* x */
    w = x - (double)xIndex[2];
    w2 = w * w;
    t = (1.0 / 6.0) * w2;
    xWeight[0] = 1.0 / 2.0 - w;
    xWeight[0] *= xWeight[0];
    xWeight[0] *= (1.0 / 24.0) * xWeight[0];
    t0 = w * (t - 11.0 / 24.0);
    t1 = 19.0 / 96.0 + w2 * (1.0 / 4.0 - t);
    xWeight[1] = t1 + t0;
    xWeight[3] = t1 - t0;
    xWeight[4] = xWeight[0] + t0 + (1.0 / 2.0) * w;
    xWeight[2] = 1.0 - xWeight[0] - xWeight[1] - xWeight[3] - xWeight[4];
    /* y */
    w = y - (double)yIndex[2];
    w2 = w * w;
    t = (1.0 / 6.0) * w2;
    yWeight[0] = 1.0 / 2.0 - w;
    yWeight[0] *= yWeight[0];
    yWeight[0] *= (1.0 / 24.0) * yWeight[0];
    t0 = w * (t - 11.0 / 24.0);
    t1 = 19.0 / 96.0 + w2 * (1.0 / 4.0 - t);
    yWeight[1] = t1 + t0;
    yWeight[3] = t1 - t0;
    yWeight[4] = yWeight[0] + t0 + (1.0 / 2.0) * w;
    yWeight[2] = 1.0 - yWeight[0] - yWeight[1] - yWeight[3] - yWeight[4];
    break;
  case 5:
    /* x */
    w = x - (double)xIndex[2];
    w2 = w * w;
    xWeight[5] = (1.0 / 120.0) * w * w2 * w2;
    w2 -= w;
    w4 = w2 * w2;
    w -= 1.0 / 2.0;
    t = w2 * (w2 - 3.0);
    xWeight[0] = (1.0 / 24.0) * (1.0 / 5.0 + w2 + w4) - xWeight[5];
    t0 = (1.0 / 24.0) * (w2 * (w2 - 5.0) + 46.0 / 5.0);
    t1 = (-1.0 / 12.0) * w * (t + 4.0);
    xWeight[2] = t0 + t1;
    xWeight[3] = t0 - t1;
    t0 = (1.0 / 16.0) * (9.0 / 5.0 - t);
    t1 = (1.0 / 24.0) * w * (w4 - w2 - 5.0);
    xWeight[1] = t0 + t1;
    xWeight[4] = t0 - t1;
    /* y */
    w = y - (double)yIndex[2];
    w2 = w * w;
    yWeight[5] = (1.0 / 120.0) * w * w2 * w2;
    w2 -= w;
    w4 = w2 * w2;
    w -= 1.0 / 2.0;
    t = w2 * (w2 - 3.0);
    yWeight[0] = (1.0 / 24.0) * (1.0 / 5.0 + w2 + w4) - yWeight[5];
    t0 = (1.0 / 24.0) * (w2 * (w2 - 5.0) + 46.0 / 5.0);
    t1 = (-1.0 / 12.0) * w * (t + 4.0);
    yWeight[2] = t0 + t1;
    yWeight[3] = t0 - t1;
    t0 = (1.0 / 16.0) * (9.0 / 5.0 - t);
    t1 = (1.0 / 24.0) * w * (w4 - w2 - 5.0);
    yWeight[1] = t0 + t1;
    yWeight[4] = t0 - t1;
    break;
  default:
    printf("Invalid spline degree\n");
    return(0.0);
  }

  /* perform interpolation */
  k = round(z);
  l = round(time);
  value = 0.0;

  for (j = 0; j <= this->_SplineDegree; j++) {
    for (i = 0; i <= this->_SplineDegree; i++) {
      value += xWeight[i] * yWeight[j] * this->_coeff(xIndex[i], yIndex[j], k, l);
    }
  }


  if (this->_clamped) {
    if (value > this->_max) return this->_max;
    if (value < this->_min) return this->_min;
  }

  return(value);
}

