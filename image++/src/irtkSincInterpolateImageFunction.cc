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

#define SINC_KERNEL_SIZE  6
#define SINC_EPSILON      0.000001
#define SINC_LUTSIZE      1000000

double *SINC_LUT = NULL;

inline double sinc(double x)
{
  // Alternatively, initialize sinc lookup table, if applicable
  if (SINC_LUT == NULL) {
    cerr << "Initializing SINC_LUT ... "; cerr.flush();

    int size;
    double alpha;

    // Allocate lookup table
    size     = round(SINC_LUTSIZE * (SINC_KERNEL_SIZE + 0.5));
    SINC_LUT = new double[size];

    // Value at zero distance
    SINC_LUT[0] = 1.;

    // Fill remaining fields (symmetric, so we only use positive distances)
    for (int i=1; i<size; i++) {
      alpha       = M_PI*double(i)/SINC_LUTSIZE;
      SINC_LUT[i] = 0.5*(1.0+cos(alpha/SINC_KERNEL_SIZE))*sin(alpha)/alpha;
    }

    cerr << "done\n";
  }

  // Return LUT value
  return SINC_LUT[round(fabs(x)*SINC_LUTSIZE)];
}

template <class VoxelType> irtkSincInterpolateImageFunction<VoxelType>::irtkSincInterpolateImageFunction()
{}

template <class VoxelType> irtkSincInterpolateImageFunction<VoxelType>::~irtkSincInterpolateImageFunction(void)
{}

template <class VoxelType> const char *irtkSincInterpolateImageFunction<VoxelType>::NameOfClass()
{
  return "irtkSincInterpolateImageFunction";
}

template <class VoxelType> void irtkSincInterpolateImageFunction<VoxelType>::Initialize()
{
  /// Initialize baseclass
  this->irtkImageFunction<VoxelType>::Initialize();

  this->_x = this->_input->GetX();
  this->_y = this->_input->GetY();
  this->_z = this->_input->GetZ();

  // Compute domain on which the sinc kernel is defined
  this->_x1 = SINC_KERNEL_SIZE;
  this->_y1 = SINC_KERNEL_SIZE;
  this->_z1 = SINC_KERNEL_SIZE;
  this->_x2 = this->_input->GetX() - (SINC_KERNEL_SIZE + 1);
  this->_y2 = this->_input->GetY() - (SINC_KERNEL_SIZE + 1);
  this->_z2 = this->_input->GetZ() - (SINC_KERNEL_SIZE + 1);

  // Compute min and max values
  this->_input->GetMinMax(&this->_min, &this->_max);
}

// Truncated sinc using Hanning window, H(dx/R)*sinc(dx), R=6 where
// sinc(dx) = sin(pi*dx)/(pi*dx), H(dx/R) = 0.5*(1+cos(pi*dx/R))
template <class VoxelType> double irtkSincInterpolateImageFunction<VoxelType>::Evaluate(double x, double y, double z, double t)
{
  int i, j, k, l;

  i = round(x);
  j = round(y);
  k = round(z);
  l = round(t);

  // Check whether there is anything to interpolate
  if ((fabs(x-i) > SINC_EPSILON) ||
      (fabs(y-j) > SINC_EPSILON) ||
      (fabs(z-k) > SINC_EPSILON )) {

    int i1, i2, i3, j1, j2, j3, k1, k2, k3;
    double wx, wy, wz, w, val, sum;

    i1  = i - SINC_KERNEL_SIZE;
    i2  = i + SINC_KERNEL_SIZE + 1;
    j1  = j - SINC_KERNEL_SIZE;
    j2  = j + SINC_KERNEL_SIZE + 1;
    k1  = k - SINC_KERNEL_SIZE;
    k2  = k + SINC_KERNEL_SIZE + 1;
    val = 0;
    sum = 0;
    for (i3 = i1; i3 < i2; i3++) {
      if ((i3 >= 0) && (i3 < this->_x)) {
        wx = sinc(i3 - x);
        for (j3 = j1; j3 < j2; j3++) {
          if ((j3 >= 0) && (j3 < this->_y)) {
            wy = sinc(j3 - y);
            for (k3 = k1; k3 < k2; k3++) {
              if ((k3 >= 0) && (k3 < this->_z)) {
                wz   = sinc(k3 - z);
                w    = wx*wy*wz;
                val += w*this->_input->Get(i3, j3, k3, l);
                sum += w;
              }
            }
          }
        }
      }
    }

    if (sum != 0) {
      val /= sum;
    } else {
      val = 0;
    }

    if (this->_clamped) {
      if (val > this->_max) return this->_max;
      if (val < this->_min) return this->_min;
    }

    return(val);

  } else {
    // Return nearest neighbour
    return this->_input->Get(i, j, k, l);
  }
}

// Truncated sinc using Hanning window, H(dx/R)*sinc(dx), R=6 where
// sinc(dx) = sin(pi*dx)/(pi*dx), H(dx/R) = 0.5*(1+cos(pi*dx/R))
template <class VoxelType> double irtkSincInterpolateImageFunction<VoxelType>::EvaluateInside(double x, double y, double z, double t)
{
  int i, j, k, l;

  i = round(x);
  j = round(y);
  k = round(z);
  l = round(t);

  // Check whether there is anything to interpolate
  if ((fabs(x-i) > SINC_EPSILON) ||
      (fabs(y-j) > SINC_EPSILON) ||
      (fabs(z-k) > SINC_EPSILON )) {

    int i1, i2, i3, j1, j2, j3, k1, k2, k3;
    double wx, wy, wz, w, val, sum;

    i1  = i - SINC_KERNEL_SIZE;
    i2  = i + SINC_KERNEL_SIZE + 1;
    j1  = j - SINC_KERNEL_SIZE;
    j2  = j + SINC_KERNEL_SIZE + 1;
    k1  = k - SINC_KERNEL_SIZE;
    k2  = k + SINC_KERNEL_SIZE + 1;
    val = 0;
    sum = 0;
    for (i3 = i1; i3 < i2; i3++) {
      wx = sinc(i3 - x);
      for (j3 = j1; j3 < j2; j3++) {
        wy = sinc(j3 - y);
        for (k3 = k1; k3 < k2; k3++) {
          wz   = sinc(k3 - z);
          w    = wx*wy*wz;
          val += w*this->_input->Get(i3, j3, k3, l);
          sum += w;
        }
      }
    }

    if (sum != 0) {
      val /= sum;
    } else {
      val = 0;
    }

    if (this->_clamped) {
      if (val > this->_max) return this->_max;
      if (val < this->_min) return this->_min;
    }

    return(val);

  } else {
    // Return nearest neighbour
    return this->_input->Get(round(t), i, j, k);
  }
}

template class irtkSincInterpolateImageFunction<irtkBytePixel>;
template class irtkSincInterpolateImageFunction<irtkGreyPixel>;
template class irtkSincInterpolateImageFunction<irtkRealPixel>;
