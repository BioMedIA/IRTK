/*=========================================================================

  Library   : Image Registration Toolkit (IRTK)
  Module    : $Id$
  Copyright : Imperial College, Department of Computing
              Visual Information Processing (VIP), 2008 onwards
  Date      : $Date$
  Version   : $Revision$
  Changes   : $Author$

=========================================================================*/

#include <irtkBiasCorrection.h>

irtkBiasCorrection::irtkBiasCorrection()
{
  // Set parameters
  _Padding   = MIN_GREY;

  // Set inputs
  _target    = NULL;
  _reference = NULL;

  // Set output
  _biasfield = NULL;
}

irtkBiasCorrection::~irtkBiasCorrection()
{
}

void irtkBiasCorrection::Initialize()
{
}

void irtkBiasCorrection::Finalize()
{
}

void irtkBiasCorrection::Run()
{
  int i, j, k, n;
  irtkRealPixel *ptr2target, *ptr2ref, *ptr2w;

  if (_reference == NULL) {
    cerr << "BiasCorrection::Run: Filter has no reference input" << endl;
    exit(1);
  }

  if (_target == NULL) {
    cerr << "BiasCorrection::Run: Filter has no target input" << endl;
    exit(1);
  }

  if (_biasfield == NULL) {
    cerr << "irtkBiasCorrection::Run: Filter has no transformation output" << endl;
    exit(1);
  }

  // Do the initial set up for all levels
  this->Initialize();

  // Compute no of unpadded voxels
  n = 0;
  ptr2target = _target->GetPointerToVoxels();
  for (i = 0; i < _target->GetNumberOfVoxels(); i++) {
    if (*ptr2target != _Padding) {
      n++;
    }
    ptr2target++;
  }

  double *x = new double[n];
  double *y = new double[n];
  double *z = new double[n];
  double *b = new double[n];
  double *w = new double[n];
  n = 0;
  ptr2target = _target->GetPointerToVoxels();
  ptr2ref    = _reference->GetPointerToVoxels();
  ptr2w      = _weights->GetPointerToVoxels();
  for (k = 0; k < _target->GetZ(); k++) {
    for (j = 0; j < _target->GetY(); j++) {
      for (i = 0; i < _target->GetX(); i++) {
        if (*ptr2target != _Padding) {
          x[n] = i;
          y[n] = j;
          z[n] = k;
          _target->ImageToWorld(x[n], y[n], z[n]);
          b[n] = *ptr2target - (double) *ptr2ref;
          w[n] = *ptr2w;
          n++;
        }
        ptr2target++;
        ptr2ref++;
        ptr2w++;
      }
    }
  }

  cout << "Computing bias field ... ";
  cout.flush();
// _biasfield->Approximate(x, y, z, b, n);
  _biasfield->WeightedLeastSquares(x, y, z, b, w, n);
  cout << "done" << endl;

  delete []x;
  delete []y;
  delete []z;
  delete []b;

  // Do the final cleaning up for all levels
  this->Finalize();
}

void irtkBiasCorrection::Apply(irtkRealImage &image)
{
  int i, j, k;
  double x, y, z, bias;

  // Default is the target image
  image = *_target;

  for (k = 0; k < image.GetZ(); k++) {
    for (j = 0; j < image.GetY(); j++) {
      for (i = 0; i < image.GetX(); i++) {
        x = i;
        y = j;
        z = k;
        image.ImageToWorld(x, y, z);
        bias = _biasfield->Bias(x, y, z);
        if (_target->Get(i, j, k) != _Padding) {
          image(i, j, k) = round(image(i, j, k) - bias);
        }
      }
    }
  }
}

void irtkBiasCorrection::ApplyToImage(irtkRealImage &image)
{
  int i, j, k;
  double x, y, z, bias;
  cerr<<"Applying bias ...";
  cerr<<_biasfield;

  for (k = 0; k < image.GetZ(); k++) {
    for (j = 0; j < image.GetY(); j++) {
      for (i = 0; i < image.GetX(); i++) {
        x = i;
        y = j;
        z = k;
        image.ImageToWorld(x, y, z);
        bias = _biasfield->Bias(x, y, z);
        if (image.Get(i, j, k) != _Padding) {
          image(i, j, k) = round(image(i, j, k) - bias);
        }
      }
    }
  }

  cerr<<"done."<<endl;
}

void irtkBiasCorrection::ApplyToImage(irtkGreyImage &image)
{
  int i, j, k;
  double x, y, z, bias;

  for (k = 0; k < image.GetZ(); k++) {
    for (j = 0; j < image.GetY(); j++) {
      for (i = 0; i < image.GetX(); i++) {
        x = i;
        y = j;
        z = k;
        image.ImageToWorld(x, y, z);
        bias = _biasfield->Bias(x, y, z);
        if (image.Get(i, j, k) != _Padding) {
          image(i, j, k) = round(image(i, j, k) / exp(bias/1000));
        }
      }
    }
  }
}


