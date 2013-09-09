/*=========================================================================

  Library   : Image Registration Toolkit (IRTK)
  Module    : $Id$
  Copyright : Imperial College, Department of Computing
              Visual Information Processing (VIP), 2011 onwards
  Date      : $Date$
  Version   : $Revision$
  Changes   : $Author$

=========================================================================*/



#include <irtkBiasCorrectionMask.h>

irtkBiasCorrectionMask::irtkBiasCorrectionMask()
{
  // mask
  _mask = NULL;
}

irtkBiasCorrectionMask::~irtkBiasCorrectionMask()
{
}

void irtkBiasCorrectionMask::SetMask( irtkRealImage* imagePtr)
{
	_mask = imagePtr;
}

void irtkBiasCorrectionMask::Run()
{
  int i, j, k, n;
  irtkRealPixel *ptr2target, *ptr2ref, *ptr2w;

  if (_reference == NULL) {
    cerr << "irtkBiasCorrectionMask::Run: Filter has no reference input" << endl;
    exit(1);
  }

  if (_target == NULL) {
    cerr << "irtkBiasCorrectionMask::Run: Filter has no target input" << endl;
    exit(1);
  }

  if (_biasfield == NULL) {
    cerr << "irtkBiasCorrectionMask::Run: Filter has no transformation output" << endl;
    exit(1);
  }

  // Do the initial set up for all levels
  this->Initialize();

  // Compute no of unpadded voxels
  n = 0;
  ptr2target = _target->GetPointerToVoxels();
  irtkRealPixel *pm = _mask->GetPointerToVoxels();
  for (i = 0; i < _target->GetNumberOfVoxels(); i++) {
    if (*ptr2target != _Padding && *pm == 1 ) {
      n++;
    }
    pm++;
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
  pm = _mask->GetPointerToVoxels();
  for (k = 0; k < _target->GetZ(); k++) {
    for (j = 0; j < _target->GetY(); j++) {
      for (i = 0; i < _target->GetX(); i++) {
        if (*ptr2target != _Padding && *pm == 1) {
          x[n] = i;
          y[n] = j;
          z[n] = k;
          _target->ImageToWorld(x[n], y[n], z[n]);
          b[n] = *ptr2target - (double) *ptr2ref;
          w[n] = *ptr2w;
          n++;
        }
        pm++;
        ptr2target++;
        ptr2ref++;
        ptr2w++;
      }
    }
  }

  cout << "Computing bias field ... ";
  cout.flush();

  _biasfield->WeightedLeastSquares(x, y, z, b, w, n);
  cout << "done" << endl;

  delete []x;
  delete []y;
  delete []z;
  delete []b;

  // Do the final cleaning up for all levels
  this->Finalize();
}

