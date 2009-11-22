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

#include <irtkRegionFilter.h>

irtkRegionFilter::irtkRegionFilter()
{

}

irtkRegionFilter::~irtkRegionFilter(void)
{
}

Bool irtkRegionFilter::RequiresBuffering(void)
{
  return True;
}

const char *irtkRegionFilter::NameOfClass()
{
  return "irtkRegionFilter";
}

void irtkRegionFilter::Run()
{
  int i, j, k, l;
  double x1, y1, z1, x2, y2, z2;

  // Do the initial set up
  this->Initialize();

  if ((_i1 < 0) || (_i1 >= _i2) ||
      (_j1 < 0) || (_j1 >= _j2) ||
      (_k1 < 0) || (_k1 >= _k2) ||
      (_l1 < 0) || (_l1 >= _l2) ||
      (_i2 > _input->GetX()) || (_j2 > _input->GetY()) || (_k2 > _input->GetZ()) || (_l2 > _input->GetT())) {
    cerr << "irtkRegionFilter::Run: Parameter out of range\n";
    exit(1);
  }

  // Initialize
  irtkImageAttributes attr = _input->GetImageAttributes();
  attr._x = _i2 - _i1;
  attr._y = _j2 - _j1;
  attr._z = _k2 - _k1;
  attr._t = _l2 - _l1;
  attr._xorigin = 0;
  attr._yorigin = 0;
  attr._zorigin = 0;
  _output->Initialize(attr);

  // Calculate position of first voxel in roi in original image
  x1 = _i1;
  y1 = _j1;
  z1 = _k1;
  _input->ImageToWorld(x1, y1, z1);

  // Calculate position of first voxel in roi in new image
  x2 = 0;
  y2 = 0;
  z2 = 0;
  _output->ImageToWorld(x2, y2, z2);

  // Shift origin of new image accordingly
  _output->PutOrigin(x1 - x2, y1 - y2, z1 - z2);

  // Copy region
  for (l = _l1; l < _l2; l++) {
    for (k = _k1; k < _k2; k++) {
      for (j = _j1; j < _j2; j++) {
        for (i = _i1; i < _i2; i++) {
          _output->PutAsDouble(i-_i1, j-_j1, k-_k1, l-_l1, _input->GetAsDouble(i, j, k, l));
        }
      }
    }
  }

  // Do the final cleaning up
  this->Finalize();
}
