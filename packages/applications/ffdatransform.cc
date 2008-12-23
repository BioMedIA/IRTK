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

#include <irtkTransformation.h>

void usage()
{
  cerr << "Usage: ffdatransform [dofin] [affine transform] [dofout]" << endl;
  exit(1);
}

int main(int argc, char **argv)
{
  int i, j, k, level;
  double x1, y1, z1, x2, y2, z2, xaxis[3], yaxis[3];

  if (argc != 4) {
    usage();
  }

  // Read non-rigid transformation
  irtkMFreeFormTransformation *mffd = new irtkMFreeFormTransformation;
  mffd->Read(argv[1]);

  // Read affine transformation
  irtkAffineTransformation *transform = new irtkAffineTransformation;
  transform->Read(argv[2]);

  for (level = 0; level < mffd->NumberOfLevels(); level++) {
    irtkBSplineFreeFormTransformation *affd = mffd->GetLocalTransformation(level);

    // Get bounding box and orientation
    affd->BoundingBox(x1, y1, z1, x2, y2, z2);
    affd->GetOrientation(xaxis, yaxis);

    // Transform orientation
    transform->Rotate(xaxis[0], xaxis[1], xaxis[2]);
    transform->Rotate(yaxis[0], yaxis[1], yaxis[2]);

    // Set orientation
    affd->PutOrientation(xaxis, yaxis);

    // Transform bounding box
    transform->Transform(x1, y1, z1);
    transform->Transform(x2, y2, z2);

    // Set bounding box
    affd->PutBoundingBox(x1, y1, z1, x2, y2, z2);

    // Get transformation matrix
    irtkMatrix m = transform->GetMatrix();

    for (k = 0; k < affd->GetZ(); k++) {
      for (j = 0; j < affd->GetY(); j++) {
        for (i = 0; i < affd->GetX(); i++) {
          affd->Get(i, j, k, x1, y1, z1);
          x2 = m(0, 0)*x1 + m(0, 1)*y1 + m(0, 2)*z1;
          y2 = m(1, 0)*x1 + m(1, 1)*y1 + m(1, 2)*z1;
          z2 = m(2, 0)*x1 + m(2, 1)*y1 + m(2, 2)*z1;
          affd->Put(i, j, k, x2, y2, z2);
        }
      }
    }
  }

  mffd->Write(argv[3]);
}







