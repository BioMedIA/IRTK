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

int main(int argc, char **argv)
{
  int i, j, k, n, t;
  double *x1, *y1, *z1, *time, *x2, *y2, *z2;

  // Read image sequence
  irtkGreyImage *image = new irtkGreyImage;
  image->Read(argv[1]);
  argv++;
  argc--;

  // Create transform
  irtkMultiLevelFreeFormTransformation *transform1 = new irtkMultiLevelFreeFormTransformation;

  // Read 4D transform
  cout << "Reading 4D transform " << argv[1] << endl;
  transform1->irtkTransformation::Read(argv[1]);
  argv++;
  argc--;

  // Cast 4D transform
  irtkBSplineFreeFormTransformation4D *ffd4D = (irtkBSplineFreeFormTransformation4D *)transform1->GetLocalTransformation(0);

  // Work out how many time frames there are ...
  n = argc - 2;
  cout << "Approximating " << n << " FFDs" << endl;

  // Create transform
  irtkMultiLevelFreeFormTransformation *transform2 = new irtkMultiLevelFreeFormTransformation;

  // Read first 3D transform
  cout << "Reading 3D transform " << argv[1] << endl;
  transform2->irtkTransformation::Read(argv[1]);
  argv++;
  argc--;

  // Read first transformation
  irtkBSplineFreeFormTransformation *ffd = (irtkBSplineFreeFormTransformation *)transform2->GetLocalTransformation(0);

  x1 = new double[n*image->GetX()*image->GetY()*image->GetZ()];
  y1 = new double[n*image->GetX()*image->GetY()*image->GetZ()];
  z1 = new double[n*image->GetX()*image->GetY()*image->GetZ()];
  x2 = new double[n*image->GetX()*image->GetY()*image->GetZ()];
  y2 = new double[n*image->GetX()*image->GetY()*image->GetZ()];
  z2 = new double[n*image->GetX()*image->GetY()*image->GetZ()];
  time = new double[n*image->GetX()*image->GetY()*image->GetZ()];

  n = 0;
  t = 1;
  for (k = 0; k < image->GetZ(); k += 2) {
    for (j = 0; j < image->GetY(); j += 2) {
      for (i = 0; i < image->GetX(); i += 2) {
        x1[n] = i;
        y1[n] = j;
        z1[n] = k;
        image->ImageToWorld(x1[n], y1[n], z1[n]);
        time[n] = image->ImageToTime(t);
        x2[n] = x1[n];
        y2[n] = y1[n];
        z2[n] = z1[n];
        ffd->LocalDisplacement(x2[n], y2[n], z2[n]);
        n++;
      }
    }
  }
  t++;

  while (argc > 2) {
    cout << "Reading 3D transform " << argv[1] << endl;
    transform2->irtkTransformation::Read(argv[1]);
    argv++;
    argc--;

    // Get transformation
    ffd = (irtkBSplineFreeFormTransformation *)transform2->GetLocalTransformation(0);

    for (k = 0; k < image->GetZ(); k += 2) {
      for (j = 0; j < image->GetY(); j += 2) {
        for (i = 0; i < image->GetX(); i += 2) {
          x1[n] = i;
          y1[n] = j;
          z1[n] = k;
          image->ImageToWorld(x1[n], y1[n], z1[n]);
          time[n] = image->ImageToTime(t);
          x2[n] = x1[n];
          y2[n] = y1[n];
          z2[n] = z1[n];
          ffd->LocalDisplacement(x2[n], y2[n], z2[n]);
          n++;
        }
      }
    }
    t++;
  }

  cout << "Approximating " << n << " points" << endl;
  cout << "Error is " << ffd4D->Approximate(x1, y1, z1, time, x2, y2, z2, n) << endl;

  irtkMultiLevelFreeFormTransformation *transform3 = new irtkMultiLevelFreeFormTransformation;
  transform3->PushLocalTransformation(ffd4D);

  // Write output
  cout << "Writing to " << argv[1] << endl;
  transform3->irtkTransformation::Write(argv[1]);

}
