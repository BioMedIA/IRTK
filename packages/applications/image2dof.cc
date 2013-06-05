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

// Default filenames
char *input_name = NULL, *input_name_x = NULL, *input_name_y = NULL, *input_name_z = NULL, *output_name = NULL;

void usage()
{
  cerr << "Usage: image2dof [image] [dx] [dy] [dz] [dof] <options> \n" << endl;
  cerr << "where <options> is one or more of the following:\n" << endl;
  cerr << "<-bspline spacing>        Fit an approximation using a BSpline FFD with the specified spacing in mm" << endl;
  exit(1);
}

int main(int argc, char **argv)
{
  int i, j, k, ok, bspline;
  double x, y, z, x1, y1, z1, x2, y2, z2, xsize, ysize, zsize, xaxis[3], yaxis[3], zaxis[3], spacing;

  if (argc < 3) {
    usage();
  }

  // Parse file names
  input_name   = argv[1];
  argc--;
  argv++;
  input_name_x = argv[1];
  argc--;
  argv++;
  input_name_y = argv[1];
  argc--;
  argv++;
  input_name_z = argv[1];
  argc--;
  argv++;
  output_name  = argv[1];
  argc--;
  argv++;

  bspline = false;
  while (argc > 1) {
    ok = false;
    if ((ok == false) && (strcmp(argv[1], "-bspline") == 0)) {
	  argc--;
	  argv++;
	  spacing = atof(argv[1]);
	  argc--;
	  argv++;
	  bspline = true;
	  ok = true;
    }
    if (ok == false) {
	  cerr << "Can not parse argument " << argv[1] << endl;
	  usage();
    }
  }

  // Read image
  irtkGreyImage image;
  image.Read(input_name);

  // Read images with displacements
  irtkRealImage dx;
  irtkRealImage dy;
  irtkRealImage dz;
  dx.Read(input_name_x);
  dy.Read(input_name_y);
  dz.Read(input_name_z);

  // Compute bounding box of data
  x1 = 0;
  y1 = 0;
  z1 = 0;
  x2 = dx.GetX()-1;
  y2 = dx.GetY()-1;
  z2 = dx.GetZ()-1;
  dx.ImageToWorld(x1, y1, z1);
  dx.ImageToWorld(x2, y2, z2);
  dx.GetOrientation(xaxis, yaxis, zaxis);
  dx.GetPixelSize(&xsize, &ysize, &zsize);

  // Create transformation
  irtkMultiLevelFreeFormTransformation *transformation = new irtkMultiLevelFreeFormTransformation;
  irtkFreeFormTransformation3D *ffd;

  if (bspline) {
	int no = dx.GetX() * dx.GetY() * dx.GetZ();
	int index;
	double dispX, dispY, dispZ, point2X, point2Y, point2Z;
	double *pointsX = new double[no];
	double *pointsY = new double[no];
	double *pointsZ = new double[no];
	double *displacementsX = new double[no];
	double *displacementsY = new double[no];
	double *displacementsZ = new double[no];

	//Create deformation
	ffd = new irtkBSplineFreeFormTransformation(x1, y1, z1, x2, y2, z2, spacing, spacing, spacing, xaxis, yaxis, zaxis);

	index = 0;
	for (k = 0; k < image.GetZ(); k++) {
	  for (j = 0; j < image.GetY(); j++) {
		for (i = 0; i < image.GetX(); i++) {
		  x = i;
		  y = j;
		  z = k;
		  image.ImageToWorld(x, y, z);
		  pointsX[index] = x;
		  pointsY[index] = y;
		  pointsZ[index] = z;

		  point2X = i + dx.Get(i, j, k);
		  point2Y = j + dy.Get(i, j, k);
		  point2Z = k + dz.Get(i, j, k);
		  image.ImageToWorld(point2X, point2Y, point2Z);

		  dispX = point2X - x;
		  dispY = point2Y - y;
		  dispZ = point2Z - z;

		  displacementsX[index] = dispX;
		  displacementsY[index] = dispY;
		  displacementsZ[index] = dispZ;

		  index++;
		}
	  }
	}

	ffd->Approximate(pointsX, pointsY, pointsZ, displacementsX, displacementsY, displacementsZ, no);

  } else {
    // Create deformation
    ffd = new irtkLinearFreeFormTransformation(x1, y1, z1, x2, y2, z2, xsize, ysize, zsize, xaxis, yaxis, zaxis);

    // Initialize point structure with transformed point positions
    for (k = 0; k < dx.GetZ(); k++) {
	  for (j = 0; j < dx.GetY(); j++) {
		for (i = 0; i < dx.GetX(); i++) {
          ffd->Put(i, j, k, dx(i, j, k), dy(i, j, k), dz(i, j, k));
        }
      }
    }
  }

  // Add deformation
  transformation->PushLocalTransformation(ffd);

  // Write transformation
  transformation->irtkTransformation::Write(output_name);
}
