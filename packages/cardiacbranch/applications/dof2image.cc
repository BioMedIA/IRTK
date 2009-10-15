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
char *image_name = NULL, *input_name = NULL, *output_name_x = NULL, *output_name_y = NULL, *output_name_z = NULL, *total_name = NULL;

void usage()
{
  cerr << "Usage: dof2image [image] [dof] [dx] [dy] [dz] <options>\n" << endl;
  cerr << "where <options> is one or more of the following:\n" << endl;
  cerr << "<-invert>                  Store the inverted displacement field" << endl;
  cerr << "<-total image>             Store the total displacement in image" << endl;
  cerr << "<-scale factor>            Scale displacement by a factor" << endl;
  cerr << "<-padding value>           Ignore padded region" << endl;
  cerr << "<-Rx1 value>               Region of interest in image" << endl;
  cerr << "<-Ry1 value>               Region of interest in image" << endl;
  cerr << "<-Rz1 value>               Region of interest in image" << endl;
  cerr << "<-Rx2 value>               Region of interest in image" << endl;
  cerr << "<-Ry2 value>               Region of interest in image" << endl;
  cerr << "<-Rz2 value>               Region of interest in image" << endl;
  exit(1);
}

int main(int argc, char **argv)
{
  int x, y, z, t, ok, invert;
  int x1, y1, z1, t1, x2, y2, z2, t2;
  irtkGreyPixel padding;
  double p1[3], p2[3], scale;

  if (argc < 6) {
    usage();
  }

  // Parse file names
  image_name  = argv[1];
  argc--;
  argv++;
  input_name  = argv[1];
  argc--;
  argv++;
  output_name_x = argv[1];
  argc--;
  argv++;
  output_name_y = argv[1];
  argc--;
  argv++;
  output_name_z = argv[1];
  argc--;
  argv++;

  // Read image
  irtkGreyImage image(image_name);

  // Read transformation
  irtkTransformation *transform = irtkTransformation::New(input_name);

  // Fix ROI
  x1 = 0;
  y1 = 0;
  z1 = 0;
  t1 = 0;
  x2 = image.GetX();
  y2 = image.GetY();
  z2 = image.GetZ();
  t2 = image.GetT();

  // Initialize parameters
  scale   = 1;
  padding = MIN_GREY;
  invert  = False;

  // Parse arguments
  while (argc > 1) {
    ok = False;
    if ((ok == False) && (strcmp(argv[1], "-padding") == 0)) {
      argc--;
      argv++;
      padding = atoi(argv[1]);
      argc--;
      argv++;
      ok = True;
    }
    if ((ok == False) && (strcmp(argv[1], "-total") == 0)) {
      argc--;
      argv++;
      total_name = argv[1];
      argc--;
      argv++;
      ok = True;
    }
    if ((ok == False) && (strcmp(argv[1], "-scale") == 0)) {
      argc--;
      argv++;
      scale = atof(argv[1]);
      argc--;
      argv++;
      ok = True;
    }
    if ((ok == False) && (strcmp(argv[1], "-invert") == 0)) {
      argc--;
      argv++;
      invert = True;
      ok = True;
    }
    if ((ok == False) && (strcmp(argv[1], "-Rx1") == 0)) {
      argc--;
      argv++;
      x1 = atoi(argv[1]);
      argc--;
      argv++;
      ok = True;
    }
    if ((ok == False) && (strcmp(argv[1], "-Rx2") == 0)) {
      argc--;
      argv++;
      x2 = atoi(argv[1]);
      argc--;
      argv++;
      ok = True;
    }
    if ((ok == False) && (strcmp(argv[1], "-Ry1") == 0)) {
      argc--;
      argv++;
      y1 = atoi(argv[1]);
      argc--;
      argv++;
      ok = True;
    }
    if ((ok == False) && (strcmp(argv[1], "-Ry2") == 0)) {
      argc--;
      argv++;
      y2 = atoi(argv[1]);
      argc--;
      argv++;
      ok = True;
    }
    if ((ok == False) && (strcmp(argv[1], "-Rz1") == 0)) {
      argc--;
      argv++;
      z1 = atoi(argv[1]);
      argc--;
      argv++;
      ok = True;
    }
    if ((ok == False) && (strcmp(argv[1], "-Rz2") == 0)) {
      argc--;
      argv++;
      z2 = atoi(argv[1]);
      argc--;
      argv++;
      ok = True;
    }
    if ((ok == False) && (strcmp(argv[1], "-Rt1") == 0)) {
      argc--;
      argv++;
      t1 = atoi(argv[1]);
      argc--;
      argv++;
      ok = True;
    }
    if ((ok == False) && (strcmp(argv[1], "-Rt2") == 0)) {
      argc--;
      argv++;
      t2 = atoi(argv[1]);
      argc--;
      argv++;
      ok = True;
    }
    if (ok == False) {
      cerr << "Can not parse argument " << argv[1] << endl;
      usage();
    }
  }

  // If there is an region of interest, use it
  if ((x1 != 0) || (x2 != image.GetX()) ||
      (y1 != 0) || (y2 != image.GetY()) ||
      (z1 != 0) || (z2 != image.GetZ()) ||
      (t1 != 0) || (t2 != image.GetZ())) {
    image = image.GetRegion(x1, y1, z1, t1, x2, y2, z2, t2);
  }

  // Create images with displacements
  irtkRealImage dx(image);
  irtkRealImage dy(image);
  irtkRealImage dz(image);

  // Initialize point structure with transformed point positions
  for (t = 0; t < image.GetT(); t++) {
    double time = image.ImageToTime(t);
    for (z = 0; z < image.GetZ(); z++) {
      for (y = 0; y < image.GetY(); y++) {
        for (x = 0; x < image.GetX(); x++) {
          if (image(x, y, z, t) > padding) {
            p1[0] = x;
            p1[1] = y;
            p1[2] = z;
            image.ImageToWorld(p1[0], p1[1], p1[2]);
            p2[0] = p1[0];
            p2[1] = p1[1];
            p2[2] = p1[2];
            if (invert == True) {
              transform->Inverse(p2[0], p2[1], p2[2], time);
            } else {
              transform->Transform(p2[0], p2[1], p2[2], time);
            }
            p2[0] -= p1[0];
            p2[1] -= p1[1];
            p2[2] -= p1[2];
            dx(x, y, z, t) = scale * p2[0];
            dy(x, y, z, t) = scale * p2[1];
            dz(x, y, z, t) = scale * p2[2];
          } else {
            dx(x, y, z, t) = 0;
            dy(x, y, z, t) = 0;
            dz(x, y, z, t) = 0;
          }
        }
      }
    }
  }

  dx.Write(output_name_x);
  dy.Write(output_name_y);
  dz.Write(output_name_z);

  if (total_name != NULL) {
    irtkRealImage total(image);
    for (t = 0; t < image.GetT(); t++) {
      for (z = 0; z < image.GetZ(); z++) {
        for (y = 0; y < image.GetY(); y++) {
          for (x = 0; x < image.GetX(); x++) {
            if (image(x, y, z, t) > padding) {
              total(x, y, z, t) = sqrt(dx(x, y, z, t) * dx(x, y, z, t) + dy(x, y, z, t) * dy(x, y, z, t) + dz(x, y, z, t) * dz(x, y, z, t));
            } else {
              total(x, y, z, t) = 0;
            }
          }
        }
      }
    }
    total.Write(total_name);
  }
}
