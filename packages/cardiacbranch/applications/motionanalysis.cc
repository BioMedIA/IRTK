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
char *image_name = NULL, *input_name = NULL, *output_name = NULL, *total_name = NULL;

typedef enum { Displacement, RadialDisplacement, CircumferentialDisplacement } irtkDeformationToImageMode;

void usage()
{
  cerr << "Usage: motionanalysis [image] [input dof] [output image]<options>\n" << endl;
  cerr << "where <options> is one or more of the following:\n" << endl;
  cerr << "<-radial>            Store the radial displacement" << endl;
  cerr << "<-displacement>      Store the displacement vector" << endl;
  cerr << "<-circumferential>   Store the circumferential displacement" << endl;
  cerr << "<-origin x y z>      Use x, y, z in world coordinates as origin" << endl;
  cerr << "<-scale factor>      Scale values by a factor" << endl;
  cerr << "<-padding value>     Ignore padded region" << endl;
  exit(1);
}

int main(int argc, char **argv)
{
  int x, y, z, ok;
  irtkGreyPixel padding;
  irtkDeformationToImageMode mode;
  double p1[3], p2[3], q1[3], q2[3], origin[3], r1, r2, scale;

  if (argc < 4) {
    usage();
  }

  // Parse file names
  image_name  = argv[1];
  argc--;
  argv++;
  input_name  = argv[1];
  argc--;
  argv++;
  output_name = argv[1];
  argc--;
  argv++;

  // Initialize parameters
  scale     = 1;
  padding   = MIN_GREY;
  mode      = Displacement;
  origin[0] = 0;
  origin[1] = 0;
  origin[2] = 0;

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
    if ((ok == False) && (strcmp(argv[1], "-displacement") == 0)) {
      argc--;
      argv++;
      mode = Displacement;
      ok = True;
    }
    if ((ok == False) && (strcmp(argv[1], "-radial") == 0)) {
      argc--;
      argv++;
      mode = RadialDisplacement;
      ok = True;
    }
    if ((ok == False) && (strcmp(argv[1], "-circumferential") == 0)) {
      argc--;
      argv++;
      mode = CircumferentialDisplacement;
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
    if ((ok == False) && (strcmp(argv[1], "-origin") == 0)) {
      argc--;
      argv++;
      origin[0] = atof(argv[1]);
      argc--;
      argv++;
      origin[1] = atof(argv[1]);
      argc--;
      argv++;
      origin[2] = atof(argv[1]);
      argc--;
      argv++;
      ok = True;
    }
    if (ok == False) {
      cerr << "Can not parse argument " << argv[1] << endl;
      usage();
    }
  }

  // Read image
  irtkGreyImage image(image_name);

  // Read transformation
  irtkTransformation *transform = irtkTransformation::New(input_name);

  // Create images with displacements
  irtkRealImage output;
  output = image;

  // Initialize point structure with transformed point positions
  for (z = 0; z < image.GetZ(); z++) {
    for (y = 0; y < image.GetY(); y++) {
      for (x = 0; x < image.GetX(); x++) {
        if (image(x, y, z) > padding) {
          if (z > 0) {
            p1[0] = x;
            p1[1] = y;
            p1[2] = z-1;
            image.ImageToWorld(p1[0], p1[1], p1[2]);
            p2[0] = p1[0];
            p2[1] = p1[1];
            p2[2] = p1[2];
            transform->Transform(p2[0], p2[1], p2[2]);
          } else {
            p1[0] = x;
            p1[1] = y;
            p1[2] = z;
            image.ImageToWorld(p1[0], p1[1], p1[2]);
            p2[0] = p1[0];
            p2[1] = p1[1];
            p2[2] = p1[2];
          }
          q1[0] = x;
          q1[1] = y;
          q1[2] = z;
          image.ImageToWorld(q1[0], q1[1], q1[2]);
          q2[0] = q1[0];
          q2[1] = q1[1];
          q2[2] = q1[2];
          transform->Transform(q2[0], q2[1], q2[2]);
          switch (mode) {
          case Displacement:
            q2[0] -= p2[0];
            q2[1] -= p2[1];
            q2[2] -= p2[2];
            output(x, y, z) = scale * sqrt(q2[0]*q2[0] + q2[1]*q2[1]);
            break;
          case RadialDisplacement:
            break;
          case CircumferentialDisplacement:
            r1 = atan2(p1[1]-origin[1], p1[0]-origin[0]);
            r2 = atan2(p2[1]-origin[1], p2[0]-origin[0]);
            break;
          }
        } else {
          output(x, y, z) = 0;
        }
      }
    }
  }
  output.Write(output_name);
}
