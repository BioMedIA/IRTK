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
char *input_name = NULL, *output_name = NULL, *dof1_name  = NULL, *dof2_name  = NULL;

void usage()
{
  cerr << "Usage: transformation2 [source] [output] <options>\n" << endl;
  cerr << "where <options> is one or more of the following:\n" << endl;
  cerr << "<-dofin1 file>     First transformation" << endl;
  cerr << "<-dofin2 file>     Second transformation" << endl;
  cerr << "<-target file>     Target image" << endl;
  cerr << "<-Rx1 value>       Region of interest" << endl;
  cerr << "<-Ry1 value>       Region of interest" << endl;
  cerr << "<-Rz1 value>       Region of interest" << endl;
  cerr << "<-Rx2 value>       Region of interest" << endl;
  cerr << "<-Ry2 value>       Region of interest" << endl;
  cerr << "<-Rz2 value>       Region of interest" << endl;
  cerr << "<-Tx1 value>       Region of interest in target image" << endl;
  cerr << "<-Ty1 value>       Region of interest in target image" << endl;
  cerr << "<-Tz1 value>       Region of interest in target image" << endl;
  cerr << "<-Tx2 value>       Region of interest in target image" << endl;
  cerr << "<-Ty2 value>       Region of interest in target image" << endl;
  cerr << "<-Tz2 value>       Region of interest in target image" << endl;
  cerr << "<-Sx1 value>       Region of interest in source image" << endl;
  cerr << "<-Sy1 value>       Region of interest in source image" << endl;
  cerr << "<-Sz1 value>       Region of interest in source image" << endl;
  cerr << "<-Sx2 value>       Region of interest in source image" << endl;
  cerr << "<-Sy2 value>       Region of interest in source image" << endl;
  cerr << "<-Sz2 value>       Region of interest in source image" << endl;
  cerr << "<-Tp value>        Target padding value" << endl;
  cerr << "<-Sp value>        Source padding value" << endl;
  cerr << "<-invert>          Invert transformation" << endl;
  cerr << "<-linear>          Linear interpolation" << endl;
  cerr << "<-bspline>         B-spline interpolation" << endl;
  cerr << "<-cspline>         Cubic spline interpolation" << endl;
  cerr << "<-sinc>            Sinc interpolation" << endl;
  exit(1);
}

int main(int argc, char **argv)
{
  int ok, invert;
  int target_padding, source_padding;
  int target_x1, target_y1, target_z1, target_x2, target_y2, target_z2;
  int source_x1, source_y1, source_z1, source_x2, source_y2, source_z2;
  irtkTransformation *transformation1 = NULL;
  irtkTransformation *transformation2 = NULL;
  irtkGreyImage *source = NULL;
  irtkGreyImage *target = NULL;
  irtkImageFunction<irtkGreyPixel> *interpolator = NULL;

  // Check command line
  if (argc < 3) {
    usage();
  }

  // Parse image
  input_name  = argv[1];
  argc--;
  argv++;
  output_name = argv[1];
  argc--;
  argv++;

  // Read image
  cout << "Reading image ... "; cout.flush();
  source = new irtkGreyImage(input_name);
  cout << "done" << endl;

  // Fix ROI
  target_x1 = -1;
  target_y1 = -1;
  target_z1 = -1;
  target_x2 = -1;
  target_y2 = -1;
  target_z2 = -1;
  source_x1 = -1;
  source_y1 = -1;
  source_z1 = -1;
  source_x2 = -1;
  source_y2 = -1;
  source_z2 = -1;

  // Other options
  invert = False;
  source_padding = 0;
  target_padding = MIN_GREY;

  while (argc > 1) {
    ok = False;
    if ((ok == False) && (strcmp(argv[1], "-dofin1") == 0)) {
      argc--;
      argv++;
      dof1_name = argv[1];
      argc--;
      argv++;
      ok = True;
    }
    if ((ok == False) && (strcmp(argv[1], "-dofin2") == 0)) {
      argc--;
      argv++;
      dof2_name = argv[1];
      argc--;
      argv++;
      ok = True;
    }
    if ((ok == False) && (strcmp(argv[1], "-target") == 0)) {
      argc--;
      argv++;
      target = new irtkGreyImage(argv[1]);
      argc--;
      argv++;
      ok = True;
    }
    if ((ok == False) && (strcmp(argv[1], "-Rx1") == 0)) {
      argc--;
      argv++;
      target_x1 = source_x1 = atoi(argv[1]);
      argc--;
      argv++;
      ok = True;
    }
    if ((ok == False) && (strcmp(argv[1], "-Rx2") == 0)) {
      argc--;
      argv++;
      target_x2 = source_x2 = atoi(argv[1]);
      argc--;
      argv++;
      ok = True;
    }
    if ((ok == False) && (strcmp(argv[1], "-Ry1") == 0)) {
      argc--;
      argv++;
      target_y1 = source_y1 = atoi(argv[1]);
      argc--;
      argv++;
      ok = True;
    }
    if ((ok == False) && (strcmp(argv[1], "-Ry2") == 0)) {
      argc--;
      argv++;
      target_y2 = source_y2 = atoi(argv[1]);
      argc--;
      argv++;
      ok = True;
    }
    if ((ok == False) && (strcmp(argv[1], "-Rz1") == 0)) {
      argc--;
      argv++;
      target_z1 = source_z1 = atoi(argv[1]);
      argc--;
      argv++;
      ok = True;
    }
    if ((ok == False) && (strcmp(argv[1], "-Rz2") == 0)) {
      argc--;
      argv++;
      target_z2 = source_z2 = atoi(argv[1]);
      argc--;
      argv++;
      ok = True;
    }
    if ((ok == False) && (strcmp(argv[1], "-Tx1") == 0)) {
      argc--;
      argv++;
      target_x1 = atoi(argv[1]);
      argc--;
      argv++;
      ok = True;
    }
    if ((ok == False) && (strcmp(argv[1], "-Tx2") == 0)) {
      argc--;
      argv++;
      target_x2 = atoi(argv[1]);
      argc--;
      argv++;
      ok = True;
    }
    if ((ok == False) && (strcmp(argv[1], "-Ty1") == 0)) {
      argc--;
      argv++;
      target_y1 = atoi(argv[1]);
      argc--;
      argv++;
      ok = True;
    }
    if ((ok == False) && (strcmp(argv[1], "-Ty2") == 0)) {
      argc--;
      argv++;
      target_y2 = atoi(argv[1]);
      argc--;
      argv++;
      ok = True;
    }
    if ((ok == False) && (strcmp(argv[1], "-Tz1") == 0)) {
      argc--;
      argv++;
      target_z1 = atoi(argv[1]);
      argc--;
      argv++;
      ok = True;
    }
    if ((ok == False) && (strcmp(argv[1], "-Tz2") == 0)) {
      argc--;
      argv++;
      target_z2 = atoi(argv[1]);
      argc--;
      argv++;
      ok = True;
    }
    if ((ok == False) && (strcmp(argv[1], "-Sx1") == 0)) {
      argc--;
      argv++;
      source_x1 = atoi(argv[1]);
      argc--;
      argv++;
      ok = True;
    }
    if ((ok == False) && (strcmp(argv[1], "-Sx2") == 0)) {
      argc--;
      argv++;
      source_x2 = atoi(argv[1]);
      argc--;
      argv++;
      ok = True;
    }
    if ((ok == False) && (strcmp(argv[1], "-Sy1") == 0)) {
      argc--;
      argv++;
      source_y1 = atoi(argv[1]);
      argc--;
      argv++;
      ok = True;
    }
    if ((ok == False) && (strcmp(argv[1], "-Sy2") == 0)) {
      argc--;
      argv++;
      source_y2 = atoi(argv[1]);
      argc--;
      argv++;
      ok = True;
    }
    if ((ok == False) && (strcmp(argv[1], "-Sz1") == 0)) {
      argc--;
      argv++;
      source_z1 = atoi(argv[1]);
      argc--;
      argv++;
      ok = True;
    }
    if ((ok == False) && (strcmp(argv[1], "-Sz2") == 0)) {
      argc--;
      argv++;
      source_z2 = atoi(argv[1]);
      argc--;
      argv++;
      ok = True;
    }
    if ((ok == False) && (strcmp(argv[1], "-Tp") == 0)) {
      argc--;
      argv++;
      target_padding = atoi(argv[1]);
      argc--;
      argv++;
      ok = True;
    }
    if ((ok == False) && (strcmp(argv[1], "-Sp") == 0)) {
      argc--;
      argv++;
      source_padding = atoi(argv[1]);
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
    if ((ok == False) && (strcmp(argv[1], "-linear") == 0)) {
      argc--;
      argv++;
      interpolator = new irtkLinearInterpolateImageFunction<irtkGreyPixel>;
      ok = True;
    }
    if ((ok == False) && (strcmp(argv[1], "-bspline") == 0)) {
      argc--;
      argv++;
      interpolator = new irtkBSplineInterpolateImageFunction<irtkGreyPixel>;
      ok = True;
    }
    if ((ok == False) && (strcmp(argv[1], "-cspline") == 0)) {
      argc--;
      argv++;
      interpolator = new irtkCSplineInterpolateImageFunction<irtkGreyPixel>;
      ok = True;
    }
    if ((ok == False) && (strcmp(argv[1], "-sinc") == 0)) {
      argc--;
      argv++;
      interpolator = new irtkSincInterpolateImageFunction<irtkGreyPixel>;
      ok = True;
    }
    if (ok == False) {
      cerr << "Can not parse argument " << argv[1] << endl;
      usage();
    }
  }

  // Create default interpolator
  if (interpolator == NULL) {
    interpolator = new irtkNearestNeighborInterpolateImageFunction<irtkGreyPixel>;
  }

  // If there is no target image use copy of source image as target image
  if (target == NULL) {
    target = new irtkGreyImage(*source);
  }

  // Compute region of interest for target image
  if (target_x1 == -1) target_x1 = 0;
  if (target_y1 == -1) target_y1 = 0;
  if (target_z1 == -1) target_z1 = 0;
  if (target_x2 == -1) target_x2 = target->GetX();
  if (target_y2 == -1) target_y2 = target->GetY();
  if (target_z2 == -1) target_z2 = target->GetZ();

  // Compute region of interest for source image
  if (source_x1 == -1) source_x1 = 0;
  if (source_y1 == -1) source_y1 = 0;
  if (source_z1 == -1) source_z1 = 0;
  if (source_x2 == -1) source_x2 = source->GetX();
  if (source_y2 == -1) source_y2 = source->GetY();
  if (source_z2 == -1) source_z2 = source->GetZ();

  // If there is an region of interest for the target image, use it
  if ((target_x1 != 0) || (target_x2 != target->GetX()) ||
      (target_y1 != 0) || (target_y2 != target->GetY()) ||
      (target_z1 != 0) || (target_z2 != target->GetZ())) {
    *target = target->GetRegion(target_x1, target_y1, target_z1,
                                target_x2, target_y2, target_z2);
  }

  // If there is an region of interest for the source image, use it
  if ((source_x1 != 0) || (source_x2 != source->GetX()) ||
      (source_y1 != 0) || (source_y2 != source->GetY()) ||
      (source_z1 != 0) || (source_z2 != source->GetZ())) {
    *source = source->GetRegion(source_x1, source_y1, source_z1,
                                source_x2, source_y2, source_z2);
  }

  if (dof1_name != NULL) {
    // Read transformation
    transformation1 = irtkTransformation::New(dof1_name);
  } else {
    // Create identity transformation
    transformation1 = new irtkRigidTransformation;
  }

  if (dof2_name != NULL) {
    // Read transformation
    transformation2 = irtkTransformation::New(dof2_name);
  } else {
    // Create identity transformation
    transformation2 = new irtkRigidTransformation;
  }


  // Create image transformation
  irtkImageTransformation2<irtkGreyPixel> *imagetransformation =
    new irtkImageTransformation2<irtkGreyPixel>;

  imagetransformation->SetInput (source);
  imagetransformation->SetTransformation(transformation1, transformation2);
  imagetransformation->SetOutput(target);
  imagetransformation->PutTargetPaddingValue(target_padding);
  imagetransformation->PutSourcePaddingValue(source_padding);
  imagetransformation->PutInterpolator(interpolator);

  // Do inverse transformation if necessary
  if (invert == True) {
    imagetransformation->InvertOn();
  } else {
    imagetransformation->InvertOff();
  }

  // Transform image
  imagetransformation->Run();

  // Write the final transformation estimate
  target->Write(output_name);
}
