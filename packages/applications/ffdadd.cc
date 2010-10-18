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

char *target_name = NULL;
char *dofin_name  = NULL;
char *dofout_name = NULL;

irtkGreyImage *target;
irtkMultiLevelFreeFormTransformation *mffd;

void usage()
{
  cerr << "Usage: ffdadd [target] [dofin] [dofout] [options]" << endl;
  cerr << "where the following options are supported:\n" << endl;
  cerr << "<-Tx1 value>  Region of interest in target image" << endl;
  cerr << "<-Ty1 value>  Region of interest in target image" << endl;
  cerr << "<-Tz1 value>  Region of interest in target image" << endl;
  cerr << "<-Tt1 value>  Region of interest in target image" << endl;
  cerr << "<-Tx2 value>  Region of interest in target image" << endl;
  cerr << "<-Ty2 value>  Region of interest in target image" << endl;
  cerr << "<-Tz2 value>  Region of interest in target image" << endl;
  cerr << "<-Tt2 value>  Region of interest in target image" << endl;
  cerr << "<-dx  value>  Control point spacing for x" << endl;
  cerr << "<-dy  value>  Control point spacing for y" << endl;
  cerr << "<-dz  value>  Control point spacing for z" << endl;
  cerr << "<-dt  value>  Control point spacing for t" << endl;
  cerr << "<-ds  value>  Control point spacing for x, y, z and t" << endl;
  cerr << "<-3D>         Generate 3D FFD (default)" << endl;
  cerr << "<-4D>         Generate 4D FFD" << endl;
  cout << endl;
  exit(1);
}

int main(int argc, char **argv)
{
  double dx, dy, dz, dt;
  int i1, j1, k1, l1, i2, j2, k2, l2, ok, dilate, fluid, ffd4D;

  // Check command line
  if (argc < 4) {
    usage();
  }

  // Parse file names
  target_name = argv[1];
  argc--;
  argv++;
  dofin_name  = argv[1];
  argc--;
  argv++;
  dofout_name = argv[1];
  argc--;
  argv++;

  // Read image image
  cout << "Reading image ... "; cout.flush();
  target = new irtkGreyImage;
  target->Read(target_name);
  cout << "done" << endl;

  // Fix ROI
  i1 = 0;
  j1 = 0;
  k1 = 0;
  l1 = 0;
  i2 = target->GetX();
  j2 = target->GetY();
  k2 = target->GetZ();
  l2 = target->GetT();

  // Fix spacing
  dx = 20;
  dy = 20;
  dz = 20;
  dt = 20;

  // Don't dilate bounding box
  dilate = false;

  // Don't dilate bounding box
  ffd4D = false;

  // No fluid
  fluid = false;

  // Parse remaining parameters
  while (argc > 1) {
    ok = false;
    if ((ok == false) && (strcmp(argv[1], "-Tx1") == 0)) {
      argc--;
      argv++;
      i1 = atoi(argv[1]);
      argc--;
      argv++;
      ok = true;
    }
    if ((ok == false) && (strcmp(argv[1], "-Tx2") == 0)) {
      argc--;
      argv++;
      i2 = atoi(argv[1]);
      argc--;
      argv++;
      ok = true;
    }
    if ((ok == false) && (strcmp(argv[1], "-Ty1") == 0)) {
      argc--;
      argv++;
      j1 = atoi(argv[1]);
      argc--;
      argv++;
      ok = true;
    }
    if ((ok == false) && (strcmp(argv[1], "-Ty2") == 0)) {
      argc--;
      argv++;
      j2 = atoi(argv[1]);
      argc--;
      argv++;
      ok = true;
    }
    if ((ok == false) && (strcmp(argv[1], "-Tz1") == 0)) {
      argc--;
      argv++;
      k1 = atoi(argv[1]);
      argc--;
      argv++;
      ok = true;
    }
    if ((ok == false) && (strcmp(argv[1], "-Tz2") == 0)) {
      argc--;
      argv++;
      k2 = atoi(argv[1]);
      argc--;
      argv++;
      ok = true;
    }
    if ((ok == false) && (strcmp(argv[1], "-Tt1") == 0)) {
      argc--;
      argv++;
      l1 = atoi(argv[1]);
      argc--;
      argv++;
      ok = true;
    }
    if ((ok == false) && (strcmp(argv[1], "-Tt2") == 0)) {
      argc--;
      argv++;
      l2 = atoi(argv[1]);
      argc--;
      argv++;
      ok = true;
    }
    if ((ok == false) && (strcmp(argv[1], "-dx") == 0)) {
      argc--;
      argv++;
      dx = atof(argv[1]);
      argc--;
      argv++;
      ok = true;
    }
    if ((ok == false) && (strcmp(argv[1], "-dy") == 0)) {
      argc--;
      argv++;
      dy = atof(argv[1]);
      argc--;
      argv++;
      ok = true;
    }
    if ((ok == false) && (strcmp(argv[1], "-dz") == 0)) {
      argc--;
      argv++;
      dz = atof(argv[1]);
      argc--;
      argv++;
      ok = true;
    }
    if ((ok == false) && (strcmp(argv[1], "-dt") == 0)) {
      argc--;
      argv++;
      dt = atof(argv[1]);
      argc--;
      argv++;
      ok = true;
    }
    if ((ok == false) && (strcmp(argv[1], "-ds") == 0)) {
      argc--;
      argv++;
      dx = atof(argv[1]);
      dy = atof(argv[1]);
      dz = atof(argv[1]);
      dt = atof(argv[1]);
      argc--;
      argv++;
      ok = true;
    }
    if ((ok == false) && (strcmp(argv[1], "-fluid") == 0)) {
      argc--;
      argv++;
      fluid = true;
      ok = true;
    }
    if ((ok == false) && (strcmp(argv[1], "-3D") == 0)) {
      argc--;
      argv++;
      ffd4D = false;
      ok = true;
    }
    if ((ok == false) && (strcmp(argv[1], "-4D") == 0)) {
      argc--;
      argv++;
      ffd4D = true;
      ok = true;
    }
    if (ok == false) {
      cerr << "Can not parse argument " << argv[1] << endl;
      usage();
    }
  }

  // Read transformation
  if (fluid == true) {
    mffd = new irtkFluidFreeFormTransformation;
  } else {
    mffd = new irtkMultiLevelFreeFormTransformation;
  }
  mffd->irtkTransformation::Read(dofin_name);

  // If there is an region of interest, use it
  if ((i1 != 0) || (i2 != target->GetX()) ||
      (j1 != 0) || (j2 != target->GetY()) ||
      (k1 != 0) || (k2 != target->GetZ()) ||
      (l1 != 0) || (l2 != target->GetT())) {
    *target = target->GetRegion(i1, j1, k1, l1, i2, j2, k2, l2);
  }

  if (ffd4D == false) {
    // Create free-form transformation
    irtkBSplineFreeFormTransformation *affd = new irtkBSplineFreeFormTransformation(*target, dx, dy, dz);

    // Add and write file
    mffd->PushLocalTransformation(affd);
  } else {
    // Create free-form transformation
    irtkBSplineFreeFormTransformation4D *affd = new irtkBSplineFreeFormTransformation4D(*target, dx, dy, dz, dt);

    // Add and write file
    mffd->PushLocalTransformation(affd);
  }

  mffd->irtkTransformation::Write(dofout_name);
}
