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
char *input_name = NULL, *output_name, *dof_name  = NULL;

typedef enum { LocalJacobian, GlobalJacobian, RelativeJacobian, TotalJacobian } JacobianMode;

void usage()
{
  cerr << "Usage: jacobian [input] [output] [ffd]\n" << endl;
  cerr << "where <options> is one or more of the following:\n" << endl;
  cerr << "<-total>                 Total jacobian (default)" << endl;
  cerr << "<-local>                 Local jacobian only " << endl;
  cerr << "<-global>                Global jacobian only " << endl;
  cerr << "<-relative>              Local jacobian divided by global Jacobian" << endl;
  cerr << "<-padding value>         Padding value" << endl;
  exit(1);
}

int main(int argc, char **argv)
{
  int i, j, k, n, m, ok, fluid, padding;
  double x, y, z, jac;
  JacobianMode jac_mode;

  // Check command line
  if (argc < 4) {
    usage();
  }

  // Parse image
  input_name  = argv[1];
  argc--;
  argv++;
  output_name = argv[1];
  argc--;
  argv++;
  dof_name = argv[1];
  argc--;
  argv++;

  // Initialize padding value
  padding = MIN_GREY;

  // Initialize jacobian mode
  jac_mode = TotalJacobian;
  fluid = false;

  while (argc > 1) {
    ok = false;
    if ((ok == false) && (strcmp(argv[1], "-padding") == 0)) {
      argc--;
      argv++;
      padding = atoi(argv[1]);
      argc--;
      argv++;
      ok = true;
    }
    if ((ok == false) && (strcmp(argv[1], "-total") == 0)) {
      argc--;
      argv++;
      jac_mode = TotalJacobian;
      ok = true;
    }
    if ((ok == false) && (strcmp(argv[1], "-local") == 0)) {
      argc--;
      argv++;
      jac_mode = LocalJacobian;
      ok = true;
    }
    if ((ok == false) && (strcmp(argv[1], "-global") == 0)) {
      argc--;
      argv++;
      jac_mode = GlobalJacobian;
      ok = true;
    }
    if ((ok == false) && (strcmp(argv[1], "-relative") == 0)) {
      argc--;
      argv++;
      jac_mode = RelativeJacobian;
      ok = true;
    }
    if ((ok == false) && (strcmp(argv[1], "-fluid") == 0)) {
      argc--;
      argv++;
      fluid = true;
      ok = true;
    }
    if (ok == false) {
      cerr << "Can not parse argument " << argv[1] << endl;
      usage();
    }
  }

  // Read image
  cout << "Reading image ... "; cout.flush();
  irtkGreyImage *image = new irtkGreyImage(input_name);
  cout << "done" << endl;

  // Read transformation
  irtkMultiLevelFreeFormTransformation *mffd;
  if (fluid == true) {
    mffd = new irtkFluidFreeFormTransformation;
  } else {
    mffd = new irtkMultiLevelFreeFormTransformation;
  }
  mffd->irtkTransformation::Read(dof_name);

  m = 0;
  n = 0;
  jac = 1;
  for (k = 0; k < image->GetZ(); k++) {
    for (j = 0; j < image->GetY(); j++) {
      for (i = 0; i < image->GetX(); i++) {
        if (image->Get(i, j, k) > padding) {
          x = i;
          y = j;
          z = k;
          image->ImageToWorld(x, y, z);
          switch (jac_mode) {
          case TotalJacobian:
            jac = mffd->irtkTransformation::Jacobian(x, y, z);
            break;
          case LocalJacobian:
            jac = mffd->irtkTransformation::LocalJacobian(x, y, z);
            break;
          case GlobalJacobian:
            jac = mffd->irtkTransformation::GlobalJacobian(x, y, z);
            break;
          case RelativeJacobian:
            jac = mffd->irtkTransformation::LocalJacobian(x, y, z) / mffd->irtkTransformation::GlobalJacobian(x, y, z);
            break;
          default:
            break;
          }
          m++;
          if (jac < 0) n++;
          image->Put(i, j, k, round(100*jac));
        }
      }
    }
  }

  if (n > 0) {
    cout << "Number of voxels with negative jacobian = " << n << " (" << (double)n/(double)m*100.0 << "%)" << endl;
  }

  // Write the final transformation estimate
  image->Write(output_name);
}
