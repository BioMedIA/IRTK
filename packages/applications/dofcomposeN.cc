/*=========================================================================

  Library   : Image Registration Toolkit (IRTK)
  Module    : $Id$
  Copyright : Imperial College, Department of Computing
              Visual Information Processing (VIP), 2011 onwards
  Date      : $Date$
  Version   : $Revision$
  Changes   : $Author$

=========================================================================*/

#include <irtkImage.h>
#include <irtkImageFunction.h>
#include <irtkTransformation.h>
#include <irtkRegistration.h>

// Default filenames
char *dofOut_name = NULL;
char **dof_name   = NULL;

#define MAX_DOFS 50

void usage()
{
  cerr << " Usage: dofcomposeN [image] [dofOut] <options>\n" << endl;
  cerr << " where <options> is one or more of the following: \n" << endl;
  cerr << " " << endl;
  cerr << " <-dofin file>      A transformation file which should be included. Multiple transformations " << endl;
  cerr << "                    can be given by repeatedly using this flag. They are processed in order " << endl;
  cerr << "                    in which they are passed. " << endl;
  cerr << " <-dofin_i file>    A transformation whose inverse should be included." << endl;
  cerr << " " << endl;
  cerr << " " << endl;
  cerr << " E.g." << endl;
  cerr << " " << endl;
  cerr << "   dofcomposeN out.dof.gz -dofin a.dof -dofin_i b.dof -dofin c.dof" << endl;
  cerr << " " << endl;
  cerr << " returns the composed transformation c b^-1 a(x) applied to each location x in the input image" << endl;
  cerr << " " << endl;
  exit(1);
}

int main(int argc, char **argv)
{
  irtkTransformation **transformation = NULL;
  int i, j, k, n, noOfDofs;
  double x1, y1, z1, x2, y2, z2;
  bool ok;

  // Check command line
  if (argc < 3) {
    usage();
  }

  // Read input image
  irtkImage *image = irtkImage::New(argv[1]);
  argc--;
  argv++;

  // Parse dofOut_name
  dofOut_name = argv[1];
  argc--;
  argv++;

  // Fix number of dofs
  noOfDofs = 0;

  dof_name = new char*[MAX_DOFS];
  bool *invert = new bool[MAX_DOFS];

  while (argc > 1) {
    ok = false;
    if ((ok == false) && (strcmp(argv[1], "-dofin") == 0)) {
      argc--;
      argv++;
      dof_name[noOfDofs] = argv[1];
      invert[noOfDofs]   = false;
      noOfDofs++;
      argc--;
      argv++;
      ok = true;
    }
    if ((ok == false) && (strcmp(argv[1], "-dofin_i") == 0)) {
      argc--;
      argv++;
      dof_name[noOfDofs] = argv[1];
      invert[noOfDofs]   = true;
      noOfDofs++;
      argc--;
      argv++;
      ok = true;
    }

    if (ok == false) {
      cerr << "Can not parse argument " << argv[1] << endl;
      usage();
    }
  }

  if (noOfDofs == 0) {
    noOfDofs = 1;
    transformation = new irtkTransformation*[noOfDofs];
    transformation[0] = new irtkRigidTransformation;
  } else {
    transformation = new irtkTransformation*[noOfDofs];
    for (n = 0; n < noOfDofs; n++) {
      transformation[n] = irtkTransformation::New(dof_name[n]);
    }
  }

  irtkFluidFreeFormTransformation *output = new irtkFluidFreeFormTransformation;

  for (n = 0; n < noOfDofs; n++) {
    // Create displacement field with a displacement vector for each voxel
    irtkLinearFreeFormTransformation *current = new irtkLinearFreeFormTransformation(*image, image->GetXSize(), image->GetYSize(), image->GetZSize());
    for (k = 0; k < image->GetZ(); k++) {
      for (j = 0; j < image->GetY(); j++) {
        for (i = 0; i < image->GetX(); i++) {
          x1 = i;
          y1 = j;
          z1 = k;
          image->ImageToWorld(x1, y1, z1);
          x2 = x1;
          y2 = y1;
          z2 = z1;
          if (invert[n] == true) {
          	transformation[n]->Inverse(x2, y2, z2);
          } else {
          	transformation[n]->Transform(x2, y2, z2);
          }
          x2 -= x1;
          y2 -= y1;
          z2 -= z1;
          current->Put(i, j, k, x2, y2, z2);
        }
      }
    }
    output->PushLocalTransformation(current);
  }

  // Write composed transformation
  output->irtkTransformation::Write(dofOut_name);
}
