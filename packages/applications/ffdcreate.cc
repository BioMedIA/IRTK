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

char *dofout_name = NULL;
char *dofin_name  = NULL;

void usage()
{
  cerr << "Usage: ffdcreate [dofout] <options>\n" << endl;
  cerr << "where <options> is one or more of the following:\n" << endl;
  cerr << "<-dofin file>   Rigid/affine transformation" << endl;
  exit(1);
}

int main(int argc, char **argv)
{
  int ok, fluid;
  irtkMultiLevelFreeFormTransformation *mffd;

  // Check command line
  if (argc < 2) {
    usage();
  }

  // Parse file names
  dofout_name = argv[1];
  argc--;
  argv++;

  // No fluid transformation
  fluid = False;

  // Parse remaining parameters
  while (argc > 1) {
    ok = False;
    if ((ok == False) && (strcmp(argv[1], "-dofin") == 0)) {
      argc--;
      argv++;
      dofin_name = argv[1];
      argc--;
      argv++;
      ok = True;
    }
    if (ok == False) {
      cerr << "Can not parse argument " << argv[1] << endl;
      usage();
    }
  }

  // Check if we want to create a fluid transformation
  if (fluid == True) {
    if (dofin_name != NULL) {
      mffd = NULL;
      cerr << "ffdcreate: Not yet implemented" << endl;
    } else {
      mffd = new irtkFluidFreeFormTransformation;
    }
  } else {
    if (dofin_name != NULL) {
      irtkTransformation *transform = irtkTransformation::New(dofin_name);
      if (strcmp(transform->NameOfClass(), "irtkRigidTransformation") == 0) {
        mffd = new irtkMultiLevelFreeFormTransformation(*((irtkRigidTransformation *)transform));
      } else {
        if (strcmp(transform->NameOfClass(), "irtkAffineTransformation") == 0) {
          mffd = new irtkMultiLevelFreeFormTransformation(*((irtkAffineTransformation *)transform));
        } else {
          cerr << "Input transformation is not of type rigid, affine " << endl;
          cerr << "or multi-level free form deformation" << endl;
          exit(1);
        }
      }
      delete transform;
    } else {
      mffd = new irtkMultiLevelFreeFormTransformation;
    }
  }

  // Write file
  mffd->irtkTransformation::Write(dofout_name);
}
