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
  cerr << "Usage: ffdmovie [dofin] [dofout] <options>\n" << endl;
  cerr << "where <options> is one or more of the following:\n" << endl;
  cerr << "<-#>   No. of transformation to interpolate (default 20)" << endl;
  exit(1);
}

int main(int argc, char **argv)
{
  int i, j, n, nFrames, ok;

  // Check command line
  if (argc < 2) {
    usage();
  }

  // Parse file names
  dofin_name = argv[1];
  argc--;
  argv++;
  dofout_name = argv[1];
  argc--;
  argv++;

  // Parse remaining parameters
  nFrames = 20;
  while (argc > 1) {
    ok = false;
    if ((ok == false) && (strcmp(argv[1], "-#") == 0)) {
      argc--;
      argv++;
      nFrames = atoi(argv[1]);
      argc--;
      argv++;
      ok = true;
    }
    if (ok == false) {
      cerr << "Can not parse argument " << argv[1] << endl;
      usage();
    }
  }

  // Create muitlevel free-form transformation
  irtkMultiLevelFreeFormTransformation *mffd = new irtkMultiLevelFreeFormTransformation;

  // Loop over levels
  for (n = 0; n < nFrames; n++) {

    // Read transformation
    mffd->irtkTransformation::Read(dofin_name);

    // Loop over levels
    for (i = 0; i < mffd->NumberOfLevels(); i++) {
      irtkFreeFormTransformation *ffd = mffd->GetLocalTransformation(i);

      // Loop over DOFs
      for (j = 0; j < ffd->NumberOfDOFs(); j++) {
        ffd->Put(j, n / (double)(nFrames - 1) * ffd->Get(j));
      }
    }

    // Write file
    char buffer[1024];
    sprintf(buffer, "%s_%d.dof", dofout_name, n);
    mffd->irtkTransformation::Write(buffer);
  }
}
