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

char *dofin_name  = NULL;
char *dofout_name = NULL;

irtkFreeFormTransformation *ffd;
irtkMultiLevelFreeFormTransformation *mffd;

void usage()
{
  cerr << "Usage: ffdsubdivide [dofin] [dofout] [options]" << endl;
  cerr << "-level <level>      MFFD level (default: current level)" << endl;
  cout << endl;
  exit(1);
}

int main(int argc, char **argv)
{
  int level, ok;

  // Check command line
  if (argc < 3) {
    usage();
  }

  dofin_name  = argv[1];
  argc--;
  argv++;
  dofout_name = argv[1];
  argc--;
  argv++;

  // Read transformation
  cout << "Reading transformation ... "; cout.flush();
  mffd = new irtkMultiLevelFreeFormTransformation;
  mffd->irtkTransformation::Read(dofin_name);
  cout << "done" << endl;

  // Default level is the last level
  level = mffd->NumberOfLevels()-1;

  // Parse remaining parameters
  while (argc > 1) {
    ok = false;
    if ((ok == false) && (strcmp(argv[1], "-level") == 0)) {
      argc--;
      argv++;
      level = atoi(argv[1]);
      argc--;
      argv++;
      ok = true;
    }
    if (ok == false) {
      cerr << "Can not parse argument " << argv[1] << endl;
      usage();
    }
  }

  // Extract current transformation level
  ffd = mffd->GetLocalTransformation(level);

  // Subdivide current transformation level
  cout << "Subdividing transformation level " << level << " ... ";
  cout.flush();
  ffd->Subdivide();
  cout << endl;

  // Write transformation
  mffd->irtkTransformation::Write(dofout_name);
}
