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

// Input transformation
char *dofin_name  = NULL;

// Output transformation
char *dofout_name = NULL;

void usage()
{
  cerr << "Usage: ffd2D [dofin] [dofout] <-level n>" << endl;
  exit(1);
}

int main(int argc, char **argv)
{
  int x, y, z, ok, level;

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
  irtkMultiLevelFreeFormTransformation *mffd = new irtkMultiLevelFreeFormTransformation;
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
      cerr << "Can't parse argument " << argv[1] << endl;
      usage();
    }
  }

  // Extract current transformation level
  irtkFreeFormTransformation3D *ffd = dynamic_cast<irtkFreeFormTransformation3D *>(mffd->GetLocalTransformation(level));

  if (ffd == NULL) {
    cerr << "Free-form transformation is not 3D" << endl;
    exit(1);
  }

  // Edit current transformation level
  cout << "Editing transformation level " << level << " ... ";
  cout.flush();

  // Get lattice dimensions and loop over lattice
  for (x = 0; x < ffd->GetX(); x++) {
    for (y = 0; y < ffd->GetY(); y++) {
      for (z = 0; z < ffd->GetZ(); z++) {
        _Status sx, sy, sz;
        ffd->GetStatus(x, y, z, sx, sy, sz);
        ffd->PutStatus(x, y, z, sx, sy, _Passive);
      }
    }
  }
  cout << "done" << endl;

  // Write transformation
  mffd->irtkTransformation::Write(dofout_name);
}

