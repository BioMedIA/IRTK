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
char *dofin_name1 = NULL, *dofin_name2 = NULL;

// Output transformation
char *dofout_name = NULL;

void usage()
{
  cerr << "Usage: ffdsubtract [dofin1] [dofin2] [dofout]\n" << endl;
  exit(1);
}

int main(int argc, char **argv)
{
  int i, j, k, level;
  double x1, y1, z1, x2, y2, z2;

  // Check command line
  if (argc != 4) {
    usage();
  }

  // Parse file names
  dofin_name1 = argv[1];
  argc--;
  argv++;
  dofin_name2 = argv[1];
  argc--;
  argv++;
  dofout_name = argv[1];
  argc--;
  argv++;

  // Read transformation no. 1
  cout << "Reading transformation ... "; cout.flush();
  irtkMultiLevelFreeFormTransformation *mffd1 = new irtkMultiLevelFreeFormTransformation;
  mffd1->irtkTransformation::Read(dofin_name1);
  cout << "done" << endl;

  // Read transformation no. 2
  cout << "Reading transformation ... "; cout.flush();
  irtkMultiLevelFreeFormTransformation *mffd2 = new irtkMultiLevelFreeFormTransformation;
  mffd2->irtkTransformation::Read(dofin_name2);
  cout << "done" << endl;

  if (mffd1->NumberOfLevels() != mffd2->NumberOfLevels()) {
    cerr << "FFDs must have the same number of levels" << endl;
    exit(1);
  }

  for (level = 0; level < mffd1->NumberOfLevels(); level++) {

    // Extract current transformation level
    irtkFreeFormTransformation3D *ffd1 = dynamic_cast<irtkFreeFormTransformation3D *>(mffd1->GetLocalTransformation(level));

    if (ffd1 == NULL) {
      cerr << "Free-form transformation is not 3D" << endl;
      exit(1);
    }
    // Extract current transformation level
    irtkFreeFormTransformation3D *ffd2 = dynamic_cast<irtkFreeFormTransformation3D *>(mffd2->GetLocalTransformation(level));

    if (ffd2 == NULL) {
      cerr << "Free-form transformation is not 3D" << endl;
      exit(1);
    }

    // Make sure that all levels have the same number of control points
    if ((ffd1->GetX() != ffd2->GetX()) ||
        (ffd1->GetY() != ffd2->GetY()) ||
        (ffd1->GetZ() != ffd2->GetZ())) {
      cerr << "Number of control points on each level must be identical"
           << endl;
      exit(1);
    }

    // Subtract levels current level to the first level
    for (i = 0; i < ffd1->GetX(); i++) {
      for (j = 0; j < ffd1->GetY(); j++) {
        for (k = 0; k < ffd1->GetZ(); k++) {

          // Get control point value on first level
          ffd1->Get(i, j, k, x1, y1, z1);

          // Get control point value on current level
          ffd2->Get(i, j, k, x2, y2, z2);

          // Add levels and store result in first level
          ffd1->Put(i, j, k, x1-x2, y1-y2, z1-z2);

        }
      }
    }
  }

  // Write transformation
  mffd1->irtkTransformation::Write(dofout_name);
}
