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
  cerr << "Usage: ffdcollapse [dofin] [dofout] [options]\n" << endl;
  cerr << "where the following options are supported:\n" << endl;
  cerr << "<-subdivide>  Subdivide levels before collapsing" << endl;
  exit(1);
}

int main(int argc, char **argv)
{
  int i, j, k, ok, subdivide, level;
  double x1, y1, z1, x2, y2, z2;

  // Check command line
  if (argc != 3) {
    usage();
  }

  // Parse file names
  dofin_name  = argv[1];
  argc--;
  argv++;
  dofout_name = argv[1];
  argc--;
  argv++;

  // Parse remaining parameters
  while (argc > 1) {
    ok = false;
    if ((ok == false) && (strcmp(argv[1], "-subdivide") == 0)) {
      argc--;
      argv++;
      subdivide = true;
      ok = true;
    }
    if (ok == false) {
      cerr << "Can not parse argument " << argv[1] << endl;
      usage();
    }
  }

  // Read transformation
  cout << "Reading transformation ... "; cout.flush();
  irtkMultiLevelFreeFormTransformation *mffd = new irtkMultiLevelFreeFormTransformation;
  mffd->irtkTransformation::Read(dofin_name);
  cout << "done" << endl;

  // Subdivide if useful
  if (subdivide == true) {
    for (level = 0; level < mffd->NumberOfLevels(); level++) {
      irtkFreeFormTransformation3D *ffd = dynamic_cast<irtkFreeFormTransformation3D *>(mffd->GetLocalTransformation(level));
    	for (i = level; i < mffd->NumberOfLevels()-1; i++) ffd->Subdivide();
    }
  }

  // Default level is the last level
  level = mffd->NumberOfLevels()-1;

  // Extract first transformation level
  irtkFreeFormTransformation3D *ffd1 = dynamic_cast<irtkFreeFormTransformation3D *>(mffd->GetLocalTransformation(0));

  if (ffd1 == NULL) {
    cerr << "Free-form transformation is not 3D" << endl;
    exit(1);
  }

  for (level = 1; level < mffd->NumberOfLevels(); level++) {

    // Extract current transformation level
    irtkFreeFormTransformation3D *ffd2 = dynamic_cast<irtkFreeFormTransformation3D *>(mffd->GetLocalTransformation(level));

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

    // Add current level to the first level
    for (i = 0; i < ffd1->GetX(); i++) {
      for (j = 0; j < ffd1->GetY(); j++) {
        for (k = 0; k < ffd1->GetZ(); k++) {

          // Get control point value on first level
          ffd1->Get(i, j, k, x1, y1, z1);

          // Get control point value on current level
          ffd2->Get(i, j, k, x2, y2, z2);

          // Add levels and store result in first level
          ffd1->Put(i, j, k, x1+x2, y1+y2, z1+z2);

        }
      }
    }
  }

  // Now remove all levels starting from the last level, except for first level
  while (mffd->NumberOfLevels() > 1) {
    irtkFreeFormTransformation *ffd = mffd->PopLocalTransformation();
    delete ffd;
  }

  // Write transformation
  mffd->irtkTransformation::Write(dofout_name);
}
