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

// Image
char *image_name = NULL;

irtkFreeFormTransformation3D *ffd;
irtkMultiLevelFreeFormTransformation *mffd;

void usage()
{
  cerr << "Usage: ffdpadding [image] [dofin] [dofout] [padding] [options]\n"
       << endl;
  cerr << "where [options] are one or more of the following:\n" << endl;
  cerr << "-level  <level>           MFFD level (default: current level)"
       << endl;
  cerr << "-status <active/passive>  Control point status (default: passive)"
       << endl;
  exit(1);
}

int main(int argc, char **argv)
{
  irtkGreyImage image;
  int ok, level, padding, index, i, j, k, x, y, z, x1, y1, z1, x2, y2, z2;

  // Check command line
  if (argc < 5) {
    usage();
  }

  // Parse file names
  image_name = argv[1];
  argc--;
  argv++;
  dofin_name  = argv[1];
  argc--;
  argv++;
  dofout_name = argv[1];
  argc--;
  argv++;
  padding = atoi(argv[1]);
  argc--;
  argv++;

  // Read image
  cout << "Reading image ... "; cout.flush();
  image.Read(image_name);
  cout << "done" << endl;

  // Read transformation
  cout << "Reading transformation ... "; cout.flush();
  mffd = new irtkMultiLevelFreeFormTransformation;
  mffd->irtkTransformation::Read(dofin_name);
  cout << "done" << endl;

  // Default level is the last level
  level = mffd->NumberOfLevels()-1;

  // Default status is passive
  _Status status = _Passive;

  // Parse remaining arguments
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
    if ((ok == false) && (strcmp(argv[1], "-status") == 0)) {
      argc--;
      argv++;
      if (strcmp(argv[1], "active") == 0) {
        status = _Active;
      } else {
        if (strcmp(argv[1], "passive") == 0) {
          status = _Passive;
        } else {
          usage();
        }
      }
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
  ffd = dynamic_cast<irtkFreeFormTransformation3D *>(mffd->GetLocalTransformation(level));

  if (ffd == NULL) {
    cerr << "Free-form transformation is not 3D" << endl;
    exit(1);
  }

  // Edit current transformation level
  cout << "Editing transformation level " << level << " ... ";
  cout.flush();

  // Calculate number of active and passive control points
  for (i = 0; i < ffd->GetX(); i++) {
    for (j = 0; j < ffd->GetY(); j++) {
      for (k = 0; k < ffd->GetZ(); k++) {
        // Convert control points to index
        index = ffd->LatticeToIndex(i, j, k);

        // Calculate bounding box of control point in voxels
        ffd->BoundingBoxImage(&image, index, x1, y1, z1, x2, y2, z2);

        ok = false;
        for (z = z1; z <= z2; z++) {
          for (y = y1; y <= y2; y++) {
            for (x = x1; x <= x2; x++) {
              if (image(x, y, z) > padding) {
                ok = true;
              }
            }
          }
        }
        if (ok == false) {
          ffd->PutStatusCP(i, j, k, status, status, status);
        }
      }
    }
  }
  cout << "done" << endl;

  // Write transformation
  mffd->irtkTransformation::Write(dofout_name);
}
