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

void usage()
{
  cerr << "Usage: ffdinfo [dofin]\n" << endl;
  exit(1);
}

int main(int argc, char **argv)
{
  int i, j, k, n, active, passive;
  double dx, dy, dz, xaxis[3], yaxis[3], zaxis[3];

  // Check command line
  if (argc != 2) {
    usage();
  }

  // Read transformation
  irtkTransformation *transform = irtkTransformation::New(argv[1]);
  irtkMultiLevelFreeFormTransformation *mffd = dynamic_cast<irtkMultiLevelFreeFormTransformation *>(transform);

  for (n = 0; n < mffd->NumberOfLevels(); n++) {

    // Extract current transformation level
    irtkFreeFormTransformation3D *ffd = dynamic_cast<irtkFreeFormTransformation3D *>(mffd->GetLocalTransformation(n));

    if (ffd == NULL) {
      cerr << "Free-form transformation is not 3D" << endl;
      exit(1);
    }

    ffd->GetOrientation(xaxis, yaxis, zaxis);
    ffd->GetSpacing(dx, dy, dz);

    // Print information about local transformation
    cout << "Local transformation no. " << n+1 << ":" << endl;
    cout << "Control points: \t" << ffd->GetX() << " x " << ffd->GetY()
         << " x " << ffd->GetZ() << endl;
    cout << "Orientation: \t\t" << xaxis[0] << " " << xaxis[1] << " "
         << xaxis[2] << endl;
    cout << "\t\t\t"  <<yaxis[0] << " " << yaxis[1] << " " << yaxis[2] << endl;
    cout << "\t\t\t" <<zaxis[0] << " " << zaxis[1] << " " << zaxis[2] << endl;
    cout << "Spacing: \t\t" << dx << " " << dy << " " << dx << endl;

    active  = 0;
    passive = 0;
    for (i = 0; i < ffd->GetX(); i++) {
      for (j = 0; j < ffd->GetY(); j++) {
        for (k = 0; k < ffd->GetZ(); k++) {
          _Status sx, sy, sz;
          ffd->GetStatus(i, j, k, sx, sy, sz);
          if (sx == _Active) {
            active++;
          } else {
            passive++;
          }
          if (sy == _Active) {
            active++;
          } else {
            passive++;
          }
          if (sz == _Active) {
            active++;
          } else {
            passive++;
          }
        }
      }
    }
    cout << "Active control points: \t " << active << " ("
         << active*100.0/(active+passive) << "%)" << endl;
    cout << "Passive control points: \t " << passive << " ("
         << passive*100.0/(active+passive) << "%)" << endl;
  }
}
