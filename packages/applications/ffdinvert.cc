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

#include <irtkTransformation.h>

char *dofin_name, *dofout_name;

void usage()
{
  cerr << "Usage: ffdinvert [dofin] [dofout]\n";
  exit(1);
}

int main(int argc, char **argv)
{
  int i, j, k, index;
  double x1, y1, z1, x2, y2, z2, error, max_error, rms_error;

  // Check command line
  if (argc != 3) {
    usage();
  }

  dofin_name  = argv[1];
  argv++;
  argc--;
  dofout_name = argv[1];
  argv++;
  argc--;

  // Read transformation
  irtkMultiLevelFreeFormTransformation *mffd = new irtkMultiLevelFreeFormTransformation;
  mffd->irtkTransformation::Read(dofin_name);

  // Check number of levels
  if (mffd->NumberOfLevels() > 1) {
    cerr << "ffdinvert: Inverting of FFDs with more than one level currently not implemented" << endl;
    exit(1);
  }

  // Invert global transformation
  mffd->Invert();

  // Extract FFD
  irtkBSplineFreeFormTransformation3D *affd1 = dynamic_cast<irtkBSplineFreeFormTransformation3D *>(mffd->PopLocalTransformation());
  if (affd1 == NULL) {
    cerr << "Local transformation is not of type irtkBSplineFreeFormTransformation3D" << endl;
    exit(1);
  }

  // Copy FFD
  irtkBSplineFreeFormTransformation3D *affd2 = new irtkBSplineFreeFormTransformation3D(*affd1);
  for (k = 0; k < affd2->GetZ(); k++) {
    for (j = 0; j < affd2->GetY(); j++) {
      for (i = 0; i < affd2->GetX(); i++) {
        affd2->Put(i, j, k, 0, 0, 0);
      }
    }
  }

  // Allocate some memory
  double *dx = new double[affd1->GetX()*affd1->GetY()*affd1->GetZ()];
  double *dy = new double[affd1->GetX()*affd1->GetY()*affd1->GetZ()];
  double *dz = new double[affd1->GetX()*affd1->GetY()*affd1->GetZ()];

  index = 0;
  error = 0;
  max_error = 0;
  rms_error = 0;
  for (k = 0; k < affd1->GetZ(); k++) {
    for (j = 0; j < affd1->GetY(); j++) {
      for (i = 0; i < affd1->GetX(); i++) {
        x1 = i;
        y1 = j;
        z1 = k;
        affd1->LatticeToWorld(x1, y1, z1);
        x2 = x1;
        y2 = y1;
        z2 = z1;
        affd1->Inverse(x2, y2, z2, 0);
        dx[index] = x2 - x1;
        dy[index] = y2 - y1;
        dz[index] = z2 - z1;
        index++;
      }
    }
  }

  // Approximate transformations
  affd2->Interpolate(dx, dy, dz);

  index = 0;
  error = 0;
  max_error = 0;
  rms_error = 0;
  for (k = 2; k < affd1->GetZ()-2; k++) {
    for (j = 2; j < affd1->GetY()-2; j++) {
      for (i = 2; i < affd1->GetX()-2; i++) {
        x1 = i;
        y1 = j;
        z1 = k;
        affd1->LatticeToWorld(x1, y1, z1);
        x2 = x1;
        y2 = y1;
        z2 = z1;
        affd2->Transform(x2, y2, z2);
        affd1->Transform(x2, y2, z2);
        error = sqrt((x1-x2)*(x1-x2)+(y1-y2)*(y1-y2)+(z1-z2)*(z1-z2));
        rms_error += error;
        if (error > max_error) max_error = error;
        index++;
      }
    }
  }
  rms_error /= double(affd1->GetX()*affd1->GetY()*affd1->GetZ());
  cerr << "RMS error for inverse FFD is = " << rms_error << endl;
  cerr << "Max error for inverse FFD is = " << max_error << endl;

  delete []dx;
  delete []dy;
  delete []dz;

  // Add transformation back
  mffd->PushLocalTransformation(affd2);

  // Write transformation collection
  mffd->irtkTransformation::Write(dofout_name);
}
