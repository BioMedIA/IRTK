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
  cerr << "Usage: ffdcompose T_1 T_2 T_out" << endl;
  cerr << "ffdcompose approximates the composition of two MFFDs T_1 and T_2 using a single level FFD T such that T ~= T_2 o T_1." << endl;
  exit(1);
}

int main(int argc, char **argv)
{
  if(argc != 4){
    usage();
  }

  char *dof1_name = argv[1];
  argc--;
  argv++;
  char *dof2_name = argv[1];
  argc--;
  argv++;
  char *dofout_name = argv[1];
  argc--;
  argv++;

  // Print out warning
  cout << "WARNING: The current implementation will ignore the global affine transformation component!" << endl;

  // Read and check the input transformations
  irtkMultiLevelFreeFormTransformation *mffd1 = dynamic_cast<irtkMultiLevelFreeFormTransformation *>(irtkTransformation::New(dof1_name));
  if(mffd1 == NULL){
    cerr << "ERROR: Transformation T_1 is not of type irtkMultiLevelFreeFormTransformation!" << endl;
    exit(1);
  }

  irtkMultiLevelFreeFormTransformation *mffd2 = dynamic_cast<irtkMultiLevelFreeFormTransformation *>(irtkTransformation::New(dof2_name));
  if(mffd2 == NULL){
    cerr << "ERROR: Transformation T_2 is not of type irtkMultiLevelFreeFormTransformation!" << endl;
    exit(1);
  }

  // The first FFD in T1
  irtkBSplineFreeFormTransformation3D *ffd1 = dynamic_cast<irtkBSplineFreeFormTransformation3D *>(mffd1->GetLocalTransformation(0));
  if(ffd1 == NULL){
    cerr << "ERROR: Transformation T_1 does not contain a valid transformation of type irtkBSplineFreeFormTransformation3D!" << endl;
    exit(1);
  }

  // Create a BSpline FFD
  irtkBSplineFreeFormTransformation3D *ffd = new irtkBSplineFreeFormTransformation3D(*ffd1);
  int nx = ffd->GetX();
  int ny = ffd->GetY();
  int nz = ffd->GetZ();

  for(int k=0; k<nz; k++){
    for(int j=0; j<ny; j++){
      for(int i=0; i<nx; i++){
	ffd->Put(i, j, k, 0, 0, 0);
      }
    }
  }

  // The deformation vectors at the control points for the composed deformation
  double *dx = new double[nx * ny * nz];
  double *dy = new double[nx * ny * nz];
  double *dz = new double[nx * ny * nz];
  int index = 0;

  for(int k=0; k<nz; k++){
    for(int j=0; j<ny; j++){
      for(int i=0; i<nx; i++){
        double x1 = i;
        double y1 = j;
        double z1 = k;
        ffd->LatticeToWorld(x1, y1, z1);
        double x2 = x1;
        double y2 = y1;
        double z2 = z1;
	mffd1->Transform(x2, y2, z2); // Here computes the composition of T_2 o T_1
	mffd2->Transform(x2, y2, z2);
        dx[index] = x2 - x1;
        dy[index] = y2 - y1;
        dz[index] = z2 - z1;
        index++;
      }
    }
  }

  // Here approximates the composition using a single-level FFD
  ffd->Interpolate(dx, dy, dz);

  // The approximatation error
  index = 0;
  double mean_error = 0;
  double max_error = -1E10;
  for(int k=0; k<nz; k++){
    for(int j=0; j<ny; j++){
      for(int i=0; i<nx; i++){
        double x1 = i;
        double y1 = j;
        double z1 = k;
        ffd->LatticeToWorld(x1, y1, z1);
        double x2 = x1;
        double y2 = y1;
        double z2 = z1;
	mffd1->Transform(x2, y2, z2); // The composition
	mffd2->Transform(x2, y2, z2);
	ffd->Transform(x1, y1, z1);   // The approximated composition

        double error = sqrt((x1-x2)*(x1-x2) + (y1-y2)*(y1-y2) + (z1-z2)*(z1-z2));
        mean_error += error;
        if(error > max_error){
	  max_error = error;
	}
        index++;
      }
    }
  }
  mean_error /= (nx * ny * nz);

  cout << "Mean error of the composed FFD is = " << mean_error << endl;
  cout << "Max error of the composed FFD is = " << max_error << endl;

  // Create a multi-level free-form transformation
  irtkMultiLevelFreeFormTransformation *mffd = new irtkMultiLevelFreeFormTransformation();
  mffd->PushLocalTransformation(ffd);

  // Write the transformation
  mffd->irtkTransformation::Write(dofout_name);

  // Free memory
  delete []dx;
  delete []dy;
  delete []dz;
}
