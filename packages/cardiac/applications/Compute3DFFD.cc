





#include <irtkRegistration2.h>

void usage()
{
  cerr << "Usage: compute3DFFD [dofin] [dofoutPrefix] [N] <options> \n" << endl;
  cerr << "takes in a 4D transformation and computes the [N] 3D transformation for time points n*1/N \n" << endl;

  exit(1);
}

int main(int argc, char **argv)
{
  double t, dt;
  int N, n, Nx, Ny, Nz, i;
  double x1, y1, z1, x2, y2, z2, dx, dy, dz, xaxis[3], yaxis[3], zaxis[3];
  double *cpX = NULL, *cpY = NULL, *cpZ = NULL, *dispX = NULL, *dispY = NULL, *dispZ = NULL, error, errorOld;

  int maxIter = 100;

  // Default filenames
  char *dofin_name  = NULL, *dofout_prefix = NULL, *dofout_current = NULL;

  // Check command line
  if (argc < 4) {
    usage();
  }

  // Parse source and target images
  dofin_name = argv[1];
  argc--;
  argv++;
  dofout_prefix = argv[1];
  argc--;
  argv++;
  N = atoi(argv[1]);
  argc--;
  argv++;

  dt = 1./double(N);
  irtkBSplineFreeFormTransformation *ffd = NULL;

  // Read target image
  cout << "Reading transformation ... "; cout.flush();
  irtkTransformation *transform = irtkTransformation::New(dofin_name);
  irtkMultiLevelFreeFormTransformation *mffd = NULL;
  irtkBSplineFreeFormTransformationPeriodic *tffd = NULL;
  if (strcmp(transform->NameOfClass(), "irtkMultiLevelFreeFormTransformation") == 0) {
    mffd = new irtkMultiLevelFreeFormTransformation(*((irtkMultiLevelFreeFormTransformation *)transform));
  } else {
    cerr << "Transformation has to be of type irtkMultiLevelFreeFormTransformation" << endl;
    exit(1);
  }
  tffd = (irtkBSplineFreeFormTransformationPeriodic *)mffd->GetLocalTransformation(0);
  cout << "done" << endl;

  // initialisation parameters
  Nx = tffd->GetX();
  Ny = tffd->GetY();
  Nz = tffd->GetZ();
  x1 = y1 = z1 = 0;
  tffd->LatticeToWorld(x1, y1, z1);
  x2 = Nx-1;
  y2 = Ny-1;
  z2 = Nz-1;
  tffd->LatticeToWorld(x2, y2, z2);
  dx = tffd->GetXSpacing();
  dy = tffd->GetYSpacing();
  dz = tffd->GetZSpacing();
  tffd->GetOrientation(xaxis, yaxis, zaxis);

  // control points
  cpX = new double[Nx*Ny*Nz];
  cpY = new double[Nx*Ny*Nz];
  cpZ = new double[Nx*Ny*Nz];
  for (i = 0; i < Nx*Ny*Nz; i++) {
	tffd->ControlPointLocation(i, cpX[i], cpY[i], cpZ[i]);
  }

  cout << "computing 3D transformations ... " << endl; cout.flush();
  for (n = 0; n < N; n++) {
cout<<n;
	t = double(n)*dt;

	// initialise 3D FFD
	ffd = new irtkBSplineFreeFormTransformation(x1, y1, z1, x2, y2, z2, dx, dy, dz, xaxis, yaxis, zaxis);

	error = 10;

	//////////////////////////////////////////
	// run Approximate function multiple times
    int iter = 0;
	do {
//cout<<", "<<iter++;
	  dispX = new double[Nx*Ny*Nz];
	  dispY = new double[Nx*Ny*Nz];
	  dispZ = new double[Nx*Ny*Nz];
	  for (i = 0; i < Nx*Ny*Nz; i++) {
	    dispX[i] = cpX[i];
	    dispY[i] = cpY[i];
	    dispZ[i] = cpZ[i];
	    tffd->Displacement(dispX[i], dispY[i], dispZ[i], t);
	  }

	  ffd->Approximate(cpX, cpY, cpZ, dispX, dispY, dispZ, Nx*Ny*Nz);

	  errorOld = error;

	  // calculate error
	  error = 0;
	  for (i = 0; i < Nx*Ny*Nz; i++) {
		error += sqrt( dispX[i]*dispX[i] + dispY[i]*dispY[i] + dispZ[i]*dispZ[i]);
	  }
	  error /= 3*Nx*Ny*Nz;

	  delete dispX;
	  delete dispY;
	  delete dispZ;

//cout<<" e: "<<error<<endl;
	} while(errorOld - error > 0.000001 && iter < maxIter); // check for convergence
	//////////////////////////////////////////
	//////////////////////////////////////////

	cout<<" -> final error: "<<error<<endl;

	// write new ffd at time t
	dofout_current = new char[200];
	sprintf(dofout_current, "%s_t%d.dof.gz", dofout_prefix, n);
	ffd->irtkTransformation::Write(dofout_current);

	delete dofout_current;
	delete ffd;
  }
  cout << "done" << endl;

  delete mffd;
}
