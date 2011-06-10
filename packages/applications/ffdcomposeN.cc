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
#include <irtkImageFunction.h>
#include <irtkTransformation.h>
#include <irtkRegistration.h>

// Default filenames
char *dofOut_name = NULL;
char **dof_name  = NULL;
char *affine_dof_name = NULL;

#define MAX_DOFS 50
#define MAX_PTS_PAREG 10000

void usage()
{
  cerr << " Usage: ffdcomposeN [dofOut] <options>\n" << endl;
  cerr << " where <options> is one or more of the following: \n" << endl;
  cerr << " " << endl;
  cerr << " <-dofin file>      A transformation file which should be included. Multiple transformations " << endl;
  cerr << "                    can be given by repeatedly using this flag. They are processed in order " << endl;
  cerr << "                    in which they are passed. " << endl;
  cerr << " <-dofin_i file>    A transformation whose inverse should be included." << endl;
  cerr << " <-affine> file     An affine transformation that should be included within the composed" << endl;
  cerr << "                    FFD. If not provided then one will be estimated." << endl;
  cerr << " " << endl;
  cerr << " " << endl;
  cerr << " E.g." << endl;
  cerr << " " << endl;
  cerr << "   ffdcomposeN out.dof.gz -dofin a.dof -dofin_i b.dof -dofin c.dof" << endl;
  cerr << " " << endl;
  cerr << " returns the composed transformation c b^-1 a(x) applied to each target location x." << endl;
  cerr << " " << endl;
  exit(1);
}

int main(int argc, char **argv)
{
  irtkTransformation **transformation = NULL;
  int i, j, k, m, xdim, ydim, zdim, noOfDofs, numberOfCPs, count;
  double x, y, z, xStore, yStore, zStore;
  bool ok;

  // Check command line
  if (argc < 3) {
    usage();
  }

  // Parse dofOut_name
  dofOut_name = argv[1];
  argc--;
  argv++;

  // Fix number of dofs
  noOfDofs = 0;

  dof_name = new char*[MAX_DOFS];
  bool *invert = new bool[MAX_DOFS];

  while (argc > 1) {
    ok = false;
    if ((ok == false) && (strcmp(argv[1], "-dofin") == 0)) {
      argc--;
      argv++;
      dof_name[noOfDofs] = argv[1];
      invert[noOfDofs]   = false;
      noOfDofs++;
      argc--;
      argv++;
      ok = true;
    }
    if ((ok == false) && (strcmp(argv[1], "-dofin_i") == 0)) {
      argc--;
      argv++;
      dof_name[noOfDofs] = argv[1];
      invert[noOfDofs]   = true;
      noOfDofs++;
      argc--;
      argv++;
      ok = true;
    }
    if ((ok == false) && (strcmp(argv[1], "-affine") == 0)) {
      argc--;
      argv++;
      affine_dof_name = argv[1];
      argc--;
      argv++;
      ok = true;
    }

    if (ok == false) {
      cerr << "Can not parse argument " << argv[1] << endl;
      usage();
    }
  }

  if (noOfDofs == 0) {
    noOfDofs = 1;
    transformation = new irtkTransformation*[noOfDofs];
    transformation[0] = new irtkRigidTransformation;
  } else {
    transformation = new irtkTransformation*[noOfDofs];
    for (m = 0; m < noOfDofs; m++) {
      transformation[m] = irtkTransformation::New(dof_name[m]);
    }
  }

  irtkMultiLevelFreeFormTransformation *mffd_out = new irtkMultiLevelFreeFormTransformation;
  mffd_out->irtkTransformation::Read(dof_name[0]);  // takes the first DOF

  // Extract FFD and get lattice dimensions
  irtkFreeFormTransformation3D *affd_out = dynamic_cast<irtkFreeFormTransformation3D *>(mffd_out->PopLocalTransformation());

  while (mffd_out->NumberOfLevels() > 0)
    (void) mffd_out->PopLocalTransformation();

  xdim = affd_out->GetX();
  ydim = affd_out->GetY();
  zdim = affd_out->GetZ();

  numberOfCPs = xdim * ydim * zdim;

  // Space to store the control point displacements.
  double *xdata = new double[numberOfCPs];
  double *ydata = new double[numberOfCPs];
  double *zdata = new double[numberOfCPs];

  count = 0;
  // Loop for each control point in the target
  for (k = 0; k < zdim; k++) {
    for (j = 0; j < ydim; j++) {
      for (i = 0; i < xdim; i++) {

        affd_out->Put(i, j, k, 0.0, 0.0, 0.0);

        x = i;
        y = j;
        z = k;

        // Transform point from lattice coordinates to target coordinates
        affd_out->LatticeToWorld(x, y, z);

        xStore = x;
        yStore = y;
        zStore = z;

        // Transform point
        for (m = 0; m < noOfDofs; m++) {
          if (invert[m])
            transformation[m]->Inverse(x, y, z);
          else
            transformation[m]->Transform(x, y, z);
        }

        // Displacement of the control point.
        xdata[count] = x - xStore;
        ydata[count] = y - yStore;
        zdata[count] = z - zStore;

        ++count;

      }
    }
  }

  // Make an identity global transformation.
  irtkAffineTransformation *trAffine = new irtkAffineTransformation;

  // Read in affine transformation if provided, otherwise estimate it.
  if (affine_dof_name != NULL) {

    trAffine->irtkTransformation::Read(affine_dof_name);

  } else {

    // Estimate an affine component for the collected displacements.
    // This will be a compromise, especially if there are lots of control points
    // outside the region of interest (e.g the brain) with a zero displacement.


    irtkPointSet targetPts;
    irtkPointSet sourcePts;

    // Collect point data.

    int noOfPoints;
    noOfPoints = affd_out->NumberOfDOFs() / 3;

    int incr = 1;
    while ((noOfPoints / incr) > MAX_PTS_PAREG) {
      incr++;
    }

    // Subsample uniformly by increments of 'incr'

    count = -1;

    // Loop over all control points.
    for (k = 0; k < zdim; k++) {
      for (j = 0; j < ydim; j++) {
        for (i = 0; i < xdim; i++) {

          count++;


          // Should we sample it or not?
          if ((count % incr) != 0)
            continue;

          // Get two copies of current image coordinates.
          x = i; y = j; z = k;
          irtkPoint p(x, y, z);
          irtkPoint q(x, y, z);

          // Transform points into target world coordinates.
          affd_out->LatticeToWorld(p);
          affd_out->LatticeToWorld(q);


          // The starting point
          targetPts.Add(p);

          // Where does p go to under all the compositions?
          q._x += xdata[count];
          q._y += ydata[count];
          q._z += zdata[count];

          sourcePts.Add(q);
        }
      }
    }

    // Estimate global affine component
    irtkPointAffineRegistration *pareg = new irtkPointAffineRegistration;
    // Set input and output for the registration filter
    irtkPointSet tmp1 = targetPts;
    irtkPointSet tmp2 = sourcePts;
    pareg->SetInput(&tmp1, &tmp2);
    pareg->SetOutput(trAffine);

    // Run registration filter
    pareg->Run();

    cout << "Estimated global affine component with following parameters:" << endl;
    trAffine->Print();
  }

  // Remove the global affine part from the estimated displacements
  count = 0;
  for (k = 0; k < zdim; k++) {
    for (j = 0; j < ydim; j++) {
      for (i = 0; i < xdim; i++) {
        x = i;
        y = j;
        z = k;

        // Transform point from lattice coordinates to target coordinates
        affd_out->LatticeToWorld(x, y, z);
        xStore = x;
        yStore = y;
        zStore = z;
        trAffine->Transform(x, y, z);

        // Subtract displacement due to affine transformation of the control point location.
        xdata[count] -= (x - xStore);
        ydata[count] -= (y - yStore);
        zdata[count] -= (z - zStore);
        ++count;
      }
    }
  }

  mffd_out->PutMatrix(trAffine->GetMatrix());

  // Interpolate the ffd and write dof
  affd_out->Interpolate(xdata, ydata, zdata);
  mffd_out->PushLocalTransformation(affd_out);
  mffd_out->irtkTransformation::Write(dofOut_name);

  delete [] xdata;
  delete [] ydata;
  delete [] zdata;
  delete [] dof_name;
  delete [] invert;

}

