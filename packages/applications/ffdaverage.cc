/*=========================================================================

  Library   : Image Registration Toolkit (IRTK)
  Module    : $Id$
  Copyright : Imperial College, Department of Computing
              Visual Information Processing (VIP), 2008 onwards
  Date      : $Date$
  Version   : $Revision$
  Changes   : $Author$

=========================================================================*/

#include <irtkTransformation.h>

char *dofout = NULL, *prefix_name = NULL, *textfile = NULL;

#define MAX_DOFS   10000
#define EPSILON    0.001

void usage()
{
  cerr << "   Usage: ffdaverage [dofout] [dofin1...N] <options>\n" << endl;

  cerr << "   Average a set of input mffds.\n" << endl;

  cerr << "   The input mffds are assumed to relate to the same target image" << endl;
  cerr << "   with the same lattice dimensions and a single level. The output " << endl;
  cerr << "   mffd will share this target and will have a single level with" << endl;
  cerr << "   the same control point structure as the input transformations.\n" << endl;

  cerr << "   By default, the rigid parameters of the output mffd will be the " << endl;
  cerr << "   identity and the scale and skew parameters are average of those " << endl;
  cerr << "   of the input mffds.\n" << endl;

  cerr << "   The input transformation control points parameters are averaged " << endl;
  cerr << "   linearly and assigned to each output control point.  Affine " << endl;
  cerr << "   correction is first done on the input transformations' local " << endl;
  cerr << "   control point values.\n" << endl ;

  cerr << "   By default, an identity mffd from the target to itself is assumed " << endl;
  cerr << "   present, i.e. the number of transformations is assumed to be one" << endl;
  cerr << "   higher than the number given as arguments.\n" << endl;

  cerr << "   Options:\n" << endl;

  cerr << "   -noaffine           	= Don't combine with affine average, i.e. " << endl;
  cerr << "                   			  resulting transformation will have identity global" << endl;
  cerr << "                   				component." << endl;
  cerr << "   -divisor      				= This is taken by default to be one higher than " << endl;
  cerr << "                   				the number of input ffds because an identity " << endl;
  cerr << "                   				transformation is assumed present.  Set manually " << endl;
  cerr << "                  					for other values." << endl;
  cerr << "   -noreset      				= Do not reset the rigid components of the global" << endl;
  cerr << "                   				matrices.\n" << endl;
  cerr << "   -dofnames file      		Name of file containing image names for <input1..inputN>\n";
  cerr << "   -prefix directory				Directory in which to find the images\n";
  cerr << "   -gaussian mean sigma  	Use Gaussian with mean and sigma for kernel-based smoothing\n";
  cerr << "   -epsilon value          Epsilon value for ignoring image in kernel-based smoothing (default = " << EPSILON << ")\n";

  exit(1);
}

void checkLevelsAndDimensions(irtkMultiLevelFreeFormTransformation *mffdFirst, irtkMultiLevelFreeFormTransformation *mffd)
{
  int x, y, z;

  x = mffd->GetLocalTransformation(0)->GetX();
  y = mffd->GetLocalTransformation(0)->GetY();
  z = mffd->GetLocalTransformation(0)->GetZ();

  if (x != mffdFirst->GetLocalTransformation(0)->GetX() ||
      y != mffdFirst->GetLocalTransformation(0)->GetY() ||
      z != mffdFirst->GetLocalTransformation(0)->GetZ()) {
    cerr << "All input transformations must have the same lattice dimensions." << endl;
    exit(1);
  }
}

void resetRigidComponents(irtkAffineTransformation *a)
{
  a->PutTranslationX(0);
  a->PutTranslationY(0);
  a->PutTranslationZ(0);
  a->PutRotationX(0);
  a->PutRotationY(0);
  a->PutRotationZ(0);
}

int main(int argc, char **argv)
{
  int ok, i, j, k, n, xdim, ydim, zdim, noOfCPs, dofCount, index, identity;
  double x, y, z, mean, sigma, identity_weight, epsilon, xAv, yAv, zAv, xFinal, yFinal, zFinal;
  int combineaffine = True;
  int resetRigid = True;

  irtkMatrix *globMats;
  irtkMatrix globMatAv(4, 4);
  irtkMatrix sumLogs(4, 4);
  irtkMatrix jac(3, 3);
  irtkMatrix globAvJac(3, 3);

  if (argc < 4) {
    usage();
  }

  // Parse arguments.
  dofout = argv[1];
  argc--;
  argv++;

  // Parse input file names and values
  char  **dofin_name   = new char *[MAX_DOFS];
  double *dofin_value  = new double[MAX_DOFS];
  double *dofin_weight = new double[MAX_DOFS];
  double  total_weight = 0;

  // Parse any remaining paramters
  dofCount = 0;

  while ((argc > 1) && (argv[1][0] != '-' )) {
    dofin_name[dofCount]   = argv[1];
    dofin_weight[dofCount] = 1;
    total_weight          += 1;
    dofCount++;
    argc--;
    argv++;
  }

  for (i = 0; i < dofCount; ++i) {
    dofin_name[i] = argv[1];
    argc--;
    argv++;
  }
  
  if (dofCount > 0){
    cout << "Read " << dofCount << " transformations from command line." << endl;
  }

  // Default: No kernel smoothing
  mean  = 0;
  sigma = 1;

  // Default: Epsilon
  epsilon = EPSILON;

  // Default: Assume identity transform should not be added
  identity = False;
  identity_weight = 0;

  // Parse options.
  while (argc > 1) {
    ok = False;
    if ((ok == False) && (strcmp(argv[1], "-noaffine") == 0)) {
      argc--;
      argv++;
      combineaffine = False;
      ok = True;
    }
    if ((ok == False) && (strcmp(argv[1], "-identity") == 0)) {
      argc--;
      argv++;
      identity = True;
      argc--;
      argv++;
      identity_weight = atof(argv[1]);
      ok = True;
    }
    if ((ok == False) && (strcmp(argv[1], "-gaussian") == 0)) {
      argc--;
      argv++;
      mean = atof(argv[1]);
      argc--;
      argv++;
      sigma = atof(argv[1]);
      argc--;
      argv++;
      ok = True;
    }
    if ((ok == False) && (strcmp(argv[1], "-dofnames") == 0)) {
      argc--;
      argv++;
      textfile = argv[1];
      argc--;
      argv++;
      ok = True;
    }
    if ((ok == False) && (strcmp(argv[1], "-prefix") == 0)) {
      argc--;
      argv++;
      prefix_name = argv[1];
      argc--;
      argv++;
      ok = True;
    }
    if ((ok == False) && (strcmp(argv[1], "-noreset") == 0)) {
      argc--;
      argv++;
      resetRigid = False;
      ok = True;
    }
    if ((ok == False) && (strcmp(argv[1], "-epsilon") == 0)) {
      argc--;
      argv++;
      epsilon = atof(argv[1]);
      argc--;
      argv++;
      ok = True;
    }
    if (ok == False) {
      cerr << "Can not parse argument " << argv[1] << endl;
      usage();
    }
  }

  if (textfile != NULL) {

    if (dofCount > 0){
      cerr << "ffdaverage : Text file specified for transformations but " << dofCount << " have already been read." << endl;
      exit(1);
    }

    ifstream in(textfile);
    if (in.is_open()) {
      while (!in.eof()) {
        dofin_name[dofCount] = new char[256];
        in >> dofin_name[dofCount] >> dofin_value[dofCount];
        if (strlen(dofin_name[dofCount]) > 0) {
          dofin_weight[dofCount] = 1.0 / sqrt(2.0) / sigma * exp(-pow((mean - dofin_value[dofCount])/sigma, 2.0)/2);

          if (dofin_weight[dofCount] > epsilon) {
            total_weight   += dofin_weight[dofCount];
            cout << dofin_value[dofCount] << " " << dofin_weight[dofCount] << " " << total_weight << endl;
          } else {
            dofin_weight[dofCount] = 0.0f;
          }

          dofCount++;
        }
      }
      in.close();
    } else {
      cout << "ffdaverage: Unable to open file " << textfile << endl;
      exit(1);
    }

    cout << "Read " << dofCount << " transformations from file " << textfile << endl;

  }

  dofin_weight[dofCount] = identity_weight;
  total_weight += identity_weight;

  for (i = 0; i < dofCount; ++i) {
    char buffer[255];

    // Read input
    if (prefix_name != NULL) {
      sprintf(buffer, "%s/%s", prefix_name, dofin_name[i]);
    } else {
      sprintf(buffer, "%s", dofin_name[i]);
    }
    dofin_name[i] = strdup(buffer);
  }

  if (dofCount < 1) {
    cerr << "ffdaverage: Must have at least one input transformation." << endl;
    exit(1);
  }

  // Read the mffds.
  // The output mffd is a copy of the first input (for now).
  irtkTransformation *transform;
  irtkMultiLevelFreeFormTransformation *mffdOut;

  irtkMultiLevelFreeFormTransformation *mffdFirst;
  irtkMultiLevelFreeFormTransformation *mffdIn;

  transform = irtkTransformation::New(dofin_name[0]);
  mffdOut   = dynamic_cast<irtkMultiLevelFreeFormTransformation *>(transform);
  mffdFirst = dynamic_cast<irtkMultiLevelFreeFormTransformation *>(transform);

  // Get space for input transformations and one extra for the identity
  // assumed from the target to itself.
  if (identity == True) {
    globMats = new irtkMatrix[dofCount + 1];
  } else {
    globMats = new irtkMatrix[dofCount];
  }

  for (i = 0; i < dofCount; ++i) {
    globMats[i].Initialize(4, 4);
  }


  ////////////////////////////////////////////////////
  // Deal with the affine part:

  globMatAv.Ident();

  // Is affine portion to be averaged?
  if (combineaffine == True) {

    // Average the affine parts separately from the displacement
    // fields.  By default, only the scales and skew are actually
    // averaged.

    for (i = 0; i < dofCount; ++i) {
      cerr << "Reading " << dofin_name[i] << " for affine averaging" << endl;
      transform = irtkTransformation::New(dofin_name[i]);
      mffdIn = dynamic_cast<irtkMultiLevelFreeFormTransformation *>(transform);

      if (resetRigid == True) {
        resetRigidComponents(mffdIn);
      }

      globMats[i] = mffdIn->irtkAffineTransformation::GetMatrix();
      delete mffdIn;
    }

    // If required, setting the divisor appropriately can ensure
    // the last identity matrix is ignored.
    if (identity == True) {
      // Put the identity matrix at the end.
      globMats[dofCount].Initialize(4, 4);
      globMats[dofCount].Ident();
      globMatAv = FrechetMean(globMats, dofin_weight, dofCount + 1, 20);
    } else {
      globMatAv = FrechetMean(globMats, dofin_weight, dofCount, 20);
    }
  }

  mffdOut->irtkAffineTransformation::PutMatrix(globMatAv);


  ////////////////////////////////////////////////////
  // Now tackle the displacement fields.

  // Copy the average Jacobian section for correcting the control
  // point values in the output transformation.
  for (i = 0; i < 3; ++i) {
    for (j = 0; j < 3; ++j) {
      globAvJac(i, j) = globMatAv(i, j);
    }
  }

  // Space for storing the output control points.
  irtkFreeFormTransformation3D *ffdOut = dynamic_cast<irtkFreeFormTransformation3D *>(mffdOut->GetLocalTransformation(0));

  if (ffdOut == NULL) {
    cerr << "Free-form transformation is not 3D" << endl;
    exit(1);
  }

  // Space for storing the control point displacements.
  double* xdata = new double[ffdOut->NumberOfDOFs()/3];
  double* ydata = new double[ffdOut->NumberOfDOFs()/3];
  double* zdata = new double[ffdOut->NumberOfDOFs()/3];

  xdim = ffdOut->GetX();
  ydim = ffdOut->GetY();
  zdim = ffdOut->GetZ();

  noOfCPs = xdim * ydim * zdim;
  for (i = 0; i < noOfCPs; ++i) {
    xdata[i] = ydata[i] = zdata[i] = 0.0;
  }

  // Correct the control point values in each of the input FFDs with
  // respect to their corresponding affine components.
  for (n = 0; n < dofCount; ++n) {

    cerr << "Reading " << dofin_name[n] << " for non-rigid averaging" << endl;
    transform = irtkTransformation::New(dofin_name[n]);
    mffdIn = dynamic_cast<irtkMultiLevelFreeFormTransformation *>(transform);

    checkLevelsAndDimensions(mffdFirst, mffdIn);

    // Re-read the original global matrices.
    globMats[n] = mffdIn->irtkAffineTransformation::GetMatrix();

    // Correcting the cp values requires the affine part to be S->T.
    globMats[n].Invert();

    irtkFreeFormTransformation3D *ffdIn = dynamic_cast<irtkFreeFormTransformation3D *>(mffdIn->GetLocalTransformation(0));

    if (ffdIn == NULL) {
      cerr << "Free-form transformation is not 3D" << endl;
      exit(1);
    }

    index = 0;
    for (k = 0; k < zdim; ++k) {
      for (j = 0; j < ydim; ++j) {
        for (i = 0; i < xdim; ++i) {
          ffdIn->Get(i, j, k, x, y, z);

          xdata[index] += dofin_weight[n] * (globMats[n](0, 0) * x + globMats[n](0, 1) * y + globMats[n](0, 2) * z);
          ydata[index] += dofin_weight[n] * (globMats[n](1, 0) * x + globMats[n](1, 1) * y + globMats[n](1, 2) * z);
          zdata[index] += dofin_weight[n] * (globMats[n](2, 0) * x + globMats[n](2, 1) * y + globMats[n](2, 2) * z);

          ++index;
        }
      }
    }

    delete mffdIn;
  }

  // Now can average the displacement fields.
  index = 0;
  for (k = 0; k < zdim; ++k) {
    for (j = 0; j < ydim; ++j) {
      for (i = 0; i < xdim; ++i) {

        // Divisor needs to be set appropriately depending on
        // whether an identity FFD is assumed to be present or
        // not.  If it is then the divisor is set by default to 1
        // + dofCount.

        xAv = xdata[index] / total_weight;
        yAv = ydata[index] / total_weight;
        zAv = zdata[index] / total_weight;
        ++index;

        if (combineaffine == True) {
          // Now introduce the affine component.
          xFinal = globAvJac(0, 0) * xAv + globAvJac(0, 1) * yAv + globAvJac(0, 2) * zAv;
          yFinal = globAvJac(1, 0) * xAv + globAvJac(1, 1) * yAv + globAvJac(1, 2) * zAv;
          zFinal = globAvJac(2, 0) * xAv + globAvJac(2, 1) * yAv + globAvJac(2, 2) * zAv;
        } else {
          xFinal = xAv;
          yFinal = yAv;
          zFinal = zAv;
        }

        ffdOut->Put(i, j, k, xFinal, yFinal, zFinal);
      }
    }
  }

  // Some reporting.
  if (combineaffine == True) {
    cout << "Affine component combined with ffd." << endl;
  } else {
    cout << "ffd has identity affine component." << endl;
  }

  if (resetRigid == True) {
    cout << "Rigid components reset to identity." << endl;

  } else {
    cout << "Rigid components used during averaging affine matrices." << endl;
  }

  mffdOut->irtkTransformation::Write(dofout);

  delete [] xdata;
  delete [] ydata;
  delete [] zdata;
}
