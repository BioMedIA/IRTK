/*=========================================================================

  Library   : Image Registration Toolkit (IRTK)
  Module    : $Id$
  Copyright : Imperial College, Department of Computing
              Visual Information Processing (VIP), 2008 onwards
  Date      : $Date$
  Version   : $Revision$
  Changes   : $Author$

=========================================================================*/

#include <irtkRegistration.h>

#include <irtkFluidFreeFormTransformation.h>

#include <irtkImageFluidRegistration.h>

#include <irtkImageFluidRegistration2D.h>

#ifdef HAS_MPI

#include <mpi.h>

int mpi_nprocs = 1, mpi_rank;

#endif

// Default filenames
char *source_name = NULL, *target_name = NULL;
char *dofin_name  = NULL, *dofout_name = NULL;
char *parameter_name = NULL;

#ifdef HAS_VTK
extern Bool interactiveVTK;
extern Bool displayVTK;
extern Bool firstVTK;
#endif

#ifdef HAS_MPI
double wtime;
#endif

void usage()
{
  cerr << "Usage: nreg [target] [source] <options> \n" << endl;
  cerr << "where <options> is one or more of the following:\n" << endl;
  cerr << "<-parameter file>    Parameter file" << endl;
  cerr << "<-dofout    file>    Final transformation estimate" << endl;
  cerr << "<-dofin     file>    Initial transformation estimate" << endl;
  cerr << "<-Rx1 value>         Region of interest in both images" << endl;
  cerr << "<-Ry1 value>         Region of interest in both images" << endl;
  cerr << "<-Rz1 value>         Region of interest in both images" << endl;
  cerr << "<-Rx2 value>         Region of interest in both images" << endl;
  cerr << "<-Ry2 value>         Region of interest in both images" << endl;
  cerr << "<-Rz2 value>         Region of interest in both images" << endl;
  cerr << "<-Tx1 value>         Region of interest in target image" << endl;
  cerr << "<-Ty1 value>         Region of interest in target image" << endl;
  cerr << "<-Tz1 value>         Region of interest in target image" << endl;
  cerr << "<-Tx2 value>         Region of interest in target image" << endl;
  cerr << "<-Ty2 value>         Region of interest in target image" << endl;
  cerr << "<-Tz2 value>         Region of interest in target image" << endl;
  cerr << "<-Sx1 value>         Region of interest in source image" << endl;
  cerr << "<-Sy1 value>         Region of interest in source image" << endl;
  cerr << "<-Sz1 value>         Region of interest in source image" << endl;
  cerr << "<-Sx2 value>         Region of interest in source image" << endl;
  cerr << "<-Sy2 value>         Region of interest in source image" << endl;
  cerr << "<-Sz2 value>         Region of interest in source image" << endl;
  cerr << "<-Tp  value>         Padding value in target image" << endl;
  cerr << "<-debug>             Enable debugging information" << endl;
  exit(1);
}

int main(int argc, char **argv)
{
  int ok;
  int target_x1, target_y1, target_z1, target_x2, target_y2, target_z2;
  int source_x1, source_y1, source_z1, source_x2, source_y2, source_z2;

  // Check command line
  if (argc < 3) {
    usage();
  }

#ifdef HAS_MPI
  MPI_Init(&argc, &argv);
  wtime = MPI_Wtime();
  MPI_Comm_size(MPI_COMM_WORLD, &mpi_nprocs);
  MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);
#endif

  // Parse source and target images
  target_name = argv[1];
  argc--;
  argv++;
  source_name = argv[1];
  argc--;
  argv++;

  // Read target image
  cout << "Reading target ... "; cout.flush();
  irtkGreyImage target(target_name);
  cout << "done" << endl;
  // Read source image
  cout << "Reading source ... "; cout.flush();
  irtkGreyImage source(source_name);
  cout << "done" << endl;

  // Fix ROI
  target_x1 = 0;
  target_y1 = 0;
  target_z1 = 0;
  target_x2 = target.GetX();
  target_y2 = target.GetY();
  target_z2 = target.GetZ();
  source_x1 = 0;
  source_y1 = 0;
  source_z1 = 0;
  source_x2 = source.GetX();
  source_y2 = source.GetY();
  source_z2 = source.GetZ();

  // Create registration filter
  irtkImageFluidRegistration *registration = NULL;
  if ((target.GetZ() == 1) && (source.GetZ() == 1)) {
    registration = new irtkImageFluidRegistration2D;
  } else {
    registration = new irtkImageFluidRegistration;
  }

  // Create initial fluid free-form deformation
  irtkFluidFreeFormTransformation *mffd = NULL;

  // Parse remaining parameters
  while (argc > 1) {
    ok = False;
    if ((ok == False) && (strcmp(argv[1], "-dofin") == 0)) {
      argc--;
      argv++;
      dofin_name = argv[1];
      argc--;
      argv++;
      ok = True;
    }
    if ((ok == False) && (strcmp(argv[1], "-dofout") == 0)) {
      argc--;
      argv++;
      dofout_name = argv[1];
      argc--;
      argv++;
      ok = True;
    }
    if ((ok == False) && (strcmp(argv[1], "-Rx1") == 0)) {
      argc--;
      argv++;
      target_x1 = source_x1 = atoi(argv[1]);
      argc--;
      argv++;
      ok = True;
    }
    if ((ok == False) && (strcmp(argv[1], "-Rx2") == 0)) {
      argc--;
      argv++;
      target_x2 = source_x2 = atoi(argv[1]);
      argc--;
      argv++;
      ok = True;
    }
    if ((ok == False) && (strcmp(argv[1], "-Ry1") == 0)) {
      argc--;
      argv++;
      target_y1 = source_y1 = atoi(argv[1]);
      argc--;
      argv++;
      ok = True;
    }
    if ((ok == False) && (strcmp(argv[1], "-Ry2") == 0)) {
      argc--;
      argv++;
      target_y2 = source_y2 = atoi(argv[1]);
      argc--;
      argv++;
      ok = True;
    }
    if ((ok == False) && (strcmp(argv[1], "-Rz1") == 0)) {
      argc--;
      argv++;
      target_z1 = source_z1 = atoi(argv[1]);
      argc--;
      argv++;
      ok = True;
    }
    if ((ok == False) && (strcmp(argv[1], "-Rz2") == 0)) {
      argc--;
      argv++;
      target_z2 = source_z2 = atoi(argv[1]);
      argc--;
      argv++;
      ok = True;
    }
    if ((ok == False) && (strcmp(argv[1], "-Tx1") == 0)) {
      argc--;
      argv++;
      target_x1 = atoi(argv[1]);
      argc--;
      argv++;
      ok = True;
    }
    if ((ok == False) && (strcmp(argv[1], "-Tx2") == 0)) {
      argc--;
      argv++;
      target_x2 = atoi(argv[1]);
      argc--;
      argv++;
      ok = True;
    }
    if ((ok == False) && (strcmp(argv[1], "-Ty1") == 0)) {
      argc--;
      argv++;
      target_y1 = atoi(argv[1]);
      argc--;
      argv++;
      ok = True;
    }
    if ((ok == False) && (strcmp(argv[1], "-Ty2") == 0)) {
      argc--;
      argv++;
      target_y2 = atoi(argv[1]);
      argc--;
      argv++;
      ok = True;
    }
    if ((ok == False) && (strcmp(argv[1], "-Tz1") == 0)) {
      argc--;
      argv++;
      target_z1 = atoi(argv[1]);
      argc--;
      argv++;
      ok = True;
    }
    if ((ok == False) && (strcmp(argv[1], "-Tz2") == 0)) {
      argc--;
      argv++;
      target_z2 = atoi(argv[1]);
      argc--;
      argv++;
      ok = True;
    }
    if ((ok == False) && (strcmp(argv[1], "-Sx1") == 0)) {
      argc--;
      argv++;
      source_x1 = atoi(argv[1]);
      argc--;
      argv++;
      ok = True;
    }
    if ((ok == False) && (strcmp(argv[1], "-Sx2") == 0)) {
      argc--;
      argv++;
      source_x2 = atoi(argv[1]);
      argc--;
      argv++;
      ok = True;
    }
    if ((ok == False) && (strcmp(argv[1], "-Sy1") == 0)) {
      argc--;
      argv++;
      source_y1 = atoi(argv[1]);
      argc--;
      argv++;
      ok = True;
    }
    if ((ok == False) && (strcmp(argv[1], "-Sy2") == 0)) {
      argc--;
      argv++;
      source_y2 = atoi(argv[1]);
      argc--;
      argv++;
      ok = True;
    }
    if ((ok == False) && (strcmp(argv[1], "-Sz1") == 0)) {
      argc--;
      argv++;
      source_z1 = atoi(argv[1]);
      argc--;
      argv++;
      ok = True;
    }
    if ((ok == False) && (strcmp(argv[1], "-Sz2") == 0)) {
      argc--;
      argv++;
      source_z2 = atoi(argv[1]);
      argc--;
      argv++;
      ok = True;
    }
    if ((ok == False) && (strcmp(argv[1], "-Tp") == 0)) {
      argc--;
      argv++;
      registration->SetTargetPadding(atoi(argv[1]));
      argc--;
      argv++;
      ok = True;
    }
    if ((ok == False) && (strcmp(argv[1], "-debug") == 0)) {
      argc--;
      argv++;
      ok = True;
      registration->SetDebugFlag(True);
    }
    if ((ok == False) && (strcmp(argv[1], "-parameter") == 0)) {
      argc--;
      argv++;
      ok = True;
      parameter_name = argv[1];
      argc--;
      argv++;
    }
    if (ok == False) {
      cerr << "Can not parse argument " << argv[1] << endl;
      usage();
    }
  }

  // If there is an region of interest, use it
  if ((target_x1 != 0) || (target_x2 != target.GetX()) ||
      (target_y1 != 0) || (target_y2 != target.GetY()) ||
      (target_z1 != 0) || (target_z2 != target.GetZ())) {
    target = target.GetRegion(target_x1, target_y1, target_z1,
                              target_x2, target_y2, target_z2);
  }

  // If there is an region of interest for the source image, use it
  if ((source_x1 != 0) || (source_x2 != source.GetX()) ||
      (source_y1 != 0) || (source_y2 != source.GetY()) ||
      (source_z1 != 0) || (source_z2 != source.GetZ())) {
    source = source.GetRegion(source_x1, source_y1, source_z1,
                              source_x2, source_y2, source_z2);
  }

  if (dofin_name != NULL) {
    irtkTransformation *transform = irtkTransformation::New(dofin_name);
    if (strcmp(transform->NameOfClass(), "irtkRigidTransformation") == 0) {
      mffd = new irtkFluidFreeFormTransformation(*((irtkRigidTransformation *)transform));
      delete transform;
    } else {
      if (strcmp(transform->NameOfClass(), "irtkAffineTransformation") == 0) {
        mffd = new irtkFluidFreeFormTransformation(*((irtkAffineTransformation *)transform));
        delete transform;
      } else {
        if (strcmp(transform->NameOfClass(), "irtkFluidFreeFormTransformation") == 0) {
          mffd = dynamic_cast<irtkFluidFreeFormTransformation *>(transform);
        } else {
          cerr << "Input transformation is not of type rigid or affine " << endl;
          cerr << "or fluid free form deformation" << endl;
        }
      }
    }
  } else {
    // Otherwise use identity transformation to start
    mffd = new irtkFluidFreeFormTransformation;
  }

  // Set input and output for the registration filter
  registration->SetInput(&target, &source);
  registration->SetOutput(mffd);

  // Read parameter if there any, otherwise make an intelligent guess
  if (parameter_name != NULL) {
    registration->irtkImageRegistration::Read(parameter_name);
  }

  // Run registration filter
  registration->Run();

#ifdef HAS_MPI
  if (mpi_rank == 0)
#endif
  {
    // Write the final transformation estimate
    if (dofout_name != NULL) {
      mffd->irtkTransformation::Write(dofout_name);
    }
  }

#ifdef HAS_MPI
  cout << "Processor #" << mpi_rank << " spent " << (MPI_Wtime()-wtime) << "s" << endl;
  MPI_Finalize();
#endif
}
