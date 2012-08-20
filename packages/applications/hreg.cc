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

#define MAX_LEVELS 100

// Default filenames
char *source_name;
char *target_name;
char *mask_name = NULL;
char *dofin_name;
char *dofout_name[MAX_LEVELS];
char *parameter_name[MAX_LEVELS];

int paddingValue, numberOfLevels, debug;

irtkImageFreeFormRegistrationMode mode;

#ifdef HAS_VTK
extern bool interactiveVTK;
extern bool displayVTK;
extern bool firstVTK;
#endif

void usage()
{
  cerr << "Usage: hreg [target] [source] [subdivisions] <options> \n" << endl;
  cerr << "where <options> is one or more of the following:\n" << endl;
  cerr << "<-parameter files>   Parameter files" << endl;
  cerr << "<-dofout    files>   Final transformation estimates" << endl;
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
  cerr << "<-ds  value>         Initial control point spacing" << endl;
  cerr << "<-dx  value>         Initial control point spacing for x" << endl;
  cerr << "<-dy  value>         Initial control point spacing for y" << endl;
  cerr << "<-dz  value>         Initial control point spacing for z" << endl;
  cerr << "<-debug>             Enable debugging information" << endl;
  cerr << "<-mask file>         Use a mask to define the ROI. The mask" << endl;
  cerr << "                     must have the same dimensions as the target." << endl;
  cerr << "                     Voxels in the mask with zero or less are " << endl;
  cerr << "                     padded in the target." << endl;
  cerr << "Note that you need to specify the number of times the FFD " << endl;
  cerr << "should be subdivided. You need to specify a filename for " << endl;
  cerr << "the parameters and final transformation estimate for each" << endl;
  cerr << "subdivision." << endl;
  exit(1);
}

void padding(irtkGreyImage image, irtkMultiLevelFreeFormTransformation *mffd)
{
  int i, j, k, x, y, z, x1, y1, z1, x2, y2, z2, ok, index;

  // Extract current transformation level
  irtkFreeFormTransformation3D *ffd = dynamic_cast<irtkFreeFormTransformation3D *>(mffd->GetLocalTransformation(mffd->NumberOfLevels()-1));

  if (ffd == NULL) {
    cerr << "Free-form transformation is not 3D" << endl;
    exit(1);
  }

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
              if (image(x, y, z) > paddingValue) {
                ok = true;
              }
            }
          }
        }
        if (ok == false) {
          ffd->PutStatusCP(i, j, k, _Passive, _Passive, _Passive);
        }
      }
    }
  }
}

void mreg(irtkGreyImage target, irtkGreyImage source,
          irtkMultiLevelFreeFormTransformation *mffd, int i)
{
  irtkImageFreeFormRegistration registration;

  // Set input and output for the registration filter
  registration.SetInput(&target, &source);
  registration.SetOutput(mffd);

  // Set padding value
  if (paddingValue > MIN_GREY) {
    registration.SetTargetPadding(paddingValue);
  }

  // Debug flag
  registration.SetDebugFlag(debug);

  // Set mode
  registration.SetMode(mode);

  // Read default parameter
  registration.irtkImageRegistration::Read(parameter_name[i]);

  // Run registration filter
  registration.Run();

  if (dofout_name[i] != NULL) {
    mffd->irtkTransformation::Write(dofout_name[i]);
  }

  // Extract current transformation level
  irtkFreeFormTransformation3D *affd = dynamic_cast<irtkFreeFormTransformation3D *>(mffd->PopLocalTransformation());

  if (affd == NULL) {
    cerr << "Free-form transformation is not 3D" << endl;
    exit(1);
  }

  affd->Subdivide();
  mffd->PushLocalTransformation(affd);
}

int main(int argc, char **argv)
{
  int i, ok, no_areg;
  int target_x1, target_y1, target_z1, target_x2, target_y2, target_z2;
  int source_x1, source_y1, source_z1, source_x2, source_y2, source_z2;
  double dx, dy, dz;
  irtkAffineTransformation transformation;

  // Check command line
  if (argc < 3) {
    usage();
  }

  // Parse source and target images
  target_name = argv[1];
  argc--;
  argv++;
  source_name = argv[1];
  argc--;
  argv++;
  numberOfLevels = atoi(argv[1]);
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

  // Fix padding
  paddingValue = MIN_GREY;

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

  // Fix spacing
  dx = 20;
  dy = 20;
  dz = 20;

  // Parse remaining parameters
  no_areg = false;
  debug   = false;
  mode    = RegisterXYZ;

  while (argc > 1) {
    ok = false;
    if ((ok == false) && (strcmp(argv[1], "-dofin") == 0)) {
      argc--;
      argv++;
      dofin_name = argv[1];
      argc--;
      argv++;
      ok = true;
    }
    if ((ok == false) && (strcmp(argv[1], "-dofout") == 0)) {
      argc--;
      argv++;
      for (i = 0; i < numberOfLevels; i++) {
        dofout_name[i] = argv[1];
        argc--;
        argv++;
      }
      ok = true;
    }
    if ((ok == false) && (strcmp(argv[1], "-Rx1") == 0)) {
      argc--;
      argv++;
      target_x1 = source_x1 = atoi(argv[1]);
      argc--;
      argv++;
      ok = true;
    }
    if ((ok == false) && (strcmp(argv[1], "-Rx2") == 0)) {
      argc--;
      argv++;
      target_x2 = source_x2 = atoi(argv[1]);
      argc--;
      argv++;
      ok = true;
    }
    if ((ok == false) && (strcmp(argv[1], "-Ry1") == 0)) {
      argc--;
      argv++;
      target_y1 = source_y1 = atoi(argv[1]);
      argc--;
      argv++;
      ok = true;
    }
    if ((ok == false) && (strcmp(argv[1], "-Ry2") == 0)) {
      argc--;
      argv++;
      target_y2 = source_y2 = atoi(argv[1]);
      argc--;
      argv++;
      ok = true;
    }
    if ((ok == false) && (strcmp(argv[1], "-Rz1") == 0)) {
      argc--;
      argv++;
      target_z1 = source_z1 = atoi(argv[1]);
      argc--;
      argv++;
      ok = true;
    }
    if ((ok == false) && (strcmp(argv[1], "-Rz2") == 0)) {
      argc--;
      argv++;
      target_z2 = source_z2 = atoi(argv[1]);
      argc--;
      argv++;
      ok = true;
    }
    if ((ok == false) && (strcmp(argv[1], "-Tx1") == 0)) {
      argc--;
      argv++;
      target_x1 = atoi(argv[1]);
      argc--;
      argv++;
      ok = true;
    }
    if ((ok == false) && (strcmp(argv[1], "-Tx2") == 0)) {
      argc--;
      argv++;
      target_x2 = atoi(argv[1]);
      argc--;
      argv++;
      ok = true;
    }
    if ((ok == false) && (strcmp(argv[1], "-Ty1") == 0)) {
      argc--;
      argv++;
      target_y1 = atoi(argv[1]);
      argc--;
      argv++;
      ok = true;
    }
    if ((ok == false) && (strcmp(argv[1], "-Ty2") == 0)) {
      argc--;
      argv++;
      target_y2 = atoi(argv[1]);
      argc--;
      argv++;
      ok = true;
    }
    if ((ok == false) && (strcmp(argv[1], "-Tz1") == 0)) {
      argc--;
      argv++;
      target_z1 = atoi(argv[1]);
      argc--;
      argv++;
      ok = true;
    }
    if ((ok == false) && (strcmp(argv[1], "-Tz2") == 0)) {
      argc--;
      argv++;
      target_z2 = atoi(argv[1]);
      argc--;
      argv++;
      ok = true;
    }
    if ((ok == false) && (strcmp(argv[1], "-Sx1") == 0)) {
      argc--;
      argv++;
      source_x1 = atoi(argv[1]);
      argc--;
      argv++;
      ok = true;
    }
    if ((ok == false) && (strcmp(argv[1], "-Sx2") == 0)) {
      argc--;
      argv++;
      source_x2 = atoi(argv[1]);
      argc--;
      argv++;
      ok = true;
    }
    if ((ok == false) && (strcmp(argv[1], "-Sy1") == 0)) {
      argc--;
      argv++;
      source_y1 = atoi(argv[1]);
      argc--;
      argv++;
      ok = true;
    }
    if ((ok == false) && (strcmp(argv[1], "-Sy2") == 0)) {
      argc--;
      argv++;
      source_y2 = atoi(argv[1]);
      argc--;
      argv++;
      ok = true;
    }
    if ((ok == false) && (strcmp(argv[1], "-Sz1") == 0)) {
      argc--;
      argv++;
      source_z1 = atoi(argv[1]);
      argc--;
      argv++;
      ok = true;
    }
    if ((ok == false) && (strcmp(argv[1], "-Sz2") == 0)) {
      argc--;
      argv++;
      source_z2 = atoi(argv[1]);
      argc--;
      argv++;
      ok = true;
    }
    if ((ok == false) && (strcmp(argv[1], "-Tp") == 0)) {
      argc--;
      argv++;
      paddingValue = atoi(argv[1]);
      argc--;
      argv++;
      ok = true;
    }
    if ((ok == false) && (strcmp(argv[1], "-y_only") == 0)) {
      argc--;
      argv++;
      mode = RegisterY;
      ok = true;
    }
    if ((ok == false) && (strcmp(argv[1], "-xy_only") == 0)) {
      argc--;
      argv++;
      mode = RegisterXY;
      ok = true;
    }
    if ((ok == false) && (strcmp(argv[1], "-dx") == 0)) {
      argc--;
      argv++;
      dx = atof(argv[1]);
      argc--;
      argv++;
      ok = true;
    }
    if ((ok == false) && (strcmp(argv[1], "-dy") == 0)) {
      argc--;
      argv++;
      dy = atof(argv[1]);
      argc--;
      argv++;
      ok = true;
    }
    if ((ok == false) && (strcmp(argv[1], "-dz") == 0)) {
      argc--;
      argv++;
      dz = atof(argv[1]);
      argc--;
      argv++;
      ok = true;
    }
    if ((ok == false) && (strcmp(argv[1], "-ds") == 0)) {
      argc--;
      argv++;
      dx = atof(argv[1]);
      dy = atof(argv[1]);
      dz = atof(argv[1]);
      argc--;
      argv++;
      ok = true;
    }
    if ((ok == false) && (strcmp(argv[1], "-parameter") == 0)) {
      argc--;
      argv++;
      ok = true;
      for (i = 0; i < numberOfLevels; i++) {
        parameter_name[i] = argv[1];
        argc--;
        argv++;
      }
    }
    if ((ok == false) && (strcmp(argv[1], "-debug") == 0)) {
      argc--;
      argv++;
      ok = true;
      debug = true;
    }
    if ((ok == false) && (strcmp(argv[1], "-mask") == 0)) {
      argc--;
      argv++;
      mask_name = argv[1];
      argc--;
      argv++;
      ok = true;
    }
    if (ok == false) {
      cerr << "Can not parse argument " << argv[1] << endl;
      usage();
    }
  }

  // Is there a mask to use?
  if (mask_name != NULL) {
    irtkGreyImage mask(mask_name);

    if (mask.GetX() != target.GetX() ||
        mask.GetY() != target.GetY() ||
        mask.GetZ() != target.GetZ()) {
      cerr << "Mask given does not match target dimensions." << endl;
      exit(1);
    }

    irtkGreyPixel *ptr2target, *ptr2mask;
    int voxels;

    voxels     = target.GetNumberOfVoxels();
    ptr2target = target.GetPointerToVoxels();
    ptr2mask   = mask.GetPointerToVoxels();

    for (i = 0; i < voxels; ++i) {
      if (*ptr2mask <= 0) {
        *ptr2target = paddingValue;
      }
      ++ptr2mask;
      ++ptr2target;
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

  // Print some information
  for (i = 0; i < numberOfLevels; i++) {
    cout << "Performing non-rigid registration using "
         << parameter_name[i] << " writing result to " << dofout_name[i]
         << endl;
  }

  // Create transformation
  irtkMultiLevelFreeFormTransformation *mffd =
    new irtkMultiLevelFreeFormTransformation;

  // Read transformation
  if (dofin_name != NULL) {
    irtkTransformation *transform = irtkTransformation::New(dofin_name);
    if (strcmp(transform->NameOfClass(), "irtkRigidTransformation") == 0) {
      mffd = new irtkMultiLevelFreeFormTransformation(*((irtkRigidTransformation *)transform));
    } else {
      if (strcmp(transform->NameOfClass(), "irtkAffineTransformation") == 0) {
        mffd = new irtkMultiLevelFreeFormTransformation(*((irtkAffineTransformation *)transform));
      } else {
        if (strcmp(transform->NameOfClass(), "irtkMultiLevelFreeFormTransformation") == 0) {
          mffd = new irtkMultiLevelFreeFormTransformation(*((irtkMultiLevelFreeFormTransformation *)transform));
        } else {
          cerr << "Input transformation is not of type rigid, affine " << endl;
          cerr << "or multi-level free form deformation" << endl;
          exit(1);
        }
      }
    }
    delete transform;
  } else {
    mffd = new irtkMultiLevelFreeFormTransformation;
  }

  // Create ffd
  irtkBSplineFreeFormTransformation *affd = new irtkBSplineFreeFormTransformation(target, dx, dy, dz);

  // Add ffd
  mffd->PushLocalTransformation(affd);

  i = 0;
  while (i < numberOfLevels) {
    // Perform padding
    padding(target, mffd);

    // Perform non-rigid registration
    mreg(target, source, mffd, i);

    // Next level
    i++;
  }
}


