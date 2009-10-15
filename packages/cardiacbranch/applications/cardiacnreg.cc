/*=========================================================================

  Library   : Image Registration Toolkit (IRTK)
  Module    : $Id: nreg.cc 2 2008-12-23 12:40:14Z dr $
  Copyright : Imperial College, Department of Computing
              Visual Information Processing (VIP), 2008 onwards
  Date      : $Date: 2008-12-23 12:40:14 +0000 (‰∫? 23 ÂçÅ‰∫åÊú?2008) $
  Version   : $Revision: 2 $
  Changes   : $Author: dr $

=========================================================================*/

#include <irtkRegistration.h>

// Default filenames
char *source_name = NULL, *target_name = NULL;
char *tsource_name = NULL, *ttarget_name=NULL;
char *dofin_name  = NULL, *dofout_name = NULL;
char *parin_name  = NULL, *parout_name = NULL;
char *mask_name = NULL;

void usage()
{
  cerr << "Usage: nreg [untaggedtarget] [untaggedsource] [taggedtarget] [taggedsource] <options> \n" << endl;
  cerr << "where <options> is one or more of the following:\n" << endl;
  cerr << "<-parin file>        Read parameter from file" << endl;
  cerr << "<-parout file>       Write parameter to file" << endl;
  cerr << "<-dofin  file>       Read transformation from file" << endl;
  cerr << "<-dofout file>       Write transformation to file" << endl;
  cerr << "<-Tp  value>         Padding value in target image" << endl;
  cerr << "<-ds  value>         Initial control point spacing" << endl;
  cerr << "<-debug>             Enable debugging information" << endl;
  cerr << "<-mask file>         Use a mask to define the ROI. The mask" << endl;
  cerr << "                     must have the same dimensions as the target." << endl;
  cerr << "                     Voxels in the mask with zero or less are " << endl;
  cerr << "                     padded in the target." << endl;
  exit(1);
}

int main(int argc, char **argv)
{
  double spacing;
  int ok, padding;
  // Check command line
  if (argc < 5) {
    usage();
  }

  // Parse source and target images
  target_name = argv[1];
  argc--;
  argv++;
  source_name = argv[1];
  argc--;
  argv++;
  ttarget_name = argv[1];
  argc--;
  argv++;
  tsource_name = argv[1];
  argc--;
  argv++;


  // Read target image
  cout << "Reading untaggedtarget ... "; cout.flush();
  irtkGreyImage target(target_name);
  cout << "done" << endl;
  // Read source image
  cout << "Reading untaggedsource ... "; cout.flush();
  irtkGreyImage source(source_name);
  cout << "done" << endl;

   // Read ttarget image
  cout << "Reading taggedtarget ... "; cout.flush();
  irtkGreyImage ttarget(ttarget_name);
  cout << "done" << endl;
  // Read tsource image
  cout << "Reading taggedsource ... "; cout.flush();
  irtkGreyImage tsource(tsource_name);
  cout << "done" << endl;


  // Create registration filter
  irtkCardiacImageFreeFormRegistration *registration = NULL;
  registration = new irtkCardiacImageFreeFormRegistration;

  // Create initial multi-level free-form deformation
  irtkMultiLevelFreeFormTransformation *mffd = NULL;

  // Default parameters
  padding   = MIN_GREY;
  spacing   = 0;

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
    
    if ((ok == False) && ((strcmp(argv[1], "-parameter") == 0) || (strcmp(argv[1], "-parin") == 0))) {
      argc--;
      argv++;
      ok = True;
      parin_name = argv[1];
      argc--;
      argv++;
    }
    if ((ok == False) && (strcmp(argv[1], "-parout") == 0)) {
      argc--;
      argv++;
      ok = True;
      parout_name = argv[1];
      argc--;
      argv++;
    }
    if ((ok == False) && (strcmp(argv[1], "-mask") == 0)) {
      argc--;
      argv++;
      mask_name = argv[1];
      argc--;
      argv++;
      ok = True;
    }
    if (ok == False) {
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
    int voxels, i;

    voxels     = target.GetNumberOfVoxels();
    ptr2target = target.GetPointerToVoxels();
    ptr2mask   = mask.GetPointerToVoxels();

    for (i = 0; i < voxels; ++i) {
      if (*ptr2mask <= 0) {
        *ptr2target = padding;
      }
      ++ptr2mask;
      ++ptr2target;
    }
  }

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
          cerr << "Input transformation is not of type rigid or affine " << endl;
          cerr << "or multi-level free form deformation" << endl;
          exit(1);
        }
      }
    }
    delete transform;
  } else {
    // Otherwise use identity transformation to start
    mffd = new irtkMultiLevelFreeFormTransformation;
  }

  // Set input and output for the registration filter
  registration->SetInput(&target, &source, &ttarget, &tsource);
  registration->SetOutput(mffd);

  // Read parameter if there any, otherwise make an intelligent guess
  if (parin_name != NULL) {
    registration->irtkImageRegistration::Read(parin_name);
  } else {
    registration->GuessParameter();
  }

  // Override parameter settings if necessary
  if (spacing > 0) {
    registration->SetDX(spacing);
    registration->SetDY(spacing);
    registration->SetDZ(spacing);
  }

  // Write parameters if necessary
  if (parout_name != NULL) {
    registration->irtkImageRegistration::Write(parout_name);
  }

  // Run registration filter
  registration->Run();

  // Write the final transformation estimate
  if (dofout_name != NULL) {
    mffd->irtkTransformation::Write(dofout_name);
  }
}
