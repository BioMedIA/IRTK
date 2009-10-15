#/*=========================================================================

  Library   : Image Registration Toolkit (IRTK)
  Module    : $Id$
  Copyright : Imperial College, Department of Computing
              Visual Information Processing (VIP), 2009 onwards
  Date      : $Date$
  Version   : $Revision$
  Changes   : $Author$

=========================================================================*/

#include <irtkRegistration.h>

#ifdef HAS_VTK

// Default filenames
char *model_name = NULL,  *image_name = NULL;
char *dofin_name  = NULL, *dofout_name = NULL;
char *parin_name  = NULL, *parout_name = NULL;

void usage()
{
  cerr << "Usage: mnreg [model] [image] <options> \n" << endl;
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
  int ok;
  
  // Check command line
  if (argc < 3) {
    usage();
  }

  // Parse model and images
  model_name = argv[1];
  argc--;
  argv++;
  image_name = argv[1];
  argc--;
  argv++;

  // Read model
  cout << "Reading model ... " << model_name << endl;
  vtkPolyDataReader *model_reader = vtkPolyDataReader::New();
  model_reader->SetFileName(model_name);
  model_reader->Modified();
  model_reader->Update();
  vtkPolyData *model = vtkPolyData::New();
  model = model_reader->GetOutput();
  model->Update();

  // Read image
  cout << "Reading image ... "; cout.flush();
  irtkGreyImage image(image_name);
  cout << "done" << endl;

  // Create registration filter
  irtkModelFreeFormRegistration *registration = new irtkModelFreeFormRegistration;

  // Create initial multi-level free-form deformation
  irtkMultiLevelFreeFormTransformation *mffd = NULL;

  // Default parameters
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
    if ((ok == False) && (strcmp(argv[1], "-debug") == 0)) {
      argc--;
      argv++;
      ok = True;
      registration->SetDebugFlag(True);
    }
    if ((ok == False) && (strcmp(argv[1], "-ds") == 0)) {
      argc--;
      argv++;
      spacing = atof(argv[1]);
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
    if (ok == False) {
      cerr << "Can not parse argument " << argv[1] << endl;
      usage();
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
  registration->SetInput(model, &image);
  registration->SetOutput(mffd);

  // Read parameter if there any, otherwise make an intelligent guess
  if (parin_name != NULL) {
    registration->irtkModelRegistration::Read(parin_name);
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
    registration->irtkModelRegistration::Write(parout_name);
  }

  // Run registration filter
  registration->Run();

  // Write the final transformation estimate
  if (dofout_name != NULL) {
    mffd->irtkTransformation::Write(dofout_name);
  }
}

#else

#include <irtkImage.h>

int main( int argc, char *argv[] )
{
  cerr << argv[0] << " this program needs to be compiled with vtk enabled." << endl;
  return 0;
}

#endif

