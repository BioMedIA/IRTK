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

// Default filenames
char *target_name = NULL, *source_name = NULL, *dof1_name = NULL, *dof2_name = NULL;

void usage()
{
  cerr << "Usage: compare [target] [source] [dof1] [dof2] <options>\n";
  cerr << "where <options> is one or more of the following:\n" << endl;
  cerr << "<-Rx1 pixel> \t Region of interest" << endl;
  cerr << "<-Ry1 pixel> \t Region of interest" << endl;
  cerr << "<-Rz1 pixel> \t Region of interest" << endl;
  cerr << "<-Rx2 pixel> \t Region of interest" << endl;
  cerr << "<-Ry2 pixel> \t Region of interest" << endl;
  cerr << "<-Rz2 pixel> \t Region of interest" << endl;
  cerr << "<-Tx1 pixel> \t Region of interest in target image" << endl;
  cerr << "<-Ty1 pixel> \t Region of interest in target image" << endl;
  cerr << "<-Tz1 pixel> \t Region of interest in target image" << endl;
  cerr << "<-Tx2 pixel> \t Region of interest in target image" << endl;
  cerr << "<-Ty2 pixel> \t Region of interest in target image" << endl;
  cerr << "<-Tz2 pixel> \t Region of interest in target image" << endl;
  cerr << "<-Sx1 pixel> \t Region of interest in source image" << endl;
  cerr << "<-Sy1 pixel> \t Region of interest in source image" << endl;
  cerr << "<-Sz1 pixel> \t Region of interest in source image" << endl;
  cerr << "<-Sx2 pixel> \t Region of interest in source image" << endl;
  cerr << "<-Sy2 pixel> \t Region of interest in source image" << endl;
  cerr << "<-Sz2 pixel> \t Region of interest in source image" << endl;
  cerr << "<-source  t> \t Threshold in source image" << endl;
  cerr << "<-target  t> \t Threshold in target image" << endl;
  cerr << "<-fast>      \t Use eight points only" << endl;
  cerr << "<-image>     \t saves an image with the error for each pixel"<<endl;
  exit(1);
}

double calculate_error(double x, double y, double z, irtkTransformation *t1,
                       irtkTransformation *t2)
{
  double x1, x2, y1, y2, z1, z2;

  x1 = x2 = x;
  y1 = y2 = y;
  z1 = z2 = z;
  t1->Transform(x1, y1, z1);
  t2->Transform(x2, y2, z2);
  return (x1 - x2) * (x1 - x2) + (y1 - y2) * (y1 - y2) + (z1 - z2) * (z1 - z2);
}

int main(int argc, char **argv)
{
  irtkTransformation *t1, *t2;
  int useImage,ok, fast, m, n, x, y, z;
  int target_threshold, source_threshold;
  int target_x1, target_y1, target_z1, target_x2, target_y2, target_z2;
  int source_x1, source_y1, source_z1, source_x2, source_y2, source_z2;
  double error, rms, max, mean, x1, y1, z1;

  char* imageName = NULL;
  useImage   = False;
  // Check command line
  if (argc < 5) {
    usage();
  }

  // Parse image
  target_name = argv[1];
  argc--;
  argv++;
  source_name = argv[1];
  argc--;
  argv++;

  // Parse gold standard transformation
  dof1_name = argv[1];
  argc--;
  argv++;
  dof2_name = argv[1];
  argc--;
  argv++;

  // Read image
  irtkGreyImage *target = new irtkGreyImage(target_name);
  irtkGreyImage *source = new irtkGreyImage(source_name);
  irtkGreyImage *image  = new irtkGreyImage;

  // Fix ROI
  target_x1 = 0;
  target_y1 = 0;
  target_z1 = 0;
  target_x2 = target->GetX();
  target_y2 = target->GetY();
  target_z2 = target->GetZ();
  source_x1 = 0;
  source_y1 = 0;
  source_z1 = 0;
  source_x2 = source->GetX();
  source_y2 = source->GetY();
  source_z2 = source->GetZ();
  target_threshold = MIN_GREY;
  source_threshold = MIN_GREY;

  // Read gold standard transformation
  t1 = irtkTransformation::New(dof1_name);

  // Read second transformation
  t2 = irtkTransformation::New(dof2_name);

  // Initialize remaining stuff
  m = 0;
  n = 0;
  fast   = False;
  while (argc > 1) {
    ok = False;
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
    if ((ok == False) && (strcmp(argv[1], "-target") == 0)) {
      argc--;
      argv++;
      target_threshold = atoi(argv[1]);
      argc--;
      argv++;
      ok = True;
    }
    if ((ok == False) && (strcmp(argv[1], "-source") == 0)) {
      argc--;
      argv++;
      source_threshold = atoi(argv[1]);
      argc--;
      argv++;
      ok = True;
    }
    if ((ok == False) && (strcmp(argv[1], "-image") == 0)) {
      argc--;
      argv++;
      useImage   = True;
      imageName  = argv[1];
      argc--;
      argv++;
      ok = True;
    }
    if ((ok == False) && (strcmp(argv[1], "-fast") == 0)) {
      argc--;
      argv++;
      fast = True;
      ok = True;
    }
    if (ok == False) {
      cerr << "Can not parse argument " << argv[1] << endl;
      usage();
    }
  }

  // If there is an region of interest for the target image, use it
  if ((target_x1 != 0) || (target_x2 != target->GetX()) ||
      (target_y1 != 0) || (target_y2 != target->GetY()) ||
      (target_z1 != 0) || (target_z2 != target->GetZ())) {
    *target = target->GetRegion(target_x1, target_y1, target_z1,
                                target_x2, target_y2, target_z2);
  }

  // If there is an region of interest for the source image, use it
  if ((source_x1 != 0) || (source_x2 != source->GetX()) ||
      (source_y1 != 0) || (source_y2 != source->GetY()) ||
      (source_z1 != 0) || (source_z2 != source->GetZ())) {
    *source = source->GetRegion(source_x1, source_y1, source_z1,
                                source_x2, source_y2, source_z2);
  }


  //Initiliaze the output image
  if (useImage==True) {
    image->Initialize(target->GetImageAttributes());
  }

  n    = 0;
  rms  = 0;
  max  = 0;
  mean = 0;
  if (fast == True) {
    for (z = 1; z <= 2; z++) {
      for (y = 1; y <= 2; y++) {
        for (x = 1; x <= 2; x++) {
          x1 = x / 3.0 * target->GetX();
          y1 = y / 3.0 * target->GetY();
          z1 = z / 3.0 * target->GetZ();
          target->ImageToWorld(x1, y1, z1);
          // Calculate difference
          error = calculate_error(x1, y1, z1, t1, t2);
          rms  += error;
          mean += sqrt(error);
          if (error > max) {
            max = error;
          }
          n++;
        }
      }
    }
  } else {

    if ((target_threshold != MIN_GREY) || (source_threshold != MIN_GREY)) {

      // Transform source image with gold standard
      irtkGreyImage tmp(*target);
      irtkImageTransformation imagetransformation;
      imagetransformation.SetInput (source, t1);
      imagetransformation.SetOutput(&tmp);
      imagetransformation.PutSourcePaddingValue(MIN_GREY);
      imagetransformation.Run();
      delete source;
      source = &tmp;

      // Calculate error, ignoring padded voxels
      for (z = 0; z < target->GetZ(); z++) {
        for (y = 0; y < target->GetY(); y++) {
          for (x = 0; x < target->GetX(); x++) {
            if ((source->Get(x, y, z) > source_threshold) &&
                (target->Get(x, y, z) > target_threshold)) {
              x1 = x;
              y1 = y;
              z1 = z;
              target->ImageToWorld(x1, y1, z1);
              // Calculate difference
              error = calculate_error(x1, y1, z1, t1, t2);
              if (useImage==True)
                image->PutAsDouble(x, y, z, sqrt(error));
              rms  += error;
              mean += sqrt(error);
              if (error > max) {
                max = error;
              }
              n++;
            }
          }
        }
      }
    } else {
      // Calculate error
      for (z = 0; z < target->GetZ(); z++) {
        for (y = 0; y < target->GetY(); y++) {
          for (x = 0; x < target->GetX(); x++) {
            x1 = x;
            y1 = y;
            z1 = z;
            target->ImageToWorld(x1, y1, z1);
            // Calculate difference
            error = calculate_error(x1, y1, z1, t1, t2);
            if (useImage==True)
              image->PutAsDouble(x, y, z, sqrt(error));
            rms  += error;
            mean += sqrt(error);
            if (error > max) {
              max = error;
            }
            n++;
          }
        }
      }
    }
  }

  // Write the image
  if (useImage==True) {
    cout<<"Output Image: "<<imageName<<" .........."<<endl;
    image->Write(imageName);
  }
  // Calculate final error
  max  = sqrt(max);
  rms  = sqrt(rms / (double)n);
  mean = mean / (double)n;

  // Print final error
  cout.setf(ios::right);
  cout.setf(ios::fixed);
  cout.precision(4);

  // Print detailed output or not
  cout << "Error is " << max  << " mm (max)" << endl;
  cout << "Error is " << rms  << " mm (rms)" << endl;
  cout << "Error is " << mean << " mm (mean)" << endl;

}















