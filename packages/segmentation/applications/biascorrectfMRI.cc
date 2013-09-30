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
#include <irtkMultiChannelImage.h>
#include <irtkEMClassification.h>
#include <irtkEMClassificationBiasCorrectionfMRI.h>
#include <irtkErosion.h>

char *output_name;
double *mean, *var, *c;
irtkRealImage * background;
// Default parameters
double threshold = 0.001;
int iterations = 50;
int padding    = 0;
double sigma = 15;

void usage()
{
  cerr << "Usage: biascprrectfMRI [image] [output] " << endl;
  exit(1);
}

int main(int argc, char **argv)
{
  int i, n, ok;

  if (argc < 3) {
    usage();
    exit(1);
  }

  // Input image
  irtkRealImage image;
  image.Read(argv[1]);
  argc--;
  argv++;

  //output segmentation
  output_name = argv[1];
  argc--;
  argv++;

  // Number of tissues
  //n = atoi(argv[1]);
  //argc--;
  //argv++;

  n=2;
  //Default settings for Gaussian Mixture parameters
  mean = new double[n];
  var  = new double[n];
  c    = new double[n];

  irtkRealPixel min, max;
  image.GetMinMax(&min,&max);

  for (i=0;i<n;i++)  mean[i] = min + i*(max-min)/(double) n;
  for (i=0;i<n;i++)  var[i] = ((max-min)/(double) n)*((max-min)/(double) n)/16;
  for (i=0;i<n;i++)  c[i] = 1.0/(double) n;


  // Parse remaining parameters
  while (argc > 1) {
    ok = false;
    if ((ok == false) && (strcmp(argv[1], "-mean") == 0)) {
      argc--;
      argv++;

      for (i=0;i<n;i++) {
        mean[i] = atof(argv[1]);
        argc--;
        argv++;
      }
      ok = true;
    }
    if ((ok == false) && (strcmp(argv[1], "-variance") == 0)) {
      argc--;
      argv++;

      for (i=0;i<n;i++) {
        var[i] = atof(argv[1])*atof(argv[1]);
        argc--;
        argv++;
      }
      ok = true;
    }
    if ((ok == false) && (strcmp(argv[1], "-c") == 0)) {
      argc--;
      argv++;

      for (i=0;i<n;i++) {
        c[i] = atof(argv[1]);
        argc--;
        argv++;
      }
      ok = true;
    }
    if ((ok == false) && (strcmp(argv[1], "-iterations") == 0)) {
      argc--;
      argv++;
      iterations = atoi(argv[1]);
      argc--;
      argv++;
      ok = true;
    }
    if ((ok == false) && (strcmp(argv[1], "-padding") == 0)) {
      argc--;
      argv++;
      padding = atoi(argv[1]);
      argc--;
      argv++;
      ok = true;
    }
    if ((ok == false) && (strcmp(argv[1], "-threshold") == 0)) {
      argc--;
      argv++;
      threshold = atof(argv[1]);
      argc--;
      argv++;
      ok = true;
    }
    if ((ok == false) && (strcmp(argv[1], "-background") == 0)) {
      argc--;
      argv++;
      background = new irtkRealImage;
      background->Read(argv[1]);
      argc--;
      argv++;
      ok = true;
    }
    if ((ok == false) && (strcmp(argv[1], "-sigma") == 0)) {
      argc--;
      argv++;
      sigma = atof(argv[1]);
      cerr << "sigma  = " << sigma <<endl;
      argc--;
      argv++;
      ok = true;
    }

    if (ok == false) {
      cerr << "Can not parse argument " << argv[1] << endl;
      usage();
    }
  }

  //Split image into low ad high intensities
  irtkEMClassification classification;
  classification.SetInput(image);
  classification.SetPadding(padding);
  classification.CreateMask();
  classification.InitialiseGMMParameters(2,mean,var,c);

  double rel_diff;
  i=1;
  do {
    cout << "Iteration = " << i << " / " << iterations << endl;
    rel_diff = classification.IterateGMM(i,false, false);
    i++;
  } while ((rel_diff>threshold)&&(i<iterations));
  
  delete mean;
  delete var;
  delete c;
  

  //create mask for the higher intensity voxels
  irtkRealImage segmentation;
  classification.ConstructSegmentationNoBG(segmentation);
  segmentation.Write("seg.nii.gz");
  //exit(1);
  irtkMultiChannelImage mch;
  mch.AddImage(segmentation);
  irtkRealImage mask = mch.ExtractLabel(0,2);
  irtkErosion<irtkRealPixel> erosion;
  erosion.SetConnectivity(CONNECTIVITY_06);
  erosion.SetInput(&mask);
  erosion.SetOutput(&mask);
  erosion.Run();
  mask.Write("mask.nii.gz");
  mch.SetImage(0,image);
  mch.SetMask(mask);
  mch.SetPadding(0);
  mch.Brainmask();
  mch.Write(0,"fmrisignal.nii.gz");

  
  n=1;
  //Default settings for Gaussian Mixture parameters
  mean = new double[n];
  var  = new double[n];
  c    = new double[n];

  //image.GetMinMax(&min,&max);

  for (i=0;i<n;i++)  mean[i] = min + i*(max-min)/(double) n;
  for (i=0;i<n;i++)  var[i] = ((max-min)/(double) n)*((max-min)/(double) n)/16;
  for (i=0;i<n;i++)  c[i] = 1.0/(double) n;

  irtkEMClassificationBiasCorrectionfMRI classificationBC(sigma);
  classificationBC.SetInput(mch.GetImage(0));
  classificationBC.SetPadding(padding);
  classificationBC.CreateMask();
  classificationBC.InitialiseGMMParameters(1,mean,var,c);
  i=1;
  do {
    cout << "Iteration = " << i << " / " << iterations << endl;
    rel_diff = classificationBC.IterateGMM(i,false, false);
    i++;
  } while ((rel_diff>threshold)&&(i<iterations));

  classificationBC.ApplyBias(image);
  image.Write(output_name);
  classificationBC.WriteGaussianParameters("parameters.txt");

  delete mean;
  delete var;
  delete c;

}

