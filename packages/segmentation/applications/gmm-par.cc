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
#include <irtkEMClassification.h>

char *output_name;
double *mean, *var, *c;
irtkRealImage * background;
// Default parameters
double treshold = 0.0001;
int iterations = 50;
int padding    = 0;

void usage()
{
  cerr << "Usage: gmm [image] [output] [n] <-mean mean1 ... meann> <-stdev s1 ... sn> <-c c1 ... cn> " << endl;
  cerr << "<-iterations> <-padding> <-background> <-treshold>" << endl;
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
  n = atoi(argv[1]);
  argc--;
  argv++;

  //Default settings for Gaussian Mixture parameters
  mean = new double[n];
  var  = new double[n];
  c    = new double[n];

  irtkRealPixel min, max;
  image.GetMinMax(&min,&max);

  for (i=0;i<n;i++)  mean[i] = min + i*(max-min)/(double) n;
  for (i=0;i<n;i++)  var[i] = ((max-min)/(double) n)*((max-min)/(double) n);
  for (i=0;i<n;i++)  c[i] = 1.0/(double) n;


  // Parse remaining parameters
  while (argc > 1) {
    ok = False;
    if ((ok == False) && (strcmp(argv[1], "-mean") == 0)) {
      argc--;
      argv++;

      for (i=0;i<n;i++) {
        mean[i] = atof(argv[1]);
        argc--;
        argv++;
      }
      ok = True;
    }
    if ((ok == False) && (strcmp(argv[1], "-variance") == 0)) {
      argc--;
      argv++;

      for (i=0;i<n;i++) {
        var[i] = atof(argv[1])*atof(argv[1]);
        argc--;
        argv++;
      }
      ok = True;
    }
    if ((ok == False) && (strcmp(argv[1], "-c") == 0)) {
      argc--;
      argv++;

      for (i=0;i<n;i++) {
        c[i] = atof(argv[1]);
        argc--;
        argv++;
      }
      ok = True;
    }
    if ((ok == False) && (strcmp(argv[1], "-iterations") == 0)) {
      argc--;
      argv++;
      iterations = atoi(argv[1]);
      argc--;
      argv++;
      ok = True;
    }
    if ((ok == False) && (strcmp(argv[1], "-padding") == 0)) {
      argc--;
      argv++;
      padding = atoi(argv[1]);
      argc--;
      argv++;
      ok = True;
    }
    if ((ok == False) && (strcmp(argv[1], "-treshold") == 0)) {
      argc--;
      argv++;
      treshold = atof(argv[1]);
      argc--;
      argv++;
      ok = True;
    }
    if ((ok == False) && (strcmp(argv[1], "-background") == 0)) {
      argc--;
      argv++;
      background = new irtkRealImage;
      background->Read(argv[1]);
      argc--;
      argv++;
      ok = True;
    }
    if (ok == False) {
      cerr << "Can not parse argument " << argv[1] << endl;
      usage();
    }
  }

  irtkEMClassification classification;
  classification.SetInput(image);
  classification.SetPadding(padding);
  classification.CreateMask();
  classification.InitialiseGMMParameters(n,mean,var,c);

  double rel_diff;
  i=1;
  do {
    cout << "Iteration = " << i << " / " << iterations << endl;
    rel_diff = classification.IterateGMM(i);
    i++;
  } while ((rel_diff>treshold)&&(i<iterations));


  classification.WriteGaussianParameters("parameters.txt");

  irtkRealImage segmentation;
  classification.ConstructSegmentationNoBG(segmentation);
  segmentation.Write(output_name);
  if (n==11) {
    classification.WriteProbMap(0,"csf.nii.gz");
    classification.WriteProbMap(1,"gray.nii.gz");
    classification.WriteProbMap(2,"caudate.nii.gz");
    classification.WriteProbMap(3,"putamen.nii.gz");
    classification.WriteProbMap(4,"nigra.nii.gz");
    classification.WriteProbMap(5,"cerebellum.nii.gz");
    classification.WriteProbMap(6,"thalamus.nii.gz");
    classification.WriteProbMap(7,"pallidum.nii.gz");
    classification.WriteProbMap(8,"brainstem.nii.gz");
    classification.WriteProbMap(9,"white.nii.gz");
    classification.WriteProbMap(10,"cerebellum-white.nii.gz");
    classification.WriteProbMap(11,"other.nii.gz");
  }
  if (n==5) {
    classification.WriteProbMap(0,"caudate.nii.gz");
    classification.WriteProbMap(1,"putamen.nii.gz");
    classification.WriteProbMap(2,"thalamus.nii.gz");
    classification.WriteProbMap(3,"pallidum.nii.gz");
    classification.WriteProbMap(4,"white.nii.gz");
  }
  if (n==2) {
    classification.WriteProbMap(0,"gray.nii.gz");
    classification.WriteProbMap(1,"white.nii.gz");
  }

  if (n==3) {
    classification.WriteProbMap(0,"csf.nii.gz");
    classification.WriteProbMap(1,"gray.nii.gz");
    classification.WriteProbMap(2,"white.nii.gz");
  }

}

