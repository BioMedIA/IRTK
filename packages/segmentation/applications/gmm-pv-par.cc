/*=========================================================================

  Library   : Image Registration Toolkit (IRTK)
  Module    : $Id: gmm-par.cc 2 2008-12-23 12:40:14Z dr $
  Copyright : Imperial College, Department of Computing
              Visual Information Processing (VIP), 2008 onwards
  Date      : $Date: 2008-12-23 12:40:14 +0000 (Tue, 23 Dec 2008) $
  Version   : $Revision: 2 $
  Changes   : $Author: dr $

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

bool equal_variance = false;
bool uniform_prior = false;
bool uniform_prior_last_iter = false;
bool pv_post = false;

void usage()
{
  cerr << "Usage: gmm-pv-par [image] [output] [n] <-mean mean1 ... meann> <-stdev s1 ... sn> <-c c1 ... cn> " << endl;
  cerr<<"<-stdev_equal s> <-uniform_prior> <uniform_prior_last_iter>"<<endl;
  cerr << "<-iterations i> <-padding p> <-background image> <-threshold t>" << endl;
  cerr<<endl;
  cerr<<"Fits mixture of Gaussians to an image with optional partial volume removal."<< endl;
  cerr<<endl;
  cerr<<"Options:"<< endl;
  cerr<<" -stdev_equal  ...  variances (s1 ... sn) are fixed to be equal for all classes"<< endl;
  cerr<<" -uniform_prior ... mixing proportions (c1 ... cn) are fixed to be equal"<< endl;
  cerr<<" -uniform_prior_last_iter ... only last iteration is calculated with uniform priors"<< endl;
  cerr<<" -pv_postprocess ... GMM fitting is followed by removal of partial volume misclassifications"<< endl;
  cerr<<"       under construction ..."<< endl;
  cerr<<""<< endl;
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
    if ((ok == False) && (strcmp(argv[1], "-stdev") == 0)) {
      argc--;
      argv++;

      for (i=0;i<n;i++) {
        var[i] = atof(argv[1])*atof(argv[1]);
        argc--;
        argv++;
      }
      ok = True;
    }

    if ((ok == False) && (strcmp(argv[1], "-stdev_equal") == 0)) {
      argc--;
      argv++;

      for (i=0;i<n;i++) {
        var[i] = atof(argv[1])*atof(argv[1]);
      }
      equal_variance = true;
      argc--;
      argv++;

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
    if ((ok == False) && (strcmp(argv[1], "-uniform_prior") == 0)) {
      argc--;
      argv++;
      uniform_prior = true;
      ok = True;
    }
    if ((ok == False) && (strcmp(argv[1], "-uniform_prior_last_iter") == 0)) {
      argc--;
      argv++;
      uniform_prior_last_iter = true;
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
    if ((ok == False) && (strcmp(argv[1], "-threshold") == 0)) {
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
    if ((ok == False) && (strcmp(argv[1], "-pv_postprocess") == 0)) {
      argc--;
      argv++;
      pv_post=true;
      ok = True;
    }
    if (ok == False) {
      cerr << "Can not parse argument " << argv[1] << endl;
      usage();
    }
  }

  //Initialise
  irtkEMClassification classification;
  classification.SetInput(image);
  classification.SetPadding(padding);
  classification.CreateMask();
  classification.InitialiseGMMParameters(n,mean,var,c);

  //Fit mixture of Gaussians to the image

  double rel_diff;
  i=1;
  do {
    cout << "Iteration = " << i << " / " << iterations << endl;
    rel_diff = classification.IterateGMM(i,equal_variance, uniform_prior);
    i++;
  } while ((rel_diff>treshold)&&(i<iterations));

  if (uniform_prior_last_iter)
    classification.EStepGMM(true);

  classification.WriteGaussianParameters("parameters.txt");

  irtkRealImage segmentation;
  classification.ConstructSegmentationNoBG(segmentation);
  segmentation.Write(output_name);
  
  if(pv_post)
  {
    classification.InitialiseAtlas();
    classification.UniformPrior();
  }
}

