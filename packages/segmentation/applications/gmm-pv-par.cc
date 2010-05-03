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
bool first = true;

bool debug=false;

// specific parameters for neonatal segmentation
irtkGreyImage * non_brain = NULL;
irtkGreyImage * smask = NULL;
double lcc_treshold = 0.5;


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
  cerr<<" -pv_neonatal ... GMM fitting is followed by removal of partial volume misclassifications"<< endl;
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
    if ((ok == False) && (strcmp(argv[1], "-pv_neonatal") == 0)) {
      argc--;
      argv++;
      pv_post=true;
      ok = True;
    }
    if ((ok == False) && (strcmp(argv[1], "-debug") == 0)) {
      argc--;
      argv++;
      debug=true;
      ok = True;
    }
    if ((ok == False) && (strcmp(argv[1], "-lcc_treshold") == 0)) {
      argc--;
      argv++;
      lcc_treshold = atof(argv[1]);
      argc--;
      argv++;
      ok = True;
    }
    if ((ok == False) && (strcmp(argv[1], "-subcortical_mask") == 0)) {
      argc--;
      argv++;
      cerr<<"Reading subcortical mask: "<<argv[1]<<endl;
      smask = new irtkGreyImage;
      smask->Read(argv[1]);
      argc--;
      argv++;
      ok = True;
    }
    if ((ok == False) && (strcmp(argv[1], "-non_brain") == 0)) {
      argc--;
      argv++;
      non_brain = new irtkGreyImage;
      non_brain->Read(argv[1]);
      argc--;
      argv++;
      ok = True;
    }


    if (ok == False) {
      cerr << "Can not parse argument " << argv[1] << endl;
      usage();
    }
  }

  //Initialise
  irtkEMClassification classification;
  classification.SetDebugFlag(debug);
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
  classification.ConstructSegmentation();
  classification.WriteSegmentation("segGMM.nii.gz");
  
  if(pv_post)
  {
    classification.InitialiseAtlas();
    classification.InitialisePVSegmentation();
    classification.UniformPrior();


     i=0;
     cerr<<"Iteration "<<i+1<<": ";
     while ((i<iterations)&&(classification.PVStep(3,4,2,5,1,lcc_treshold,first,non_brain)))
   {
     cerr<<"Iteration "<<i+2<<": ";
     first = false;
     if (debug)
     {
       char name[100];
       sprintf(name, "seg%d.nii.gz",i+10);
       classification.WriteSegmentation(name);
     }
     i++;
   }
    classification.WriteSegmentation("seg-post.nii.gz");
    classification.ConstructSegmentationBrainNonBrain(3,4,2,5,1,smask);
    classification.WriteSegmentation(output_name);
    classification.WritePVProbMap(0,"bg.nii.gz");
    classification.WritePVProbMap(1,"cortex.nii.gz");
    classification.WritePVProbMap(2,"wm1.nii.gz");
    classification.WritePVProbMap(3,"wm2.nii.gz");
    classification.WritePVProbMap(4,"csf.nii.gz");
  }
}

