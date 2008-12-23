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
#include <irtkEMClassificationTemplateBiasCorrection.h>
#include <irtkMultiChannelImage.h>

char *corrected_name=NULL, *weights_name=NULL, *bias_name=NULL, *parameters_name=NULL;
double cutoff=0.5, voxelsize=5, sigma = 0;
int cP=3, padding=0, iterations=5;
bool logtransformed=false;

void usage()
{
  cerr << "Usage: bc [image] [reference] [corrected] <-bias> <-iterations> <-padding> <-cutoff> <-voxelsize> <-weights> <-cP> <-sigma> <-logtransformed>" << endl;
  exit(1);
}

int main(int argc, char **argv)
{
  int ok;

  if (argc < 3) {
    usage();
    exit(1);
  }

  // Input image
  irtkRealImage image;
  irtkRealImage uncorrected;
  image.Read(argv[1]);
  uncorrected.Read(argv[1]);
  argc--;
  argv++;

  //reference image
  irtkRealImage reference;
  reference.Read(argv[1]);
  argc--;
  argv++;

  // corrected image
  corrected_name = argv[1];
  argc--;
  argv++;

  // Parse remaining parameters
  while (argc > 1) {
    ok = False;
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
    if ((ok == False) && (strcmp(argv[1], "-cutoff") == 0)) {
      argc--;
      argv++;
      cutoff = atof(argv[1]);
      argc--;
      argv++;
      ok = True;
    }

    if ((ok == False) && (strcmp(argv[1], "-voxelsize") == 0)) {
      argc--;
      argv++;
      voxelsize = atof(argv[1]);
      argc--;
      argv++;
      ok = True;
    }
    if ((ok == False) && (strcmp(argv[1], "-sigma") == 0)) {
      argc--;
      argv++;
      sigma = atof(argv[1]);
      argc--;
      argv++;
      ok = True;
    }
    if ((ok == False) && (strcmp(argv[1], "-cP") == 0)) {
      argc--;
      argv++;
      cP = atoi(argv[1]);
      argc--;
      argv++;
      ok = True;
    }
    if ((ok == False) && (strcmp(argv[1], "-weights") == 0)) {
      argc--;
      argv++;
      weights_name = argv[1];
      argc--;
      argv++;
      ok = True;
    }
    if ((ok == False) && (strcmp(argv[1], "-bias") == 0)) {
      argc--;
      argv++;
      bias_name = argv[1];
      argc--;
      argv++;
      ok = True;
    }
    if ((ok == False) && (strcmp(argv[1], "-logtransformed") == 0)) {
      argc--;
      argv++;
      logtransformed = true;
      ok = True;
    }
    if (ok == False) {
      cerr << "Can not parse argument " << argv[1] << endl;
      usage();
    }
  }

  irtkEMClassificationTemplateBiasCorrection classification(image,reference,cP,padding,voxelsize);
  classification.Initialise();
  classification.WriteProbMap(0,"inliers.nii.gz");
  classification.WriteProbMap(1,"outliers.nii.gz");
  classification.WriteInput("difference.nii.gz");

  int i = 1;
  while (i<=iterations) {

    cerr<<endl<<endl<<endl<<"Loop "<<i<<"/"<<iterations<<endl<<endl;
    classification.IterateGMM(i);
    i++;
  }
  /*
     irtkMultiChannelImage corrected;
     corrected.SetPadding(0);
     corrected.AddImage(image);
     corrected.AddImage(image);
     corrected.Log(0);
     classification.ApplyBias(corrected.GetImage(0));
     corrected.Exp(0);
     corrected.AdjustMean(0,1);
     corrected.Write(0,corrected_name);
  */
  classification.CorrectTarget();
  classification.WriteTarget(corrected_name);
  classification.WriteBias("bias");
  classification.WriteProbMap(0,"tissue0.nii.gz");
  classification.WriteProbMap(1,"tissue1.nii.gz");
}

