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
#include <irtkHistogram.h>
#include <irtkMultiChannelImage.h>
#include <irtkResampling.h>
#include <irtkResamplingWithPadding.h>
#include <irtkEMClassificationBiasCorrection.h>

char *corrected_name=NULL, *weights_name=NULL, *bias_name=NULL, *parameters_name=NULL;
double cutoff=0.5, voxelsize=5, sigma = 0;
int cP=3, padding=0, iterations=5;
bool logtransformed=false;

void usage()
{
  cerr << "Usage: bc [image] [corrected] <-bias> <-iterations> <-padding> <-cutoff> <-voxelsize> <-weights> <-cP> <-sigma> <-logtransformed>" << endl;
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

  // corrected image
  corrected_name = argv[1];
  argc--;
  argv++;

  // Parse remaining parameters
  while (argc > 1) {
    ok = false;
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
    if ((ok == false) && (strcmp(argv[1], "-cutoff") == 0)) {
      argc--;
      argv++;
      cutoff = atof(argv[1]);
      argc--;
      argv++;
      ok = true;
    }

    if ((ok == false) && (strcmp(argv[1], "-voxelsize") == 0)) {
      argc--;
      argv++;
      voxelsize = atof(argv[1]);
      argc--;
      argv++;
      ok = true;
    }
    if ((ok == false) && (strcmp(argv[1], "-sigma") == 0)) {
      argc--;
      argv++;
      sigma = atof(argv[1]);
      argc--;
      argv++;
      ok = true;
    }
    if ((ok == false) && (strcmp(argv[1], "-cP") == 0)) {
      argc--;
      argv++;
      cP = atoi(argv[1]);
      argc--;
      argv++;
      ok = true;
    }
    if ((ok == false) && (strcmp(argv[1], "-weights") == 0)) {
      argc--;
      argv++;
      weights_name = argv[1];
      argc--;
      argv++;
      ok = true;
    }
    if ((ok == false) && (strcmp(argv[1], "-bias") == 0)) {
      argc--;
      argv++;
      bias_name = argv[1];
      argc--;
      argv++;
      ok = true;
    }
    if ((ok == false) && (strcmp(argv[1], "-logtransformed") == 0)) {
      argc--;
      argv++;
      logtransformed = true;
      ok = true;
    }
    if (ok == false) {
      cerr << "Can not parse argument " << argv[1] << endl;
      usage();
    }
  }

  irtkMultiChannelImage mch;
  mch.SetPadding(padding);
  mch.AddImage(image);
  if (!logtransformed) mch.Log(0);

  int i;
  irtkRealImage i3;
  cerr << "Resampling with padding ... ";
  irtkResamplingWithPadding<irtkRealPixel> resampling(voxelsize, voxelsize, voxelsize, padding);
  resampling.SetInput(&(mch.GetImage(0)));
  resampling.SetOutput(&i3);
  resampling.Run();
  cerr << "done." << endl;

  i3.Write("resampled.nii.gz");

  cerr<<"Initializing biasfield ...";
  irtkBSplineBiasField *biasfield = new irtkBSplineBiasField(i3, cP, cP, cP);
  cerr<<"done."<<endl;

  cerr << "Creating classification object ... ";
  irtkEMClassificationBiasCorrection classification;
  classification.SetInput(i3);
  classification.SetPadding(padding);
  classification.SetBiasField(biasfield);
  cerr << "done." << endl;

  cerr <<"Estimating GMM parameters ... ";
  irtkHistogram h(100);
  irtkRealPixel imin, imax;
  i3.GetMinMaxPad(&imin, &imax,padding);
  cerr<<" min = "<<imin<<", max = "<<imax<<" ... ";
  h.PutMin(imin);
  h.PutMax(imax);
  irtkRealPixel *ptr=i3.GetPointerToVoxels();
  for (i=0; i<i3.GetNumberOfVoxels(); i++) {
    if (ptr[i]!=padding) h.AddSample(ptr[3]);
    ptr++;
  }
  double mean, variance;
  mean = h.Mean();
  variance=h.Variance();
  cerr<<"mean="<<mean<<" variance="<<sqrt(variance)<<" ... done."<<endl;
  double *m=new double[2];
  double *s=new double[2];
  double *c=new double[2];
  m[0]=m[1]=mean;
  s[0]=variance/4;
  s[1]=variance*4;
  c[0]=c[1]=0.5;

  classification.InitialiseGMMParameters(2,m,s,c);
  classification.WriteProbMap(0,"inliers.nii.gz");
  classification.WriteProbMap(1,"outliers.nii.gz");

  cerr<<endl;

  delete[] m;
  delete[] s;
  delete[] c;


  i=1;
  while (i<=iterations) {

    cerr<<endl<<endl<<endl<<"Loop "<<i<<"/"<<iterations<<endl<<endl;
    classification.IterateGMM(i, false, false);
    i++;
  }
  biasfield->Write("bias");
  cerr<<"Applying bias ...";
  classification.ApplyBias((mch.GetImage(0)));
  cerr<<"done"<<endl;
  mch.Exp(0);
  mch.Write(0,corrected_name);

}

