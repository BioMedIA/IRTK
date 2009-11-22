
#include <irtkImage.h>
#include <irtkEMClassificationTemplateBiasCorrection.h>
#include <irtkMultiChannelImage.h>

char *corrected_name=NULL, *bias_name=NULL;
double nomatch = false, voxelsize=3, sigma = 0;
int padding=0, iterations=5;
double spacing = 75;
double variance_ratio = 4;
double c0 = 0.5;


void usage()
{
  cerr << "Usage: bc2 [image] [reference] [corrected] "
		  "<-bias> <-iterations> <-padding> <-nomatch> <-voxelsize> <-spacing>" << endl;
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
  while (argc > 1){
    ok = False;
    if ((ok == False) && (strcmp(argv[1], "-iterations") == 0)){
      argc--;
      argv++;
      iterations = atoi(argv[1]);
      argc--;
      argv++;
      ok = True;
    }
    if ((ok == False) && (strcmp(argv[1], "-padding") == 0)){
      argc--;
      argv++;
      padding = atoi(argv[1]);
      argc--;
      argv++;
      ok = True;
    }
    if ((ok == False) && (strcmp(argv[1], "-nomatch") == 0)){
      argc--;
      argv++;
      nomatch = true;
      ok = True;
    }

    if ((ok == False) && (strcmp(argv[1], "-voxelsize") == 0)){
      argc--;
      argv++;
      voxelsize = atof(argv[1]);
      argc--;
      argv++;
      ok = True;
    }
   if ((ok == False) && (strcmp(argv[1], "-spacing") == 0)){
      argc--;
      argv++;
      spacing = atof(argv[1]);
      argc--;
      argv++;
      ok = True;
    }
    if ((ok == False) && (strcmp(argv[1], "-bias") == 0)){
      argc--;
      argv++;
      bias_name = argv[1];
      argc--;
      argv++;
      ok = True;
    }
    if (ok == False){
      cerr << "Can not parse argument " << argv[1] << endl;
      usage();
    }
  }

  irtkEMClassificationTemplateBiasCorrection classification(image,reference,spacing,padding,voxelsize);
  classification.Initialise(nomatch);
  //classification.WriteProbMap(0,"inliers.nii.gz");
  //classification.WriteProbMap(1,"outliers.nii.gz");
  //classification.WriteInput("difference.nii.gz");


  int i=1;
  double rel_diff;
  double f=0,f0=0;

  do
  {
    cout << "Iteration = " << i << " / " << iterations << endl;
    rel_diff = classification.IterateGMM(i);

    char buffer1[100];
    char buffer2[100];
    char buffer3[100];
    char buffer4[100];
    char buffer5[100];
    char buffer6[100];
 
    classification.CorrectTarget();
    sprintf(buffer1, "c%d.nii.gz",i);
    sprintf(buffer2, "b%d",i);
    sprintf(buffer3, "pm0%d.nii.gz",i);
    sprintf(buffer4, "pm1%d.nii.gz",i);
    sprintf(buffer5, "in%d.nii.gz",i);
    sprintf(buffer6, "w%d.nii.gz",i);

    classification.WriteTarget(buffer1);
    classification.WriteBias(buffer2);
    classification.WriteProbMap(0,buffer3);
    classification.WriteProbMap(1,buffer4);
    classification.WriteInput(buffer5);
    classification.WriteWeights(buffer6);
 

    i++;

  }
  while((rel_diff>0.001)&&(i<iterations));
/*
  i = 1;
  while(i<=iterations)
  {

    cerr<<endl<<endl<<endl<<"Loop "<<i<<"/"<<iterations<<endl<<endl;
    classification.IterateGMM(i);


    char buffer1[100];
    char buffer2[100];
    char buffer3[100];
    char buffer4[100];
    char buffer5[100];
    char buffer6[100];
    char buffer7[100];
    classification.CorrectTarget();
    sprintf(buffer1, "c%d.nii.gz",i);
    sprintf(buffer2, "b%d",i);
    sprintf(buffer3, "pm0%d.nii.gz",i);
    sprintf(buffer4, "pm1%d.nii.gz",i);
    sprintf(buffer5, "in%d.nii.gz",i);
    sprintf(buffer6, "w%d.nii.gz",i);
    //sprintf(buffer7, "rm%d.nii.gz",i);
    classification.WriteTarget(buffer1);
    classification.WriteBias(buffer2);
    classification.WriteProbMap(0,buffer3);
    classification.WriteProbMap(1,buffer4);
    classification.WriteInput(buffer5);
    classification.WriteWeights(buffer6);
    //classification.WriteMatchedReference(buffer7);

    i++;
  }
*/
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
   if (bias_name!=0) classification.WriteBias(bias_name);
   //classification.WriteProbMap(0,"tissue0.nii.gz");
   //classification.WriteProbMap(1,"tissue1.nii.gz");
}

