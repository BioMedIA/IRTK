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

#include <irtkEMClassificationBiasCorrection.h>

char *output_segmentation, *output_bias, *output_corrected;

void usage()
{
  cerr << "Usage: emsbc [image] [n] [atlas 1 ... atlas n] [segmented image] [corrected image] <-sigma>" << endl;
  exit(1);
}

int main(int argc, char **argv)
{
  int i, n, ok, padding, iterations;
  bool nobg=false;
  double sigma;

  if (argc < 4) {
    usage();
    exit(1);
  }

  // Input image
  irtkRealImage image;
  image.Read(argv[1]);
  cerr<<"Input image: "<<argv[1]<<endl;
  argc--;
  argv++;

  // Number of tissues
  n = atoi(argv[1]);
  argc--;
  argv++;
  cerr<<"Number of structures: "<<atoi(argv[1])<<endl;

  // Probabilistic atlas
  irtkRealImage **atlas = new irtkRealImage*[n];
  irtkRealImage *background=NULL;

  // Read atlas for each tissue
  for (i = 0; i < n; i++) {
    atlas[i] = new irtkRealImage;
    atlas[i]->Read(argv[1]);
    cerr<<"Prior "<<i<<": "<<argv[1]<<endl;
    argc--;
    argv++;
  }

  // File name for segmentation
  output_segmentation = argv[1];
  cerr<<"Output: "<<argv[1]<<endl;
  argc--;
  argv++;

  // File name for bias corrected image
  output_corrected = argv[1];
  cerr<<"Corrected: "<<argv[1]<<endl;
  argc--;
  argv++;

  // Default parameters
  iterations = 50;
  padding    = MIN_GREY;
  sigma = 60;

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
    if ((ok == false) && (strcmp(argv[1], "-bias") == 0)) {
      argc--;
      argv++;
      output_bias = argv[1];
      argc--;
      argv++;
      ok = true;
    }
    if ((ok == false) && (strcmp(argv[1], "-background") == 0)) {
      argc--;
      argv++;
      background = new irtkRealImage;
      background->Read(argv[1]);
      cerr << "background  = " << argv[1] <<endl;
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
    if ((ok == false) && (strcmp(argv[1], "-nobg") == 0)) {
      argc--;
      argv++;
      nobg=true;
      cerr << "Not adding background tissue."<<endl;
      ok = true;
    }

    if (ok == false) {
      cerr << "Can not parse argument " << argv[1] << endl;
      usage();
    }

  }

   // Create classification object
  irtkEMClassificationBiasCorrection *classification;
  if (nobg) classification= new irtkEMClassificationBiasCorrection(n, atlas, sigma);
  else classification= new irtkEMClassificationBiasCorrection(n, atlas, background,sigma);
  classification->SetInput(image);
  classification->SetPadding(padding);
  classification->CreateMask();

  classification->MStep();
  classification->Print();
  //classification.EStep();
  //classification.WStep();
  //classification.BStep();
  /*
    for (i = 0; i < iterations; i++){
      cout << "Iteration = " << i+1 << " / " << iterations << endl;
      classification.Iterate(i);
    }
  */
  double rel_diff;
  i=0;
  do {
    cout << "Iteration = " << i+1 << " / " << iterations << endl;
    rel_diff = classification->Iterate(i);
    i++;
  } while ((rel_diff>0.001)&&(i<iterations));  //(i<iterations); 


  // Bias corrected image
  //classification.ConstructBiasCorrectedImage(image);
  //image.Write(output_biascorrection);

  // Save segmentation
  classification->ConstructSegmentation(image);
  image.Write(output_segmentation);

  // Save bias field
  //if (output_biasfield != NULL) {
    //biasfield->Write(output_biasfield);
  //}

  if (n==12) {
    classification->WriteProbMap(0,"csf.nii.gz");
    classification->WriteProbMap(1,"gray.nii.gz");
    classification->WriteProbMap(2,"caudate.nii.gz");
    classification->WriteProbMap(3,"putamen.nii.gz");
    classification->WriteProbMap(4,"nigra.nii.gz");
    classification->WriteProbMap(5,"cerebellum.nii.gz");
    classification->WriteProbMap(6,"thalamus.nii.gz");
    classification->WriteProbMap(7,"pallidum.nii.gz");
    classification->WriteProbMap(8,"brainstem.nii.gz");
    classification->WriteProbMap(9,"wm.nii.gz");
    classification->WriteProbMap(10,"callosum.nii.gz");
    classification->WriteProbMap(11,"cerebellum-white.nii.gz");
    classification->WriteProbMap(12,"other.nii.gz");
  }
  if (n==11) {
    classification->WriteProbMap(0,"csf.nii.gz");
    classification->WriteProbMap(1,"gray.nii.gz");
    classification->WriteProbMap(2,"caudate.nii.gz");
    classification->WriteProbMap(3,"putamen.nii.gz");
    classification->WriteProbMap(4,"nigra.nii.gz");
    classification->WriteProbMap(5,"cerebellum.nii.gz");
    classification->WriteProbMap(6,"thalamus.nii.gz");
    classification->WriteProbMap(7,"pallidum.nii.gz");
    classification->WriteProbMap(8,"brainstem.nii.gz");
    classification->WriteProbMap(9,"wm.nii.gz");
    classification->WriteProbMap(10,"cerebellum-white.nii.gz");
    classification->WriteProbMap(11,"other.nii.gz");
  }
  if (n==7) {
    classification->WriteProbMap(0,"csf.nii.gz");
    classification->WriteProbMap(1,"gray.nii.gz");
    classification->WriteProbMap(2,"caudate.nii.gz");
    classification->WriteProbMap(3,"putamen.nii.gz");
    classification->WriteProbMap(4,"thalamus.nii.gz");
    classification->WriteProbMap(5,"pallidum.nii.gz");
    classification->WriteProbMap(6,"white.nii.gz");
    classification->WriteProbMap(7,"other.nii.gz");
  }
  if (n==3) {
    classification->WriteProbMap(0,"csf.nii.gz");
    classification->WriteProbMap(1,"gray.nii.gz");
    classification->WriteProbMap(2,"white.nii.gz");
  }

  classification->WriteGaussianParameters("parameters.txt");

}

