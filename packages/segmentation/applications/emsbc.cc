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

#include <irtkGaussian.h>

char *output_segmentation, *output_biascorrection, *output_biasfield;

void usage()
{
  cerr << "Usage: emsbc [image] [n] [atlas 1 ... atlas n] [segmented image] [corrected image] <options>" << endl;
  exit(1);
}

int main(int argc, char **argv)
{
  int i, n, ok, padding, iterations;

  if (argc < 4) {
    usage();
    exit(1);
  }

  // Input image
  irtkRealImage image;
  image.Read(argv[1]);
  argc--;
  argv++;

  // Number of tissues
  n = atoi(argv[1]);
  argc--;
  argv++;

  // Probabilistic atlas
  irtkRealImage **atlas = new irtkRealImage*[n];
  irtkRealImage *background=NULL;

  // Read atlas for each tissue
  for (i = 0; i < n; i++) {
    atlas[i] = new irtkRealImage;
    atlas[i]->Read(argv[1]);
    argc--;
    argv++;
  }

  // File name for segmentation
  output_segmentation = argv[1];
  argc--;
  argv++;

  // File name for bias corrected image
  output_biascorrection = argv[1];
  argc--;
  argv++;

  // Default parameters
  iterations = 50;
  padding    = MIN_GREY;

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
    if ((ok == false) && (strcmp(argv[1], "-biasfield") == 0)) {
      argc--;
      argv++;
      output_biasfield = argv[1];
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

    if (ok == false) {
      cerr << "Can not parse argument " << argv[1] << endl;
      usage();
    }

  }

  // Create bias field
  //irtkBSplineBiasField *biasfield = new irtkBSplineBiasField(image, 40, 40, 40);
  irtkBSplineBiasField *biasfield = new irtkBSplineBiasField(image, 3, 3, 3);

  // Create classification object
  irtkEMClassificationBiasCorrection classification(n, atlas, background);
  classification.SetInput(image);
  classification.SetPadding(padding);
  classification.SetBiasField(biasfield);
  classification.Initialise();
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
    rel_diff = classification.Iterate(i);
    i++;
  } while ((rel_diff>0.001)&&(i<iterations));


  // Bias corrected image
  classification.ConstructBiasCorrectedImage(image);
  image.Write(output_biascorrection);

  // Save segmentation
  classification.ConstructSegmentation(image);
  image.Write(output_segmentation);

  // Save bias field
  if (output_biasfield != NULL) {
    biasfield->Write(output_biasfield);
  }

  if (n==11) {
    classification.WriteProbMap(0,"csf.hdr");
    classification.WriteProbMap(1,"gray.hdr");
    classification.WriteProbMap(2,"caudate.hdr");
    classification.WriteProbMap(3,"putamen.hdr");
    classification.WriteProbMap(4,"nigra.hdr");
    classification.WriteProbMap(5,"cerebellum.hdr");
    classification.WriteProbMap(6,"thalamus.hdr");
    classification.WriteProbMap(7,"pallidum.hdr");
    classification.WriteProbMap(8,"brainstem.hdr");
    classification.WriteProbMap(9,"white.hdr");
    classification.WriteProbMap(10,"cerebellum-white.hdr");
    classification.WriteProbMap(11,"other.hdr");
  }
  if (n==7) {
    classification.WriteProbMap(0,"csf.hdr");
    classification.WriteProbMap(1,"gray.hdr");
    classification.WriteProbMap(2,"caudate.hdr");
    classification.WriteProbMap(3,"putamen.hdr");
    classification.WriteProbMap(4,"thalamus.hdr");
    classification.WriteProbMap(5,"pallidum.hdr");
    classification.WriteProbMap(6,"white.hdr");
    classification.WriteProbMap(7,"other.hdr");
  }
  if (n==3) {
    classification.WriteProbMap(0,"csf.hdr");
    classification.WriteProbMap(1,"gray.hdr");
    classification.WriteProbMap(2,"white.hdr");
  }

  classification.WriteGaussianParameters("parameters.txt");

}

