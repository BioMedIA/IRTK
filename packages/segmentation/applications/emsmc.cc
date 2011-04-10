/*=========================================================================

  Library   : Image Registration Toolkit (IRTK)
  Module    : $Id: ems.cc 37 2009-04-27 15:28:58Z mm3 $
  Copyright : Imperial College, Department of Computing
              Visual Information Processing (VIP), 2008 onwards
  Date      : $Date: 2009-04-27 16:28:58 +0100 (一, 27 四月 2009) $
  Version   : $Revision: 37 $
  Changes   : $Author: mm3 $

=========================================================================*/

#include <irtkEMClassification.h>

char *output_name;

void usage()
{
  cerr << "Usage: ems [image] [n] [atlas 1 ... atlas n] [c1 ... cn cb] [output]" << endl;
  cerr << "n number of atlas" << endl;
  cerr << "ci number of components of atlas i" << endl;
  cerr << "cb number of components of background" << endl;
  exit(1);
}

int main(int argc, char **argv)
{
  int i, n, m, ok, padding, iterations;
  bool nobg=false,outputprob=false;

  if (argc < 5) {
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
  int *c = new int[n+1];
  irtkRealImage *background=NULL;

  // Read atlas for each tissue
  for (i = 0; i < n; i++) {
    atlas[i] = new irtkRealImage;
    atlas[i]->Read(argv[1]);
    cerr << "Image " << i <<" = " << argv[1] <<endl;
    argc--;
    argv++;
  }

  for (i = 0; i < n+1; i++) {
    c[i] = atoi(argv[1]);
    argc--;
    argv++;
  }

  // File name for output
  output_name = argv[1];
  argc--;
  argv++;

  // Default parameters
  iterations = 50;
  padding    = -1;//MIN_GREY;

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
      cerr << "padding  = " << padding <<endl;
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
    if ((ok == false) && (strcmp(argv[1], "-nobg") == 0)) {
      argc--;
      argv++;
      nobg=true;
      cerr << "Not adding background tissue."<<endl;
      ok = true;
    }
	if ((ok == false) && (strcmp(argv[1], "-outputprob") == 0)) {
      argc--;
      argv++;
      outputprob=true;
      cerr << "Output probability segmentations."<<endl;
      ok = true;
    }
    if (ok == false) {
      cerr << "Can not parse argument " << argv[1] << endl;
      usage();
    }
  }

  irtkEMClassificationMultiComp *classification;
  if (nobg) classification= new irtkEMClassificationMultiComp(n, atlas, c);
  else classification= new irtkEMClassificationMultiComp(n, atlas, background, c);

  classification->SetPadding(padding);
  classification->SetInput(image);
  classification->Initialise();

  double rel_diff;
  i=0;
  do {
    cout << "Iteration = " << i+1 << " / " << iterations << endl;
    rel_diff = classification->Iterate(i);
    i++;
  } while ((rel_diff>0.0001)&&(i<iterations));

  if(outputprob){
	  for (i = 0; i < n+1; i++) {
	    irtkRealImage prob;
		classification->GetProbMap(i,prob);
		char buffer[255];
		sprintf(buffer,"%s%.2d.nii.gz",output_name,i);
		prob.Write(buffer);
	  }
  }else{
	  irtkRealImage segmentation;
	  classification->ConstructSegmentation(segmentation);
	  segmentation.Write(output_name);
  }
  
  char buffer[255];
  sprintf(buffer,"%sparameters.txt",output_name);
  classification->WriteGaussianParameters(buffer,1);



  delete classification;
  delete []atlas;
  delete []c;


}

