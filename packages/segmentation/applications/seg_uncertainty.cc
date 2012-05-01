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


char *output_name,*parameter_name;

void usage()
{
  cerr << "Usage: seg_uncertainty [segmentation result] [n] [prob 1 ... prob n prob background] [parameterfile] [output]" << endl;
  exit(1);
}

int main(int argc, char **argv)
{
  int i,j,k,n;

  if (argc < 6) {
    usage();
    exit(1);
  }

  irtkGreyImage segmentation;
  segmentation.Read(argv[1]);
  argc--;
  argv++;

  irtkRealImage output;
  output.Initialize(segmentation.GetImageAttributes());

  // Number of tissues
  n = atoi(argv[1]);
  argc--;
  argv++;

  // Probabilistic atlas
  irtkRealImage **atlas = new irtkRealImage*[n+1];


  // Read atlas for each tissue
  for (i = 0; i < n+1; i++) {
    atlas[i] = new irtkRealImage;
    atlas[i]->Read(argv[1]);
    cerr << "Image " << i <<" = " << argv[1] <<endl;
    argc--;
    argv++;
  }

  // File name for output
  parameter_name = argv[1];
  argc--;
  argv++;

  // File name for output
  output_name = argv[1];
  argc--;
  argv++;

  // file in file out
  ifstream fin(parameter_name);

  //Get parameters
  double *mi,*sigma;
  mi = new double[n+1];
  sigma = new double[n+1];
  for (i = 0; i < n+1; i++){
	fin >> mi[i];
	fin >> sigma[i];
	sigma[i] = sqrt(sigma[i]);
  }
  fin.close();

  //Evaluate uncertainty
  int value;
  for (i=0;i<segmentation.GetX();i++){
	  for (j=0;j<segmentation.GetY();j++){
		  for (k=0;k<segmentation.GetZ();k++){
			  value = segmentation.GetAsDouble(i,j,k);
			  if(value == 0){
				  value = n+1;
			  }
			  output.PutAsDouble(i,j,k,
				  atlas[value - 1]->GetAsDouble(i,j,k)
				  /log(sigma[value - 1]+1.00001));
		  }
	  }
  }

  //Output uncertainty
  output.Write(output_name);

  //Clear
  for (i = 0; i < n+1; i++){
	  delete atlas[i];
  }
  delete []atlas;
  delete []mi;
  delete []sigma;

}

