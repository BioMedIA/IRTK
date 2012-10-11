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

#include <irtkFileToImage.h>

char *input_name = NULL, *output_name = NULL;

void usage()
{
  cerr << "Usage: info [image]\n";
  cerr << "<-time> [filename] Output overall time to a file" << endl;
  exit(1);
}

int main(int argc, char **argv)
{
  double min, max;
  int ok;
  irtkBaseImage *image;

  if (argc < 2) {
    usage();
  }

  cout << "Information from irtkFileToImage" << endl;
  irtkFileToImage *reader = irtkFileToImage::New(argv[1]);
  argc--;
  argv++;
  while (argc > 1) {
      ok = false;
      if ((ok == false) && (strcmp(argv[1], "-time") == 0)) {
          argc--;
          argv++;
          output_name = argv[1];
          argc--;
          argv++;
          ok = true;
      }
      if (ok == false) {
          cerr << "Can not parse argument " << argv[1] << endl;
          usage();
      }
  }
  reader->Print();
  
  cout << "Information from irtkBaseImage" << endl;
  image = reader->GetOutput();
  image->Print();
  image->GetMinMaxAsDouble(&min, &max);
  cout << min << " " << max << endl;
  cout << "Image to world matrix" << endl;
  image->GetImageToWorldMatrix().Print();
  cout << "World to image matrix" << endl;
  image->GetWorldToImageMatrix().Print();

  if(output_name){
	  ofstream fout(output_name,ios::app);
	  fout << image->GetTSize()*image->GetT() << " ";
	  fout.close();
  }

  return 0;
}
