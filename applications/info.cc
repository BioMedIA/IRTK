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

char *input_name = NULL;

void usage()
{
  cerr << "Usage: info [image]\n";
  exit(1);
}

int main(int argc, char **argv)
{
	double min, max;
  irtkBaseImage *image;

  if (argc < 2) {
    usage();
  }

  cout << "Information from irtkFileToImage" << endl;
  irtkFileToImage *reader = irtkFileToImage::New(argv[1]);
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

  return 0;
}
