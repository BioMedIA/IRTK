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
  cerr << "Usage: binarize [in] [out] [threshold1] [threshold2] [value1] [value2]" << endl;
  exit(1);
}

int main(int argc, char **argv)
{
  int i, j, k, l;
  irtkBaseImage *input;
  double threshold1, threshold2, valueIn, value1, value2;

  if (argc != 7) {
    usage();
  }

  input_name  = argv[1];
  argc--;
  argv++;
  output_name = argv[1];
  argc--;
  argv++;
  threshold1 = atof(argv[1]);
  argc--;
  argv++;
  threshold2 = atof(argv[1]);
  argc--;
  argv++;
  value1 = atof(argv[1]);
  argc--;
  argv++;
  value2 = atof(argv[1]);
  argc--;
  argv++;

  // Read input
  irtkFileToImage *reader = irtkFileToImage::New(input_name);
  input = reader->GetOutput();

  for (l = 0; l < input->GetT(); ++l){
    for (k = 0; k < input->GetZ(); ++k){
      for (j = 0; j < input->GetY(); ++j){
        for (i = 0; i < input->GetX(); ++i){
          valueIn = input->GetAsDouble(i, j, k, l);
          if (valueIn >=threshold1 && valueIn <= threshold2){
            input->PutAsDouble(i, j, k, l, value1);
          } else {
            input->PutAsDouble(i, j, k, l, value2);
          }
        }
      }
    }
  }

  // Write image
  input->Write(output_name);

  return 0;
}

