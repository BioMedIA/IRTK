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

char *output_name = NULL;

void usage()
{
  cerr << "Usage: average [image] [output]\n";
  exit(1);
}

int main(int argc, char **argv)
{
  int x, y, z, t;
  double average;
  irtkRealPixel std;
  irtkRealImage *image;
  irtkRealImage *output;

  if (argc < 3) {
    usage();
  }

  image = new irtkRealImage(argv[1]);
  output = new irtkRealImage(*image);
  argc--;
  argv++;
  output_name = argv[1];
  argc--;
  argv++;

  average = image->GetAverage();
  std = image->GetSD();

  for(t = 0; t < output->GetT(); t++){
      for(z = 0; z < output->GetZ(); z++){
          for(y = 0; y < output->GetY(); y++){
              for(x = 0; x < output->GetX(); x++){
                  if(output->GetAsDouble(x,y,z,t) > average + 3.0*std){
                      output->PutAsDouble(x,y,z,t,average + 3.0*std);
                  }else if(output->GetAsDouble(x,y,z,t) < average - 3.0*std){
                      output->PutAsDouble(x,y,z,t,average - 3.0*std);
                  }
              }
          }
      }
  }

  output->Write(output_name);

  delete image;
  delete output;
}
