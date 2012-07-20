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
  cerr << "Usage: autocontrast [image] [output]\n";
  cerr << "<options>\n" << endl;
  cerr << "where <options> is one or more of the following:\n" << endl;
  cerr << "<-Tp value>            Target padding value all pixel < value = value" << endl;
  cerr << "<-Tp2 maxvalue>        Target padding maxium value all pixel > maxvalue = maxvalue" << endl;
  exit(1);
}

int main(int argc, char **argv)
{
  int x, y, z, t, ok;
  double average,padding1,padding2;
  irtkRealPixel std;
  irtkRealImage *image;
  irtkRealImage *output;

  if (argc < 3) {
    usage();
  }

  padding1 = MIN_GREY;
  padding2 = MAX_GREY;

  image = new irtkRealImage(argv[1]);
  output = new irtkRealImage(*image);
  argc--;
  argv++;
  output_name = argv[1];
  argc--;
  argv++;

  while (argc > 1) {
      ok = false;
      if ((ok == false) && (strcmp(argv[1], "-Tp") == 0)) {
          argc--;
          argv++;
          padding1 = atoi(argv[1]);
          argc--;
          argv++;
          ok = true;
      }
      if ((ok == false) && (strcmp(argv[1], "-Tp2") == 0)) {
          argc--;
          argv++;
          padding2 = atoi(argv[1]);
          argc--;
          argv++;
          ok = true;
      }
      if (ok == false) {
          cerr << "Can not parse argument " << argv[1] << endl;
          usage();
      }
  }

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
                  if(output->GetAsDouble(x,y,z,t) < padding1)
                      output->PutAsDouble(x,y,z,t,padding1);
                  if(output->GetAsDouble(x,y,z,t) > padding2)
                      output->PutAsDouble(x,y,z,t,padding2);
              }
          }
      }
  }

  output->Write(output_name);

  delete image;
  delete output;
}
