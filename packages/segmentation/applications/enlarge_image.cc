/*=========================================================================

  Library   : packages/segmentation
  Module    : $RCSfile: gmm.cc,v $
  Authors   : Maria Murgasova and Daniel Rueckert
  Copyright : Imperial College, Department of Computing
              Visual Information Processing (VIP), 2006
  Purpose   :
  Date      : $Date: 2007/05/24 15:22:19 $

=========================================================================*/

#include <irtkImage.h>

char *output_name;
irtkRealImage image;


void usage()
{
  cerr << "Usage: enlarge_image [image] [output] <-x voxels> <-y voxels> <-z voxels>"<<endl;
  exit(1);
}

int main(int argc, char **argv)
{

  if (argc < 3) {
    usage();
    exit(1);
  }

  int ex=0,ey=0,ez=0;
  int ok;

  image.Read(argv[1]);
  argc--;
  argv++;

  output_name=argv[1];
  argc--;
  argv++;

  // Parse remaining parameters
  while (argc > 1){
    ok = false;
    if ((ok == false) && (strcmp(argv[1], "-x") == 0)){
      argc--;
      argv++;
      ex = atoi(argv[1]);
      argc--;
      argv++;
      ok = true;
    }
    if ((ok == false) && (strcmp(argv[1], "-y") == 0)){
      argc--;
      argv++;
      ey = atoi(argv[1]);
      argc--;
      argv++;
      ok = true;
    }
    if ((ok == false) && (strcmp(argv[1], "-z") == 0)){
      argc--;
      argv++;
      ez = atoi(argv[1]);
      argc--;
      argv++;
      ok = true;
    }
     if (ok == false){
      cerr << "Can not parse argument " << argv[1] << endl;
      usage();
    }
  } 


  irtkRealImage new_image; 
  irtkImageAttributes attr = image.GetImageAttributes();
  attr._x += 2*ex;
  attr._y += 2*ey;
  attr._z += 2*ez;
  new_image.Initialize(attr);
  int i,j,k;
  for(i=0; i<image.GetX();i++)
    for(j=0; j<image.GetY();j++)
      for(k=0; k<image.GetZ();k++)
      {
        new_image.Put(i+ex,j+ey,k+ez,image.Get(i,j,k));
      }
  new_image.Write(output_name);
 
}



