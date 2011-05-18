/*=========================================================================

  Library   : packages/segmentation
  Module    : $RCSfile: gmm.cc,v $
  Authors   : Maria Murgasova and Daniel Rueckert
  Copyright : Imperial College, Department of Computing
              Visual Information Processing (VIP), 2006
  Purpose   :
  Date      : $Date$

=========================================================================*/

#include <irtkImage.h>

char *output_name;
irtkRealImage image;


void usage()
{
  cerr << "Usage: enlarge_image [image] [output] <-x voxels> <-y voxels> <-z voxels> <-percent> <-value value>"<<endl;
  exit(1);
}

int main(int argc, char **argv)
{

  if (argc < 3) {
    usage();
    exit(1);
  }

  int ex=0,ey=0,ez=0,percent=0,value=0;
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
	if ((ok == false) && (strcmp(argv[1], "-percent") == 0)){
      argc--;
      argv++;
      percent = 1;
      ok = true;
    }
	if ((ok == false) && (strcmp(argv[1], "-value") == 0)){
      argc--;
      argv++;
      value = atoi(argv[1]);
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
  if(percent == 1){
	  ex = round(double(image.GetX()*ex)/100.0);
	  ey = round(double(image.GetY()*ey)/100.0);
	  ez = round(double(image.GetZ()*ez)/100.0);
  }
  attr._x += 2*ex;
  attr._y += 2*ey;
  attr._z += 2*ez;
  new_image.Initialize(attr);
  int i,j,k,l;
  for(l=0;l<new_image.GetT();l++)
	  for(i=0; i<new_image.GetX();i++)
		  for(j=0; j<new_image.GetY();j++)
			  for(k=0; k<new_image.GetZ();k++){
			    new_image.Put(i,j,k,l,value);
			  }

  for(l=0;l<image.GetT();l++)
	  for(i=0; i<image.GetX();i++)
		  for(j=0; j<image.GetY();j++)
			  for(k=0; k<image.GetZ();k++)
			  {
				new_image.Put(i+ex,j+ey,k+ez,l,image.Get(i,j,k,l));
			  }
  new_image.Write(output_name);
 
}



