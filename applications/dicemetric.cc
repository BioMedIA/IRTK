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
char *outputname = NULL;

void usage()
{
  cerr << "Usage: dicemetric [segmentationimageA] [segmentationimageB]" << endl;
  cerr << "<-output file>       Result output file" << endl;
  cerr << "<-minvalue value>         Min counting value in image default 1" << endl;
  cerr << "<-maxvalue value>         Max counting value in image default 1" << endl;
  cerr << "<-minZ value>             Min Z slice number to evaluate default 0" << endl;
  cerr << "<-maxZ value>             Max Z slice number to evaluate default max" << endl;
  exit(1);
}

int main(int argc, char **argv)
{
  double min,max;
  int i,j,k,l,a,b,ok,minz,maxz;
  int countAnd, countA, countB, fraction;

  if (argc < 3) {
    usage();
  }
  min = 1; 
  max = 1;

  irtkGreyImage imageA(argv[1]);
  argc--;
  argv++;
  irtkGreyImage imageB(argv[1]);
  argc--;
  argv++;
  minz = 0;
  maxz = imageA.GetZ();
  if(imageA.GetX() != imageB.GetX() || imageA.GetY() != imageB.GetY() 
	  || imageA.GetZ() != imageB.GetZ() || imageA.GetT() != imageB.GetT()){
	cerr << "Image sizes of A and B do not correspond!" << endl;
  }
  while (argc > 1) {
    ok = false;
    if ((ok == false) && (strcmp(argv[1], "-output") == 0)) {
      argc--;
      argv++;
      outputname = argv[1];
	  argc--;
	  argv++;
      ok = true;
    }
	if ((ok == false) && (strcmp(argv[1], "-minvalue") == 0)) {
      argc--;
      argv++;
      min = atoi(argv[1]);
	  argc--;
	  argv++;
      ok = true;
    }
	if ((ok == false) && (strcmp(argv[1], "-maxvalue") == 0)) {
      argc--;
      argv++;
      max = atoi(argv[1]);
	  argc--;
	  argv++;
      ok = true;
    }
	if ((ok == false) && (strcmp(argv[1], "-minZ") == 0)) {
      argc--;
      argv++;
      minz = atoi(argv[1]);
	  argc--;
	  argv++;
      ok = true;
    }
	if ((ok == false) && (strcmp(argv[1], "-maxZ") == 0)) {
      argc--;
      argv++;
      maxz = atoi(argv[1]);
	  argc--;
	  argv++;
      ok = true;
    }
    if (ok == false) {
      cerr << "Can not parse argument " << argv[1] << endl;
      usage();
    }
  }
  for(l = 0; l < imageA.GetT(); l++ ){
      countAnd = 0; countA = 0; countB = 0;
	  for(k = minz; k < maxz; k++){
		  for(j = 0; j < imageA.GetY(); j++){
			  for(i = 0; i < imageA.GetX(); i++){
				  a = 0;b = 0;
				  if(imageA.GetAsDouble(i,j,k,l) >= min && imageA.GetAsDouble(i,j,k,l) <= max){
					countA++;
					a = 1;
				  }
				  if(imageB.GetAsDouble(i,j,k,l) >= min && imageB.GetAsDouble(i,j,k,l) <= max){
					countB++;
					b = 1;
				  }
				  if(a && b){
					countAnd++;
				  }
			  }
		  }

	  }
	  fraction = (countA+countB)/2;
	  if(fraction != 0)
		  cout << "Dice metric of frame "<< l << " is " << double(countAnd)/double(fraction) << endl;
	  else
		  cout << "Dice metric of frame "<< l << " is 1"<< endl;
	  if(outputname){
		  cerr << "Writing Results: " << outputname << endl;
		  ofstream fout(outputname,ios::app);
		  if(fraction != 0)
			  fout << double(countAnd)/double(fraction);
		  else
			  fout << 1;
		  fout.close();
	  }
	  fraction = countA+countB-countAnd;
	  if(fraction != 0)
		cout << "Overlap metric of frame "<< l << " is " << double(countAnd)/double(fraction) << endl;
	  else
		cout << "Overlap metric of frame "<< l << " is 1"<< endl;
	  if(outputname){
		  cerr << "Writing Results: " << outputname << endl;
		  ofstream fout(outputname,ios::app);
		  if(fraction != 0)
			  fout << " " << double(countAnd)/double(fraction) <<endl;
		  else
			  fout << " " << 1 <<endl;
		  fout.close();
	  }
  }

  return 0;
}
