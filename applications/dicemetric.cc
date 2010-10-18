/*=========================================================================

  Library   : Image Registration Toolkit (IRTK)
  Module    : $Id: padding.cc 2 2008-12-23 12:40:14Z dr $
  Copyright : Imperial College, Department of Computing
              Visual Information Processing (VIP), 2008 onwards
  Date      : $Date: 2008-12-23 12:40:14 +0000 (�? 23 十二�?2008) $
  Version   : $Revision: 2 $
  Changes   : $Author: dr $

=========================================================================*/

#include <irtkImage.h>
char *outputname = NULL;

void usage()
{
  cerr << "Usage: dicemetric [segmentationimageA] [segmentationimageB]" << endl;
  cerr << "<-output file>       Result output file" << endl;
  cerr << "<-minvalueA value>        Min counting value in imageA default 1" << endl;
  cerr << "<-minvalueB value>        Min counting value in imageB default 1" << endl;
  cerr << "<-maxvalueA value>        Max counting value in imageA default 1" << endl;
  cerr << "<-maxvalueB value>        Max counting value in imageB default 1" << endl;
  cerr << "<-minZ value>             Min Z slice number to evaluate default 0" << endl;
  cerr << "<-maxZ value>             Max Z slice number to evaluate default max" << endl;
  exit(1);
}

int main(int argc, char **argv)
{
  double minA,minB,maxA,maxB;
  int i,j,k,l,a,b,ok,minz,maxz;
  int countAnd, countA, countB, fraction;

  if (argc < 3) {
    usage();
  }
  minA = 1; minB = 1; maxA = 1; maxB = 1;

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
	if ((ok == false) && (strcmp(argv[1], "-minvalueA") == 0)) {
      argc--;
      argv++;
      minA = atoi(argv[1]);
	  argc--;
	  argv++;
      ok = true;
    }
	if ((ok == false) && (strcmp(argv[1], "-minvalueB") == 0)) {
      argc--;
      argv++;
      minB = atoi(argv[1]);
	  argc--;
	  argv++;
      ok = true;
    }
	if ((ok == false) && (strcmp(argv[1], "-maxvalueA") == 0)) {
      argc--;
      argv++;
      maxA = atoi(argv[1]);
	  argc--;
	  argv++;
      ok = true;
    }
	if ((ok == false) && (strcmp(argv[1], "-maxvalueB") == 0)) {
      argc--;
      argv++;
      maxB = atoi(argv[1]);
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
				  if(imageA.GetAsDouble(i,j,k,l) >= minA && imageA.GetAsDouble(i,j,k,l) <= maxA){
					countA++;
					a = 1;
				  }
				  if(imageB.GetAsDouble(i,j,k,l) >= minB && imageB.GetAsDouble(i,j,k,l) <= maxB){
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
