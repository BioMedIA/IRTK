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
  cerr << "Dicemetric of two images from label min to lable max" << endl;
  cerr << "<-output file>       Result output file" << endl;
  cerr << "<-minvalue value>        Min counting value of label default 1" << endl;
  cerr << "<-maxvalue value>        Max counting value of label default 1" << endl;
  cerr << "<-minZ value>             Min Z slice number to evaluate default 0" << endl;
  cerr << "<-maxZ value>             Max Z slice number to evaluate default max" << endl;
  exit(1);
}

int main(int argc, char **argv)
{
  double min,max;
  int i,j,k,l,ok,minz,maxz;
  int *countAnd, *countA, *countB, fraction;

  if (argc < 3) {
    usage();
  }
  min = 1; max = 1;

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

  // Initialize
  countA = new int[int(max-min+1)];
  countB = new int[int(max-min+1)];
  countAnd = new int[int(max-min+1)];

  for(i = 0; i < max-min+1; i++){
      countA[i] = 0;
      countB[i] = 0;
      countAnd[i] = 0;
  }

  for(l = 0; l < imageA.GetT(); l++ ){
	  for(k = minz; k < maxz; k++){
		  for(j = 0; j < imageA.GetY(); j++){
			  for(i = 0; i < imageA.GetX(); i++){
				  if(imageA.GetAsDouble(i,j,k,l) >= min && imageA.GetAsDouble(i,j,k,l) <= max){
					countA[int(imageA.GetAsDouble(i,j,k,l) - min)]++;
				  }
				  if(imageB.GetAsDouble(i,j,k,l) >= min && imageB.GetAsDouble(i,j,k,l) <= max){
					countB[int(imageB.GetAsDouble(i,j,k,l) - min)]++;
				  }
				  if(imageA.GetAsDouble(i,j,k,l) ==  imageB.GetAsDouble(i,j,k,l)){
					countAnd[int(imageA.GetAsDouble(i,j,k,l) - min)]++;
				  }
			  }
		  }

	  }

      //output
      for(i = 0; i < max-min+1; i++){
          fraction = (countA[i]+countB[i])/2;
          if(fraction != 0)
              cout << "Dice metric of frame "<< l << " label " << i+min << " is " << double(countAnd[i])/double(fraction) << endl;
          else
              cout << "Dice metric of frame "<< l << " label " << i+min << " is 1"<< endl;
          if(outputname){
              cerr << "Writing Results: " << outputname << endl;
              ofstream fout(outputname,ios::app);
              if(fraction != 0)
                  fout << double(countAnd[i])/double(fraction);
              else
                  fout << 1;
              fout.close();
          }
          fraction = countA[i]+countB[i]-countAnd[i];
          if(fraction != 0)
              cout << "Overlap metric of frame "<< l << " label " << i+min << " is " << double(countAnd[i])/double(fraction) << endl;
          else
              cout << "Overlap metric of frame "<< l << " label " << i+min << " is 1"<< endl;
          if(outputname){
              cerr << "Writing Results: " << outputname << endl;
              ofstream fout(outputname,ios::app);
              if(fraction != 0)
                  fout << " " << double(countAnd[i])/double(fraction) <<endl;
              else
                  fout << " " << 1 <<endl;
              fout.close();
          }
      }
  }

  return 0;
}
