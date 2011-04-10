/*=========================================================================

  Library   : Image Registration Toolkit (IRTK)
  Module    : $Id: padding.cc 2 2008-12-23 12:40:14Z dr $
  Copyright : Imperial College, Department of Computing
              Visual Information Processing (VIP), 2008 onwards
  Date      : $Date: 2008-12-23 12:40:14 +0000 (‰∫? 23 ÂçÅ‰∫åÊú?2008) $
  Version   : $Revision: 2 $
  Changes   : $Author: dr $

=========================================================================*/

#include <irtkImage.h>

void usage()
{
  cerr << "Usage: temporalalign [imageA] [imageB] [output]" << endl;
  cerr << "<-St1>                   Start time for imageA, must" << endl;
  cerr << "<-St2>                   Start time for imageB, must" << endl;
  cerr << "<-Et1>                   End time for imageA, must" << endl;
  cerr << "<-Et2>                   End time for imageB, must" << endl;
  exit(1);
}

int main(int argc, char **argv)
{
  int x, y, z, t, ok;
  int i, l;
  double st[2],et[2],gap[2];
  double *temporal[2];
  char* outputname;

  if (argc < 8) {
    usage();
  }

  irtkGreyImage imageA(argv[1]);
  argc--;
  argv++;
  irtkGreyImage imageB(argv[1]);
  argc--;
  argv++;
  outputname = argv[1];
  argc--;
  argv++;

  //generate output image
  irtkImageAttributes atr = imageB.GetImageAttributes();
  atr._t = imageA.GetT();
  irtkGreyImage output(atr);

  //check number of frames.  
  if(imageA.GetT() == 1 || imageB.GetT() == 1){
	  cerr << "Temproal alignment not needed: number of frame is one!" << endl;
	  exit(1);
  }

  while (argc > 1) {
    ok = false;
    if ((ok == false) && (strcmp(argv[1], "-St1") == 0)) {
      argc--;
      argv++;
      st[0] = atoi(argv[1]);
	  argc--;
      argv++;
      ok = true;
    }
	if ((ok == false) && (strcmp(argv[1], "-St2") == 0)) {
      argc--;
      argv++;
      st[1] = atoi(argv[1]);
	  argc--;
      argv++;
      ok = true;
    }
	if ((ok == false) && (strcmp(argv[1], "-Et1") == 0)) {
      argc--;
      argv++;
      et[0] = atoi(argv[1]);
	  argc--;
      argv++;
      ok = true;
    }
	if ((ok == false) && (strcmp(argv[1], "-Et2") == 0)) {
      argc--;
      argv++;
      et[1] = atoi(argv[1]);
	  argc--;
      argv++;
      ok = true;
    }
    if (ok == false) {
      cerr << "Can not parse argument " << argv[1] << endl;
      usage();
    }
  }

  //evaluate temporal info for image A and image B
  gap[0] = (et[0] - st[0])/(double(imageA.GetT())-1.0);
  gap[1] = (et[1] - st[1])/(double(imageB.GetT())-1.0);

  temporal[0] = new double[imageA.GetT()];
  temporal[1] = new double[imageB.GetT()];

  //generate temporal info
  for(i = 0; i < imageA.GetT(); i ++){
	  temporal[0][i] = st[0] + i * gap[0];
  }
  for(i = 0; i < imageB.GetT(); i ++){
	  temporal[1][i] = st[1] + i * gap[1];
  }

  for(t = 0; t < imageA.GetT(); t++){
	  //find closest time frame (l) in image B
	  double distance = 10000;
	  l = 0;
	  for(i = 0; i < imageB.GetT(); i++){
		  if(fabs(temporal[1][i] - temporal[0][t]) < distance){
			distance = fabs(temporal[1][i] - temporal[0][t]);
			l = i;
		  }
	  }
	  if(distance <= gap[1] || distance <= gap[0]){
		  for (z = 0; z < imageB.GetZ(); z++) {
			  for (y = 0; y < imageB.GetY(); y++) {
				  for (x = 0; x < imageB.GetX(); x++) {
					  output.PutAsDouble(x,y,z,t,imageB.GetAsDouble(x,y,z,l));
				  }
			  }
		  }
	  }else{
		  for (z = 0; z < imageB.GetZ(); z++) {
			  for (y = 0; y < imageB.GetY(); y++) {
				  for (x = 0; x < imageB.GetX(); x++) {
					  output.PutAsDouble(x,y,z,t,0);
				  }
			  }
		  }
	  }
  }
  output.Write(outputname);

  return 0;
}
