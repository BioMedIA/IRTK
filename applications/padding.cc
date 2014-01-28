/*=========================================================================

  Library   : Image Registration Toolkit (IRTK)
  Module    : $Id$
  Copyright : Imperial College, Department of Computing
              Visual Information Processing (VIP), 2008 onwards
  Date      : $Date$
  Version   : $Revision$
  Changes   : $Author$

Copyright (c) 1999-2014 and onwards, Imperial College London
All rights reserved.
See LICENSE for details

=========================================================================*/

#include <irtkImage.h>
#include <irtkFileToImage.h>

void usage()
{
  cerr << "Usage: padding [imageA] [imageB] [output] [Value in imageB] [Padding Value in output]" << endl;
  cerr << "<-invert>         not equal Value in imageB do padding" << endl;
  exit(1);
}

int main(int argc, char **argv)
{
  int x, y, z, t, ok, invert;
  double i, j, k;
  irtkGreyPixel threshold;
  irtkRealPixel padding;
  char* outputname;

  if (argc < 6) {
    usage();
  }
  invert = 0;

  irtkFileToImage *readerA = irtkFileToImage::New(argv[1]);
  irtkBaseImage *imageA = readerA->GetOutput();
  argc--;
  argv++;
  irtkFileToImage *readerB = irtkFileToImage::New(argv[1]);
  irtkBaseImage *imageB = readerB->GetOutput();
  argc--;
  argv++;
  outputname = argv[1];
  argc--;
  argv++;
  threshold = atoi(argv[1]);
  argc--;
  argv++;
  padding = atoi(argv[1]);
  argc--;
  argv++;
  while (argc > 1) {
    ok = false;
    if ((ok == false) && (strcmp(argv[1], "-invert") == 0)) {
      argc--;
      argv++;
      invert = 1;
      ok = true;
    }
    if (ok == false) {
      cerr << "Can not parse argument " << argv[1] << endl;
      usage();
    }
  }

  for (t = 0; t < imageA->GetT(); t++) {
    for (z = 0; z < imageA->GetZ(); z++) {
      for (y = 0; y < imageA->GetY(); y++) {
        for (x = 0; x < imageA->GetX(); x++) {
          i = x; j = y; k = z;
          imageA->ImageToWorld(i,j,k);
          imageB->WorldToImage(i,j,k);
          i = round(i); j = round(j); k = round(k);
          if(i >= 0 && i < imageB->GetX()
              && j >= 0 && j < imageB->GetY()
              && k >= 0 && k < imageB->GetZ()
              && t >= 0 && t < imageB->GetT()){
            if(invert){
              if (imageB->GetAsDouble(i, j, k, t) != threshold) imageA->PutAsDouble(x, y, z, t, padding);
            }else{
              if (imageB->GetAsDouble(i, j, k, t) == threshold) imageA->PutAsDouble(x, y, z, t, padding);
            }
          }
        }
      }
    }
  }

  imageA->Write(outputname);

  return 0;
}
