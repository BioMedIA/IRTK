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

char *output_name = NULL;

void usage()
{
  cerr << "Usage: average [image] [outputfile]\n";
  exit(1);
}

int main(int argc, char **argv)
{
  double average;
  irtkRealPixel std;
  irtkRealImage *image;

  if (argc < 3) {
    usage();
  }

  image = new irtkRealImage(argv[1]);
  argc--;
  argv++;

  average = image->GetAverage();
  std = image->GetSD();

  ofstream fout(argv[1]);
  argc--;
  argv++;
  fout << average << " " << std << endl;
  fout.close();

  delete image;
}
