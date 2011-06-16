/*=========================================================================
 
  Library   : Image Registration Toolkit (IRTK)
  Module    : $Id: info.cc 8 2009-03-02 16:12:58Z dr $
  Copyright : Imperial College, Department of Computing
              Visual Information Processing (VIP), 2008 onwards
  Date      : $Date: 2009-03-02 16:12:58 +0000 (一, 02 三月 2009) $
  Version   : $Revision: 8 $
  Changes   : $Author: dr $
 
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
