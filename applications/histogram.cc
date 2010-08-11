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

#include <irtkHistogram_1D.h>

void usage()
{
  cerr << "Creates a histogram for a specific image" << endl;
  cerr << "Usage: [input] [output] <options>\n" << endl;
  cerr << "where <options> is one or more of the following:\n" << endl;
  cerr << "<-Rx1 pixel>       Region of interest" << endl;
  cerr << "<-Ry1 pixel>       Region of interest" << endl;
  cerr << "<-Rz1 pixel>       Region of interest" << endl;
  cerr << "<-Rx2 pixel>       Region of interest" << endl;
  cerr << "<-Ry2 pixel>       Region of interest" << endl;
  cerr << "<-Rz2 pixel>       Region of interest" << endl;

  exit(1);
}



int main(int argc, char **argv)
{
  int ok, x1, y1, z1, x2, y2, z2;

  char *input_name = NULL;
  char *output_name = NULL;

  // Check command line
  if (argc < 2) {
    usage();
  }

  // Parse image
  input_name  = argv[1];
  argc--;
  argv++;

  // Parse output name
  output_name  = argv[1];
  argc--;
  argv++;

  // Read image
  cout << "Reading image ... "<<input_name <<endl;
  irtkGreyImage input;
  input.Read(input_name);
  cout << "done" << endl;

  // Fix ROI
  x1 = 0;
  y1 = 0;
  z1 = 0;
  x2 = input.GetX();
  y2 = input.GetY();
  z2 = input.GetZ();

  while (argc > 1) {
    ok = False;

    if ((ok == False) && (strcmp(argv[1], "-Rx1") == 0)) {
      argc--;
      argv++;
      x1 = atoi(argv[1]);
      argc--;
      argv++;
      ok = True;
    }
    if ((ok == False) && (strcmp(argv[1], "-Rx2") == 0)) {
      argc--;
      argv++;
      x2 = atoi(argv[1]);
      argc--;
      argv++;
      ok = True;
    }
    if ((ok == False) && (strcmp(argv[1], "-Ry1") == 0)) {
      argc--;
      argv++;
      y1 = atoi(argv[1]);
      argc--;
      argv++;
      ok = True;
    }
    if ((ok == False) && (strcmp(argv[1], "-Ry2") == 0)) {
      argc--;
      argv++;
      y2 = atoi(argv[1]);
      argc--;
      argv++;
      ok = True;
    }
    if ((ok == False) && (strcmp(argv[1], "-Rz1") == 0)) {
      argc--;
      argv++;
      z1 = atoi(argv[1]);
      argc--;
      argv++;
      ok = True;
    }
    if ((ok == False) && (strcmp(argv[1], "-Rz2") == 0)) {
      argc--;
      argv++;
      z2 = atoi(argv[1]);
      argc--;
      argv++;
      ok = True;
    }
    if (ok == False) {
      cerr << "Can not parse argument " << argv[1] << endl;
      usage();
    }
  }


  input = input.GetRegion(x1, y1, z1, x2, y2, z2);


  irtkGreyPixel min ;
  irtkGreyPixel max ;

  input.GetMinMax(&min,&max);


  cout<<"Creating the histogram"<<endl;
  // Create the histogram
  irtkHistogram_1D<int> *hist = new irtkHistogram_1D<int>(double(min), double(max), 1);

  for (int x=0; x < input.GetX(); ++x) {
    for (int y = 0; y < input.GetY(); y++) {
      for (int z = 0; z < input.GetZ(); z++) {
        hist->AddSample(input.Get(x,y,z));
      }
    }
  }

  // Write the histogram
  hist->Write(output_name);

  delete hist;
  return 0;

}
