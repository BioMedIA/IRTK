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

#include <irtkDilation.h>
#include <irtkErosion.h>

char *input_name = NULL, *output_name = NULL;

void usage()
{
  cerr << "Usage: closing [in] [out] <options>\n";
  cerr << "Where <options> are one or more of the following:\n";
  cerr << "\t<-iterations n>    Number of iterations\n";
  exit(1);
}

int main(int argc, char **argv)
{
  bool ok;
  int i, iterations;
  irtkGreyImage image;

  // Check command line
  if (argc < 3) {
    usage();
  }

  // Read input and output names
  input_name  = argv[1];
  argc--;
  argv++;
  output_name = argv[1];
  argc--;
  argv++;

  // Read input
  image.Read(input_name);

  // Parse remaining parameters
  iterations = 1;
  while (argc > 1) {
    ok = false;
    if ((ok == false) && (strcmp(argv[1], "-iterations") == 0)) {
      argc--;
      argv++;
      iterations = atoi(argv[1]);
      argc--;
      argv++;
      ok = true;
    }
    if (ok == false) {
      cerr << "Unknown option: " << argv[1] << endl;
      usage();
    }
  }

  cout << "Closing ... "; cout.flush();
  irtkDilation<irtkGreyPixel> dilation;
  dilation.SetInput(&image);
  dilation.SetOutput(&image);
  irtkErosion<irtkGreyPixel> erosion;
  erosion.SetInput(&image);
  erosion.SetOutput(&image);
  for (i = 0; i < iterations; i++) {
    dilation.Run();
  } 
  for (i = 0; i < iterations; i++) {
    erosion.Run();
  }
  cout << "done" << endl;

  // Save result
  image.Write(output_name);

  return 0;
}
