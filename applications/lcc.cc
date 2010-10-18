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

#include <irtkLargestConnectedComponent.h>

char *input_name = NULL, *output_name = NULL;

void usage()
{
  cerr << "Usage: lcc [in] [out] <options>\n";
  cerr << "Where <options> are one or more of the following:\n";
  cerr << "\t<-2D>             2D largest connected component\n";
  cerr << "\t<-3D>             3D largest connected component (default)\n";
  cerr << "\t<-label value>    Value of label for which to search for largest connected component\n";
  exit(1);
}

int main(int argc, char **argv)
{
  bool ok;
  int label;
  irtkByteImage image;

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

  irtkLargestConnectedComponent<irtkBytePixel> lcc;

  // Parse remaining parameters
  label = 1;
  while (argc > 1) {
    ok = false;
    if ((ok == false) && (strcmp(argv[1], "-label") == 0)) {
      argc--;
      argv++;
      label = atoi(argv[1]);
      argc--;
      argv++;
      ok = true;
    }
    if ((ok == false) && (strcmp(argv[1], "-2D") == 0)) {
      argc--;
      argv++;
      lcc.SetMode2D(true);
      ok = true;
    }
    if ((ok == false) && (strcmp(argv[1], "-3D") == 0)) {
      argc--;
      argv++;
      lcc.SetMode2D(false);
      ok = true;
    }
    if (ok == false) {
      cerr << "Unknown option: " << argv[1] << endl;
      usage();
    }
  }

  cout << "Searching for largest connected component ... "; cout.flush();
  lcc.SetInput (&image);
  lcc.SetOutput(&image);
  lcc.SetClusterLabel(label);
  lcc.Run();
  cout << "done" << endl;

  // Save result
  image.Write(output_name);

  return 0;
}
