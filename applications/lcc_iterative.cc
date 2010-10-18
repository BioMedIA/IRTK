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

#include <irtkLargestConnectedComponentIterative.h>

char *input_name = NULL, *output_name = NULL;

void usage()
{
  cerr << "Usage: lcc_iterative [in] [out] <options>\n";
  cerr << "Where <options> are one or more of the following:\n";
  cerr << "\t<-label value>    Value of the target label for which we seek the\n";
  cerr << "\t                  largest connected component.\n";
  cerr << "\t<-3D>             3D largest connected component (Default).\n";
  cerr << "\t<-2D>             2D largest connected component.\n";
  cerr << "\t<-allClusters>    Label all clusters matching target label with a\n";
  cerr << "\t                  sequence starting at 1, 2, ...\n";

  exit(1);
}

int main(int argc, char **argv)
{
  bool ok;
  int label;
  irtkGreyImage image;
  irtkGreyImage output;

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
  output.Read(input_name);

  irtkLargestConnectedComponentIterative<irtkGreyPixel> lcc;

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
    if ((ok == false) && (strcmp(argv[1], "-allClusters") == 0)) {
      argc--;
      argv++;
      lcc.SetAllClustersMode(true);
      ok = true;
    }
    if (ok == false) {
      cerr << "Unknown option: " << argv[1] << endl;
      usage();
    }
  }

  lcc.SetInput (&image);
  lcc.SetOutput(&output);
  lcc.SetTargetLabel(label);

  if (lcc.GetAllClustersMode() == false) {
    cout << "Searching for largest connected component with target label ("<< label <<") ... ";
  } else {
    cout << "Labeling all clusters with the target label (" << label <<") ... ";
  }
  cout.flush();

  lcc.Run();

  cout << "done" << endl;

  // Save result
  output.Write(output_name);

  return 0;
}
