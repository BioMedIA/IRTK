/*

  Library   : Image Registration Toolkit (IRTK)
  Module    : $Id$
  pa100 $
  Copyright : Imperial College, Department of Computing
              Visual Information Processing (VIP), 2008 onwards
  Date      : $Date$
  $
  Version   : $Revision$
  Changes   : $Author$

*/


#include <irtkImage.h>

#include <irtkModeFilter.h>

char *input_name = NULL, *output_name = NULL;

void usage()
{
  cerr << "Usage: modefilter [in] [out] <options>\n";
  cerr << "" << endl;
  cerr << "Apply a mode filter over local window for each voxel in a label image." << endl;
  cerr << "The window consists of the immediate neighbours of each voxel according to the " << endl;
  cerr << "type of connecivity specified." << endl;
  cerr << "<options> are one or more of the following:\n";
  cerr << "\t<-iterations n>    Number of iterations\n";
  cerr << "\t<-connectivity n>  Type of voxel neighbourhood connectivity. "<< endl;
  cerr << "\t                   Valid choices are 6, 18 or 26 (default)\n";
  exit(1);
}

int main(int argc, char **argv)
{
  bool ok;
  int i, iterations, connectivity;
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

  connectivity = 26;

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
    if ((ok == false) && (strcmp(argv[1], "-connectivity") == 0)) {
      argc--;
      argv++;
      connectivity = atoi(argv[1]);
      argc--;
      argv++;
      ok = true;
    }
    if (ok == false) {
      cerr << "Unknown option: " << argv[1] << endl;
      usage();
    }
  }

  cout << "Running filter ... "; cout.flush();
  irtkModeFilter<irtkGreyPixel> modefilter;

  cout << "Setting connectivity to " << connectivity << endl;
  switch (connectivity){
  case 6:
    modefilter.SetConnectivity(CONNECTIVITY_06);
    break;
  case 18:
    modefilter.SetConnectivity(CONNECTIVITY_18);
    break;
  case 26:
    modefilter.SetConnectivity(CONNECTIVITY_26);
    break;
  default:
    cerr << "Invalid connectivity ("<< connectivity <<"). Exiting." << endl;
    usage();
    break;
  }

  modefilter.SetInput(&image);
  modefilter.SetOutput(&image);
  for (i = 0; i < iterations; i++){
    modefilter.Run();
  }

  cout << "done" << endl;

  // Save result
  image.Write(output_name);

  return 0;
}


