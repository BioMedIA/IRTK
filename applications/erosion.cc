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

#include <irtkErosion.h>

char *input_name = NULL, *output_name = NULL;

void usage()
{
  cerr << "Usage: erosion [in] [out] <options>\n";
  cerr << "Where <options> are one or more of the following:\n";
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

  cout << "Eroding ... "; cout.flush();
  irtkErosion<irtkGreyPixel> erosion;

  cout << "Setting connectivity to " << connectivity << endl;
  switch (connectivity){
  case 6:
    erosion.SetConnectivity(CONNECTIVITY_06);
    break;
  case 18:
  	erosion.SetConnectivity(CONNECTIVITY_18);
  	break;
  case 26:
  	erosion.SetConnectivity(CONNECTIVITY_26);
  	break;
  default:
  	cerr << "Invalid connectivity ("<< connectivity <<"). Exiting." << endl;
  	usage();
  }

  erosion.SetInput(&image);
  erosion.SetOutput(&image);
  for (i = 0; i < iterations; i++) erosion.Run();
  cout << "done" << endl;

  // Save result
  image.Write(output_name);

  return 0;
}
