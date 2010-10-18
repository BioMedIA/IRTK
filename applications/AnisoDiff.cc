/*=========================================================================

  Date      : $Date: 03.03.2009$
  Changes   : $Author: Laurent Risser $

=========================================================================*/

#include <irtkImage.h>
#include <irtkAnisoDiffusion.h>

char *input_name = NULL, *output_name = NULL;

void usage()
{
  cerr << "Usage: AnisoDiff [in] [out] <options>\n";
  cerr << "Where <options> are one or more of the following:\n";
  cerr << "\t<-iterations n>     Number of iterations (default=5)\n";
  cerr << "\t<-SemiImplicit n>   1: semi implicit ADI scheme / 0: explicit scheme  (default=1)\n";
  cerr << "\t<-TimeDependent n>  1: 4D scheme / 0: 3D scheme (default=0)\n";
  cerr << "\t<-ax n>             Threshold on gray level variations along x (default=3.)\n";
  cerr << "\t<-ay n>             Threshold on gray level variations along y (default=3.)\n";
  cerr << "\t<-az n>             Threshold on gray level variations along z (default=3.)\n";
  cerr << "\t<-at n>             Threshold on gray level variations along t (default=3.)\n";
  cerr << "\t<-dTau n>           Virtual time step of the PDE (default=1.)\n";
  exit(1);
}

int main(int argc, char **argv)
{
  bool ok;
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

  //allocate anisoDiffusion class
  anisoDiffusion<irtkGreyPixel> anisodiff;
  anisodiff.SetInput(&image);
  anisodiff.SetOutput(&image);

  // Parse remaining parameters
  while (argc > 1) {
    ok = false;
    if ((ok == false) && (strcmp(argv[1], "-iterations") == 0)) {
      argc--;
      argv++;
      anisodiff.ITERATIONS_NB = atoi(argv[1]);
      argc--;
      argv++;
      ok = true;
    }
    if ((ok == false) && (strcmp(argv[1], "-SemiImplicit") == 0)) {
      argc--;
      argv++;
      anisodiff.SemiImplicit = atoi(argv[1]);
      argc--;
      argv++;
      ok = true;
    }
    if ((ok == false) && (strcmp(argv[1], "-TimeDependent") == 0)) {
      argc--;
      argv++;
      anisodiff.TimeDependent = atoi(argv[1]);
      argc--;
      argv++;
      ok = true;
    }
    if ((ok == false) && (strcmp(argv[1], "-ax") == 0)) {
      argc--;
      argv++;
      anisodiff.ax = atof(argv[1]);
      argc--;
      argv++;
      ok = true;
    }
    if ((ok == false) && (strcmp(argv[1], "-ay") == 0)) {
      argc--;
      argv++;
      anisodiff.ay = atof(argv[1]);
      argc--;
      argv++;
      ok = true;
    }
    if ((ok == false) && (strcmp(argv[1], "-az") == 0)) {
      argc--;
      argv++;
      anisodiff.az = atof(argv[1]);
      argc--;
      argv++;
      ok = true;
    }
    if ((ok == false) && (strcmp(argv[1], "-at") == 0)) {
      argc--;
      argv++;
      anisodiff.at = atof(argv[1]);
      argc--;
      argv++;
      ok = true;
    }
    if ((ok == false) && (strcmp(argv[1], "-dTau") == 0)) {
      argc--;
      argv++;
      anisodiff.dTau = atof(argv[1]);
      argc--;
      argv++;
      ok = true;
    }
    if (ok == false) {
      cerr << "Unknown option: " << argv[1] << endl;
      usage();
    }
  }

  //run process
  cout << "Anisotropic diffusion ... \n"; cout.flush();
  anisodiff.Run();
  cout << "done" << endl;

  // Save result
  image.Write(output_name);

  return 0;
}
