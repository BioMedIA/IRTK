/*=========================================================================

  Date      : $Date: 29.06.2009$
  Changes   : $Authors: Laurent Risser, Francois-Xavier Vialard$

=========================================================================*/

#include <irtkImage.h>
#include <irtkLargeDeformationGradientLagrange.h>

void usage(){
  cerr << "Usage: largedeformations [Template] [Target] [Deformed Template]<options>\n";
  cerr << "Where <options> are one or more of the following:\n";
  cerr << "  Primary options:\n";
  cerr << "    <-iterations n>     Number of iterations (default=10)\n";
  cerr << "    <-subdivisions n>   Number of subdivisons (default=10)\n";
  cerr << "    <-sigma n>          Std. dev. of the Gaussian kernel (default=3.0 voxels)\n";
  cerr << "    <-MaxVeloUpdt n>    Maximum velocity update at each iteration (default=3.0 voxels)\n";
  cerr << "  Inputs and Outputs:\n";
  cerr << "    <-PrefixInputVF n>  Prefix of the files containing an initial velocity field (default=\"Null\")\n";
  cerr << "    <-PrefixOutputVF n> Prefix of the files where the final velocity field is saved (default=\"Null\")\n";
  cerr << "  Secondary options:\n";
  cerr << "    <-margins n>        Margin of the image where the calculations are reduced  (default=0 voxels)\n";
  cerr << "    <-reparamFreq n>    Frequency of the velocity field reparameterizations  (default=10 iterations)\n";
  cerr << "    <-DeltaVox n>       Size of a voxel in each direction (default=1)\n";
  cerr << "    <-alpha n>          Weight of the velocity field smoothness for the energy  (default=0.01)\n";
  cerr << "    <-gamma n>          Weight of the velocity field deformations for the energy (default=0.1)\n";
  cerr << "    <-epsilon n>        Threshold on the energy gradient convergence  (default=0.1)\n";

  exit(1);
}

int main(int argc, char **argv){
  char *template_name = NULL, *target_name = NULL, *deformed_template_name = NULL;
  Bool ok;
  irtkGreyImage template_image;
  irtkGreyImage target_image;

  // Check command line
  if (argc < 4) {
    usage();
  }

  // Read input and output names
  template_name  = argv[1];
  argc--;
  argv++;
  target_name = argv[1];
  argc--;
  argv++;
  deformed_template_name = argv[1];
  argc--;
  argv++;
  
  // Read inputs
  template_image.Read(template_name);

  //allocate LargeDefGradLagrange class
  LargeDefGradLagrange<irtkGreyPixel> LargeDefGradLagrange;
  LargeDefGradLagrange.SetInput(&template_image);
  LargeDefGradLagrange.SetOutput(&template_image);
  LargeDefGradLagrange.target_image.Read(target_name);

  
  // Parse remaining parameters
  while (argc > 1) {
    ok = False;
    if ((ok == False) && (strcmp(argv[1], "-iterations") == 0)) {
      argc--;
      argv++;
      LargeDefGradLagrange.iteration_nb = atoi(argv[1]);
      argc--;
      argv++;
      ok = True;
    }
    if ((ok == False) && (strcmp(argv[1], "-subdivisions") == 0)) {
      argc--;
      argv++;
      LargeDefGradLagrange.NbTimeSubdiv = atoi(argv[1]);
      argc--;
      argv++;
      ok = True;
    }
    if ((ok == False) && (strcmp(argv[1], "-margins") == 0)) {
      argc--;
      argv++;
      LargeDefGradLagrange.Margin = atoi(argv[1]);
      argc--;
      argv++;
      ok = True;
    }
    if ((ok == False) && (strcmp(argv[1], "-sigma") == 0)) {
      argc--;
      argv++;
      LargeDefGradLagrange.sigma = atof(argv[1]);
      argc--;
      argv++;
      ok = True;
    }
    if ((ok == False) && (strcmp(argv[1], "-DeltaVox") == 0)) {
      argc--;
      argv++;
      LargeDefGradLagrange.DeltaVox = atof(argv[1]);
      argc--;
      argv++;
      ok = True;
    }
    if ((ok == False) && (strcmp(argv[1], "-MaxVeloUpdt") == 0)) {
      argc--;
      argv++;
      LargeDefGradLagrange.MaxVelocityUpdate = atof(argv[1]);
      argc--;
      argv++;
      ok = True;
    }
    if ((ok == False) && (strcmp(argv[1], "-epsilon") == 0)) {
      argc--;
      argv++;
      LargeDefGradLagrange.epsilon = atof(argv[1]);
      argc--;
      argv++;
      ok = True;
    }
    if ((ok == False) && (strcmp(argv[1], "-reparamFreq") == 0)) {
      argc--;
      argv++;
      LargeDefGradLagrange.reparametrization = atoi(argv[1]);
      argc--;
      argv++;
      ok = True;
    }
    if ((ok == False) && (strcmp(argv[1], "-alpha") == 0)) {
      argc--;
      argv++;
      LargeDefGradLagrange.alpha = atof(argv[1]);
      argc--;
      argv++;
      ok = True;
    }
    if ((ok == False) && (strcmp(argv[1], "-gamma") == 0)) {
      argc--;
      argv++;
      LargeDefGradLagrange.gamma = atof(argv[1]);
      argc--;
      argv++;
      ok = True;
    }
    if ((ok == False) && (strcmp(argv[1], "-PrefixInputVF") == 0)) {
      argc--;
      argv++;
      strcpy(LargeDefGradLagrange.PrefixInputVF,argv[1]);
      argc--;
      argv++;
      ok = True;
    }
    if ((ok == False) && (strcmp(argv[1], "-PrefixOutputVF") == 0)) {
      argc--;
      argv++;
      strcpy(LargeDefGradLagrange.PrefixOutputVF,argv[1]);
      argc--;
      argv++;
      ok = True;
    }
    if (ok == False) usage();
  }

  //run process
  cout << "Large Deformation registration using Beg 05's technique ... \n"; cout.flush();
  LargeDefGradLagrange.Run();
  cout << "done" << endl;

  // Save result
  template_image.Write(deformed_template_name);

  return 0;
}
