/*=========================================================================

  Date      : $Date: 29.06.2009$
  Changes   : $Authors: Laurent Risser, Francois-Xavier Vialard$

=========================================================================*/

#include <irtkImage.h>
#include <irtkLargeDeformationGradientLagrange.h>

void usage()
{
  cerr << "Usage: largedeformations [Template] [Target] [Deformed Template]<options>\n";
  cerr << "Where <options> are one or more of the following:\n";
  cerr << "\t<-iterations n>     Number of iterations (default=3)\n";
  cerr << "\t<-subdivisions n>   Number of subdivisons (default=10)\n";
  cerr << "\t<-sigma n>          Std. dev. of the Gaussian kernel (default=4)\n";
  cerr << "\t<-DeltaVox n>       Step between two voxels (default=3)\n";
  cerr << "\t<-MaxVeloUpdt n>    Maximum velocity update in voxels at each iteration (default=10)\n";
  cerr << "\t<-epsilon n>        Threshold on the energy gradient convergence  (default=0.1)\n";
  cerr << "\t<-WgtVelocity n>    Weight of the velocity field to evaluate energy gradients (default=0.01)\n";
  cerr << "\t<-PrefixInputVF n>  Prefix of the files containing an initial velocity field (default=\"Null\")\n";
  cerr << "\t<-PrefixOutputVF n> Prefix of the files where the final velocity field is saved (default=\"Null\")\n";
  exit(1);
}

int main(int argc, char **argv)
{
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
    if ((ok == False) && (strcmp(argv[1], "-WgtVelocity") == 0)) {
      argc--;
      argv++;
      LargeDefGradLagrange.WeightVelocityField = atof(argv[1]);
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
