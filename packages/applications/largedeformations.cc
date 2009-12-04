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
  cerr << "    <-MaxVeloUpdt n>    Maximum velocity update at each iteration (default=2.0 voxels)\n";
  cerr << "  Inputs and Outputs:\n";
  cerr << "    <-PrefixInputVF n>  Prefix of the files containing an initial velocity field (default=\"Null\")\n";
  cerr << "    <-PrefixOutputVF n> Prefix of the files where the final velocity field is saved (default=\"Null\")\n";
  cerr << "  Kernels ((alpha and gamma) by default or another option):\n";
  cerr << "    <-alpha n>          Velocity field smoothness weight in the kernel(default=0.01)\n";
  cerr << "    <-gamma n>          Velocity field deformations weight in the kernel (default=0.1)\n";
  cerr << "    <-sigma n>          Gaussian kernel of std. dev. Sigma\n";
  cerr << "    <-AnisoGauss n>     Anisotropic Gaussian kernel of weight 'Wgt' and std. dev. 'SigmaX SigmaY SigmaZ'\n";
  cerr << "    <-MultiGauss2 n>    Kernel = Sum of 2 weighted anisotropic Gaussians - W1 SX1 SY1 SZ1  W2 SX2 SY2 SZ2\n";
  cerr << "    <-MultiGauss3 n>    Kernel = Sum of 3 weighted anisotropic Gaussians - W1 ... SZ3\n";
  cerr << "    <-MultiGauss4 n>    Kernel = Sum of 4 weighted anisotropic Gaussians - W1 ... SZ4\n";
  cerr << "    <-ChainGauss2 n>    Chain of 2 deformations using aniso. Gaussians kernels - W1 ... SZ2\n";
  cerr << "    <-ChainGauss3 n>    Chain of 3 deformations using aniso. Gaussians kernels - W1 ... SZ3\n";
  cerr << "    <-ChainGauss4 n>    Chain of 4 deformations using aniso. Gaussians kernels - W1 ... SZ4\n";
  cerr << "  Secondary options:\n";
  cerr << "    <-margins n>        Margin of the image where the calculations are reduced  (default=0 voxels)\n";
  cerr << "    <-WghtVeloField n>  Weight of the velocity field in the energy (default=1.) \n";
  cerr << "    <-reparamFreq n>    Frequency of the velocity field reparameterizations  (default=200 iterations)\n";
  cerr << "    <-epsilon n>        Threshold on the energy gradient convergence  (default=0.1)\n";
  
  exit(1);
}

int main(int argc, char **argv){
  char *template_name = NULL, *target_name = NULL, *deformed_template_name = NULL;
  Bool ok;
  irtkGreyImage template_image;
  irtkGreyImage target_image;
  LargeDefGradLagrange<irtkGreyPixel> LargeDefGradLagrange;

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
  LargeDefGradLagrange.SetInput(&template_image);
  LargeDefGradLagrange.SetOutput(&template_image);
  LargeDefGradLagrange.target_image.Read(target_name);

  
  // Parse remaining parameters
  while (argc > 1) {
    //1 - Primary options
    ok = False;
    if ((ok == False) && (strcmp(argv[1], "-iterations") == 0)) {
      argc--; argv++;
      LargeDefGradLagrange.iteration_nb = atoi(argv[1]);
      argc--; argv++;
      ok = True;
    }
    if ((ok == False) && (strcmp(argv[1], "-subdivisions") == 0)) {
      argc--; argv++;
      LargeDefGradLagrange.NbTimeSubdiv = atoi(argv[1]);
      argc--; argv++;
      ok = True;
    }
    if ((ok == False) && (strcmp(argv[1], "-MaxVeloUpdt") == 0)) {
      argc--; argv++;
      LargeDefGradLagrange.MaxVelocityUpdate = atof(argv[1]);
      argc--; argv++;
      ok = True;
    }
    //2 - Inputs and Outputs
    if ((ok == False) && (strcmp(argv[1], "-PrefixInputVF") == 0)) {
      argc--; argv++;
      strcpy(LargeDefGradLagrange.PrefixInputVF,argv[1]);
      argc--; argv++;
      ok = True;
    }
    if ((ok == False) && (strcmp(argv[1], "-PrefixOutputVF") == 0)) {
      argc--; argv++;
      strcpy(LargeDefGradLagrange.PrefixOutputVF,argv[1]);
      argc--; argv++;
      ok = True;
    }
    //3 - Kernels
    if ((ok == False) && (strcmp(argv[1], "-alpha") == 0)) {
      argc--; argv++;
      LargeDefGradLagrange.alpha = atof(argv[1]);
      argc--; argv++;
      ok = True;
    }
    if ((ok == False) && (strcmp(argv[1], "-gamma") == 0)) {
      argc--; argv++;
      LargeDefGradLagrange.gamma = atof(argv[1]);
      argc--; argv++;
      ok = True;
    }
    if ((ok == False) && (strcmp(argv[1], "-sigma") == 0)) {
      argc--; argv++;
      LargeDefGradLagrange.sigmaX1 = atof(argv[1]);
      LargeDefGradLagrange.sigmaY1 = atof(argv[1]);
      LargeDefGradLagrange.sigmaZ1 = atof(argv[1]);
      argc--; argv++;
      ok = True;
    }
    if ((ok == False) && (strcmp(argv[1], "-AnisoGauss") == 0)) {
      argc--; argv++;
      LargeDefGradLagrange.weight1 = atof(argv[1]);
      argc--; argv++;
      LargeDefGradLagrange.sigmaX1 = atof(argv[1]);
      argc--; argv++;
      LargeDefGradLagrange.sigmaY1 = atof(argv[1]);
      argc--; argv++;
      LargeDefGradLagrange.sigmaZ1 = atof(argv[1]);
      argc--; argv++;
      ok = True;
    }
    if ((ok == False) && (strcmp(argv[1], "-MultiGauss2") == 0)) {
      argc--; argv++;
      LargeDefGradLagrange.weight1 = atof(argv[1]);
      argc--; argv++;
      LargeDefGradLagrange.sigmaX1 = atof(argv[1]);
      argc--; argv++;
      LargeDefGradLagrange.sigmaY1 = atof(argv[1]);
      argc--; argv++;
      LargeDefGradLagrange.sigmaZ1 = atof(argv[1]);
      argc--; argv++;
      LargeDefGradLagrange.weight2 = atof(argv[1]);
      argc--; argv++;
      LargeDefGradLagrange.sigmaX2 = atof(argv[1]);
      argc--; argv++;
      LargeDefGradLagrange.sigmaY2 = atof(argv[1]);
      argc--; argv++;
      LargeDefGradLagrange.sigmaZ2 = atof(argv[1]);
      argc--; argv++;
      ok = True;
    }
    if ((ok == False) && (strcmp(argv[1], "-MultiGauss3") == 0)) {
      argc--; argv++;
      LargeDefGradLagrange.weight1 = atof(argv[1]);
      argc--; argv++;
      LargeDefGradLagrange.sigmaX1 = atof(argv[1]);
      argc--; argv++;
      LargeDefGradLagrange.sigmaY1 = atof(argv[1]);
      argc--; argv++;
      LargeDefGradLagrange.sigmaZ1 = atof(argv[1]);
      argc--; argv++;
      LargeDefGradLagrange.weight2 = atof(argv[1]);
      argc--; argv++;
      LargeDefGradLagrange.sigmaX2 = atof(argv[1]);
      argc--; argv++;
      LargeDefGradLagrange.sigmaY2 = atof(argv[1]);
      argc--; argv++;
      LargeDefGradLagrange.sigmaZ2 = atof(argv[1]);
      argc--; argv++;
      LargeDefGradLagrange.weight3 = atof(argv[1]);
      argc--; argv++;
      LargeDefGradLagrange.sigmaX3 = atof(argv[1]);
      argc--; argv++;
      LargeDefGradLagrange.sigmaY3 = atof(argv[1]);
      argc--; argv++;
      LargeDefGradLagrange.sigmaZ3 = atof(argv[1]);
      argc--; argv++;
      ok = True;
    }
    if ((ok == False) && (strcmp(argv[1], "-MultiGauss4") == 0)) {
      argc--; argv++;
      LargeDefGradLagrange.weight1 = atof(argv[1]);
      argc--; argv++;
      LargeDefGradLagrange.sigmaX1 = atof(argv[1]);
      argc--; argv++;
      LargeDefGradLagrange.sigmaY1 = atof(argv[1]);
      argc--; argv++;
      LargeDefGradLagrange.sigmaZ1 = atof(argv[1]);
      argc--; argv++;
      LargeDefGradLagrange.weight2 = atof(argv[1]);
      argc--; argv++;
      LargeDefGradLagrange.sigmaX2 = atof(argv[1]);
      argc--; argv++;
      LargeDefGradLagrange.sigmaY2 = atof(argv[1]);
      argc--; argv++;
      LargeDefGradLagrange.sigmaZ2 = atof(argv[1]);
      argc--; argv++;
      LargeDefGradLagrange.weight3 = atof(argv[1]);
      argc--; argv++;
      LargeDefGradLagrange.sigmaX3 = atof(argv[1]);
      argc--; argv++;
      LargeDefGradLagrange.sigmaY3 = atof(argv[1]);
      argc--; argv++;
      LargeDefGradLagrange.sigmaZ3 = atof(argv[1]);
      argc--; argv++;
      LargeDefGradLagrange.weight4 = atof(argv[1]);
      argc--; argv++;
      LargeDefGradLagrange.sigmaX4 = atof(argv[1]);
      argc--; argv++;
      LargeDefGradLagrange.sigmaY4 = atof(argv[1]);
      argc--; argv++;
      LargeDefGradLagrange.sigmaZ4 = atof(argv[1]);
      argc--; argv++;
      ok = True;
    }
    if ((ok == False) && (strcmp(argv[1], "-ChainGauss2") == 0)) {
      argc--; argv++;
      LargeDefGradLagrange.weight1 = atof(argv[1]);
      argc--; argv++;
      LargeDefGradLagrange.sigmaX1 = atof(argv[1]);
      argc--; argv++;
      LargeDefGradLagrange.sigmaY1 = atof(argv[1]);
      argc--; argv++;
      LargeDefGradLagrange.sigmaZ1 = atof(argv[1]);
      argc--; argv++;
      LargeDefGradLagrange.weight2 = atof(argv[1]);
      argc--; argv++;
      LargeDefGradLagrange.sigmaX2 = atof(argv[1]);
      argc--; argv++;
      LargeDefGradLagrange.sigmaY2 = atof(argv[1]);
      argc--; argv++;
      LargeDefGradLagrange.sigmaZ2 = atof(argv[1]);
      argc--; argv++;
      LargeDefGradLagrange.NbKernels=2;
      ok = True;
    }
    if ((ok == False) && (strcmp(argv[1], "-ChainGauss3") == 0)) {
      argc--; argv++;
      LargeDefGradLagrange.weight1 = atof(argv[1]);
      argc--; argv++;
      LargeDefGradLagrange.sigmaX1 = atof(argv[1]);
      argc--; argv++;
      LargeDefGradLagrange.sigmaY1 = atof(argv[1]);
      argc--; argv++;
      LargeDefGradLagrange.sigmaZ1 = atof(argv[1]);
      argc--; argv++;
      LargeDefGradLagrange.weight2 = atof(argv[1]);
      argc--; argv++;
      LargeDefGradLagrange.sigmaX2 = atof(argv[1]);
      argc--; argv++;
      LargeDefGradLagrange.sigmaY2 = atof(argv[1]);
      argc--; argv++;
      LargeDefGradLagrange.sigmaZ2 = atof(argv[1]);
      argc--; argv++;
      LargeDefGradLagrange.weight3 = atof(argv[1]);
      argc--; argv++;
      LargeDefGradLagrange.sigmaX3 = atof(argv[1]);
      argc--; argv++;
      LargeDefGradLagrange.sigmaY3 = atof(argv[1]);
      argc--; argv++;
      LargeDefGradLagrange.sigmaZ3 = atof(argv[1]);
      argc--; argv++;
      LargeDefGradLagrange.NbKernels=3;
      ok = True;
    }
    if ((ok == False) && (strcmp(argv[1], "-ChainGauss4") == 0)) {
      argc--; argv++;
      LargeDefGradLagrange.weight1 = atof(argv[1]);
      argc--; argv++;
      LargeDefGradLagrange.sigmaX1 = atof(argv[1]);
      argc--; argv++;
      LargeDefGradLagrange.sigmaY1 = atof(argv[1]);
      argc--; argv++;
      LargeDefGradLagrange.sigmaZ1 = atof(argv[1]);
      argc--; argv++;
      LargeDefGradLagrange.weight2 = atof(argv[1]);
      argc--; argv++;
      LargeDefGradLagrange.sigmaX2 = atof(argv[1]);
      argc--; argv++;
      LargeDefGradLagrange.sigmaY2 = atof(argv[1]);
      argc--; argv++;
      LargeDefGradLagrange.sigmaZ2 = atof(argv[1]);
      argc--; argv++;
      LargeDefGradLagrange.weight3 = atof(argv[1]);
      argc--; argv++;
      LargeDefGradLagrange.sigmaX3 = atof(argv[1]);
      argc--; argv++;
      LargeDefGradLagrange.sigmaY3 = atof(argv[1]);
      argc--; argv++;
      LargeDefGradLagrange.sigmaZ3 = atof(argv[1]);
      argc--; argv++;
      LargeDefGradLagrange.weight4 = atof(argv[1]);
      argc--; argv++;
      LargeDefGradLagrange.sigmaX4 = atof(argv[1]);
      argc--; argv++;
      LargeDefGradLagrange.sigmaY4 = atof(argv[1]);
      argc--; argv++;
      LargeDefGradLagrange.sigmaZ4 = atof(argv[1]);
      argc--; argv++;
      LargeDefGradLagrange.NbKernels=4;
      ok = True;
    }

    //4 - Secondary options
    if ((ok == False) && (strcmp(argv[1], "-margins") == 0)) {
      argc--; argv++;
      LargeDefGradLagrange.Margin = atoi(argv[1]);
      argc--; argv++;
      ok = True;
    }
    if ((ok == False) && (strcmp(argv[1], "-WghtVeloField") == 0)) {
      argc--; argv++;
      LargeDefGradLagrange.WghtVelField = atof(argv[1]);
      argc--; argv++;
      ok = True;
    }
    if ((ok == False) && (strcmp(argv[1], "-reparamFreq") == 0)) {
      argc--; argv++;
      LargeDefGradLagrange.reparametrization = atoi(argv[1]);
      argc--; argv++;
      ok = True;
    }
    if ((ok == False) && (strcmp(argv[1], "-epsilon") == 0)) {
      argc--; argv++;
      LargeDefGradLagrange.epsilon = atof(argv[1]);
      argc--; argv++;
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
