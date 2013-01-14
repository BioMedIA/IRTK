/*=========================================================================
  Date      : $Date: 09.02.2010$
=========================================================================*/

#include <irtkImage.h>
#include <irtkLargeDeformationShooting.h>

void usage(){
  cerr << "Usage: largedeformationsShoot [Template] [Target] <options>\n";
  cerr << "Where <options> are one or more of the following:\n";
  cerr << "  Primary options:\n";
  cerr << "    <-iterations n>             Number of iterations (default=10)\n";
  cerr << "    <-subdivisions n>           Number of subdivisons between t=0 and t=1 (default=10)\n";
  cerr << "    <-MaxVeloUpdt n>            Size of the maximum updates of the vector field (Default=0.5 voxels)\n";
  cerr << "    <-alpha n>                  Weight of the norm in the cost function (Default=0.001)\n";
  cerr << "  Kernels (Default: -sigma 3):\n";
  cerr << "    <-Gauss n>                  Gaussian kernel of std. dev. Sigma (in voxels)\n";
  cerr << "    <-M_Gauss n>                Sum of Gaussian kernels (max 4) -- n = k W1 S1 ... Wk Sk   (k=[#kernels], W.=weight, S.=Sigma)\n" ;
  cerr << "  Inputs and Outputs:\n";
  cerr << "    <-PrefixOutputs n>          Prefix of the output files (default=\"Outputs\")\n";
  cerr << "    <-InputInitMomentum n>      Initial Momentum to initiate the gradient descent (default=\"Null\")\n";
  cerr << "    <-InputIniMoTxt n>          Initial Momentum in an ascii file (instead of nifti) (default=\"Null\")\n";
  cerr << "    <-OutIniMoTxt n>            Outputs the initial momentum in an ascii file (default=\"Null\")\n";
  cerr << "    <-OutVeloField n>           Outputs the 3D+t velocity field (default=\"Null\")\n";
  cerr << "    <-OutDistEnSim n>           Outputs the distance, enrgy and similarity measure (default=\"Null\")\n";
  cerr << "    <-OutDeformation n>         Outputs the deformation (default=\"Null\")\n";
  cerr << "  Secondary options:\n";
  cerr << "    <-indicatorLimiter n>       UpWind -> 0, MinMod -> 1 (default), SuperBee -> 2  \n";
  cerr << "    <-indicatorRungeKutta n>    Euler -> 0 (default), Runge-Kutta -> 1\n";
  cerr << "    <-margins n>                Margin of the image where the calculations are reduced  (default=3 voxels)\n";
  cerr << "    <-GreyLevAlign n>           Grey level linear alignment of each channel -- n = [Padding Src] [Padding Trg]\n";
  exit(1);
}

int main(int argc, char **argv){
  bool ok;
  EulerianShooting Shoot;
  int temp;
  
  // Check command line
  if (argc < 3) 
  {
    usage();
  }

  // read mandatory parameters
  strcpy(Shoot.SourceImageName,argv[1]);
  argc--;  argv++;
  strcpy(Shoot.TargetImageName,argv[1]);
  argc--;  argv++;
  
  // Parse remaining parameters
  while (argc > 1) {
    //1 - Primary options
    ok = false;
    if ((ok == false) && (strcmp(argv[1], "-subdivisions") == 0)) {
      argc--; argv++;
      Shoot.NbTimes = atoi(argv[1]);
      if (Shoot.NbTimes<2) Shoot.NbTimes=2;
      argc--; argv++;
      ok = true;
    }
	if ((ok == false) && (strcmp(argv[1], "-PrefixOutputs") == 0)) {
      argc--; argv++;
      strcpy(Shoot.PrefixOutputs,argv[1]);
      argc--; argv++;
      ok = true;
    }
    if ((ok == false) && (strcmp(argv[1], "-iterations") == 0)) {
      argc--; argv++;
      Shoot.NbIter = atoi(argv[1]);
      if (Shoot.NbIter<0) Shoot.NbIter=0;
      argc--; argv++;
      ok = true;
    }
    if ((ok == false) && (strcmp(argv[1], "-InputInitMomentum") == 0)) {
      argc--; argv++;
      strcpy(Shoot.InputInitialMomentumName,argv[1]);
      Shoot.indicatorInitialMomentum = 1;
      argc--; argv++;
      ok = true;
    }

    if ((ok == false) && (strcmp(argv[1], "-InputIniMoTxt") == 0)) {
      argc--; argv++;
      strcpy(Shoot.InputInitialMomentumName,argv[1]);
      Shoot.indicatorInitialMomentum = 2;
      argc--; argv++;
      ok = true;
    }
    
    if ((ok == false) && (strcmp(argv[1], "-OutIniMoTxt") == 0)) {
      argc--; argv++;
      Shoot.OutIniMoTxt = 1;
      ok = true;
    }
    if ((ok == false) && (strcmp(argv[1], "-OutVeloField") == 0)) {
      argc--; argv++;
      Shoot.OutVeloField = 1;
      ok = true;
    }
    if ((ok == false) && (strcmp(argv[1], "-OutDistEnSim") == 0)) {
      argc--; argv++;
      Shoot.OutDistEnSim = 1;
      ok = true;
    }
    if ((ok == false) && (strcmp(argv[1], "-OutDeformation") == 0)) {
      argc--; argv++;
      Shoot.OutDeformation = 1;
      ok = true;
    }


    
    
    
    if ((ok == false) && (strcmp(argv[1], "-Gauss") == 0)) {
      argc--; argv++;
      Shoot.sigmaX1 = atof(argv[1]);
      Shoot.sigmaY1 = atof(argv[1]);
      Shoot.sigmaZ1 = atof(argv[1]);
      argc--; argv++;
      ok = true;
    }
    if ((ok == false) && (strcmp(argv[1], "-M_Gauss") == 0)) {
      argc--; argv++;
      temp= atoi(argv[1]);
      argc--; argv++;
      if (temp>=1){
        Shoot.weight1 = atof(argv[1]);
        argc--; argv++;
        Shoot.sigmaX1 = atof(argv[1]); Shoot.sigmaY1 = atof(argv[1]); Shoot.sigmaZ1 = atof(argv[1]);
        argc--; argv++;
      }
      if (temp>=2){
        Shoot.weight2 = atof(argv[1]);
        argc--; argv++;
        Shoot.sigmaX2 = atof(argv[1]); Shoot.sigmaY2 = atof(argv[1]); Shoot.sigmaZ2 = atof(argv[1]);
        argc--; argv++;
      }
      if (temp>=3){
        Shoot.weight3 = atof(argv[1]);
        argc--; argv++;
        Shoot.sigmaX3 = atof(argv[1]); Shoot.sigmaY3 = atof(argv[1]); Shoot.sigmaZ3 = atof(argv[1]);
        argc--; argv++;
      }
      if (temp>=4){
        Shoot.weight4 = atof(argv[1]);
        argc--; argv++;
        Shoot.sigmaX4 = atof(argv[1]); Shoot.sigmaY4 = atof(argv[1]); Shoot.sigmaZ4 = atof(argv[1]);
        argc--; argv++;
      }
      ok = true;
    }
    if ((ok == false) && (strcmp(argv[1], "-alpha") == 0)) {
      argc--; argv++;
      Shoot.alpha = atof(argv[1]);
      argc--; argv++;
      ok = true;
    }
    if ((ok == false) && (strcmp(argv[1], "-MaxVeloUpdt") == 0)) {
      argc--; argv++;
	  Shoot.MaxUpdate = atof(argv[1]);
      argc--; argv++;
      ok = true;
    }
	if ((ok == false) && (strcmp(argv[1], "-indicatorLimiter") == 0)) {
      argc--; argv++;
	  Shoot.indicatorLimiter = atoi(argv[1]);
      argc--; argv++;
      ok = true;
    }
	if ((ok == false) && (strcmp(argv[1], "-indicatorRungeKutta") == 0)) {
      argc--; argv++;
	  Shoot.indicatorRungeKutta = atoi(argv[1]);
      argc--; argv++;
      ok = true;
    }
	if ((ok == false) && (strcmp(argv[1], "-margins") == 0)) {
      argc--; argv++;
      Shoot.Margin = atoi(argv[1]);
      argc--; argv++;
      ok = true;
    }
	if ((ok == false) && (strcmp(argv[1], "-GreyLevAlign") == 0)) {
      argc--; argv++;
      Shoot.GreyLevAlign = 1;
      Shoot.GLA_Padding_Src = atof(argv[1]);
      argc--; argv++;
      Shoot.GLA_Padding_Trg = atof(argv[1]);
      argc--; argv++;
      ok = true;
    }
	if (ok == false) 
	{
		usage();
	}
  }
  
  //check the inputs
  printf("\nPARAMETERS:\n");
  printf("  Source image:           %s\n",Shoot.SourceImageName);
  printf("  Target image:           %s\n",Shoot.TargetImageName);
  printf("  Subdivisions number:    %d\n",Shoot.NbTimes);
  printf("  Gaussian Std. dev.:     %lf\n\n",Shoot.sigmaX1);
  
  
  //run process
  cout << "Large deformation registration using geodesic shooting... \n"; cout.flush();
  Shoot.Run();
  cout << "done" << endl;

  return 0;
}


