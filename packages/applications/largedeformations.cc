/*=========================================================================
 
 Date      : $Date: 29.04.2010$
 Changes   : $Authors: Laurent Risser, Francois-Xavier Vialard$
 
 =========================================================================*/

#include <irtkLargeDeformationGradientLagrange.h>

void usage(){
	cerr << "Usage: largedeformations [Source] [Target] <options>\n";
	cerr << "Where <options> are one or more of the following:\n";
	cerr << "  Primary options:\n";
	cerr << "    <-iterations n>      Number of iterations (default=10)\n";
	cerr << "    <-subdivisions n>    Number of subdivisons (default=10)\n";
	cerr << "    <-MaxVeloUpdt n>     Maximum velocity update at each iteration (default=0.4 voxels)\n";
	cerr << "  Inputs and Outputs:\n";
	cerr << "    <-PrefixInputs n>    Prefix of the files containing an initial velocity field (default=\"Null\")\n";
	cerr << "    <-PrefixOutputs n>   Prefix of the files containing the outputs (default=\"Outputs\")\n";
	cerr << "    <-AddChannel W S T>  Add a channel -> W=weight (wgt of ref channel is 1) S=Source T=Target\n";
	cerr << "    <-Mask n>            Definition of a mask (default=\"Null\")\n";
	cerr << "  Kernels (Default: -Gauss 1):\n";
	cerr << "    <-Gauss n>           Gaussian kernel -- n = S  (the std. dev.)\n";
	cerr << "    <-M_Gauss n>         Sum of Gaussian kernels (max 7) -- n = k W1 S1 ... Wk Sk   (k=[#kernels], W.=weight)\n" ;
	cerr << "    <-A_Gauss n>         Anisotropic Gaussian kernel -- n = SX SY SZ    (the std. dev.)\n";
	cerr << "    <-M_A_Gauss n>       Sum of anisotropic Gaussian kernels (max 7) -- n = k W1 SX1 SY1 SZ1 ... Wk SXk SYk SZk\n" ;
	cerr << "  Secondary options:\n";
	cerr << "    <-symmetric>         Symmetric registration\n";
	cerr << "    <-epsilon n>         Threshold on the normalized max update of the velicty field (default=0.2)\n";
	cerr << "    <-GreyLevAlign n>    Grey level linear alignment of each channel (Inputs: Padding Src - Padding Trg)\n";
	cerr << "    <-margins n>         Margin of the image where the calculations are reduced  (default=0 voxels)\n";
	cerr << "    <-WghtVeloField n>   Weight of the velocity field in the energy (default=1.) \n";
	cerr << "    <-RefMaxGrad n>      Value to manage the convergence. Automatically configured if <0 (default=-1.)\n";
	cerr << "    <-MapSrcImag n>      Compose the deformation with a f. mapping of the Src image (default=\"Null\")\n";
	cerr << "    <-MapTrgImag n>      Compose the deformation with a b. mapping of the Trg image (default=\"Null\")\n";
	cerr << "  Special Outputs:\n";
	cerr << "    <-SplitKernels>      Split the contribution of each kernel\n";
	cerr << "    <-AOD>               Amplitude of the deformations from each voxel of the source image\n";
	cerr << "    <-DetJacobian>       Determinant of the Jacobian at each voxel\n";
	cerr << "    <-FinalDefVec>       Vector field of the estimated deformation from [Source] to [Target]\n";
	cerr << "    <-FinalDefInvVec>    Vector field of the estimated deformation from [Target] to [Source]\n";
	cerr << "    <-InitMomentum>      Estimated initial momentum\n";
	cerr << "    <-ShowSSD>           Show the Sum of the Squared Differences at t=1 ieration after iteration\n";
	cerr << "  \n";
	cerr << "  Measure the inverse typical amplitude of the updates between [Source] and [Target] \n";
	cerr << "  with a Gaussian kernel of stddev SigmaX SigmaY SigmaZ: <-TypicInvAmp SigmaX SigmaY SigmaZ>. \n";
	cerr << "  This measure is 1/g(...) in the MICCAI paper and 1/f(...) in the WBIR paper.\n";
	cerr << "  Rem: Only this option will be taken into account and no further computation done.\n";
	
	exit(1);
}

int main(int argc, char **argv){
	LargeDefGradLagrange LargeDef;
	bool ok;
	
	// Check command line
	if (argc < 3) {
		usage();
	}
	
	// Read the name of input and output images (= 1st channel)
	strcpy(LargeDef.SourceFiles[LargeDef.NbChannels],argv[1]);
	argc--;  argv++;
	strcpy(LargeDef.TargetFiles[LargeDef.NbChannels],argv[1]);
	argc--;  argv++;
	LargeDef.NbChannels++;
	
	// Parse remaining parameters
	while (argc > 1) {
		//1 - Primary options
		ok = false;
		if ((ok == false) && (strcmp(argv[1], "-iterations") == 0)) {
			argc--; argv++;
			LargeDef.iteration_nb = atoi(argv[1]);
			argc--; argv++;
			ok = true;
		}
		if ((ok == false) && (strcmp(argv[1], "-subdivisions") == 0)) {
			argc--; argv++;
			LargeDef.NbTimeSubdiv = atoi(argv[1]);
			argc--; argv++;
			ok = true;
		}
		if ((ok == false) && (strcmp(argv[1], "-MaxVeloUpdt") == 0)) {
			argc--; argv++;
			LargeDef.MaxVelocityUpdate = atof(argv[1]);
			argc--; argv++;
			ok = true;
		}
		//2 - Inputs and Outputs
		if ((ok == false) && (strcmp(argv[1], "-PrefixInputs") == 0)) {
			argc--; argv++;
			strcpy(LargeDef.PrefixInputs,argv[1]);
			argc--; argv++;
			ok = true;
		}
		if ((ok == false) && (strcmp(argv[1], "-PrefixOutputs") == 0)) {
			argc--; argv++;
			strcpy(LargeDef.PrefixOutputs,argv[1]);
			argc--; argv++;
			ok = true;
		}
		if ((ok == false) && (strcmp(argv[1], "-AddChannel") == 0) && (LargeDef.NbChannels<100)) {
			argc--; argv++;
			LargeDef.weightChannel[LargeDef.NbChannels] = atof(argv[1]);
			argc--; argv++;
			strcpy(LargeDef.SourceFiles[LargeDef.NbChannels],argv[1]);
			argc--; argv++;
			strcpy(LargeDef.TargetFiles[LargeDef.NbChannels],argv[1]);
			argc--; argv++;
			LargeDef.NbChannels++;
			ok = true;
			if (LargeDef.NbChannels==100) cout << "\n \n MAXIMUM NUMBER OF 100 CHANNELS IS REACHED !!!\n \n ";
		}
		if ((ok == false) && (strcmp(argv[1], "-Mask") == 0)) {
			argc--; argv++;
			strcpy(LargeDef.MaskFile,argv[1]);
			argc--; argv++;
			ok = true;
		}
		//3 - Kernels
		if ((ok == false) && (strcmp(argv[1], "-Gauss") == 0)) {
			LargeDef.weight1 = 100.;
			argc--; argv++;
			LargeDef.sigmaX1 = atof(argv[1]);
			LargeDef.sigmaY1 = atof(argv[1]);
			LargeDef.sigmaZ1 = atof(argv[1]);
			argc--; argv++;
			LargeDef.NbKernels=1;
			ok = true;
		}
		if ((ok == false) && (strcmp(argv[1], "-A_Gauss") == 0)) {
			LargeDef.weight1 = 100.;
			argc--; argv++;
			LargeDef.sigmaX1 = atof(argv[1]);
			argc--; argv++;
			LargeDef.sigmaY1 = atof(argv[1]);
			argc--; argv++;
			LargeDef.sigmaZ1 = atof(argv[1]);
			argc--; argv++;
			LargeDef.NbKernels=1;
			ok = true;
		}
		if ((ok == false) && (strcmp(argv[1], "-M_Gauss") == 0)) {
			argc--; argv++;
			LargeDef.NbKernels= atoi(argv[1]);
			argc--; argv++;
			if (LargeDef.NbKernels>=1){
				LargeDef.weight1 = atof(argv[1]);
				argc--; argv++;
				LargeDef.sigmaX1 = atof(argv[1]);
				LargeDef.sigmaY1 = atof(argv[1]);
				LargeDef.sigmaZ1 = atof(argv[1]);
				argc--; argv++;
			}
			if (LargeDef.NbKernels>=2){
				LargeDef.weight2 = atof(argv[1]);
				argc--; argv++;
				LargeDef.sigmaX2 = atof(argv[1]);
				LargeDef.sigmaY2 = atof(argv[1]);
				LargeDef.sigmaZ2 = atof(argv[1]);
				argc--; argv++;
			}
			if (LargeDef.NbKernels>=3){
				LargeDef.weight3 = atof(argv[1]);
				argc--; argv++;
				LargeDef.sigmaX3 = atof(argv[1]);
				LargeDef.sigmaY3 = atof(argv[1]);
				LargeDef.sigmaZ3 = atof(argv[1]);
				argc--; argv++;
			}
			if (LargeDef.NbKernels>=4){
				LargeDef.weight4 = atof(argv[1]);
				argc--; argv++;
				LargeDef.sigmaX4 = atof(argv[1]);
				LargeDef.sigmaY4 = atof(argv[1]);
				LargeDef.sigmaZ4 = atof(argv[1]);
				argc--; argv++;
			}
			if (LargeDef.NbKernels>=5){
				LargeDef.weight5 = atof(argv[1]);
				argc--; argv++;
				LargeDef.sigmaX5 = atof(argv[1]);
				LargeDef.sigmaY5 = atof(argv[1]);
				LargeDef.sigmaZ5 = atof(argv[1]);
				argc--; argv++;
			}
			if (LargeDef.NbKernels>=6){
				LargeDef.weight6 = atof(argv[1]);
				argc--; argv++;
				LargeDef.sigmaX6 = atof(argv[1]);
				LargeDef.sigmaY6 = atof(argv[1]);
				LargeDef.sigmaZ6 = atof(argv[1]);
				argc--; argv++;
			}
			if (LargeDef.NbKernels>=7){
				LargeDef.weight7 = atof(argv[1]);
				argc--; argv++;
				LargeDef.sigmaX7 = atof(argv[1]);
				LargeDef.sigmaY7 = atof(argv[1]);
				LargeDef.sigmaZ7 = atof(argv[1]);
				argc--; argv++;
			}
			ok = true;
		}
		if ((ok == false) && (strcmp(argv[1], "-M_A_Gauss") == 0)) {
			argc--; argv++;
			LargeDef.NbKernels= atoi(argv[1]);
			argc--; argv++;
			if (LargeDef.NbKernels>=1){
				LargeDef.weight1 = atof(argv[1]);
				argc--; argv++;
				LargeDef.sigmaX1 = atof(argv[1]);
				argc--; argv++;
				LargeDef.sigmaY1 = atof(argv[1]);
				argc--; argv++;
				LargeDef.sigmaZ1 = atof(argv[1]);
				argc--; argv++;
			}
			if (LargeDef.NbKernels>=2){
				LargeDef.weight2 = atof(argv[1]);
				argc--; argv++;
				LargeDef.sigmaX2 = atof(argv[1]);
				argc--; argv++;
				LargeDef.sigmaY2 = atof(argv[1]);
				argc--; argv++;
				LargeDef.sigmaZ2 = atof(argv[1]);
				argc--; argv++;
			}
			if (LargeDef.NbKernels>=3){
				LargeDef.weight3 = atof(argv[1]);
				argc--; argv++;
				LargeDef.sigmaX3 = atof(argv[1]);
				argc--; argv++;
				LargeDef.sigmaY3 = atof(argv[1]);
				argc--; argv++;
				LargeDef.sigmaZ3 = atof(argv[1]);
				argc--; argv++;
			}
			if (LargeDef.NbKernels>=4){
				LargeDef.weight4 = atof(argv[1]);
				argc--; argv++;
				LargeDef.sigmaX4 = atof(argv[1]);
				argc--; argv++;
				LargeDef.sigmaY4 = atof(argv[1]);
				argc--; argv++;
				LargeDef.sigmaZ4 = atof(argv[1]);
				argc--; argv++;
			}
			if (LargeDef.NbKernels>=5){
				LargeDef.weight5 = atof(argv[1]);
				argc--; argv++;
				LargeDef.sigmaX5 = atof(argv[1]);
				argc--; argv++;
				LargeDef.sigmaY5 = atof(argv[1]);
				argc--; argv++;
				LargeDef.sigmaZ5 = atof(argv[1]);
				argc--; argv++;
			}
			if (LargeDef.NbKernels>=6){
				LargeDef.weight6 = atof(argv[1]);
				argc--; argv++;
				LargeDef.sigmaX6 = atof(argv[1]);
				argc--; argv++;
				LargeDef.sigmaY6 = atof(argv[1]);
				argc--; argv++;
				LargeDef.sigmaZ6 = atof(argv[1]);
				argc--; argv++;
			}
			if (LargeDef.NbKernels>=7){
				LargeDef.weight7 = atof(argv[1]);
				argc--; argv++;
				LargeDef.sigmaX7 = atof(argv[1]);
				argc--; argv++;
				LargeDef.sigmaY7 = atof(argv[1]);
				argc--; argv++;
				LargeDef.sigmaZ7 = atof(argv[1]);
				argc--; argv++;
			}
			ok = true;
		}
		
		
		
		//4 - Secondary options
		if ((ok == false) && (strcmp(argv[1], "-GreyLevAlign") == 0)) {
			argc--; argv++;
			LargeDef.GreyLevAlign = 1;
			LargeDef.GLA_Padding_Src = atof(argv[1]);
			argc--; argv++;
			LargeDef.GLA_Padding_Trg = atof(argv[1]);
			argc--; argv++;
			ok = true;
		}
		if ((ok == false) && (strcmp(argv[1], "-SplitKernels") == 0)) {
			argc--; argv++;
			LargeDef.SplitKernels = 1;
			ok = true;
		}
		if ((ok == false) && (strcmp(argv[1], "-symmetric") == 0)) {
			argc--; argv++;
			LargeDef.symmetric = 1;
			ok = true;
		}
		if ((ok == false) && (strcmp(argv[1], "-margins") == 0)) {
			argc--; argv++;
			LargeDef.Margin = atoi(argv[1]);
			argc--; argv++;
			ok = true;
		}
		if ((ok == false) && (strcmp(argv[1], "-WghtVeloField") == 0)) {
			argc--; argv++;
			LargeDef.WghtVelField = atof(argv[1]);
			argc--; argv++;
			ok = true;
		}
		if ((ok == false) && (strcmp(argv[1], "-RefMaxGrad") == 0)) {
			argc--; argv++;
			LargeDef.RefMaxGrad = atof(argv[1]);
			argc--; argv++;
			ok = true;
		}
		if ((ok == false) && (strcmp(argv[1], "-epsilon") == 0)) {
			argc--; argv++;
			LargeDef.epsilon = atof(argv[1]);
			argc--; argv++;
			ok = true;
		}
		if ((ok == false) && (strcmp(argv[1], "-MapSrcImag") == 0)) {
			argc--; argv++;
			strcpy(LargeDef.MappingSrc,argv[1]);
			argc--; argv++;
			ok = true;
		}
		if ((ok == false) && (strcmp(argv[1], "-MapTrgImag") == 0)) {
			argc--; argv++;
			strcpy(LargeDef.MappingTrg,argv[1]);
			argc--; argv++;
			ok = true;
		}
		if ((ok == false) && (strcmp(argv[1], "-TypicInvAmp") == 0)) {
			argc--; argv++;
			LargeDef.sigmaX1 = atof(argv[1]);
			argc--; argv++;
			LargeDef.sigmaY1 = atof(argv[1]);
			argc--; argv++;
			LargeDef.sigmaZ1 = atof(argv[1]);
			argc--; argv++;
			LargeDef.NbKernels=1;
			LargeDef.MeasureTypicAmp=1;
			ok = true;
		}
		//5) Special outputs
		if ((ok == false) && (strcmp(argv[1], "-AOD") == 0)) {
			argc--; argv++;
			LargeDef.FlowLength = 1;
			ok = true;
		}
		if ((ok == false) && (strcmp(argv[1], "-DetJacobian") == 0)) {
			argc--; argv++;
			LargeDef.DetJacobian = 1;
			ok = true;
		}
		if ((ok == false) && (strcmp(argv[1], "-FinalDefVec") == 0)) {
			argc--; argv++;
			LargeDef.FinalDefVec = 1;
			ok = true;
		}
		if ((ok == false) && (strcmp(argv[1], "-InitMomentum") == 0)) {
			argc--; argv++;
			LargeDef.CptInitMomentum = 1;
			ok = true;
		}
		if ((ok == false) && (strcmp(argv[1], "-FinalDefInvVec") == 0)) {
			argc--; argv++;
			LargeDef.FinalDefInvVec = 1;
			ok = true;
		}
		if ((ok == false) && (strcmp(argv[1], "-ShowSSD") == 0)) {
			argc--; argv++;
			LargeDef.ShowSSD = 1;
			ok = true;
		}
		if (ok == false) usage();
	}
	
	//run process
	cout << "Large Deformation registration using Beg 05's technique ... \n"; cout.flush();
	LargeDef.Run();
	cout << "done" << endl;
	
	return 0;
}
