/*=========================================================================

  Library   : Image Registration Toolkit (IRTK)
  Module    : $Id$
  Copyright : Imperial College, Department of Computing
              Visual Information Processing (VIP), 2011 onwards
  Date      : $Date$
  Version   : $Revision$
  Changes   : $Author$

=========================================================================*/



/*
 * 	This class implements the Paper:
 * 	Neonatal brain segmentation using second order neighboorhood information, Ledig C. et al
 * 	PaPI Workshop 2012 in conjunction with MICCAI 2012, Nice *
 *
 */

#include <irtkImage.h>
#include <irtkEMClassification.h>
#include <irtkPolynomialBiasField.h>
#include <irtkEMClassification2ndOrderMRF.h>
#include <irtkGaussian.h>
#include <irtkMatrix.h>
#include <iostream>
#include <fstream>
#include <sstream>

#include <irtkConvolutionWithGaussianDerivative.h>

char *output_segmentation, *output_biascorrection, *output_biasfield, *connections, *output_pvsegmentation, *mask, *posteriors, *connections_2nd;

void readConnectivity2nd(vector< vector< pair< pair<int,int>, double > > > &G2nd, char* filename );
void readConnectivityMatrix( irtkMatrix* mat, char* filename);

void usage()
{
	cerr << "Usage: emSeg2ndOrderMRF [image] [n] [atlas 1 ... atlas n] [segmented image (output)] <options>" << endl;
	cerr << "where <options> is one or more of the following:\n" << endl;
	cerr << "<-mrf file>" 						<<	"\t\t\t\t\tuse MRF, file contains a (n+1)x(n+1) (+1 because of background, n if -nobackground is set) connection matrix defining with entry (i,j) if class i is adjacent (i,j)=1, distant (i,j)=2 or identic (i==j) to class j. " << endl;
	cerr << "<-mrf2 file>" 						<<	"\t\t\t\t\tuse 2nd order neighborhood for MRF energy, file contains a sparse representation of weights. Each line has the form j k l G_jkl, where class k and l are the neigbors of j. By default G_jkl=0." << endl;
	cerr << "<-biasfielddegree number>" 		<< 	"\t\t\tpolynomial degree (in one dimension) of biasfield (default = 0)" << endl;
	cerr << "<-biasfield file>"    				<<	"\t\t\t\tsave final bias field (for log transformed intensities) to file. " << endl;
	cerr << "<-correctedimage file>"			<< 	"\t\t\tsave bias corrected image (biasfielddegree should be >0)]." << endl;
	cerr << "<-outputpv file>"    				<<	"\t\t\t\tsave hard segmentation containing partial volume classes to file " << endl;
	cerr << "<-norelax>"    					<<	"\t\t\t\t\tdo not relax priors after first convergence " << endl;
	cerr << "<-mask filename>"					<<	"\t\t\t\tmask" << endl;
	cerr << "<-nolog>"							<<	"\t\t\t\t\tdo not log transform data" << endl;
	cerr << "<-nobackground>"					<<	"\t\t\t\t\tdo not add background class" << endl;
	cerr << "<-normalize>"						<<	"\t\t\t\t\tnormalize atlas" << endl;
	cerr << "<-posteriors path>"				<<	"\t\t\t\tsave posteriors to path" << endl;
	cerr << "<-pv number1 number2>"  			<<	"\t\t\t\tadd partial volume class between class number1 and number2" << endl;
	cerr << "<-mrfweights number1 number2>"  	<<	"\t\t\tspecify MRF weights of adjacent (num1) and distant (num2) classes, default (0.15/0.75)" << endl;
	cerr << "<-relaxfactor number1>"  			<<	"\t\t\t\tspecify priors relaxation factor in [0,1], 0=no relaxation, default 0.5" << endl;
	cerr << "<-relaxruns number1>"  			<<	"\t\t\t\tnumber of relaxation runs, default 1" << endl;
	cerr << "<-relaxruns_2nd number1>"  		<<	"\t\t\tnumber of relaxation runs dependent on 2nd order connectivity, default 0" << endl;
	cerr << "<-relaxfactor_2nd number1>"  		<<	"\t\t\trelaxation factor 0=no relaxation, default=0.5" << endl;
	cerr << "<-refine number1>"					<<	"\t\t\t\trefine final segmentation where number1 is the index of the background label. default=no refinement" << endl;
	cerr << "<-addbackground number>"			<<	"\t\t\t\tadds background label, prior is the undefined part of the other priors with number is the max prior for a voxel" << endl;
	exit(1);
}

int main(int argc, char **argv)
{
  int i, n, ok, padding, maxIterations;

  // Default parameters
  int biasfield_degree = 0;
  bool nolog = false;
  bool normalize = false;
  output_biasfield = NULL;
  connections = NULL;
  connections_2nd = NULL;

  maxIterations = 50;
  padding    = MIN_GREY;

  bool nomrf = true;

  bool addbg = false;
  double max_prior = 1;

  bool refine_seg = false;

  int refine_label = -1;

  int relax_runs = 1;
  double relax_factor = 0.5;

  int relax_runs_2nd = 0;
  double relax_factor_2nd = 0.5;

  double mrf_weight_adjacent = 0.15;
  double mrf_weight_distant = 0.75;

  int no_mrf_correction = 1;
  int no_background = 0;
  int no_relax = 0;

  vector< pair<int,int> > pv_classes;

    // File name for bias corrected image
  output_biascorrection = NULL;

  output_pvsegmentation = NULL;

  posteriors = NULL;
  mask = NULL;

  if (argc < 4) {
    usage();
    exit(1);
  }

  // Input image
  irtkRealImage image;
  image.Read(argv[1]);

  image.Print();
  argc--;
  argv++;

  // Number of tissues
  n = atoi(argv[1]);
  argc--;
  argv++;

  // Probabilistic atlas
  irtkRealImage **atlas = new irtkRealImage*[n];
  irtkRealImage *background=NULL;

  // Read atlas for each tissue
  for (i = 0; i < n; i++) {
    atlas[i] = new irtkRealImage;
    atlas[i]->Read(argv[1]);
    atlas[i]->Print();
    cerr << "Image " << i <<" = " << argv[1] <<endl;


    // resize atlases to image size
    // hack! This is not good and should not happen!
    if( atlas[i]->GetX() != image.GetX() ||  atlas[i]->GetY() != image.GetY() ||  atlas[i]->GetZ() != image.GetZ() )
    {
      cerr << "WARNING: atlas has different dimension from image!" << endl;
  	  irtkRealImage* atlas_new = new irtkRealImage();
  	  atlas_new->Initialize(image.GetImageAttributes());

	  for( int j = 0; j < image.GetX(); ++j )
	  {
		  for( int k = 0; k < image.GetY(); ++k )
		  {
		  	  for( int l = 0; l < image.GetZ(); ++l )
		  	  {
		  		  double x = j;
		  		  double y = k;
		  		  double z = l;
		  		  atlas_new->ImageToWorld(x,y,z);
		  		  atlas[i]->WorldToImage(x,y,z);
  				  if( x < atlas[i]->GetX() && y < atlas[i]->GetY() && z < atlas[i]->GetZ() \
  						  && x > 0 && y > 0 && z > 0 )
  					atlas_new->Put(j,k,l,atlas[i]->Get(x,y,z));
  				  else
					atlas_new->Put(j,k,l,0);
  			  }
  		  }
  	  }
  	  delete atlas[i];
  	  atlas[i] = atlas_new;
    }

    argc--;
    argv++;
  }

  // File name for segmentation
  output_segmentation = argv[1];
  argc--;
  argv++;

  // Parse remaining parameters
  while (argc > 1) {
    ok = false;

    if ((ok == false) && (strcmp(argv[1], "-correctedimage") == 0)) {
      argc--;
      argv++;
      output_biascorrection = argv[1];
      argc--;
      argv++;
      ok = true;
    }
    if ((ok == false) && (strcmp(argv[1], "-mrf") == 0)) {
      argc--;
      argv++;
      connections = argv[1];
      nomrf = false;
      argc--;
      argv++;
      ok = true;
    }
    if ((ok == false) && (strcmp(argv[1], "-mrf2") == 0)) {
      argc--;
      argv++;
      cout << "Using 2nd_order connectivity information" << endl;
      connections_2nd = argv[1];
      argc--;
      argv++;
      ok = true;
    }
    if ((ok == false) && (strcmp(argv[1], "-nolog") == 0)) {
      argc--;
      argv++;
      nolog = true;
      ok = true;
    }
    if ((ok == false) && (strcmp(argv[1], "-outputpv") == 0)) {
      argc--;
      argv++;
      output_pvsegmentation = argv[1];
      argc--;
      argv++;
      ok = true;
    }
    if ((ok == false) && (strcmp(argv[1], "-relaxfactor") == 0)) {
      argc--;
      argv++;
      relax_factor = atof(argv[1]);
      argc--;
      argv++;
      cout << "Using relaxation factor " << relax_factor << endl;
      ok = true;
    }
    if ((ok == false) && (strcmp(argv[1], "-relaxruns") == 0)) {
      argc--;
      argv++;
      relax_runs = atof(argv[1]);
      argc--;
      argv++;
      cout << "Relaxing priors " << relax_runs << " times" << endl;
      ok = true;
    }
    if ((ok == false) && (strcmp(argv[1], "-refine") == 0)) {
      argc--;
      argv++;
      refine_seg = true;
      refine_label = atoi(argv[1]);
      argc--;
      argv++;
      cout << "Refining final segmentation with background label: " << refine_label << endl;
      ok = true;
    }
    if ((ok == false) && (strcmp(argv[1], "-relaxfactor_2nd") == 0)) {
      argc--;
      argv++;
      relax_factor_2nd = atof(argv[1]);
      argc--;
      argv++;
      cout << "Using relaxation factor for 2nd order " << relax_factor_2nd << endl;
      ok = true;
    }
    if ((ok == false) && (strcmp(argv[1], "-relaxruns_2nd") == 0)) {
      argc--;
      argv++;
      relax_runs_2nd = atof(argv[1]);
      argc--;
      argv++;
      cout << "Relaxing priors " << relax_runs << " times dependent on 2nd order information" << endl;
      ok = true;
    }
    if ((ok == false) && (strcmp(argv[1], "-mrfweights") == 0)) {
      argc--;
      argv++;
      mrf_weight_adjacent = atof(argv[1]);
      argc--;
      argv++;
      mrf_weight_distant = atof(argv[1]);
      argc--;
      argv++;
      cout << "Using MRF weights " << mrf_weight_adjacent << "/" << mrf_weight_distant << endl;
      ok = true;
    }
    if ((ok == false) && (strcmp(argv[1], "-posteriors") == 0)) {
      argc--;
      argv++;
      posteriors = argv[1];
      cout << "posteriors will be saved at: " << posteriors << endl;
      argc--;
      argv++;
      ok = true;
    }
    if ((ok == false) && (strcmp(argv[1], "-iterations") == 0)) {
      argc--;
      argv++;
      maxIterations = atoi(argv[1]);
      argc--;
      argv++;
      ok = true;
    }
    if ((ok == false) && (strcmp(argv[1], "-norelax") == 0)) {
      cout << "no prior relaxation" << endl;
      argc--;
      argv++;
      no_relax = true;
      ok = true;
    }
    if ((ok == false) && (strcmp(argv[1], "-normalize") == 0)) {
      cout << "normalizing atlas" << endl;
      argc--;
      argv++;
      normalize = true;
      ok = true;
    }
    if ((ok == false) && (strcmp(argv[1], "-nobackground") == 0)) {
      cout << "no background" << endl;
      argc--;
      argv++;
      no_background = true;
      ok = true;
    }
    if ((ok == false) && (strcmp(argv[1], "-mask") == 0)) {
    	cout << "using manual mask! padding has no effect!" << endl;
        argc--;
        argv++;
        mask = argv[1];
        argc--;
        argv++;
        ok = true;
    }
    if ((ok == false) && (strcmp(argv[1], "-pv") == 0)) {
      argc--;
      argv++;
      int a, b;
      a = atoi(argv[1]);
      argc--;
      argv++;
      b = atoi(argv[1]);
      argc--;
      argv++;
      pv_classes.push_back(make_pair(a,b));
      cout << "will add partial volume class between " << a << " and " << b << endl;
      ok = true;
    }
    if ((ok == false) && (strcmp(argv[1], "-padding") == 0)) {
      argc--;
      argv++;
      padding = atoi(argv[1]);
      argc--;
      argv++;
      ok = true;
    }
    if ((ok == false) && (strcmp(argv[1], "-biasfield") == 0)) {
      argc--;
      argv++;
      output_biasfield = argv[1];
      cout << "Output biasfield to: " << output_biasfield << endl;
      argc--;
      argv++;
      ok = true;
    }
    if ((ok == false) && (strcmp(argv[1], "-biasfielddegree") == 0)) {
      argc--;
      argv++;
      biasfield_degree = atoi(argv[1]);
      cout << "Degree of biasfield polynomial: " << biasfield_degree << endl;
      argc--;
      argv++;
      ok = true;
    }
    if ((ok == false) && (strcmp(argv[1], "-addbackground") == 0)) {
      argc--;
      argv++;
      background = new irtkRealImage;
      *background = *atlas[0];
      cout << "Adding background" << endl;
      cout << "Normalizing priors" << endl;
      normalize = true;
      addbg = true;
      no_background = false;
      ok = true;
    }
    if ((ok == false) && (strcmp(argv[1], "-nomrf") == 0)) {
      argc--;
      argv++;
      nomrf = true;
      cout << "NO MRF " << endl;
      ok = true;
    }
    if (ok == false) {
      cerr << "Can not parse argument " << argv[1] << endl;
      usage();
    }

  }

  if( connections_2nd != NULL && connections == NULL )
  {
	  cerr << "MRF must be definited (-mrf filename) if 2nd order MRF is used!" << endl;
	  exit(1);
  }

  if( normalize )
  {
	  cout << "Normalising atlas...";
	  // normalize atlas
	  for( int i = 0; i < atlas[0]->GetNumberOfVoxels(); ++i )
	  {
		  double denom = 0;
		  for( int l = 0; l < n; ++l )
		  {
			  irtkRealPixel* ptr=atlas[l]->GetPointerToVoxels();
			  denom += ptr[i];
		  }
		  if( denom > max_prior ) max_prior = denom;
	  }
	  for( int i = 0; i < atlas[0]->GetNumberOfVoxels(); ++i )
	  {
		  double denom = 0;
		  for( int l = 0; l < n; ++l )
		  {
			  irtkRealPixel* ptr=atlas[l]->GetPointerToVoxels();
			  denom += ptr[i];
		  }
		  if (addbg)
		  {
			  if (denom > max_prior ) cout << denom << endl;
			  irtkRealPixel* bgptr = background->GetPointerToVoxels();
			  bgptr[i] = (max_prior - denom) / max_prior;
			  if( bgptr[i] < 0.2 )
			  {
				  bgptr[i] = 0;
			  }
			  else
			  {
				  denom = max_prior;
			  }
		  }
		  if ( denom )
		  {
			  for( int l = 0; l < n; ++l )
			  {
				  irtkRealPixel* ptr=atlas[l]->GetPointerToVoxels();
				  ptr[i] /= denom;
			  }
		  }
		  else
		  {
			  for( int l = 0; l < n; ++l )
			  {
				  irtkRealPixel* ptr=atlas[l]->GetPointerToVoxels();
				  if( l == 0 )
					  ptr[i] = 1.0;
				  else
					  ptr[i] = 0.0;
			  }
		  }
	  }
	  cout << "done!" << endl;
  }

  // logtransform image
  // be careful log transformation might transform intensities to the actual padding! --> change padding value
  int min_log_voxel = 1000;
  irtkRealImage padding_map = image;

  for( int i = 0; i < image.GetX(); ++i)
  {
	  for( int j = 0; j < image.GetY(); ++j)
	  {
		  for( int k = 0; k < image.GetZ(); ++k)
		  {
			  if(image.Get(i,j,k) == padding )
			  {
				  padding_map.Put(i,j,k, 1);
			  }
			  else
			  {
				  padding_map.Put(i,j,k, 0);
			  }
			  if( !nolog )
			  {
				  // no log transformation possible for 0 intensity values, tread like padding voxels
				  if( image.Get(i,j,k) == 0 )
				  {
					  image.Put(i,j,k, padding);
					  padding_map.Put(i,j,k,1);
				  }

				  if( image.Get(i,j,k) > 0 && image.Get(i,j,k) != padding )
				  {
					  image.Put(i,j,k, log(image.Get(i,j,k)));
					  min_log_voxel = min(min_log_voxel,  static_cast<int>( image.Get(i,j,k)) );
				  }
			  }
		  }
	  }
  }
  int original_padding = padding;

  if ( !nolog )
  {
	  padding = min_log_voxel-1;

	  for( int i = 0; i < image.GetX(); ++i)
	  {
		  for( int j = 0; j < image.GetY(); ++j)
		  {
			  for( int k = 0; k < image.GetZ(); ++k)
			  {
				  if( padding_map.Get(i,j,k) == 1 )
				  {
					  image.Put(i,j,k, padding);
				  }
			  }
		  }
	  }
  }

   // Create classification object
  irtkEMClassification2ndOrderMRF* classification;//(n, atlas);
  if( no_background )
  {
	  classification = new irtkEMClassification2ndOrderMRF(n, atlas);
  }
  else
  {
	  classification = new irtkEMClassification2ndOrderMRF(n, atlas, background);
	  n = n+1;
  }

  irtkMatrix* G;
  if( connections != NULL )
  {
	  G = new irtkMatrix(n,n);
	  readConnectivityMatrix(G, connections);
	  //G->Print();
	  classification->SetInput(image, *G);
	  no_mrf_correction = 0;
  }
  else
  {
	  // no MRF correction if theres no MRF
	  no_mrf_correction = 1;
	  G = new irtkMatrix(1,1);
	  classification->SetInput(image, *G);
  }

  vector< vector< pair< pair<int,int>, double > > > G2nd;
  G2nd.resize(n);

  if( connections_2nd != NULL )
  {
	  readConnectivity2nd(G2nd, connections_2nd);
	  classification->SetMRF_2nd_order(G2nd);
  }

  classification->SetPadding(padding);
  if( !nolog )
  {
	  classification->SetLogTransformed(true);
  }

  classification->CreateMask();
  classification->Initialise();
  classification->SetMRFWeights(mrf_weight_adjacent, mrf_weight_distant);
  classification->SetRelaxationFactor( relax_factor );
  classification->SetRelaxationFactor_2nd( relax_factor_2nd );

  irtkRealImage maskImage;

  if ( mask != NULL )
  {
	  maskImage.Read(mask);
	  irtkRealPixel* ptr = maskImage.GetPointerToVoxels();

	  for( int i = 0; i < maskImage.GetNumberOfVoxels(); ++i )
	  {
		  if( ptr[i] > 0 ) *ptr = 1;
		  else *ptr = 0;
	  }
	  classification->SetMask(maskImage);
  }

  double rel_diff = 1.0;
  bool stop = false;
  int curr_biasfield_degree = 1;
  int iter = 0;
  int improvePhase = 0;

  // Create bias field
  irtkPolynomialBiasField *biasfield = new irtkPolynomialBiasField(image, curr_biasfield_degree);
  classification->SetBiasField(biasfield);

  bool BFupdate = false;
  bool MRFupdate = false;
  vector<int> pv_positions;
  int number_current_iterations = 0;
  int gyri_sulci = 0;
  bool PVon = false;

  while( !stop && iter < maxIterations )
  {
	  if( !MRFupdate || nomrf )
		  classification->EStep();
	  else if (connections_2nd)
		  classification->EStepMRF_2nd_order();
	  else
	  	  classification->EStepMRF();

	  if( BFupdate ) {
		  classification->WStep();
		  classification->BStep();
	  }

	  if( !PVon )
	  {
		  classification->MStep();
	  }
	  else
	  {
		  //	Alternative cardoso implementation, not necessarily better!!
		  classification->MStepPV();
	  }
	  classification->Print();
	  rel_diff = classification->LogLikelihood();

	  if( rel_diff > 0.05 )
	  {
		  gyri_sulci++;
	  }

	  if( rel_diff < 0.005  && iter < maxIterations && number_current_iterations )
	  {
		  switch( improvePhase )
		  {
		  case 0:
			  if( curr_biasfield_degree < biasfield_degree && number_current_iterations)
			  {
				  curr_biasfield_degree++;
				  delete biasfield;
				  biasfield = new irtkPolynomialBiasField(image, curr_biasfield_degree);
				  classification->SetBiasField(biasfield);

				  BFupdate = true;
				  if( curr_biasfield_degree == biasfield_degree )
				  {
					  improvePhase++;
				  }
				  break;
			  }
			  else
			  {
				  improvePhase++;
			  }
		  case 1:
			  if( !MRFupdate && !no_mrf_correction && !nomrf )
			  {
				  MRFupdate = true;
				  improvePhase++;
				  break;
			  }
			  else
			  {
				  improvePhase++;
			  }
		  case 2:
			  if( relax_runs_2nd )
			  {
				  // relax priors
				  cout << "relaxing priors (2nd order) now..." << endl;
				  classification->RStep_2nd_order();
				  cout << "done!" << endl;
				  relax_runs_2nd--;
				  if( relax_runs_2nd == 0 )
				  {
						  improvePhase++;
				  }
				  break;
			  }
			  else
			  {
				  improvePhase++;
			  }
		  case 3:
			  if( pv_classes.size() && !PVon )
			  {
				  for( unsigned int i = 0; i < pv_classes.size(); ++i )
				  {
					  int a, b;
					  pair<int, int> tmp = pv_classes[i];
					  a = tmp.first;
					  b = tmp.second;

					  // add partial volume classes
					  cout << "adding partial volume classes between " << a << " " << b << endl;
					  pv_positions.push_back( classification->AddPartialVolumeClass( a, b) );
					  cout << "New PV Class at position: " << pv_positions[i] << endl;
					  cout << "done!" << endl;
				  }
				  PVon = true;
				  improvePhase++;
				  break;
			  }
			  else
			  {
				  improvePhase++;
			  }
		  case 4:
			  if( !no_relax && relax_runs )
			  {
				  // relax priors
				  cout << "relaxing priors now..." << endl;
				  classification->RStep();
				  cout << "done!" << endl;
				  relax_runs--;
				  if( relax_runs == 0 )
				  {
						  improvePhase++;
				  }
				  break;
			  }
			  else
			  {
				  improvePhase++;
			  }
		  case 5:
			  stop = true;
			  break;
		  }
		  number_current_iterations = 0;

	  }
	  else
	  {
		  number_current_iterations++;
	  }
	  iter++;
  }

  irtkRealImage corrected_image = image;
  irtkRealImage biasfield_image = image;

  // Bias corrected image
  classification->GetBiasCorrectedImage(corrected_image);

  // Save bias field
  if (output_biasfield != NULL) {
	cout << "Writing biasfield " << output_biasfield;
    classification->GetBiasField( biasfield_image );
    biasfield_image.Write(output_biasfield);
    cout << "    done!!!!" << endl;
  }

  if( output_biascorrection != NULL )
  {
	  for( int k = 0; k < corrected_image.GetZ(); ++k)
	  {
		  for( int j = 0; j < corrected_image.GetY(); ++j)
		  {
			  for( int i = 0; i < corrected_image.GetX(); ++i)
			  {
				  if( !nolog )
				  {
					  if( corrected_image.Get(i,j,k) == padding )
					  {
						  corrected_image.Put(i,j,k, original_padding);
					  }
					  else
					  {
						  corrected_image.Put(i,j,k, exp(corrected_image.Get(i,j,k) ));
					  }
				  }
			  }
		  }
	  }
	  corrected_image.Write(output_biascorrection);
  }

  // Save segmentation
  if( refine_seg )
  {
	  classification->RefineSegmentation(refine_label);
  }

  classification->ConstructSegmentation(corrected_image);
  corrected_image.Write(output_segmentation);

  if( pv_classes.size() )
  {
	  classification->removePVclasses();
  }
  if( output_pvsegmentation != NULL )
  {
	  classification->ConstructSegmentation(corrected_image);
	  corrected_image.Write(output_pvsegmentation);
  }

  irtkImageAttributes attr = image.GetImageAttributes();
  irtkRealImage probMap(attr);

  if( posteriors != NULL )
  {
	  for( int i = 0; i < n; ++i )
	  {
		  classification->GetProbMapMasked(i,probMap);
		  if( posteriors[strlen(posteriors)-1] == '\n' ) posteriors[strlen(posteriors)-1] = '\0';
		  stringstream str;
		  str << posteriors << "posteriors_" << i << ".nii.gz";
		  cout << str.str().c_str() << endl;
		  probMap.GetFrame(0).Write(str.str().c_str());
	  }
  }

  //classification->WriteGaussianParameters("parameters.txt");
  delete G;
  delete classification;
}

void readConnectivity2nd(vector< vector< pair< pair<int,int>, double > > > &G2nd, char* filename )
{
	  // Open file
	  ifstream from;
	  from.open(filename);

	  int j,k,l,w;
	  while( !from.eof() )
	  {
		  from >> j >> k >> l >> w;
		  if( from.eof() ) break;
		  cout << "2nd order MRF: Penalising neighbors of " << j << ": " << k << ", " << l << " with weight: " << w << endl;
		  G2nd[j].push_back(make_pair(make_pair(k,l), w));
	  }
	  from.close();
}

void readConnectivityMatrix( irtkMatrix* mat, char* filename)
{
	  // Open file
	  ifstream from;
	  from.open(filename);

	  int value;
	  for( int i = 0; i < mat->Rows(); ++i )
	  {
		  for( int j = 0; j < mat->Cols(); ++j )
		  {
			  from >> value;
			  // distant = 2, close = 1, same = 0
			  if( i == j )
				  mat->Put(i,j,0);
			  else if( value == 1 )
				  mat->Put(i,j,1);
			  else
				  mat->Put(i,j,2);
		  }
	  }
	  from.close();
}
