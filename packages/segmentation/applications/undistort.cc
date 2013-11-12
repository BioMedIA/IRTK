/*=========================================================================

  Library   : Image Registration Toolkit (IRTK)
  Module    : $Id: reconstructionb0.cc 998 2013-10-15 15:24:16Z mm3 $
  Copyright : Imperial College, Department of Computing
              Visual Information Processing (VIP), 2011 onwards
  Date      : $Date: 2013-10-15 16:24:16 +0100 (Tue, 15 Oct 2013) $
  Version   : $Revision: 998 $
  Changes   : $Author: mm3 $

=========================================================================*/

#include <irtkImage.h>
#include <irtkTransformation.h>
#include <irtkReconstruction.h>
#include <irtkReconstructionb0.h>
#include <vector>
#include <string>
using namespace std;

//Application to perform reconstruction of volumetric MRI from thick slices.

void usage()
{
  cerr << "Usage: undistort [input] [output] [fieldmap] <options>\n" << endl;
  cerr << "Options:\n" << endl;
  cerr << "<-x>        phase encoding direction is x [default: y]" << endl;
  cerr << "<-y>        phase encoding direction is y [default: y]" << endl;
  cerr << "<-wfs wfs>  water fat shift in pixels [default: 16.895]" << endl;
  cerr << "<-res res>  acquired resolution in phase-encoding direction [default: 2.32mm]" << endl;
  cerr << "<-3T>       3T scan [default: 1.5mm]" << endl;
  cerr << "<-shim_fm sx sy sz> shim values for fieldmap [default: 0 0 0]" << endl;  
  cerr << "<-shim_image sx sy sz> shim values for image [default: 0 0 0]" << endl;  
  cerr << "<-sinc>     Use sinc interpolation [default: linear]" << endl;
  cerr << "<-debug>    Save intermediate results [default: off]" << endl;
  cerr << endl;

  exit(1);
}

int main(int argc, char **argv)
{
  
  //utility variables
  int ok;
  
  //declare variables for input
  /// Name for output volume
  char * output_name = NULL;
  //phase encode is y
  bool swap = true;
  //water fat shift
  double wfs = 16.895;
  //acquired resolution in phase encoding direction
  double res = 2.32;
  //is it 3T? default 1.5T
  bool is3T = false;
  //shim values for fieldmap
  double fx=0,fy=0,fz=0;
  //shim values for image
  double sx=0,sy=0,sz=0;
  //use sinc interpolation
  bool sinc = false;
  //save intermediate results
  bool debug = false;
  
  int i,j,k,t;
  double x,y,z;
  double factor,sfactor;


  if (argc < 3)
    usage();

  //read template
  irtkRealImage image(argv[1]);
  cout<<"Image ... "<<argv[1]<<endl;
  cout.flush();
  argc--;
  argv++;

  //read output name
  output_name = argv[1];
  argc--;
  argv++;
  cout<<"output name ... "<<output_name<<endl;
  cout.flush();

  //read template
  irtkRealImage fieldmap(argv[1]);
  cout<<"Image ... "<<argv[1]<<endl;
  cout.flush();
  argc--;
  argv++;

  // Parse options.
  while (argc > 1){
    ok = false;
    
    if ((ok == false) && (strcmp(argv[1], "-x") == 0)){
      argc--;
      argv++;
      swap=false;
      cout<< "Phase endoding direction is x."<<endl;
      cout.flush();
      ok = true;
    }

    if ((ok == false) && (strcmp(argv[1], "-y") == 0)){
      argc--;
      argv++;
      swap=false;
      cout<< "Phase endoding direction is y."<<endl;
      cout.flush();
      ok = true;
    }
    
    if ((ok == false) && (strcmp(argv[1], "-wfs") == 0)){
      argc--;
      argv++;
      wfs=atof(argv[1]);
      cout<< "WFS is "<<wfs<<endl;
      cout.flush();
      argc--;
      argv++;
      ok = true;
    }

    if ((ok == false) && (strcmp(argv[1], "-res") == 0)){
      argc--;
      argv++;
      res=atof(argv[1]);
      cout<< "Res is "<<res<<endl;
      cout.flush();
      argc--;
      argv++;
      ok = true;
    }

    if ((ok == false) && (strcmp(argv[1], "-3T") == 0)){
      argc--;
      argv++;
      is3T=true;
      cout<< "Scanner is 3T."<<endl;
      cout.flush();
      ok = true;
    }

    if ((ok == false) && (strcmp(argv[1], "-shim_fm") == 0)){
      argc--;
      argv++;
      fx=atof(argv[1]);
      argc--;
      argv++;
      fy=atof(argv[1]);
      argc--;
      argv++;
      fz=atof(argv[1]);
      argc--;
      argv++;
      cout<< "Shim for fieldmap is "<<fx<<" "<<fy<<" "<<fz<<" "<<endl;
      cout.flush();
      ok = true;
    }

    if ((ok == false) && (strcmp(argv[1], "-shim_image") == 0)){
      argc--;
      argv++;
      sx=atof(argv[1]);
      argc--;
      argv++;
      sy=atof(argv[1]);
      argc--;
      argv++;
      sz=atof(argv[1]);
      argc--;
      argv++;
      cout<< "Shim for image is "<<sx<<" "<<sy<<" "<<sz<<" "<<endl;
      cout.flush();
      ok = true;
    }

    if ((ok == false) && (strcmp(argv[1], "-debug") == 0)){
      argc--;
      argv++;
      debug=true;
      cout<< "Debug on."<<endl;
      cout.flush();
      ok = true;
    }


    if ((ok == false) && (strcmp(argv[1], "-sinc") == 0)){
      argc--;
      argv++;
      sinc=true;
      cout<< "Use sinc."<<endl;
      cout.flush();
      ok = true;
    }

    if (ok == false){
      cerr << "Can not parse argument " << argv[1] << endl;
      usage();
    }
  }

  if(debug)
    fieldmap.Write("fieldmap.nii.gz");
  
  //differece (in Hz) between resonance frequency of water and fat at 1.5T 
  //omega = 3.5 ppm x 1.5T x gama (gama=2.675e8 Hz/T)
  double omega = 223.52;
  
  //double omega for 3T
  if(is3T)
    omega*=2;
  
  //to convert fieldmap from Hz to distortion in mm
  //d=f*wfs*res/omega [mm]
  factor = wfs*res/omega;
  
  //correct shimming from fieldmap  
  sfactor = factor*267/(2*3.1415);
  //differential shim
  //sx-=fx;
  //sy-=fy;
  //sz-=fz;
  
  //resample fieldmap and add image shim
  irtkImageAttributes attr = image.GetImageAttributes();
  attr._t=1;
  irtkRealImage resfieldmap(attr);
  resfieldmap=0;
  irtkImageFunction *interpolatorLin = new irtkLinearInterpolateImageFunction;
  interpolatorLin->SetInput(&fieldmap);
  interpolatorLin->Initialize();

  for (k = 0; k < resfieldmap.GetZ(); k++) {
    for (j = 0; j < resfieldmap.GetY(); j++) {
      for (i = 0; i < resfieldmap.GetX(); i++) {
        x = i;
        y = j;
        z = k;
        resfieldmap.ImageToWorld(x, y, z);
	resfieldmap(i,j,k)= sfactor * ((fy-sy)*x-(fx-sx)*y+(fz-sz)*z);
	
        fieldmap.WorldToImage(x,y,z);
	if ((x > -0.5) && (x < fieldmap.GetX()-0.5) && 
	    (y > -0.5) && (y < fieldmap.GetY()-0.5) &&
            (z > -0.5) && (z < fieldmap.GetZ()-0.5))
	  {
	    resfieldmap(i,j,k)+= factor * interpolatorLin->Evaluate(x,y,z);
	  }
      }
    }
  }
  if(debug)
    resfieldmap.Write("fieldmap_resampled.nii.gz");
  
  cout<<fx-sx<<" "<<fy-sy<<" "<<fz-sz<<endl;
  cout<<factor<<endl;
  cout<<sfactor<<endl;
  
  irtkRealImage output = image;
  output=0;

  irtkImageFunction *interpolatorSinc;
  if(sinc)
    interpolatorSinc = new irtkSincInterpolateImageFunction;
  else
    interpolatorSinc = new irtkLinearInterpolateImageFunction;
  
  interpolatorSinc->SetInput(&image);
  interpolatorSinc->Initialize();
  
  attr = image.GetImageAttributes();
  for (t = 0; t < image.GetT(); t++) {
    for (k = 0; k < image.GetZ(); k++) {
      for (j = 0; j < image.GetY(); j++) {
        for (i = 0; i < image.GetX(); i++) {
          x = i;
          y = j;
          z = k;
	  //move it by fieldmap converted to voxels (reconstructed reslution)
	  y+=resfieldmap(i,j,k)/attr._dy;
	  if ((x > -0.5) && (x < image.GetX()-0.5) && 
	      (y > -0.5) && (y < image.GetY()-0.5) &&
              (z > -0.5) && (z < image.GetZ()-0.5))
	  {
	    output(i,j,k,t) = interpolatorSinc->Evaluate(x,y,z,t);
	  }
        }
      }
    }
  }

  output.Write(output_name);
  //The end of main()

  /*adjustment by jacobian does not seem to bing good result - perhaps need smoothing?
  irtkRealImage outputjac = image;
  outputjac=0;

  irtkImageFunction *interpolator = new irtkLinearInterpolateImageFunction;
  interpolator->SetInput(&resfieldmap);
  interpolator->Initialize();
  double lower, higher;
  int lj,hj;
  attr = image.GetImageAttributes();
  for (t = 0; t < image.GetT(); t++) {
    for (k = 0; k < image.GetZ(); k++) {
      for (j = 0; j < image.GetY(); j++) {
        for (i = 0; i < image.GetX(); i++) {
          x = i;
          y = j;
          z = k;
	  //move it by fieldmap converted to voxels (reconstructed reslution)
	  //y+=resfieldmap(i,j,k)/attr._dy;
	  if ((x > 0) && (x < image.GetX()-1) && 
	      (y > 0) && (y < image.GetY()-1) &&
              (z > 0) && (z < image.GetZ()-1))
	  {
	    lower = y-0.5 + interpolator->Evaluate(x,y-0.5,z)/attr._dy;
	    higher = y+0.5 +interpolator->Evaluate(x,y+0.5,z)/attr._dy;
	    //outputjac(i,j,k,t) = higher-lower;
	    //output(i,j,k,t)*=(higher-lower);
	    lj=round(lower);
	    hj=round(higher);
	    if(lj>hj)
	      output(i,j,k,t)=0;
	    else
	    {
	      if(lj==hj)
	        output(i,j,k,t)=image(i,hj,k,t)*(higher-lower);
	        //outputjac(i,j,k,t)=(higher-lower);
	      else
	      {
		output(i,j,k,t)=(lj+0.5-lower)*image(i,lj,k,t);
		
		for(int jj=(lj+1);jj<hj;jj++)
		  output(i,j,k,t)+=image(i,jj,k,t);
		output(i,j,k,t)+=(higher-hj+0.5)*image(i,hj,k,t);
	      }
	    }
	  }
        }
      }
    }
  }

  outputjac.Write("jacobian.nii.gz");
  output.Write("output_jac.nii.gz");
*/
}  
