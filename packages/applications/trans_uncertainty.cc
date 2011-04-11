/*=========================================================================

  Library   : target Registration Toolkit (IRTK)
  Module    : $Id: hard_target.cc 2 2008-12-23 12:40:14Z dr $
  Copyright : Imperial College, Department of Computing
              Visual Information Processing (VIP), 2008 onwards
  Date      : $Date: 2008-12-23 12:40:14 +0000 (二, 23 十二月 2008) $
  Version   : $Revision: 2 $
  Changes   : $Author: dr $

=========================================================================*/

#include <irtkRegistration.h>

#include <irtkRician.h>


char *output_name,*parameter_name;

void usage()
{
  cerr << "Usage: trans_uncertainty [target] [transformed_source] [transformation_jacobian] [output]" << endl;
  exit(1);
}

int main(int argc, char **argv)
{
  int i, j, k;
  float meant, sigmat, pt, pti, meanti, sigmati;
  irtkRician Rician;
  irtkGaussian gaussian;
  //irtkMutualInformationSimilarityMetric metric;

  if (argc < 5) {
    usage();
    exit(1);
  }

  irtkGreyImage target;
  target.Read(argv[1]);
  argc--;
  argv++;

  irtkRealImage output;
  output.Initialize(target.GetImageAttributes());

  irtkGreyImage difference;
  difference.Initialize(target.GetImageAttributes());

  // source
  irtkGreyImage source;
  source.Read(argv[1]);
  argc--;
  argv++;

  // jacobian
  irtkRealImage jacobian;
  jacobian.Read(argv[1]);
  argc--;
  argv++;

  // File name for output
  output_name = argv[1];
  argc--;
  argv++;

  //Equalize images
  double min,max;
  target.GetMinMaxAsDouble(&min,&max);
  irtkImageHistogram_1D<irtkGreyPixel> histogram;
  histogram.Evaluate(&target,-1);
  histogram.Equalize(min,max);
  histogram.BackProject(&target);
  histogram.Evaluate(&source,-1);
  histogram.Equalize(min,max);
  histogram.BackProject(&source);

  //calculate difference
  meant = target.GetAverage();
  meanti = source.GetAverage();
  for (i=0;i<target.GetX();i++){
	  for (j=0;j<target.GetY();j++){
		  for (k=0;k<target.GetZ();k++){
			  if(source.GetAsDouble(i,j,k)>0){
			  difference.PutAsDouble(i,j,k,
				  target.GetAsDouble(i,j,k)
				  - source.GetAsDouble(i,j,k));
			  }else{
			  difference.PutAsDouble(i,j,k,meant - meanti);
			  }
		  }
	  }
  }

  difference.Write("uncertainty\\difference.nii.gz");

  //Evaluate uncertainty
  meant = jacobian.GetAverage(0);
  sigmat = jacobian.GetSD(0);
  sigmat = sigmat*sigmat;

  Rician.Initialise(meant,sigmat);
  Rician.Approximate();

  meanti = difference.GetAverage(0);
  sigmati = difference.GetSD(0);
  sigmati = sigmati*sigmati;

  gaussian.Initialise(meanti,sigmati);

  meant = Rician.Evaluate(meant);
  meanti = gaussian.Evaluate(meanti);

  for (i=0;i<target.GetX();i++){
	  for (j=0;j<target.GetY();j++){
		  for (k=0;k<target.GetZ();k++){
			  pt = jacobian.GetAsDouble(i,j,k);
			  pt = Rician.Evaluate(pt)/meant;
			  pti = difference.GetAsDouble(i,j,k);
			  pti = gaussian.Evaluate(pti)/meanti;
			  output.PutAsDouble(i,j,k,
				  pt*pti);
		  }
	  }
  }

  //Output uncertainty
  output.Write(output_name);

}

