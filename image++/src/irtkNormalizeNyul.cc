#include <irtkNormalizeNyul.h>

/* "This program implements intensity normalization algorithm published in "
       "Nyul, L.G.; Udupa, J.K.; Xuan Zhang, "
       "\"New variants of a method of MRI scale standardization,\""
       "Medical Imaging, IEEE Transactions on , vol.19, no.2, pp.143-150, Feb. 2000 "
       "http://dx.doi.org/10.1109/42.836373 "

    Source code adapted from an implementation by Vladimir Fonov in EZminc (https://github.com/vfonov/EZminc)

       */


irtkNormalizeNyul::irtkNormalizeNyul(irtkRealImage source, irtkRealImage target){
	_source = source;
	_target = target;
	_source_padding = 0;
	_target_padding = 0;
}


void irtkNormalizeNyul::SetMask(irtkRealImage source_mask, irtkRealImage target_mask){
	irtkRealPixel * iPtr, *mPtr;
	iPtr = _source.GetPointerToVoxels();
	mPtr = source_mask.GetPointerToVoxels();
	for(int i = 0; i < _source.GetNumberOfVoxels(); i++){
		if(*mPtr<1){
			*iPtr = _source_padding;
		}
		iPtr++;
		mPtr++;
	}
	iPtr = _target.GetPointerToVoxels();
	mPtr = source_mask.GetPointerToVoxels();
	for(int i = 0; i < _target.GetNumberOfVoxels(); i++){
		if(*mPtr<1){
			*iPtr = _target_padding;
		}
		iPtr++;
		mPtr++;
	}
}

irtkRealImage irtkNormalizeNyul::GetOutput(){
	return _source;
}

irtkRealImage irtkNormalizeNyul::GetTarget(){
	return _target;
}

void irtkNormalizeNyul::SetPadding(int source_padding, int target_padding){
	_source_padding = source_padding;
	_target_padding = target_padding;
}

void irtkNormalizeNyul::Run(){
	   int verbose=0;
//	  int normalize=0;
//	  int bimodalT=0;
//	  int debug=0;
	  int hist_bins=4000;
	  int steps=10;
//	  int clobber=0;
//	  int fix_zero_padding=0;

	  double cut_off=0.01;
	  float src_min,src_max;
	  float trg_min,trg_max;

	  irtkImageHistogram_1D<irtkRealPixel> src_hist_s, trg_hist_s;

	  src_hist_s.PutNumberOfBins(hist_bins);
	  trg_hist_s.PutNumberOfBins(hist_bins);

	  src_hist_s.Evaluate(&_source,_source_padding);
	  trg_hist_s.Evaluate(&_target,_target_padding);

	  src_min = src_hist_s.GetMin();
	  src_max = src_hist_s.GetMax();

	  trg_min = trg_hist_s.GetMin();
	  trg_max = trg_hist_s.GetMax();

	  double step  = (1.0 - 2.0*cut_off/100.0) / steps ;

	  std::vector<double> src_levels_s;
	  std::vector<double> trg_levels_s;

	  //provide mapping for background
	  src_levels_s.push_back(src_min);

	  trg_levels_s.push_back(trg_min);

	  if(verbose){
		std::cout<<"[ min ] "<<src_levels_s[0]<<" => "<<trg_levels_s[0]<<std::endl;
		std::cout << "nr tgt samples: " << trg_hist_s.NumberOfSamples() << ". nr src samples: " << src_hist_s.NumberOfSamples() << endl;
		std::cout << "nr tgt bins: " << trg_hist_s.NumberOfBins() << ". nr src bins: " << src_hist_s.NumberOfBins() << endl;
		for(int i = 0 ; i < src_hist_s.NumberOfBins(); i++){

		}
	  }



	  for(int i=0;i<=steps;i++)
	  {
		double pct =  i * step + cut_off/100;

		double src_lev_s=src_hist_s.CDFToVal(pct);
		double trg_lev_s=trg_hist_s.CDFToVal(pct);

		if(verbose)
			std::cout << "step " << i << " perc: " << pct << " src perc: " << src_lev_s << " trg perc: " << trg_lev_s << endl;

		if(trg_lev_s-trg_levels_s[trg_levels_s.size()-1]<(trg_max-trg_min)/100000.0)
		{
		  std::cerr<<"Warning: "<<pct*100<<" percentile collapses in target, skipping"<<std::endl;
		} else {

		  src_levels_s.push_back(src_lev_s);
		  trg_levels_s.push_back(trg_lev_s);

		  if(verbose)
			std::cout<<"[ "<<pct*100.0<<" ] "<<src_levels_s[src_levels_s.size()-1]<<" => "<<trg_levels_s[trg_levels_s.size()-1]<<std::endl;
		}
	  }
	  //provide mapping for upper range

	  src_levels_s.push_back(src_max);
	  trg_levels_s.push_back(trg_max);


	if(verbose)
		std::cout<<"[ max ] "<<src_levels_s[src_levels_s.size()-1]<<" => "<<trg_levels_s[trg_levels_s.size()-1]<<std::endl;

	if(verbose)
		std::cout<<"Recalculating intensities..."<<std::flush;

	irtkRealPixel * ptr = _source.GetPointerToVoxels();

	for(int i=0;i<_source.GetNumberOfVoxels();i++)
	{
	  //use LUT to map the intensities
	  unsigned int bin;
	  double input=*ptr;

	  double output=input;
	  for(bin=0;bin<src_levels_s.size();bin++)
		if(input <= src_levels_s[bin]) break;

	  if(input <= this->_source_padding)
		  output = this->_source_padding;
	  else if(bin==0) // first bin ?
		  output=trg_levels_s[0];
	  else if(bin>=(src_levels_s.size()-1))
		  output=trg_levels_s[trg_levels_s.size()-1];
	  else
		  output=(input-src_levels_s[bin-1])
		  /(src_levels_s[bin]-src_levels_s[bin-1])
		  *(trg_levels_s[bin]-trg_levels_s[bin-1])+trg_levels_s[bin-1];

//	  src.c_buf()[i]=output;
	  *ptr = output;
	  ptr++;
	}

	if(verbose)
	  std::cout<<"Done!"<<std::endl;
}

