/*=========================================================================

  Library   : Image Registration Toolkit (IRTK)
  Module    : $Id: irtkEMClassificationMultiComp.cc 2 2008-12-23 12:40:14Z dr $
  Copyright : Imperial College, Department of Computing
              Visual Information Processing (VIP), 2008 onwards
  Date      : $Date: 2008-12-23 12:40:14 +0000 (‰∫? 23 ÂçÅ‰∫åÊú?2008) $
  Version   : $Revision: 2 $
  Changes   : $Author: dr $

=========================================================================*/

#include <irtkEMClassification.h>

irtkEMClassificationMultiComp::irtkEMClassificationMultiComp() : irtkEMClassification()
{
  cerr<<"irtkEMClassificationMultiComp() must have atlas ";
  exit(1);
}

irtkEMClassificationMultiComp::~irtkEMClassificationMultiComp()
{
	delete []_number_of_components;
	delete []_ne;
	delete []_ns;
}

irtkEMClassificationMultiComp::irtkEMClassificationMultiComp(int noAtlas, irtkRealImage **atlas, int *noComp)
{
	int i,j,k;
	//caculate notissues
	int noTissues = 0;
	for (i=0;i<noAtlas;i++){
		noTissues += noComp[i];
	}

	_atlas.AddProbabilityMaps(noAtlas, atlas);

	//process atlas
	irtkRealImage **output = new irtkRealImage*[noTissues];
	for (i = 0; i < noAtlas; i++) {
		for (j = 0; j < noComp[i]; j++){
			output[i+j] = new irtkRealImage(*atlas[i]);
		}
	}
	_output.AddProbabilityMaps(noTissues, output);
	_atlas.NormalizeAtlas();
	_output.NormalizeAtlas();

	_padding = MIN_GREY;
	_number_of_atlas = noAtlas;
	_number_of_tissues = noTissues;
	_number_of_components = noComp;
	_number_of_voxels = 0;

	_number_of_components = new int[_number_of_atlas];
	_mi = new double[_number_of_tissues];
	_sigma = new double[_number_of_tissues];
	_c = new double[_number_of_tissues];
	_ne = new int[_number_of_atlas];
	_ns = new int[_number_of_atlas];
	for (j = 0; j < _number_of_atlas; j++){
		_number_of_components[j] = noComp[j];
	}
	for (j = 0; j < _number_of_atlas; j++) {
		_ns[j] = 0; _ne[j] = 0;
		for (k = 0; k < j+1; k++){
			_ne[j] += _number_of_components[k];
		}
		for (k = 0; k < j; k++){
			_ns[j] += _number_of_components[k];
		}
	}

	_f = 0;
	_G = NULL;
	_background_tissue = -1;
	_debug = false;
	delete []output;
}

irtkEMClassificationMultiComp::irtkEMClassificationMultiComp(int noAtlas, irtkRealImage **atlas, irtkRealImage *background, int *noComp)
{
	int i,j,k,newed;
	newed = 0;
	if(noComp == NULL){
		noComp = new int[noAtlas + 1];
		for (i = 0; i < noAtlas + 1; i++)
			noComp[i] = 1;
		newed = 1;
	}
	//caculate notissues
	int noTissues = 0;
	for (i=0;i<noAtlas+1;i++){
		noTissues += noComp[i];
	}

	cout << "Number of atlas: "<<noAtlas<<endl;
	cout << "Number of components: "<<noTissues<<endl;

	_atlas.AddProbabilityMaps(noAtlas, atlas);

	//process atlas
	irtkRealImage **output = new irtkRealImage*[noTissues];
	k = 0;
	for (i = 0; i < noAtlas; i++) {
		for (j = 0; j < noComp[i]; j++){
			output[k+j] = new irtkRealImage(*atlas[i]);
		}
		k += noComp[i];
	}

	if (background == NULL) {
		cerr<<"adding background ..."<<endl;
		_atlas.AddBackground();
		for (j = 0; j < noComp[noAtlas]; j++){
			output[k+j] = new irtkRealImage(_atlas.GetImage(noAtlas));
		}
		_output.AddProbabilityMaps(noTissues, output);
		cerr<<"done"<<endl;
	} else {
		cerr<<"normalize ...";
		_atlas.NormalizeAtlas(*background);
		for (j = 0; j < noComp[noAtlas]; j++){
			output[k+j] = new irtkRealImage(_atlas.GetImage(noAtlas));
		}
		_output.AddProbabilityMaps(noTissues, output);
		cerr<<"done ..."<<endl;
	}
	_background_tissue = noTissues;
	_padding = MIN_GREY;
	_number_of_atlas = noAtlas+1;
	_number_of_tissues = noTissues;

	_number_of_components = new int[_number_of_atlas];
	_mi = new double[_number_of_tissues];
	_sigma = new double[_number_of_tissues];
	_c = new double[_number_of_tissues];
	_ne = new int[_number_of_atlas];
	_ns = new int[_number_of_atlas];
	for (j = 0; j < _number_of_atlas; j++){
		if(newed == 1){
			_number_of_components[j] = 1;
		}else{
			_number_of_components[j] = noComp[j];
		}
	}
	for (j = 0; j < _number_of_atlas; j++) {
		_ns[j] = 0; _ne[j] = 0;
		for (k = 0; k < j+1; k++){
			_ne[j] += _number_of_components[k];
		}
		for (k = 0; k < j; k++){
			_ns[j] += _number_of_components[k];
		}
	}


	_f = 0;
	_number_of_voxels = 0;
	_G = NULL;
	_debug = false;
	delete []output;
	if(newed == 1){
		delete []noComp;
	}
}

void irtkEMClassificationMultiComp::Initialise()
{
	int i,j,n;
	
	//Evaluate Guassian filters
	this->MStep();

	// A easy approach to split
	for (i = 0; i < _number_of_atlas; i++) {
		n = _number_of_components[i];
		for (j = _ns[i]; j < _ne[i]; j++){
			_mi[j] = _mi[j] + double(2*(j-_ns[i])+1-n)*sqrt(_sigma[j])/double(n);
			_sigma[j] = _sigma[j]/n;
		}
	}
	// Update approach to split
	/*int *current_no,split_no,to_no,max_sigma;
	
	current_no = new int[_number_of_atlas];
	//current number of components
	for (i = 0; i < _number_of_atlas; i++){
		current_no[i] = 1;
	}
	
	//for
	for (i = 0; i < _number_of_atlas; i++){
		//check if current atlas needs split
		n = _number_of_components[i];
		while (current_no[i] < n){
			//find out which one to split
			max_sigma = 0;
			for (j = _ns[i]; j < _ns[i] + current_no[i]; j++){
				if(max_sigma < _sigma[j]){
					max_sigma = _sigma[j];
					split_no = j;
				}
			}
			//find out which one to split to
			to_no = _ns[i] + current_no[i];
			//split by sigma
			_mi[split_no] = _mi[split_no] - sqrt(_sigma[split_no])/2;
			_sigma[split_no] = _sigma[split_no]/2;
			_mi[to_no] = _mi[to_no] + sqrt(_sigma[split_no])/2;
			_sigma[to_no] = _sigma[split_no];
			//estimate new gmm based on atlas > 1
			this->InitialEStep(i,current_no[i]+1);
			this->InitialMStep(i,current_no[i]+1);
			//add current number
			current_no[i] ++;
		}
	}//end
	delete []current_no;
	*/

	//Update output
	this->EStep();
	Print();	
}

void irtkEMClassificationMultiComp::MStep()
{
  int i, j, k;
  double fraction;
  vector<double> mi_num(_number_of_tissues);
  vector<double> sigma_num(_number_of_tissues);
  vector<double> denom(_number_of_tissues);

  for (k = 0; k < _number_of_tissues; k++) {
    mi_num[k] = 0;
  }

  for (k = 0; k < _number_of_tissues; k++) {
    sigma_num[k] = 0;
  }

  for (k = 0; k < _number_of_tissues; k++) {
    denom[k] = 0;
  }

  _output.First();
  irtkRealPixel *ptr = _input.GetPointerToVoxels();
  irtkRealPixel *pm = _mask.GetPointerToVoxels();

  for (i = 0; i < _input.GetNumberOfVoxels(); i++) {
    // Check for backgound
    if (*pm == 1) {
      for (k = 0; k < _number_of_tissues; k++) {
        mi_num[k] += _output.GetValue(k) * *ptr;
        denom[k]  += _output.GetValue(k);
      }
    }
    ptr++;
    pm++;
    _output.Next();
  }

  for (k = 0; k < _number_of_tissues; k++) {
    if (denom[k] != 0) {
      _mi[k] = mi_num[k] / denom[k];
    } else {
	  _mi[k] = -1;
      //cerr << "Division by zero while computing tissue mean!" << endl;
      //exit(1);
    }
  }

  for (j = 0; j < _number_of_atlas; j++) {
	  fraction = 0;
	  for (k = _ns[j]; k < _ne[j]; k++){
		  fraction += denom[k];
	  }
	  for (k = _ns[j]; k < _ne[j]; k++){
		  if(_mi[k]!=-1)
			  _c[k] = denom[k]/fraction;
		  else
			  _c[k] = 0;
	  }
  }

  _output.First();
  ptr = _input.GetPointerToVoxels();
  pm = _mask.GetPointerToVoxels();

  for (i = 0; i < _input.GetNumberOfVoxels(); i++) {
    // Check for backgound
    if (*pm == 1) {
      for (k = 0; k <_number_of_tissues; k++) {
		  if(_mi[k]!=-1)
			 sigma_num[k] += (_output.GetValue(k) * (*ptr - _mi[k]) * (*ptr - _mi[k]));
      }
    }
    ptr++;
    pm++;
    _output.Next();
  }

  for (k = 0; k <_number_of_tissues; k++) {
	  if(_mi[k]!=-1)
		_sigma[k] = sigma_num[k] / denom[k];
	  else
		_sigma[k] = -1;
  }
}

void irtkEMClassificationMultiComp::EStep()
{
  int i, j, k, x, y, z, t;
  double v,distance;

  irtkGaussian* G = new irtkGaussian[_number_of_tissues];
  double *center_x = new double[_number_of_tissues];
  double *center_z = new double[_number_of_tissues];
  double *center_y = new double[_number_of_tissues];
  double *center_weight = new double[_number_of_tissues];

  for (k = 0; k < _number_of_tissues; k++) {
    G[k].Initialise( _mi[k], _sigma[k]);
	center_x[k] = 0; center_y[k] = 0; 
	center_z[k] = 0; center_weight[k] = 0;
  }

  _atlas.First();
  _output.First();
  irtkRealPixel *ptr = _input.GetPointerToVoxels();
  irtkRealPixel *pm = _mask.GetPointerToVoxels();
  for (i=0; i< _input.GetNumberOfVoxels(); i++) {
	  double* gv = new double[_number_of_tissues];
	  double* numerator = new double[_number_of_tissues];
	  if (*pm == 1) {
		  v = *ptr;
		  for (j = 0; j < _number_of_atlas; j++){
			  for (k = _ns[j]; k < _ne[j]; k++) {
				  if(_mi[k]!=-1){
					  gv[k] = G[k].Evaluate(v);
					  numerator[k] = gv[k] * _atlas.GetValue(j);
					  center_weight[k] += numerator[k];
					  _input.GetImageAttributes().IndexToLattice(i,&x,&y,&z,&t);
					  center_x[k] += x*numerator[k]; 
					  center_y[k] += y*numerator[k]; 
					  center_z[k] += z*numerator[k];
				  }
			  }
		  }
	  }
	  ptr++;
	  pm++;
	  _atlas.Next();
	  _output.Next();
	  delete[] gv;
	  delete[] numerator;
  }
  for (k = 0; k < _number_of_tissues; k++) {
	  if(center_weight[k] != 0){
		  center_x[k] = center_x[k]/center_weight[k];
		  center_y[k] = center_y[k]/center_weight[k];
		  center_z[k] = center_z[k]/center_weight[k];
	  }
  }

  _atlas.First();
  _output.First();
  ptr = _input.GetPointerToVoxels();
  pm = _mask.GetPointerToVoxels();
  for (i=0; i< _input.GetNumberOfVoxels(); i++) {
    double* gv = new double[_number_of_tissues];
    double* numerator = new double[_number_of_tissues];
    double denominator=0;
    if (*pm == 1) {
      v = *ptr;
	  for (j = 0; j < _number_of_atlas; j++){
		  for (k = _ns[j]; k < _ne[j]; k++) {
			  if(_mi[k]!=-1){
				  gv[k] = G[k].Evaluate(v);
				  _input.GetImageAttributes().IndexToLattice(i,&x,&y,&z,&t);
				  distance = sqrt(pow(x-center_x[k],2)+pow(y-center_y[k],2)+pow(z-center_z[k],2)) + 1;
				  numerator[k] = gv[k] * _atlas.GetValue(j) * _c[k] / distance;
				  denominator += numerator[k];
			  }else{
				  numerator[k] = 0;
			  }
		  }
	  }
      for (k = 0; k < _number_of_tissues; k++) {
        if (denominator != 0) {
          double value = numerator[k]/denominator;
          _output.SetValue(k, value);
          if ((value < 0) || (value > 1)) {
            cerr << "Probability value out of range = " << value << endl;
            cerr << value << " " << k << " " << i << " " << _atlas.GetValue(k) <<  " " << _sigma[k] << endl;
            exit(1);
          }
        } else {
          cerr << "Division by 0 while computing probabilities" << endl;
          cerr<<"tissue="<<k;
          if (k==0) _output.SetValue(k, 1);
          else     _output.SetValue(k, 0);


          //exit(1);
        }
      }
    } else {
		for (j = 0; j < _number_of_atlas - 1; j++) {
			for (k = _ns[j]; k < _ne[j]; k++) {
				_output.SetValue(k, 0);
			}
		}
		for (k = _ns[_number_of_atlas - 1]; k < _ne[_number_of_atlas - 1]; k++) {
			_output.SetValue(k, 1);
		}
    }
    ptr++;
    pm++;
    _atlas.Next();
    _output.Next();
    delete[] gv;
    delete[] numerator;
  }
  delete[] G;
  delete[] center_x;
  delete[] center_y;
  delete[] center_z;
  delete[] center_weight;

}

void irtkEMClassificationMultiComp::InitialMStep(int j, int n)
{
  int i, k;
  double fraction;
  vector<double> mi_num(n);
  vector<double> sigma_num(n);
  vector<double> denom(n);

  for (k = 0; k < n; k++) {
    mi_num[k] = 0;
  }

  for (k = 0; k < n; k++) {
    sigma_num[k] = 0;
  }

  for (k = 0; k < n; k++) {
    denom[k] = 0;
  }

  _output.First();
  _atlas.First();
  irtkRealPixel *ptr = _input.GetPointerToVoxels();
  irtkRealPixel *pm = _mask.GetPointerToVoxels();

  for (i = 0; i < _input.GetNumberOfVoxels(); i++) {
    // Check for backgound
    if (*pm == 1 && _atlas.GetValue(j) > 0) {
      for (k = 0; k < n; k++) {
        mi_num[k] += _output.GetValue(k+_ns[j]) * *ptr;
        denom[k]  += _output.GetValue(k+_ns[j]);
      }
    }
    ptr++;
    pm++;
    _output.Next();
	_atlas.Next();
  }

  for (k = _ns[j]; k < _ns[j] + n; k++) {
    if (denom[k - _ns[j]] != 0) {
      _mi[k] = mi_num[k - _ns[j]] / denom[k - _ns[j]];
    } else {
      cerr << "Division by zero while computing tissue mean!" << endl;
      exit(1);
    }
  }

  fraction = 0;
  for (k = _ns[j]; k < _ns[j] + n; k++){
	  fraction += denom[k-_ns[j]];
  }
  for (k = _ns[j]; k < _ns[j] + n; k++){
	  _c[k] = denom[k-_ns[j]]/fraction;
  }
  
  _output.First();
  _atlas.First();
  ptr = _input.GetPointerToVoxels();
  pm = _mask.GetPointerToVoxels();

  for (i = 0; i < _input.GetNumberOfVoxels(); i++) {
    // Check for backgound
    if (*pm == 1 && _atlas.GetValue(j) > 0) {
      for (k = _ns[j]; k < _ns[j] + n; k++) {
        sigma_num[k - _ns[j]] += (_output.GetValue(k) * (*ptr - _mi[k]) * (*ptr - _mi[k]));
      }
    }
    ptr++;
    pm++;
    _output.Next();
	_atlas.Next();
  }

  for (k = _ns[j]; k < _ns[j] + n; k++) {
    _sigma[k] = sigma_num[k-_ns[j]] / denom[k-_ns[j]];
  }
}

void irtkEMClassificationMultiComp::InitialEStep(int j, int n)
{
  int i, k;
  double x;

  irtkGaussian* G = new irtkGaussian[n];

  for (k = 0; k < n; k++) {
    G[k].Initialise( _mi[_ns[j]+k], _sigma[_ns[j]+k]);
  }

  _atlas.First();
  _output.First();
  irtkRealPixel *ptr = _input.GetPointerToVoxels();
  irtkRealPixel *pm = _mask.GetPointerToVoxels();
  for (i=0; i< _input.GetNumberOfVoxels(); i++) {
    double* gv = new double[n];
    double* numerator = new double[n];
    double denominator = 0;
    if (*pm == 1 && _atlas.GetValue(j) > 0) {
      x = *ptr;
	  for (k = 0; k < n; k++) {
		  gv[k] = G[k].Evaluate(x);
		  numerator[k] = gv[k] * _c[k + _ns[j]];
		  denominator += numerator[k];
	  }
      for (k = _ns[j]; k < _ns[j]+n; k++) {
        if (denominator != 0) {
          double value = numerator[k - _ns[j]]/denominator;
          _output.SetValue(k, value);
          if ((value < 0) || (value > 1)) {
            cerr << "Probability value out of range = " << value << endl;
            cerr << value << " " << k << " " << i << " " << _atlas.GetValue(k) <<  " " << _sigma[k] << endl;
            exit(1);
          }
        } else {
          cerr << "Division by 0 while computing probabilities" << endl;
          cerr<<"tissue="<<k;
          if (k==0) _output.SetValue(k, 1);
          else     _output.SetValue(k, 0);


          //exit(1);
        }
      }
    } else {
      for (k = _ns[j]; k < _ns[j]+n; k++) {
        _output.SetValue(k, 0);
      }
    }
    ptr++;
    pm++;
    _atlas.Next();
    _output.Next();
    delete[] gv;
    delete[] numerator;
  }
  delete[] G;

}

double irtkEMClassificationMultiComp::Iterate(int)
{
  this->MStep();
  this->EStep();
  Print();
  PrintGMM();
  return LogLikelihood();
}

double irtkEMClassificationMultiComp::LogLikelihood()
{
  int i, k;
  double temp, f;
  cerr<< "Log likelihood: ";
  irtkGaussian* G = new irtkGaussian[_number_of_tissues];
  double* gv = new double[_number_of_tissues];

  for (k = 0; k < _number_of_tissues; k++) {
    G[k].Initialise( _mi[k], _sigma[k]);
  }

  irtkRealPixel *ptr = _input.GetPointerToVoxels();
  irtkRealPixel *pm = _mask.GetPointerToVoxels();
  _output.First();
  f = 0;
  for (i = 0; i < _input.GetNumberOfVoxels(); i++) {
    if (*pm == 1) {
      temp = 0;
      for (k=0; k < _number_of_tissues; k++) {
        // Estimation of gaussian probability of intensity (*ptr) for tissue k
        gv[k] = G[k].Evaluate(*ptr);
        // Probability that current voxel is from tissue k
        temp += gv[k] *	_output.GetValue(k);
      }

      if ((temp > 1) || (temp < 0)) {
        cerr << "Could not compute likelihood, probability out of range = " << temp << endl;
        exit(1);
      }
      f += log(temp);
    }
    ptr++;
    pm++;
    _output.Next();
  }

  f = -f;
  double diff, rel_diff;
  diff = _f-f;

  if (_f == 0) rel_diff = 1;
  else rel_diff = diff/_f;

  _f=f;

  cerr << "f= "<< f << " diff = " << diff << " rel_diff = " << rel_diff <<endl;
  delete[] G;
  delete[] gv;

  return rel_diff;
}

void irtkEMClassificationMultiComp::ConstructSegmentation(irtkRealImage &segmentation)
{
  int i, j, k;
  double max,value;
  irtkRealPixel *ptr;

  cerr<<"Constructing segmentation"<<endl;

  // Initialize pointers of probability maps
  _output.First();

  // Initialize segmentation to same size as input
  segmentation = _input;
  ptr = segmentation.GetPointerToVoxels();

  for (i = 0; i < _input.GetNumberOfVoxels(); i++) {

    if (*ptr != _padding) {
      max  = 0;
      *ptr = 0;
      for (j = 0; j < _number_of_atlas; j++) {
		  value = 0;
		  for (k = _ns[j]; k < _ne[j]; k++){
			  value += _output.GetValue(k);
		  }
		  if (value > max) {
			  max  = value;
			  *ptr = j+1;
			  if ( (j+1) == _number_of_atlas) *ptr=0;
		  }
      }
    }else *ptr = 0;
    ptr++;
    _output.Next();
  }
}
void irtkEMClassificationMultiComp::ConstructSegmentationWithPadding(irtkRealImage &segmentation)
{
  int i, j, k;
  double max,value;
  irtkRealPixel *ptr;

  cerr<<"Constructing segmentation with Padding"<<endl;
   // Initialize pointers of probability maps
  _output.First();

  // Initialize segmentation to same size as input
  segmentation = _input;
  ptr = segmentation.GetPointerToVoxels();

  for (i = 0; i < _input.GetNumberOfVoxels(); i++) {

    if (*ptr != _padding) {
      max  = 0;
      *ptr = 0;
	  for (j = 0; j < _number_of_atlas; j++) {
		  value = 0;
		  for (k = _ns[j]; k < _ne[j]; k++){
			  value += _output.GetValue(k);
		  }
		  if (value > max) {
			  max  = value;
			  *ptr = j+1;
		  }
      }
    } else *ptr=_padding;
    ptr++;
    _output.Next();
  }
}
void irtkEMClassificationMultiComp::GetProbMap(int i,irtkRealImage& image){
	int j,k;
	irtkRealPixel *ptr;

	if  (i < _number_of_atlas) {
		_output.First();
		// Initialize segmentation to same size as input
		image = _input;

		ptr = image.GetPointerToVoxels();

		for (j = 0; j < _input.GetNumberOfVoxels(); j++) {
			*ptr = 0;
			for (k = _ns[i]; k < _ne[i]; k++) {
				*ptr += _output.GetValue(k);
			}
			ptr++;
			_output.Next();
		}
	} else {
		cerr << "irtkProbabilisticAtlas::Write: No such probability map" << endl;
		exit(1);
	}
}
void irtkEMClassificationMultiComp::ConstructSegmentation()
{
  int i, j, k;
  double max,value;
  irtkRealPixel *ptr;

  cerr<<"Constructing segmentation with Padding"<<endl;
   // Initialize pointers of probability maps
  _output.First();

  // Initialize segmentation to same size as input
  _segmentation = _input;
  ptr = _segmentation.GetPointerToVoxels();

  for (i = 0; i < _input.GetNumberOfVoxels(); i++) {
    if (*ptr != _padding) {
      max  = 0;
      *ptr = 0;
      for (j = 0; j < _number_of_atlas; j++) {
		  value = 0;
		  for (k = _ns[j]; k < _ne[j]; k++){
			  value += _output.GetValue(k);
		  }
		  if (value > max) {
			  max  = value;
			  *ptr = j+1;
		  }
      }
    }else *ptr = 0;
    ptr++;
    _output.Next();
  }
}
