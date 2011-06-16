/*=========================================================================

  Library   : Image Registration Toolkit (IRTK)
  Module    : $Id$
  Copyright : Imperial College, Department of Computing
              Visual Information Processing (VIP), 2008 onwards
  Date      : $Date$
  Version   : $Revision$
  Changes   : $Author$

=========================================================================*/

#include <irtkEMClassification.h>

irtkEMClassification::irtkEMClassification()
{
  _padding = MIN_GREY;
  _mi = NULL;
  _sigma = NULL;
  _c = NULL;
  _f = 0;
  _G = NULL;
  _number_of_voxels = 0;
  _background_tissue = -1;
  _debug = false;

}

irtkEMClassification::irtkEMClassification(int noTissues, irtkRealImage **atlas, irtkRealImage *background)
{
  _atlas.AddProbabilityMaps(noTissues, atlas);
  _output.AddProbabilityMaps(noTissues, atlas);
  if (background == NULL) {
    cerr<<"adding background ..."<<endl;
    _atlas.AddBackground();
    _output.AddBackground();
    cerr<<"done"<<endl;
  } else {
    cerr<<"normalize ...";
    _atlas.NormalizeAtlas(*background);
    _output.NormalizeAtlas(*background);
    cerr<<"done ..."<<endl;
  }
  _background_tissue = noTissues;
  _padding = MIN_GREY;
  _number_of_tissues = noTissues+1;
  _mi = new double[_number_of_tissues+1];
  _sigma = new double[_number_of_tissues+1];
  _c = new double[_number_of_tissues+1];
  _f = 0;
  _number_of_voxels = 0;
  _G = NULL;
  _debug = false;
}

irtkEMClassification::irtkEMClassification(int noTissues, irtkRealImage **atlas)
{
  _atlas.AddProbabilityMaps(noTissues, atlas);
  _output.AddProbabilityMaps(noTissues, atlas);
  _atlas.NormalizeAtlas();
  _output.NormalizeAtlas();

  _padding = MIN_GREY;
  _number_of_tissues = noTissues;
  _number_of_voxels = 0;
  _mi = new double[_number_of_tissues];
  _sigma = new double[_number_of_tissues];
  _c = new double[_number_of_tissues];
  _f = 0;
  _G = NULL;
  _background_tissue = -1;
  _debug = false;
}

irtkEMClassification::~irtkEMClassification()
{
  if (_G!=NULL) delete []_G;
  delete []_mi;
  delete []_sigma;
  delete []_c;
}

void irtkEMClassification::SetInput(const irtkRealImage &image)
{
  _input = image;
  _estimate = image;
  _weights = image;
  _weightsB = image;
  _weightsR = image;
  _number_of_voxels=_input.GetNumberOfVoxels();
  CreateMask();
}

void irtkEMClassification::CreateMask()
{
  _mask=_input;
  irtkRealPixel *p=_mask.GetPointerToVoxels();
  for (int i=0; i<_mask.GetNumberOfVoxels(); i++) {
    if (*p!=_padding) *p=1;
    else *p=0;
    p++;
  }
}

void irtkEMClassification::SetMask(irtkRealImage &mask)
{
  _mask=mask;
}

void irtkEMClassification::Initialise()
{
  this->MStep();
  Print();
}

void irtkEMClassification::InitialiseGMM()
{
  this->MStepGMM();
  Print();
}


void irtkEMClassification::InitialisePosteriors()
{
  int i;

  if (_number_of_tissues == 0) {
    cerr<<"Can't initialize atlas, no tissues given"<<endl;
    exit(1);
  }

  if (_number_of_voxels == 0) {
    cerr<<"Can't initialize atlas, no input given"<<endl;
    exit(1);
  }

  if (_output.GetNumberOfTissues() == 0) {
    cerr<<"Initializing atlas ...";
    for (i = 0; i < _number_of_tissues; i++) _output.AddImage(_input);
    cerr<<" done."<<endl;
  } else {
    if (_output.GetNumberOfTissues()!=_number_of_tissues) {
      cerr<<"Error: Number of tissues mismatch:"<<_output.GetNumberOfTissues()<<endl;
      exit(1);
    }
  }
  EStepGMM();

}

void irtkEMClassification::InitialisePVSegmentation()
{
  int i;

  if (_number_of_tissues == 0) {
    cerr<<"Can't initialize atlas, no tissues given"<<endl;
    exit(1);
  }

  if (_number_of_voxels == 0) {
    cerr<<"Can't initialize atlas, no input given"<<endl;
    exit(1);
  }

  if (_pv_output.GetNumberOfTissues() == 0) {
    cerr<<"Initializing PV ...";
    for (i = 0; i < _number_of_tissues; i++) _pv_output.AddImage(_input);
    cerr<<" done."<<endl;
  } else {
/*
    if (_output.GetNumberOfTissues()!=_number_of_tissues) {
      cerr<<"Error: Number of tissues mismatch:"<<_output.GetNumberOfTissues()<<endl;
      exit(1);
    }
*/   cerr<<"Warning: Attempt to initialise PV segmentation. Current number of tissues: "<<_pv_output.GetNumberOfTissues()<<endl;
  }

}

void irtkEMClassification::InitialiseAtlas()
{
  int i;

  if (_number_of_tissues == 0) {
    cerr<<"Can't initialize atlas, no tissues given"<<endl;
    exit(1);
  }

  if (_number_of_voxels == 0) {
    cerr<<"Can't initialize atlas, no input given"<<endl;
    exit(1);
  }

  if (_atlas.GetNumberOfTissues() == 0) {
    cerr<<"Initializing atlas ...";
    for (i = 0; i < _number_of_tissues; i++) _atlas.AddImage(_input);
    cerr<<" done."<<endl;
  } else {
    if (_atlas.GetNumberOfTissues()!=_number_of_tissues) {
      cerr<<"Error: Number of tissues mismatch:"<<_output.GetNumberOfTissues()<<endl;
      exit(1);
    }
  }
  EStepGMM();

}

void irtkEMClassification::UniformPrior()
{
  int i, k;
  
  _atlas.First();
  irtkRealPixel *pm = _mask.GetPointerToVoxels();
  for (i=0; i< _input.GetNumberOfVoxels(); i++) {
    if (*pm == 1) {
      for (k = 0; k < _number_of_tissues; k++) {
        _atlas.SetValue(k,1.0/_number_of_tissues);
        }
      }
      pm++;
      _atlas.Next();
    }
}

void irtkEMClassification::InitialiseGMMParameters3()
{
  int i;

  cerr <<"Estimating GMM parameters ... ";
  irtkHistogram h(100);
  irtkRealPixel imin, imax;
  _input.GetMinMaxPad(&imin, &imax,_padding);
  cerr<<" min = "<<imin<<", max = "<<imax<<" ... ";
  h.PutMin(imin);
  h.PutMax(imax);
  irtkRealPixel *ptr=_input.GetPointerToVoxels();
  for (i=0; i<_input.GetNumberOfVoxels(); i++) {
    if (*ptr!=_padding) h.AddSample(*ptr);
    ptr++;
  }
  double mean, variance;
  mean = h.Mean();
  variance=h.Variance();
  cerr<<"mean="<<mean<<" variance="<<sqrt(variance)<<" ... done."<<endl;

  _number_of_tissues=3;
  _mi=new double[3];
  _sigma=new double[3];
  _c=new double[3];
  _mi[1]=mean;
  _mi[0]=mean-sqrt(variance)/2;
  _mi[2]=mean+sqrt(variance)/2;
  _sigma[0]=variance/4;
  _sigma[1]=variance/4;
  _sigma[2]=variance/4;
  _c[0]=0.10;
  _c[1]=0.65;
  _c[2]=0.25;
}

void irtkEMClassification::InitialiseGMMParameters(int n)
{
  int i;

  cerr <<"Estimating GMM parameters ... ";
  irtkRealPixel imin, imax;
  _input.GetMinMaxPad(&imin, &imax,_padding);
  cerr<<" min = "<<imin<<", max = "<<imax<<" ... ";
  _number_of_tissues=n;
  _mi=new double[n];
  _sigma=new double[n];
  _c=new double[n];
  for (i=0; i<n; i++) {
    _mi[i]=imin+i*(imax-imin)/(n-1);
    _sigma[i]=((imax-imin)/(n-1))*((imax-imin)/(n-1));
    _c[i]=1.0/n;
  }
  PrintGMM();
}


void irtkEMClassification::InitialiseGMMParameters(int n, double *m, double *s, double *c)
{
  int i;

  _number_of_tissues = n;
  _mi = new double[_number_of_tissues];
  _sigma = new double[_number_of_tissues];
  _c = new double[_number_of_tissues];

  for ( i=0; i<n; i++) {
    _mi[i]=m[i];
    _sigma[i]=s[i];
    _c[i]=c[i];
  }

  PrintGMM();
  InitialisePosteriors();

}

void irtkEMClassification::MStep()
{
  int i, k;
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
      cerr << "Division by zero while computing tissue mean!" << endl;
      exit(1);
    }
  }

  _output.First();
  ptr = _input.GetPointerToVoxels();
  pm = _mask.GetPointerToVoxels();

  for (i = 0; i < _input.GetNumberOfVoxels(); i++) {
    // Check for backgound
    if (*pm == 1) {
      for (k = 0; k <_number_of_tissues; k++) {
        sigma_num[k] += (_output.GetValue(k) * (*ptr - _mi[k]) * (*ptr - _mi[k]));
      }
    }
    ptr++;
    pm++;
    _output.Next();
  }

  for (k = 0; k <_number_of_tissues; k++) {
    _sigma[k] = sigma_num[k] / denom[k];
  }
}

void irtkEMClassification::EStep()
{
  int i, k;
  double x;

  irtkGaussian* G = new irtkGaussian[_number_of_tissues];

  for (k = 0; k < _number_of_tissues; k++) {
    G[k].Initialise( _mi[k], _sigma[k]);
  }

  _atlas.First();
  _output.First();
  irtkRealPixel *ptr = _input.GetPointerToVoxels();
  irtkRealPixel *pm = _mask.GetPointerToVoxels();
  for (i=0; i< _input.GetNumberOfVoxels(); i++) {
    double* gv = new double[_number_of_tissues];
    double* numerator = new double[_number_of_tissues];
    double denominator=0, temp=0;
    if (*pm == 1) {
      x = *ptr;
      for (k = 0; k < _number_of_tissues; k++) {
        temp = G[k].Evaluate(x);
        gv[k] = temp;
        temp = temp * _atlas.GetValue(k);
        numerator[k] = temp;
        denominator += temp;
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
        }
      }
    } else {
      for (k = 0; k < _number_of_tissues - 1; k++) {
        _output.SetValue(k, 0);
      }
      _output.SetValue(_number_of_tissues - 1, 1);
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

void irtkEMClassification::WStep()
{
  int i,k;
  double num, den;
  cerr<<"Calculating weights ...";
  irtkRealPixel *pi=_input.GetPointerToVoxels();
  irtkRealPixel *pw=_weights.GetPointerToVoxels();
  irtkRealPixel *pwB=_weightsB.GetPointerToVoxels();
  irtkRealPixel *pe=_estimate.GetPointerToVoxels();
  irtkRealPixel *pm=_mask.GetPointerToVoxels();
  _output.First();
  _atlas.First();

  for (i=0; i< _input.GetNumberOfVoxels(); i++) {
    if ((*pm == 1)&&(*pi != _padding)) {
      num=0;
      den=0;
      for (k=0; k<_number_of_tissues; k++) {
        num += _output.GetValue(k)*_mi[k]/_sigma[k];
        den += _output.GetValue(k)/_sigma[k];
      }
      *pw=den;
      *pwB=_output.GetValue(0)/_sigma[0];
      *pe=num/den;

      /*if(_background_tissue >=0)
      {
        if(_atlas.GetValue(_background_tissue) == 1) *pe=_padding;
      }
      */
    } else {
      *pw=_padding;
      *pe=_padding;
    }

    pi++;
    pm++;
    pw++;
    pwB++;
    pe++;
    _output.Next();
    _atlas.Next();
  }
  _estimate.Write("_e.nii.gz");
  _weights.Write("_weights.nii.gz");
  _weightsR.Write("_weightsR.nii.gz");
  _weightsB.Write("_weightsB.nii.gz");
  cerr<<"done."<<endl;
}

void irtkEMClassification::GetMean(double *mean){
	int i;
	for(i=0;i<_number_of_tissues;i++){
		mean[i] = _mi[i];
	}
}

void irtkEMClassification::GetVariance(double *variance){
	int i;
	for(i=0;i<_number_of_tissues;i++){
		variance[i] = sqrt(_sigma[i]);
	}	
}

void irtkEMClassification::BrainmaskInput()
{
  int i;

  cerr<<"brainmasking ...";
  irtkRealPixel *pi=_input.GetPointerToVoxels();
  _atlas.First();
  if (_background_tissue >=0) {
    for (i=0; i< _input.GetNumberOfVoxels(); i++) {
      if (_atlas.GetValue(_background_tissue) == 1) *pi=_padding;

      pi++;
      _atlas.Next();
    }
    cerr<<"done."<<endl;

  } else {
    cerr<<"no background tissue."<<endl;
  }
  CreateMask();
  _mask.Write("_m.nii.gz");

  _input.Write("_i.nii.gz");
}

void irtkEMClassification::MStepGMM(bool uniform_prior)
{
  int i, k;
  vector<double> mi_num(_number_of_tissues);
  vector<double> sigma_num(_number_of_tissues);
  vector<double> denom(_number_of_tissues);
  vector<double> num_vox(_number_of_tissues);

  for (k = 0; k < _number_of_tissues; k++) {
    mi_num[k] = 0;
  }

  for (k = 0; k < _number_of_tissues; k++) {
    sigma_num[k] = 0;
  }

  for (k = 0; k < _number_of_tissues; k++) {
    denom[k] = 0;
  }

  for (k = 0; k < _number_of_tissues; k++) {
    num_vox[k] = 0;
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
        num_vox[k]+= 1;
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
      cerr <<"Tissue "<< k <<": Division by zero while computing tissue mean!" << endl;
      exit(1);
    }
     if (uniform_prior) _c[k]=1.0/_number_of_tissues;
     else _c[k]=denom[k]/num_vox[k];
  }

  _output.First();
  ptr = _input.GetPointerToVoxels();
  pm = _mask.GetPointerToVoxels();

  for (i = 0; i < _input.GetNumberOfVoxels(); i++) {
    // Check for backgound
    if (*pm == 1) {
      for (k = 0; k <_number_of_tissues; k++) {
        sigma_num[k] += (_output.GetValue(k) * (*ptr - _mi[k]) * (*ptr - _mi[k]));
      }
    }
    ptr++;
    pm++;
    _output.Next();
  }

  for (k = 0; k <_number_of_tissues; k++) {
    _sigma[k] = sigma_num[k] / denom[k];
	if(_sigma[k]<1)
		_sigma[k] = 1;
  }
}

void irtkEMClassification::MStepVarGMM(bool uniform_prior)
{
  int i, k;
  vector<double> mi_num(_number_of_tissues);
  double sigma_num;
  vector<double> denom(_number_of_tissues);
  vector<double> num_vox(_number_of_tissues);

  for (k = 0; k < _number_of_tissues; k++) {
    mi_num[k] = 0;
  }

    sigma_num = 0;

  for (k = 0; k < _number_of_tissues; k++) {
    denom[k] = 0;
  }

  for (k = 0; k < _number_of_tissues; k++) {
    num_vox[k] = 0;
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
        num_vox[k]+= 1;
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
      cerr << "Division by zero while computing tissue mean!" << endl;
      //exit(1);
    }
     if (uniform_prior) _c[k]=1.0/_number_of_tissues;
     else _c[k]=denom[k]/num_vox[k];
  }

  _output.First();
  ptr = _input.GetPointerToVoxels();
  pm = _mask.GetPointerToVoxels();

  for (i = 0; i < _input.GetNumberOfVoxels(); i++) {
    // Check for backgound
    if (*pm == 1) {
      for (k = 0; k <_number_of_tissues; k++) {
        sigma_num += (_output.GetValue(k) * (*ptr - _mi[k]) * (*ptr - _mi[k]));
      }
    }
    ptr++;
    pm++;
    _output.Next();
  }

  double sum =0;
  for (k = 0; k <_number_of_tissues; k++) sum += denom[k];
  for (k = 0; k <_number_of_tissues; k++) {
    if (sum>0) _sigma[k] = sigma_num / sum;
  }
}



void irtkEMClassification::EStepGMM(bool uniform_prior)
{
  int i, k;
  double x;
  double *gv = new double[_number_of_tissues];
  double *numerator = new double[_number_of_tissues];
  irtkGaussian *G = new irtkGaussian[_number_of_tissues];

  for (k = 0; k < _number_of_tissues; k++) {
    G[k].Initialise( _mi[k], _sigma[k]);
  }

  _output.First();
  irtkRealPixel *ptr = _input.GetPointerToVoxels();
  irtkRealPixel *pm = _mask.GetPointerToVoxels();

  for (i=0; i< _input.GetNumberOfVoxels(); i++) {
    double denominator=0, temp=0;
    if (*pm == 1) {
      x = *ptr;
      for (k = 0; k < _number_of_tissues; k++) {
        temp = G[k].Evaluate(x);
        gv[k] = temp;
        if (!uniform_prior) temp = temp * _c[k];
        numerator[k] = temp;
        denominator += temp;
      }
      for (k = 0; k < _number_of_tissues; k++) {
        if (denominator != 0) {
          double value = numerator[k]/denominator;
          _output.SetValue(k, value);
          if ((value < 0) || (value > 1)) {
            cerr << "Probability value out of range = " << value << endl;
            exit(1);
          }
        } else {
          cerr << "Division by 0 while computing probabilities" << endl;
          cerr<<"tissue="<<k;
          exit(1);
        }
      }
    } else {
      for (k = 0; k < _number_of_tissues - 1; k++) {
        _output.SetValue(k, 0);
      }
      _output.SetValue(_number_of_tissues - 1, 1);
    }
    ptr++;
    pm++;
    _output.Next();
  }

  delete[] G;
  delete[] gv;
  delete[] numerator;
}

void irtkEMClassification::Print()
{
  int k;

  cerr << "mean: " <<endl;
  for (k = 0;  k <_number_of_tissues; k++) {
    cerr << "tissue " << k << ": ";
    cerr << _mi[k] << endl;
  }

  cerr << "sigma: " << endl;
  for (k = 0; k < _number_of_tissues; k++) {
    cerr << "tissue " << k << ":";
    cerr << sqrt(_sigma[k]) << endl;
  }
}

void irtkEMClassification::PrintGMM()
{
  int k;
  Print();
  cerr << "c: " << endl;
  for (k = 0; k < _number_of_tissues; k++) {
    cerr << "tissue " << k << ":";
    cerr << _c[k] << endl;
  }

}

double irtkEMClassification::Iterate(int)
{
  this->EStep();
  this->MStep();
  Print();
  return LogLikelihood();
}

double irtkEMClassification::IterateGMM(int iteration, bool equal_var, bool uniform_prior)
{
  if (iteration > 1) this->EStepGMM();
  if (equal_var) this->MStepVarGMM(uniform_prior);
  else this->MStepGMM(uniform_prior);
  PrintGMM();

  return LogLikelihoodGMM();
}

double irtkEMClassification::LogLikelihood()
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
        temp += gv[k] * _output.GetValue(k);
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

double irtkEMClassification::LogLikelihoodGMM()
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
        temp += gv[k] * _c[k];
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


void irtkEMClassification::ConstructSegmentation(irtkRealImage &segmentation)
{
  int i, j;
  double max;
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
      for (j = 0; j < _number_of_tissues; j++) {
        if (_output.GetValue(j) > max) {
          max  = _output.GetValue(j);
          *ptr = j+1;
          if ( (j+1) == _number_of_tissues) *ptr=0;
        }
      }
    }
    ptr++;
    _output.Next();
  }
}
void irtkEMClassification::ConstructSegmentationWithPadding(irtkRealImage &segmentation)
{
  int i, j;
  double max;
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
      for (j = 0; j < _number_of_tissues; j++) {
        if (_output.GetValue(j) > max) {
          max  = _output.GetValue(j);
          *ptr = j+1;
          //if ( (j+1) == _number_of_tissues) *ptr=0;
        }
      }
    } else *ptr=_padding;
    ptr++;
    _output.Next();
  }
}
void irtkEMClassification::GetProbMap(int i,irtkRealImage& image){
	if  (i < _number_of_tissues) {
		image = _output.GetImage(i);
	} else {
		cerr << "irtkProbabilisticAtlas::Write: No such probability map" << endl;
		exit(1);
	}
}
void irtkEMClassification::ConstructSegmentation()
{
  int i, j;
  double max;
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
      for (j = 0; j < _number_of_tissues; j++) {
        if (_output.GetValue(j) > max) {
          max  = _output.GetValue(j);
          *ptr = j+1;
          //if ( (j+1) == _number_of_tissues) *ptr=0;
        }
      }
    } else *ptr=_padding;
    ptr++;
    _output.Next();
  }
}

void irtkEMClassification::ConstructSegmentationFromPV()
{
  int i, j;
  double max;
  irtkRealPixel *ptr,*ptrm;

  if (_debug) cerr<<"Updating segmentation with PV"<<endl;

  // Initialize pointers of probability maps
  _output.First();
  _pv_output.First();
  ptr = _segmentation.GetPointerToVoxels();
  ptrm = _mask.GetPointerToVoxels();


  for (i = 0; i < _input.GetNumberOfVoxels(); i++) 
  {
    if (*ptrm>_padding) 
    {
      max  = 0;
      *ptr = 0;
      for (j = 0; j < _number_of_tissues; j++) 
        if (_pv_output.GetValue(j) > max) 
        {
          max  = _pv_output.GetValue(j);
          *ptr = j+1;
        }
    }
    ptr++;
    ptrm++;
    _pv_output.Next();
  }
}

void irtkEMClassification::ConstructSegmentationNoBG(irtkRealImage &segmentation)
{
  int i, j;
  double max;
  irtkRealPixel *ptr;

  cerr<<"Constructing segmentation NoBG"<<endl;

  // Initialize pointers of probability maps
  _output.First();

  // Initialize segmentation to same size as input
  segmentation = _input;
  ptr = segmentation.GetPointerToVoxels();

  for (i = 0; i < _input.GetNumberOfVoxels(); i++) {
    if (*ptr>_padding) {
      max  = 0;
      *ptr = 0;
      for (j = 0; j < _number_of_tissues; j++) {
        if (_output.GetValue(j) > max) {
          max  = _output.GetValue(j);
          *ptr = j+1;
        }
      }
    }
    ptr++;
    _output.Next();
  }
}

void irtkEMClassification::WriteProbMap(int i, const char *filename)
{
  _output.Write(i, filename);
}

void irtkEMClassification::WriteGaussianParameters(const char *file_name, int flag)
{
  cerr << "Writing GaussianDistributionParameters: " << file_name << endl;

  ofstream fileOut(file_name);

  if (!fileOut) {
    cerr << "Can't open file " << file_name << endl;
    exit(1);
  }

  int k,l,m;

  if(flag){
	  // out put without names
	  for (k=0; k<_number_of_tissues; k++) {
			  fileOut << _mi[k] << " " << _sigma[k] << endl;
	  }
  }else{
	  // out put with names
	  fileOut << "mi: " <<endl;
	  for (k=0; k<_number_of_tissues; k++) {
		  fileOut << "Tissue " << k << ": (";
		  for (l=0; l < 1/*_input.GetNumberOfChannels()*/; l++) {
			  fileOut << _mi[k];//.Get(l);
			  if (l == 0/*_input.GetNumberOfChannels() - 1*/) fileOut << ")" << endl;
			  else fileOut << ", ";
		  }
	  }

	  fileOut << "sigma: " << endl;
	  for (k=0; k<_number_of_tissues; k++) {
		  fileOut << "Tissue " << k << ": (";
		  //<<endl << "(";

		  for (l=0; l < 1/*_input.GetNumberOfChannels()*/; l++) {
			  //fileOut << "(";
			  for (m = 0; m < 1/*_input.GetNumberOfChannels()*/; m++) {
				  double s = _sigma[k];//.Get(m,l);
				  if ( s >= 0) fileOut << sqrt(s);
				  else fileOut << -sqrt(-s);
				  //if (m < 1/*_input.GetNumberOfChannels() - 1*/) fileOut << ", ";
				  //else 
				  //fileOut <<  ")" << endl;
			  }
		  }
		  //if (l == 0/*_input.GetNumberOfChannels() - 1*/) fileOut << ")" << endl;
		  fileOut <<  ")" << endl;
	  }
  }
}

void irtkEMClassification::WriteWeights(const char *filename)
{
  irtkRealImage w(_weights);
  irtkRealPixel *pw = w.GetPointerToVoxels();
  int i;
  double sigma_min=_sigma[0];

  for (i=1; i<_number_of_tissues; i++) if (sigma_min > _sigma[i]) sigma_min = _sigma[i];
  cerr<<"sigma min = "<<sigma_min<<endl;

  for (i=0; i<w.GetNumberOfVoxels(); i++) {
    if (*pw != _padding) *pw=(*pw) * sigma_min * 100;
    pw++;
  }

  w.Write(filename);
}

double irtkEMClassification::PointLogLikelihoodGMM(double x)
{
  int k;
  double temp=0;

  for (k=0; k < _number_of_tissues; k++) temp+= _G[k].Evaluate(x) * _c[k];
  if (-log(temp)> 1000000) exit(1);

  if ((temp > 1) || (temp < 0)) {
    cerr << "Could not compute likelihood, probability out of range = " << temp << endl;
    exit(1);
  }
  return -log(temp);
}

double irtkEMClassification::PointLogLikelihoodGMM(double x, double)
{
  cerr<<"May be we do not want this one ;-)"<<endl;
  return PointLogLikelihoodGMM(x);
}

double irtkEMClassification::PointLogLikelihoodGMMnomatch(double x, double)
{
  cerr<<"May be we do not want this one ;-)"<<endl;
  exit(1);
  return PointLogLikelihoodGMM(x);
}



void irtkEMClassification::GInit()
{
  if (_G!=NULL) delete[] _G;
  _G = new irtkGaussian[_number_of_tissues];

  for (int k = 0; k < _number_of_tissues; k++) {
    _G[k].Initialise( _mi[k], _sigma[k]);
  }
}

bool irtkEMClassification::PVStep(int wm1Label, int wm2Label, int cortexLabel, int csfLabel, int backgroundLabel, double lcc_treshold, bool first, irtkGreyImage *)
{
  int i,j,k,l,n;
  int labelId[5], labelCount[4];
  double lambda = 0, wm1Prior,wm2Prior,cortexPrior,csfPrior,cortexNew,csfNew,backgroundPrior;
  bool reduceWM = 0, reduceGM = 0;
  int increaseBG = 0;
  bool change=true;
  int offset[26][3] = {{1,0,0},{0,1,0},{0,0,1},{-1,0,0},{0,-1,0},{0,0,-1},{1,1,0},{1,0,1},{0,1,1},{1,-1,0},{1,0,-1},{0,1,-1},{-1,1,0},{-1,0,1},{0,-1,1},{-1,-1,0},{-1,0,-1},{0,-1,-1},{1,1,1},{1,1,-1},{1,-1,1},{-1,1,1},{1,-1,-1},{-1,1,-1},{-1,-1,1},{-1,-1,-1}};
  int neighbour_num = 26;
  double PVtissues[10];


  cerr<<"PVStep ... ";//<<endl;

  labelId[0]=wm1Label;
  labelId[1]=cortexLabel;
  labelId[2]=csfLabel;
  labelId[3]=backgroundLabel;
  labelId[4]=wm2Label;

  if (first)
  {
    if (_debug) cerr<<"First iteration."<<endl;
    _brain = _mask;
    irtkGreyImage bg(_segmentation);
    irtkMeanShift mshbg(bg);
    mshbg.SetOutput(&bg);
    mshbg.Lcc(1);
    bg.Write("bglcc.nii.gz");
    for(k=0;k<_segmentation.GetZ();k++)
      for(j=0;j<_segmentation.GetY();j++)
        for(i=0;i<_segmentation.GetX();i++)
        {
          if ( _segmentation.Get(i,j,k)==csfLabel)
          {
            _atlas.SetValue(i,j,k,backgroundLabel-1,0);
            _atlas.SetValue(i,j,k,cortexLabel-1,0);
            _atlas.SetValue(i,j,k,wm1Label-1,0);
            _atlas.SetValue(i,j,k,wm2Label-1,0);
            _atlas.SetValue(i,j,k,csfLabel-1,1);
            _brain.Put(i,j,k,0);
          }
          if (( _segmentation.Get(i,j,k)==backgroundLabel)&&(bg.Get(i,j,k)==1))
          {
            _atlas.SetValue(i,j,k,backgroundLabel-1,1);
            _atlas.SetValue(i,j,k,cortexLabel-1,0);
            _atlas.SetValue(i,j,k,wm1Label-1,0);
            _atlas.SetValue(i,j,k,wm2Label-1,0);
            _atlas.SetValue(i,j,k,csfLabel-1,0);
            _brain.Put(i,j,k,0);
          }
       }
       if (_debug) _brain.Write("brain1.nii.gz");
  }
  

  //largest connected component on WM, brain (WM and GM) and background
  irtkGreyImage wm(_segmentation);
  irtkGreyImage gm(_segmentation);
  for(k=0;k<_segmentation.GetZ();k++)
    for(j=0;j<_segmentation.GetY();j++)
      for(i=0;i<_segmentation.GetX();i++)
        if (( _segmentation.Get(i,j,k)==wm1Label)||( _segmentation.Get(i,j,k) == wm2Label )) 
        {
          wm.Put(i,j,k,1);
          gm.Put(i,j,k,1);
        }
        else 
        {
          wm.Put(i,j,k,0);
          if ( _segmentation.Get(i,j,k)==cortexLabel) gm.Put(i,j,k,1);
          else gm.Put(i,j,k,0);
        }
 
   irtkDilation<irtkGreyPixel> dilation;
   irtkErosion<irtkGreyPixel> erosion;


   if (_debug) wm.Write("wm.nii.gz");
   irtkMeanShift msh(wm);
   irtkGreyImage wmlcc(wm);
   msh.SetOutput(&wmlcc);
   msh.LccS(1, lcc_treshold);
   if (_debug) wmlcc.Write("wm-lcc.nii.gz");

   irtkGreyImage gmlcc(gm);
   if (_debug) gm.Write("gm2.nii.gz");
   erosion.SetInput(&gm);
   erosion.SetOutput(&gm);
   erosion.Run();
   if (_debug) gm.Write("gm-e-lcc.nii.gz");
   irtkMeanShift msh2(gm);
   msh2.SetOutput(&gmlcc);
   msh2.Lcc(1);
   dilation.SetInput(&gmlcc);
   dilation.SetOutput(&gmlcc);
   dilation.Run();
   if (_debug) gmlcc.Write("gm-lcc.nii.gz");

  //DistanceTransform(3,4);

  change = false;
  for(k=0;k<_segmentation.GetZ();k++)
    for(j=0;j<_segmentation.GetY();j++)
      for(i=0;i<_segmentation.GetX();i++)
      {

        if(_mask.GetAsDouble(i,j,k)==1)
        {

          //find which tissues to reduce

          reduceWM = 0;
          reduceGM = 0;
          increaseBG = 0;

          if((_segmentation.GetAsDouble(i,j,k)==wm1Label) || (_segmentation.GetAsDouble(i,j,k)==wm2Label))
            if (wmlcc.Get(i,j,k)==0) 
            {
              reduceWM = 1;
              change=true;
            }
          if((_segmentation.GetAsDouble(i,j,k)==wm1Label) || (_segmentation.GetAsDouble(i,j,k)==wm2Label) || (_segmentation.GetAsDouble(i,j,k)==cortexLabel))
            if (gmlcc.Get(i,j,k)==0)
            { 
              reduceWM = 1;
              reduceGM = 1;
              //increaseBG = 1;
              change=true;
             }

          //Calculate number of labels in 26-neighbourhood
          for(l=0;l<=4;l++) labelCount[l]=0;
          for(n=0;n<neighbour_num;n++)
            for(l=0;l<=4;l++) 
              if(((i+offset[n][0])>=0) && ((i+offset[n][0])<_segmentation.GetX()) && ((j+offset[n][1])>=0) && ((j+offset[n][1])<_segmentation.GetY()) && ((k+offset[n][2])>=0) && ((k+offset[n][2])<_segmentation.GetZ())) 
                if (_segmentation.Get(i+offset[n][0],j+offset[n][1],k+offset[n][2])==labelId[l]) labelCount[l]++;

          double sum = 0;
          for(n=0;n<_number_of_tissues;n++)
            for(l=0;l<n;l++) 
            {
              PVtissues[round(n*(n-1)/2) + l] = labelCount[n]*labelCount[l];
              sum += labelCount[n]*labelCount[l];
            }
          for(n=0;n<_number_of_tissues;n++)
            for(l=0;l<n;l++) 
            {
              PVtissues[round(n*(n-1)/2) + l] /= sum;
            }
          

          //background for WM voxels
          if((_segmentation.GetAsDouble(i,j,k)==wm1Label)||(_segmentation.GetAsDouble(i,j,k)==wm2Label))
            if(labelCount[3]>0) 
            {
              reduceWM = 1;
              //increaseBG = 1;
            if((_segmentation.GetAsDouble(i,j,k)==wm1Label)||(_segmentation.GetAsDouble(i,j,k)==wm2Label))
              change = true;
            }

           //identify csf-gm boundary csf index = 2 and gm index = 1 => csf-gm PV index =  2
          if((_segmentation.GetAsDouble(i,j,k)==wm1Label)||(_segmentation.GetAsDouble(i,j,k)==wm2Label))
             if (PVtissues[2]>0.33) 
             { 
               reduceWM = 1;
               change = true;
             }

           //identify csf-gm boundary csf index = 2 and gm index = 1 => csf-gm PV index =  2, with degree of certanty about wm
           if((_segmentation.GetAsDouble(i,j,k)==wm1Label)||(_segmentation.GetAsDouble(i,j,k)==wm2Label))
             if ((PVtissues[2]>0)&&((labelCount[0]+labelCount[3])<13) )
             { 
               reduceWM = 1;               
               change = true;
             }

           //identify csf-bg boundary csf index = 2 and bg index = 3 => csf-bg PV index =  5 
          //also consider wm-bg boundary as csf can be misclassified for wm. wm index = 0,5 and bg index = 3 => bg-wm PV index =  3,9
           if(_segmentation.GetAsDouble(i,j,k)==cortexLabel)
             if ((PVtissues[3]>0.33)||(PVtissues[9]>0.33)||(PVtissues[5]>0.33))
             { 
               reduceWM = 1;
               reduceGM = 1;
               change = true;
             }

           //identify csf-wm1 boundary csf index = 2 and wm2 index = 0 => csf-wm1 PV index =  1
           if((_segmentation.GetAsDouble(i,j,k)==wm2Label))
             if (PVtissues[1]>0.33) 
             { 
               reduceWM = 1;
               reduceGM = 1;
               change = true;
             }

          cortexPrior = _atlas.GetValue(i,j,k,cortexLabel-1);
          wm1Prior = _atlas.GetValue(i,j,k,wm1Label-1);
          wm2Prior = _atlas.GetValue(i,j,k,wm2Label-1);
          csfPrior = _atlas.GetValue(i,j,k,csfLabel-1);
          backgroundPrior = _atlas.GetValue(i,j,k,backgroundLabel-1);


          if ((reduceWM)&&(!reduceGM))
          {
            _atlas.SetValue(i,j,k,wm1Label-1,lambda*wm1Prior);
            _atlas.SetValue(i,j,k,wm2Label-1,lambda*wm2Prior);
            if ((cortexPrior>0)||(csfPrior>0))
            {
             cortexNew = cortexPrior+ (1-lambda)*(wm1Prior+wm2Prior)*cortexPrior/(cortexPrior+csfPrior);
             csfNew = csfPrior + (1-lambda)*(wm1Prior+wm2Prior)*csfPrior/(cortexPrior+csfPrior);
            }
            else
            {
	      cortexNew = cortexPrior+ (1-lambda)*(wm1Prior+wm2Prior)*0.5;
              csfNew = csfPrior+ (1-lambda)*(wm1Prior+wm2Prior)*0.5;
            }
            _atlas.SetValue(i,j,k,cortexLabel-1,cortexNew);
            _atlas.SetValue(i,j,k,csfLabel-1,csfNew);
            _brain.Put(i,j,k,0);
           }

            if ((reduceWM)&&(reduceGM))
            {
              _atlas.SetValue(i,j,k,wm1Label-1,lambda*wm1Prior);
              _atlas.SetValue(i,j,k,wm2Label-1,lambda*wm2Prior);
              _atlas.SetValue(i,j,k,cortexLabel-1,lambda*cortexPrior);
              csfNew = csfPrior + (1-lambda)*(wm1Prior+wm2Prior+cortexPrior);
              _atlas.SetValue(i,j,k,csfLabel-1,csfNew);
              _brain.Put(i,j,k,0);
            }
          }
        }

    if (_debug) _atlas.Write(0,"csfPVmod.nii.gz");
    if (_debug) _atlas.Write(1,"gmPVmod.nii.gz");
    if (_debug) _atlas.Write(2,"wm1PVmod.nii.gz");
    if (_debug) _atlas.Write(3,"wm2PVmod.nii.gz");
    if (_debug) _atlas.Write(4,"bgPVmod.nii.gz");
    if (_debug) _brain.Write("brain2.nii.gz");
    
    EStep();
    ConstructSegmentationWithPadding(_segmentation);
    if (change) cerr<<"there was change."<<endl;
    else cerr<<"no more change."<<endl;
    return change;

}

void irtkEMClassification::WritePVProbMap(int i, const char *filename)
{
  _pv_output.Write(i, filename);
}


void irtkEMClassification::ConstructPVSegmentation()
{
  int i,j,k,l,n;
  double temp;

//Please note that in this approach we trust the segmentation a lot so we do not assume that there is noise (this is reasonable for modern MR)
//Single voxel surounded by other tissues will therefore be considered as correctly classified rather than noisy

  if (_debug) cerr<<"Constructing PV segmentation"<<endl;

  vector<int> labelCount(_number_of_tissues);
  int offset[26][3] = {{1,0,0},{0,1,0},{0,0,1},{-1,0,0},{0,-1,0},{0,0,-1},{1,1,0},{1,0,1},{0,1,1},{1,-1,0},{1,0,-1},{0,1,-1},{-1,1,0},{-1,0,1},{0,-1,1},{-1,-1,0},{-1,0,-1},{0,-1,-1},{1,1,1},{1,1,-1},{1,-1,1},{-1,1,1},{1,-1,-1},{-1,1,-1},{-1,-1,1},{-1,-1,-1}};
  int neighbour_num = 26;
  int pure_treshold = 19;

  int label1,label2;


  for(k=0;k<_segmentation.GetZ();k++)
    for(j=0;j<_segmentation.GetY();j++)
      for(i=0;i<_segmentation.GetX();i++)
      {
        if(_mask.GetAsDouble(i,j,k)==1)
        {
          //initialise number of labels in neighbourhood
          for(l=0;l<_number_of_tissues;l++) labelCount[l]=0;

          //count number of labels in neighbourhood
          for(n=0;n<neighbour_num;n++)
            for(l=0;l<_number_of_tissues;l++) 
              if(((i+offset[n][0])>=0) && ((i+offset[n][0])<_segmentation.GetX()) && ((j+offset[n][1])>=0) && ((j+offset[n][1])<_segmentation.GetY()) && ((k+offset[n][2])>=0) && ((k+offset[n][2])<_segmentation.GetZ())) 
                if (_segmentation.Get(i+offset[n][0],j+offset[n][1],k+offset[n][2]) == (l+1)) labelCount[l]++;
           //also count the current voxels
           labelCount[round(_segmentation.Get(i,j,k)-1)]++;

          //initialise pv segmentation
          for(l=0;l<_number_of_tissues;l++) _pv_output.SetValue(i,j,k,l,0);

          //pure tissue
          if(labelCount[round(_segmentation.Get(i,j,k)-1)]>=pure_treshold)
            _pv_output.SetValue(i,j,k,round(_segmentation.Get(i,j,k)-1),1);
          //mixed tissue
          else
          {
           //find 2 main tissues in neighbourhood
            label1=0; label2=1;
            if(labelCount[label2]>labelCount[label1])
            {label1=1;label2=0;}
            for(l=2;l<_number_of_tissues;l++) 
              if (labelCount[l]>=labelCount[label1])
              { 
                label2=label1;
                label1=l;
              }
              else if (labelCount[l]>=labelCount[label2]) label2 = l;
            //Calculate mixing proportions
            temp = (_mi[label2]-_input.Get(i,j,k))/(_mi[label2]-_mi[label1]);
           // cerr<<temp<<" ";
            if (temp<0) temp=0;
            if (temp>1) temp=1;
            _pv_output.SetValue(i,j,k,label1,temp);
            _pv_output.SetValue(i,j,k,label2,1-temp);
            //cerr<<temp<<" ";
          }
        }
      }
  if (_debug) _pv_output.Write(0, "csfPV.nii.gz");
  if (_debug) _pv_output.Write(1, "gmPV.nii.gz");
  if (_debug) _pv_output.Write(2, "wm1PV.nii.gz");
  if (_debug) _pv_output.Write(3, "wm2PV.nii.gz");
  if (_debug) _pv_output.Write(4, "bgPV.nii.gz");

}

void irtkEMClassification::ConstructSegmentationBrainNonBrain(int wm1Label, int wm2Label, int cortexLabel, int csfLabel, int backgroundLabel, irtkGreyImage * smask)
{
  cerr<<"Constructing segmentation Brain-NonBrain"<<endl;
  int i,j,k,l,n;
  int labelId[5], labelCount[5], nonBrainCount;
  //int offset[6][3] = {{1,0,0},{0,1,0},{0,0,1},{-1,0,0},{0,-1,0},{0,0,-1}};
  //int neighbour_num = 6;
  int offset[26][3] = {{1,0,0},{0,1,0},{0,0,1},{-1,0,0},{0,-1,0},{0,0,-1},{1,1,0},{1,0,1},{0,1,1},{1,-1,0},{1,0,-1},{0,1,-1},{-1,1,0},{-1,0,1},{0,-1,1},{-1,-1,0},{-1,0,-1},{0,-1,-1},{1,1,1},{1,1,-1},{1,-1,1},{-1,1,1},{1,-1,-1},{-1,1,-1},{-1,-1,1},{-1,-1,-1}};
  int neighbour_num = 26;
  int pure_treshold = 19;
  int nonbrain_treshold = 19;
  int nontissue_treshold = 7;
  double distance_treshold = 3;
  double bg_distance_treshold = 3;

  labelId[0]=wm1Label;
  labelId[1]=cortexLabel;
  labelId[2]=csfLabel;
  labelId[3]=backgroundLabel;
  labelId[4]=wm2Label;
  double alpha;
  _brain.Write("brain.nii.gz");
  irtkGreyImage braincount(_brain);

  if(smask != NULL) _smask = *smask;
  else _smask = _mask;

  irtkGreyImage background(_segmentation);
  irtkGreyPixel * ptr = background.GetPointerToVoxels();
  for (i=0; i<background.GetNumberOfVoxels();i++)
  {
    if(*ptr <= 1) *ptr = 1;
    else *ptr = 0;
    ptr++;
  }
  background.Write("background.nii.gz");
  irtkMeanShift msh(background);
  msh.SetOutput(&background);
  msh.Lcc(1);
  background.Write("background-lcc.nii.gz");
  irtkGreyImage pvseg(_segmentation);
  _segmentation = background;
  DistanceTransform(1,1);
  background = _distance;
  background.Write("background-distance.nii.gz");

  _segmentation = pvseg;
  
  DistanceTransform(3,4);
  
 //cerr<<_segmentation.GetX()<<" "<<_segmentation.GetY()<<" "<<_segmentation.GetZ()<<endl;
  for(k=0;k<_segmentation.GetZ();k++)
    for(j=0;j<_segmentation.GetY();j++)
      for(i=0;i<_segmentation.GetX();i++)
      {
        if(_mask.GetAsDouble(i,j,k)==1)
        {
          for(l=0;l<=4;l++) labelCount[l]=0;
          nonBrainCount = 0;

          for(n=0;n<neighbour_num;n++)
            for(l=0;l<=4;l++) 
              if(((i+offset[n][0])>=0) && ((i+offset[n][0])<_segmentation.GetX()) && ((j+offset[n][1])>=0) && ((j+offset[n][1])<_segmentation.GetY()) && ((k+offset[n][2])>=0) && ((k+offset[n][2])<_segmentation.GetZ())) 
                if (_segmentation.Get(i+offset[n][0],j+offset[n][1],k+offset[n][2])==labelId[l]) labelCount[l]++;

          for(n=0;n<neighbour_num;n++)
            if(((i+offset[n][0])>=0) && ((i+offset[n][0])<_segmentation.GetX()) && ((j+offset[n][1])>=0) && ((j+offset[n][1])<_segmentation.GetY()) && ((k+offset[n][2])>=0) && ((k+offset[n][2])<_segmentation.GetZ())) 
              if (_brain.Get(i+offset[n][0],j+offset[n][1],k+offset[n][2])==0) nonBrainCount++;

              braincount.Put(i,j,k,nonBrainCount);


              //csf-bg
              if ((nonBrainCount>nonbrain_treshold) || ((_segmentation.Get(i,j,k) == backgroundLabel)&&(_brain.Get(i,j,k)==0)) || ((_distance.Get(i,j,k)>distance_treshold)&&(background(i,j,k)<bg_distance_treshold)))
              {
                alpha = (_input.Get(i,j,k)-_mi[backgroundLabel-1])/(_mi[csfLabel-1]-_mi[backgroundLabel-1]);
                if (labelCount[3]>=pure_treshold) alpha = 0;
                if (labelCount[2]>=pure_treshold) alpha = 1;
                if (alpha < 0) alpha = 0;
                if (alpha > 1) alpha = 1;
                _pv_output.SetValue(i,j,k,backgroundLabel-1,1-alpha);
                _pv_output.SetValue(i,j,k,csfLabel-1,alpha);
                _pv_output.SetValue(i,j,k,cortexLabel-1,0);
                _pv_output.SetValue(i,j,k,wm1Label-1,0);
                _pv_output.SetValue(i,j,k,wm2Label-1,0);
              }
              else
              {
                //gm-wm1-wm2
		//it is brain tissue and distance from wm is within limits and there is WM tissue in neighbourhood 
                if ((_brain.Get(i,j,k)==1)&&(_distance.Get(i,j,k)<=distance_treshold)&&((labelCount[0]+labelCount[4])>nontissue_treshold))
                {
                  //gm-wm1
                  alpha = (_input.Get(i,j,k)-_mi[cortexLabel-1])/(_mi[wm1Label-1]-_mi[cortexLabel-1]);
                  if (labelCount[1]>=pure_treshold) alpha = 0;
                  if (labelCount[0]>=pure_treshold) alpha = 1;
                  if (alpha < 0) alpha = 0;
                  if (alpha <= 1)
                  {
                    _pv_output.SetValue(i,j,k,cortexLabel-1,1-alpha);
                    _pv_output.SetValue(i,j,k,wm1Label-1,alpha);
                    _pv_output.SetValue(i,j,k,backgroundLabel-1,0);
                    _pv_output.SetValue(i,j,k,csfLabel-1,0);
                    _pv_output.SetValue(i,j,k,wm2Label-1,0);
                  }
                  //wm1-wm2
                  else
                  {
                    alpha = (_input.Get(i,j,k)-_mi[wm1Label-1])/(_mi[wm2Label-1]-_mi[wm1Label-1]);
                    if (labelCount[4]>=pure_treshold) alpha = 1;
                    if (alpha > 1) alpha = 1;
                    _pv_output.SetValue(i,j,k,wm1Label-1,1-alpha);
                    _pv_output.SetValue(i,j,k,wm2Label-1,alpha);
                    _pv_output.SetValue(i,j,k,backgroundLabel-1,0);
                    _pv_output.SetValue(i,j,k,csfLabel-1,0);
                    _pv_output.SetValue(i,j,k,cortexLabel-1,0);
                   }
                 }
                 //csf-gm
		 //if it is non-brain or distance from WM is outside limit or there is not WM tissue in neigbourhood
                 //if (_brain.Get(i,j,k)==0) //&&(nonBrainCount<=nonbrain_treshold))
                 else
                 {
		   //but only in the proximity of non-brain tissue csf is considered
		   if((nonBrainCount>nontissue_treshold)||(_brain.Get(i,j,k)==0))
		   {
                     alpha = (_input.Get(i,j,k)-_mi[cortexLabel-1])/(_mi[csfLabel-1]-_mi[cortexLabel-1]);
                     if (alpha < 0) alpha = 0;
                     if (alpha > 1) alpha = 1;
                     _pv_output.SetValue(i,j,k,cortexLabel-1,1-alpha);
                     _pv_output.SetValue(i,j,k,csfLabel-1,alpha);
                     _pv_output.SetValue(i,j,k,backgroundLabel-1,0);
                     _pv_output.SetValue(i,j,k,wm1Label-1,0);
                     _pv_output.SetValue(i,j,k,wm2Label-1,0);
		   }
		   //only GM can be present here
		   else
		   {
                     _pv_output.SetValue(i,j,k,cortexLabel-1,1);
                     _pv_output.SetValue(i,j,k,csfLabel-1,0);
                     _pv_output.SetValue(i,j,k,backgroundLabel-1,0);
                     _pv_output.SetValue(i,j,k,wm1Label-1,0);
                     _pv_output.SetValue(i,j,k,wm2Label-1,0);
		   }
                 }
              }
    }
  }
  if (_debug) braincount.Write("nonbraincount.nii.gz");
  if (_debug) _distance.Write("distance.nii.gz");
  ConstructSegmentationFromPV();
}


void irtkEMClassification::DistanceTransform( int label1, int label2)
{
  int i,j,k,x,y,z;
  double dx,dy,dz;
  double dv[3];
  double d,value;
  _distance=_mask;
  int offset[26][3] = {{1,0,0},{0,1,0},{0,0,1},{-1,0,0},{0,-1,0},{0,0,-1},{1,1,0},{1,0,1},{0,1,1},{1,-1,0},{1,0,-1},{0,1,-1},{-1,1,0},{-1,0,1},{0,-1,1},{-1,-1,0},{-1,0,-1},{0,-1,-1},{1,1,1},{1,1,-1},{1,-1,1},{-1,1,1},{1,-1,-1},{-1,1,-1},{-1,-1,1},{-1,-1,-1}};
  int neighbour_num = 26;

  _input.GetPixelSize(&dx,&dy,&dz);
  dv[0]=dx;
  dv[1]=dy;
  dv[2]=dz;

  irtkRealPixel *pd = _distance.GetPointerToVoxels();
  irtkRealPixel *pm = _mask.GetPointerToVoxels();
  irtkGreyPixel *psm = _smask.GetPointerToVoxels();
  irtkRealPixel *ps = _segmentation.GetPointerToVoxels();

  //initialize distance map
  for(i=0;i<_distance.GetNumberOfVoxels();i++)
  {
    if((*pm==0)||(*psm==1)) *pd = -1;
    else  *pd = 1000;
    if ((*pd>=0)&&((*ps == label1)||(*ps==label2))) *pd =0;
    pd++;
    pm++;
    psm++;
    ps++;
  }
  if (_debug) _distance.Write("d1.nii.gz");
  if (_debug) _smask.Write("smask.nii.gz");

  //initialise queue for calculating distance
  queue<irtkPoint> voxels;

  for(k=0;k<_distance.GetZ();k++)
    for(j=0;j<_distance.GetY();j++)
      for(i=0;i<_distance.GetX();i++)
        if(_distance.GetAsDouble(i,j,k)==0)
        {
          irtkPoint p(i,j,k);
          voxels.push(p);
        }
  while (!voxels.empty())
  {
    x=voxels.front()._x;
    y=voxels.front()._y;
    z=voxels.front()._z;
    value = _distance.Get(x,y,z);

    for (i=0; i<neighbour_num; i++)
    {
      d=0;
      for(j=0;j<3;j++) d+=offset[i][j]*offset[i][j]*dv[j]*dv[j];
      d=sqrt(d);
      if (((x+offset[i][0])>=0)&&((x+offset[i][0])<_distance.GetX()) && ((y+offset[i][1])>=0) && ((y+offset[i][1])<_distance.GetY()) && ((z+offset[i][2])>=0) && ((z+offset[i][2])<_distance.GetZ()))
        if (_distance.Get(x+offset[i][0],y+offset[i][1],z+offset[i][2]) > (value+d))
        {
          //cerr<<x<<" "<<y<<" "<<z<<endl;
          _distance.Put(x+offset[i][0],y+offset[i][1],z+offset[i][2],value+d);
          irtkPoint p(x+offset[i][0],y+offset[i][1],z+offset[i][2]);
          voxels.push(p);
        }
    }
    voxels.pop();
  }

  if (_debug) _distance.Write("d2.nii.gz");
}

void irtkEMClassification::WritePVSegmentation()
{
  for (int i=0; i<_number_of_tissues;i++)
  {
    char name[100];
    sprintf(name, "tissue%d.nii.gz",i);
    _pv_output.Write(i,name);
  }
}

