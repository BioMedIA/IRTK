/*=========================================================================

  Library   : Image Registration Toolkit (IRTK)
  Module    : $Id$
  Copyright : Imperial College, Department of Computing
              Visual Information Processing (VIP), 2008 onwards
  Date      : $Date$
  Version   : $Revision$
  Changes   : $Author$

=========================================================================*/

#include <irtkImage.h>

#include <irtkEMClassification.h>

#include <irtkGaussian.h>

#include <irtkHistogram.h>


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

}

irtkEMClassification::irtkEMClassification(int noTissues, irtkRealImage **atlas, irtkRealImage *background)
{
  _atlas.AddProbabilityMaps(noTissues, atlas);
  _output.AddProbabilityMaps(noTissues, atlas);
  if (background == NULL) {
    _atlas.AddBackground();
    _output.AddBackground();
  } else {
    _atlas.NormalizeAtlas(*background);
    _output.NormalizeAtlas(*background);
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
  _mask.Write("_mask.nii.gz");
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
  InitialiseAtlas();

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


          //exit(1);
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
  double num, den, dn[2];
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
void irtkEMClassification::BrainmaskInput()
{
  int i,k;
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

void irtkEMClassification::MStepGMM()
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
    _c[k]=denom[k]/num_vox[k];
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

void irtkEMClassification::EStepGMM()
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
        temp = temp * _c[k];
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

double irtkEMClassification::Iterate(int iteration)
{
  this->EStep();
  this->MStep();
  return LogLikelihood();
}

double irtkEMClassification::IterateGMM(int iteration)
{
  if (iteration > 1) this->EStepGMM();
  this->MStepGMM();
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
  /*  if (_number_of_tissues != 3)
    {
      cerr<<"Expects 3 tissue classes."<<endl;
      return;
    }
  */
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

void irtkEMClassification::WriteProbMap(int i, char *filename)
{
  _output.Write(i, filename);
}

void irtkEMClassification::WriteGaussianParameters(char *file_name)
{
  cerr << "Writing GaussianDistributionParameters: " << file_name << endl;

  ofstream fileOut(file_name);

  if (!fileOut) {
    cerr << "Can't open file " << file_name << endl;
    exit(1);
  }

  int k,l,m;

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
    fileOut << "Tissue " << k << ":" <<endl << "(";

    for (l=0; l < 1/*_input.GetNumberOfChannels()*/; l++) {
      fileOut << "(";
      for (m = 0; m < 1/*_input.GetNumberOfChannels()*/; m++) {
        double s = _sigma[k];//.Get(m,l);
        if ( s >= 0) fileOut << sqrt(s);
        else fileOut << -sqrt(-s);
        if (m < 1/*_input.GetNumberOfChannels() - 1*/) fileOut << ", ";
        else fileOut <<  ")" << endl;
      }
    }
    if (l == 0/*_input.GetNumberOfChannels() - 1*/) fileOut << ")" << endl;
  }
}

void irtkEMClassification::WriteWeights(char *filename)
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

double irtkEMClassification::PointLogLikelihoodGMM(double x, double y)
{
  cerr<<"May be we do not want this one ;-)"<<endl;
  return PointLogLikelihoodGMM(x);
}

double irtkEMClassification::PointLogLikelihoodGMMnomatch(double x, double y)
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


