/*=========================================================================

  Library   : Image Registration Toolkit (IRTK)
  Module    : $Id$
  Copyright : Imperial College, Department of Computing
              Visual Information Processing (VIP), 2008 onwards
  Date      : $Date$
  Version   : $Revision$
  Changes   : $Author$

=========================================================================*/

#include <irtkMultiChannelImage.h>
#include <irtkImage.h>
#include <irtkHistogram.h>
#include <irtkGaussianBlurring.h>
#include <irtkGaussianBlurringWithPadding.h>

///Constructor

void irtkMultiChannelImage::AddImage(const irtkRealImage &image)
{
  if (_images.size() == 0) _number_of_voxels = image.GetNumberOfVoxels();
  else
    if (_number_of_voxels != image.GetNumberOfVoxels()) {
      cerr << "Image sizes mismatch" << endl;
      cerr<<_number_of_voxels<<" "<<image.GetNumberOfVoxels()<<endl;
      exit(1);
    }

  _images.push_back(image);
  _pointers.push_back(image.GetPointerToVoxels());

}

void irtkMultiChannelImage::SetImage(int channel, irtkRealImage &image)
{
  if (channel < (int) _images.size()) {
    _images[channel]=image;
    _pointers[channel]=image.GetPointerToVoxels();
  } else {
    cerr << "channel "<<channel<<"does not exist" << endl;
    exit(1);
  }
}

void irtkMultiChannelImage::SetMask(irtkRealImage &mask)
{
    _mask=mask;
}


irtkRealImage& irtkMultiChannelImage::GetImage(int channel)
{
  return _images[channel];
}

void irtkMultiChannelImage::First()
{
  int i;
  for (i=0; i<(int)_pointers.size(); i++) _pointers[i] = _images[i].GetPointerToVoxels();
}

void irtkMultiChannelImage::Next()
{
  int i;
  for (i=0; i<(int)_pointers.size(); i++) _pointers[i]++;
}

irtkRealPixel irtkMultiChannelImage::GetValue(int channel)
{
  if ((channel >= 0) && (channel < (int)_images.size())) return *_pointers[channel];
  else {
    cerr << "Channel identificator " << channel <<" out of range." <<endl;
    exit(1);
  }
}

irtkRealPixel irtkMultiChannelImage::GetValue(int x, int y, int z, int channel)
{
  if ((channel >= 0) && (channel < (int)_images.size())) return _images[channel].Get(x,y,z);
  else {
    cerr << "Channel identificator " << channel <<" out of range." <<endl;
    exit(1);
  }
}


void irtkMultiChannelImage::SetValue(int channel, irtkRealPixel value)
{
  if ((channel >= 0) && (channel < (int)_images.size())) *_pointers[channel] = value;
  else {
    cerr << "Channel identificator out of range." <<endl;
    exit(1);
  }
}

void irtkMultiChannelImage::SetValue(int x, int y, int z, int channel, irtkRealPixel value)
{
  if ((channel >= 0) && (channel < (int)_images.size()))  _images[channel].Put( x, y, z, value);
  else {
    cerr << "Channel identificator " << channel <<" out of range." <<endl;
    exit(1);
  }
}

int irtkMultiChannelImage::GetNumberOfVoxels()
{
  return _number_of_voxels;
}

int irtkMultiChannelImage::GetX()
{
  return _images[0].GetX();
}

int irtkMultiChannelImage::GetY()
{
  return _images[0].GetY();
}

int irtkMultiChannelImage::GetZ()
{
  return _images[0].GetZ();
}


void irtkMultiChannelImage::Write(int channel, const char *file_name)
{
  if ((channel >= 0) && (channel < (int)_images.size())) _images[channel].Write(file_name);
  else {
    cerr << "Channel identificator out of range." <<endl;
    exit(1);
  }
}

int irtkMultiChannelImage::GetNumberOfChannels()
{
  return _images.size();
}

void irtkMultiChannelImage::GetMinMax(int channel, int &min, int &max)
{
  irtkRealPixel mn, mx;
  _images[channel].GetMinMax(&mn, &mx);
  min = (int) mn;
  max = (int) mx;
}


irtkRealImage irtkMultiChannelImage::Average()
{
  int i,j;
  irtkRealPixel *ptr;
  irtkRealImage res(_images[0]);
  double temp;

  cerr << "Averaging" << "...";

  ptr = res.GetPointerToVoxels();
  First();

  for (i=0; i < _number_of_voxels; i++) {
    temp = *ptr;
    if (GetValue(0)!=_padding) {
      for (j=1; j < (int)_images.size(); j++) temp += GetValue(j);
      *ptr = temp / _images.size();
    } else *ptr = _padding;

    ptr++;
    Next();
  }

  cerr << "done." << endl;
  return res;
}

irtkRealImage irtkMultiChannelImage::Add(bool quiet)
{
  int i,j;
  irtkRealPixel *ptr;
  irtkRealImage res(_images[0]);
  double temp;
  bool ok;

  if (!quiet)
     cerr << "Adding" << "...";

  ptr = res.GetPointerToVoxels();
  First();

  for (i=0; i < _number_of_voxels; i++) {
    temp = *ptr;
    ok = true;
    if (GetValue(0)!=_padding) {
      for (j=1; j < (int)_images.size(); j++) temp += GetValue(j);
      *ptr = temp;
    } else *ptr = _padding;

    ptr++;
    Next();
  }

  if (!quiet)
    cerr << "done." << endl;
  return res;
}


irtkRealImage irtkMultiChannelImage::Multiply(bool quiet)
{
  int i,j;
  irtkRealPixel *ptr;
  irtkRealImage res(_images[0]);
  double temp;

  if (!quiet)
    cerr << "Multiplying" << "...";

  ptr = res.GetPointerToVoxels();
  First();

  for (i=0; i < _number_of_voxels; i++) {
    temp = *ptr;
    if (GetValue(0)!=_padding) {
      for (j=1; j < (int)_images.size(); j++) temp *= GetValue(j);
      *ptr = temp;
    } else *ptr = _padding;

    ptr++;
    Next();
  }

  if (!quiet)
    cerr << "done." << endl;
  return res;
}

double irtkMultiChannelImage::ComputeVolume(int channel, int label)
{
  //todo: to find out a real volume of 1 voxel

  int i;
  double volume=0, dx, dy, dz, voxel;
  if ((channel < 0) && (channel >=(int) _images.size())) {
    cerr << "Channel identificator out of range." <<endl;
    exit(1);
  }

  First();

  for (i = 0; i < _number_of_voxels; i++) {
    if (GetValue(channel)==label) volume ++;
    Next();
  }
  _images[channel].GetPixelSize(&dx, &dy, &dz);
  cerr << "Voxel size = " << dx << " " << dy << " " << dz <<endl;
  voxel = dx * dy * dz/1000;
  cerr << "Voxel volume = " << voxel <<endl;
  cerr << "Volume = " << volume << "voxels" <<endl;
  cerr << "Volume = " << volume*voxel << " cm^3" <<endl;
  return volume*voxel;
}

irtkRealImage irtkMultiChannelImage::Subtract(bool quiet)
{
  int i,j;
  irtkRealPixel *ptr;
  irtkRealImage res(_images[0]);
  double temp;

  if (!quiet)
    cerr << "Subtracting" << "...";

  ptr = res.GetPointerToVoxels();
  First();

  for (i=0; i < _number_of_voxels; i++) {
    temp = *ptr;
    if (GetValue(0)!=_padding) {
      for (j=1; j < (int)_images.size(); j++) temp -= GetValue(j);
      *ptr = temp;
    } else *ptr = _padding;

    ptr++;
    Next();
  }

  if (!quiet)
    cerr << "done." << endl;
  return res;
}

irtkRealImage irtkMultiChannelImage::Divide(double scale, bool quiet)
{
  int i,j;
  irtkRealPixel *ptr;
  irtkRealImage res(_images[0]);
  double temp;
  
  if (!quiet)
    cerr << "Dividing" << "...";

  ptr = res.GetPointerToVoxels();
  First();

  for (i=0; i < _number_of_voxels; i++) {
    temp = *ptr;
    if (GetValue(0)!=_padding) {
      for (j=1; j < (int)_images.size(); j++)
        if (GetValue(j)!=0) temp /= GetValue(j);
      *ptr = temp*scale;
    } else *ptr = _padding;

    ptr++;
    Next();
  }

  if (!quiet)
   cerr << "done." << endl;
  return res;
}

void irtkMultiChannelImage::Log(int channel, double scale, bool quiet)
{
  int i;

  if (!quiet)
    cerr << "Log" << "...";

  First();
  for (i=0; i < _number_of_voxels; i++) {
    if ((GetValue(channel)>0)&&(GetValue(channel)!=_padding)) SetValue(channel, scale*log(GetValue(channel)));
    else SetValue(channel, _padding);
    Next();
  }

  if (!quiet)
    cerr << "done." << endl;
}

void irtkMultiChannelImage::Exp(int channel, double scale, bool quiet)
{
  int i;

  if (!quiet)
    cerr << "Exp" << "...";

  First();
  for (i=0; i < _number_of_voxels; i++) {
    if (GetValue(channel)!=_padding) SetValue(channel, exp(GetValue(channel)/scale));
    Next();
  }

  if (!quiet)
    cerr << "done." << endl;
}

irtkRealImage irtkMultiChannelImage::Max()
{
  int i,j;
  irtkRealPixel *ptr;
  irtkRealImage maxImg(_images[0]);
  double temp;

  cerr << "Max" << "...";
  ptr = maxImg.GetPointerToVoxels();
  First();

  for (i=0; i < _number_of_voxels; i++) {
    temp = GetValue(0);
    if (GetValue(0)!=_padding) {
      for (j=1; j < (int)_images.size(); j++) if (GetValue(j) >= temp) temp = GetValue(j);
      *ptr = temp;
    } else *ptr = _padding;

    ptr++;
    Next();
  }

  cerr << "done." << endl;
  return maxImg;
}

irtkRealImage irtkMultiChannelImage::CreateMask()
{
  int i,j;
  irtkRealPixel *ptr;
  double temp;

  cerr << "Creating mask" << "...";
  _mask = _images[0];
  ptr = _mask.GetPointerToVoxels();
  First();

  for (i=0; i < _number_of_voxels; i++) {
    temp = 1;
    for (j=0; j < (int)_images.size(); j++) if (GetValue(j) == _padding) temp = 0;
    *ptr = temp;
    ptr++;
    Next();
  }

  cerr << "done." << endl;
  return _mask;
}

irtkRealImage irtkMultiChannelImage::ExtractLabel(int channel, int label)
{
  int i,j;
  irtkRealPixel *ptr;
  irtkRealImage labels;
  double temp;

  cerr << "Extracting label "<<label << " ...";
  labels = _images[channel];
  ptr = labels.GetPointerToVoxels();
 // First();
  for (i=0; i < _number_of_voxels; i++) {
    temp=0;
    if (*ptr == label) temp = 1;
    *ptr = temp;
    ptr++;
    //Next();
  }

  cerr << "done." << endl;
  return labels;
}

void irtkMultiChannelImage::Brainmask(bool quiet)
{
  int i,j;
  irtkRealPixel *ptr;
  double temp;

  if (!quiet) 
    cerr << "brainmasking all images" << "...";
  ptr = _mask.GetPointerToVoxels();
  First();

  for (i=0; i < _number_of_voxels; i++) {
    temp = 1;
    for (j=0; j < (int)_images.size(); j++) if (*ptr == 0) SetValue(j,_padding);
    ptr++;
    Next();
  }

  if (!quiet) 
    cerr << "done." << endl;

}

irtkRealImage irtkMultiChannelImage::Brainmask(int channel)
{
  int i,j;
  irtkRealPixel *ptr;
  double temp;

  //if (!quiet) 
    cout << "brainmasking" << "...";
  ptr = _mask.GetPointerToVoxels();
  First();

  for (i=0; i < _number_of_voxels; i++) {
    temp = 1;
    if (*ptr == 0) SetValue(channel,_padding);
    ptr++;
    Next();
  }

  //if (!quiet) 
    cout << "done." << endl;
  return _images[channel];

}



void irtkMultiChannelImage::SetPadding(int padding)
{
  _padding=padding;
}

void irtkMultiChannelImage::HistogramMatching(int tchannel, int rchannel, double sigma)
{
  //irtkRealImage im,r;
  //im = _images[tchannel];
  //r = _images[rchannel];

  int i,n;

  if (sigma>0) {
    cerr<<"Blurring."<<endl;
    irtkGaussianBlurringWithPadding<irtkRealPixel> gaussianBlurring(sigma, _padding);
    gaussianBlurring.SetInput (&_images[tchannel]);
    gaussianBlurring.SetOutput(&_images[tchannel]);
    gaussianBlurring.Run();

    gaussianBlurring.SetInput (&_images[rchannel]);
    gaussianBlurring.SetOutput(&_images[rchannel]);
    gaussianBlurring.Run();

    Brainmask();
    Write(tchannel,"blurred_image.nii.gz");
    Write(rchannel,"blurred_reference.nii.gz");
  }

  cerr<<"Matching intensities."<<endl;

  irtkRealPixel rmin, rmax, tmin, tmax;
  _images[tchannel].GetMinMax(&tmin, &tmax);
  _images[rchannel].GetMinMax(&rmin, &rmax);

  irtkHistogram href((int)rmax);
  irtkHistogram htar((int)tmax);
  double *v = new double [round(rmax)];

  int htmax = (int)tmax;
  int hrmax = (int)rmax;


  href.PutMin(0);
  href.PutMax(rmax);
  htar.PutMin(0);
  htar.PutMax(tmax);

  n=GetNumberOfVoxels();

  cerr<<"Creating histograms."<<endl;
  First();
  for (i=0; i<n; i++) {
    if ( GetValue(tchannel)!=_padding ) htar.AddSample(GetValue(tchannel));
    if ( GetValue(rchannel)!=_padding ) href.AddSample(GetValue(rchannel));
    Next();
  }

  cerr<<"Calculating transformation."<<endl;
  double R=0, T=0;
  int k=0;
  double nt=htar.NumberOfSamples(), nr=href.NumberOfSamples();

  for (i=0; i<hrmax; i++) {
    R+=href(i);
    while ((R/nr)>=(T/nt)) {
      T+=htar(k);
      if (k<htmax-1) k++;
    }
    v[i]=htar.BinToVal(k);
  }

  cerr<<endl;

  /*  irtkGreyImage graph(htmax+1,hrmax+1,1);
    for(i=0; i<hrmax; i++)
    {
      graph.Put(v[i],i,0,1);
    }
    graph.Write("graph.nii.gz");
  */

  cerr<<"adjusting intensities of reference to the image."<<endl;

  First();
  for (i=0; i<n; i++) {
    if (GetValue(rchannel)!=_padding) {
      SetValue(rchannel,v[href.ValToBin(GetValue(rchannel))]);
    }
    Next();
  }

  delete [] v;

  //Write(rchannel,"rm.nii.gz");

}

void irtkMultiChannelImage::HistogramEqualization(int channel)
{
  int i,n;
  cerr<<"Equalizing histogram."<<endl;

  irtkRealPixel mn, mx;
  _images[channel].GetMinMax(&mn, &mx);

  irtkHistogram h((int)mx);
  double *v = new double [round(mx)];

  int hmx = (int)mx;

  h.PutMin(0);
  h.PutMax(hmx);

  n=GetNumberOfVoxels();

  cerr<<"Creating histograms."<<endl;
  First();
  for (i=0; i<n; i++) {
    if ( GetValue(channel)!=_padding ) h.AddSample(GetValue(channel));
    Next();
  }

  cerr<<"Calculating transformation."<<endl;
  double R=0, T=0;
  int k=0;
  double ns=h.NumberOfSamples();

  for (i=0; i<hmx; i++) {
    R+=h(i);
    while ((R/ns)>=(T/ns)) {
      T+=ns/mx;
      if (k<hmx-1) k++;
    }
    v[i]=k;
    cerr<<k<<" ";
  }

  cerr<<endl;

  cerr<<"Equalizing intensities of the image."<<endl;

  First();
  for (i=0; i<n; i++) {
    if (GetValue(channel)!=_padding) {
      int value = (int)GetValue(channel);
      SetValue(channel,v[value]);
    }
    Next();
  }

  delete [] v;

  //Write(channel,"equalized.nii.gz");

}

void irtkMultiChannelImage::AdjustMean(int tchannel, int rchannel)
{
  int n=GetNumberOfVoxels();
  int i;
  double sum=0;
  int count=0;

  cerr<<"Calculating mean difference."<<endl;
  First();
  for (i=0; i<n; i++) {
    if (( GetValue(tchannel)!=_padding )&&(GetValue(rchannel)!=_padding)) {
      sum += (GetValue(tchannel) - GetValue(rchannel));
      count++;
    }
    Next();
  }

  double mean = sum/count;
  cerr<<"Sum: "<<sum<<endl;
  cerr<<"Count: "<<count<<endl;
  cerr<<"Mean difference: "<<mean<<endl;

  cerr<<"Adjusting image to reference."<<endl;
  First();
  for (i=0; i<n; i++) {
    if ( GetValue(tchannel)!=_padding ) {
      SetValue(tchannel,GetValue(tchannel) - mean);
    }
    Next();
  }

}


double irtkMultiChannelImage::Average(int channel, bool quiet)
{
  int n=GetNumberOfVoxels();
  int i;
  double sum=0;
  int count=0;

  if (!quiet)
    cerr<<"Calculating intensity average."<<endl;
  First();
  for (i=0; i<n; i++) {
    if ( GetValue(channel)!=_padding ) 
    {
      sum += GetValue(channel);
      count++;
    }
    Next();
  }

  double mean = sum/count;
  //cerr<<"Sum: "<<sum<<endl;
  //cerr<<"Count: "<<count<<endl;
  
  if (!quiet)
    cerr<<"Average: "<<mean<<endl;

  return mean;
}


double irtkMultiChannelImage::StDev(int channel, bool quiet)
{
  int n=GetNumberOfVoxels();
  int i;
  double sum=0;
  int count=0;
  double mean = Average(channel,true);
  
  if (!quiet)
    cerr<<"Calculating intensity variance."<<endl;
  First();
  for (i=0; i<n; i++) {
    if ( GetValue(channel)!=_padding ) 
    {
      sum += (GetValue(channel)-mean)*(GetValue(channel)-mean);
      count++;
    }
    Next();
  }

  double stdev = sqrt(sum/count);
  //cerr<<"Sum: "<<sum<<endl;
  //cerr<<"Count: "<<count<<endl;
  
  if (!quiet)
    cerr<<"StDev: "<<stdev<<endl;

  return stdev;
}

void irtkMultiChannelImage::Sqrt(int channel)
{
  int i;

  cerr << "Log" << "...";

  First();
  for (i=0; i < _number_of_voxels; i++) {
    if ((GetValue(channel)>0)&&(GetValue(channel)!=_padding)) SetValue(channel, sqrt(GetValue(channel)));
    else SetValue(channel, _padding);
    Next();
  }

  cerr << "done." << endl;
}



