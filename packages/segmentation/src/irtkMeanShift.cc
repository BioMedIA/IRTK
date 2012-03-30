
// queue::push/pop

#include <irtkMeanShift.h>


irtkMeanShift::irtkMeanShift(irtkGreyImage& image, int padding, int nBins)
{
  _image = image;
  _orig_image=image;
  _nBins=nBins;
  _padding=padding;
  _brain=NULL;
  _density = new double[_nBins];
  for(int i=0;i<_nBins;i++) _density[i]=0;
  _bg=-1;
  _gm=-1;
  _wm=-1;
  _bin_width=1;

}

irtkMeanShift::~irtkMeanShift()
{
  delete[] _density;
}

irtkRealImage irtkMeanShift::ReturnMask()
{
  irtkRealImage mask = _image;
  irtkRealPixel *ptr = mask.GetPointerToVoxels();
  int n = _image.GetNumberOfVoxels();
  for(int i=0;i<n;i++)
  {
    if(*ptr!=_padding)
      *ptr=1;
    else
      *ptr=0;
    ptr++;
  }
  return mask;
}


void irtkMeanShift::SetOutput(irtkGreyImage *output)
{
  _output = output;
}

double irtkMeanShift::ValueToBin(double value)
{
  return _nBins*(value-_imin)/(_imax+1-_imin);
}

double irtkMeanShift::BinToValue(int bin)
{
  return _imin+bin*(_imax+1-_imin)/_nBins;
}

double irtkMeanShift::GenerateDensity(double cut_off)
{
  int i,j;
  irtkGreyPixel *ptr=_image.GetPointerToVoxels();
  int n = _image.GetNumberOfVoxels();
  int voxels=0;

  _image.GetMinMax(&_imin,&_imax);
  cout<<_imin<<" "<<_imax<<endl;
  _bin_width = (_imax-_imin+1)/_nBins;

  cout<<"generating density..."<<endl;
  for(i=0;i<n;i++)
  {
    if(*ptr>_padding)
    {
      j = ValueToBin(*ptr);
      _density[j]++;
      voxels++;
    }
    ptr++;
  }
  cout<<"done"<<endl;


  cout<<"Cutting off 2% of highest intensities"<<endl;
  double sum=0,bs=0;
  for(i=0; i<_nBins; i++) sum+=_density[i];
  i=_nBins-1;
  while(bs<cut_off*sum)
  {
    bs+=_density[i];
    i--;
  }

  _limit = BinToValue(i);

  cout<<"limit="<<_limit<<endl<<endl;
  return _limit;
}


void irtkMeanShift::AddPoint(int x, int y, int z)
{
  if((x>=0)&&(y>=0)&&(z>=0)&&(x<_image.GetX())&&(y<_image.GetY())&&(z<_image.GetZ()))
  {
    bool mask=true;
    irtkPoint p;
    if(_brain!=NULL) if(_brain->Get(x,y,z)==1) mask=false;

    if((_map.Get(x,y,z)==1)&&(_image.Get(x,y,z)<_treshold)&&(mask))
    {
      p._x=x;
      p._y=y;
      p._z=z;

      _q.push(p);

      _map.Put(x,y,z,0);
    }
  }
}

void irtkMeanShift::AddPoint(int x, int y, int z, int label)
{
  if((x>=0)&&(y>=0)&&(z>=0)&&(x<_image.GetX())&&(y<_image.GetY())&&(z<_image.GetZ()))
  {
    irtkPoint p;

    if((_map.Get(x,y,z)==0)&&(_image.Get(x,y,z)==label))
    {
      p._x=x;
      p._y=y;
      p._z=z;

      _q.push(p);

      _map.Put(x,y,z,1);
      _clusterSize++;
    }
  }
}

double irtkMeanShift::msh(double y, double h)
{
  double points, sum, y0;
  irtkGreyPixel* ptr;
  int n = _image.GetNumberOfVoxels();

  //cout<<"msh: position = "<<y<<", bandwidth = "<<h<<endl;

  if(h<=_bin_width) return y;

  do
  {
    ptr = _image.GetPointerToVoxels();
    sum=0;
    points=0;

    for(int i = 0; i < n; i++)
    {
      if (*ptr > _padding)
      {
        if ( abs(*ptr-y) <= h)
        {
          points++;
	  sum += *ptr;
        }
      }
      ptr++;
    }
    y0=y;
    y=sum/points;
  }
  while (abs(y-y0) >= 1);

  return y;
}

double irtkMeanShift::findGMvar()
{
  double k;
  int mean = ValueToBin(_gm);
  double sigma=0;

  for (int i=mean; i>=0; i--)
  {
    k=sqrt(2*log(_density[mean]/(double)_density[i]));
    sigma=abs(BinToValue(i)-_gm)/k;
    cout<<"Pos "<<BinToValue(i)<<": sigma = "<<sigma<<endl;
  }

  return sigma;
}

double irtkMeanShift::findMax(double tr1, double tr2)
{
  double mx=0, pos=-1;
  int j1,j2,i0,i;
  //j1 = nBins*(tr1-imin)/(imax+1-imin);
  //j2 = nBins*(tr2-imin)/(imax+1-imin);
  j1=ValueToBin(tr1);
  j2=ValueToBin(tr2);

  //cout<<imin<<" "<<imax<<endl;
  //cout<<j1<<" "<<j2<<endl;

  i0 = j1;
  for (i=j1; i<=j2; i++)
  {
    //cout<<_density[i]<<" ";
    if(mx<_density[i])
    {
     i0=i;
     mx=_density[i];
     //cout<<"density["<<i<<"]="<<_density[i]<<", mx="<<mx<<endl;
    }
  }
  //cout<<i0<<endl;
  pos=BinToValue(i0);//imin+i0*(imax+1-imin)/nBins;
  //cout<<pos<<endl;

  cout<<"Maximum for interval <"<<tr1<<", "<<tr2<<"> is "<<pos<<endl;
  return pos;
}

double irtkMeanShift::findMin(double tr1, double tr2)
{
  double mx=100000000, pos=-1;
  int j1,j2,i0,i;
  //j1 = nBins*(tr1-imin)/(imax+1-imin);
  //j2 = nBins*(tr2-imin)/(imax+1-imin);

  //cout<<imin<<" "<<imax<<endl;
  //cout<<j1<<" "<<j2<<endl;

  j1=ValueToBin(tr1);
  j2=ValueToBin(tr2);
  i0=j1;

  for (i=j1; i<=j2; i++)
  {
    if(mx>_density[i])
    {
     i0=i;
     mx=_density[i];
     //cout<<"density["<<i<<"]="<<density[i]<<", mx="<<mx<<endl;
    }
  }
  pos=BinToValue(i0);//imin+i0*(imax+1-imin)/nBins;
  //pos=imin+i0*(imax+1-imin)/nBins;

  cout<<"Minimum for interval <"<<tr1<<", "<<tr2<<"> is "<<pos<<endl;
  return pos;
}

double irtkMeanShift::split(double pos1, double pos2, double bw, double h1, double h2)
{
  double tr=0.02*_nBins;
  double m1, m2, m3;
  bool change = true;

  cout<<"split: "<<pos1<<" "<<pos2<<" "<<bw<<" "<<h1<<" "<<h2<<endl;

  while(((pos2-pos1) > _bin_width)&&(change)&&( bw>_bin_width ))
  {
    change=false;
    cout<<"pos1="<<pos1<<" pos2="<<pos2<<endl;
    tr=_bin_width;
    if (tr<1) tr=1;
    m1=msh(pos1,bw);
    m2=msh(pos2,bw);
    cout<<"m1="<<m1<<", m2="<<m2<<endl;
    
    if((m2-m1)< tr)
    {
      h2=bw;
      bw=(h1+h2)/2;
      if(bw<=h2-1)
      {
        cout<<"reducing bandwidth: bw="<<bw<<endl;
        cout<<"h1="<<h1<<", h2="<<h2<<endl;
        change=true;
      }
    }
    else
    {
      cout<<"changing positions:"<<endl;
      m3=msh((pos1+pos2)/2, bw);
      cout<<"m3="<<m3<<endl;
      if((m3-m1)<tr)
      {
        pos1=(pos1+pos2)/2;
        change=true;
        cout<<"pos1="<<pos1<<endl;
      }
      else
      {
        if((m2-m3)<tr)
	{
	  pos2=(pos1+pos2)/2;
	  change=true;
        cout<<"pos2="<<pos2<<endl;
	}
	else
	{
	  h1=bw;
	  bw=(h1+h2)/2;
	  if(bw>=h1+1)
	  {
            cout<<"increasing bandwidth: bw="<<bw<<endl;
            cout<<"h1="<<h1<<", h2="<<h2<<endl;

	    change=true;
	  }
          else 
          {
            //cout<<"bandwidth fixed:"<<bw<<endl;
            cout<<"change=false, exiting split."<<endl;
            
          }
	}
      }
      ///////end changing positions
    }
  }

  cout<<"treshold="<<pos1<<endl;

  return pos1;
}


void irtkMeanShift::SetTreshold(double treshold)
{
  cout<<"treshold = "<<treshold<<endl;

  _treshold = treshold;
}

void irtkMeanShift::SetTreshold()
{

  double pos1=0, pos2=_limit;
  double bw=2*_nBins;
  double h1=0, h2=bw;

  cout<<"calculating treshold ... ";

  _limit1 = split(pos1, pos2, bw, h1, h2);
  _bg=findMax(0,_limit1);
  _gm=findMax(_limit1,_limit);
  _treshold = findMin(_bg,_gm);
  //_treshold = _limit1;
}


int irtkMeanShift::Lcc(int label, bool add_second)
{
 int i,j,k;
 int lcc_size = 0;
 int lcc2_size = 0;
 int lcc_x = 0, lcc_y = 0, lcc_z = 0, lcc2_x = 0, lcc2_y = 0, lcc2_z = 0;

 //cout<<"Finding Lcc"<<endl;
  _map=_image;
  //_image.Write("image.nii.gz");

  irtkGreyPixel* ptr=_map.GetPointerToVoxels();
  int n = _image.GetNumberOfVoxels();
  for(i=0;i<n;i++)
  {
    *ptr=0;
    ptr++;
  }

  for (i=0; i<_image.GetX();i++)
    for (j=0; j<_image.GetY();j++)
      for (k=0; k<_image.GetZ();k++)
      {
        _clusterSize = 0;
        if ((_map.Get(i,j,k)==0)&&(_image.Get(i,j,k)==label)) Grow(i,j,k,label);
        if (_clusterSize > lcc_size)
        {
          lcc2_size = lcc_size;
          lcc2_x = lcc_x;
          lcc2_y = lcc_y;
          lcc2_z = lcc_z;

          lcc_size = _clusterSize;
          lcc_x = i;
          lcc_y = j;
          lcc_z = k;
        }
      }

  ptr=_map.GetPointerToVoxels();
  n = _image.GetNumberOfVoxels();
  for(i=0;i<n;i++)
  {
    *ptr=0;
    ptr++;
  }
  
  Grow(lcc_x,lcc_y,lcc_z,label);

  if((add_second)&&(lcc2_size > 0.5*lcc_size))
  {
    cout<<"Adding second largest cluster too. ";
    Grow(lcc2_x,lcc2_y,lcc2_z,label);
  }
  //_map.Write("lcc.nii.gz");
  *_output = _map;

  return lcc_size;
}

int irtkMeanShift::LccS(int label, double treshold)
{
 int i,j,k;
 int lcc_size = 0;
 int lcc_x, lcc_y, lcc_z;
 queue<irtkPoint> seed;
 queue<int> size;

 //cout<<"Finding Lcc and all cluster of 70% of the size of Lcc"<<endl;
  _map=_image;
  //_image.Write("image.nii.gz");

  irtkGreyPixel* ptr=_map.GetPointerToVoxels();
  int n = _image.GetNumberOfVoxels();
  for(i=0;i<n;i++)
  {
    *ptr=0;
    ptr++;
  }

  for (i=0; i<_image.GetX();i++)
    for (j=0; j<_image.GetY();j++)
      for (k=0; k<_image.GetZ();k++)
      {
        _clusterSize = 0;
        if ((_map.Get(i,j,k)==0)&&(_image.Get(i,j,k)==label)) Grow(i,j,k,label);
        if (_clusterSize > lcc_size)
        {
          lcc_size = _clusterSize;
          lcc_x = i;
          lcc_y = j;
          lcc_z = k;
        }

        if (_clusterSize > 0)
        {
          size.push(_clusterSize);
          irtkPoint p(i,j,k);
          seed.push(p);
        }
      }

  ptr=_map.GetPointerToVoxels();
  n = _image.GetNumberOfVoxels();
  for(i=0;i<n;i++)
  {
    *ptr=0;
    ptr++;
  }
  
  while (!size.empty())
  {
    if ( size.front()> treshold*lcc_size ) Grow(seed.front()._x,seed.front()._y,seed.front()._z,label);
    seed.pop();
    size.pop();
  }
  //_map.Write("lcc.nii.gz");
  *_output = _map;

  return lcc_size;
}


void irtkMeanShift::Grow(int x, int y, int z, int label)
{
  AddPoint(x,y,z,label);

  while (!_q.empty())
  {
    x=_q.front()._x;
    y=_q.front()._y;
    z=_q.front()._z;

    _q.pop();

    AddPoint(x-1,y,z,label);
    AddPoint(x,y-1,z,label);
    AddPoint(x,y,z-1,label);
    AddPoint(x+1,y,z,label);
    AddPoint(x,y+1,z,label);
    AddPoint(x,y,z+1,label);
  }

}


void irtkMeanShift::RegionGrowing()
{
 int i; 
 cout<<"Removing background"<<endl;
  _map=_image;

  irtkGreyPixel* ptr=_map.GetPointerToVoxels();
  int n = _image.GetNumberOfVoxels();

  for(i=0;i<n;i++)
  {
    *ptr=1;
    ptr++;
  }

  AddPoint(0,0,0);
  AddPoint(_image.GetX()-1,0,0);
  AddPoint(0,_image.GetY()-1,0);
  AddPoint(0,0,_image.GetZ()-1);
  AddPoint(_image.GetX()-1,_image.GetY()-1,0);
  AddPoint(_image.GetX()-1,0,_image.GetZ()-1);
  AddPoint(0,_image.GetY()-1,_image.GetZ()-1);
  AddPoint(_image.GetX()-1,_image.GetY()-1,_image.GetZ()-1);

  int x,y,z;

  while (!_q.empty())
  {
    x=_q.front()._x;
    y=_q.front()._y;
    z=_q.front()._z;

    _q.pop();

    AddPoint(x-1,y,z);
    AddPoint(x,y-1,z);
    AddPoint(x,y,z-1);
    AddPoint(x+1,y,z);
    AddPoint(x,y+1,z);
    AddPoint(x,y,z+1);
  }
}


void irtkMeanShift::RemoveBackground()
{
 int i; 
 cout<<"Removing background"<<endl;
  _map=_image;

  irtkGreyPixel* ptr=_map.GetPointerToVoxels();
  int n = _image.GetNumberOfVoxels();

  for(i=0;i<n;i++)
  {
    *ptr=1;
    ptr++;
  }

  AddPoint(0,0,0);
  AddPoint(_image.GetX()-1,0,0);
  AddPoint(0,_image.GetY()-1,0);
  AddPoint(0,0,_image.GetZ()-1);
  AddPoint(_image.GetX()-1,_image.GetY()-1,0);
  AddPoint(_image.GetX()-1,0,_image.GetZ()-1);
  AddPoint(0,_image.GetY()-1,_image.GetZ()-1);
  AddPoint(_image.GetX()-1,_image.GetY()-1,_image.GetZ()-1);

  int x,y,z;

  while (!_q.empty())
  {
    x=_q.front()._x;
    y=_q.front()._y;
    z=_q.front()._z;

    _q.pop();

    AddPoint(x-1,y,z);
    AddPoint(x,y-1,z);
    AddPoint(x,y,z-1);
    AddPoint(x+1,y,z);
    AddPoint(x,y+1,z);
    AddPoint(x,y,z+1);
  }

  cout<< "dilating and eroding ... ";
  _brain = new irtkGreyImage(_map);

  int iterations = 3;
  irtkDilation<irtkGreyPixel> dilation;
  dilation.SetInput(_brain);
  dilation.SetOutput(_brain);
  for (i = 0; i < iterations; i++) dilation.Run();

  irtkErosion<irtkGreyPixel> erosion;
  erosion.SetInput(_brain);
  erosion.SetOutput(_brain);
  for (i = 0; i < iterations; i++) erosion.Run();

  cout<<"recalculating ... ";

  ptr=_map.GetPointerToVoxels();
  for(i=0;i<n;i++)
  {
    *ptr=1;
    ptr++;
  }

  AddPoint(0,0,0);
  AddPoint(_image.GetX()-1,0,0);
  AddPoint(0,_image.GetY()-1,0);
  AddPoint(0,0,_image.GetZ()-1);
  AddPoint(_image.GetX()-1,_image.GetY()-1,0);
  AddPoint(_image.GetX()-1,0,_image.GetZ()-1);
  AddPoint(0,_image.GetY()-1,_image.GetZ()-1);
  AddPoint(_image.GetX()-1,_image.GetY()-1,_image.GetZ()-1);

  while (!_q.empty())
  {
    x=_q.front()._x;
    y=_q.front()._y;
    z=_q.front()._z;

    _q.pop();

    AddPoint(x-1,y,z);
    AddPoint(x,y-1,z);
    AddPoint(x,y,z-1);
    AddPoint(x+1,y,z);
    AddPoint(x,y+1,z);
    AddPoint(x,y,z+1);
  }

   cout<<"eroding ... ";

  *_brain = _map;

  erosion.SetInput(_brain);
  erosion.SetOutput(_brain);
  for (i = 0; i < iterations; i++) erosion.Run();

  cout<<"final recalculation ...";

  ptr=_map.GetPointerToVoxels();
  for(i=0;i<n;i++)
  {
    *ptr=1;
    ptr++;
  }

  AddPoint(0,0,0);
  AddPoint(_image.GetX()-1,0,0);
  AddPoint(0,_image.GetY()-1,0);
  AddPoint(0,0,_image.GetZ()-1);
  AddPoint(_image.GetX()-1,_image.GetY()-1,0);
  AddPoint(_image.GetX()-1,0,_image.GetZ()-1);
  AddPoint(0,_image.GetY()-1,_image.GetZ()-1);
  AddPoint(_image.GetX()-1,_image.GetY()-1,_image.GetZ()-1);

  while (!_q.empty())
  {
    x=_q.front()._x;
    y=_q.front()._y;
    z=_q.front()._z;

    _q.pop();

    AddPoint(x-1,y,z);
    AddPoint(x,y-1,z);
    AddPoint(x,y,z-1);
    AddPoint(x+1,y,z);
    AddPoint(x,y+1,z);
    AddPoint(x,y,z+1);
  }

  delete _brain;
  cout<<"done."<<endl;

  ptr=_map.GetPointerToVoxels();
  irtkGreyPixel *ptr_b=_image.GetPointerToVoxels();
  for(i=0;i<n;i++)
  {
    if(*ptr==0) *ptr_b = _padding;
    ptr++;
    ptr_b++;
  }


}

void irtkMeanShift::Write(char *output_name)
{
  _image.Write(output_name);
}

void irtkMeanShift::WriteMap(char *output_name)
{
  _map.Write(output_name);
}

void irtkMeanShift::FindWMGMmeans()
{
  cout<<"Finding WM and GM mean"<<endl;
  cout<<"Make sure that background has been removed when calling this method!!!"<<endl;

  _gm=findMax(0, _limit);
  _limit2 = split(_gm, _limit, 2*_nBins, 0, 2*_nBins);
  _wm=findMax(_limit2,_limit);
  _split2=findMin(_gm,_wm);
  cout<<"GM mean = "<<_gm<<endl;
  cout<<"WM mean = "<<_wm<<endl;

/*
  _limit2 = split(0, _limit, 2*_nBins, 0, 2*_nBins);


  _gm=findMax(0, _limit2);
  _wm=findMax(_limit2,_limit);

  if(_bg>-1) _split1=findMin(_bg,_gm);
  _split2=findMin(_gm,_wm);

  cout<<"GM mean = "<<_gm<<endl;
  cout<<"WM mean = "<<_wm<<endl;

  irtkGreyImage gmwm(_image);
  irtkGreyPixel *ptr=gmwm.GetPointerToVoxels();
  int n=gmwm.GetNumberOfVoxels();
  for(int i=0;i<n;i++)
  {
    if(*ptr>_padding)
    {
      if(*ptr < _split2) *ptr=0;
      else *ptr=1;
    }
    ptr++;
  }
  gmwm.Write("gm-wm-split.nii.gz");

*/
}




