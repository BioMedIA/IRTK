/*=========================================================================

  Date: 14.09.2009
  Author: Laurent Risser
  +++ Works very well +++

=========================================================================*/


#include <irtkImageFastFourierTransform.h>

void four1NR(float * data, unsigned long nn, int isign)
{
  unsigned long n,mmax,m,j,istep,i;
  double wtemp,wr,wpr,wpi,wi,theta;
  float tempr,tempi;

  n=nn << 1;
  j=1;
  for (i=1;i<n;i+=2) {
    if (j > i) {
      tempr=data[j]; data[j]=data[i]; data[i]=tempr;
      tempr=data[j+1]; data[j+1]=data[i+1]; data[i+1]=tempr;
    }
    m=n >> 1;
    while (m >= 2 && j > m) {
      j -= m;
      m >>= 1;
    }
    j += m;
  }
  mmax=2;
  while (n > mmax) {
    istep=mmax << 1;
    theta=isign*(6.28318530717959/mmax);
    wtemp=sin(0.5*theta);
    wpr = -2.0*wtemp*wtemp;
    wpi=sin(theta);
    wr=1.0;
    wi=0.0;
    for (m=1;m<mmax;m+=2) {
      for (i=m;i<=n;i+=istep) {
        j=i+mmax;
        tempr=wr*data[j]-wi*data[j+1];
        tempi=wr*data[j+1]+wi*data[j];
        data[j]=data[i]-tempr;
        data[j+1]=data[i+1]-tempi;
        data[i] += tempr;
        data[i+1] += tempi;
      }
      wr=(wtemp=wr)*wpr-wi*wpi+wr;
      wi=wi*wpr+wtemp*wpi+wi;
    }
    mmax=istep;
  }
}



void DirectFFT(irtkGenericImage<float> * RealSignal,irtkGenericImage<float> * ImaginarySignal){
  int SizeX,SizeY,SizeZ;
  float SqrtSizeX,SqrtSizeY,SqrtSizeZ;
  int x,y,z;
  float * dataX;
  float * dataY;
  float * dataZ;
  cout << "yo\n";
  //1) extract the size of the images
  SizeX=RealSignal->GetX();
  SizeY=RealSignal->GetY();
  SizeZ=RealSignal->GetZ();
  
  SqrtSizeX=static_cast<float>(sqrt(static_cast<double>(SizeX)));
  SqrtSizeY=static_cast<float>(sqrt(static_cast<double>(SizeY)));
  SqrtSizeZ=static_cast<float>(sqrt(static_cast<double>(SizeZ)));

  
  //2) perform the fft along x axis
  dataX = new float [SizeX*2+1];
  for (z = 0; z < SizeZ; z++) for (y = 0; y < SizeY; y++){
    for (x = 0; x < SizeX; x++){
      dataX[2*x+1]=RealSignal->Get(x, y, z, 0);
      dataX[2*x+2]=ImaginarySignal->Get(x, y, z, 0);
    }
    four1NR(dataX, (unsigned long)SizeX, 1);
    for (x = 0; x < SizeX; x++){
      RealSignal->Put(x, y, z, 0,dataX[2*x+1]/SqrtSizeX);
      ImaginarySignal->Put(x, y, z, 0,dataX[2*x+2]/SqrtSizeX);
    }
  }
  delete dataX;
  
  //3) perform the fft along y axis
  dataY = new float [SizeY*2+1];
  for (z = 0; z < SizeZ; z++) for (x = 0; x < SizeX; x++){
    for (y = 0; y < SizeY; y++){
      dataY[2*y+1]=RealSignal->Get(x, y, z, 0);
      dataY[2*y+2]=ImaginarySignal->Get(x, y, z, 0);
    }
    four1NR(dataY, (unsigned long)SizeY, 1);
    for (y = 0; y < SizeY; y++){
      RealSignal->Put(x, y, z, 0,dataY[2*y+1]/SqrtSizeY);
      ImaginarySignal->Put(x, y, z, 0,dataY[2*y+2]/SqrtSizeY);
    }
  }
  delete dataY;
  
  
  //4) perform the fft along z axis
  dataZ = new float [SizeZ*2+1];
  for (y = 0; y < SizeY; y++) for (x = 0; x < SizeX; x++){
    for (z = 0; z < SizeZ; z++){
      dataZ[2*z+1]=RealSignal->Get(x, y, z, 0);
      dataZ[2*z+2]=ImaginarySignal->Get(x, y, z, 0);
    }
    four1NR(dataZ, (unsigned long)SizeZ, 1);
    for (z = 0; z < SizeZ; z++){
      RealSignal->Put(x, y, z, 0,dataZ[2*z+1]/SqrtSizeZ);
      ImaginarySignal->Put(x, y, z, 0,dataZ[2*z+2]/SqrtSizeZ);
    }
  }
  delete dataZ;
}


void InverseFFT(irtkGenericImage<float> * RealSignal,irtkGenericImage<float> * ImaginarySignal){
  int SizeX,SizeY,SizeZ;
  float SqrtSizeX,SqrtSizeY,SqrtSizeZ;
  int x,y,z;
  float * dataX;
  float * dataY;
  float * dataZ;
  
  //1) extract the size of the images
  SizeX=RealSignal->GetX();
  SizeY=RealSignal->GetY();
  SizeZ=RealSignal->GetZ();
  
  SqrtSizeX=static_cast<float>(sqrt(static_cast<double>(SizeX)));
  SqrtSizeY=static_cast<float>(sqrt(static_cast<double>(SizeY)));
  SqrtSizeZ=static_cast<float>(sqrt(static_cast<double>(SizeZ)));
  
  
  //2) perform the ifft along z axis
  dataZ = new float [SizeZ*2+1];
  for (y = 0; y < SizeY; y++) for (x = 0; x < SizeX; x++){
    for (z = 0; z < SizeZ; z++){
      dataZ[2*z+1]=RealSignal->Get(x, y, z, 0);
      dataZ[2*z+2]=ImaginarySignal->Get(x, y, z, 0);
    }
    four1NR(dataZ, (unsigned long)SizeZ, -1);
    for (z = 0; z < SizeZ; z++){
      RealSignal->Put(x, y, z, 0,dataZ[2*z+1]/SqrtSizeZ);
      ImaginarySignal->Put(x, y, z, 0,dataZ[2*z+2]/SqrtSizeZ);
    }
  }
  delete dataZ;
  
  //3) perform the ifft along y axis
  dataY = new float [SizeY*2+1];
  for (z = 0; z < SizeZ; z++) for (x = 0; x < SizeX; x++){
    for (y = 0; y < SizeY; y++){
      dataY[2*y+1]=RealSignal->Get(x, y, z, 0);
      dataY[2*y+2]=ImaginarySignal->Get(x, y, z, 0);
    }
    four1NR(dataY, (unsigned long)SizeY, -1);
    for (y = 0; y < SizeY; y++){
      RealSignal->Put(x, y, z, 0,dataY[2*y+1]/SqrtSizeY);
      ImaginarySignal->Put(x, y, z, 0,dataY[2*y+2]/SqrtSizeY);
    }
  }
  delete dataY;
  
  //4) perform the ifft along x axis
  dataX = new float [SizeX*2+1];
  for (z = 0; z < SizeZ; z++) for (y = 0; y < SizeY; y++){
    for (x = 0; x < SizeX; x++){
      dataX[2*x+1]=RealSignal->Get(x, y, z, 0);
      dataX[2*x+2]=ImaginarySignal->Get(x, y, z, 0);
    }
    four1NR(dataX, (unsigned long)SizeX, -1);
    for (x = 0; x < SizeX; x++){
      RealSignal->Put(x, y, z, 0,dataX[2*x+1]/SqrtSizeX);
      ImaginarySignal->Put(x, y, z, 0,dataX[2*x+2]/SqrtSizeX);
    }
  }
  delete dataX;
}

void ConvolutionInFourier(irtkGenericImage<float> * RealPartSignal,irtkGenericImage<float> * ImaginaryPartSignal,irtkGenericImage<float> * RealPartFilter,irtkGenericImage<float> * ImaginaryPartFilter){
  int x,y,z;
  float a,b,c,d;
  float CoefMult;
  
  
  //1) FFT
  DirectFFT(RealPartSignal,ImaginaryPartSignal);
  DirectFFT(RealPartFilter,ImaginaryPartFilter);
  
  //2) filtering in Fourier spaces
  CoefMult=(float)(sqrt(float(RealPartSignal->GetX())*sqrt((float)RealPartSignal->GetY())*sqrt((float)RealPartSignal->GetZ())));
  
  for (z = 0; z < RealPartSignal->GetZ(); z++) for (y = 0; y < RealPartSignal->GetY(); y++) for (x = 0; x < RealPartSignal->GetX(); x++){
    a=RealPartSignal->Get(x, y, z, 0);
    b=ImaginaryPartSignal->Get(x, y, z, 0);
    c=RealPartFilter->Get(x, y, z, 0)*CoefMult;
    d=ImaginaryPartFilter->Get(x, y, z, 0)*CoefMult;
    
    RealPartSignal->Put(x, y, z, 0, a*c-b*d);
    ImaginaryPartSignal->Put(x, y, z, 0, c*b+a*d);
  }
  
  //3) IFFT
  InverseFFT(RealPartSignal,ImaginaryPartSignal);
  InverseFFT(RealPartFilter,ImaginaryPartFilter);
}


void DeconvolutionInFourier(irtkGenericImage<float> * RealPartSignal,irtkGenericImage<float> * ImaginaryPartSignal,irtkGenericImage<float> * RealPartFilter,irtkGenericImage<float> * ImaginaryPartFilter){
  int x,y,z;
  float a,b,c,d;
  float CoefMult;
  
  
  //1) FFT
  DirectFFT(RealPartSignal,ImaginaryPartSignal);
  DirectFFT(RealPartFilter,ImaginaryPartFilter);
  
  //2) filtering in Fourier spaces
  CoefMult=(float)(sqrt((float)RealPartSignal->GetX())*sqrt((float)RealPartSignal->GetY())*sqrt((float)RealPartSignal->GetZ()));
  
  for (z = 0; z < RealPartSignal->GetZ(); z++) for (y = 0; y < RealPartSignal->GetY(); y++) for (x = 0; x < RealPartSignal->GetX(); x++){
    a=RealPartSignal->Get(x, y, z, 0);
    b=ImaginaryPartSignal->Get(x, y, z, 0);
    c=RealPartFilter->Get(x, y, z, 0)*CoefMult;
    d=ImaginaryPartFilter->Get(x, y, z, 0)*CoefMult;
    
    RealPartSignal->Put(x, y, z, 0, (a*c+b*d)/(c*c+d*d));
    ImaginaryPartSignal->Put(x, y, z, 0, (c*b-a*d)/(c*c+d*d));
  }
  
  //3) IFFT
  InverseFFT(RealPartSignal,ImaginaryPartSignal);
  InverseFFT(RealPartFilter,ImaginaryPartFilter);
}

