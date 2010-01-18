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
  CoefMult=(float)(sqrt((double)RealPartSignal->GetX())*sqrt((double)RealPartSignal->GetY())*sqrt((double)RealPartSignal->GetZ()));
  
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


void ConvolutionInFourierNoFilterTransfo(irtkGenericImage<float> * RealPartSignal,irtkGenericImage<float> * ImaginaryPartSignal,irtkGenericImage<float> * RealPartFilterTransformedFrSpace,irtkGenericImage<float> * ImaginaryPartFilterTransformedFrSpace){
  int x,y,z;
  float a,b,c,d;
  float CoefMult;
  
  
  //1) FFT
  DirectFFT(RealPartSignal,ImaginaryPartSignal);
  
  //2) filtering in Fourier spaces
  CoefMult=(float)(sqrt((double)RealPartSignal->GetX())*sqrt((double)RealPartSignal->GetY())*sqrt((double)RealPartSignal->GetZ()));
  
  for (z = 0; z < RealPartSignal->GetZ(); z++) for (y = 0; y < RealPartSignal->GetY(); y++) for (x = 0; x < RealPartSignal->GetX(); x++){
    a=RealPartSignal->Get(x, y, z, 0);
    b=ImaginaryPartSignal->Get(x, y, z, 0);
    c=RealPartFilterTransformedFrSpace->Get(x, y, z, 0)*CoefMult;
    d=ImaginaryPartFilterTransformedFrSpace->Get(x, y, z, 0)*CoefMult;
    
    RealPartSignal->Put(x, y, z, 0, a*c-b*d);
    ImaginaryPartSignal->Put(x, y, z, 0, c*b+a*d);
  }
  
  //3) IFFT
  InverseFFT(RealPartSignal,ImaginaryPartSignal);
}







void DeconvolutionInFourier(irtkGenericImage<float> * RealPartSignal,irtkGenericImage<float> * ImaginaryPartSignal,irtkGenericImage<float> * RealPartFilter,irtkGenericImage<float> * ImaginaryPartFilter){
  int x,y,z;
  float a,b,c,d;
  float CoefMult;
  
  
  //1) FFT
  DirectFFT(RealPartSignal,ImaginaryPartSignal);
  DirectFFT(RealPartFilter,ImaginaryPartFilter);
  
  //2) filtering in Fourier spaces
  CoefMult=(float)(sqrt((double)RealPartSignal->GetX())*sqrt((double)RealPartSignal->GetY())*sqrt((double)RealPartSignal->GetZ()));
  
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





void DeconvolutionInFourierNoFilterTransfo(irtkGenericImage<float> * RealPartSignal,irtkGenericImage<float> * ImaginaryPartSignal,irtkGenericImage<float> * RealPartFilterTransformedFrSpace,irtkGenericImage<float> * ImaginaryPartFilterTransformedFrSpace){
  int x,y,z;
  float a,b,c,d;
  float CoefMult;
  
  
  //1) FFT
  DirectFFT(RealPartSignal,ImaginaryPartSignal);
  
  //2) filtering in Fourier spaces
  CoefMult=(float)(sqrt((double)RealPartSignal->GetX())*sqrt((double)RealPartSignal->GetY())*sqrt((double)RealPartSignal->GetZ()));
  
  for (z = 0; z < RealPartSignal->GetZ(); z++) for (y = 0; y < RealPartSignal->GetY(); y++) for (x = 0; x < RealPartSignal->GetX(); x++){
    a=RealPartSignal->Get(x, y, z, 0);
    b=ImaginaryPartSignal->Get(x, y, z, 0);
    c=RealPartFilterTransformedFrSpace->Get(x, y, z, 0)*CoefMult;
    d=ImaginaryPartFilterTransformedFrSpace->Get(x, y, z, 0)*CoefMult;
    
    RealPartSignal->Put(x, y, z, 0, (a*c+b*d)/(c*c+d*d));
    ImaginaryPartSignal->Put(x, y, z, 0, (c*b-a*d)/(c*c+d*d));
  }
  
  //3) IFFT
  InverseFFT(RealPartSignal,ImaginaryPartSignal);
}






void MakeGaussianFilter(float sigma,irtkGenericImage<float> * RealPartFilter,irtkGenericImage<float> * ImaginaryPartFilter){
  int x,y,z;
  double TempDouble;
  int NX,NY,NZ;
  float SumLoc;
  
  //variables
  NX=RealPartFilter->GetX();
  NY=RealPartFilter->GetY();
  NZ=RealPartFilter->GetZ();
  
  //background at 0
  for (z=0;z<NZ;z++) for (y=0;y<NY;y++) for (x=0;x<NX;x++) RealPartFilter->Put(x,y,z,0,0);
  for (z=0;z<NZ;z++) for (y=0;y<NY;y++) for (x=0;x<NX;x++) ImaginaryPartFilter->Put(x,y,z,0,0);

  //gaussian filter
  for (z=0;z<NZ/2;z++) for (y=0;y<NY/2;y++) for (x=0;x<NX/2;x++){
    TempDouble=exp(-(x*x+y*y+z*z)/(2.*sigma*sigma));
    RealPartFilter->Put(x,y,z,0,(float)(TempDouble));
  }
  for (z=0;z<NZ/2;z++) for (y=0;y<NY/2;y++) for (x=NX/2;x<NX;x++){
    TempDouble=exp(-((NX-x)*(NX-x)+y*y+z*z)/(2.*sigma*sigma));
    RealPartFilter->Put(x,y,z,0,(float)(TempDouble));
  }
  for (z=0;z<NZ/2;z++) for (y=NY/2;y<NY;y++) for (x=0;x<NX/2;x++){
    TempDouble=exp(-(x*x+(NY-y)*(NY-y)+z*z)/(2.*sigma*sigma));
    RealPartFilter->Put(x,y,z,0,(float)(TempDouble));
  }
  for (z=0;z<NZ/2;z++) for (y=NY/2;y<NY;y++) for (x=NX/2;x<NX;x++){
    TempDouble=exp(-((NX-x)*(NX-x)+(NY-y)*(NY-y)+z*z)/(2.*sigma*sigma));
    RealPartFilter->Put(x,y,z,0,(float)(TempDouble));
  }
  for (z=NZ/2;z<NZ;z++) for (y=0;y<NY/2;y++) for (x=0;x<NX/2;x++){
    TempDouble=exp(-(x*x+y*y+(NZ-z)*(NZ-z))/(2.*sigma*sigma));
    RealPartFilter->Put(x,y,z,0,(float)(TempDouble));
  }
  for (z=NZ/2;z<NZ;z++) for (y=0;y<NY/2;y++) for (x=NX/2;x<NX;x++){
    TempDouble=exp(-((NX-x)*(NX-x)+y*y+(NZ-z)*(NZ-z))/(2.*sigma*sigma));
    RealPartFilter->Put(x,y,z,0,(float)(TempDouble));
  }
  for (z=NZ/2;z<NZ;z++) for (y=NY/2;y<NY;y++) for (x=0;x<NX/2;x++){
    TempDouble=exp(-(x*x+(NY-y)*(NY-y)+(NZ-z)*(NZ-z))/(2.*sigma*sigma));
    RealPartFilter->Put(x,y,z,0,(float)(TempDouble));
  }
  for (z=NZ/2;z<NZ;z++) for (y=NY/2;y<NY;y++) for (x=NX/2;x<NX;x++){
    TempDouble=exp(-((NX-x)*(NX-x)+(NY-y)*(NY-y)+(NZ-z)*(NZ-z))/(2.*sigma*sigma));
    RealPartFilter->Put(x,y,z,0,(float)(TempDouble));
  }
  
  //normalization
  SumLoc=0.;
  for (z=0;z<NZ;z++) for (y=0;y<NY;y++) for (x=0;x<NX;x++) SumLoc+=RealPartFilter->Get(x,y,z,0);
  for (z=0;z<NZ;z++) for (y=0;y<NY;y++) for (x=0;x<NX;x++) RealPartFilter->Put(x,y,z,0,RealPartFilter->Get(x,y,z,0)/SumLoc);
}





void MakeAnisotropicGaussianFilter(float weight,float sigmaX,float sigmaY,float sigmaZ,irtkGenericImage<float> * RealPartFilter,irtkGenericImage<float> * ImaginaryPartFilter){
  int x,y,z;
  double TempDouble;
  int NX,NY,NZ;
  float SumLoc;
  
  //variables
  NX=RealPartFilter->GetX();
  NY=RealPartFilter->GetY();
  NZ=RealPartFilter->GetZ();
  
  //background at 0
  for (z=0;z<NZ;z++) for (y=0;y<NY;y++) for (x=0;x<NX;x++) RealPartFilter->Put(x,y,z,0,0);
  for (z=0;z<NZ;z++) for (y=0;y<NY;y++) for (x=0;x<NX;x++) ImaginaryPartFilter->Put(x,y,z,0,0);
  
  //gaussian filter
  for (z=0;z<NZ/2;z++) for (y=0;y<NY/2;y++) for (x=0;x<NX/2;x++){
    TempDouble=exp( -(float)(x*x)/(2.*sigmaX*sigmaX) -(float)(y*y)/(2.*sigmaY*sigmaY) -(float)(z*z)/(2.*sigmaZ*sigmaZ) );
    RealPartFilter->Put(x,y,z,0,(float)(TempDouble));
  }
  for (z=0;z<NZ/2;z++) for (y=0;y<NY/2;y++) for (x=NX/2;x<NX;x++){
    TempDouble=exp( -(float)((NX-x)*(NX-x))/(2.*sigmaX*sigmaX) -(float)(y*y)/(2.*sigmaY*sigmaY) -(float)(z*z)/(2.*sigmaZ*sigmaZ) );
    RealPartFilter->Put(x,y,z,0,(float)(TempDouble));
  }
  for (z=0;z<NZ/2;z++) for (y=NY/2;y<NY;y++) for (x=0;x<NX/2;x++){
    TempDouble=exp( -(float)(x*x)/(2.*sigmaX*sigmaX) -(float)((NY-y)*(NY-y))/(2.*sigmaY*sigmaY) -(float)(z*z)/(2.*sigmaZ*sigmaZ));
    RealPartFilter->Put(x,y,z,0,(float)(TempDouble));
  }
  for (z=0;z<NZ/2;z++) for (y=NY/2;y<NY;y++) for (x=NX/2;x<NX;x++){
    TempDouble=exp( -(float)((NX-x)*(NX-x))/(2.*sigmaX*sigmaX) -(float)((NY-y)*(NY-y))/(2.*sigmaY*sigmaY) -(float)(z*z)/(2.*sigmaZ*sigmaZ));
    RealPartFilter->Put(x,y,z,0,(float)(TempDouble));
  }
  for (z=NZ/2;z<NZ;z++) for (y=0;y<NY/2;y++) for (x=0;x<NX/2;x++){
    TempDouble=exp(-(float)(x*x)/(2.*sigmaX*sigmaX) -(float)(y*y)/(2.*sigmaY*sigmaY) -(float)((NZ-z)*(NZ-z))/(2.*sigmaZ*sigmaZ));
    RealPartFilter->Put(x,y,z,0,(float)(TempDouble));
  }
  for (z=NZ/2;z<NZ;z++) for (y=0;y<NY/2;y++) for (x=NX/2;x<NX;x++){
    TempDouble=exp(-(float)((NX-x)*(NX-x))/(2.*sigmaX*sigmaX) -(float)(y*y)/(2.*sigmaY*sigmaY) -(float)((NZ-z)*(NZ-z))/(2.*sigmaZ*sigmaZ));
    RealPartFilter->Put(x,y,z,0,(float)(TempDouble));
  }
  for (z=NZ/2;z<NZ;z++) for (y=NY/2;y<NY;y++) for (x=0;x<NX/2;x++){
    TempDouble=exp(-(float)(x*x)/(2.*sigmaX*sigmaX) -(float)((NY-y)*(NY-y))/(2.*sigmaY*sigmaY) -(float)((NZ-z)*(NZ-z))/(2.*sigmaZ*sigmaZ));
    RealPartFilter->Put(x,y,z,0,(float)(TempDouble));
  }
  for (z=NZ/2;z<NZ;z++) for (y=NY/2;y<NY;y++) for (x=NX/2;x<NX;x++){
    TempDouble=exp(-(float)((NX-x)*(NX-x))/(2.*sigmaX*sigmaX) -(float)((NY-y)*(NY-y))/(2.*sigmaY*sigmaY) -(float)((NZ-z)*(NZ-z))/(2.*sigmaZ*sigmaZ));
    RealPartFilter->Put(x,y,z,0,(float)(TempDouble));
  }

  //normalization
  SumLoc=0.;
  for (z=0;z<NZ;z++) for (y=0;y<NY;y++) for (x=0;x<NX;x++) SumLoc+=RealPartFilter->Get(x,y,z,0);
  for (z=0;z<NZ;z++) for (y=0;y<NY;y++) for (x=0;x<NX;x++) RealPartFilter->Put(x,y,z,0,weight*RealPartFilter->Get(x,y,z,0)/SumLoc);
}






void MakeSumOf2AnisotropicGaussianFilters(float weight1,float sigmaX1,float sigmaY1,float sigmaZ1,float weight2,float sigmaX2,float sigmaY2,float sigmaZ2,irtkGenericImage<float> * RealPartFilter,irtkGenericImage<float> * ImaginaryPartFilter){
  int x,y,z;
  int NX,NY,NZ;
  float SumLoc;
  irtkGenericImage<float> ImageTemp;
  float sigmaX,sigmaY,sigmaZ;
      
  //variables and allocation
  NX=RealPartFilter->GetX();
  NY=RealPartFilter->GetY();
  NZ=RealPartFilter->GetZ();
  ImageTemp = irtkGenericImage<float>(NX,NY,NZ,1);
  
  //background at 0
  for (z=0;z<NZ;z++) for (y=0;y<NY;y++) for (x=0;x<NX;x++) RealPartFilter->Put(x,y,z,0,0);
  for (z=0;z<NZ;z++) for (y=0;y<NY;y++) for (x=0;x<NX;x++) ImaginaryPartFilter->Put(x,y,z,0,0);

  //gaussian filter n1
  sigmaX=sigmaX1;
  sigmaY=sigmaY1;
  sigmaZ=sigmaZ1;
  for (z=0;z<NZ/2;z++) for (y=0;y<NY/2;y++) for (x=0;x<NX/2;x++)
        ImageTemp.Put(x,y,z,0,(float)(exp( -(float)(x*x)/(2.*sigmaX*sigmaX) -(float)(y*y)/(2.*sigmaY*sigmaY) -(float)(z*z)/(2.*sigmaZ*sigmaZ))));
  for (z=0;z<NZ/2;z++) for (y=0;y<NY/2;y++) for (x=NX/2;x<NX;x++)
        ImageTemp.Put(x,y,z,0,(float)(exp( -(float)((NX-x)*(NX-x))/(2.*sigmaX*sigmaX) -(float)(y*y)/(2.*sigmaY*sigmaY) -(float)(z*z)/(2.*sigmaZ*sigmaZ))));
  for (z=0;z<NZ/2;z++) for (y=NY/2;y<NY;y++) for (x=0;x<NX/2;x++)
        ImageTemp.Put(x,y,z,0,(float)(exp( -(float)(x*x)/(2.*sigmaX*sigmaX) -(float)((NY-y)*(NY-y))/(2.*sigmaY*sigmaY) -(float)(z*z)/(2.*sigmaZ*sigmaZ))));
  for (z=0;z<NZ/2;z++) for (y=NY/2;y<NY;y++) for (x=NX/2;x<NX;x++)
        ImageTemp.Put(x,y,z,0,(float)(exp( -(float)((NX-x)*(NX-x))/(2.*sigmaX*sigmaX) -(float)((NY-y)*(NY-y))/(2.*sigmaY*sigmaY) -(float)(z*z)/(2.*sigmaZ*sigmaZ))));
  for (z=NZ/2;z<NZ;z++) for (y=0;y<NY/2;y++) for (x=0;x<NX/2;x++)
        ImageTemp.Put(x,y,z,0,(float)(exp(-(float)(x*x)/(2.*sigmaX*sigmaX) -(float)(y*y)/(2.*sigmaY*sigmaY) -(float)((NZ-z)*(NZ-z))/(2.*sigmaZ*sigmaZ))));
  for (z=NZ/2;z<NZ;z++) for (y=0;y<NY/2;y++) for (x=NX/2;x<NX;x++)
        ImageTemp.Put(x,y,z,0,(float)(exp(-(float)((NX-x)*(NX-x))/(2.*sigmaX*sigmaX) -(float)(y*y)/(2.*sigmaY*sigmaY) -(float)((NZ-z)*(NZ-z))/(2.*sigmaZ*sigmaZ))));
  for (z=NZ/2;z<NZ;z++) for (y=NY/2;y<NY;y++) for (x=0;x<NX/2;x++)
        ImageTemp.Put(x,y,z,0,(float)(exp(-(float)(x*x)/(2.*sigmaX*sigmaX) -(float)((NY-y)*(NY-y))/(2.*sigmaY*sigmaY) -(float)((NZ-z)*(NZ-z))/(2.*sigmaZ*sigmaZ))));
  for (z=NZ/2;z<NZ;z++) for (y=NY/2;y<NY;y++) for (x=NX/2;x<NX;x++)
        ImageTemp.Put(x,y,z,0,(float)(exp(-(float)((NX-x)*(NX-x))/(2.*sigmaX*sigmaX) -(float)((NY-y)*(NY-y))/(2.*sigmaY*sigmaY) -(float)((NZ-z)*(NZ-z))/(2.*sigmaZ*sigmaZ))));
  
  //normalization of filter n1 and copy in RealPartFilter
  SumLoc=0.;
  for (z=0;z<NZ;z++) for (y=0;y<NY;y++) for (x=0;x<NX;x++) SumLoc+=ImageTemp.Get(x,y,z,0);
  for (z=0;z<NZ;z++) for (y=0;y<NY;y++) for (x=0;x<NX;x++) RealPartFilter->Put(x,y,z,0,RealPartFilter->Get(x,y,z,0)+weight1*ImageTemp.Get(x,y,z,0)/SumLoc);

  //gaussian filter n2
  sigmaX=sigmaX2;
  sigmaY=sigmaY2;
  sigmaZ=sigmaZ2;
  for (z=0;z<NZ/2;z++) for (y=0;y<NY/2;y++) for (x=0;x<NX/2;x++)
        ImageTemp.Put(x,y,z,0,(float)(exp( -(float)(x*x)/(2.*sigmaX*sigmaX) -(float)(y*y)/(2.*sigmaY*sigmaY) -(float)(z*z)/(2.*sigmaZ*sigmaZ))));
  for (z=0;z<NZ/2;z++) for (y=0;y<NY/2;y++) for (x=NX/2;x<NX;x++)
        ImageTemp.Put(x,y,z,0,(float)(exp( -(float)((NX-x)*(NX-x))/(2.*sigmaX*sigmaX) -(float)(y*y)/(2.*sigmaY*sigmaY) -(float)(z*z)/(2.*sigmaZ*sigmaZ))));
  for (z=0;z<NZ/2;z++) for (y=NY/2;y<NY;y++) for (x=0;x<NX/2;x++)
        ImageTemp.Put(x,y,z,0,(float)(exp( -(float)(x*x)/(2.*sigmaX*sigmaX) -(float)((NY-y)*(NY-y))/(2.*sigmaY*sigmaY) -(float)(z*z)/(2.*sigmaZ*sigmaZ))));
  for (z=0;z<NZ/2;z++) for (y=NY/2;y<NY;y++) for (x=NX/2;x<NX;x++)
        ImageTemp.Put(x,y,z,0,(float)(exp( -(float)((NX-x)*(NX-x))/(2.*sigmaX*sigmaX) -(float)((NY-y)*(NY-y))/(2.*sigmaY*sigmaY) -(float)(z*z)/(2.*sigmaZ*sigmaZ))));
  for (z=NZ/2;z<NZ;z++) for (y=0;y<NY/2;y++) for (x=0;x<NX/2;x++)
        ImageTemp.Put(x,y,z,0,(float)(exp(-(float)(x*x)/(2.*sigmaX*sigmaX) -(float)(y*y)/(2.*sigmaY*sigmaY) -(float)((NZ-z)*(NZ-z))/(2.*sigmaZ*sigmaZ))));
  for (z=NZ/2;z<NZ;z++) for (y=0;y<NY/2;y++) for (x=NX/2;x<NX;x++)
        ImageTemp.Put(x,y,z,0,(float)(exp(-(float)((NX-x)*(NX-x))/(2.*sigmaX*sigmaX) -(float)(y*y)/(2.*sigmaY*sigmaY) -(float)((NZ-z)*(NZ-z))/(2.*sigmaZ*sigmaZ))));
  for (z=NZ/2;z<NZ;z++) for (y=NY/2;y<NY;y++) for (x=0;x<NX/2;x++)
        ImageTemp.Put(x,y,z,0,(float)(exp(-(float)(x*x)/(2.*sigmaX*sigmaX) -(float)((NY-y)*(NY-y))/(2.*sigmaY*sigmaY) -(float)((NZ-z)*(NZ-z))/(2.*sigmaZ*sigmaZ))));
  for (z=NZ/2;z<NZ;z++) for (y=NY/2;y<NY;y++) for (x=NX/2;x<NX;x++)
        ImageTemp.Put(x,y,z,0,(float)(exp(-(float)((NX-x)*(NX-x))/(2.*sigmaX*sigmaX) -(float)((NY-y)*(NY-y))/(2.*sigmaY*sigmaY) -(float)((NZ-z)*(NZ-z))/(2.*sigmaZ*sigmaZ))));
  
  //normalization of filter n2 and copy in RealPartFilter
  SumLoc=0.;
  for (z=0;z<NZ;z++) for (y=0;y<NY;y++) for (x=0;x<NX;x++) SumLoc+=ImageTemp.Get(x,y,z,0);
  for (z=0;z<NZ;z++) for (y=0;y<NY;y++) for (x=0;x<NX;x++) RealPartFilter->Put(x,y,z,0,RealPartFilter->Get(x,y,z,0)+weight2*ImageTemp.Get(x,y,z,0)/SumLoc);


}


void MakeSumOf3AnisotropicGaussianFilters(float weight1,float sigmaX1,float sigmaY1,float sigmaZ1,float weight2,float sigmaX2,float sigmaY2,float sigmaZ2,float weight3,float sigmaX3,float sigmaY3,float sigmaZ3,irtkGenericImage<float> * RealPartFilter,irtkGenericImage<float> * ImaginaryPartFilter){
  int x,y,z;
  int NX,NY,NZ;
  float SumLoc;
  irtkGenericImage<float> ImageTemp;
  float sigmaX,sigmaY,sigmaZ;
  
  //init with the 2 first gaussian filters
  MakeSumOf2AnisotropicGaussianFilters(weight1,sigmaX1,sigmaY1,sigmaZ1,weight2,sigmaX2,sigmaY2,sigmaZ2,RealPartFilter,ImaginaryPartFilter);

  //variables and allocation
  NX=RealPartFilter->GetX();
  NY=RealPartFilter->GetY();
  NZ=RealPartFilter->GetZ();
  ImageTemp = irtkGenericImage<float>(NX,NY,NZ,1);
  
  //gaussian filter n3
  sigmaX=sigmaX3;
  sigmaY=sigmaY3;
  sigmaZ=sigmaZ3;
  for (z=0;z<NZ/2;z++) for (y=0;y<NY/2;y++) for (x=0;x<NX/2;x++)
        ImageTemp.Put(x,y,z,0,(float)(exp( -(float)(x*x)/(2.*sigmaX*sigmaX) -(float)(y*y)/(2.*sigmaY*sigmaY) -(float)(z*z)/(2.*sigmaZ*sigmaZ))));
  for (z=0;z<NZ/2;z++) for (y=0;y<NY/2;y++) for (x=NX/2;x<NX;x++)
        ImageTemp.Put(x,y,z,0,(float)(exp( -(float)((NX-x)*(NX-x))/(2.*sigmaX*sigmaX) -(float)(y*y)/(2.*sigmaY*sigmaY) -(float)(z*z)/(2.*sigmaZ*sigmaZ))));
  for (z=0;z<NZ/2;z++) for (y=NY/2;y<NY;y++) for (x=0;x<NX/2;x++)
        ImageTemp.Put(x,y,z,0,(float)(exp( -(float)(x*x)/(2.*sigmaX*sigmaX) -(float)((NY-y)*(NY-y))/(2.*sigmaY*sigmaY) -(float)(z*z)/(2.*sigmaZ*sigmaZ))));
  for (z=0;z<NZ/2;z++) for (y=NY/2;y<NY;y++) for (x=NX/2;x<NX;x++)
        ImageTemp.Put(x,y,z,0,(float)(exp( -(float)((NX-x)*(NX-x))/(2.*sigmaX*sigmaX) -(float)((NY-y)*(NY-y))/(2.*sigmaY*sigmaY) -(float)(z*z)/(2.*sigmaZ*sigmaZ))));
  for (z=NZ/2;z<NZ;z++) for (y=0;y<NY/2;y++) for (x=0;x<NX/2;x++)
        ImageTemp.Put(x,y,z,0,(float)(exp(-(float)(x*x)/(2.*sigmaX*sigmaX) -(float)(y*y)/(2.*sigmaY*sigmaY) -(float)((NZ-z)*(NZ-z))/(2.*sigmaZ*sigmaZ))));
  for (z=NZ/2;z<NZ;z++) for (y=0;y<NY/2;y++) for (x=NX/2;x<NX;x++)
        ImageTemp.Put(x,y,z,0,(float)(exp(-(float)((NX-x)*(NX-x))/(2.*sigmaX*sigmaX) -(float)(y*y)/(2.*sigmaY*sigmaY) -(float)((NZ-z)*(NZ-z))/(2.*sigmaZ*sigmaZ))));
  for (z=NZ/2;z<NZ;z++) for (y=NY/2;y<NY;y++) for (x=0;x<NX/2;x++)
        ImageTemp.Put(x,y,z,0,(float)(exp(-(float)(x*x)/(2.*sigmaX*sigmaX) -(float)((NY-y)*(NY-y))/(2.*sigmaY*sigmaY) -(float)((NZ-z)*(NZ-z))/(2.*sigmaZ*sigmaZ))));
  for (z=NZ/2;z<NZ;z++) for (y=NY/2;y<NY;y++) for (x=NX/2;x<NX;x++)
        ImageTemp.Put(x,y,z,0,(float)(exp(-(float)((NX-x)*(NX-x))/(2.*sigmaX*sigmaX) -(float)((NY-y)*(NY-y))/(2.*sigmaY*sigmaY) -(float)((NZ-z)*(NZ-z))/(2.*sigmaZ*sigmaZ))));
  
  //normalization of filter n3 and copy in RealPartFilter
  SumLoc=0.;
  for (z=0;z<NZ;z++) for (y=0;y<NY;y++) for (x=0;x<NX;x++) SumLoc+=ImageTemp.Get(x,y,z,0);
  for (z=0;z<NZ;z++) for (y=0;y<NY;y++) for (x=0;x<NX;x++) RealPartFilter->Put(x,y,z,0,RealPartFilter->Get(x,y,z,0)+weight3*ImageTemp.Get(x,y,z,0)/SumLoc);
}


void MakeSumOf4AnisotropicGaussianFilters(float weight1,float sigmaX1,float sigmaY1,float sigmaZ1,float weight2,float sigmaX2,float sigmaY2,float sigmaZ2,float weight3,float sigmaX3,float sigmaY3,float sigmaZ3,float weight4,float sigmaX4,float sigmaY4,float sigmaZ4,irtkGenericImage<float> * RealPartFilter,irtkGenericImage<float> * ImaginaryPartFilter){
  int x,y,z;
  int NX,NY,NZ;
  float SumLoc;
  irtkGenericImage<float> ImageTemp;
  float sigmaX,sigmaY,sigmaZ;
  
  //init with the 3 first gaussian filters
  MakeSumOf3AnisotropicGaussianFilters(weight1,sigmaX1,sigmaY1,sigmaZ1,weight2,sigmaX2,sigmaY2,sigmaZ2,weight3,sigmaX3,sigmaY3,sigmaZ3,RealPartFilter,ImaginaryPartFilter);

  //variables and allocation
  NX=RealPartFilter->GetX();
  NY=RealPartFilter->GetY();
  NZ=RealPartFilter->GetZ();
  ImageTemp = irtkGenericImage<float>(NX,NY,NZ,1);
  
  //gaussian filter n4
  sigmaX=sigmaX4;
  sigmaY=sigmaY4;
  sigmaZ=sigmaZ4;
  for (z=0;z<NZ/2;z++) for (y=0;y<NY/2;y++) for (x=0;x<NX/2;x++)
        ImageTemp.Put(x,y,z,0,(float)(exp( -(float)(x*x)/(2.*sigmaX*sigmaX) -(float)(y*y)/(2.*sigmaY*sigmaY) -(float)(z*z)/(2.*sigmaZ*sigmaZ))));
  for (z=0;z<NZ/2;z++) for (y=0;y<NY/2;y++) for (x=NX/2;x<NX;x++)
        ImageTemp.Put(x,y,z,0,(float)(exp( -(float)((NX-x)*(NX-x))/(2.*sigmaX*sigmaX) -(float)(y*y)/(2.*sigmaY*sigmaY) -(float)(z*z)/(2.*sigmaZ*sigmaZ))));
  for (z=0;z<NZ/2;z++) for (y=NY/2;y<NY;y++) for (x=0;x<NX/2;x++)
        ImageTemp.Put(x,y,z,0,(float)(exp( -(float)(x*x)/(2.*sigmaX*sigmaX) -(float)((NY-y)*(NY-y))/(2.*sigmaY*sigmaY) -(float)(z*z)/(2.*sigmaZ*sigmaZ))));
  for (z=0;z<NZ/2;z++) for (y=NY/2;y<NY;y++) for (x=NX/2;x<NX;x++)
        ImageTemp.Put(x,y,z,0,(float)(exp( -(float)((NX-x)*(NX-x))/(2.*sigmaX*sigmaX) -(float)((NY-y)*(NY-y))/(2.*sigmaY*sigmaY) -(float)(z*z)/(2.*sigmaZ*sigmaZ))));
  for (z=NZ/2;z<NZ;z++) for (y=0;y<NY/2;y++) for (x=0;x<NX/2;x++)
        ImageTemp.Put(x,y,z,0,(float)(exp(-(float)(x*x)/(2.*sigmaX*sigmaX) -(float)(y*y)/(2.*sigmaY*sigmaY) -(float)((NZ-z)*(NZ-z))/(2.*sigmaZ*sigmaZ))));
  for (z=NZ/2;z<NZ;z++) for (y=0;y<NY/2;y++) for (x=NX/2;x<NX;x++)
        ImageTemp.Put(x,y,z,0,(float)(exp(-(float)((NX-x)*(NX-x))/(2.*sigmaX*sigmaX) -(float)(y*y)/(2.*sigmaY*sigmaY) -(float)((NZ-z)*(NZ-z))/(2.*sigmaZ*sigmaZ))));
  for (z=NZ/2;z<NZ;z++) for (y=NY/2;y<NY;y++) for (x=0;x<NX/2;x++)
        ImageTemp.Put(x,y,z,0,(float)(exp(-(float)(x*x)/(2.*sigmaX*sigmaX) -(float)((NY-y)*(NY-y))/(2.*sigmaY*sigmaY) -(float)((NZ-z)*(NZ-z))/(2.*sigmaZ*sigmaZ))));
  for (z=NZ/2;z<NZ;z++) for (y=NY/2;y<NY;y++) for (x=NX/2;x<NX;x++)
        ImageTemp.Put(x,y,z,0,(float)(exp(-(float)((NX-x)*(NX-x))/(2.*sigmaX*sigmaX) -(float)((NY-y)*(NY-y))/(2.*sigmaY*sigmaY) -(float)((NZ-z)*(NZ-z))/(2.*sigmaZ*sigmaZ))));
  
  //normalization of filter n4 and copy in RealPartFilter
  SumLoc=0.;
  for (z=0;z<NZ;z++) for (y=0;y<NY;y++) for (x=0;x<NX;x++) SumLoc+=ImageTemp.Get(x,y,z,0);
  for (z=0;z<NZ;z++) for (y=0;y<NY;y++) for (x=0;x<NX;x++) RealPartFilter->Put(x,y,z,0,RealPartFilter->Get(x,y,z,0)+weight4*ImageTemp.Get(x,y,z,0)/SumLoc);
}


