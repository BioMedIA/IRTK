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

#include <irtkNonLocalMedianFilter.h>

#include <irtkGradientImage.h>

#include <nr.h>

template <class VoxelType> irtkNonLocalMedianFilter<VoxelType>::irtkNonLocalMedianFilter(
    int Sigma, irtkGenericImage<irtkGreyPixel>* input2, irtkGenericImage<irtkRealPixel>* input3)
{
  _Sigma = Sigma;
  _localweight = NULL;
  _localneighbor = NULL;

  if(input2 != NULL)
      _input2 = input2;
  else
      _input2 = NULL;

  if(input3 != NULL)
      _input3 = input3;
  else
      _input3 = NULL;
}

template <class VoxelType> irtkNonLocalMedianFilter<VoxelType>::~irtkNonLocalMedianFilter(void)
{
}

template <class VoxelType> bool irtkNonLocalMedianFilter<VoxelType>::RequiresBuffering(void)
{
  return false;
}

template <class VoxelType> const char *irtkNonLocalMedianFilter<VoxelType>::NameOfClass()
{
  return "irtkNonLocalMedianFilter";
}

template <class VoxelType> double irtkNonLocalMedianFilter<VoxelType>::Run(int x, int y, int z, int t){
    double distancev,currentv,sumofweight,dx,dy,dz,ds;
    int x1,y1,z1,x2,y2,z2,i,j,k,index;

    // if edge use orignal size or use smaller window size
    if(_edge->Get(x,y,z,t) > 0){
        // intialize range
        x1 = x - _Sigma/2;
        x2 = x + _Sigma/2;
        y1 = y - _Sigma/2;
        y2 = y + _Sigma/2;
        z1 = z - _Sigma/2;
        z2 = z + _Sigma/2;
    }else{
        // intialize range
        x1 = x - _Sigma/4;
        x2 = x + _Sigma/4;
        y1 = y - _Sigma/4;
        y2 = y + _Sigma/4;
        z1 = z - _Sigma/4;
        z2 = z + _Sigma/4;
    }

    currentv = _input->GetAsDouble(x,y,z,t);
    
    index = 1;
    sumofweight = 0;
    dx = _input->GetXSize();
    dy = _input->GetYSize();
    dz = _input->GetZSize();
    // do the job find neighbor and weight
    for(k = z1; k <= z2; k++){
        for(j = y1; j <= y2; j++){
            for(i = x1; i<= x2; i++){
                // check if within image
                if(i>=0 && i < _input->GetX()
                    && j>=0 && j < _input->GetY()
                    && k>=0 && k < _input->GetZ()){
                        distancev = (i-x)*(i-x)*dx*dx
                            + (j-y)*(j-y)*dy*dy + (k-z)*(k-z)*dz*dz;
                    if(_input2 == NULL)
                        distancev += (_input->GetAsDouble(i,j,k,t)-_input->GetAsDouble(x,y,z,t))
                        * (_input->GetAsDouble(i,j,k,t)-_input->GetAsDouble(x,y,z,t));
                    else
                        distancev += (_input2->GetAsDouble(i,j,k)-_input2->GetAsDouble(x,y,z))
                        * (_input2->GetAsDouble(i,j,k)-_input2->GetAsDouble(x,y,z));


                    _localneighbor[index] = _input->GetAsDouble(i,j,k,t);

                    _localweight[index] = this->EvaluateWeight(distancev);

                    if(_input3 != NULL)
                        _localweight[index] = _localweight[index]*_input3->GetAsDouble(i,j,k);

                    //accumulate sum of weight
                    sumofweight += _localweight[index];
                    //index
                    index++;
                }
            }
        }
    }

    // normalize weight
    for(i = 1; i < index; i++){
        _localweight[i] /= sumofweight;
    }

    if(index > 1)
        return weightedmedian(index,0,currentv,_localneighbor,_localweight);
    else
        return _input->GetAsDouble(x,y,z,t);

}

template <class VoxelType> void irtkNonLocalMedianFilter<VoxelType>::Initialize()
{

    //Run
    this->irtkImageToImage<VoxelType>::Initialize();

    int x,y,z,t;
    double threshold;

    if(_Sigma%2 == 0) _Sigma ++;
    // the neighborhood box must be odd
   _localweight = new float[_Sigma*_Sigma*_Sigma+1];
   _localneighbor = new float[_Sigma*_Sigma*_Sigma*2+2];

   if(_input2 != NULL){
       if(_input2->GetX() != _input->GetX() || _input2->GetY() != _input->GetY() 
           || _input2->GetZ() != _input->GetZ()){
               cerr << "median fileter image dimension does not correspond!" << endl;
               exit(1);
       }
   }

   // edge
   _edge = new irtkGenericImage<VoxelType>(_input->GetImageAttributes());
   irtkGradientImage<VoxelType> gradient;
   gradient.SetInput (_input);
   gradient.SetOutput(_edge);
   gradient.Run();

   // Calculate edge using threshold
   for (t = 0; t < _edge->GetT(); t++) {
       threshold = 0;
       for (z = 0; z < _edge->GetZ(); z++) {
           for (y = 0; y < _edge->GetY(); y++) {
               for (x = 0; x < _input->GetX(); x++) {
                   if(_edge->GetAsDouble(x,y,z,t)>threshold){
                       threshold = _edge->GetAsDouble(x,y,z,t);
                   }
               }
           }
       }
       threshold = threshold/2.0;
       for (z = 0; z < _edge->GetZ(); z++) {
           for (y = 0; y < _edge->GetY(); y++) {
               for (x = 0; x < _input->GetX(); x++) {
                   if(_edge->GetAsDouble(x,y,z,t) <= threshold){
                       _edge->PutAsDouble(x,y,z,t,0);
                   }
               }
           }
       }
   }
}

template <class VoxelType> void irtkNonLocalMedianFilter<VoxelType>::Run()
{
    //Finalize
    this->irtkImageToImage<VoxelType>::Run();

}

template <class VoxelType> void irtkNonLocalMedianFilter<VoxelType>::Finalize()
{
    //Finalize
    this->irtkImageToImage<VoxelType>::Finalize();

    delete []_localweight;
    delete []_localneighbor;
    _localweight = NULL;
    _localneighbor = NULL;

    delete _edge;
}

template class irtkNonLocalMedianFilter<unsigned char>;
template class irtkNonLocalMedianFilter<short>;
template class irtkNonLocalMedianFilter<unsigned short>;
template class irtkNonLocalMedianFilter<float>;
template class irtkNonLocalMedianFilter<double>;