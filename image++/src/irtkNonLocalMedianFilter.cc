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
    int Sigma, irtkGenericImage<irtkGreyPixel>* input2, irtkGenericImage<irtkRealPixel>* input3
    , irtkGenericImage<VoxelType>* input4)
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

    if(input4 != NULL)
        _input4 = input4;
    else
        _input4 = NULL;

    _Lambda = 0;
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
    double distancev,currentv,sumofweight;
    int x1,y1,z1,x2,y2,z2,i,j,k,i1,j1,k1,index;

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

    if(_input->GetZ() == 1){
        z1 = z;
        z2 = z;
    }

    if(_input4 == NULL)
       currentv = _input->GetAsDouble(x,y,z,t);
    else
       currentv = _input4->GetAsDouble(x,y,z,t);

    index = 1;
    sumofweight = 0;
    // do the job find neighbor and weight
    for(k = z1; k <= z2; k++){
        for(j = y1; j <= y2; j++){
            for(i = x1; i<= x2; i++){
                // check if within image
                // mirror the image around boundary according to Deqing Sun's matlab code;
                i1 = i; j1 = j; k1 = k;
                if(i < 0) 
                    i1 = abs(i);
                else if(i>=_input->GetX()) 
                    i1 = 2*_input->GetX() - i - 2;

                if(j < 0)
                    j1 = abs(j);
                else if(j>=_input->GetY()) 
                    j1 = 2*_input->GetY() - j - 2;

                if(k < 0)
                    k1 = abs(k);
                else if(k>=_input->GetZ()) 
                    k1 = 2*_input->GetZ() - k - 2;

                if(i1>=0 && i1 < _input->GetX()
                    && j1>=0 && j1 < _input->GetY()
                    && k1>=0 && k1 < _input->GetZ()){
                        // now done to busniss
                        distancev = (i-x)*(i-x)*_dx*_dx
                            + (j-y)*(j-y)*_dy*_dy + (k-z)*(k-z)*_dz*_dz;
                        if(_input2 != NULL)
                            distancev += (_input2->GetAsDouble(i1,j1,k1)-_input2->GetAsDouble(x,y,z))
                            * (_input2->GetAsDouble(i1,j1,k1)-_input2->GetAsDouble(x,y,z))*_ds;

                        if(distancev > 0){
                            _localneighbor[index] = _input->GetAsDouble(i1,j1,k1,t);

                            _localweight[index] = this->EvaluateWeight(distancev);

                            if(_input3 != NULL)
                                _localweight[index] = _localweight[index]*_input3->GetAsDouble(i1,j1,k1);

                            //accumulate sum of weight
                            sumofweight += _localweight[index];
                            //index
                            index++;
                        }
                }
            }
        }
    }

    // normalize weight
    for(i = 1; i < index; i++){
        _localweight[i] /= sumofweight;
    }

    if(index > 1)
        return weightedmedian(index,_Lambda,currentv,_localneighbor,_localweight);
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

    if(_input2 == NULL){
        VoxelType min,max;
        _input->GetMinMax(&min,&max);
        _ds = max-min;
    }
    else{
        irtkGreyPixel min,max;
        _input2->GetMinMax(&min,&max);
        _ds = max-min;
    }
    _ds = 256*256/_ds/_ds;

    _dx = _input->GetXSize();
    _dy = _input->GetYSize();
    _dz = _input->GetZSize();
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