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
        this->_input2 = input2;
    else
        this->_input2 = NULL;

    if(input3 != NULL)
        this->_input3 = input3;
    else
        this->_input3 = NULL;

    if(input4 != NULL)
        this->_input4 = input4;
    else
        this->_input4 = NULL;

    _Lambda = 0;
}

template <class VoxelType> irtkNonLocalMedianFilter<VoxelType>::~irtkNonLocalMedianFilter(void)
{
}

template <class VoxelType> bool irtkNonLocalMedianFilter<VoxelType>::RequiresBuffering(void)
{
    return true;
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

    if(this->_input->GetZ() == 1){
        z1 = z;
        z2 = z;
    }

    if(this->_input4 == NULL)
       currentv = this->_input->GetAsDouble(x,y,z,t);
    else
       currentv = this->_input4->GetAsDouble(x,y,z,t);

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
                else if(i>=this->_input->GetX())
                    i1 = 2*this->_input->GetX() - i - 2;

                if(j < 0)
                    j1 = abs(j);
                else if(j>=this->_input->GetY())
                    j1 = 2*this->_input->GetY() - j - 2;

                if(k < 0)
                    k1 = abs(k);
                else if(k>=this->_input->GetZ())
                    k1 = 2*this->_input->GetZ() - k - 2;

                if(i1>=0 && i1 < this->_input->GetX()
                    && j1>=0 && j1 < this->_input->GetY()
                    && k1>=0 && k1 < this->_input->GetZ()){
                        // Now down to business
                        distancev = (i-x)*(i-x)*_dx*_dx
                            + (j-y)*(j-y)*_dy*_dy + (k-z)*(k-z)*_dz*_dz;
                        if(this->_input2 != NULL)
                            distancev += (this->_input2->GetAsDouble(i1,j1,k1)-this->_input2->GetAsDouble(x,y,z))
                            * (this->_input2->GetAsDouble(i1,j1,k1)-this->_input2->GetAsDouble(x,y,z))*_ds;

                        distancev = this->EvaluateWeight(distancev);

                        if(distancev > 0 && distancev < 1){
                            _localneighbor[index] = this->_input->GetAsDouble(i1,j1,k1,t);

                            _localweight[index] = distancev;

                            if(this->_input3 != NULL)
                                _localweight[index] = _localweight[index]*this->_input3->GetAsDouble(i1,j1,k1);

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
        return this->_input->GetAsDouble(x,y,z,t);

}

template <class VoxelType> void irtkNonLocalMedianFilter<VoxelType>::Initialize()
{

    //Run
    this->irtkImageToImage<VoxelType>::Initialize();

    int x,y,z,t,index;
    double threshold;

    if(_Sigma%2 == 0) _Sigma ++;
    // the neighborhood box must be odd
    _localweight = new float[_Sigma*_Sigma*_Sigma+1];
    _localneighbor = new float[_Sigma*_Sigma*_Sigma*2+2];

    if(this->_input2 != NULL){
        if(this->_input2->GetX() != this->_input->GetX() || this->_input2->GetY() != this->_input->GetY()
            || this->_input2->GetZ() != this->_input->GetZ()){
                cerr << "median fileter image dimension does not correspond!" << endl;
                exit(1);
        }
    }

    // edge
    _edge = new irtkGenericImage<VoxelType>(this->_input->GetImageAttributes());
    irtkGradientImage<VoxelType> gradient;
    gradient.SetInput (this->_input);
    gradient.SetOutput(_edge);
    gradient.Run();

    // Calculate edge using threshold
    for (t = 0; t < _edge->GetT(); t++) {
        threshold = 0;
        index = 0;
        for (z = 0; z < _edge->GetZ(); z++) {
            for (y = 0; y < _edge->GetY(); y++) {
                for (x = 0; x < this->_input->GetX(); x++) {
                     threshold += _edge->GetAsDouble(x,y,z,t);
                     index ++;
                }
            }
        }
        threshold = threshold/index;
        if(threshold > 0){
            for (z = 0; z < _edge->GetZ(); z++) {
                for (y = 0; y < _edge->GetY(); y++) {
                    for (x = 0; x < this->_input->GetX(); x++) {
                        if(_edge->GetAsDouble(x,y,z,t) <= threshold){
                            _edge->PutAsDouble(x,y,z,t,0);
                        }    
                    }
                }
            }
        }
    }

    if(this->_input2 == NULL){
        VoxelType min,max;
        this->_input->GetMinMax(&min,&max);
        _ds = max-min;
    }
    else{   
        irtkGreyPixel min,max;
        this->_input2->GetMinMax(&min,&max);
        _ds = max-min;
    }
    _ds = 128*128/_ds/_ds;

    _dx = this->_input->GetXSize();
    _dy = this->_input->GetYSize();
    _dz = this->_input->GetZSize();
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
