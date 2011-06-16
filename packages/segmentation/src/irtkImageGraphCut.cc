/*=========================================================================

Library   : Image Registration Toolkit (IRTK)
Module    : $Id$
Copyright : Imperial College, Department of Computing
Visual Information Processing (VIP), 2008 onwards
Date      : $Date$
Version   : $Revision$
Changes   : $Author$

=========================================================================*/

#include <irtkImageGraphCut.h>

template <class VoxelType> irtkImageGraphCut<VoxelType>::irtkImageGraphCut()
{
    // Set in- and outputs
    _input  = NULL;
    _output = NULL;
    /// source weight image
    _weight = NULL;
    _mode = 4;
    _numberOfImages = 0;
    _numberOfWeights = 0.0;
    _labels = 2;
}

template <class VoxelType> irtkImageGraphCut<VoxelType>::~irtkImageGraphCut()
{
    // Set in- and outputs
    if(_input)
        delete []_input;
    if(_weight)
        delete []_weight;
    _input  = NULL;
    _output = NULL;
    /// source weight image
    _weight = NULL;
}

template <class VoxelType> void irtkImageGraphCut<VoxelType>::SetInput(int n, irtkGenericImage<VoxelType> **image,int l, irtkRealImage **weight)
{
    _numberOfImages = n;
    _numberOfWeights = 1.0;
    _input = new irtkGenericImage<VoxelType>*[n];
    for(int i = 0; i < _numberOfImages; i++){
        if (image[i] != NULL) {
            _input[i] = image[i];
        } else {
            cerr << "irtkImageGraphCut::SetInput: "<< i <<"th Input is not an image\n";
            exit(1);
        }
    }
    _labels = l;
    _weight = new irtkRealImage*[l];
    for(int i = 0; i < _labels; i++){
        if (weight[i] != NULL) {
            _weight[i] = weight[i];
        } else {
            cerr << "irtkImageGraphCut::SetInput: "<< i <<"th Weight is not an image\n";
            exit(1);
        }
    }
    _dx = _input[0]->GetXSize();
    _dy = _input[0]->GetYSize();
    _dz = _input[0]->GetZSize();
    _dt = _input[0]->GetTSize();
}

template <class VoxelType> void irtkImageGraphCut<VoxelType>::SetInput(irtkGenericImage<VoxelType> *image,int l, irtkRealImage **weight)
{
    _numberOfImages = 1;
    _numberOfWeights = 1.0;
    _input = new irtkGenericImage<VoxelType>*[1];
    for(int i = 0; i < _numberOfImages; i++){
        if (image != NULL) {
            _input[i] = image;
        } else {
            cerr << "irtkImageGraphCut::SetInput: "<< i <<"th Input is not an image\n";
            exit(1);
        }
    }
    _labels = l;
    _weight = new irtkRealImage*[l];
    for(int i = 0; i < _labels; i++){
        if (weight[i] != NULL) {
            _weight[i] = weight[i];
        } else {
            cerr << "irtkImageGraphCut::SetInput: "<< i <<"th Weight is not an image\n";
            exit(1);
        }
    }
    _dx = _input[0]->GetXSize();
    _dy = _input[0]->GetYSize();
    _dz = _input[0]->GetZSize();
    _dt = _input[0]->GetTSize();
}

template <class VoxelType> void irtkImageGraphCut<VoxelType>::SetOutput(irtkGenericImage<VoxelType> *image)
{
    if (image != NULL) {
        _output = image;
    } else {
        cerr << "irtkImageGraphCut::SetOutput: Output is not an image\n";
        exit(1);
    }
}

template <class VoxelType> void irtkImageGraphCut<VoxelType>::Initialize()
{
    // Check inputs and outputs
    if (_input == NULL && _output == NULL) {
        cerr << this->NameOfClass() << "::Run: has no input nor output" << endl;
        exit(1);
    }
    for(int i = 0; i < _numberOfImages; i++){
        if (_input[i] == NULL) {
            cerr << this->NameOfClass() << "::Run: has no input: "<<i << endl;
            exit(1);
        }
    }
    // Check inputs and outputs
    if (_weight == NULL) {
        cerr << this->NameOfClass() << "::Run: has no weight" << endl;
        exit(1);
    }
    for(int i = 0; i < _labels; i++){
        if (_weight[i] == NULL) {
            cerr << this->NameOfClass() << "::Run: has no weight: "<<i << endl;
            exit(1);
        }
    }

    if (_input[0]->GetX() != _output->GetX() || _input[0]->GetY() != _output->GetY() 
        || _input[0]->GetZ() != _output->GetZ() || _input[0]->GetT() != _output->GetT()
        || _input[0]->GetX() != _weight[0]->GetX() || _input[0]->GetY() != _weight[0]->GetY() 
        || _input[0]->GetZ() != _weight[0]->GetZ() || _input[0]->GetT() != _weight[0]->GetT())
    {
        cerr << this->NameOfClass() << "::Run: Image Demensions does not correspond" << endl;
        exit(1);
    }

    for(int i = 0; i < _numberOfImages; i++){
        if (_input[0]->GetX() != _input[i]->GetX() || _input[0]->GetY() != _input[i]->GetY() 
            || _input[0]->GetZ() != _input[i]->GetZ() || _input[0]->GetT() != _input[i]->GetT())
        {
            cerr << this->NameOfClass() << "::Run: Image: "<<i<<" Demensions does not correspond" << endl;
            exit(1);
        }
    }

    for(int i = 0; i < _labels; i++){
        if (_weight[0]->GetX() != _weight[i]->GetX() || _weight[0]->GetY() != _weight[i]->GetY() 
            || _weight[0]->GetZ() != _weight[i]->GetZ() || _weight[0]->GetT() != _weight[i]->GetT())
        {
            cerr << this->NameOfClass() << "::Run: Image: "<<i<<" Demensions does not correspond" << endl;
            exit(1);
        }
    }
}

template <class VoxelType> void irtkImageGraphCut<VoxelType>::Finalize()
{

}

template <class VoxelType> void irtkImageGraphCut<VoxelType>::AddBoundaryTerm(GCoptimizationGeneralGraph *graph, int count, 
    int i,int j, int k, int l,
    int xoff, int yoff, int zoff, int toff,double divide)
{
    double tmpweight,weight;
    weight = 0;
    int n;
    for(n=0;n<_numberOfImages;n++){
        tmpweight = abs(_input[n]->GetAsDouble(i,j,k,l) - _input[n]->GetAsDouble(i+xoff,j+yoff,k+zoff,l+toff));
        weight += 1/(log(10+tmpweight)/log(double(10)));
    }
    weight = weight / divide / _numberOfImages;
    graph->setNeighbors(count,_input[0]->GetImageAttributes().LatticeToIndex(i+xoff,j+yoff,k+zoff,l+toff), round(weight*1000.0));
}

template <class VoxelType> void irtkImageGraphCut<VoxelType>::Run(double lambda, int connect)
{
    int i,j,k,l,count,n;
    count = 0;

    // Do the initial set up
    this->Initialize();

    try{
        GCoptimizationGeneralGraph *graph = new GCoptimizationGeneralGraph(_input[0]->GetNumberOfVoxels(),_labels);
        int *smooth = new int[_labels*_labels];
        for ( int l1 = 0; l1 < _labels; l1++ )
            for (int l2 = 0; l2 < _labels; l2++ )
                if(l1 == l2){
                    smooth[l1+l2*_labels] = 0;
                }else{
                    smooth[l1+l2*_labels] = 1;
                }

        graph->setSmoothCost(smooth);

        int *datacost = new int[_input[0]->GetNumberOfVoxels()*_labels];
        for (l = 0; l < _input[0]->GetT(); l++) {
            for (k = 0; k < _input[0]->GetZ(); k++) {
                for (j = 0; j < _input[0]->GetY(); j++) {
                    for (i = 0; i < _input[0]->GetX(); i++) {
                        count = _input[0]->GetImageAttributes().LatticeToIndex(i,j,k,l);
                        //evaluate weights
                        for(n = 0; n < _labels; n++){
                            datacost[count*_labels+n] = round(lambda*1000.0*(1.0-_weight[n]->GetAsDouble(i,j,k,l)));
                        }
                    }
                }
            }
        }
        graph->setDataCost(datacost);

        for (l = 0; l < _input[0]->GetT(); l++) {
            for (k = 0; k < _input[0]->GetZ(); k++) {
                for (j = 0; j < _input[0]->GetY(); j++) {
                    for (i = 0; i < _input[0]->GetX(); i++) {
                        count = _input[0]->GetImageAttributes().LatticeToIndex(i,j,k,l);

                        //add edges (boundary term)
                        if(i < _input[0]->GetX() - 1 && this->_mode > 0){
                            //evaluate boundary weight from inputs
                            AddBoundaryTerm(graph,count,i,j,k,l,1,0,0,0,_dx);
                        }
                        if(j < _input[0]->GetY() - 1 && this->_mode > 1){
                            AddBoundaryTerm(graph,count,i,j,k,l,0,1,0,0,_dy);
                        }
                        if(k < _input[0]->GetZ() - 1 && this->_mode > 2){
                            AddBoundaryTerm(graph,count,i,j,k,l,0,0,1,0,_dz);
                        }
                        if(l < _input[0]->GetT() - 1 && this->_mode > 3){
                            AddBoundaryTerm(graph,count,i,j,k,l,0,0,0,1,_dt);
                        }
                        //add more edges if connection is set to cubic
                        if(connect)
                        {
                            //3.110
                            if(j < _input[0]->GetY() - 1 && i < _input[0]->GetX() - 1 && this->_mode > 1)
                                AddBoundaryTerm(graph,count,i,j,k,l,1,1,0,0,sqrt(_dx*_dx + _dy*_dy));
                            //4.-110
                            if(j < _input[0]->GetY() - 1 && i > 0 && this->_mode > 1)
                                AddBoundaryTerm(graph,count,i,j,k,l,-1,1,0,0,sqrt(_dx*_dx + _dy*_dy));

                            //6.111
                            if(k < _input[0]->GetZ() - 1&& j < _input[0]->GetY() - 1 && i < _input[0]->GetX() - 1&& this->_mode > 2)
                                AddBoundaryTerm(graph,count,i,j,k,l,1,1,1,0,sqrt(_dx*_dx + _dy*_dy + _dz*_dz));
                            //7.101
                            if(k < _input[0]->GetZ() - 1 && i < _input[0]->GetX() - 1 && this->_mode > 2)
                                AddBoundaryTerm(graph,count,i,j,k,l,1,0,1,0,sqrt(_dx*_dx + _dz*_dz));
                            //8.1-11
                            if(k < _input[0]->GetZ() - 1&& j > 0 && i < _input[0]->GetX() - 1&& this->_mode > 2)
                                AddBoundaryTerm(graph,count,i,j,k,l,1,-1,1,0,sqrt(_dx*_dx + _dy*_dy + _dz*_dz));
                            //9.011
                            if(k < _input[0]->GetZ() - 1 && j < _input[0]->GetY() - 1&& this->_mode > 2)
                                AddBoundaryTerm(graph,count,i,j,k,l,0,1,1,0,sqrt(_dy*_dy + _dz*_dz));
                            //10.0-11
                            if(k < _input[0]->GetZ() - 1 && j >0 && this->_mode > 2)
                                AddBoundaryTerm(graph,count,i,j,k,l,0,-1,1,0,sqrt(_dy*_dy + _dz*_dz));
                            //11.-111
                            if(k < _input[0]->GetZ() - 1&& j < _input[0]->GetY() - 1 && i > 0 && this->_mode > 2)
                                AddBoundaryTerm(graph,count,i,j,k,l,-1,1,1,0,sqrt(_dx*_dx + _dy*_dy + _dz*_dz));
                            //12.-101
                            if(k < _input[0]->GetZ() - 1 && i > 0 && this->_mode > 2)
                                AddBoundaryTerm(graph,count,i,j,k,l,-1,0,1,0,sqrt(_dx*_dx + _dz*_dz));
                            //13.-1-11
                            if(k < _input[0]->GetZ() - 1&& j > 0 && i > 0 && this->_mode > 2)
                                AddBoundaryTerm(graph,count,i,j,k,l,-1,-1,1,0,sqrt(_dx*_dx + _dy*_dy + _dz*_dz));
                        }
                    }
                }
            }
        }

        int term;
        graph->setLabelCost(1);
        graph->expansion();
        for (l = 0; l < _output->GetT(); l++) {
            for (k = 0; k < _output->GetZ(); k++) {
                for (j = 0; j < _output->GetY(); j++) {
                    for (i = 0; i < _output->GetX(); i++) {
                        term = 0;
                        term = graph->whatLabel(_output->GetImageAttributes().LatticeToIndex(i,j,k,l));
                        _output->PutAsDouble(i,j,k,l,term);
                    }
                }
            }
        }
        delete graph;
        delete []smooth;
        delete []datacost;
    }
    catch (GCException e){
        e.Report();
    }

    // Do the final cleaning up
    this->Finalize();
}

template class irtkImageGraphCut<unsigned char>;
template class irtkImageGraphCut<short>;
template class irtkImageGraphCut<unsigned short>;
template class irtkImageGraphCut<float>;
template class irtkImageGraphCut<double>;
