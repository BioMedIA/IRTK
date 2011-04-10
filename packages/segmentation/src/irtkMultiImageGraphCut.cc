/*=========================================================================

Library   : Image Registration Toolkit (IRTK)
Module    : $Id: irtkMultiImageGraphCut.cc 101 2009-11-04 16:27:54Z dr $
Copyright : Imperial College, Department of Computing
Visual Information Processing (VIP), 2008 onwards
Date      : $Date: 2009-11-04 16:27:54 +0000 (‰∏? 04 ÂçÅ‰∏ÄÊú?2009) $
Version   : $Revision: 101 $
Changes   : $Author: dr $

=========================================================================*/

#include <irtkMultiImageGraphCut.h>
#include <irtkGradientImage.h>

void mdoIt2(char * errMsg)
{
	cout << "Error: " << errMsg << endl;
}

template <class VoxelType> irtkMultiImageGraphCut<VoxelType>::irtkMultiImageGraphCut()
{
	// Set in- and outputs
	_input  = NULL;
	_output = NULL;
	/// source weight image
	_sourceweight = NULL;
	/// sink weight image
	_sinkweight = NULL;
	_mode = 4;
	_numberOfImages = 0;
	_dx = NULL;
	_dy = NULL;
	_dz = NULL;
	_imageoffset = NULL;
}

template <class VoxelType> irtkMultiImageGraphCut<VoxelType>::~irtkMultiImageGraphCut()
{
	if(_dx) delete []_dx;
	if(_dy) delete []_dy;
	if(_dz) delete []_dz;
	_dx = NULL;
	_dy = NULL;
	_dz = NULL;
	//if(_output) delete []_output;
	if(_sourceweight) delete []_sourceweight;
	if(_sinkweight) delete []_sinkweight;
	if(_imageoffset) delete []_imageoffset;
	if(_input) delete []_input;
	_input  = NULL;
	_output = NULL;
	_sourceweight = NULL;
	_sinkweight = NULL;
}

template <class VoxelType> void irtkMultiImageGraphCut<VoxelType>::SetInput(int n, irtkGenericImage<VoxelType> **image,irtkGenericImage<VoxelType> **sourceweight,irtkGenericImage<VoxelType> **sinkweight)
{
	_numberOfImages = n;
	_input = new irtkGenericImage<VoxelType>*[n];
	_dx = new double[n];
	_dy = new double[n];
	_dz = new double[n];
	_imageoffset = new double[n];
	_totalVoxel = 0;
	for(int i = 0; i < _numberOfImages; i++){
		if (image[i] != NULL) {
			_input[i] = image[i];
			_dx[i] = _input[i]->GetXSize();
			_dy[i] = _input[i]->GetYSize();
			_dz[i] = _input[i]->GetZSize();
		} else {
			cerr << "irtkMultiImageGraphCut::SetInput: "<< i <<"th Input is not an image\n";
			exit(1);
		}
		_imageoffset[i] = _totalVoxel;
		_totalVoxel += _input[i]->GetNumberOfVoxels();
	}
	_sourceweight = new irtkGenericImage<VoxelType>*[n];
	for(int i = 0; i < _numberOfImages; i++){
		if (sourceweight[i] != NULL) {
			_sourceweight[i] = sourceweight[i];
		} else {
			cerr << "irtkMultiImageGraphCut::SetInput: "<< i <<"th Sourceweight is not an image\n";
			exit(1);
		}
	}
	_sinkweight = new irtkGenericImage<VoxelType>*[n];
	for(int i = 0; i < _numberOfImages; i++){
		if (sinkweight[i] != NULL) {
			_sinkweight[i] = sinkweight[i];
		} else {
			cerr << "irtkMultiImageGraphCut::SetInput: "<< i <<"th SinkWeight is not an image\n";
			exit(1);
		}
	}
	_dt = 1;
}

template <class VoxelType> void irtkMultiImageGraphCut<VoxelType>::SetOutput(irtkGenericImage<VoxelType> **image)
{
	if (image != NULL) {
		_output = image;
	} else {
		cerr << "irtkMultiImageGraphCut::SetOutput: Output is not an image\n";
		exit(1);
	}
}

template <class VoxelType> void irtkMultiImageGraphCut<VoxelType>::Initialize()
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
	if (_sinkweight == NULL || _sourceweight == NULL) {
		cerr << this->NameOfClass() << "::Run: has no weight" << endl;
		exit(1);
	}

	for(int i = 0; i < _numberOfImages; i++){
		if (_sinkweight[i] == NULL || _sourceweight[i] == NULL) {
			cerr << this->NameOfClass() << "::Run: has no weight: "<<i << endl;
			exit(1);
		}
	}

	for(int i = 0; i < _numberOfImages; i++){
		if (_input[i]->GetX() != _output[i]->GetX() || _input[i]->GetY() != _output[i]->GetY() 
			|| _input[i]->GetZ() != _output[i]->GetZ() || _input[i]->GetT() != _output[i]->GetT()
			|| _input[i]->GetX() != _sinkweight[i]->GetX() || _input[i]->GetY() != _sinkweight[i]->GetY() 
			|| _input[i]->GetZ() != _sinkweight[i]->GetZ() || _input[i]->GetT() != _sinkweight[i]->GetT()
			|| _input[i]->GetX() != _sourceweight[i]->GetX() || _input[i]->GetY() != _sourceweight[i]->GetY() 
			|| _input[i]->GetZ() != _sourceweight[i]->GetZ() || _input[i]->GetT() != _sourceweight[i]->GetT())
		{
			cerr << this->NameOfClass() << "::Run: Image Demensions does not correspond" << endl;
			exit(1);
		}
	}
}

template <class VoxelType> void irtkMultiImageGraphCut<VoxelType>::Finalize()
{

}

template <class VoxelType> void irtkMultiImageGraphCut<VoxelType>::AddBoundaryTerm(Graph<double, double, double>& graph, int count, 
	int i,int j, int k, int l,  int n,
	int xoff, int yoff, int zoff, int toff,double divide)
{
	double tmpweight,weight;
	weight = 0;
	tmpweight = abs(_input[n]->GetAsDouble(i,j,k,l) - _input[n]->GetAsDouble(i+xoff,j+yoff,k+zoff,l+toff));
	weight += 1/(log(1+pow(tmpweight,2.0))/log(10.0)+0.000001);
	weight = weight / divide;
	graph.add_edge(count,_imageoffset[n]+_input[n]->GetImageAttributes().LatticeToIndex(i+xoff,j+yoff,k+zoff,l+toff), weight, weight);
}

template <class VoxelType> void irtkMultiImageGraphCut<VoxelType>::AddImageTerm(Graph<double, double, double>& graph, int count, 
	int count2,double divide)
{
	double weight;
	weight = 1.0 / divide;
	graph.add_edge(count,count2, weight, weight);
}

template <class VoxelType> void irtkMultiImageGraphCut<VoxelType>::Run(double lambda, int connect)
{
	int i,j,k,l,n,count;
	double sinkweight,sourceweight,weight,tmpweight;
	//irtkGenericImage<VoxelType> **tmpedge;
	//VoxelType *intsd,*edgesd;
	count = 0;
	sinkweight = 0;
	sourceweight = 0;
	weight = 0;
	tmpweight = 0;

	// Do the initial set up
	this->Initialize();

	void (*errFunction)(char *) = mdoIt2;
	Graph<double, double, double> graph (100, 600, errFunction);

	graph.add_node(_totalVoxel);
	for (n = 0; n < _numberOfImages; n++){
		for (l = 0; l < _input[n]->GetT(); l++) {
			for (k = 0; k < _input[n]->GetZ(); k++) {
				for (j = 0; j < _input[n]->GetY(); j++) {
					for (i = 0; i < _input[n]->GetX(); i++) {
						count = _input[n]->GetImageAttributes().LatticeToIndex(i,j,k,l);
						//evaluate weights
						sinkweight = _sinkweight[n]->GetAsDouble(i,j,k,l);
						sourceweight = _sourceweight[n]->GetAsDouble(i,j,k,l);

						//add edges (boundary term)
						if(i < _input[n]->GetX() - 1 && this->_mode > 0){
							//evaluate boundary weight from inputs
							AddBoundaryTerm(graph,_imageoffset[n]+count,i,j,k,l,n,1,0,0,0,_dx[n]);
						}
						if(j < _input[n]->GetY() - 1 && this->_mode > 1){
							AddBoundaryTerm(graph,_imageoffset[n]+count,i,j,k,l,n,0,1,0,0,_dy[n]);
						}
						if(k < _input[n]->GetZ() - 1 && this->_mode > 2){
							AddBoundaryTerm(graph,_imageoffset[n]+count,i,j,k,l,n,0,0,1,0,_dz[n]);
						}
						if(l < _input[n]->GetT() - 1 && this->_mode > 3){
							AddBoundaryTerm(graph,_imageoffset[n]+count,i,j,k,l,n,0,0,0,1,_dt);
						}
						//add more edges if connection is set to cubic
						if(connect)
						{
							//3.110
							if(j < _input[n]->GetY() - 1 && i < _input[n]->GetX() - 1 && this->_mode > 1)
								AddBoundaryTerm(graph,_imageoffset[n]+count,i,j,k,l,n,1,1,0,0,sqrt(_dx[n]*_dx[n] + _dy[n]*_dy[n]));
							//4.-110
							if(j < _input[n]->GetY() - 1 && i > 0 && this->_mode > 1)
								AddBoundaryTerm(graph,_imageoffset[n]+count,i,j,k,l,n,-1,1,0,0,sqrt(_dx[n]*_dx[n] + _dy[n]*_dy[n]));

							//6.111
							if(k < _input[n]->GetZ() - 1&& j < _input[n]->GetY() - 1 && i < _input[n]->GetX() - 1&& this->_mode > 2)
								AddBoundaryTerm(graph,_imageoffset[n]+count,i,j,k,l,n,1,1,1,0,sqrt(_dx[n]*_dx[n] + _dy[n]*_dy[n] + _dz[n]*_dz[n]));
							//7.101
							if(k < _input[n]->GetZ() - 1 && i < _input[n]->GetX() - 1 && this->_mode > 2)
								AddBoundaryTerm(graph,_imageoffset[n]+count,i,j,k,l,n,1,0,1,0,sqrt(_dx[n]*_dx[n] + _dz[n]*_dz[n]));
							//8.1-11
							if(k < _input[n]->GetZ() - 1&& j > 0 && i < _input[n]->GetX() - 1&& this->_mode > 2)
								AddBoundaryTerm(graph,_imageoffset[n]+count,i,j,k,l,n,1,-1,1,0,sqrt(_dx[n]*_dx[n] + _dy[n]*_dy[n] + _dz[n]*_dz[n]));
							//9.011
							if(k < _input[n]->GetZ() - 1 && j < _input[n]->GetY() - 1&& this->_mode > 2)
								AddBoundaryTerm(graph,_imageoffset[n]+count,i,j,k,l,n,0,1,1,0,sqrt(_dy[n]*_dy[n] + _dz[n]*_dz[n]));
							//10.0-11
							if(k < _input[n]->GetZ() - 1 && j >0 && this->_mode > 2)
								AddBoundaryTerm(graph,_imageoffset[n]+count,i,j,k,l,n,0,-1,1,0,sqrt(_dy[n]*_dy[n] + _dz[n]*_dz[n]));
							//11.-111
							if(k < _input[n]->GetZ() - 1&& j < _input[n]->GetY() - 1 && i > 0 && this->_mode > 2)
								AddBoundaryTerm(graph,_imageoffset[n]+count,i,j,k,l,n,-1,1,1,0,sqrt(_dx[n]*_dx[n] + _dy[n]*_dy[n] + _dz[n]*_dz[n]));
							//12.-101
							if(k < _input[n]->GetZ() - 1 && i > 0 && this->_mode > 2)
								AddBoundaryTerm(graph,_imageoffset[n]+count,i,j,k,l,n,-1,0,1,0,sqrt(_dx[n]*_dx[n] + _dz[n]*_dz[n]));
							//13.-1-11
							if(k < _input[n]->GetZ() - 1&& j > 0 && i > 0 && this->_mode > 2)
								AddBoundaryTerm(graph,_imageoffset[n]+count,i,j,k,l,n,-1,-1,1,0,sqrt(_dx[n]*_dx[n] + _dy[n]*_dy[n] + _dz[n]*_dz[n]));
						}

						//add weight (region term)
						graph.add_tweights(_imageoffset[n]+count,lambda*sourceweight,lambda*sinkweight);
					}
				}
			}
		}
	}

	//Add multiple image terms
	for (n = 0; n < _numberOfImages; n++){
		for (l = 0; l < _input[n]->GetT(); l++) {
			for (k = 0; k < _input[n]->GetZ(); k++) {
				for (j = 0; j < _input[n]->GetY(); j++) {
					for (i = 0; i < _input[n]->GetX(); i++) {
						count = _input[n]->GetImageAttributes().LatticeToIndex(i,j,k,l);
						double x[3],y[3],z[3],sdx,sdy,sdz,t;
						for (int m = 0; m < 3; m ++){
							x[m] = i + 0.5*(m-1);
							y[m] = j + 0.5*(m-1);
							z[m] = k + 0.5*(m-1);
							_input[n]->ImageToWorld(x[m],y[m],z[m]);
						}
						t = l;
						_input[n]->ImageToTime(t);
						sdx = _input[n]->GetXSize();
						sdy = _input[n]->GetYSize();
						sdz = _input[n]->GetZSize();
						for(int m = 0; m < _numberOfImages; m++){
							if(n != m ){
								double tx[3],ty[3],tz[3],tdx,tdy,tdz;
								int count2;
								for(int p = 0; p < 3; p++){
									tx[p] = x[p];
									ty[p] = y[p];
									tz[p] = z[p];
									_input[m]->WorldToImage(tx[p],ty[p],tz[p]);
								}
								_input[m]->TimeToImage(t);
								tdx = _input[m]->GetXSize();
								tdy = _input[m]->GetYSize();
								tdz = _input[m]->GetZSize();
								int ti,tj,tk;
								ti = round(tx[1]); tj = round(ty[1]); tk = round(tz[1]);
								if(t >= 0 && t < _input[m]->GetT()
									&&ti >= 0 && ti < _input[m]->GetX()
									&& tj >= 0 && tj < _input[m]->GetY()
									&& tk >= 0 && tk < _input[m]->GetZ()){
										//evaluate dice metric
										double dx,dy,dz,ds;
										double px[4],py[4],pz[4];
										px[0] = tx[0]; px[1] = tx[2];
										px[2] = double(ti) - 0.5; px[3] = double(ti) + 0.5;
										py[0] = ty[0]; py[1] = ty[2];
										py[2] = double(tj) - 0.5; py[3] = double(tj) + 0.5;
										pz[0] = tz[0]; pz[1] = tz[2];
										pz[2] = double(tk) - 0.5; pz[3] = double(tk) + 0.5;
										sort(px,px+4);
										sort(py,py+4);
										sort(pz,pz+4);
										dx = (px[2]-px[1])/(px[3]-px[0]);
										dy = (py[2]-py[1])/(py[3]-py[0]);
										dz = (pz[2]-pz[1])/(pz[3]-pz[0]);
										ds = dx*dy*dz;
										sdx*sdy*sdz>tdx*tdy*tdz?
											ds = ds*(sdx*sdy*sdz)/(tdx*tdy*tdz):
										    ds = ds*(tdx*tdy*tdz)/(sdx*sdy*sdz);
								count2 = _input[m]->GetImageAttributes().LatticeToIndex(ti,tj,tk,t);
								AddImageTerm(graph,_imageoffset[n]+count,_imageoffset[m]+count2,_dt/ds);
								}
							}
						}
					}
				}
			}
		}
	}

	double flow;
	int term;
	flow = graph.maxflow();
	for (n = 0; n < _numberOfImages; n++){
		for (l = 0; l < _output[n]->GetT(); l++) {
			for (k = 0; k < _output[n]->GetZ(); k++) {
				for (j = 0; j < _output[n]->GetY(); j++) {
					for (i = 0; i < _output[n]->GetX(); i++) {
						term = 0;
						term = graph.what_segment(_imageoffset[n] + _output[n]->GetImageAttributes().LatticeToIndex(i,j,k,l));
						_output[n]->PutAsDouble(i,j,k,l,term);
					}
				}
			}
		}
	}

	// Do the final cleaning up
	this->Finalize();
}

template class irtkMultiImageGraphCut<unsigned char>;
template class irtkMultiImageGraphCut<short>;
template class irtkMultiImageGraphCut<unsigned short>;
template class irtkMultiImageGraphCut<float>;
template class irtkMultiImageGraphCut<double>;
