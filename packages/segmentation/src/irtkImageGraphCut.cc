/*=========================================================================

  Library   : Image Registration Toolkit (IRTK)
  Module    : $Id: irtkImageGraphCut.cc 101 2009-11-04 16:27:54Z dr $
  Copyright : Imperial College, Department of Computing
              Visual Information Processing (VIP), 2008 onwards
  Date      : $Date: 2009-11-04 16:27:54 +0000 (‰∏? 04 ÂçÅ‰∏ÄÊú?2009) $
  Version   : $Revision: 101 $
  Changes   : $Author: dr $

=========================================================================*/

#include <irtkImageGraphCut.h>
#include <irtkGradientImage.h>

void doIt2(char * errMsg)
{
cout << "Error: " << errMsg << endl;
}

template <class VoxelType> irtkImageGraphCut<VoxelType>::irtkImageGraphCut()
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
  _numberOfWeights = 0.0;
}

template <class VoxelType> irtkImageGraphCut<VoxelType>::~irtkImageGraphCut()
{
  // Set in- and outputs
  //if(_input)
	  //delete []_input;
  _input  = NULL;
  _output = NULL;
  /// source weight image
  _sourceweight = NULL;
  /// sink weight image
  _sinkweight = NULL;
}

template <class VoxelType> void irtkImageGraphCut<VoxelType>::SetInput(int n, irtkGenericImage<VoxelType> **image,irtkGenericImage<VoxelType> *sourceweight,irtkGenericImage<VoxelType> *sinkweight)
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
  if (sourceweight != NULL) {
    _sourceweight = sourceweight;
  } else {
    cerr << "irtkImageGraphCut::SetInput: Input is not an image\n";
    exit(1);
  }
  if (sinkweight != NULL) {
    _sinkweight = sinkweight;
  } else {
    cerr << "irtkImageGraphCut::SetInput: Input is not an image\n";
    exit(1);
  }
  _dx = _input[0]->GetXSize();
  _dy = _input[0]->GetYSize();
  _dz = _input[0]->GetZSize();
  _dt = _input[0]->GetTSize();
}

template <class VoxelType> void irtkImageGraphCut<VoxelType>::SetInput(irtkGenericImage<VoxelType> *image,irtkGenericImage<VoxelType> *sourceweight,irtkGenericImage<VoxelType> *sinkweight)
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
  if (sourceweight != NULL) {
    _sourceweight = sourceweight;
  } else {
    cerr << "irtkImageGraphCut::SetInput: Input is not an image\n";
    exit(1);
  }
  if (sinkweight != NULL) {
    _sinkweight = sinkweight;
  } else {
    cerr << "irtkImageGraphCut::SetInput: Input is not an image\n";
    exit(1);
  }
  _dx = _input[0]->GetXSize();
  _dy = _input[0]->GetYSize();
  _dz = _input[0]->GetZSize();
  _dt = _input[0]->GetTSize();
}

template <class VoxelType> void irtkImageGraphCut<VoxelType>::AddWeight(irtkGenericImage<VoxelType> *sourceweight,irtkGenericImage<VoxelType> *sinkweight, double weight)
{
	int i,j,k,l;
	_numberOfWeights += weight;
  if (sourceweight == NULL) {
    cerr << "irtkImageGraphCut::SetInput: Input is not an image\n";
    exit(1);
  }
  if (sinkweight == NULL) {
    cerr << "irtkImageGraphCut::SetInput: Input is not an image\n";
    exit(1);
  }
 
  if (_sourceweight->GetX() != sourceweight->GetX() || _sourceweight->GetY() != sourceweight->GetY() 
	  || _sourceweight->GetZ() != sourceweight->GetZ() || _sourceweight->GetT() != sourceweight->GetT())
  {
	  cerr << this->NameOfClass() << "::Run: SourceWeight: Demensions does not correspond" << endl;
	  exit(1);
  }
  if (_sinkweight->GetX() != sinkweight->GetX() || _sinkweight->GetY() != sinkweight->GetY() 
	  || _sinkweight->GetZ() != sinkweight->GetZ() || _sinkweight->GetT() != sinkweight->GetT())
  {
	  cerr << this->NameOfClass() << "::Run: SinkWeight: Demensions does not correspond" << endl;
	  exit(1);
  }
  for (l = 0; l < sourceweight->GetT(); l++) {
	  for (k = 0; k < sourceweight->GetZ(); k++) {
		for (j = 0; j < sourceweight->GetY(); j++) {
			for (i = 0; i < sourceweight->GetX(); i++) {
				_sourceweight->PutAsDouble(i,j,k,l,((_numberOfWeights - weight)*_sourceweight->GetAsDouble(i,j,k,l) 
					+ weight*sourceweight->GetAsDouble(i,j,k,l))/_numberOfWeights);
			}
		}
	  }
  }
  for (l = 0; l < sinkweight->GetT(); l++) {
	  for (k = 0; k < sinkweight->GetZ(); k++) {
		for (j = 0; j < sinkweight->GetY(); j++) {
			for (i = 0; i < sinkweight->GetX(); i++) {
				_sinkweight->PutAsDouble(i,j,k,l,((_numberOfWeights - weight)*_sinkweight->GetAsDouble(i,j,k,l) 
					+ weight*sinkweight->GetAsDouble(i,j,k,l))/_numberOfWeights);
			}
		}
	  }
  }
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
  if (_sinkweight == NULL && _sourceweight == NULL) {
    cerr << this->NameOfClass() << "::Run: has no weight" << endl;
    exit(1);
  }

  if (_input[0]->GetX() != _output->GetX() || _input[0]->GetY() != _output->GetY() 
	  || _input[0]->GetZ() != _output->GetZ() || _input[0]->GetT() != _output->GetT()
	  || _input[0]->GetX() != _sinkweight->GetX() || _input[0]->GetY() != _sinkweight->GetY() 
	  || _input[0]->GetZ() != _sinkweight->GetZ() || _input[0]->GetT() != _sinkweight->GetT()
	  || _input[0]->GetX() != _sourceweight->GetX() || _input[0]->GetY() != _sourceweight->GetY() 
	  || _input[0]->GetZ() != _sourceweight->GetZ() || _input[0]->GetT() != _sourceweight->GetT())
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
}

template <class VoxelType> void irtkImageGraphCut<VoxelType>::Finalize()
{
  
}

template <class VoxelType> void irtkImageGraphCut<VoxelType>::AddBoundaryTerm(Graph<double, double, double>& graph, int count, 
																			  int i,int j, int k, int l,
																			  int xoff, int yoff, int zoff, int toff,double divide)
{
	double tmpweight,weight;
	weight = 0;
	int n;
	for(n=0;n<_numberOfImages;n++){
		tmpweight = abs(_input[n]->GetAsDouble(i,j,k,l) - _input[n]->GetAsDouble(i+xoff,j+yoff,k+zoff,l+toff));
		weight += 1/(log(1+pow(tmpweight,2.0))/log(10.0)+0.000001);
	}
	weight = weight / divide / _numberOfImages;
	graph.add_edge(count,_input[0]->GetImageAttributes().LatticeToIndex(i+xoff,j+yoff,k+zoff,l+toff), weight, weight);
}

template <class VoxelType> void irtkImageGraphCut<VoxelType>::Run(double lambda, int connect)
{
  int i,j,k,l,count;
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

  void (*errFunction)(char *) = doIt2;
  Graph<double, double, double> graph (100, 600, errFunction);

  graph.add_node(_input[0]->GetNumberOfVoxels());

  for (l = 0; l < _input[0]->GetT(); l++) {
	  for (k = 0; k < _input[0]->GetZ(); k++) {
		for (j = 0; j < _input[0]->GetY(); j++) {
			for (i = 0; i < _input[0]->GetX(); i++) {
				count = _input[0]->GetImageAttributes().LatticeToIndex(i,j,k,l);
				//evaluate weights
				sinkweight = _sinkweight->GetAsDouble(i,j,k,l);
				sourceweight = _sourceweight->GetAsDouble(i,j,k,l);

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

				//add weight (region term)
				graph.add_tweights(count,lambda*sourceweight,lambda*sinkweight);
			}
		}
	}
  }

  double flow;
  int term;
  flow = graph.maxflow();
  for (l = 0; l < _output->GetT(); l++) {
	  for (k = 0; k < _output->GetZ(); k++) {
		for (j = 0; j < _output->GetY(); j++) {
			for (i = 0; i < _output->GetX(); i++) {
				term = 0;
				term = graph.what_segment(_output->GetImageAttributes().LatticeToIndex(i,j,k,l));
				_output->PutAsDouble(i,j,k,l,term);
			}
		}
	}
  }

  // Do the final cleaning up
  this->Finalize();
}

template class irtkImageGraphCut<unsigned char>;
template class irtkImageGraphCut<short>;
template class irtkImageGraphCut<unsigned short>;
template class irtkImageGraphCut<float>;
template class irtkImageGraphCut<double>;
