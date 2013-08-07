#include <irtkGraphCutSegmentation_4D.h>
#include <iostream>
#include <irtkImage.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <string>
#include <sstream>
#include <irtkGradientImage.h>
#include <nr.h> 

irtkGraphCutSegmentation_4D::irtkGraphCutSegmentation_4D(int numTissues, irtkRealImage foregroundAtlas)
{

	_numTissues = numTissues;
	_foregroundAtlas = foregroundAtlas;
	_gaussians = new irtkGaussian[numTissues];
	_padding = 0;
	useMask = false;
}

irtkGraphCutSegmentation_4D::~irtkGraphCutSegmentation_4D()
{

	delete _segmentation;

}

void irtkGraphCutSegmentation_4D::SetTissuePriors(irtkRealImage ** priors){
	_tissuePriors = priors;
}

void irtkGraphCutSegmentation_4D::setMask(irtkGreyImage mask){
	_mask = mask;
	useMask = true;
}

void irtkGraphCutSegmentation_4D::SetInput(const irtkRealImage &input, irtkGenericImage<float> gradMagnitude)
{
	_input = input;
	_xDim = _input.GetX();
	_yDim = _input.GetY();
	_zDim = _input.GetZ();
	_segmentation = new irtkGenericImage<irtkGreyPixel>(_xDim, _yDim, _zDim, _input.GetT());
	_segmentation->Initialize(_input.GetImageAttributes());
	_gradMagnitude = gradMagnitude;
	int count = 0;
	int numberOfValues = 0;
	irtkRealPixel *ptr = _input.GetPointerToVoxels();
	for(int i = 0; i < _input.GetNumberOfVoxels(); i++){
				if(*ptr > _padding)
					numberOfValues++;
				ptr++;
	}
	float *gradMagnData = new float[1+numberOfValues];
	irtkRealPixel * gradPtr = _gradMagnitude.GetPointerToVoxels();
	irtkRealPixel * forePtr = _foregroundAtlas.GetPointerToVoxels();
	irtkGreyPixel * segPtr = _segmentation->GetPointerToVoxels();
	irtkRealPixel * inPtr = _input.GetPointerToVoxels();
	for(int i = 0; i < _segmentation->GetNumberOfVoxels(); i++){
		if(*inPtr <= _padding){
			*segPtr = _padding;
		}
		else{
			int tissue = 0;
			if(*forePtr > 0.5){
				tissue = 1;
			}
			*segPtr = tissue + 1;
			gradMagnData[1+count] = *gradPtr;
			count++;
		}
		gradPtr++;
		forePtr++;
		segPtr++;
		inPtr++;
	}
	sort(numberOfValues, gradMagnData);
	int index = 1 + (numberOfValues-1)/2;
	int median = gradMagnData[index];
	ptr = _gradMagnitude.GetPointerToVoxels();
	irtkRealPixel *inpPtr = _input.GetPointerToVoxels();
	count = 0;
	for(int i = 0; i < _gradMagnitude.GetNumberOfVoxels(); i++){
		if(*inpPtr > _padding){
			gradMagnData[count+1] = abs(*ptr - median);
			count++;
		}
		ptr++;
		inpPtr++;
	}

	sort(numberOfValues, gradMagnData);
	int MAD = gradMagnData[index];
	if(MAD == 0)
		MAD = 1;
	_sigma = 1.4826 * MAD;
	cout << "sigma " << _sigma << endl;
	delete gradMagnData;
}



void irtkGraphCutSegmentation_4D::GenerateGaussians(int tissueNum, double my, double s){
	_gaussians[tissueNum].Initialise(my, s);
}


void irtkGraphCutSegmentation_4D::SetMog(MOG m){
	mog = m;
}


void irtkGraphCutSegmentation_4D::SetMi(int m1, int m2, int m3)
{
	mi1 = m1;
	mi2 = m2;
	mi3 = m3;
}

void irtkGraphCutSegmentation_4D::setParameters(double l, double g, double c, double sigmaG)
{
	
	_sigmaG = .02;
	_mul=sigmaG;
	_c = c;
	_lambda = l;
	_gamma = g;
}

double irtkGraphCutSegmentation_4D::GetNeighbourWeight(irtkRealPixel *ptr1, irtkRealPixel *ptr2, irtkRealPixel *grad, double sigmaG, double tMul)
{
	double mul = 1;
	if(sigmaG > 1){
		sigmaG -= 1;
	}
	mul *= tMul;
	double rho = log(1 + 0.5 * pow(fabs(*ptr1 - *ptr2)/_sigma, 2 ));
	double wI = 1 / (1 + rho);
	double wB = exp(-*grad/sigmaG);
	double c;
	c = _c;
	double weight = c*wI + (1-c)*wB;
	weight *= mul;
	return weight;
}
double irtkGraphCutSegmentation_4D::GetTerminalWeight(irtkRealPixel *ptr, int label, int x, int y, int z)

{
	return GetTerminalWeight(ptr, label, x, y, z, 0);
}


double irtkGraphCutSegmentation_4D::GetTerminalWeight(irtkRealPixel *ptr, int label, int x, int y, int z, int t)
{
	double gamma = _gamma ;
	double weight = 0;
	double wa = 0;
	double foreP = _foregroundAtlas.GetAsDouble(x, y, z, t)/100;
	double p = 0;
	if(label-1 == _numTissues-1){
		p = 0;
		p = mog.evaluate(*ptr, x, y, z, t);
		weight = -log(p);
		double atlasVal = 1 - foreP;
		wa = -log(atlasVal);
	}
	else{
		p = _gaussians[label-1].Evaluate(*ptr);
		weight = -log(p);
		double atlasVal = foreP;
		wa = -log(atlasVal);
	}

	weight *= gamma;
	wa *= (1-gamma);
	weight += wa;
	weight = 1/weight;
	double lambda;
	lambda = _lambda;
	weight *= lambda;
	return weight;
}

int irtkGraphCutSegmentation_4D::GetWrongClassified()
{
	return _wrongClass;
}

void irtkGraphCutSegmentation_4D::Iterate()
{

	int iter = 0;
	int v = 0;
	do{
		v = 0;
		for(int i = 1; i <= _numTissues; i++){
			for(int j = 1; j <= _numTissues; j++){
				if(i < j && i <= _numTissues - 1){
					GenerateGraph(i, j);
					int res = UpdateSegmentation(i, j);
					delete _graph;
					delete _nodes;
					if(v == 0)
						v = res;
				}
			}
		}
		iter++;
	}while((v != 0) && iter < 10);
	cout << "iterations: " << iter << endl;
}

void doIt(char * errMsg)
{
	cout << "Error: " << errMsg << endl;
}

void irtkGraphCutSegmentation_4D::GenerateGraph(int label1, int label2)
{
	numNodes = 2;
	_nodes = new irtkGenericImage<irtkRealPixel>(_xDim, _yDim, _zDim, _input.GetT());
	_nodes->Initialize(_input.GetImageAttributes());
	irtkRealPixel *ptrNodes = _nodes->GetPointerToVoxels();
	for(int i = 0; i < _nodes->GetNumberOfVoxels(); i++){
		*ptrNodes=_padding;
		ptrNodes++;
	}
	void (*errFunction)(char *) = doIt;
	_graph = new Graph<double, double, double>(numNodes, numNodes * 6, errFunction);
	irtkGreyPixel *ptr;
	numNodes = 0;
	irtkRealPixel * inPtr = _input.GetPointerToVoxels();
	ptr = _segmentation->GetPointerToVoxels();
	ptrNodes = _nodes->GetPointerToVoxels();
	for(int i = 0; i < _input.GetNumberOfVoxels(); i++){
		if(*inPtr > _padding){
			if(*ptr == label1 || *ptr == label2){
				int nodeID = _graph->add_node(1);
				*ptrNodes = nodeID;
				numNodes++;
			}
		}
		inPtr++;
		ptr++;
		ptrNodes++;
	}
	for(int y = 0; y < _yDim; y++){
		for(int x = 0; x < _xDim; x++){
			for(int z = 0; z < _zDim; z++){
				for(int t = 0; t < _input.GetT(); t++){
					if((useMask && _mask.Get(x,y,z, t) > 0 && _nodes->Get(x, y, z, t) > _padding) || (! useMask && _nodes->Get(x, y, z, t) > _padding)){
						double tWeightL1 = GetTerminalWeight(_input.GetPointerToVoxels(x, y, z, t), label1,  x, y, z);
						double tWeightL2 = GetTerminalWeight(_input.GetPointerToVoxels(x, y, z, t), label2,  x, y, z);
						if(y<_yDim-1)
							if(_nodes->Get(x, y+1, z, t) > _padding){
								double weight = GetNeighbourWeight(_input.GetPointerToVoxels(x, y, z, t), _input.GetPointerToVoxels(x, y+1, z, t), _gradMagnitude.GetPointerToVoxels(x, y, z, t), _sigmaG, 1);
								_graph->add_edge(_nodes->Get(x, y, z, t), _nodes->Get(x, y+1, z, t), weight, weight);
							}
						if(x<_xDim-1)
							if(_nodes->Get(x+1, y, z, t) > _padding){
								double weight = GetNeighbourWeight(_input.GetPointerToVoxels(x, y, z, t), _input.GetPointerToVoxels(x+1, y, z, t), _gradMagnitude.GetPointerToVoxels(x, y, z, t), _sigmaG, 1);
								_graph->add_edge(_nodes->Get(x, y, z, t), _nodes->Get(x+1, y, z, t), weight, weight);
							} 
						if(z<_zDim-1)
							if(_nodes->Get(x, y, z+1, t) > _padding){
								double weight = GetNeighbourWeight(_input.GetPointerToVoxels(x, y, z, t), _input.GetPointerToVoxels(x, y, z+1, t), _gradMagnitude.GetPointerToVoxels(x, y, z, t), _sigmaG, 1);
								_graph->add_edge(_nodes->Get(x, y, z, t), _nodes->Get(x, y, z+1, t), weight, weight);
							}
						if(t <_input.GetT()-1)
							if(_nodes->Get(x, y, z, t+1) > _padding){
							 
								double weight = GetNeighbourWeight(_input.GetPointerToVoxels(x, y, z, t), _input.GetPointerToVoxels(x, y, z, t+1), _gradMagnitude.GetPointerToVoxels(x, y, z, t), _sigmaG+1, 0.0);
								_graph->add_edge(_nodes->Get(x, y, z, t), _nodes->Get(x, y, z, t+1), weight, weight);
							}
						_graph->add_tweights(_nodes->Get(x, y, z, t), tWeightL1, tWeightL2);
					}
				}
			}
		}
	}
}


int irtkGraphCutSegmentation_4D::UpdateSegmentation(int label1, int label2)
{
	double flow;
	flow = _graph->maxflow();
	int changed = 0;
	int wrongClass = 0;
	irtkGreyPixel * segPtr = _segmentation->GetPointerToVoxels();
	irtkRealPixel * inPtr =  _input.GetPointerToVoxels();
	irtkRealPixel * nodesPtr = _nodes->GetPointerToVoxels();
	for(int i = 0; i < _input.GetNumberOfVoxels(); i++){
		int segment;
		int newLabel;
		if(*nodesPtr <= _padding){
			segment = *segPtr;
			newLabel = segment;
		}
		else{
			if(*inPtr <= _padding){
				newLabel = *inPtr;
			}
			else{
				segment = _graph->what_segment(*nodesPtr);
				if(segment == 0){
					newLabel = label1;
				}
				else{
					newLabel = label2;
				}
			}
		}
	//	int corrVal = 0;
		if(*nodesPtr > _padding){
			if(*segPtr != newLabel){
				*segPtr = newLabel;
				changed++;
			}
		}
		segPtr++;
		inPtr++;
		nodesPtr++;
	}
	_wrongClass = wrongClass;
	return changed;
}

irtkRealImage irtkGraphCutSegmentation_4D::GetSegmentation()
{
	return *_segmentation;
}


////////////////////////////////////////////////////////////////////////////////////////////////////////




//MOG (Mixture of Gaussians) members:


/////////////////////////////////////////////////////////////////////////////////////////////////////////
MOG::MOG(){}

MOG::MOG(int numT, int numS){
	numStructures = numS;
	numTissues = numT;
	tissueProb = new irtkRealImage*[numTissues];
	structureProb = new irtkRealImage*[numStructures];
	tissueGaussian = new irtkGaussian*[numTissues];
	structureGaussian = new irtkGaussian*[numStructures];
}

void MOG::setTissueProb(int num, irtkRealImage *prob, double my, double sigma){
	tissueProb[num] = prob;
	tissueGaussian[num] = new irtkGaussian;
	tissueGaussian[num]->Initialise(my, sigma);
}

void MOG::setStructureProb(int num, irtkRealImage *prob, double my, double sigma){
	structureProb[num] = prob;
	structureGaussian[num] = new irtkGaussian;
	structureGaussian[num]->Initialise(my, sigma);
}

double MOG::evaluate(double val, int x, int y, int z){
	return evaluate(val, x, y, z, 0);
}

double MOG::evaluate(double val, int x, int y, int z, int t){

	double prob = 0;
	double pS = 0, pT = 0;
	double gammaStruct = 0;
	for(int i = 0; i < numStructures; i++){
		double gamma = structureProb[i]->GetAsDouble(x, y, z, t)/255;
		pS += gamma * structureGaussian[i]->Evaluate(val);
		gammaStruct += gamma;
	}
	if(gammaStruct != 1){
		for(int i = 0; i < numTissues; i++){
			pT += tissueProb[i]->GetAsDouble(x, y, z) *  tissueGaussian[i]->Evaluate(val);
		}
		if(pT == 0){
			for(int i = 0; i < numTissues; i++){
						pT += tissueGaussian[i]->Evaluate(val)/3;
					}

			}
	}


	prob = (1-gammaStruct) * pT + (gammaStruct) * pS;
	return prob;
}

