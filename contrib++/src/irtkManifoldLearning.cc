/*=========================================================================

  Library   : Image Registration Toolkit (IRTK)
  Module    : $Id$
  Copyright : Imperial College, Department of Computing
              Visual Information Processing (VIP), 2008 onwards
  Date      : $Date$
  Version   : $Revision$
  Changes   : $Author$

Copyright (c) 1999-2014 and onwards, Imperial College London
All rights reserved.
See LICENSE for details

=========================================================================*/

#include <irtkManifoldLearning.h>
#include <irtkEigenAnalysis.h>
#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <math.h>
#include <sstream>
#include <algorithm>
#include <sys/types.h>
#include <errno.h>
#include <sstream>

 
irtkManifoldLearning::irtkManifoldLearning(){
	_distances = NULL;
}

double irtkManifoldLearning::GetDistance(int i, int j){
	if(_distances == NULL){
		EstablishDistances();
	}
	return _distances[i][j];
}

void irtkManifoldLearning::EstablishDistances(){
	if(_distances == NULL){
		_distances = new double*[_nrSubjects];
		for(int i = 0; i < _nrSubjects; i++){
			_distances[i] = new double[_nrSubjects];
		}
	}
	for(int i = 0; i < _nrSubjects; i++){
		_distances[i][i] = 0;
		for(int j = i+1; j < _nrSubjects; j++){
			double dist = 0;
			for(int f = 0; f < _nrFeatures; f++){
				dist += (_features[i][f]-_features[j][f])*(_features[i][f]-_features[j][f]);
			}
			dist = sqrt(dist);
			_distances[i][j] = _distances[j][i] = dist;
		}
	}
}

void irtkManifoldLearning::GetNeighbours(int node, int nrNeighbours, int * neighbours){
	double * distTemp = new double[_nrSubjects];

	if(_distances == NULL){
		EstablishDistances();
	}
	for(int i = 0; i < _nrSubjects; i++){
		distTemp[i] = _distances[node][i];
//		cout << _distances[node][i];
	}

	int max = 65000;
//	neighbours = new int[nrNeighbours];
	for(int k = 0; k < nrNeighbours; k++){
		double min = max;
		int nextNeighbour = -1;
		for(int i = 0; i < _nrSubjects; i++){
			if(i != node && distTemp[i] < min){
				min = distTemp[i];
				nextNeighbour = i;				
			}		
		}

		distTemp[nextNeighbour] = max;
		cout << nextNeighbour << endl;
		neighbours[k] = nextNeighbour;		
	}
	delete distTemp;	
}

double ** irtkManifoldLearning::GetEmbedding(){
	
	

	return _features;
}

void irtkManifoldLearning::Initialize(double ** input){
	_w.Initialize(_nrSubjects, _nrSubjects);
	for(int i = 0; i < _nrSubjects; i++){
		for(int j = 0; j < _nrSubjects; j++){
			_w(i,j) = input[i][j];		
		}
	
	}
	

}

void irtkManifoldLearning::Initialize(irtkMatrix input){

	_w.Initialize(_nrSubjects, _nrSubjects);
	for(int i = 0; i < _nrSubjects; i++){
		for(int j = 0; j < _nrSubjects; j++){
			_w(i,j) = input.Get(i, j);
		}

	}


}

void irtkManifoldLearning::Initialize(string csvFilename){
	
	
	ifstream input;
	input.open(csvFilename.c_str());
	string value;
	int rows = 0;
	while(input.good()){
		getline(input, value, '\n');
		rows++;	
	}
	rows--;
	if(rows != _nrSubjects){
		cerr << "Dimensions in .csv file do not agree with nrSubjects" << endl;
	}
	_w.Initialize(rows, rows);
	ifstream input2;
	input2.open(csvFilename.c_str());
	int rowctr = 0;
	int colctr = 0;
	while(input2.good()){
		getline(input2, value, '\n');
		int pos = value.find(",");
		int posOld = -1;
		while(colctr < rows && rowctr < rows){		
			string stringdata = value.substr(posOld+1,pos-posOld-1);
			posOld = pos;			
			pos = value.find(",", posOld+1);
			double data = atof(stringdata.c_str());
		 	_w(rowctr, colctr) = data;
		 	colctr++;			
		}
		colctr = 0;
		rowctr ++;
	}
	
}

void irtkManifoldLearning::WriteFeatures(string filename){
	ofstream output;
	output.open(filename.c_str());
	for(int i = 0; i < _nrSubjects; i++){
		for(int j = 0; j < _nrFeatures; j++){
			output << _features[i][j];
			if(j == _nrFeatures-1){
				output << endl;
			}
			else{
				output << ",";
			}
		}
	}

}

 
