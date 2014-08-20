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

#include <irtkSpectralClustering.h>
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

irtkSpectralClustering::irtkSpectralClustering(int nrSubjects, int nrFeatures){	
	_nrSubjects = nrSubjects;
	_nrFeatures = nrFeatures;
	_features = new double*[nrSubjects];
	for(int i = 0; i < _nrSubjects; i++){
		_features[i] = new double[_nrFeatures];
	}
	

}



void irtkSpectralClustering::DoSpectralEmbedding(){
		
	int rows = _nrSubjects;
	irtkMatrix wTrans;
	wTrans = _w;
	wTrans.Transpose();	
	_w = _w + wTrans;
	irtkMatrix dd, d, lsym;
	lsym.Initialize(rows, rows);
	for(int i = 0; i < rows; i++){
	 	for(int j = 0; j < rows; j++){
	 		lsym(i,j) = -_w(i,j);
	 	}
	}	
	for(int i = 0; i < rows; i++){
		double dVal = 0;
	 	for(int j = 0; j < rows; j++){	 			
	 			dVal += _w(i,j);
	 	}
	 	double ddVal = 1/sqrt(dVal);
	 	lsym(i,i) += dVal;
	 	for(int k = 0; k < rows; k++){
	 		lsym(i,k) *= ddVal;
	 		lsym(k,i) *= ddVal;	 		
	 	}	 	
	 }
	 irtkEigenAnalysis ea(rows);
	 for(int i = 0; i < rows; i++){
	 	for(int j = 0; j < rows; j++){
	 		ea.Matrix(i,j) = lsym(i,j);
	 	}
	 }
	 ea.IncrSortEigenStuff();
//	 cout << endl;
	 for(int i = 0; i < rows; i++){
	 	for(int j = 0; j < _nrFeatures; j++){
	 		_features[i][j] = ea.Eigenvector( i, j+1);
//	 		cout << ea.Eigenvector( i, j+1) << ",";
	 	}
//	 	cout << endl;
	 }
	 
}
