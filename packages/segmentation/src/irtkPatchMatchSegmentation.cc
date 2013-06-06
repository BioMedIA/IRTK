/*=========================================================================

Library   : Image Registration Toolkit (IRTK)
Module    : $Id$
Copyright : Imperial College, Department of Computing
Visual Information Processing (VIP), 2008 onwards
Date      : $Date$
Version   : $Revision$
Changes   : $Author$

=========================================================================*/

#include <irtkSegmentationFunction.h>

irtkPatchMatchSegmentation::irtkPatchMatchSegmentation
	(irtkGreyImage *target, irtkGreyImage **source, 
	irtkRealImage *targetdistance, irtkRealImage **sourcedistance,
	int radius, int nimages, int nneighbour)
	:irtkPatchMatch(target,source,radius,nimages,nneighbour,1)
{
	//construct distance
	targetlabeldistance = targetdistance;
	sourcelabeldistance = sourcedistance;
	if(targetdistance != NULL)
		sweight = 1.0/targetdistance->GetT();
	else
		sweight = 0;

	distancenorm = NULL;

	//normalize the distances with respect to image
	if(targetlabeldistance != NULL){
		this->normalizedistances();
	}

	this->initialguess();
}

irtkPatchMatchSegmentation::~irtkPatchMatchSegmentation()
{
	for(int i = 0; i < target->GetNumberOfVoxels(); i++){
		delete []nnfs[i];
	}
	delete []nnfs;

	delete targetgradient;
	for(int n = 0; n < nimages; n++){
		delete sourcesgradient[n];
		delete search[n];
	}
	delete []sourcesgradient;
	delete []search;

	delete []votevalues;
	delete []voteweights;
	if(decimated != NULL){
		delete decimated;
		delete decimatedgradient;
		delete blured;
		for(int n = 0; n < nimages; n++){
			delete bluredgradient[n];
		}
		delete []bluredgradient;
	}

	if(distancenorm != NULL){
		delete []distancenorm;
	}
}

void irtkPatchMatchSegmentation::normalizedistances(){
	distancenorm = new double[targetlabeldistance->GetT()];
	double min,max;
	for(int t = 0; t < targetlabeldistance->GetT(); t++){
		targetlabeldistance->GetMinMax(&max,&min);
		for(int k = 0; k < targetlabeldistance->GetZ(); k++){
			for(int j = 0; j < targetlabeldistance->GetY(); j++){
				for(int i = 0; i < targetlabeldistance->GetX(); i++){
					if(targetlabeldistance->GetAsDouble(i,j,k,t) < min){
						min = targetlabeldistance->GetAsDouble(i,j,k,t);
					}
					if(targetlabeldistance->GetAsDouble(i,j,k,t) > max){
						max = targetlabeldistance->GetAsDouble(i,j,k,t);
					}
				}
			}
		}
		min = max - min;
		distancenorm[t] = maxdistance/min;

		cout << "distance norm for distance map " << t << " is " << distancenorm[t] << endl;

	}
}

void irtkPatchMatchSegmentation::initialguess(){
	//populate NNF field with initial guess
	cout << "populating NNF field with initial total random guess" << endl;
	// initialize random seed
	srand ( time(NULL) );
	int index = 0;
	double x,y,z;
	for(int k = 0; k < target->GetZ(); k++){
		for(int j = 0; j < target->GetY(); j++){
			for(int i = 0; i < target->GetX(); i++){
				for(int n = 0; n < nneighbour; n++){
					//TODO need to change this if source dimention != target dimention
					nnfs[index][n].n = rand()%nimages;
					x = rand()%target->GetX();
					y = rand()%target->GetY();
					z = rand()%target->GetZ();
					//x = i;
					//y = j;
					//z = k;
					target->ImageToWorld(x,y,z);
					sources[nnfs[index][n].n]->WorldToImage(x,y,z);

					if(x < 0) x = 0;
					if(y < 0) y = 0;
					if(z < 0) z = 0;
					if(x > sources[nnfs[index][n].n]->GetX() - 1)
						x = sources[nnfs[index][n].n]->GetX() - 1;
					if(y > sources[nnfs[index][n].n]->GetY() - 1)
						y = sources[nnfs[index][n].n]->GetY() - 1;
					if(z > sources[nnfs[index][n].n]->GetZ() - 1)
						z = sources[nnfs[index][n].n]->GetZ() - 1;

					nnfs[index][n].x = round(x);
					nnfs[index][n].y = round(y);
					nnfs[index][n].z = round(z);
					nnfs[index][n].weight = maxdistance;
				}
				index++;
			}
		}

	}
	cout << "done" << endl;
}

void irtkPatchMatchSegmentation::setWeight(double weight){
	if(targetlabeldistance != NULL)
		sweight = weight/targetlabeldistance->GetT();
}

/// calculate distance between patches
double irtkPatchMatchSegmentation::distance(int x1, int y1, int z1, int x2, int y2, int z2, int n){
	int i,j,k,i1,j1,k1,i2,j2,k2,t,g,tmpradius,increase;
	double dif = 0, tmp, count;
	short value1, value2, values;
	count = 0;
	increase = 1;
	tmpradius = this->radius;

	if(target->Get(x1,y1,z1) < 0){
		return maxdistance - 1;
	}

	for(k = - tmpradius; k <= tmpradius; k+=increase){
		k1 = k + z1;
		k2 = k + z2;
		for(j = -tmpradius; j <= tmpradius; j+=increase){
			j1 = j + y1;
			j2 = j + y2;
			for(i = -tmpradius; i <= tmpradius; i+=increase){

				i1 = i + x1;
				i2 = i + x2;

				if(i1 < 0 || (i1 > target->GetX() - 1) 
					|| j1 < 0 || (j1 > target->GetY() - 1) 
					|| k1 < 0 || (k1 > target->GetZ() - 1) ){
						//do not count
				}else if(i2 < 0 || (i2 > sources[n]->GetX() - 1)
					|| j2 < 0 || (j2 > sources[n]->GetY() - 1)
					|| k2 < 0 || (k2 > sources[n]->GetZ() - 1)){
						if(isdecimated == false){
							dif += maxdistance*4;
							count += 4;
						}
						if(decimated != NULL){
							dif += maxdistance*3;
							count += 3;
						}
						if(targetlabeldistance != NULL){
							dif += maxdistance*sweight*targetlabeldistance->GetT();
							count += sweight*targetlabeldistance->GetT();
						}
				}else{
					//distance between resolved image and atlases
					if(isdecimated == false){
						//reconstructed image not just decimated image
						value1 = target->Get(i1,j1,k1);
						if(value1 >= 0){		
							values = sources[n]->Get(i2,j2,k2);
							tmp = double(value1 - values);
							dif += sqrt(tmp*tmp);
							count++;
							//distance between gradient
							for(g = 0; g < 3; g++){
								value1 = targetgradient->Get(i1,j1,k1,g);
								value2 = sourcesgradient[n]->Get(i2,j2,k2,g);
								tmp = double(value1 - value2);
								dif += sqrt(tmp*tmp);
								count++;
							}
						}
					}
					//distance between decimated and atlases
					if(decimated != NULL){
						//check if decimated image exists
						value1 = decimated->Get(i1,j1,k1);
						if(value1 >= 0){
							values = blured[n]->Get(i2,j2,k2);
							tmp = double(value1 - values);
							dif += sqrt(tmp*tmp);
							count++;
							//distance between gradient
							for(g = 0; g < 2; g++){
								value1 = decimatedgradient->Get(i1,j1,k1,g);
								value2 = bluredgradient[n]->Get(i2,j2,k2,g);
								tmp = double(value1 - value2);
								dif += sqrt(tmp*tmp);
								count++;
							}
						}
					}
					//distance between label distances
					if(targetlabeldistance != NULL && sweight > 0){
						//check if decimated image exists
						for(t = 0; t < targetlabeldistance->GetT(); t++){
							value1 = targetlabeldistance->Get(i1,j1,k1,t);
							values = sourcelabeldistance[n]->Get(i2,j2,k2,t);
							tmp = double(value1 - values);
							dif += sweight*sqrt(tmp*tmp)*distancenorm[t];
							count += sweight;
						}
					}
				}	
			}
		}
	}

	if(count < 1)
		return maxdistance - 1;
	else
		return dif/count;
}

/// expectation maximization
void irtkPatchMatchSegmentation::EMstep(){
	/// Expectation create the vote matrix
	int index = 0, x, y, z, n, o;
	double weight;

	for(int i = 0; i < target->GetNumberOfVoxels(); i++){
		votevalues[i] = 0;
		voteweights[i] = 0;
	}

	for(int k = 0; k < target->GetZ(); k++){
		for(int j = 0; j < target->GetY(); j++){
			for(int i = 0; i < target->GetX(); i++){
				for(int o = 0; o < nneighbour; o++){
					x = nnfs[index][o].x;
					y = nnfs[index][o].y;
					z = nnfs[index][o].z;
					n = nnfs[index][o].n;
					weight = nnfs[index][o].weight;

					weight = exp(-(weight/maxdistance*32));

					this->votepatch(i,j,k,x,y,z,n,weight);
				}
				index++;
			}
		}
	}
	/// Maximization non local mean fusion
	index = 0;
	for(int k = 0; k < target->GetZ(); k++){
		for(int j = 0; j < target->GetY(); j++){
			for(int i = 0; i < target->GetX(); i++){
				if(voteweights[index] > 0 && target->Get(i,j,k) >= 0){
					target->PutAsDouble(i,j,k,votevalues[index]/voteweights[index]);
				}else{
					target->PutAsDouble(i,j,k,-1);
				}
				index++;
			}
		}
	}

	isdecimated = false;

	if(debug == true){
		char buffer[255];
		sprintf(buffer, "tmptarget%d.nii.gz", rand()%1000);
		target->Write(buffer);
	}
}