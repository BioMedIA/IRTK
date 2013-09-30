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

irtkMAPatchMatchSuperResolution::irtkMAPatchMatchSuperResolution
	(irtkGreyImage *target, irtkGreyImage **source, int radius, int nimages, int nneighbour)
	:irtkMAPatchMatch(target, source, radius, nimages, nneighbour){

		this->decimated = NULL;
		this->isdecimated = false;

		this->votevalues = new double[target->GetNumberOfVoxels()];
		this->voteweights = new double[target->GetNumberOfVoxels()];

		this->initialguess();
}

irtkMAPatchMatchSuperResolution::~irtkMAPatchMatchSuperResolution(){
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
}

void irtkMAPatchMatchSuperResolution::runEMIteration(int maxiterations){

	double localweight;
	double previousweight;
	int iteration = 0;

	localweight = this->initialize();

	cout << "minimizing NNF..." << endl;

	//debug tmpflow
	if(debug == true){
		irtkImageAttributes atr = target->GetImageAttributes();
		atr._t = 5;
		irtkRealImage *flowtmp = new irtkRealImage(atr);

		int index = 0;
		for(int k = 0; k < target->GetZ(); k++){
			for(int j = 0; j < target->GetY(); j++){
				for(int i = 0; i < target->GetX(); i++){
					flowtmp->Put(i,j,k,0,this->nnfs[index][0].x);
					flowtmp->Put(i,j,k,1,this->nnfs[index][0].y);
					flowtmp->Put(i,j,k,2,this->nnfs[index][0].z);
					flowtmp->Put(i,j,k,3,this->nnfs[index][0].n);
					flowtmp->Put(i,j,k,4,this->nnfs[index][0].weight);
					index++;
				}
			}
		}

		char buffer[255];
		sprintf(buffer, "flow_initial%d.nii.gz", rand()%1000);
		flowtmp->Write(buffer);
		delete flowtmp;
	}

	while((iteration == 0 || (previousweight-localweight)>0.0001 ) && iteration < maxiterations){
		cout << "iteration: " << iteration << " ";
		previousweight = localweight;
		if(debug == true && iteration == 0 && isdecimated == true){
			localweight = this->minimizeflowwithdebug();
		}else{
			localweight = this->minimizeflow();
		}
		cout << "total distance " << localweight << endl;
		iteration++;
	}
	cout << "done" << endl;

	cout << "EM step... ";
	this->EMstep();
	cout << "done" << endl;
}

/// Generate HD image
void irtkMAPatchMatchSuperResolution::generateImage(){
	this->initialize();
	this->EMstep();
}

/// initialize field's weight
double irtkMAPatchMatchSuperResolution::initialize(){
	cout << "update gradients" << endl;
	// Compute spatial gradient of source image
	irtkGradientImageFilter<short> gradient(irtkGradientImageFilter<short>::GRADIENT_VECTOR);
	irtkGreyImage tmp = *target;
	gradient.SetInput (&tmp);
	gradient.SetOutput(targetgradient);
	gradient.SetPadding(-1);
	gradient.Run();

	if(debug == true){
		target->Write("targetincode.nii.gz");
		targetgradient->Write("targetgradient.nii.gz");
	}

	if(decimated != NULL){
		tmp = *decimated;
		gradient.SetInput (&tmp);
		gradient.SetOutput(decimatedgradient);
		gradient.SetPadding(-1);
		gradient.Run();

		if(debug == true){
			decimatedgradient->Write("decimatedgradient.nii.gz");
		}
	}

	cout << "initialize NNF's weight...";
	cout.flush();
	double totalweight = 0;
	/// recalculate of the fields
	int index = 0, iteration;
	for(int k = 0; k < target->GetZ(); k++){
		for(int j = 0; j < target->GetY(); j++){
			for(int i = 0; i < target->GetX(); i++){
				for(int n = 0; n < nneighbour; n++){
					double difference = this->distance(i,j,k,
						nnfs[index][n].x,nnfs[index][n].y,
						nnfs[index][n].z,nnfs[index][n].n);

					nnfs[index][n].weight = difference;

					// if weight == maxdistance, the link is not good try to find a better link
					iteration = 0;
					while(nnfs[index][n].weight >= maxdistance && iteration < 10){
						randomlink(i,j,k, n, index);
						iteration++;
					}
					totalweight += nnfs[index][n].weight;
				}
				index++;
			}
		}
	}
	cout << "done" << endl;

	return totalweight;
}

/// vote patch
void irtkMAPatchMatchSuperResolution::votepatch(int xt, int yt, int zt, int xs, int ys, int zs, int n, double weight){
	int index1;
	for(int k = - radius; k <= radius; k++){
		for(int j = - radius; j <= radius; j++){
			for(int i = - radius; i <= radius; i++){

				if(xt + i < 0 || xt + i > target->GetX() - 1
					|| xs + i < 0 || xs + i > sources[n]->GetX() - 1
					|| yt + j < 0 || yt + j > target->GetY() - 1
					|| ys + j < 0 || ys + j > sources[n]->GetY() - 1
					|| zt + k < 0 || zt + k > target->GetZ() - 1
					|| zs + k < 0 || zs + k > sources[n]->GetZ() - 1){
						//do nothing

				}else{
					index1 = target->VoxelToIndex(xt+i,yt+j,zt+k);

					if(sources[n]->GetAsDouble(xs+i,ys+j,zs+k) >= 0){

						votevalues[index1] += weight*sources[n]->GetAsDouble(xs+i,ys+j,zs+k);
						voteweights[index1] += weight;
					}
				}
			}
		}
	}
}


/// expectation maximization
void irtkMAPatchMatchSuperResolution::EMstep(){
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
				if(voteweights[index] > 0){
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

/// find minum flow
double irtkMAPatchMatchSuperResolution::minimizeflow(){
	int i,j,k, x,y,z,fcount, bcount, rcount, count, index,index1,index2,index3;
	double sumweight = 0;
	/// Forward propergation
	fcount = 0;
	index = 0;
	//x-1
	index1 = -1;
	//y-1
	index2 = -target->GetX();
	//z-1
	index3 = -target->GetX()*target->GetY();
	for(int k = 0; k < target->GetZ(); k++){
		for(int j = 0; j < target->GetY(); j++){
			for(int i = 0; i < target->GetX(); i++){
				count = 0;
				//propergate from x - 1
				if(i > 0){
					x = i -1;
					y = j;
					z = k;
					count += this->propergate(x,y,z,i,j,k,1,0,0,index1,index);
				}
				//propergate from y - 1
				if(j > 0){
					x = i;
					y = j - 1;
					z = k;
					count += this->propergate(x,y,z,i,j,k,0,1,0,index2,index);
				}
				//propergate from z - 1
				if(k > 0){
					x = i;
					y = j;
					z = k - 1;
					count += this->propergate(x,y,z,i,j,k,0,0,1,index3,index);
				}

				if(count > 0)
					fcount ++;

				index++;
				index1++;
				index2++;
				index3++;
			}
		}
	}

	/// Random search
	rcount = 0;
	index = 0;
	for(int k = 0; k < target->GetZ(); k++){
		for(int j = 0; j < target->GetY(); j++){
			for(int i = 0; i < target->GetX(); i++){
				for(int o = 0; o < nneighbour; o++){
					if(this->randomlink(i,j,k,o,index) > 0){
						rcount++;
					}
				}
				index++;
			}
		}
	}

	/// Backward propergation
	bcount = 0;
	index = target->GetNumberOfVoxels() - 1;
	//x-1
	index1 = index+1;
	//y-1
	index2 = index+target->GetX();
	//z-1
	index3 = index+target->GetX()*target->GetY();
	for(int k = target->GetZ() - 1; k >= 0; k--){
		for(int j = target->GetY() - 1; j >= 0; j--){
			for(int i = target->GetX() - 1; i >=0 ; i--){
				count = 0;
				//propergate from x + 1
				if(i < target->GetX() - 1){
					x = i + 1;
					y = j;
					z = k;
					count += this->propergate(x,y,z,i,j,k,-1,0,0,index1,index);
				}
				//propergate from y + 1
				if(j < target->GetY() - 1){
					x = i;
					y = j + 1;
					z = k;
					count += this->propergate(x,y,z,i,j,k,0,-1,0,index2,index);
				}
				//propergate from z + 1
				if(k < target->GetZ() - 1){
					x = i;
					y = j;
					z = k + 1;
					count += this->propergate(x,y,z,i,j,k,0,0,-1,index3,index);
				}

				if(count > 0)
					bcount ++;

				index--;
				index1--;
				index2--;
				index3--;
			}
		}
	}

	/// Random search
	index = 0;
	for(int k = 0; k < target->GetZ(); k++){
		for(int j = 0; j < target->GetY(); j++){
			for(int i = 0; i < target->GetX(); i++){
				for(int o = 0; o < nneighbour; o++){
					if(this->randomlink(i,j,k,o,index) > 0){
						rcount++;
					}
					if(this->nnfs[index][o].weight < maxdistance - 1){
						sumweight += this->nnfs[index][o].weight;
					}
				}
				index++;
			}
		}
	}

	cout << "number of fields changed: " << fcount << " " << bcount << " " << rcount << " ";

	if(debug == true){
		/// Debug write flow to a image using the dynamicid

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

		irtkImageAttributes atr = target->GetImageAttributes();
		atr._t = 6;
		irtkRealImage *flowtmp = new irtkRealImage(atr);

		index = 0;
		for(int k = 0; k < target->GetZ(); k++){
			for(int j = 0; j < target->GetY(); j++){
				for(int i = 0; i < target->GetX(); i++){
					flowtmp->Put(i,j,k,0,this->nnfs[index][0].x);
					flowtmp->Put(i,j,k,1,this->nnfs[index][0].y);
					flowtmp->Put(i,j,k,2,this->nnfs[index][0].z);
					flowtmp->Put(i,j,k,3,this->nnfs[index][0].n);
					flowtmp->Put(i,j,k,4,this->nnfs[index][0].weight);
					if(voteweights[index] > 0){
						flowtmp->PutAsDouble(i,j,k,5,votevalues[index]/voteweights[index]);
					}else{
						flowtmp->PutAsDouble(i,j,k,5,-1);
					}
					index++;
				}
			}
		}

		char buffer[255];
		sprintf(buffer, "flow%d_%d_%d.nii.gz", fcount, bcount, rcount);
		flowtmp->Write(buffer);
		delete flowtmp;
	}

	return log(sumweight);
}

/// vote weight matrix
void irtkMAPatchMatchSuperResolution::voteweight(int zd, int mode){
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

	/// Debug write flow to a image using the dynamicid

	irtkImageAttributes atr = target->GetImageAttributes();
	atr._t = 6;
	irtkRealImage *flowtmp = new irtkRealImage(atr);

	index = 0;
	for(int k = 0; k < target->GetZ(); k++){
		for(int j = 0; j < target->GetY(); j++){
			for(int i = 0; i < target->GetX(); i++){
				flowtmp->Put(i,j,k,0,this->nnfs[index][0].x);
				flowtmp->Put(i,j,k,1,this->nnfs[index][0].y);
				flowtmp->Put(i,j,k,2,this->nnfs[index][0].z);
				flowtmp->Put(i,j,k,3,this->nnfs[index][0].n);	
				if(voteweights[index] > 0){
					if(zd >= k || mode > 0){
						flowtmp->Put(i,j,k,4,this->nnfs[index][0].weight);
						flowtmp->PutAsDouble(i,j,k,5,votevalues[index]/voteweights[index]);
					}else{
						double tmpvalue = (this->nnfs[index][0].n
							+this->nnfs[index][0].weight)/2;
						if(tmpvalue > this->maxdistance - 1)
							tmpvalue = this->maxdistance - 1;
						flowtmp->Put(i,j,k,4,tmpvalue);
						tmpvalue = votevalues[index]/voteweights[index]/3
							+this->nnfs[index][0].n+this->nnfs[index][0].weight;
						if(tmpvalue > this->maxdistance - 1)
							tmpvalue = this->maxdistance - 1;
						flowtmp->PutAsDouble(i,j,k,5,tmpvalue);
					}
				}else{
					flowtmp->Put(i,j,k,4,this->nnfs[index][0].weight);
					flowtmp->PutAsDouble(i,j,k,5,-1);
				}
				index++;
			}
		}
	}

	char buffer[255];
	sprintf(buffer, "flowdebug%d_%d.nii.gz", mode, zd);
	flowtmp->Write(buffer);
	delete flowtmp;
}

/// find minum flow
double irtkMAPatchMatchSuperResolution::minimizeflowwithdebug(){
	int i,j,k, x,y,z,fcount, bcount, rcount, count, index,index1,index2,index3;
	double sumweight = 0;
	/// Forward propergation
	fcount = 0;
	index = 0;
	//x-1
	index1 = -1;
	//y-1
	index2 = -target->GetX();
	//z-1
	index3 = -target->GetX()*target->GetY();
	for(int k = 0; k < target->GetZ(); k++){
		for(int j = 0; j < target->GetY(); j++){
			for(int i = 0; i < target->GetX(); i++){
				count = 0;
				//propergate from x - 1
				if(i > 0){
					x = i -1;
					y = j;
					z = k;
					count += this->propergate(x,y,z,i,j,k,1,0,0,index1,index);
				}
				//propergate from y - 1
				if(j > 0){
					x = i;
					y = j - 1;
					z = k;
					count += this->propergate(x,y,z,i,j,k,0,1,0,index2,index);
				}
				//propergate from z - 1
				if(k > 0){
					x = i;
					y = j;
					z = k - 1;
					count += this->propergate(x,y,z,i,j,k,0,0,1,index3,index);
				}

				if(count > 0)
					fcount ++;

				index++;
				index1++;
				index2++;
				index3++;
			}
		}

		if(k%4 == 0)
			this->voteweight(k,0);

	}

	/// Random search
	rcount = 0;
	index = 0;
	for(int k = 0; k < target->GetZ(); k++){
		for(int j = 0; j < target->GetY(); j++){
			for(int i = 0; i < target->GetX(); i++){
				for(int o = 0; o < nneighbour; o++){
					if(this->randomlink(i,j,k,o,index) > 0){
						rcount++;
					}
				}
				index++;
			}
		}
		if(k%4 == 0)
			this->voteweight(k,1);
	}

	/// Backward propergation
	bcount = 0;
	index = target->GetNumberOfVoxels() - 1;
	//x-1
	index1 = index+1;
	//y-1
	index2 = index+target->GetX();
	//z-1
	index3 = index+target->GetX()*target->GetY();
	for(int k = target->GetZ() - 1; k >= 0; k--){
		for(int j = target->GetY() - 1; j >= 0; j--){
			for(int i = target->GetX() - 1; i >=0 ; i--){
				count = 0;
				//propergate from x + 1
				if(i < target->GetX() - 1){
					x = i + 1;
					y = j;
					z = k;
					count += this->propergate(x,y,z,i,j,k,-1,0,0,index1,index);
				}
				//propergate from y + 1
				if(j < target->GetY() - 1){
					x = i;
					y = j + 1;
					z = k;
					count += this->propergate(x,y,z,i,j,k,0,-1,0,index2,index);
				}
				//propergate from z + 1
				if(k < target->GetZ() - 1){
					x = i;
					y = j;
					z = k + 1;
					count += this->propergate(x,y,z,i,j,k,0,0,-1,index3,index);
				}

				if(count > 0)
					bcount ++;

				index--;
				index1--;
				index2--;
				index3--;
			}
		}
		if(k%4 == 0)
			this->voteweight(target->GetZ() - 1 - k,2);
	}

	/// Random search
	index = 0;
	for(int k = 0; k < target->GetZ(); k++){
		for(int j = 0; j < target->GetY(); j++){
			for(int i = 0; i < target->GetX(); i++){
				for(int o = 0; o < nneighbour; o++){
					if(this->randomlink(i,j,k,o,index) > 0){
						rcount++;
					}
					if(this->nnfs[index][o].weight < maxdistance - 1){
						sumweight += this->nnfs[index][o].weight;
					}
				}
				index++;
			}
		}
		if(k%4 == 0)
			this->voteweight(k,3);
	}

	cout << "number of fields changed: " << fcount << " " << bcount << " " << rcount << " ";

	return log(sumweight);
}

double irtkMAPatchMatchSuperResolution::distance(int x1, int y1, int z1, int x2, int y2, int z2, int n){
	int i,j,k,i1,j1,k1,i2,j2,k2,t,count,g,tmpradius,increase;
	double dif = 0, tmp;
	short value1, value2, values;
	count = 0;
	increase = 1;
	tmpradius = this->radius;

	//check is search
	if(x2 < 0 || x2 > sources[n]->GetX() - 1
		|| y2 < 0 || y2 > sources[n]->GetY() - 1
		|| z2 < 0 || z2 > sources[n]->GetZ() - 1
		|| search[n]->Get(x2,y2,z2) < 1)
		return maxdistance;

	if(isdecimated == false && target->Get(x1,y1,z1) < 0){
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

				if(i1 < 0|| i2 < 0 ||
					(i1 > target->GetX() - 1) 
					|| (i2 > sources[n]->GetX() - 1)
					|| j1 < 0|| j2 < 0
					|| (j1 > target->GetY() - 1) 
					|| (j2 > sources[n]->GetY() - 1)
					|| k1 < 0|| k2 < 0
					|| (k1 > target->GetZ() - 1) 
					|| (k2 > sources[n]->GetZ() - 1)){
						if(isdecimated == false){
							dif += maxdistance*4;
							count += 4;
						}
						if(decimated != NULL){
							dif += maxdistance*3;
							count += 3;
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
				}	
			}
		}
	}

	if(count < 1)
		return maxdistance - 1;
	else
		return dif/count;
}

void irtkMAPatchMatchSuperResolution::setDecimatedImage(irtkGreyImage *input)
{
	double x, y, z, offsetx, offsety, offsetz;
	decimated = new irtkGreyImage(target->GetImageAttributes());
	decimatedgradient = new irtkGreyImage();

	for(int k = 0; k < decimated->GetZ(); k++){
		for(int j = 0; j < decimated->GetY(); j++){
			for(int i = 0; i < decimated->GetX(); i++){
				decimated->Put(i,j,k,-1);
				target->Put(i,j,k,-1);
			}
		}
	}

	for(int k = 0; k < input->GetZ(); k++){
		for(int j = 0; j < input->GetY(); j++){
			for(int i = 0; i < input->GetX(); i++){
				x = i;
				y = j;
				z = k;

				input->ImageToWorld(x,y,z);
				decimated->WorldToImage(x,y,z);

				//find nneighbor
				decimated->Put(round(x),round(y),round(z),input->Get(i,j,k));
				target->Put(round(x),round(y),round(z),input->Get(i,j,k));
			}
		}
	}

	blured = new irtkGreyImage*[nimages];

	for(int n = 0; n < nimages; n++){
		blured[n] = new irtkGreyImage(*sources[n]);
		//FWHM 2.35 sigma
		irtkGaussianBlurringWithPadding<irtkGreyPixel> gaussianBlurring(input->GetZSize()/2.94,0);
		gaussianBlurring.SetInput (blured[n]);
		gaussianBlurring.SetOutput(blured[n]);
		gaussianBlurring.RunZ();
	}

	//create blured image
	cout << "calculating blured gradient image" << endl;

	irtkGradientImageFilter<short> gradient(irtkGradientImageFilter<short>::GRADIENT_VECTOR);

	bluredgradient = new irtkGreyImage*[nimages];
	for(int n = 0; n < nimages; n++){
		bluredgradient[n] = new irtkGreyImage();
		irtkGreyImage tmp = *blured[n];
		gradient.SetInput (&tmp);
		gradient.SetOutput(bluredgradient[n]);
		gradient.SetPadding(-1);
		gradient.Run();
	}
	cout << "done" << endl;

	isdecimated = true;

	if(debug == true){
		decimated->Write("decimatedimage.nii.gz");
	}
}