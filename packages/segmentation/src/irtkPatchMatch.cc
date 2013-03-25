/*=========================================================================

Library   : Image Registration Toolkit (IRTK)
Module    : $Id$
Copyright : Imperial College, Department of Computing
Visual Information Processing (VIP), 2008 onwards
Date      : $Date$
Version   : $Revision$
Changes   : $Author$

=========================================================================*/

#include<irtkPatchMatch.h>

irtkPatchMatch::irtkPatchMatch(irtkGreyImage *target, irtkGreyImage **source, int radius, int nimages, int nneighbour, int slevels){

	this->target = target;
	this->sources = source;
	this->radius = radius;
	this->nimages = nimages;
	this->nneighbour = nneighbour;
	this->randomrate = 0.5;
	this->decimated = NULL;
	this->isdecimated = false;
	this->slevels = slevels;
	this->debug = false;
	target->GetMinMax(&minlabel, &maxlabel);
	maxdistance = maxlabel - minlabel;

	//create NNF field
	this->nnfs = new NearstNeighbor*[target->GetNumberOfVoxels()];
	for(int i = 0; i < target->GetNumberOfVoxels(); i++){
		nnfs[i] = new NearstNeighbor[this->nneighbour];
	}

	this->votevalues = new double[target->GetNumberOfVoxels()];
	this->voteweights = new double[target->GetNumberOfVoxels()];

	//populate NNF field with initial guess
	cout << "populating NNF field with initial guess" << endl;
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
					x = i;
					y = j;
					z = k;
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

	//calculate gradient image
	cout << "calculating atlas gradient image" << endl;

	irtkGradientImageFilter<short> gradient(irtkGradientImageFilter<short>::GRADIENT_VECTOR);

	targetgradient = new irtkGreyImage();

	sourcesgradient = new irtkGreyImage*[nimages];
	search = new irtkGreyImage*[nimages];
	for(int n = 0; n < nimages; n++){
		sourcesgradient[n] = new irtkGreyImage();
		irtkGreyImage tmp = *sources[n];
		gradient.SetInput (&tmp);
		gradient.SetOutput(sourcesgradient[n]);
		gradient.SetPadding(-1);
		gradient.Run();

		search[n] = new irtkGreyImage(*source[n]);

	}

	this->createsearchimages();

	cout << "done" << endl;
}

irtkPatchMatch::~irtkPatchMatch(){
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

void irtkPatchMatch::runEMIteration(int maxiterations){

	double localweight;
	double previousweight;
	int iteration = 0;

	localweight = this->initialize();

	cout << "minimizing NNF..." << endl;

	//debug tmpflow
	if(debug == true){
		irtkRealImage *flowtmp = new irtkRealImage(target->GetX(),target->GetY(),target->GetZ(),5);

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

	while(((previousweight-localweight)>0.0001 || iteration == 0) && iteration < maxiterations){
		cout << "iteration: " << iteration << " ";
		previousweight = localweight;
		localweight = this->minimizeflow();
		cout << "total distance " << localweight << endl;
		iteration++;
	}
	cout << "done" << endl;

	cout << "EM step... ";
	this->EMstep();
	cout << "done" << endl;
}

/// Generate HD image
void irtkPatchMatch::generateImage(){
	this->initialize();
	this->EMstep();
}

/// initialize field's weight
double irtkPatchMatch::initialize(){
	cout << "update gradients" << endl;
	// Compute spatial gradient of source image
	irtkGradientImageFilter<short> gradient(irtkGradientImageFilter<short>::GRADIENT_VECTOR);
	irtkGreyImage tmp = *target;
	gradient.SetInput (&tmp);
	gradient.SetOutput(targetgradient);
	gradient.SetPadding(-1);
	gradient.Run();

	if(debug == true){
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

int irtkPatchMatch::randomlink(int i, int j, int k, int o, int index){
	int x, y, z, n, ti, tj, tk, tn, wi, wj, wk, count;
	if(index < 0)
		index = target->VoxelToIndex(i,j,k);
	double dp;

	count = 0;


	if(nnfs[index][o].weight > 0){

		if(randomrate > 0){
			wi = target->GetX()*randomrate;
			wj = target->GetY()*randomrate;
			wk = target->GetZ()*randomrate;

			wi = min(wi,min(wj,wk));

			if(wi < 1)
				wi = 1;
		}else{
			wi = 0;
		}

		ti = nnfs[index][o].x;
		tj = nnfs[index][o].y;
		tk = nnfs[index][o].z;
		tn = nnfs[index][o].n;

		while(wi>0) {

			n = tn;

			x = ti + rand()%(2*wi+1)-wi;
			y = tj + rand()%(2*wi+1)-wi;
			z = tk + rand()%(2*wi+1)-wi;

			if(x < 0) x = 0;
			if(y < 0) y = 0;
			if(z < 0) z = 0;
			if(x > sources[n]->GetX() - 1)
				x = sources[n]->GetX() - 1;
			if(y > sources[n]->GetY() - 1)
				y = sources[n]->GetY() - 1;
			if(z > sources[n]->GetZ() - 1)
				z = sources[n]->GetZ() - 1;		

			dp = this->distance(i,j,k,x,y,z,n);

			if(dp < nnfs[index][o].weight){
				nnfs[index][o].x = x;
				nnfs[index][o].y = y;
				nnfs[index][o].z = z;
				nnfs[index][o].n = n;
				nnfs[index][o].weight = dp;
				count ++;
			}

			wi/=2;
		}

		if(randomrate > 0){
			wi = target->GetX()*randomrate;
			wj = target->GetY()*randomrate;
			wk = target->GetZ()*randomrate;

			wi = min(wi,min(wj,wk));

			if(wi < 1)
				wi = 1;
		}else{
			wi = 0;
		}

		while(wi>0) {

			n = rand()%nimages;

			x = ti + rand()%(2*wi+1)-wi;
			y = tj + rand()%(2*wi+1)-wi;
			z = tk + rand()%(2*wi+1)-wi;

			if(x < 0) x = 0;
			if(y < 0) y = 0;
			if(z < 0) z = 0;
			if(x > sources[n]->GetX() - 1)
				x = sources[n]->GetX() - 1;
			if(y > sources[n]->GetY() - 1)
				y = sources[n]->GetY() - 1;
			if(z > sources[n]->GetZ() - 1)
				z = sources[n]->GetZ() - 1;		

			/*if(i >= radius && i < target->GetX() - radius
			&& j >= radius && j < target->GetY() - radius
			&& k >= radius && k < target->GetZ() - radius
			&& x >= radius && x < sources[n]->GetX() - radius
			&& y >= radius && y < sources[n]->GetY() - radius
			&& z >= radius && z < sources[n]->GetZ() - radius)
			dp = this->distancefast(i,j,k,x,y,z,n);

			else*/
			dp = this->distance(i,j,k,x,y,z,n);

			if(dp < nnfs[index][o].weight){
				nnfs[index][o].x = x;
				nnfs[index][o].y = y;
				nnfs[index][o].z = z;
				nnfs[index][o].n = n;
				nnfs[index][o].weight = dp;
				count ++;
			}

			wi/=2;
		}
	}
	return count;
}

/// vote patch
void irtkPatchMatch::votepatch(int xt, int yt, int zt, int xs, int ys, int zs, int n, double weight){
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
/// vote label
void irtkPatchMatch::votelabel(int xt, int yt, int zt, int xs, int ys, int zs, int n, double weight){
	int index1;
	for(int k = - radius; k <= radius; k++){
		for(int j = - radius; j <= radius; j++){
			for(int i = - radius; i <= radius; i++){

				if(xt + i < 0 || xt + i > target->GetX() - 1
					|| xs + i < 0 || xs + i > labels[n]->GetX() - 1
					|| yt + j < 0 || yt + j > target->GetY() - 1
					|| ys + j < 0 || ys + j > labels[n]->GetY() - 1
					|| zt + k < 0 || zt + k > target->GetZ() - 1
					|| zs + k < 0 || zs + k > labels[n]->GetZ() - 1){

						// do nothing

				}else{

					index1 = target->VoxelToIndex(xt+i,yt+j,zt+k);

					votelabels[index1][labels[n]->Get(xs+i,ys+j,zs+k) - minlabel] += weight;

				}
			}
		}
	}
}

/// expectation maximization
void irtkPatchMatch::EMstep(){
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

void irtkPatchMatch::generateLabels(irtkGreyImage *label, irtkGreyImage **labels){
	/// Expectation create the vote matrix
	int index = 0, x, y, z, n, o,NumberOfLabels;
	double weight;
	this->labels = labels;

	label->Initialize(target->GetImageAttributes());

	//initialization
	labels[0]->GetMinMax(&minlabel,&maxlabel);
	NumberOfLabels = maxlabel - minlabel + 1;
	votelabels = new double*[label->GetNumberOfVoxels()];
	for(int i = 0; i < label->GetNumberOfVoxels(); i++){
		votelabels[i] = new double[NumberOfLabels];
		for(int j = 0; j < NumberOfLabels; j++){
			votelabels[i][j] = 0;
		}
	}

	for(int k = 0; k < label->GetZ(); k++){
		for(int j = 0; j < label->GetY(); j++){
			for(int i = 0; i < label->GetX(); i++){
				for(int o = 0; o < nneighbour; o++){
					x = nnfs[index][o].x;
					y = nnfs[index][o].y;
					z = nnfs[index][o].z;
					n = nnfs[index][o].n;
					weight = nnfs[index][o].weight;

					weight = 16/(weight + 1);
					if(weight > 1)
						weight = 1;
					else
						weight = weight * weight;

					this->votelabel(i,j,k,x,y,z,n,weight);
				}
				index++;
			}
		}
	}

	/// Maximization non local mean fusion
	index = 0;
	for(int k = 0; k < label->GetZ(); k++){
		for(int j = 0; j < label->GetY(); j++){
			for(int i = 0; i < label->GetX(); i++){
				double maxiumweight = 0;
				short maxiumproblabel = minlabel;
				for(int t = 0; t < NumberOfLabels; t++){
					if(votelabels[index][t] > maxiumweight){
						maxiumweight = votelabels[index][t];
						maxiumproblabel = minlabel + t;
					}
				}
				label->PutAsDouble(i,j,k,maxiumproblabel);
				index++;
			}
		}
	}
}

/// find minum flow
double irtkPatchMatch::minimizeflow(){
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
				//propergate from x - 1
				if(i < target->GetX() - 1){
					x = i + 1;
					y = j;
					z = k;
					count += this->propergate(x,y,z,i,j,k,-1,0,0,index1,index);
				}
				//propergate from y - 1
				if(j < target->GetY() - 1){
					x = i;
					y = j + 1;
					z = k;
					count += this->propergate(x,y,z,i,j,k,0,-1,0,index2,index);
				}
				//propergate from z - 1
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
					sumweight += this->nnfs[index][o].weight;
				}
				index++;
			}
		}
	}

	cout << "number of fields changed: " << fcount << " " << bcount << " " << rcount << " ";

	if(debug == true){
		/// Debug write flow to a image using the dynamicid
		irtkRealImage *flowtmp = new irtkRealImage(target->GetX(),target->GetY(),target->GetZ(),5);

		index = 0;
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
		sprintf(buffer, "flow%d_%d_%d.nii.gz", fcount, bcount, rcount);
		flowtmp->Write(buffer);
		delete flowtmp;
	}

	return log(sumweight);
}

/// propergate from one to another using the offsets
int irtkPatchMatch::propergate(int x, int y, int z, int i, int j, int k, int offsetx, int offsety, int offestz, int index1, int index2){
	int count, m;
	double dp;
	// get the index of the nnfs
	if(index1 < 0)
		index1 = target->VoxelToIndex(x,y,z);
	if(index2 < 0)
		index2 = target->VoxelToIndex(i,j,k);
	count = 0;
	for(int n = 0; n < nneighbour; n++){
		if(nnfs[index2][n].weight > 0){
			x = nnfs[index1][n].x + offsetx;
			y = nnfs[index1][n].y + offsety;
			z = nnfs[index1][n].z + offestz;
			m = nnfs[index1][n].n;



			/*if(i >= radius && i < target->GetX() - radius
			&& j >= radius && j < target->GetY() - radius
			&& k >= radius && k < target->GetZ() - radius
			&& x >= radius && x < sources[m]->GetX() - radius
			&& y >= radius && y < sources[m]->GetY() - radius
			&& z >= radius && z < sources[m]->GetZ() - radius)
			dp = this->distancefast(i,j,k,x,y,z,m);

			else*/
			dp = this->distance(i,j,k,x,y,z,m);

			// check if new distance is smaller than any of the old distances, if so, replace
			for(int o = 0; o < nneighbour; o++){
				if(dp < nnfs[index2][o].weight){
					nnfs[index2][o].x = x;
					nnfs[index2][o].y = y;
					nnfs[index2][o].z = z;
					nnfs[index2][o].n = m;
					nnfs[index2][o].weight = dp;
					count ++;
					break;
				}
			}
		}
	}
	return count;
}

/// calculate distance between patches
double irtkPatchMatch::distance(int x1, int y1, int z1, int x2, int y2, int z2, int n){
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

	for(t = 0; t < slevels; t++){
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
		increase = increase * 2;
		tmpradius = tmpradius * 2;
	}

	if(count == 0)
		return maxdistance - 1;
	else
		return dif/count;
}

void irtkPatchMatch::createsearchimages(){
	int i,j,k,n;
	for(n = 0; n < nimages; n++){
		for(k = 0; k < sources[n]->GetZ(); k++){
			for(j = 0; j < sources[n]->GetY(); j++){
				for(i = 0; i < sources[n]->GetX(); i++){
					search[n]->Put(i,j,k,this->checkissearch(i,j,k,n));
				}
			}
		}
	}
}

/// check if it is search
int irtkPatchMatch::checkissearch(int x2, int y2, int z2, int n){
	int i,j,k,i2,j2,k2;
	short values;
	for(k = - this->radius; k <= this->radius; k++){
		k2 = k + z2;
		for(j = -this->radius; j <= this->radius; j++){
			j2 = j + y2;
			for(i = -this->radius; i <= this->radius; i++){

				i2 = i + x2;

				if(i2 < 0 || (i2 > sources[n]->GetX() - 1)
					|| j2 < 0 || (j2 > sources[n]->GetY() - 1)
					|| k2 < 0 || (k2 > sources[n]->GetZ() - 1)){

				}else{
					values = sources[n]->Get(i2,j2,k2);
					if(values < 0){
						return 0;
					}
				}
			}
		}	
	}

	return 1;
}

void irtkPatchMatch::upSampleNNF(irtkPatchMatch *reference)
{
	int i,j,k, n, o, p, index1,index2;
	double x, y, z, weight, offsetx, offsety, offsetz;
	NearstNeighbor** targetnnfs = reference->getNearstNeighbor();

	index1 = 0;
	for(k = 0; k < reference->getTarget()->GetZ(); k++){
		for(j = 0; j < reference->getTarget()->GetY(); j++){
			for(i = 0; i < reference->getTarget()->GetX(); i++){
				x = i;
				y = j;
				z = k;
				reference->getTarget()->ImageToWorld(x,y,z);
				target->WorldToImage(x,y,z);
				if(x < 0) x = 0;
				if(x > target->GetX() - 1) x = target->GetX() - 1;
				if(y < 0) y = 0;
				if(y > target->GetY() - 1) y = target->GetY() - 1;
				if(z < 0) z = 0;
				if(z > target->GetZ() - 1) z = target->GetZ() - 1;

				index2 = target->VoxelToIndex(round(x),round(y),round(z));
				// remember the offset in the image coordinate
				offsetx = x - round(x);
				offsety = y - round(y);
				offsetz = z - round(z);

				for(n = 0; n < nneighbour; n++){

					x = nnfs[index2][n].x + offsetx;
					y = nnfs[index2][n].y + offsety;
					z = nnfs[index2][n].z + offsetz;
					p = nnfs[index2][n].n;
					weight = nnfs[index2][n].weight;

					sources[p]->ImageToWorld(x,y,z);
					reference->getSource(p)->WorldToImage(x,y,z);

					if(x < 0) x = 0;
					if(x > reference->getSource(p)->GetX() - 1) 
						x = reference->getSource(p)->GetX() - 1;
					if(y < 0) y = 0;
					if(y > reference->getSource(p)->GetY() - 1) 
						y = reference->getSource(p)->GetY() - 1;
					if(z < 0) z = 0;
					if(z > reference->getSource(p)->GetZ() - 1) 
						z = reference->getSource(p)->GetZ() - 1;

					targetnnfs[index1][n].x = round(x);
					targetnnfs[index1][n].y = round(y);
					targetnnfs[index1][n].z = round(z);
					targetnnfs[index1][n].n = p;
					targetnnfs[index1][n].weight = weight;

				}

				index1++;
			}
		}
	}
}

void irtkPatchMatch::setDecimatedImage(irtkGreyImage *input)
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