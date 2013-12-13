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

irtkMAPatchMatch::irtkMAPatchMatch(irtkGreyImage *target, irtkGreyImage **source, int radius, int nimages, int nneighbour){

	this->target = target;
	this->sources = source;
	this->radius = radius;
	this->nimages = nimages;
	this->nneighbour = nneighbour;
	this->randomrate = 0.5;
	this->debug = false;
	short minvalue, maxvalue;
	target->GetMinMax(&minvalue, &maxvalue);
	maxdistance = maxvalue - minvalue;

	//create NNF field
	this->nnfs = new NearstNeighbor*[target->GetNumberOfVoxels()];
	for(int i = 0; i < target->GetNumberOfVoxels(); i++){
		nnfs[i] = new NearstNeighbor[this->nneighbour];
	}

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

	cout.flush();

	this->initialguess();
}

void irtkMAPatchMatch::initialguess(){
	//populate NNF field with initial guess
	cout << "populating NNF field with spatial initial guess" << endl;
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
}

irtkMAPatchMatch::~irtkMAPatchMatch(){
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
}

void irtkMAPatchMatch::run(int maxiterations){

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
		if(debug == true && iteration == 0){
			localweight = this->minimizeflowwithdebug();
		}else{
			localweight = this->minimizeflow();
		}
		cout << "total distance " << localweight << endl;
		iteration++;
	}
	cout << "done" << endl;

}

/// initialize field's weight
double irtkMAPatchMatch::initialize(){
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

int irtkMAPatchMatch::randomlink(int i, int j, int k, int o, int index){
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

			if(target->GetZ() > 1)
				wi = min(wi,min(wj,wk));
			else
				wi = min(wi,wj);

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

			if(target->GetZ() > 1)
				z = tk + rand()%(2*wi+1)-wi;
			else
				z = tk;

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


/// find minum flow
double irtkMAPatchMatch::minimizeflow(){
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

		irtkImageAttributes atr = target->GetImageAttributes();
		atr._t = 5;
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

/// find minum flow
double irtkMAPatchMatch::minimizeflowwithdebug(){
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

/// vote weight matrix
void irtkMAPatchMatch::voteweight(int zd, int mode){
	/// Expectation create the vote matrix
	int index = 0, x, y, z, n, o;

	/// Debug write flow to a image using the dynamicid

	irtkImageAttributes atr = target->GetImageAttributes();
	atr._t = 5;
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
				index++;
			}
		}
	}

	char buffer[255];
	sprintf(buffer, "flowdebug%d_%d.nii.gz", mode, zd);
	flowtmp->Write(buffer);
	delete flowtmp;
}

void irtkMAPatchMatch::outputmap(char* name){
	/// Debug write flow to a image using the dynamicid
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

	flowtmp->Write(name);
	delete flowtmp;
}

/// propergate from one to another using the offsets
int irtkMAPatchMatch::propergate(int x, int y, int z, int i, int j, int k, int offsetx, int offsety, int offestz, int index1, int index2){
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

//wenzhe question for cuda impelementation
//1. how to include and check if cuda is on //ifdef HAS_CUDA?
//2. how to design this, host memory includes target, targetgradient, atlases, atlasesgradient, the mapping and weighting
//3. how to return the result double
/// calculate distance between patches



double irtkMAPatchMatch::distance(int x1, int y1, int z1, int x2, int y2, int z2, int n){
	int i,j,k,i1,j1,k1,i2,j2,k2,t,count,g,tmpradius,tmpradiusx, tmpradiusy, tmpradiusz, increase;
	double dif = 0, tmp;
	short value1, value2, values;
	count = 0;
	increase = 1;
	tmpradius = this->radius;

	tmpradiusx = tmpradius;
	tmpradiusy = tmpradius;
	if(target->GetZ() > 1)
		tmpradiusz = tmpradius;
	else
		tmpradiusz = 0;

	//check is search
	if(x2 < 0 || x2 > sources[n]->GetX() - 1
		|| y2 < 0 || y2 > sources[n]->GetY() - 1
		|| z2 < 0 || z2 > sources[n]->GetZ() - 1
		|| search[n]->Get(x2,y2,z2) < 1)
		return maxdistance;

	if(target->Get(x1,y1,z1) < 0){
		return maxdistance - 1;
	}

	for(k = - tmpradiusz; k <= tmpradiusz; k+=increase){
		k1 = k + z1;
		k2 = k + z2;
		for(j = -tmpradiusy; j <= tmpradiusy; j+=increase){
			j1 = j + y1;
			j2 = j + y2;
			for(i = -tmpradiusx; i <= tmpradiusx; i+=increase){
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
						dif += maxdistance*4;
						count += 4;
				}else{
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
			}
		}
	}

	if(count < 1)
		return maxdistance - 1;
	else
		return dif/count;
}

void irtkMAPatchMatch::createsearchimages(){
	int i,j,k,n;
	for(n = 0; n < nimages; n++){
		for(k = 0; k < sources[n]->GetZ(); k++){
			for(j = 0; j < sources[n]->GetY(); j++){
				for(i = 0; i < sources[n]->GetX(); i++){
					//TODO
					//search[n]->Put(i,j,k,this->checkissearch(i,j,k,n));
					search[n]->Put(i,j,k,1);
				}
			}
		}
	}
}

/// check if it is search
int irtkMAPatchMatch::checkissearch(int x2, int y2, int z2, int n){
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

void irtkMAPatchMatch::upSampleNNF(irtkMAPatchMatch *reference)
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