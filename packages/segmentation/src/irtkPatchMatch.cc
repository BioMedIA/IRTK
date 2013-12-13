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

irtkPatchMatch::irtkPatchMatch(irtkGreyImage *target, irtkGreyImage *source, int radius, int nneighbour){

	this->target = target;
	this->source = source;
	this->radius = radius;
	this->nneighbour = nneighbour;
	this->randomrate = 0.5;
	this->debug = false;
	short minvalue, maxvalue;
	target->GetMinMax(&minvalue, &maxvalue);
	maxdistance = maxvalue - minvalue;

	if(target->GetX() != source->GetX()
		|| target->GetY() != source->GetY()
		|| target->GetZ() != source->GetZ()
		|| target->GetT() != source->GetT()){
			cerr << "Image dimention does not match " << endl;
			cerr << "Target image: " << endl;
			target->Print();
			cerr << "Source image: " << endl;
			source->Print();
			exit(1);
	}

	//create NNF field
	this->nnfs = new NearstNeighborFlow*[target->GetNumberOfVoxels()];
	for(int i = 0; i < target->GetNumberOfVoxels(); i++){
		nnfs[i] = new NearstNeighborFlow[this->nneighbour];
	}

	cout.flush();

	this->initialguess();
}

void irtkPatchMatch::initialguess(){
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
					nnfs[index][n].x = 0;
					nnfs[index][n].y = 0;
					nnfs[index][n].z = 0;
					nnfs[index][n].weight = maxdistance;
				}
				index++;
			}
		}
	}
	cout << "done" << endl;
}

irtkPatchMatch::~irtkPatchMatch(){
	for(int i = 0; i < target->GetNumberOfVoxels(); i++){
		delete []nnfs[i];
	}
	delete []nnfs;
}

void irtkPatchMatch::run(int maxiterations){

	double localweight;
	double previousweight;
	int iteration = 0;

	localweight = this->initialize();

	cout << "minimizing NNF..." << endl;

	//debug tmpflow
	if(debug == true){
		irtkImageAttributes atr = target->GetImageAttributes();
		atr._t = 4;
		irtkRealImage *flowtmp = new irtkRealImage(atr);

		int index = 0;
		for(int k = 0; k < target->GetZ(); k++){
			for(int j = 0; j < target->GetY(); j++){
				for(int i = 0; i < target->GetX(); i++){
					flowtmp->Put(i,j,k,0,this->nnfs[index][0].x);
					flowtmp->Put(i,j,k,1,this->nnfs[index][0].y);
					flowtmp->Put(i,j,k,2,this->nnfs[index][0].z);
					flowtmp->Put(i,j,k,3,this->nnfs[index][0].weight);
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
double irtkPatchMatch::initialize(){

	if(debug == true){
		target->Write("targetincode.nii.gz");
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
						i+nnfs[index][n].x,j+nnfs[index][n].y,k+nnfs[index][n].z);

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
	int x, y, z, ti, tj, tk, wi, wj, wk, count;
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

		ti = nnfs[index][o].x + i;
		tj = nnfs[index][o].y + j;
		tk = nnfs[index][o].z + k;

		while(wi>0) {

			x = ti + rand()%(2*wi+1)-wi;
			y = tj + rand()%(2*wi+1)-wi;
			if(target->GetZ() > 1)
				z = tk + rand()%(2*wi+1)-wi;
			else
				z = tk;

			if(x < 0) x = 0;
			if(y < 0) y = 0;
			if(z < 0) z = 0;
			if(x > source->GetX() - 1)
				x = source->GetX() - 1;
			if(y > source->GetY() - 1)
				y = source->GetY() - 1;
			if(z > source->GetZ() - 1)
				z = source->GetZ() - 1;		

			dp = this->distance(i,j,k,x,y,z);

			if(dp < nnfs[index][o].weight){
				nnfs[index][o].x = x - i;
				nnfs[index][o].y = y - j;
				nnfs[index][o].z = z - k;
				nnfs[index][o].weight = dp;
				count ++;
			}
			wi/=2;
		}
	}
	return count;
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
		int index = 0, x, y, z, o;

		irtkImageAttributes atr = target->GetImageAttributes();
		atr._t = 4;
		irtkRealImage *flowtmp = new irtkRealImage(atr);

		index = 0;
		for(int k = 0; k < target->GetZ(); k++){
			for(int j = 0; j < target->GetY(); j++){
				for(int i = 0; i < target->GetX(); i++){
					flowtmp->Put(i,j,k,0,this->nnfs[index][0].x);
					flowtmp->Put(i,j,k,1,this->nnfs[index][0].y);
					flowtmp->Put(i,j,k,2,this->nnfs[index][0].z);
					flowtmp->Put(i,j,k,3,this->nnfs[index][0].weight);
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
double irtkPatchMatch::minimizeflowwithdebug(){
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
void irtkPatchMatch::voteweight(int zd, int mode){
	/// Expectation create the vote matrix
	int index = 0, x, y, z, o;

	/// Debug write flow to a image using the dynamicid

	irtkImageAttributes atr = target->GetImageAttributes();
	atr._t = 4;
	irtkRealImage *flowtmp = new irtkRealImage(atr);

	index = 0;
	for(int k = 0; k < target->GetZ(); k++){
		for(int j = 0; j < target->GetY(); j++){
			for(int i = 0; i < target->GetX(); i++){
				flowtmp->Put(i,j,k,0,this->nnfs[index][0].x);
				flowtmp->Put(i,j,k,1,this->nnfs[index][0].y);
				flowtmp->Put(i,j,k,2,this->nnfs[index][0].z);
				flowtmp->Put(i,j,k,3,this->nnfs[index][0].weight);
				index++;
			}
		}
	}

	char buffer[255];
	sprintf(buffer, "flowdebug%d_%d.nii.gz", mode, zd);
	flowtmp->Write(buffer);
	delete flowtmp;
}

void irtkPatchMatch::outputgraph(char* name){
	int index1, index2;
	double distance;
	ofstream fout(name,ios::app);	  
	// Output links between neighbors in the image
	// Target
	for(int k = 0; k < target->GetZ(); k++){
		for(int j = 0; j < target->GetY(); j++){
			for(int i = 0; i < target->GetX(); i++){
				index1 = target->VoxelToIndex(i,j,k);
				if(i > 0){
					index2 = target->VoxelToIndex(i - 1, j, k);
					distance = this->selfdistance(i,j,k,i-1,j,k);
					fout << index1 << " " << index2 << " " << distance << endl;
				}
				if(j > 0){
					index2 = target->VoxelToIndex(i, j - 1, k);
					distance = this->selfdistance(i,j,k,i,j-1,k);
					fout << index1 << " " << index2 << " " << distance << endl;
				}
				if(k > 0){
					index2 = target->VoxelToIndex(i, j, k - 1);
					distance = this->selfdistance(i,j,k,i,j,k-1);
					fout << index1 << " " << index2 << " " << distance << endl;
				}
				if(i < target->GetX() - 1){
					index2 = target->VoxelToIndex(i + 1, j, k);
					distance = this->selfdistance(i,j,k,i+1,j,k);
					fout << index1 << " " << index2 << " " << distance << endl;
				}
				if(j < target->GetY() - 1){
					index2 = target->VoxelToIndex(i, j + 1, k);
					distance = this->selfdistance(i,j,k,i,j+1,k);
					fout << index1 << " " << index2 << " " << distance << endl;
				}
				if(k < target->GetZ() - 1){
					index2 = target->VoxelToIndex(i, j, k + 1);
					distance = this->selfdistance(i,j,k,i,j,k+1);
					fout << index1 << " " << index2 << " " << distance << endl;
				}
			}
		}
	}

	// Source
	for(int k = 0; k < source->GetZ(); k++){
		for(int j = 0; j < source->GetY(); j++){
			for(int i = 0; i < source->GetX(); i++){
				index1 = source->VoxelToIndex(i,j,k) + target->GetNumberOfVoxels()/target->GetT();
				if(i > 0){
					index2 = source->VoxelToIndex(i - 1, j, k) + target->GetNumberOfVoxels()/target->GetT();
					distance = this->selfdistance(i,j,k,i-1,j,k,1);
					fout << index1 << " " << index2 << " " << distance << endl;
				}
				if(j > 0){
					index2 = source->VoxelToIndex(i, j - 1, k) + target->GetNumberOfVoxels()/target->GetT();
					distance = this->selfdistance(i,j,k,i,j-1,k,1);
					fout << index1 << " " << index2 << " " << distance << endl;
				}
				if(k > 0){
					index2 = source->VoxelToIndex(i, j, k - 1) + target->GetNumberOfVoxels()/target->GetT();
					distance = this->selfdistance(i,j,k,i,j,k-1,1);
					fout << index1 << " " << index2 << " " << distance << endl;
				}
				if(i < source->GetX() - 1){
					index2 = source->VoxelToIndex(i + 1, j, k) + target->GetNumberOfVoxels()/target->GetT();
					distance = this->selfdistance(i,j,k,i+1,j,k,1);
					fout << index1 << " " << index2 << " " << distance << endl;
				}
				if(j < source->GetY() - 1){
					index2 = source->VoxelToIndex(i, j + 1, k) + target->GetNumberOfVoxels()/target->GetT();
					distance = this->selfdistance(i,j,k,i,j+1,k,1);
					fout << index1 << " " << index2 << " " << distance << endl;
				}
				if(k < source->GetZ() - 1){
					index2 = source->VoxelToIndex(i, j, k + 1) + target->GetNumberOfVoxels()/target->GetT();
					distance = this->selfdistance(i,j,k,i,j,k+1,1);
					fout << index1 << " " << index2 << " " << distance << endl;
				}
			}
		}
	}

	// Output links between images from nnfs
	index1 = 0;
	for(int k = 0; k < target->GetZ(); k++){
		for(int j = 0; j < target->GetY(); j++){
			for(int i = 0; i < target->GetX(); i++){
				index2 = source->VoxelToIndex(this->nnfs[index1][0].x, this->nnfs[index1][0].y,
					this->nnfs[index1][0].z) + target->GetNumberOfVoxels()/target->GetT();
				distance = this->nnfs[index1][0].weight;
				fout << index1 << " " << index2 << " " << distance << endl;
				index1++;
			}
		}
	}

	fout.close();
}

void irtkPatchMatch::outputmap(char* name){
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
				flowtmp->Put(i,j,k,3,this->nnfs[index][0].weight);
				flowtmp->Put(i,j,k,4,this->source->GetAsDouble(
					i+this->nnfs[index][0].x,j+this->nnfs[index][0].y,k+this->nnfs[index][0].z));
				index++;
			}
		}
	}

	flowtmp->Write(name);
	delete flowtmp;
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
			x = i + nnfs[index1][n].x + offsetx;
			y = j + nnfs[index1][n].y + offsety;
			z = k + nnfs[index1][n].z + offestz;

			dp = this->distance(i,j,k,x,y,z);

			// check if new distance is smaller than any of the old distances, if so, replace
			for(int o = 0; o < nneighbour; o++){
				if(dp < nnfs[index2][o].weight){
					nnfs[index2][o].x = nnfs[index1][n].x;
					nnfs[index2][o].y = nnfs[index1][n].y;
					nnfs[index2][o].z = nnfs[index1][n].z;
					nnfs[index2][o].weight = dp;
					count ++;
					break;
				}
			}
		}
	}
	return count;
}

double irtkPatchMatch::distance2D(int x1, int y1, int x2, int y2){
	int i,j,k,o,i1,j1,i2,j2,t,count,g,tmpradius,increase;
	double dif = 0, tmp;
	short value1, value2, values;
	count = 0;
	increase = 1;
	tmpradius = this->radius;

	//check is search
	if(x2 < 0 || x2 > source->GetX() - 1
		|| y2 < 0 || y2 > source->GetY() - 1
		|| source->Get(x2,y2,0) < 0)
		return maxdistance;

	if(target->Get(x1,y1,0) < 0){
		return maxdistance - 1;
	}

	for(o = 0; o < target->GetT(); o++){
		k = 0;
		for(j = -tmpradius; j <= tmpradius; j+=increase){
			j1 = j + y1;
			j2 = j + y2;
			for(i = -tmpradius; i <= tmpradius; i+=increase){
				i1 = i + x1;
				i2 = i + x2;
				if(i1 < 0|| i2 < 0 ||
					(i1 > target->GetX() - 1) 
					|| (i2 > source->GetX() - 1)
					|| j1 < 0|| j2 < 0
					|| (j1 > target->GetY() - 1) 
					|| (j2 > source->GetY() - 1)){
						dif += maxdistance;
						count ++;
				}else{
					value1 = target->Get(i1,j1,0,o);
					if(value1 >= 0){		
						values = source->Get(i2,j2,0,o);
						tmp = double(value1 - values);
						dif += sqrt(tmp*tmp);
						count++;
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

double irtkPatchMatch::distance(int x1, int y1, int z1, int x2, int y2, int z2){
	if(target->GetZ() > 1)
		return this->distance3D(x1, y1, z1, x2, y2, z2);
	else
		return this->distance2D(x1, y1, x2, y2);
}

double irtkPatchMatch::distance3D(int x1, int y1, int z1, int x2, int y2, int z2){
	int i,j,k,o,i1,j1,k1,i2,j2,k2,t,count,g,tmpradius,increase;
	double dif = 0, tmp;
	short value1, value2, values;
	count = 0;
	increase = 1;
	tmpradius = this->radius;

	//check is search
	if(x2 < 0 || x2 > source->GetX() - 1
		|| y2 < 0 || y2 > source->GetY() - 1
		|| z2 < 0 || z2 > source->GetZ() - 1
		|| source->Get(x2,y2,z2) < 0)
		return maxdistance;

	if(target->Get(x1,y1,z1) < 0){
		return maxdistance - 1;
	}

	for(o = 0; o < target->GetT(); o++){
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
						|| (i2 > source->GetX() - 1)
						|| j1 < 0|| j2 < 0
						|| (j1 > target->GetY() - 1) 
						|| (j2 > source->GetY() - 1)
						|| k1 < 0|| k2 < 0
						|| (k1 > target->GetZ() - 1) 
						|| (k2 > source->GetZ() - 1)){
							dif += maxdistance;
							count ++;
					}else{
						value1 = target->Get(i1,j1,k1,o);
						if(value1 >= 0){		
							values = source->Get(i2,j2,k2,o);
							tmp = double(value1 - values);
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

double irtkPatchMatch::selfdistance2D(int x1, int y1, int x2, int y2, int mode){
	int i,j,k,o,i1,j1,i2,j2,t,count,g,tmpradius,increase;
	double dif = 0, tmp;
	short value1, value2, values;
	count = 0;
	increase = 1;
	tmpradius = this->radius;

	if(mode == 0){
		//check is search
		if(x2 < 0 || x2 > target->GetX() - 1
			|| y2 < 0 || y2 > target->GetY() - 1
			|| target->Get(x2,y2,0) < 0)
			return maxdistance;

		if(target->Get(x1,y1,0) < 0){
			return maxdistance - 1;
		}

		for(o = 0; o < target->GetT(); o++){
			k = 0;
			for(j = -tmpradius; j <= tmpradius; j+=increase){
				j1 = j + y1;
				j2 = j + y2;
				for(i = -tmpradius; i <= tmpradius; i+=increase){
					i1 = i + x1;
					i2 = i + x2;
					if(i1 < 0|| i2 < 0 ||
						(i1 > target->GetX() - 1) 
						|| (i2 > target->GetX() - 1)
						|| j1 < 0|| j2 < 0
						|| (j1 > target->GetY() - 1) 
						|| (j2 > target->GetY() - 1)){
							dif += maxdistance;
							count ++;
					}else{
						value1 = target->Get(i1,j1,0,o);
						if(value1 >= 0){		
							values = target->Get(i2,j2,0,o);
							tmp = double(value1 - values);
							dif += sqrt(tmp*tmp);
							count++;
						}
					}	
				}
			}
		}
	}else{
		//check is search
		if(x2 < 0 || x2 > source->GetX() - 1
			|| y2 < 0 || y2 > source->GetY() - 1
			|| source->Get(x2,y2,0) < 0)
			return maxdistance;

		if(source->Get(x1,y1,0) < 0){
			return maxdistance - 1;
		}

		for(o = 0; o < source->GetT(); o++){
			k = 0;
			for(j = -tmpradius; j <= tmpradius; j+=increase){
				j1 = j + y1;
				j2 = j + y2;
				for(i = -tmpradius; i <= tmpradius; i+=increase){
					i1 = i + x1;
					i2 = i + x2;
					if(i1 < 0|| i2 < 0 ||
						(i1 > source->GetX() - 1) 
						|| (i2 > source->GetX() - 1)
						|| j1 < 0|| j2 < 0
						|| (j1 > source->GetY() - 1) 
						|| (j2 > source->GetY() - 1)){
							dif += maxdistance;
							count ++;
					}else{
						value1 = source->Get(i1,j1,0,o);
						if(value1 >= 0){		
							values = source->Get(i2,j2,0,o);
							tmp = double(value1 - values);
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

double irtkPatchMatch::selfdistance(int x1, int y1, int z1, int x2, int y2, int z2, int mode){
	if(target->GetZ() > 1)
		return this->selfdistance3D(x1, y1, z1, x2, y2, z2, mode);
	else
		return this->selfdistance2D(x1, y1, x2, y2, mode);
}

double irtkPatchMatch::selfdistance3D(int x1, int y1, int z1, int x2, int y2, int z2, int mode){
	int i,j,k,o,i1,j1,k1,i2,j2,k2,t,count,g,tmpradius,increase;
	double dif = 0, tmp;
	short value1, value2, values;
	count = 0;
	increase = 1;
	tmpradius = this->radius;

	if(mode == 0){
		//check is search
		if(x2 < 0 || x2 > target->GetX() - 1
			|| y2 < 0 || y2 > target->GetY() - 1
			|| z2 < 0 || z2 > target->GetZ() - 1
			|| target->Get(x2,y2,z2) < 0)
			return maxdistance;

		if(target->Get(x1,y1,z1) < 0){
			return maxdistance - 1;
		}

		for(o = 0; o < target->GetT(); o++){
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
							|| (i2 > target->GetX() - 1)
							|| j1 < 0|| j2 < 0
							|| (j1 > target->GetY() - 1) 
							|| (j2 > target->GetY() - 1)
							|| k1 < 0|| k2 < 0
							|| (k1 > target->GetZ() - 1) 
							|| (k2 > target->GetZ() - 1)){
								dif += maxdistance;
								count ++;
						}else{
							value1 = target->Get(i1,j1,k1,o);
							if(value1 >= 0){		
								values = target->Get(i2,j2,k2,o);
								tmp = double(value1 - values);
								dif += sqrt(tmp*tmp);
								count++;
							}
						}	
					}
				}
			}
		}
	}else{
		//check is search
		if(x2 < 0 || x2 > source->GetX() - 1
			|| y2 < 0 || y2 > source->GetY() - 1
			|| z2 < 0 || z2 > source->GetZ() - 1
			|| source->Get(x2,y2,z2) < 0)
			return maxdistance;

		if(source->Get(x1,y1,z1) < 0){
			return maxdistance - 1;
		}

		for(o = 0; o < source->GetT(); o++){
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
							(i1 > source->GetX() - 1) 
							|| (i2 > source->GetX() - 1)
							|| j1 < 0|| j2 < 0
							|| (j1 > source->GetY() - 1) 
							|| (j2 > source->GetY() - 1)
							|| k1 < 0|| k2 < 0
							|| (k1 > source->GetZ() - 1) 
							|| (k2 > source->GetZ() - 1)){
								dif += maxdistance;
								count ++;
						}else{
							value1 = source->Get(i1,j1,k1,o);
							if(value1 >= 0){		
								values = source->Get(i2,j2,k2,o);
								tmp = double(value1 - values);
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