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
#include <irtkMeanFilter.h>

irtkMAPatchMatchSegmentation::irtkMAPatchMatchSegmentation
	(irtkGreyImage *target, irtkGreyImage **source, 
	irtkRealImage *targetdistance, irtkRealImage **sourcedistance,
	irtkGreyImage *label, irtkGreyImage **labels,
	int radius, int nimages, int nneighbour)
	:irtkMAPatchMatch(target,source,radius,nimages,nneighbour)
{
	//construct distance
	this->targetlabeldistance = targetdistance;
	this->sourcelabeldistance = sourcedistance;
	this->label = label;
	this->labels = labels;

	this->context_target = new irtkGreyImage(this->target->GetImageAttributes());

	irtkMeanFilter<irtkGreyPixel> *filter = new irtkMeanFilter<irtkGreyPixel>();
	filter->SetInput(target);
	filter->SetOutput(context_target);
	filter->SetkernelRadius(radius);
	filter->Run();

	this->context_sources = new irtkGreyImage*[nimages];
	for(int n = 0; n < nimages; n++){
		context_sources[n] = new irtkGreyImage(this->sources[n]->GetImageAttributes());
		filter->SetInput(sources[n]);
		filter->SetOutput(context_sources[n]);
		filter->Run();
	}

	delete filter;

	if(label == NULL || labels == NULL){
		cerr << "labels are empty, please input labels" << endl;
		exit(1);
	}

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

	this->radius_x = round(this->radius / target->GetXSize());
	this->radius_y = round(this->radius / target->GetYSize());
	this->radius_z = round(this->radius / target->GetZSize());

	if(this->radius_x < 1) this->radius_x = 1;
	if(this->radius_y < 1) this->radius_y = 1;
	if(this->radius_z < 1) this->radius_z = 1;

}

/// check if it is search
int irtkMAPatchMatchSegmentation::checkissearch(int x2, int y2, int z2, int n){
	int i,j,k,i2,j2,k2;
	short values;
	for(k = - this->radius_z; k <= this->radius_z; k++){
		k2 = k + z2;
		for(j = -this->radius_y; j <= this->radius_y; j++){
			j2 = j + y2;
			for(i = -this->radius_x; i <= this->radius_x; i++){

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

irtkMAPatchMatchSegmentation::~irtkMAPatchMatchSegmentation()
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

	delete context_target;
	delete []context_sources;

	if(distancenorm != NULL){
		delete []distancenorm;
	}
}

void irtkMAPatchMatchSegmentation::normalizedistances(){
	distancenorm = new double[targetlabeldistance->GetT()];
	double min,max;
	irtkRealPixel min_tmp, max_tmp;
	for(int t = 0; t < targetlabeldistance->GetT(); t++){
		targetlabeldistance->GetMinMax(&max_tmp,&min_tmp);
		min = static_cast<double>(min_tmp);
		max = static_cast<double>(max_tmp);
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

/// vote label
void irtkMAPatchMatchSegmentation::votelabel(int xt, int yt, int zt, int xs, int ys, int zs, int n, double weight){

	if( xs < 0 || xs > labels[n]->GetX() - 1
		|| ys < 0 || ys > labels[n]->GetY() - 1
		|| zs < 0 || zs > labels[n]->GetZ() - 1){

			// do nothing

	}else{

		votelabels[labels[n]->Get(xs,ys,zs) - minlabel] += weight;

	}
}

void irtkMAPatchMatchSegmentation::generateLabels(){
	/// Expectation create the vote matrix
	int index,x, y, z, n, o,NumberOfLabels;
	double weight;

	label->Initialize(target->GetImageAttributes());

	//initialization
	labels[0]->GetMinMax(&minlabel,&maxlabel);
	NumberOfLabels = maxlabel - minlabel + 1;
	votelabels = new double[NumberOfLabels];

	for(int k = 0; k < label->GetZ(); k++){
		for(int j = 0; j < label->GetY(); j++){
			for(int i = 0; i < label->GetX(); i++){

				//find maxium label
				double maxiumweight = 0;
				short maxiumproblabel = minlabel;
				if(target->Get(i,j,k) > -1){

					// initialize labels
					for(int l = 0; l < NumberOfLabels; l++){
						votelabels[l] = 0;
					}

					//vote labels with neighbours
					for(int offsetz = -radius_z; offsetz <= radius_z; offsetz++){
						for(int offsety = -radius_y; offsety <= radius_y; offsety++){
							for(int offsetx = -radius_x; offsetx <= radius_x; offsetx++){
								if(k+offsetz < 0 || k+offsetz >= label->GetZ()
									|| j+offsety < 0 || j+offsety >= label->GetY()
									|| i+offsetx < 0 || i+offsetx >= label->GetX()){
								}else{
									index = label->VoxelToIndex(i+offsetx,j+offsety,k+offsetz);
									for(int o = 0; o < nneighbour; o++){
										x = nnfs[index][o].x - offsetx;
										y = nnfs[index][o].y - offsety;
										z = nnfs[index][o].z - offsetz;
										n = nnfs[index][o].n;
										weight = nnfs[index][o].weight;

										weight = exp(-(weight/maxdistance*32));

										this->votelabel(i,j,k,x,y,z,n,weight);
									}
								}
							}
						}
					}

					for(int t = 0; t < NumberOfLabels; t++){
						if(votelabels[t] > maxiumweight){
							maxiumweight = votelabels[t];
							maxiumproblabel = minlabel + t;
						}
					}
				}

				//put label
				label->PutAsDouble(i,j,k,maxiumproblabel);
			}
		}
	}

	delete []votelabels;
}

/// find minum flow
double irtkMAPatchMatchSegmentation::minimizeflow(){

	if(debug == true)
		cout << "irtkMAPatchMatchSegmentation::minimizeflow" << endl;

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
		atr._t = 6;
		irtkRealImage *flowtmp = new irtkRealImage(atr);

		this->generateLabels();

		index = 0;
		for(int k = 0; k < target->GetZ(); k++){
			for(int j = 0; j < target->GetY(); j++){
				for(int i = 0; i < target->GetX(); i++){
					flowtmp->Put(i,j,k,0,this->nnfs[index][0].x);
					flowtmp->Put(i,j,k,1,this->nnfs[index][0].y);
					flowtmp->Put(i,j,k,2,this->nnfs[index][0].z);
					flowtmp->Put(i,j,k,3,this->nnfs[index][0].n);
					flowtmp->Put(i,j,k,4,this->nnfs[index][0].weight);
					flowtmp->Put(i,j,k,5,this->label->Get(i,j,k));
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
void irtkMAPatchMatchSegmentation::voteweight(int zd, int mode){

	if(debug == true)
		cout << "irtkMAPatchMatchSegmentation::voteweight" << endl;

	/// Expectation create the vote matrix
	int index = 0, x, y, z, n, o;
	double weight;

	this->generateLabels();

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
				flowtmp->Put(i,j,k,4,this->nnfs[index][0].weight);
				flowtmp->PutAsDouble(i,j,k,5,label->Get(i,j,k));
				index++;
			}
		}
	}

	char buffer[255];
	sprintf(buffer, "flowdebug%d_%d.nii.gz", mode, zd);
	flowtmp->Write(buffer);
	delete flowtmp;
}

void irtkMAPatchMatchSegmentation::initialguess(){
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
					//x = rand()%target->GetX();
					//y = rand()%target->GetY();
					//z = rand()%target->GetZ();
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

void irtkMAPatchMatchSegmentation::setWeight(double weight){
	if(targetlabeldistance != NULL)
		sweight = weight/targetlabeldistance->GetT();
}

void irtkMAPatchMatchSegmentation::setDirectionalRadius(double x, double y, double z){
	this->radius_x = round(x / target->GetXSize());
	this->radius_y = round(y / target->GetYSize());
	this->radius_z = round(z / target->GetZSize());

	if(this->radius_x < 1) this->radius_x = 1;
	if(this->radius_y < 1) this->radius_y = 1;
	if(this->radius_z < 1) this->radius_z = 1;
}

/// calculate distance between patches
double irtkMAPatchMatchSegmentation::distance(int x1, int y1, int z1, int x2, int y2, int z2, int n){

	int i,j,k,i1,j1,k1,i2,j2,k2,t,g;
	int tmpradius_x, tmpradius_y, tmpradius_z;
	int increase_x, increase_y, increase_z;
	double dif = 0, tmp, count;
	short value1, value2, values;
	count = 0;
	increase_x = increase_y = increase_z = 1;
	tmpradius_x = radius_x;
	tmpradius_y = radius_y;
	tmpradius_z = radius_z;

	if(target->GetZ() <= 1)
		tmpradius_z = 0;

	if(target->Get(x1,y1,z1) < 0){
		return maxdistance - 1;
	}

	for(k = - tmpradius_z; k <= tmpradius_z; k+=increase_z){
		k1 = k + z1;
		k2 = k + z2;
		for(j = -tmpradius_y; j <= tmpradius_y; j+=increase_y){
			j1 = j + y1;
			j2 = j + y2;
			for(i = -tmpradius_x; i <= tmpradius_x; i+=increase_x){

				i1 = i + x1;
				i2 = i + x2;

				if(i1 < 0 || (i1 > target->GetX() - 1) 
					|| j1 < 0 || (j1 > target->GetY() - 1) 
					|| k1 < 0 || (k1 > target->GetZ() - 1) ){
						//do not count
				}else if(i2 < 0 || (i2 > sources[n]->GetX() - 1)
					|| j2 < 0 || (j2 > sources[n]->GetY() - 1)
					|| k2 < 0 || (k2 > sources[n]->GetZ() - 1)){
						dif += maxdistance*4;
						count += 4;
						if(targetlabeldistance != NULL){
							dif += maxdistance*sweight*targetlabeldistance->GetT();
							count += sweight*targetlabeldistance->GetT();
						}
				}else{
					//distance between resolved image and atlases
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
	}

	//now contextual informations
	//tmpradius_x = radius_x*8;
	//tmpradius_y = radius_y*8;
	//tmpradius_z = radius_z*8;
	//increase_x = radius_x*2;
	//increase_y = radius_y*2;
	//increase_z = radius_z*2;
	//
	//short mean_value_target, mean_value_source;
	//mean_value_target = context_target->Get(x1,y1,z1);
	//mean_value_source = context_sources[n]->Get(x2,y2,z2);

	//for(k = - tmpradius_z; k <= tmpradius_z; k+=increase_z){
	//	k1 = k + z1;
	//	k2 = k + z2;
	//	for(j = -tmpradius_y; j <= tmpradius_y; j+=increase_y){
	//		j1 = j + y1;
	//		j2 = j + y2;
	//		for(i = -tmpradius_x; i <= tmpradius_x; i+=increase_x){
	//			i1 = i + x1;
	//			i2 = i + x2;
	//			if(i1 < 0 || (i1 > target->GetX() - 1) 
	//				|| j1 < 0 || (j1 > target->GetY() - 1) 
	//				|| k1 < 0 || (k1 > target->GetZ() - 1) ){
	//					//do not count
	//			}else if(i2 < 0 || (i2 > sources[n]->GetX() - 1)
	//				|| j2 < 0 || (j2 > sources[n]->GetY() - 1)
	//				|| k2 < 0 || (k2 > sources[n]->GetZ() - 1)){
	//					dif += maxdistance;
	//					count += 1;
	//			}else{
	//				//distance between resolved image and atlases
	//				//reconstructed image not just decimated image
	//				value1 = context_target->Get(i1,j1,k1);
	//				if(value1 >= 0){		
	//					value1 -= mean_value_target;
	//					values = context_sources[n]->Get(i2,j2,k2) - mean_value_source;
	//					tmp = double(value1 - values);
	//					dif += sqrt(tmp*tmp);
	//					count++;
	//				}
	//			}	
	//		}
	//	}
	//}

	if(count < 1)
		return maxdistance - 1;
	else
		return dif/count;
}
