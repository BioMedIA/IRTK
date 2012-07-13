/*
 * irtkBatchBasedSegmentation.cc
 *
 *  Created on: Jun 21, 2011
 *      Author: rw1008
 */

#include<irtkPatchBasedSegmentation.h>
#include<math.h>






irtkPatchBasedSegmentation::irtkPatchBasedSegmentation(irtkRealImage image, irtkRealImage ** atlases, irtkGreyImage ** labels, int nAtlases, int patchSize, int neighbourhoodsize){
	_image = image;
	_atlases = atlases;
	_labels = labels;
	irtkGreyPixel min_lab, max_lab;
	_labels[0]->GetMinMax(&min_lab, &max_lab);
	_availableLabels = new int[int(max_lab)];
	_nrLabels = int(max_lab)+1;
	_patchSize = patchSize;
	_neighbourhoodSize = neighbourhoodsize;
	_nAtlases = nAtlases;
	_segmentations = new irtkGreyImage*[nAtlases];
	_verbose = 0;
	_padding = -1;
	_wMin = new double[_nAtlases];
	_maxVal = 65000;
	for(int i = 0 ; i < _nAtlases; i++){
		_segmentations[i] = new irtkGreyImage;
		_segmentations[i]->Initialize(_labels[0]->GetImageAttributes());
		_wMin[i] = _maxVal;

	}
	_useMask = false;
	_patchMeasuresAvailable = false;
	_winnerTakesAll = false;


}

irtkPatchBasedSegmentation::~irtkPatchBasedSegmentation(){
	delete [] _availableLabels;
	for(int i = 0; i < _nAtlases; i++){
		delete _segmentations[i];
	}
	delete [] _segmentations;
}

int irtkPatchBasedSegmentation::GetNeighbourhoodSize(){
	return _neighbourhoodSize;
}



double irtkPatchBasedSegmentation::EvaluatePatch(int x_image, int y_image, int z_image, int x_atlas, int y_atlas, int z_atlas, int atlasNo){
	double w = 0;
	int N = 0;
	for(int x = -_patchSize; x <= _patchSize; x++){
		for(int y = -_patchSize; y <= _patchSize; y++){
			for(int z = -_patchSize; z <= _patchSize; z++){
					if(!((x_atlas + x ) < 0 || ( x_atlas + x ) >= _image.GetX() ||
						(y_atlas + y )< 0 || ( y_atlas + y ) >= _image.GetY() ||
						(z_atlas + z )< 0 || (z_atlas + z )>= _image.GetZ() ||
						(x_image + x )< 0 || (x_image + x )>= _image.GetX() ||
						(y_image + y )< 0 || (y_image + y ) >= _image.GetY() ||
						(z_image + z )< 0 || (z_image + z )>= _image.GetZ() )){
							int val_atlas = _atlases[atlasNo]->Get(x_atlas+x, y_atlas+y, z_atlas+z);
							int val_image = _image.Get(x_image+x, y_image+y, z_image+z);
							if(val_image > 0 && val_atlas > 0 ){
								w += (val_atlas - val_image) * (val_atlas - val_image);
								N++ ;
							}
					}
				}
			}
	}
	double ret_val = w / N;
	return ret_val;
}



void irtkPatchBasedSegmentation::EvaluateNeighbourhood(int x_image, int y_image, int z_image, int atlasNo, double *** tmpVals){

	int xCtr = 0;
	double my_i = _patchMean[0]->Get(x_image, y_image, z_image);
	double sigma_i = _patchStd[0]->Get(x_image, y_image, z_image);
	irtkRealImage i;
	int ctr = 0;
	for(int x = x_image-_neighbourhoodSize; x <= x_image+_neighbourhoodSize; x++){
		int yCtr = 0;
		for(int y = y_image-_neighbourhoodSize; y <= y_image+_neighbourhoodSize; y++){
			int zCtr = 0;
			for(int z = z_image-_neighbourhoodSize; z <= z_image+_neighbourhoodSize; z++){
					if(!( x < 0 ||   x  >= _image.GetX() ||
							y < 0 || y >= _image.GetY() ||
							z < 0 || z >= _image.GetZ() )){
						double my_a = _patchMean[atlasNo]->Get(x, y, z);
						double sigma_a = _patchStd[atlasNo]->Get(x, y, z);
						double ss = (2*my_i*my_a / (my_i*my_i + my_a*my_a)) * (2*sigma_i*sigma_a / (sigma_i*sigma_i + sigma_a*sigma_a));
						double val = _maxVal;
						if(ss>.95)
							if(IsForeground(x, y, z)){
								val = EvaluatePatch(x_image, y_image, z_image, x, y, z, atlasNo);
								tmpVals[xCtr][yCtr][zCtr] = val;
							}
						if(val > 0 && val < _wMin[atlasNo]){
							_wMin[atlasNo] = val;
						}
					}

					zCtr++;
				}
			yCtr++;
			}
		xCtr++;
	}
}

void irtkPatchBasedSegmentation::Run(){
	irtkGreyImage seg;
	if(! _patchMeasuresAvailable){
		EstablishPatchMeasures();
	}
	int nSize = _neighbourhoodSize*2+1;
	for(int x = 0; x < _image.GetX(); x++){
		for(int y = 0; y < _image.GetY(); y++){
			for(int z = 0; z < _image.GetZ(); z++){
				if(IsForeground(x, y, z)){
					for(int a = 0; a < _nAtlases; a++){
						_wMin[a] = _maxVal;
					}
					double **** tmpVals;
					tmpVals = new double***[_nAtlases];
					for(int i = 0; i < _nAtlases; i++){
						tmpVals[i] = new double**[nSize];
						for(int j = 0; j < nSize; j++){
							tmpVals[i][j] = new double*[nSize];
							for(int k = 0; k < nSize; k++){
								tmpVals[i][j][k] = new double[nSize];
								for(int l = 0; l < nSize; l++){
									tmpVals[i][j][k][l] = _maxVal;
								}
							}
						}
					}
					for(int i = 0; i < _nAtlases; i++){
						EvaluateNeighbourhood(x, y, z, i, tmpVals[i]);
					}

					double h = _maxVal;
					for(int a = 0; a < _nAtlases; a++){
						if(_wMin[a] < h)
							h = _wMin[a];
					}
					for(int a = 0; a < _nAtlases; a++){
						for(int i = 0; i < nSize; i++){
							for(int j = 0; j < nSize; j++){
								for(int k = 0; k < nSize; k++){
									double valNew = tmpVals[a][i][j][k];
									valNew = exp(-valNew/(h));
									tmpVals[a][i][j][k] = valNew;
								}
							}
						}
					}
					for(int i = 0; i < _nAtlases; i++){
						FindLabel(i, x, y, z, tmpVals[i]);
					}
					for(int i = 0; i < _nAtlases; i++){
						for(int j = 0; j < nSize; j++){
							for(int k = 0; k < nSize; k++){
								delete [] tmpVals[i][j][k];
							}
							delete [] tmpVals[i][j];
						}
						delete [] tmpVals[i];
					}
					delete [] tmpVals;
				}

			}
		}
	}
}

bool irtkPatchBasedSegmentation::IsForeground(int x, int y, int z){
	bool isForeground = false;
	if(_image.Get(x, y, z) > 0)
		isForeground = true;
	if(_useMask)
		if(_mask.Get(x, y, z) <= 0)
			isForeground = false;
	return isForeground;
}



void irtkPatchBasedSegmentation::FindLabel(int atlasNo, int x_image, int y_image, int z_image, double *** tmpVals){

	double probLabels [_nrLabels];
	for(int i = 0; i < _nrLabels; i++){
		probLabels[i] = 0;
	}
	double overallVal = 0;
	int nSize = _neighbourhoodSize*2 + 1;
	for(int i = 0; i < nSize; i++){
		for(int j = 0; j < nSize; j++){
			for(int k = 0; k < nSize; k++){
				int x_atlas = x_image - _neighbourhoodSize + i;
				int y_atlas = y_image - _neighbourhoodSize + j;
				int z_atlas = z_image - _neighbourhoodSize + k;
				int label = _padding;
				if(!( x_atlas < 0 ||   x_atlas  >= _labels[atlasNo]->GetX() ||
					  y_atlas < 0 || y_atlas >= _labels[atlasNo]->GetY() ||
					  z_atlas < 0 || z_atlas >= _labels[atlasNo]->GetZ() )){
					  label = _labels[atlasNo]->Get(x_atlas, y_atlas, z_atlas);
				}
				if(label > _padding){
					probLabels[label] += tmpVals[i][j][k];
					if(_winnerTakesAll && probLabels[label] < tmpVals[i][j][k]){
						cout << "Winner takes all" << endl;
						probLabels[label] = tmpVals[i][j][k];
					}
					overallVal += tmpVals[i][j][k];
				}
			}
		}
	}
	double maxVal = 0;
	double maxLabel = -1;
	_verbose = 0;
	for(int i = 0; i < _nrLabels; i++){
		probLabels[i] /= overallVal;
		if(_verbose)
			cout << "label " << i << ". Prob " << probLabels[i] << endl;
		double combVal = probLabels[i];
		if(combVal> maxVal){
			maxVal = combVal;
			maxLabel = i;
		}

	}

	_segmentations[atlasNo]->Put(x_image, y_image, z_image, maxLabel);
}



void irtkPatchBasedSegmentation::EstablishPatchMeasures(){
	cout << "Establish patch measures..." << endl;
	_patchMean = new irtkGreyImage*[_nAtlases+1];
	_patchStd = new irtkGreyImage*[_nAtlases+1];
	for(int i = 0; i < _nAtlases+1; i++){
		_patchMean[i] = new irtkGreyImage;
		_patchMean[i]->Initialize(_image.GetImageAttributes());
		_patchStd[i] = new irtkGreyImage;
		_patchStd[i]->Initialize(_image.GetImageAttributes());
	}
	for(int x = _patchSize; x < _image.GetX()-(_patchSize+1); x++){
		for(int y = _patchSize; y < _image.GetY()-(_patchSize+1); y++){
			for(int z = _patchSize; z < _image.GetZ()-(_patchSize+1); z++){
				if(_image.Get(x, y, z) > 0){
					double mean = 0;
					double std = 0;
					int ctr = 0;
					for(int x1 = x-_patchSize; x1<=x+_patchSize; x1++){
						for(int y1 = y-_patchSize; y1<=y+_patchSize; y1++){
							for(int z1 = z-_patchSize; z1<=z+_patchSize; z1++){

								mean += _image.Get(x1, y1, z1);
								ctr++;
							}
						}
					}
					if(mean>0)
						mean /= ctr;
					for(int x1 = x-_patchSize; x1<=x+_patchSize; x1++){
						for(int y1 = y-_patchSize; y1<=y+_patchSize; y1++){
							for(int z1 = z-_patchSize; z1<=z+_patchSize; z1++){
								std += (_image.Get(x1, y1, z1)-mean)*(_image.Get(x1, y1, z1)-mean);
							}
						}
					}
					if(std>0){
						std /= ctr;
						std = sqrt(std);
					}
					_patchMean[0]->Put(x, y, z, mean);
					_patchStd[0]->Put(x, y, z, std);
					for(int i = 1; i < _nAtlases+1; i++){
						double mean = 0;
						double std = 0;
						int ctr = 0;
						for(int x1 = x-_patchSize; x1<=x+_patchSize; x1++){
							for(int y1 = y-_patchSize; y1<=y+_patchSize; y1++){
								for(int z1 = z-_patchSize; z1<=z+_patchSize; z1++){
									mean += _atlases[i-1]->Get(x1, y1, z1);
									ctr++;
								}
							}
						}
						if(mean<0)
							mean = 0;
						if(mean>0)
							mean /= ctr;
						for(int x1 = x-_patchSize; x1<=x+_patchSize; x1++){
							for(int y1 = y-_patchSize; y1<=y+_patchSize; y1++){
								for(int z1 = z-_patchSize; z1<=z+_patchSize; z1++){
									std += (_atlases[i-1]->Get(x1, y1, z1)-mean)*(_atlases[i-1]->Get(x1, y1, z1)-mean);
								}
							}
						}
						if(std>0){
							std /= ctr;
							std = sqrt(std);
						}
						_patchMean[i]->Put(x, y, z, mean);
						_patchStd[i]->Put(x, y, z, std);
					}
				}
			}
		}
	}
	cout << "done" << endl;
	_patchMeasuresAvailable = true;
}

void irtkPatchBasedSegmentation::GetConsensusSegmentation(){

	_segmentationOutput.Initialize(_image.GetImageAttributes());
	for(int x = 0; x < _image.GetX(); x++){
		for(int y = 0; y < _image.GetY(); y++){
			for(int z = 0; z < _image.GetZ(); z++){
				if(_image.Get(x, y, z)>0){
					int labels[_nrLabels];
					for(int i = 0; i < _nrLabels; i++){
						labels[i] = 0;
					}

					for(int i = 0; i < _nAtlases; i++){

						labels[_segmentations[i]->Get(x, y, z)] += 1;
					}
//
					int maxLabel = -1;
					int max = 0;
					for(int i = 0; i < _nrLabels; i++){
						if(labels[i]>max){
							maxLabel = i;
							max = labels[i];
						}
					}
					_segmentationOutput.Put(x, y, z, 0, maxLabel);
				}
			}
		}
	}
}

void irtkPatchBasedSegmentation::GetProbabilisticSegmentation(){
	cout << "PatchBasedSegmentation::GetProbabilisticSegmentation() -- assumes binary segmentation!" << endl;
	_segmentationOutput.Initialize(_image.GetImageAttributes());
	for(int x = 0; x < _image.GetX(); x++){
		for(int y = 0; y < _image.GetY(); y++){
			for(int z = 0; z < _image.GetZ(); z++){
				if(_image.Get(x, y, z)>0){
					int count = 0;
					for(int i = 0; i < _nAtlases; i++){

						if(_segmentations[i]->Get(x, y, z) > 0)
							count += 100;
					}
					count /= _nAtlases;

					_segmentationOutput.Put(x, y, z, 0, count);
				}
			}
		}
	}
}

void irtkPatchBasedSegmentation::WriteSegmentation(char * fileName){
		_segmentationOutput.Write(fileName);
}

void irtkPatchBasedSegmentation::SetPatchMeasures(irtkGreyImage ** patchMean, irtkGreyImage ** patchStd){
	_patchMean = patchMean;
	_patchStd = patchStd;
	_patchMeasuresAvailable = true;
}


