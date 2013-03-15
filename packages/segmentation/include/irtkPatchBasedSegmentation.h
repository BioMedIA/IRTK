/*
 * irtkBatchBasedSegmentation.cc
 *
 *  Created on: Jun 21, 2011
 *      Author: rw1008
 */

#ifndef IRTKPATCHBASEDSEGMENTAION_H_
#define IRTKPATCHBASEDSEGMENTAION_H_

#include <irtkImage.h>
#include <irtkRegistration.h>
#include <irtkGaussian.h>


class irtkPatchBasedSegmentation{

private:

	irtkRealImage _image, _hdimage, **_atlases;
	irtkGreyImage **_labels;
	irtkGreyImage **_segmentations;
	irtkGreyImage ** _patchMean, ** _patchStd;
	irtkGreyImage  _segmentationOutput;

	double * priorMean, * priorStd;
	double _maskMean, _maskStd;
	double *_atlasMean, *_atlasStd;

	int _patchSize, _neighbourhoodSize;
	int _simMeasure;
	int _maxVal, _padding;
	//_patchSize: nr of voxels around voxel of interest
	double _sigma, _beta, _epsilon, * _wMin;
	//int _x_image, _y_image, _z_image;
	double EvaluatePatch(int x_image, int y_image, int z_image, int x_atlas, int y_atlas, int z_atlas, int atlasNo);
	double *** _allSimilarities;

	int * _availableLabels;

	int _verbose;

	int _nAtlases;

	int _minCtr;
	bool _useMask, _patchMeasuresAvailable;
	bool _winnerTakesAll;
	irtkGreyImage _mask;
	void FindImage(int, int, int, double ****);


public:
	irtkPatchBasedSegmentation(irtkRealImage, irtkRealImage **, irtkGreyImage **, int, int, int);
	~irtkPatchBasedSegmentation();
	bool IsForeground(int, int, int);
	void EvaluateNeighbourhood(int, int, int, int, double ***);
	int GetLabel(int, int, int, int);
	void FindLabel(int, int, int, int, double ***);
	void Run();
	void WriteSegmentation(char *);
	void EstablishPatchMeasures();
	int GetNeighbourhoodSize();
	void AddMinVal(double);
	double EvalMinVal();
	irtkRealImage _hImage;
	int _nrLabels;
	irtkGreyImage GetImage();
	void Output(int);
	void SetMask(irtkGreyImage mask);
	irtkGreyImage GetSegmentation();

	void GetConsensusSegmentation();
	irtkRealImage GetConsensusImage();
	void GetProbabilisticSegmentation();
	void GetCardiacSegmentation(double ed, double es);
	void SetPatchMeasures(irtkGreyImage ** patchMean, irtkGreyImage ** patchStd);
	void SetPatchSizes(int patchSize, int neighborhoodSize);
	void SetPadding(int padding);
	void WinnerTakesAll();

};

inline void irtkPatchBasedSegmentation::WinnerTakesAll(){
	 _winnerTakesAll = true;
}

inline void irtkPatchBasedSegmentation::SetPadding(int padding){
	_padding = padding;
}

inline void irtkPatchBasedSegmentation::SetPatchSizes(int patchSize, int neighborhoodSize){
	_patchSize = patchSize;
	_neighbourhoodSize = neighborhoodSize;
}



inline irtkGreyImage irtkPatchBasedSegmentation::GetSegmentation(){

		return _segmentationOutput;
}

inline void irtkPatchBasedSegmentation::SetMask(irtkGreyImage mask){
	_useMask = true;
	_mask = mask;
	_maskMean = 0;
	_maskStd = 0;
	_atlasMean = new double[_nAtlases];
	_atlasStd = new double[_nAtlases];
	irtkGreyPixel * mPtr = _mask.GetPointerToVoxels();
	irtkRealPixel * iPtr = _image.GetPointerToVoxels();
	irtkRealPixel ** aPtrs = new irtkRealPixel*[_nAtlases];
	for(int i = 0; i < _nAtlases; i++){
		aPtrs[i] = _atlases[i]->GetPointerToVoxels();
		_atlasMean[i] = 0;
		_atlasStd[i] = 0;
	}
	int ctr = 0;
	for(int i = 0; i < _image.GetNumberOfVoxels(); i++){
		if(*mPtr>0){
			ctr++;
			_maskMean += *iPtr;

			for(int j = 0; j < _nAtlases; j++){

				_atlasMean[j] += *aPtrs[j];

			}

		}
		for(int j = 0; j < _nAtlases; j++){
			aPtrs[j]++;
		}
		iPtr++;
		mPtr++;
	}
	_maskMean /= ctr;

	mPtr = _mask.GetPointerToVoxels();
	iPtr = _image.GetPointerToVoxels();
	for(int i = 0; i < _nAtlases; i++){
		aPtrs[i] = _atlases[i]->GetPointerToVoxels();
		_atlasMean[i] /= ctr;
	}

	for(int i = 0; i < _image.GetNumberOfVoxels(); i++){
		if(*mPtr>0){
			_maskStd += (*iPtr-_maskMean)*(*iPtr-_maskMean);
			for(int j = 0; j < _nAtlases; j++){
				_atlasStd[j] += (*aPtrs[j]-_atlasMean[j])*(*aPtrs[j]-_atlasMean[j]);

			}

		}
		for(int j = 0; j < _nAtlases; j++){
			aPtrs[j]++;
		}
		iPtr++;
		mPtr++;
	}

	_maskStd /= (ctr-1);
	_maskStd = sqrt(_maskStd);
	for(int i = 0; i < _nAtlases; i++){
		_atlasStd[i] /= (ctr-1);
		_atlasStd[i] = sqrt(_atlasStd[i]);
	}

}



inline irtkGreyImage irtkPatchBasedSegmentation::GetImage(){
	return _image;
}



inline int irtkPatchBasedSegmentation::GetLabel(int atlasN, int x, int y, int z){

	int label = _labels[atlasN]->Get(x, y, z);

	return label;
}


#endif /*IRTKPATCHBASEDSEGMENTAION_H_*/
