#include <irtkImage.h>
#include <irtkRegistration.h>
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
#include <irtkPairwiseSimilarity.h>

#ifdef HAS_TBB
#include "tbb/task_scheduler_init.h"
#include "tbb/parallel_for.h"
#include "tbb/blocked_range.h"

using namespace tbb;

class MultiThreadedSimilarity {

	bool _useMasks;
	irtkGreyImage ** _images, ** _regions;
	int i, _nrRegions, _similarityType, _nrRows;
	double *** _results;
	bool completeImage;
	vector <string> _files;
	bool _transformImages, _twoSets, _singleRegionMask;
	int _padding;

public:
	MultiThreadedSimilarity(bool useMasks , irtkGreyImage ** images, irtkGreyImage ** masks, int _i, double *** results, int nrRegions, int simType, bool twoSets, int nrRows, bool singleRegionMask, int padding) {

		_useMasks = useMasks;
		_images = images;
		_regions = masks;
		i = _i;
		_nrRegions = nrRegions;
		completeImage = false;
		_results = results;
		_transformImages = false;
		_similarityType = simType;
		_twoSets = twoSets;
		_nrRows = nrRows;
		_singleRegionMask = singleRegionMask;
		_padding = padding;
	}
	void operator()(const blocked_range<int> &r) const {
		for (int j = r.begin(); j!= r.end(); j++) {
			if(j > i || _twoSets) {

				irtkSimilarityMetric **metrics;
				metrics = new irtkSimilarityMetric*[_nrRegions];
				irtkGreyImage *i1, *i2;
				i1 = _images[i];
				if(_transformImages) {

				}
				else {
					i2 = _images[j];
				}
				cout << i << " " << j << endl;

				irtkGreyPixel * i1Ptr = i1->GetPointerToVoxels();
				irtkGreyPixel * i2Ptr = i2->GetPointerToVoxels();
				int numberOfBins = 256;

				for(int k = 0; k < _nrRegions; k++) {

					switch(_similarityType) {
//					if(_similarityType == 1){
						case 1 :
						irtkGreyPixel i1_min, i1_max, i1_nbinss;
						irtkGreyPixel i2_min, i2_max, i2_nbinss;

						i1->GetMinMax(&i1_min, &i1_max);
						i2->GetMinMax(&i2_min, &i2_max);

						i2_nbinss = irtkCalculateNumberOfBins(i2, numberOfBins,
								i2_min, i2_max);
						i1_nbinss = irtkCalculateNumberOfBins(i1, numberOfBins,
								i1_min, i1_max);
						metrics[k] = new irtkNormalisedMutualInformationSimilarityMetric(i1_nbinss, i2_nbinss);
//					}
						break;
//					else{
						case 2 :
						metrics[k] = new irtkSSDSimilarityMetric();
						break;
						case 3 :
						metrics[k] = new irtkCrossCorrelationSimilarityMetric();
						break;
						default:
						metrics[k] = new irtkSSDSimilarityMetric();
					}
					metrics[k]->Reset();
				}

				cout << "region: finished. nrRegions: " << _nrRegions << endl;
				irtkGreyPixel **regionPtr;
				int nrRegionPointers = _nrRegions;
				if(_singleRegionMask) {
					nrRegionPointers = 1;
				}
				if(_useMasks) {
					regionPtr = new irtkGreyPixel*[nrRegionPointers];
					for(int k = 0; k < nrRegionPointers; k++) {
						regionPtr[k] = new irtkGreyPixel;
						regionPtr[k] = _regions[k]->GetPointerToVoxels();
						//		cout << "region: " << k << endl;
					}
				}
				else {
					regionPtr = new irtkGreyPixel*[1];
					regionPtr[0] = _regions[0]->GetPointerToVoxels();
				}

				int valCtr = 0;
				for(int ii = 0; ii < i1->GetNumberOfVoxels(); ii++) {
					int val1 = *i1Ptr;
					int val2 = *i2Ptr;
					if(val1 > _padding && val2 > _padding) {
						for(int r = 0; r < _nrRegions; r++) {
							bool voxelInRegion = false;
							if(!_singleRegionMask) {
								if(*regionPtr[r]==0) {
									voxelInRegion = true;
								}
							}
							else {
								if(*regionPtr[0]==r) {
									voxelInRegion = true;
								}
							}
							if(voxelInRegion) {
								//								if(valCtr == 11)
								//if(valCtr <= 20)
								//cout << ii << " " << val1 << ", " << val2 << endl;
								//									cout << *regionPtr[0] << endl;
								metrics[r]->Add(val1, val2);
								valCtr += 1;
							}

						}
					}
					for(int r = 0; r < nrRegionPointers; r++) {
						regionPtr[r]++;
					}
					i1Ptr++;
					i2Ptr++;
				}
				cout << "Vals: " << valCtr << endl;
				int noValue = 0;
				for(int r = 0; r < _nrRegions; r++) {
					double results = -1;

					results = metrics[r]->Evaluate();
					//cout << results << " " << i << " " << j << endl;
					if(_twoSets) {
						_results[r][i][j-_nrRows] = results;
					}
					else {
						_results[r][i][j] = results;
						_results[r][j][i] = results;
					}
				}

				for(int r = 0; r < _nrRegions; r++) {
					delete metrics[r];
				}
				delete [] metrics;
				delete [] regionPtr;

			}
		}
	}
};

#endif

irtkPairwiseSimilarity::irtkPairwiseSimilarity()
{
	_nrRegions = -1;
	_transformImages = false;
	_padding = 0;
}
;
void irtkPairwiseSimilarity::Initialize(int nrSubjects, int nrSubjects2, bool useMasks,
		irtkGreyImage ** regions, int simType, int nrRegions)
{

	_regions = regions;
	_nrSubjects = nrSubjects;
	_images = new irtkGreyImage*[_nrSubjects + nrSubjects2];
	_nrRows = _nrSubjects;
	if (nrSubjects2 > 0) {
		_twoSets = true;
		_nrCols = nrSubjects2;
	}
	else {
		_twoSets = false;
		_nrCols = _nrSubjects;
	}
//	short min, max;
//	_regions[0]->GetMinMax(&min, &max);
//	_nrRegions = max+1;
	_nrRegions = nrRegions;
	_singleRegionMask = false;
	//check if single mask with several non-overlapping regions
	if (_nrRegions == 1) {
		irtkGreyPixel maskMin, maskMax;
		_regions[0]->GetMinMax(&maskMin, &maskMax);
		if (maskMax > 0) {
			_singleRegionMask = true;
			_nrRegions = maskMax + 1;
		}
	}
	cout << _nrRegions << endl;
	_useMasks = useMasks;
	_results = new double**[_nrRegions];
	for (int i = 0; i < _nrRegions; i++) {
		_results[i] = new double*[_nrRows];
		for (int j = 0; j < _nrRows; j++) {
			_results[i][j] = new double[_nrCols];
			for (int k = 0; k < _nrCols; k++) {
				_results[i][j][k] = 0;
			}
		}
	}
	_similarityType = simType;

}

void irtkPairwiseSimilarity::Initialize(int nrSubjects, bool useMasks)
{

	_nrSubjects = nrSubjects;
	_nrRegions = 1;
	_useMasks = useMasks;
	SetupRegions();
	_results = new double**[_nrRegions];
	for (int i = 0; i < _nrRegions; i++) {
		_results[i] = new double*[_nrSubjects];
		for (int j = 0; j < _nrSubjects; j++) {
			_results[i][j] = new double[_nrSubjects];
			for (int k = 0; k < _nrSubjects; k++) {
				_results[i][j][k] = 0;
			}
		}
	}

}

void irtkPairwiseSimilarity::LoadImages(irtkGreyImage ** images)
{
	_images = images;
	if (_nrRegions == -1)
		SetupRegions();
}

void irtkPairwiseSimilarity::LoadImages(string imagefilename, string imageDir)
{
	vector < string > files = vector<string>();
	string line;
	ifstream from(imagefilename.c_str());
	if (!from) {
		cerr << "Can't open file "
				<< endl;
		exit(1);
	}

	if (from.is_open()) {
		while (!from.eof()) {
			getline(from, line);
			string fileName = line;
			files.push_back(fileName);
		}
		from.close();
	}

	if ((int)files.size() != _nrSubjects && !(_twoSets && _nrRows + _nrCols == (int)files.size())) {
		cerr << "Nr of files does not agree with nrSubjects provided" << endl;
	}
	int nrImages = _nrRows;
	if (_twoSets)
		nrImages += _nrCols;
	for (int i = 0; i < nrImages; i++) {
		cout << i << " " << files[i] << endl;
		string name = imageDir + files[i];
		_images[i] = new irtkGreyImage();
		_images[i]->Read(name.c_str());
	}
	if (_nrRegions == -1)
		SetupRegions();
}

void irtkPairwiseSimilarity::SetTemplate(string filename)
{
	_template = new irtkGreyImage;
	_template->Read(filename.c_str());
}

void irtkPairwiseSimilarity::SetupRegions()
{
//	_useMasks = false;
	if (_useMasks) {
		_regions = new irtkGreyImage*[_nrRegions];
		for (int r = 0; r < _nrRegions; r++) {
			_regions[r] = new irtkGreyImage;
			_regions[r]->Initialize(_images[0]->GetImageAttributes());
			for (int x = 0; x < _regions[0]->GetX(); x++) {
				for (int y = 0; y < _regions[0]->GetY(); y++) {
					for (int z = 0; z < _regions[0]->GetZ(); z++) {
						_regions[r]->Put(x, y, z, 0, 1);
					}
				}
			}

		}
	}
	else {
		_regions = new irtkGreyImage*[1];
		_regions[0] = new irtkGreyImage;
		_regions[0]->Initialize(_template->GetImageAttributes());
		double spacing = 20;
		int ctr = 1;
		for (int x = spacing / 2; x < _template->GetX() - spacing / 2; x += spacing) {
			for (int y = spacing / 2; y < _template->GetY() - spacing / 2; y += spacing) {
				for (int z = spacing / 2; z < _template->GetZ() - spacing / 2; z += spacing) {
					int inside = 0;
					int insideCtr = 0;
					for (int x1 = x - spacing / 2; x1 <= x + spacing / 2; x1++) {
						for (int y1 = y - spacing / 2; y1 <= y + spacing / 2; y1++) {
							for (int z1 = z - spacing / 2; z1 <= z + spacing / 2; z1++) {
								insideCtr++;
								if (x1 < _template->GetX() && y1 < _template->GetY() && z1 < _template->GetZ())
									if (_template->Get(x1, y1, z1) > 0) {
										inside++;
									}

							}
						}
					}
					if (inside > spacing * spacing * spacing * .75) {
						for (int x1 = x - spacing / 2; x1 <= x + spacing / 2; x1++) {
							for (int y1 = y - spacing / 2; y1 <= y + spacing / 2; y1++) {
								for (int z1 = z - spacing / 2; z1 <= z + spacing / 2; z1++) {
									if (x1 < _template->GetX() && y1 < _template->GetY() && z1 < _template->GetZ())
										_regions[0]->Put(x1, y1, z1, ctr);
								}
							}
						}
						ctr++;
					}

				}
			}
		}
		_nrRegions = ctr;
		cout << "ctr " << ctr << endl;
		_regions[0]->Write("/vol/vipdata/users/rw1008/ADNI/test_regions_20mm.nii.gz");

	}
}

void irtkPairwiseSimilarity::WriteSimilarities(string filename)
{
	cout << "Writing " << filename << endl;
	for (int r = 0; r < _nrRegions; r++) {
		ofstream output;
		ostringstream stm;
		stm << r;
		//	string filename = dir + "/" + namePrefix + "_" + stm.str() + ".csv";
		output.open(filename.c_str());
		if (_twoSets)
			for (int i = 0; i < _nrRows; i++) {
				for (int j = 0; j < _nrCols; j++) {
					output << _results[r][i][j];
					if (j == _nrCols - 1) {
						output << endl;
					}
					else {
						output << ",";
					}
				}
			}
		else
			for (int i = 0; i < _nrSubjects; i++) {
				for (int j = 0; j < _nrSubjects; j++) {
					output << _results[r][i][j];
					if (j == _nrSubjects - 1) {
						output << endl;
					}
					else {
						output << ",";
					}
				}
			}
		output.close();
	}
}

void irtkPairwiseSimilarity::WriteSimilarities(string dir, string namePrefix)
{

	for (int r = 0; r < _nrRegions; r++) {
		ofstream output;
		ostringstream stm;
		stm << r;
		string filename = dir + "/" + namePrefix + "_" + stm.str() + ".csv";
		cout << "Write " << filename << endl;
		output.open(filename.c_str());
		if (_twoSets)
			for (int i = 0; i < _nrRows; i++) {
				for (int j = 0; j < _nrCols; j++) {
					output << _results[r][i][j];
					if (j == _nrCols - 1) {
						output << endl;
					}
					else {
						output << ",";
					}
				}
			}
		else
			for (int i = 0; i < _nrSubjects; i++) {
				for (int j = 0; j < _nrSubjects; j++) {
					output << _results[r][i][j];
					if (j == _nrSubjects - 1) {
						output << endl;
					}
					else {
						output << ",";
					}
				}
			}
		output.close();
	}
}

void irtkPairwiseSimilarity::GetSimilarities()
{
	for (int i = 0; i < _nrRows; i++) {
		cout << i << ": ";
#ifdef HAS_TBB
//			cout << "Nr of Regions: " << _nrRegions << endl;
		task_scheduler_init init;
		MultiThreadedSimilarity evaluate(_useMasks , _images, _regions,i, _results, _nrRegions,_similarityType, _twoSets, _nrRows,_singleRegionMask,_padding);
		int blocks = -1;
		if(i==0)
		blocks = _nrSubjects;
		else
		blocks = 1;
//			blocks = 10;
//			blocks = 8;
//			blocks = _nrSubjects;
		if(_twoSets) {
			if(i==0)
			blocks = _nrRows+_nrCols;
			else
			blocks = 24;

			parallel_for(blocked_range<int>(_nrRows, _nrRows+_nrCols, int(blocks)), evaluate);
		}
		else
		parallel_for(blocked_range<int>(1, _nrRows, int(blocks)), evaluate);
#else
		for (int j = 0; j < _nrRows + _nrCols; j++) {
			if (j > i || _twoSets) {
				//	cout << _twoSets << " " << j << " " << _nrRows << " " << _nrSubjects << endl;
				if (!_twoSets || (j > (_nrRows - 1))) {
					irtkSimilarityMetric **metrics;
					metrics = new irtkSimilarityMetric*[_nrRegions];
					irtkGreyImage *i1, *i2;
					i1 = _images[i];
					if (_transformImages) {

					}
					else {
						i2 = _images[j];
					}
					//cout << i << " " << j << endl;

					irtkGreyPixel * i1Ptr = i1->GetPointerToVoxels();
					irtkGreyPixel * i2Ptr = i2->GetPointerToVoxels();
					int numberOfBins = 256;

					for (int k = 0; k < _nrRegions; k++) {

						switch (_similarityType)
						{
						//					if(_similarityType == 1){
						case 1:
							irtkGreyPixel i1_min, i1_max, i1_nbinss;
							irtkGreyPixel i2_min, i2_max, i2_nbinss;

							i1->GetMinMax(&i1_min, &i1_max);
							i2->GetMinMax(&i2_min, &i2_max);

							i2_nbinss = irtkCalculateNumberOfBins(i2, numberOfBins,
									i2_min, i2_max);
							i1_nbinss = irtkCalculateNumberOfBins(i1, numberOfBins,
									i1_min, i1_max);
							metrics[k] = new irtkNormalisedMutualInformationSimilarityMetric(i1_nbinss,
									i2_nbinss);
							//					}
							break;
							//					else{
						case 2:
							metrics[k] = new irtkSSDSimilarityMetric();
							break;
						case 3:
							metrics[k] = new irtkCrossCorrelationSimilarityMetric();
							break;
						default:
							metrics[k] = new irtkSSDSimilarityMetric();
						}
						metrics[k]->Reset();
					}

//						cout << "region: finished. nrRegions: " << _nrRegions << endl;
					irtkGreyPixel **regionPtr;
					int nrRegionPointers = _nrRegions;
					if (_singleRegionMask) {
						nrRegionPointers = 1;
					}
					if (_useMasks) {
						regionPtr = new irtkGreyPixel*[nrRegionPointers];
						for (int k = 0; k < nrRegionPointers; k++) {
							regionPtr[k] = new irtkGreyPixel;
							regionPtr[k] = _regions[k]->GetPointerToVoxels();
							//		cout << "region: " << k << endl;
						}
					}
					else {
						regionPtr = new irtkGreyPixel*[1];
						regionPtr[0] = _regions[0]->GetPointerToVoxels();
					}
//						cout << "a" << endl;
					int valCtr = 0;
					for (int ii = 0; ii < i1->GetNumberOfVoxels(); ii++) {
						int val1 = *i1Ptr;
						int val2 = *i2Ptr;
						if (val1 > 0 && val2 > 0) {
							for (int r = 0; r < _nrRegions; r++) {
								bool voxelInRegion = false;
								if (!_singleRegionMask) {
									if (*regionPtr[r] == 0) {
										voxelInRegion = true;
									}
								}
								else {
									if (*regionPtr[0] == r) {
										voxelInRegion = true;
									}
								}
								if (voxelInRegion) {
									//								if(valCtr == 11)
									//if(valCtr <= 20)
									//cout << ii << " " << val1 << ", " << val2 << endl;
									//									cout << *regionPtr[0] << endl;
									metrics[r]->Add(val1, val2);
									valCtr += 1;
								}

							}
						}
						for (int r = 0; r < nrRegionPointers; r++) {
							regionPtr[r]++;
						}
						i1Ptr++;
						i2Ptr++;
					}
//						cout << "Vals: " << valCtr << endl;
					for (int r = 0; r < _nrRegions; r++) {
						double results = -1;
//							cout << "a" << endl;
						results = metrics[r]->Evaluate();
//							cout << "res: " << results << endl;

						cout << results << endl;
						if (_twoSets) {
//								cout << r << " " << i << " " << j << " " << _nrRows << endl;
							_results[r][i][j - _nrRows] = results;
//								cout << "b" << endl;
						}
						else {
							_results[r][i][j] = results;
							_results[r][j][i] = results;
						}
					}

					for (int r = 0; r < _nrRegions; r++) {
						delete metrics[r];
					}
					delete[] metrics;
					delete[] regionPtr;

				}
			}
		}

#endif
	}

}

void irtkPairwiseSimilarity::SetPadding(int padding)
{
	_padding = padding;
}

double irtkPairwiseSimilarity::GetSimilarity(int region, int row, int col)
{
	return _results[region][row][col];
}

