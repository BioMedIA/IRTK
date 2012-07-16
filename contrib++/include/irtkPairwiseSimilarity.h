#include <irtkImage.h>

#ifndef IRTKPAIRWISESIMILARITY_H_
#define IRTKPAIRWISESIMILARITY_H_

#define NMI 1;
#define SSD 2;


class irtkPairwiseSimilarity
{
	protected:
	
	private:
		irtkGreyImage ** _images;
		irtkGreyImage ** _regions;
		irtkGreyImage * _template;
		double *** _results;
//		double ** _regions;
		int _nrSubjects;
		int _nrCols;
		int _nrRows;
		int _nrRegions;
		bool _useMasks;
		bool _twoSets;
		bool _transformImages;
		int _similarityType;
		bool _singleRegionMask;
		int _padding;

		
		
		
	public:
	void SetupRegions();
		irtkPairwiseSimilarity();
		void Initialize(int, bool	);
		void Initialize(int, int, bool, irtkGreyImage **	,int, int);
		void LoadImages(string imagefilename, string);
		void LoadImages(irtkGreyImage **);
		void GetSimilarities();
		int IsRegion(int, int, int);
		void WriteSimilarities(string dir, string namePrefix);
		void WriteSimilarities(string);
		void SetTemplate(string filename);
		double GetSimilarity(int region, int row, int col);
		void SetPadding(int padding);
		
		
		
	
		
		
};



#endif /*PAIRWISESIMILARITY_H_*/
