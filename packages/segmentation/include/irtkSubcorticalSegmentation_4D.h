#ifndef IRTKSUBCORTICALSEGMENTATION_4D_H_
#define IRTKSUBCORTICALSEGMENTATION_4D_H_

#include <string>

class irtkSubcorticalSegmentation_4D{

	private:
		string imageName;
		string imageDir;
		string atlasDir;
		string weightImageDir;
		string outputDir;
		string tissuePriorDir;
		irtkGreyImage mask;
		irtkRealImage **tissuePriors;
		int numTissues;
		int structure;
		irtkRealImage **atlasI;
		string **names;
		int **structValues;
		irtkRealImage input;
		bool useMask;
		double _my[83];
		double _sigma[83];
		double tissueMy[3];
		double tissueS[3];
		double counter[3];
		irtkProbabilisticAtlas atlas;
		void readStructureNames(char *);

	public:
		irtkSubcorticalSegmentation_4D();
		irtkSubcorticalSegmentation_4D(string, char *);
		void init();
		void init(irtkRealImage , irtkRealImage **, irtkRealImage **, int, string, string **, int **, string);
		void readParameters(char *);
		
		void generateGaussians(double);
		void doSegmentation(double, double, double, double);
		irtkRealImage getBackground(irtkRealImage, irtkRealImage **, int);
		bool Read(const char *, const char*);
		bool Read(string line);
		void setMask(const char *);
};



#endif /*IRTKSUBCORTICALSEGMENTATION_H_*/
