#ifndef IRTKGRAPHCUTSEGMENTATION_4D_H_
#define IRTKGRAPHCUTSEGMENTATION_4D_H_

#include <irtkImage.h>
#include <irtkProbabilisticAtlas.h>
#include <irtkGaussian.h>
#include <irtkImage.h>
#include <graph.h>
#include <string>

#ifndef MOG_H_
#define MOG_H_

class MOG{

	private:
		int numTissues;
		int numStructures;
		irtkRealImage **tissueProb;
		irtkRealImage **structureProb;
		irtkGaussian **tissueGaussian;
		irtkGaussian **structureGaussian;

	public:
		MOG(int, int);
		MOG();
		void setTissueProb(int, irtkRealImage*, double, double);
		void setStructureProb(int, irtkRealImage*, double, double);
		double evaluate(double, int, int, int);
		double evaluate(double, int, int, int, int);

};

#endif /*MOG_H_*/



class irtkGraphCutSegmentation_4D : public irtkObject
{
 

protected:
	irtkGenericImage<irtkRealPixel> _input;
	irtkRealImage _foregroundAtlas;
	irtkGaussian* _gaussians;
	irtkGaussian* _tissueProb;
	irtkGenericImage<irtkGreyPixel> *_segmentation;
	irtkRealImage segmentation;
	irtkGreyImage _mask;
	int _numTissues;
	irtkRealPixel _padding;
	bool useMask;

private:
	int _xDim;
	int _yDim;
	int _zDim;
	int numNodes;
	int mi1;
	int mi2;
	int mi3;
	int _wrongClass;
	double _lambda;
	double _gamma;
	double _sigma;
	double _sigmaG;
	double _c;
	double _mul;
	Graph<double, double, double> *_graph;
	irtkGenericImage<irtkRealPixel> *_nodes;
	irtkGenericImage<irtkRealPixel> _gradMagnitude;
	irtkRealImage **_tissuePriors;
	void GenerateGraph(int, int);
	int UpdateSegmentation(int, int);
	MOG mog;

public:
	double GetNeighbourWeight(irtkRealPixel *, irtkRealPixel *, irtkRealPixel *, double, double);
	double GetTerminalWeight(irtkRealPixel *, int, int, int, int, int);
	double GetTerminalWeight(irtkRealPixel *, int, int, int, int);
	void SetMog(MOG);
	void SetTissuePriors(irtkRealImage **);
	irtkRealImage GetSegmentation();
	void setParameters(double, double, double, double);
	irtkGraphCutSegmentation_4D(int, irtkRealImage);
	void SetMi(int, int, int);
	void Iterate();
	void SetInput(const irtkRealImage &, irtkGenericImage<float> );
	virtual ~irtkGraphCutSegmentation_4D();
	void GenerateGaussians();
	void GenerateTissueGaussians(double, double, double, double, double, double);
	void GenerateGaussians(int, double, double);
	int GetWrongClassified();
	void setMask(irtkGreyImage);
};

class GetGaussian : public irtkObject
{
public:
static double getGaussian(double);

};

#endif /*IRTKGraphCutSegmentation_4D_H_*/



