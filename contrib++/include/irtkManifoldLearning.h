#ifndef IRTKMANIFOLDLEARNING_H_
#define IRTKMANIFOLDLEARNING_H_

#include <irtkEigenAnalysis.h>
#include <string>

class irtkManifoldLearning
{
	protected:
	
	
		irtkMatrix _w;
		int _nrFeatures;
		int _nrSubjects;
		double ** _features;
		double ** _distances;
		void EstablishDistances();
		irtkManifoldLearning();
	
	
	
	public:
		void Initialize(string csvFilename);
		void Initialize(double ** input);
		void Initialize(irtkMatrix input);
		double ** GetEmbedding();
		void WriteFeatures(string);
		double GetDistance(int, int);
		void GetNeighbours(int node, int nrNeighbours, int * neighbourList);
};

#endif /*IRTKMANIFOLDLEARNING_H_*/
