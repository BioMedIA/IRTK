#ifndef SPECTRALCLUSTERING_H_
#define SPECTRALCLUSTERING_H_

#include <irtkEigenAnalysis.h>
#include <irtkManifoldLearning.h>
#include <string>

class irtkSpectralClustering : public irtkManifoldLearning
{
	protected:
	
	private:

	
	public:
		irtkSpectralClustering(int nrSubjects, int nrFeatures);
		void DoSpectralEmbedding();

};

#endif /*SPECTRALCLUSTERING_H_*/
