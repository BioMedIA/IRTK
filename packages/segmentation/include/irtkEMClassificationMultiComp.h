/*=========================================================================

Library   : Image Registration Toolkit (IRTK)
Module    : $Id$
Copyright : Imperial College, Department of Computing
Visual Information Processing (VIP), 2008 onwards
Date      : $Date$
Version   : $Revision$
Changes   : $Author$

=========================================================================*/

#ifndef _IRTKEMCLASSIFICATIONMULTICOMP_H

#define _IRTKEMCLASSIFICATIONMULTICOMP_H

/*

EM maximisation algorithm with Multiple Component estimation

*/

class irtkEMClassificationMultiComp : public irtkEMClassification
{
protected:
	/// Number of components of every atlas
	int *_number_of_components;
	int _number_of_atlas;
	int *_ns;
	int *_ne;

public:

	/// Constructor
	irtkEMClassificationMultiComp(int noAtlas, irtkRealImage **atlas, int *noComp);

	/// Constructor
	irtkEMClassificationMultiComp();

	/// Constructor
	~irtkEMClassificationMultiComp();

	/// Constructor
	irtkEMClassificationMultiComp(int noAtlas, irtkRealImage **atlas, irtkRealImage *background, int *noComp);

	/// Initialize multiple components filter
	virtual void Initialise();

	/// Estimates posterior probabilities
	virtual void EStep();

	/// Estimates parameters
	virtual void MStep();

	/// Susbended
	/// Estimates posterior probabilities
	virtual void InitialEStep(int j, int n);

	/// Susbended
	/// Estimates parameters
	virtual void InitialMStep(int j, int n);

	/// Computes log likelihood for current parameters
	virtual double LogLikelihood();

	/// Execute one iteration and return log likelihood
	virtual double Iterate(int iteration);

	/// Construct segmentation based on current posterior probabilities
	virtual void ConstructSegmentation(irtkRealImage &);
	/// Construct segmentation based on current posterior probabilities
	virtual void ConstructSegmentation();
	/// Construct segmentation based on current posterior probabilities
	virtual void ConstructSegmentationWithPadding(irtkRealImage &);
	/// Get ProbMap
	virtual void GetProbMap(int i,irtkRealImage& image);

};

#endif
