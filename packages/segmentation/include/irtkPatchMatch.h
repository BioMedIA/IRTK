/*=========================================================================

Library   : Image Registration Toolkit (IRTK)
Module    : $Id$
Copyright : Imperial College, Department of Computing
Visual Information Processing (VIP), 2008 onwards
Date      : $Date$
Version   : $Revision$
Changes   : $Author$

=========================================================================*/

#ifndef IRTKPATCHMATCH_H_
#define IRTKPATCHMATCH_H_

struct NearstNeighbor{
	int x;
	int y;
	int z;
	int n;
	double weight;
};

class irtkPatchMatch{

protected:

	/// target image
	irtkGreyImage *target;
	/// target gradient
	irtkGreyImage *targetgradient;
	/// decimated image
	irtkGreyImage *decimated;
	/// decimated gradient
	irtkGreyImage *decimatedgradient;
	/// label image
	irtkGreyImage **labels;
	/// source image
	irtkGreyImage **sources;
	/// source gradient
	irtkGreyImage **sourcesgradient;
	/// source gradient
	irtkGreyImage **search;
	/// blured source image
	irtkGreyImage **blured;
	/// blured source gradient
	irtkGreyImage **bluredgradient;
	/// is decimated
	bool isdecimated;
	/// NNF
	NearstNeighbor **nnfs;
	/// vote
	double *votevalues;
	/// vote
	double *voteweights;
	/// vote label
	double *votelabels;
	/// Random rate
	double randomrate;
	/// Min Max labels
	short minlabel;
	short maxlabel;
	int maxdistance;
	/// Debug
	int debug;
	/// Radius
	int radius;
	/// Similarity levels
	int slevels;
	/// Number of source images
	int nimages;
	/// Number of neighborhoods
	int nneighbour;
	/// initialize field's weight
	virtual double initialize();
	/// initial guess of the mapping
	virtual void initialguess();
	/// find minum flow
	double minimizeflow();
	/// calculate distance between patches
	virtual double distance(int x1, int y1, int z1, int x2, int y2, int z2, int n);
	/// random search the space to find a better link
	int randomlink(int i, int j, int k, int n, int index = -1);
	/// createsearchimage
	void createsearchimages();
	/// test if it is search
	int checkissearch(int i, int j, int k, int n);
	/// propergate from one to another using the offset
	int propergate(int x, int y, int z, int i, int j, int k, int offsetx, int offsety, int offestz, int index1 = -1, int index2 = -1);
	/// vote from one patch to another
	void votepatch(int i, int j, int k, int x, int y, int z, int n, double weight);
	/// vote label from one patch to another
	void votelabel(int i, int j, int k, int x, int y, int z, int n, double weight);
	/// expectation maximization
	virtual void EMstep();
public:
	/// constructor
	irtkPatchMatch(irtkGreyImage *target, irtkGreyImage **source, int radius = 2, int nimages = 1, int nneighbour = 1, int slevels = 1);
	/// destructor
	virtual ~irtkPatchMatch();
	/// get nnf
	NearstNeighbor** getNearstNeighbor();
	/// upsample NNF to refrence image
	void upSampleNNF(irtkPatchMatch *reference);
	/// calculate NNF
	void runEMIteration(int maxiteration);
	/// get target image
	irtkGreyImage * getTarget();
	/// get source image
	irtkGreyImage * getSource(int n);
	/// set source images
	void setSources(irtkGreyImage **inputs);
	/// set source images
	void setDecimatedImage(irtkGreyImage *decimated);
	/// set random weight
	void setRandomrate(double rate);
	/// set random weight
	void setDebug(int debug);
	/// Use the nearest neighbor to create label images
	void generateLabels(irtkGreyImage *label, irtkGreyImage **labels);
	/// Use the nearest neighbor to create label images
	void generateImage();
	/// Generate the output image and output it
	void outputmap(char* name);
};

inline irtkGreyImage * irtkPatchMatch::getTarget()
{
	return target;
}

inline irtkGreyImage * irtkPatchMatch::getSource(int n)
{
	return sources[n];
}

inline void irtkPatchMatch::setSources(irtkGreyImage **inputs)
{
	sources = inputs;
}

inline void irtkPatchMatch::setRandomrate(double rate){
	randomrate = rate;
}

inline void irtkPatchMatch::setDebug(int d){
	debug = d;
}

inline NearstNeighbor** irtkPatchMatch::getNearstNeighbor()
{
	return nnfs;
}

#endif
