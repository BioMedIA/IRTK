/*=========================================================================

Library   : Image Registration Toolkit (IRTK)
Module    : $Id$
Copyright : Imperial College, Department of Computing
Visual Information Processing (VIP), 2008 onwards
Date      : $Date$
Version   : $Revision$
Changes   : $Author$

=========================================================================*/

#ifndef irtkMAPatchMatch_H_
#define irtkMAPatchMatch_H_

struct NearstNeighbor{
	int x;
	int y;
	int z;
	int n;
	double weight;
};

class irtkMAPatchMatch{

protected:

	/// target image
	irtkGreyImage *target;
	/// target gradient
	irtkGreyImage *targetgradient;
	/// source image
	irtkGreyImage **sources;
	/// source gradient
	irtkGreyImage **sourcesgradient;
	/// source gradient
	irtkGreyImage **search;
	/// NNF
	NearstNeighbor **nnfs;
	/// Random rate
	double randomrate;
	/// Debug
	int debug;
	/// Radius
	int radius;
	/// Number of source images
	int nimages;
	/// Number of neighborhoods
	int nneighbour;
	/// Max distance between patches
	int maxdistance;
	/// initialize field's weight
	virtual double initialize();
	/// initial guess of the mapping
	virtual void initialguess();
	/// find minum flow
	virtual double minimizeflow();
	/// find minum flow
	virtual double minimizeflowwithdebug();
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
	/// vote weight matrix
	virtual void voteweight(int zd, int mode);
public:
	/// constructor
	irtkMAPatchMatch(irtkGreyImage *target, irtkGreyImage **source, int radius = 2, int nimages = 1, int nneighbour = 1);
	/// destructor
	virtual ~irtkMAPatchMatch();
	/// get the NNF
	NearstNeighbor** getNearstNeighbor();
	/// upsample NNF to refrence image
	void upSampleNNF(irtkMAPatchMatch *reference);
	/// get target image
	irtkGreyImage * getTarget();
	/// get source image
	irtkGreyImage * getSource(int n);
	/// set source images
	void setSources(irtkGreyImage **inputs);
	/// set random weight
	void setRandomrate(double rate);
	/// calculate NNF
	void run(int maxiteration);
	/// set random weight
	void setDebug(int debug);
	/// Generate the output image and output it
	void outputmap(char* name);
};

inline irtkGreyImage * irtkMAPatchMatch::getTarget()
{
	return target;
}

inline irtkGreyImage * irtkMAPatchMatch::getSource(int n)
{
	return sources[n];
}

inline void irtkMAPatchMatch::setSources(irtkGreyImage **inputs)
{
	sources = inputs;
}

inline void irtkMAPatchMatch::setRandomrate(double rate){
	randomrate = rate;
}

inline void irtkMAPatchMatch::setDebug(int d){
	debug = d;
}

inline NearstNeighbor** irtkMAPatchMatch::getNearstNeighbor()
{
	return nnfs;
}

#endif
