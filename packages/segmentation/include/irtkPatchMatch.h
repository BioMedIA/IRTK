/*=========================================================================

Library   : Image Registration Toolkit (IRTK)
Module    : $Id$
Copyright : Imperial College, Department of Computing
Visual Information Processing (VIP), 2008 onwards
Date      : $Date$
Version   : $Revision$
Changes   : $Author$

=========================================================================*/

#ifndef irtkPatchMatch_H_
#define irtkPatchMatch_H_

struct NearstNeighborFlow{
	int x;
	int y;
	int z;
	double weight;
};

class irtkPatchMatch{

protected:

	/// target image
	irtkGreyImage *target;
	/// source image
	irtkGreyImage *source;
	/// NNF
	NearstNeighborFlow **nnfs;
	/// Random rate
	double randomrate;
	/// Debug
	int debug;
	/// Radius
	int radius;
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
	/// calculate distance between patches of different image
	virtual double distance(int x1, int y1, int z1, int x2, int y2, int z2);
	/// calculate distance between patches of different image
	virtual double distance3D(int x1, int y1, int z1, int x2, int y2, int z2);
	/// calculate distance between patches of different image
	virtual double distance2D(int x1, int y1, int x2, int y2);
	/// calculate distance between patches of same image
	virtual double selfdistance(int x1, int y1, int z1, int x2, int y2, int z2, int mode = 0);
	/// calculate distance between patches of same image
	virtual double selfdistance3D(int x1, int y1, int z1, int x2, int y2, int z2, int mode = 0);
	/// calculate distance between patches of same image
	virtual double selfdistance2D(int x1, int y1, int x2, int y2, int mode = 0);
	/// random search the space to find a better link
	int randomlink(int i, int j, int k, int n, int index = -1);
	/// propergate from one to another using the offset
	int propergate(int x, int y, int z, int i, int j, int k, int offsetx, int offsety, int offestz, int index1 = -1, int index2 = -1);
	/// vote weight matrix
	virtual void voteweight(int zd, int mode);
public:
	/// constructor
	irtkPatchMatch(irtkGreyImage *target, irtkGreyImage *source, int radius = 2, int nneighbour = 1);
	/// destructor
	virtual ~irtkPatchMatch();
	/// get the NNF
	NearstNeighborFlow** getNearstNeighbor();
	/// get target image
	irtkGreyImage * getTarget();
	/// get source image
	irtkGreyImage * getSource();
	/// set random weight
	void setRandomrate(double rate);
	/// calculate NNF
	void run(int maxiteration);
	/// set debug info
	void setDebug(int debug);
	/// Generate the output image and output it
	void outputmap(char* name);
	/// Generate the distance graph between two images
	void outputgraph(char* name);
};

inline irtkGreyImage * irtkPatchMatch::getTarget()
{
	return target;
}

inline irtkGreyImage * irtkPatchMatch::getSource()
{
	return source;
}

inline void irtkPatchMatch::setRandomrate(double rate){
	randomrate = rate;
}

inline void irtkPatchMatch::setDebug(int d){
	debug = d;
}

inline NearstNeighborFlow** irtkPatchMatch::getNearstNeighbor()
{
	return nnfs;
}

#endif
