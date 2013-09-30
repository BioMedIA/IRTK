/*=========================================================================

Library   : Image Registration Toolkit (IRTK)
Module    : $Id$
Copyright : Imperial College, Department of Computing
Visual Information Processing (VIP), 2008 onwards
Date      : $Date$
Version   : $Revision$
Changes   : $Author$

=========================================================================*/

#ifndef irtkMAPatchMatchSuperResolution_H_
#define irtkMAPatchMatchSuperResolution_H_

class irtkMAPatchMatchSuperResolution  : public irtkMAPatchMatch{

protected:

	/// decimated image
	irtkGreyImage *decimated;
	/// decimated gradient
	irtkGreyImage *decimatedgradient;
	/// blured source image
	irtkGreyImage **blured;
	/// blured source gradient
	irtkGreyImage **bluredgradient;
	/// is decimated
	bool isdecimated;
	/// vote
	double *votevalues;
	/// vote
	double *voteweights;
	/// initialize field's weight
	virtual double initialize();
	/// find minum flow
	virtual double minimizeflow();
	/// find minum flow
	virtual double minimizeflowwithdebug();
	/// calculate distance between patches
	virtual double distance(int x1, int y1, int z1, int x2, int y2, int z2, int n);
	/// vote from one patch to another
	void votepatch(int i, int j, int k, int x, int y, int z, int n, double weight);
	/// vote weight matrix
	virtual void voteweight(int zd, int mode);
	/// expectation maximization
	virtual void EMstep();
public:
	/// constructor
	irtkMAPatchMatchSuperResolution(irtkGreyImage *target, irtkGreyImage **source, int radius = 2, int nimages = 1, int nneighbour = 1);
	/// destructor
	virtual ~irtkMAPatchMatchSuperResolution();
	/// calculate NNF
	void runEMIteration(int maxiteration);
	/// set source images
	void setDecimatedImage(irtkGreyImage *decimated);
	/// Use the nearest neighbor to create label images
	void generateImage();
};

#endif
