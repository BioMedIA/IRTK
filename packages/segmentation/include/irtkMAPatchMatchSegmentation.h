/*=========================================================================

Library   : Image Registration Toolkit (IRTK)
Module    : $Id$
Copyright : Imperial College, Department of Computing
Visual Information Processing (VIP), 2008 onwards
Date      : $Date$
Version   : $Revision$
Changes   : $Author$

=========================================================================*/

#ifndef irtkMAPatchMatchSEGMENTATION_H_
#define irtkMAPatchMatchSEGMENTATION_H_

class irtkMAPatchMatchSegmentation  : public irtkMAPatchMatch{

protected:
	/// target gradient
	irtkRealImage *targetlabeldistance;
	/// source gradient
	irtkRealImage **sourcelabeldistance;
	/// spatial weight
	double sweight;
	/// normalizer for distance and intensity
	double *distancenorm;
	/// label image
	irtkGreyImage *label;
	/// label image
	irtkGreyImage **labels;
	/// vote
	double *votelabels;
	/// Min Max labels
	short minlabel;
	short maxlabel;
	/// initial guess of the mapping
	virtual void initialguess();
	/// calculate distance between patches
	virtual double distance(int x1, int y1, int z1, int x2, int y2, int z2, int n);
	/// find minum flow
	virtual double minimizeflow();
	/// vote weight matrix
	virtual void voteweight(int zd, int mode);
	/// vote label from one patch to another
	void votelabel(int i, int j, int k, int x, int y, int z, int n, double weight);
	/// normalize the distance with respect to the image
	void normalizedistances();
public:
	/// constructor
	irtkMAPatchMatchSegmentation(irtkGreyImage *target, irtkGreyImage **source, 
		irtkRealImage *targetlabeldistance, irtkRealImage **sourcelabeldistance,
		irtkGreyImage *label, irtkGreyImage **labels,
		int radius = 2, int nimages = 1, int nneighbour = 1);

	virtual ~irtkMAPatchMatchSegmentation();

	/// set weight
	virtual void setWeight(double weight);
	/// Use the nearest neighbor to create label images
	void generateLabels();
};

#endif
