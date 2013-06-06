/*=========================================================================

Library   : Image Registration Toolkit (IRTK)
Module    : $Id$
Copyright : Imperial College, Department of Computing
Visual Information Processing (VIP), 2008 onwards
Date      : $Date$
Version   : $Revision$
Changes   : $Author$

=========================================================================*/

#ifndef IRTKPATCHMATCHSEGMENTATION_H_
#define IRTKPATCHMATCHSEGMENTATION_H_

class irtkPatchMatchSegmentation  : public irtkPatchMatch{

protected:
	/// target gradient
	irtkRealImage *targetlabeldistance;
	/// source gradient
	irtkRealImage **sourcelabeldistance;
	/// spatial weight
	double sweight;
	/// normalizer for distance and intensity
	double *distancenorm;
	/// initial guess of the mapping
	virtual void initialguess();
	/// calculate distance between patches
	virtual double distance(int x1, int y1, int z1, int x2, int y2, int z2, int n);
	/// expectation maximization
	virtual void EMstep();
	/// normalize the distance with respect to the image
	void normalizedistances();
public:
	/// constructor
	irtkPatchMatchSegmentation(irtkGreyImage *target, irtkGreyImage **source, 
		irtkRealImage *targetlabeldistance, irtkRealImage **sourcelabeldistance,
		int radius = 2, int nimages = 1, int nneighbour = 1);

	virtual ~irtkPatchMatchSegmentation();

	/// set weight
	void setWeight(double weight);
};

#endif
