/*=========================================================================

Library   : Image Registration Toolkit (IRTK)
Module    : $Id: irtkShapeBasedInterpolateImageFunction.h 8 2009-03-02 16:12:58Z dr $
Copyright : Imperial College, Department of Computing
Visual Information Processing (VIP), 2008 onwards
Date      : $Date: 2009-03-02 16:12:58 +0000 (一, 02 三月 2009) $
Version   : $Revision: 8 $
Changes   : $Author: dr $

=========================================================================*/

#ifndef _IRTKSHAPEBASEDINTERPOLATEIMAGEFUNCTION_H

#define _IRTKSHAPEBASEDINTERPOLATEIMAGEFUNCTION_H

/**
* Class for nearest neighbor interpolation of images
*
* This class defines and implements the nearest neighbor interpolation  of
* images.
*/

class irtkShapeBasedInterpolateImageFunction : public irtkInterpolateImageFunction
{

private:
	/// threshold map for input image
	irtkRealImage _tinput;

	/// distance map for input image
	irtkRealImage _dmap;

	/// Isotropic resampled distance map for input image with linear interpolation
	irtkRealImage _rdmap;

	/// Isotropic resampled Input image after shape based interpolation
	irtkRealImage _rinput;

	/// Isotropic resampled cache image for refine procedure
	irtkRealImage _rcdmap;

	/// Dimension of input image in X-direction
	int _x;

	/// Dimension of input image in Y-direction
	int _y;

	/// Dimension of input image in Z-direction
	int _z;

	/// Dimension of interpolated input image in X-direction
	int _rx;

	/// Dimension of interpolated  input image in Y-direction
	int _ry;

	/// Dimension of interpolated  input image in Z-direction
	int _rz;

	/// Offsets for fast pixel access
	int _offset1, _offset2, _offset3, _offset4;
	int _offset5, _offset6, _offset7, _offset8;

public:

	/// Constructor
	irtkShapeBasedInterpolateImageFunction();

	/// Destructor
	~irtkShapeBasedInterpolateImageFunction();

	/// Returns the name of the class
	virtual const char *NameOfClass();

	/// Initialize
	virtual void Initialize();

	/// Initialize second step, fix the union property.
	virtual void Refine();

	/// Evaluate
	virtual double Evaluate(double, double, double, double = 0);

	/** Evaluate the filter at an arbitrary image location (in pixels) without
	*  handling boundary conditions. This version is faster than the method
	*  above, but is only defined inside the image domain. */
	virtual double EvaluateInside(double, double, double, double = 0);

	/// Evaluate
	virtual double EvaluateLinear(double, double, double, double = 0);

	/** Evaluate the filter at an arbitrary image location (in pixels) without
	*  handling boundary conditions. This version is faster than the method
	*  above, but is only defined inside the image domain. */
	virtual double EvaluateInsideLinear(double, double, double, double = 0);

};

#endif
