/*=========================================================================

  Library   : Image Registration Toolkit (IRTK)
  Module    : $Id$
  Copyright : Imperial College, Department of Computing
              Visual Information Processing (VIP), 2009 onwards
  Date      : $Date$
  Version   : $Revision$
  Changes   : $Author$

=========================================================================*/

#include <irtkImage.h>

#include <irtkRegistration.h>

#ifdef HAS_VTK

#define EPSILON 0.001

irtkModelCorrelationSimilarityMetric::irtkModelCorrelationSimilarityMetric(irtkGreyImage *image, int n, double ds) : irtkModelSimilarityMetric(image)
{
	// Setup the interpolator
	_interpolator = new irtkLinearInterpolateImageFunction;

	// Setup interpolation for the source image
	_interpolator->SetInput(_image);
	_interpolator->Initialize();

	// Calculate the image domain in which we can interpolate
	_interpolator->Inside(_x1, _y1, _z1, _x2, _y2, _z2);

	// Set up profile
	_n  = n;
	_ds = ds;
	_profile = new double[2 * n + 1];
}

irtkModelCorrelationSimilarityMetric::~irtkModelCorrelationSimilarityMetric()
{
	delete _profile;
	delete _interpolator;
}

double irtkModelCorrelationSimilarityMetric::Evaluate(double *point,	double *normal, double *profile)
{
	int i;
	double x, y, z, mean, var, ncc;

	// Transform point to image coordinates
	_image->WorldToImage(point[0], point[1], point[2]);

	if ((point[0] > _x1) && (point[0] < _x2) && (point[1] > _y1) && (point[1]
			< _y2) && (point[2] > _z1) && (point[2] < _z2)) {

		// Compute means and variances
		mean  = 0;
		var   = 0;
		_mean = 0;
		_var  = 0;

		for (i = 0; i < 2 * _n + 1; i++) {
			// Compute profile value
			x = point[0] + (i - _n) * _ds * normal[0];
			y = point[1] + (i - _n) * _ds * normal[1];
			z = point[2] + (i - _n) * _ds * normal[2];
			_profile[i] = _interpolator->Evaluate(x, y, z);

			// Update image mean and variance
			_mean += _profile[i];
			_var  += _profile[i] * _profile[i];

			// Update profile mean and variance
			mean += profile[i];
			var  += profile[i] * profile[i];
		}

		// Compute means and variances
		mean  =  mean / (2.0 * _n + 1);
		var   =  var  / (2.0 * _n + 1);
		var   =  var  - (mean * mean);
		_mean = _mean / (2.0 * _n + 1);
		_var  = _var  / (2.0 * _n + 1);
		_var  = _var  - (_mean * _mean);

		// Compute normalized cross-correlation
		ncc = 0;
		for (i = 0; i < 2 * _n + 1; i++) {
			ncc += (_profile[i] - _mean) * (profile[i] - mean);
		}
		ncc /= 2 * _n + 1;
		return ncc / (sqrt(var + EPSILON) * sqrt(_var + EPSILON));
	} else {
		return 0;
	}
}

#endif

