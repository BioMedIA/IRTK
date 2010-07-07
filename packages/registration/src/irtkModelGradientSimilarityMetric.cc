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

#define EPSILON 0.00001

irtkModelGradientSimilarityMetric::irtkModelGradientSimilarityMetric(
		irtkGreyImage *image) :
	irtkModelSimilarityMetric(image) {
	int x, y, z;
	double norm, dx, dy, dz;

	// Allocate gradient image
	irtkImageAttributes attr = image->GetImageAttributes();
	attr._t = 4;
	_gradient = new irtkRealImage(attr);

	// Compute gradient image
	for (z = 0; z < image->GetZ(); ++z) {
		for (y = 0; y < image->GetY(); ++y) {
			for (x = 0; x < image->GetX(); ++x) {
				if ((x > 0) && (x < image->GetX() - 1)) {
					dx = image->Get(x - 1, y, z) - image->Get(x + 1, y, z);
				} else {
					dx = 0;
				}
				if ((y > 0) && (y < image->GetY() - 1)) {
					dy = image->Get(x, y - 1, z) - image->Get(x, y + 1, z);
				} else {
					dy = 0;
				}
				if ((z > 0) && (z < image->GetZ() - 1)) {
					dz = image->Get(x, y, z - 1) - image->Get(x, y, z + 1);
				} else {
					dz = 0;
				}

				// Gradient magnitude
				norm = sqrt(dx * dx + dy * dy + dz * dz) + EPSILON;
				_gradient->Put(x, y, z, 0, norm);

				// Gradient direction
				_gradient->Put(x, y, z, 1, dx / norm);
				_gradient->Put(x, y, z, 2, dy / norm);
				_gradient->Put(x, y, z, 3, dz / norm);
			}
		}
	}

	// Setup the interpolator
	_interpolator = new irtkLinearInterpolateImageFunction;

	// Setup interpolation for the source image
	_interpolator->SetInput(_gradient);
	_interpolator->Initialize();

	// Calculate the image domain in which we can interpolate
	_interpolator->Inside(_x1, _y1, _z1, _x2, _y2, _z2);
}

irtkModelGradientSimilarityMetric::~irtkModelGradientSimilarityMetric() {
	delete _gradient;
	delete _interpolator;
}

double irtkModelGradientSimilarityMetric::Evaluate(double *point, double *normal, double *)
{
	double mag, dir[3];

	// Transform point to image coordinates
	_image->WorldToImage(point[0], point[1], point[2]);

	if ((point[0] > _x1) && (point[0] < _x2) && (point[1] > _y1) && (point[1]
			< _y2) && (point[2] > _z1) && (point[2] < _z2)) {
		if (normal != NULL) {
			// Compute gradient magnitude
			mag = _interpolator->EvaluateInside(point[0], point[1], point[2], 0);
			// Compute gradient direction
			dir[0] = _interpolator->EvaluateInside(point[0], point[1], point[2], 1);
			dir[1] = _interpolator->EvaluateInside(point[0], point[1], point[2], 2);
			dir[2] = _interpolator->EvaluateInside(point[0], point[1], point[2], 3);
			return mag * (dir[0] * normal[0] + dir[1] * normal[1] + dir[2]
					* normal[2]);
		} else {
			// Compute gradient magnitude
			mag = _interpolator->EvaluateInside(point[0], point[1], point[2], 0);

			return mag;
		}
	} else {
		return 0;
	}
}

#endif

