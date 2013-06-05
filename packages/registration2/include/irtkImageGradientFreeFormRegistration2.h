/*=========================================================================

  Library   : Image Registration Toolkit (IRTK)
  Module    : $Id$
  Copyright : Imperial College, Department of Computing
              Visual Information Processing (VIP), 2009 onwards
  Date      : $Date$
  Version   : $Revision$
  Changes   : $Author$

=========================================================================*/

#ifndef _IRTKIMAGEGRADIENTFREEFORMREGISTRATION2_H

#define _IRTKIMAGEGRADIENTFREEFORMREGISTRATION2_H

class irtkImageGradientFreeFormRegistration2 : public irtkImageFreeFormRegistration2
{

protected:

  /// 4D "Target" image : Normalised gradients of the original target (3 images in target space)
  irtkGenericImage<double> _normalisedGradientTarget;

  /// 4D "Source" image : Normalised gradients of the original source (3 images in source space)
  irtkGenericImage<double> _normalisedGradientSource;

  /// Transformed 4D "Source" image
  irtkGenericImage<double> _transformedNormalisedGradientSource;

  /// 4D "Source" image X, Y and Z gradients
  irtkGenericImage<double> _normalisedGradientSourceGradient[3];

  /// Transformed 4D "Source" image X, Y and Z gradients
  irtkGenericImage<double> _transformedNormalisedGradientSourceGradient[3];

public:

  /// Constructor
  irtkImageGradientFreeFormRegistration2();

  /// Destructor
  virtual ~irtkImageGradientFreeFormRegistration2();

  /// Initial set up for the registration
  virtual void Initialize(int);

  /// Set output for the registration filter
  virtual void SetOutput(irtkTransformation *);

  /// Spherical linear interpolation. Image result is stored in the first vector.
  virtual void Slerp(double *, double *, double *, double, double);

  /// Spherical linear interpolation of gradient. Image result is stored in the first vector.
  virtual void SlerpGradient(double *, double *, double *, double *, double *, double, double);

  /// Update state of the registration based on current transformation estimate (4D source image)
  virtual void UpdateSource();

  /// Update state of the registration based on current transformation estimate (4D source image and 4D source image gradient)
  virtual void UpdateSourceAndGradient();

  /// Evaluate the similarity measure
  virtual double Evaluate();

  /// Evaluate the gradient of the similarity measure
  virtual double EvaluateGradient(double *);

  /// Evaluate the gradient of the similarity measure for the current transformation for 2D images
  virtual void EvaluateGradient2D(double *);

  /// Evaluate the gradient of the similarity measure for the current transformation for 3D images
  virtual void EvaluateGradient3D(double *);

  /// Evaluate similarity measure: NGD
  virtual double EvaluateNGD();

  /// Evaluate gradient of similarity measure: NGD
  virtual void EvaluateGradientNGD();

  /// Evaluate similarity measure: NGP
  virtual double EvaluateNGP();

  /// Evaluate gradient of similarity measure: NGP
  virtual void EvaluateGradientNGP();

  /// Evaluate similarity measure: NGS
  virtual double EvaluateNGS();

  /// Evaluate gradient of similarity measure: NGS
  virtual void EvaluateGradientNGS();

  /// Returns the name of the class
  virtual const char *NameOfClass();
};

inline void irtkImageGradientFreeFormRegistration2::Slerp(double *result, double *v1, double *v2, double w1, double w2)
{
  double dot, angle, v1norm, v2norm, w1slerp, w2slerp;

  if ((*v1 < -1) && (*v2 < -1)) {
	result[0] = -2;
	result[1] = -2;
	result[2] = -2;
  } else if (*v1 < -1) {
	result[0] = w2 * v2[0];
	result[1] = w2 * v2[1];
	result[2] = w2 * v2[2];
  } else if (*v2 < -1) {
	result[0] = w1 * v1[0];
	result[1] = w1 * v1[1];
	result[2] = w1 * v1[2];
  } else {

	result[0] = w1 * v1[0] + w2 * v2[0];
	result[1] = w1 * v1[1] + w2 * v2[1];
	result[2] = w1 * v1[2] + w2 * v2[2];

	v1norm = sqrt((v1[0]*v1[0]) + (v1[1]*v1[1]) + (v1[2]*v1[2]));
	v2norm = sqrt((v2[0]*v2[0]) + (v2[1]*v2[1]) + (v2[2]*v2[2]));

	//Make sure both norms are greater than zero
	if ((v1norm > 0) && (v2norm > 0)) {
	  dot = ((v1[0] * v2[0]) + (v1[1] * v2[1]) + (v1[2] * v2[2])) / (v1norm * v2norm);

	  //If the vectors are not parallel, reweight according to spherical weighting
	  if (fabs(dot) < 0.999) {
		angle = acos(dot);
		w1slerp = sin(w1 * angle) / sin(angle);
		w2slerp = sin(w2 * angle) / sin(angle);

		result[0] = (w1 * v1norm + w2 * v2norm) * (w1slerp * v1[0]/v1norm + w2slerp * v2[0]/v2norm);
		result[1] = (w1 * v1norm + w2 * v2norm) * (w1slerp * v1[1]/v1norm + w2slerp * v2[1]/v2norm);
		result[2] = (w1 * v1norm + w2 * v2norm) * (w1slerp * v1[2]/v1norm + w2slerp * v2[2]/v2norm);
	  }
    }
  }
}

inline void irtkImageGradientFreeFormRegistration2::SlerpGradient(double *result, double *v1_grad, double *v2_grad, double *v1, double *v2, double w1, double w2)
{
  double dot, angle, v1norm, v2norm, w1slerp, w2slerp;

  if ((*v1 < -1) && (*v2 < -1)) {
	result[0] = 0;
	result[1] = 0;
	result[2] = 0;
  } else if (*v1 < -1) {
	result[0] = w2 * v2_grad[0];
	result[1] = w2 * v2_grad[1];
	result[2] = w2 * v2_grad[2];
  } else if (*v2 < -1) {
	result[0] = w1 * v1_grad[0];
	result[1] = w1 * v1_grad[1];
	result[2] = w1 * v1_grad[2];
  } else {
	result[0] = w1 * v1_grad[0] + w2 * v2_grad[0];
	result[1] = w1 * v1_grad[1] + w2 * v2_grad[1];
	result[2] = w1 * v1_grad[2] + w2 * v2_grad[2];

	v1norm = sqrt((v1[0]*v1[0]) + (v1[1]*v1[1]) + (v1[2]*v1[2]));
	v2norm = sqrt((v2[0]*v2[0]) + (v2[1]*v2[1]) + (v2[2]*v2[2]));

	//Make sure both norms are greater than zero
	if ((v1norm > 0) && (v2norm > 0)) {
	  dot = ((v1[0] * v2[0]) + (v1[1] * v2[1]) + (v1[2] * v2[2])) / (v1norm * v2norm);

	  //If the vectors are not parallel, reweight according to spherical weighting
	  if (fabs(dot) < 0.999) {
		angle = acos(dot);
		w1slerp = sin(w1 * angle) / sin(angle);
		w2slerp = sin(w2 * angle) / sin(angle);

		result[0] = (w1 * v1norm + w2 * v2norm) * (w1slerp * v1_grad[0]/v1norm + w2slerp * v2_grad[0]/v2norm);
		result[1] = (w1 * v1norm + w2 * v2norm) * (w1slerp * v1_grad[1]/v1norm + w2slerp * v2_grad[1]/v2norm);
		result[2] = (w1 * v1norm + w2 * v2norm) * (w1slerp * v1_grad[2]/v1norm + w2slerp * v2_grad[2]/v2norm);
	  }
    }
  }
}

inline void irtkImageGradientFreeFormRegistration2::SetOutput(irtkTransformation *transformation)
{
  if ((strcmp(transformation->NameOfClass(),
	             "irtkMultiLevelFreeFormTransformation") != 0) && (strcmp(transformation->NameOfClass(),
	                     "irtkFluidFreeFormTransformation") != 0)){
	cerr << "irtkImageFreeFormRegistration::SetOutput: Transformation must be "
	     << "irtkMultiLevelFreeFormTransformation or irtkFluidFreeFormTransformation" << endl;
	exit(0);
  }

  _transformation = transformation;
}

inline const char *irtkImageGradientFreeFormRegistration2::NameOfClass()
{
  return "irtkImageGradientFreeFormRegistration2";
}

#endif
