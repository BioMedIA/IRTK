/*
 * irtkPolynomialBiasField.h
 *
 *  Created on: 7 Jul 2011
 *      Author: cl6311
 */

#ifndef IRTKPOLYNOMIALBIASFIELD_H_
#define IRTKPOLYNOMIALBIASFIELD_H_

#define IRTKBIASFIELD_POLYNOMIAL 123

#include <irtkBiasField.h>

class irtkPolynomialBiasField : public irtkBiasField
{
private:
	int _dop;
	double* _coeff;
	int _numOfCoefficients;
	irtkRealImage *_mask;

public:
	irtkPolynomialBiasField();

	/**
	 * @param dop max degree of polynomial
	 */
	irtkPolynomialBiasField(const irtkGreyImage &image, int dop);
	~irtkPolynomialBiasField();

	virtual void SetMask( irtkRealImage *);

	/// Calculate weighted least square fit of polynomial to data
	virtual void WeightedLeastSquares(double *x1, double *y1, double *z1, double *bias, double *weights, int no);

	double Bias(double, double, double);

	/// Returns a string with the name of the instantiated class
	const char *NameOfClass();

	double Approximate(double *, double *, double *, double *, int);
	void Interpolate(double* dbias);

	/// Subdivide FFD
	void Subdivide();

	/// Reads FFD from file
	void Read (char *);

	/// Writes FFD to file
	virtual void Write(char *);

	/// Print info
	virtual void Print();

private:
	double evaluatePolynomial(double x, double y, double z);
	int getNumberOfCoefficients(int dop);
};

inline const char *irtkPolynomialBiasField::NameOfClass()
{
  return "irtkPolynomialBiasField";
}
#endif /* IRTKPOLYNOMIALBIASFIELD_H_ */
