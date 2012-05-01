/*=========================================================================

  Library   : Image Registration Toolkit (IRTK)
  Module    : $Id$
  Copyright : Imperial College, Department of Computing
              Visual Information Processing (VIP), 2008 onwards
  Date      : $Date$
  Version   : $Revision$
  Changes   : $Author$

=========================================================================*/

#include <irtkRician.h>

#include <irtkGeometry.h>

void irtkRician::Initialise(const double &mi, const double &sigma)
{
	//intialise
	_mi = mi;
	_sigma = sigma;
	_norm = 1.0 / _sigma;
	return;
}

void irtkRician::Approximate(){
	//Koay, C.G. and Basser, P.J. Analytically exact correction scheme for signal extraction from noisy magnitude MR signals 2003
	double difference,thita,ratio,temp, newthita,square_thita;
	int iteration;
	thita = sqrt(M_PI/(4.0-M_PI));
	ratio = _mi/sqrt(_sigma);
	thita = ratio - thita;
	if(thita < 0){
		cerr << "ratio is smaller than lower boundary, the function will not always converge" << endl;
		temp = 0.429; // min value of the lower boundary of temp
		_mi = sqrt(pow(_mi,2)+(temp-2)*_sigma);
		_sigma = _sigma/temp;
		_norm = 1.0 / _sigma;
		return;
	}
	difference = thita;
	iteration = 0;
	// the fix point equation
	while(difference > 0.0001 && iteration < 500){
		square_thita = pow(thita,2);
		temp = 2 + square_thita - (M_PI/8)*exp(-square_thita/2)
			*pow(((2+square_thita)*bessi0(square_thita/4)
			+square_thita*bessi1(square_thita/4)),2);
		newthita = sqrt(temp*(1+pow(ratio,2))-2);
		difference = fabs(newthita - thita);
		thita = newthita;
		iteration++;
	}
	// final evaluation
	square_thita = pow(thita,2);
	temp = 2 + square_thita - (M_PI/8)*exp(-square_thita/2)
		*pow(((2+square_thita)*bessi0(square_thita/4)
		+square_thita*bessi1(square_thita/4)),2);
	_mi = sqrt(pow(_mi,2)+(temp-2)*_sigma);
	_sigma = _sigma/temp;
	//calculate norm
	_norm = 1.0 / _sigma;
	return;
}

double irtkRician::Evaluate(const double &x)
{
	if(x>0)
		if(_mi>0)
			return _norm * x * exp(-(x*x + _mi*_mi) / (2.0 * _sigma))*Bessel(x*_mi/_sigma);
		else
			//_mi = 0 reduce to Rayleigh distribution
			return _norm * x * exp(-(x*x) / (2.0 * _sigma));
	else
		return 0;
}

double irtkRician::Bessel(const double &v)
{
  return bessi0(v);
}

