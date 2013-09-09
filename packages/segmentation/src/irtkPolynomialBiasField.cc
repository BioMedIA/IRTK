/*=========================================================================

  Library   : Image Registration Toolkit (IRTK)
  Module    : $Id$
  Copyright : Imperial College, Department of Computing
              Visual Information Processing (VIP), 2011 onwards
  Date      : $Date$
  Version   : $Revision$
  Changes   : $Author$

=========================================================================*/



#include <irtkBiasField.h>
#include <irtkMatrix.h>
#include <irtkVector.h>
#include <irtkPolynomialBiasField.h>


irtkPolynomialBiasField::irtkPolynomialBiasField()
{
}

irtkPolynomialBiasField::irtkPolynomialBiasField(const irtkGreyImage &image, int dop)
{
	_dop = dop;
	_numOfCoefficients = getNumberOfCoefficients(dop);
	_coeff = new double[_numOfCoefficients];
	memset( _coeff, 0, sizeof(double) * _numOfCoefficients );
}

void irtkPolynomialBiasField::SetMask( irtkRealImage* imagePtr)
{
	_mask = imagePtr;
}

irtkPolynomialBiasField::~irtkPolynomialBiasField()
{
	delete _coeff;
}

// implementation using symmetry of A'WA, should be by a factor of 2 faster than standard implementation (commented out at eof)
void irtkPolynomialBiasField::WeightedLeastSquares(double *x1, double *y1, double *z1, double *bias, double *weights, int no)
{
	// just consider each eachs voxel...
	int each = 1;
	no /= each;

	irtkMatrix A(no, _numOfCoefficients);
	irtkVector vecB(no);

	cout << "numOfCoefficients: " << _numOfCoefficients << endl;
	cout << "num of voxels: " << no << endl;

	double* AtWA = new double[_numOfCoefficients * _numOfCoefficients];
	memset(AtWA, 0, sizeof(double) * _numOfCoefficients * _numOfCoefficients );

	double* Basis = new double[_numOfCoefficients];

	for( int rr = 0; rr < no; ++rr )
	{
		int r = rr * each;
		double weight = weights[rr*each];
		vecB.Put(rr,bias[r]);

		double x = x1[r];
		double y = y1[r];
		double z = z1[r];

		int c = 0;

		double cur_x = 1.0;
		double cur_y = 1.0;

		for( int xd = 0; xd <= _dop; ++xd )
		{
			cur_y = 1.0;
			for( int yd = 0; yd <= _dop-xd; ++yd )
			{
				double tmp = cur_x*cur_y;
				for( int zd = 0; zd <= _dop-xd-yd; ++zd )
				{
					A.Put(rr,c, tmp);
					Basis[c] = tmp;
					c++;
					tmp *= z;
				}
				cur_y *= y;
			}
			cur_x *= x;
		}

		double* Aptr = (double*) AtWA;
		double* Basisptr1 = (double*) Basis;
		double* Basisptr2;

		for( int j2 = 0; j2 < _numOfCoefficients; j2++, Basisptr1++)
		{
	         Basisptr2= &Basis[j2];
	         Aptr= &AtWA[j2+j2*_numOfCoefficients];
	         for(int i2=j2; i2<_numOfCoefficients; i2++, Aptr++, Basisptr2++){
	             (*Aptr)+=(*Basisptr2)*(weight)*(*Basisptr1);
	         }
		}
	}

	irtkMatrix AtW = A;
	AtW.Transpose();

	for( int r = 0; r < _numOfCoefficients; ++r)
	{
		for( int c= 0; c < no; ++c )
		{
			AtW.Put(r,c, AtW(r,c) * weights[c*each] );
		}
	}

	irtkMatrix leftSide(_numOfCoefficients, _numOfCoefficients);
    for(int j2=0; j2<_numOfCoefficients; j2++)
    {
    	for(int i2=j2; i2<_numOfCoefficients; i2++)
    	{
    		leftSide.Put(i2,j2,(double)(AtWA[i2+j2*_numOfCoefficients]));
    		leftSide.Put(j2,i2,(double)(AtWA[i2+j2*_numOfCoefficients]));
    	}
    }

 	irtkVector rightSide = AtW * vecB;
 	leftSide.Invert();

	irtkVector vecC = leftSide * rightSide;

	for( int r = 0; r < _numOfCoefficients; ++r )
	{
		_coeff[r] = vecC.Get(r);
	}
}

double irtkPolynomialBiasField::evaluatePolynomial(double x, double y, double z)
{
	double res = 0;
	int n = 0;

	double cur_x = 1.0;
	double cur_y = 1.0;

	for( int xd = 0; xd <= _dop; ++xd )
	{
		cur_y = 1.0;
		for( int yd = 0; yd <= _dop-xd; ++yd )
		{
			double tmp = cur_x * cur_y;
			for( int zd = 0; zd <= _dop-xd-yd; ++zd )
			{
				res +=   tmp * _coeff[n];
				n++;
				tmp *= z;
			}
			cur_y *= y;
		}
		cur_x *= x;
	}
	return res;
}


double irtkPolynomialBiasField::Bias(double x, double y, double z)
{
	return evaluatePolynomial(x, y, z);
}

void irtkPolynomialBiasField::Interpolate(double* dbias)
{

}

void irtkPolynomialBiasField::Subdivide()
{

}

void irtkPolynomialBiasField::Read(char *name)
{
  unsigned int magic_no;
  unsigned int trans_type;

  // Open file
  irtkCifstream from;
  from.Open(name);

  // Read magic no. for transformations
  from.ReadAsUInt(&magic_no, 1, 0);

  // Read transformation type
  from.ReadAsUInt(&trans_type, 1);

  if ((magic_no != IRTKBIASFIELD_MAGIC) && (trans_type != IRTKBIASFIELD_POLYNOMIAL)) {
    cerr << "irtkPolynomialBiasField::Read: File format not recognized" << endl;
    exit(1);
  }

  if( _coeff ) delete _coeff;

  // Write data
  from.ReadAsInt(&_dop, 1);
  from.ReadAsInt(&_numOfCoefficients, 1);
  _coeff = new double[_numOfCoefficients];

  from.ReadAsDouble(_coeff, _numOfCoefficients);

  // Close file stream
  from.Close();
}


void irtkPolynomialBiasField::Write(char *name)
{
  // Open file
  irtkCofstream to;
  to.Open(name);

  // Write magic no. for transformations
  unsigned int magic_no = IRTKBIASFIELD_MAGIC;
  to.WriteAsUInt(&magic_no, 1, 0);

  // Write transformation type
  unsigned int trans_type = IRTKBIASFIELD_POLYNOMIAL;
  to.WriteAsUInt(&trans_type, 1);

  // Write data
  to.WriteAsInt(&_dop, 1);
  to.WriteAsInt(&_numOfCoefficients, 1);
  to.WriteAsDouble(_coeff, _numOfCoefficients);

  // Close file stream
  to.Close();
}

void irtkPolynomialBiasField::Print()
{
  cerr<<endl<<"Polynomial Bias Field info:" << endl;
  // Write no. of control points
  cout << "Control points: " << _x << " x " << _y << " x " << _z << endl;
  cout << "Spacing: " << _dx << " x " << _dy << " x " << _dz << endl;
  cout << "Origin: " << _origin._x << " " << _origin._y << " " << _origin._z << " " << endl;
  cout << "Orientation: " << _xaxis[0] << " " << _xaxis[1] << " "
       << _xaxis[2] << " " << _yaxis[0] << " " << _yaxis[1] << " "
       << _yaxis[2] << " " << _zaxis[0] << " " << _zaxis[1] << " " << _zaxis[2] << endl;
  cout << _numOfCoefficients << " Coefficients:" << endl;

  int n;
  for (n=0; n<_numOfCoefficients; n++) {
	  cout << _coeff[n] << " ";
  }
  cout<<endl;
}


double irtkPolynomialBiasField::Approximate(double *x1, double *y1, double *z1, double *bias, int no)
{
  return .0;
}

int irtkPolynomialBiasField::getNumberOfCoefficients(int dop)
{
	int n = 0;
	for( int x = dop; x >=0; --x )
	{
		for( int y = dop-x; y >=0; --y )
		{
			for( int z = dop-x-y; z >=0; --z )
			{
				n++;
			}
		}
	}
	return n;
}



//void irtkPolynomialBiasField::WeightedLeastSquares(double *x1, double *y1, double *z1, double *bias, double *weights, int no)
//{
//	// just consider each eachs voxel...
//	int each = 1;
//	no /= each;
//
//	irtkMatrix A(no, _numOfCoefficients);
//	irtkVector vecB(no);
//
//	cout << "numOfCoefficients: " << _numOfCoefficients << endl;
//	cout << "num of voxels: " << no << endl;
//
//	for( int rr = 0; rr < no; ++rr )
//	{
//		int r = rr * each;
//		vecB.Put(rr,bias[r]);
//
//		double x = x1[r];
//		double y = y1[r];
//		double z = z1[r];
//
//		int c = 0;
//
//		double cur_x = pow(x, _dop);
//		for( int xd = _dop; xd >= 0; --xd)
//		{
//			double cur_y = pow(y,_dop-xd);
//			for( int yd = _dop-xd; yd >= 0; --yd)
//			{
//				double cur_z = cur_x*cur_y * pow(z,_dop-xd-yd);
//				for( int zd = _dop-xd-yd; zd >= 0; --zd)
//				{
//					A.Put(rr,c, cur_z);
//					c++;
//					if( z != 0) cur_z /= z;
//				}
//				if( y != 0) cur_y /= y;
//			}
//			if( x != 0) cur_x /= x;
//		}
//	}
//
//	irtkMatrix AtW = A;
//	AtW.Transpose();
//
//	for( int r = 0; r < _numOfCoefficients; ++r)
//	{
//		for( int c= 0; c < no; ++c )
//		{
//			AtW.Put(r,c, AtW(r,c) * weights[c*each] );
//		}
//	}
//
// 	irtkVector rightSide = AtW * vecB;
// 	irtkMatrix leftSide = AtW * A;
//
// 	leftSide.Invert();
//
//	irtkVector vecC = leftSide * rightSide;
//
//	for( int r = 0; r < _numOfCoefficients; ++r )
//	{
//		_coeff[r] = vecC.Get(r);
//	}
//}
