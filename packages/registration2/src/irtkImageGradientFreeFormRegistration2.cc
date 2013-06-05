/*=========================================================================

  Library   : Image Registration Toolkit (IRTK)
  Module    : $Id$
  Copyright : Imperial College, Department of Computing
              Visual Information Processing (VIP), 2009 onwards
  Date      : $Date$
  Version   : $Revision$
  Changes   : $Author$

=========================================================================*/

#include <irtkRegistration2.h>

#include <irtkGradientImageFilter.h>

#ifdef WIN32
#include <time.h>
#else
#include <sys/resource.h>
#endif

#define EPS_NEIGH 3.0

irtkImageGradientFreeFormRegistration2::irtkImageGradientFreeFormRegistration2() : irtkImageFreeFormRegistration2()
{
  _SimilarityMeasure = NGD;
}

irtkImageGradientFreeFormRegistration2::~irtkImageGradientFreeFormRegistration2()
{}

void irtkImageGradientFreeFormRegistration2::Initialize(int level)
{
  int i, j, k, p, q, r, n;
  double norm, eps_target, eps_source;

  // Print debugging information
  this->Debug("irtkImageGradientFreeFormRegistration2::Initialize(int)");

  // Initialize base class
  this->irtkImageFreeFormRegistration2::Initialize(level);

  //Update real transformed source and real transformed source gradient
  this->irtkImageFreeFormRegistration2::UpdateSourceAndGradient();

  //Compute normalised target and source gradients
  irtkGradientImageFilter<double> gradient(irtkGradientImageFilter<double>::GRADIENT_VECTOR);
  irtkGenericImage<double> tmp, targetGradient;

  irtkImageAttributes attr = _target->GetImageAttributes();
  attr._t = 3;
  _normalisedGradientTarget.Initialize(attr);

  tmp = *_target;
  gradient.SetInput (&tmp);
  gradient.SetOutput(&targetGradient);
  gradient.SetPadding(_TargetPadding);
  gradient.Run();

  //Compute normalised gradient for target image
  cout << "Computing normalised gradient target ... "; cout.flush();

  if (_SimilarityMeasure == NGS) {

	//Compute epsilon for target image using the whole image
	eps_target = 0;
	n = 0;
    for (k = 0; k < _target->GetZ(); k++) {
	  for (j = 0; j < _target->GetY(); j++) {
	    for (i = 0; i < _target->GetX(); i++) {
	      if (_target->GetZ() == 1) {
	    	eps_target += sqrt((targetGradient(i, j, 0, 0) * targetGradient(i, j, 0, 0)) + (targetGradient(i, j, 0, 1) * targetGradient(i, j, 0, 1)) + (targetGradient(i, j, 0, 2) * targetGradient(i, j, 0, 2)));
	      } else {
			eps_target += sqrt((targetGradient(i, j, k, 0) * targetGradient(i, j, k, 0)) + (targetGradient(i, j, k, 1) * targetGradient(i, j, k, 1)) + (targetGradient(i, j, k, 2) * targetGradient(i, j, k, 2)));
	      }
	      n++;
	    }
	  }
    }
    if (n > 0) eps_target /= double(n);

    if (eps_target == 0) {
      cerr << "Target image has no structure!" << endl;
      exit(1);
    }

    for (k = 0; k < _target->GetZ(); k++) {
	  for (j = 0; j < _target->GetY(); j++) {
		for (i = 0; i < _target->GetX(); i++) {
		  norm = sqrt((eps_target * eps_target) + (targetGradient(i, j, k, 0) * targetGradient(i, j, k, 0)) + (targetGradient(i, j, k, 1) * targetGradient(i, j, k, 1)) + (targetGradient(i, j, k, 2) * targetGradient(i, j, k, 2)));
		  _normalisedGradientTarget(i, j, k, 0) = targetGradient(i, j, k, 0) / norm;
		  _normalisedGradientTarget(i, j, k, 1) = targetGradient(i, j, k, 1) / norm;
		  _normalisedGradientTarget(i, j, k, 2) = targetGradient(i, j, k, 2) / norm;
		}
	  }
    }

  } else {

    for (k = 0; k < _target->GetZ(); k++) {
      for (j = 0; j < _target->GetY(); j++) {
        for (i = 0; i < _target->GetX(); i++) {
		  eps_target = 1e-4;
		  n = 0;

	      //Compute epsilon for target image using neighbourhood patches
		  if (_target->GetZ() == 1) {
			for (q = j-floor(EPS_NEIGH/2.0); q <= j+floor(EPS_NEIGH/2.0); q++) {
			  for (p = i-floor(EPS_NEIGH/2.0); p <= i+floor(EPS_NEIGH/2.0); p++) {
			    if ((p >= 0) && (p < _target->GetX()) && (q >= 0) && (q < _target->GetY())) {
				  eps_target += sqrt((targetGradient(p, q, 0, 0) * targetGradient(p, q, 0, 0)) + (targetGradient(p, q, 0, 1) * targetGradient(p, q, 0, 1)) + (targetGradient(p, q, 0, 2) * targetGradient(p, q, 0, 2)));
				  n++;
			    }
			  }
			}
	      } else {
			for (r = k-floor(EPS_NEIGH/2.0); r <= k+floor(EPS_NEIGH/2.0); r++) {
			  for (q = j-floor(EPS_NEIGH/2.0); q <= j+floor(EPS_NEIGH/2.0); q++) {
				for (p = i-floor(EPS_NEIGH/2.0); p <= i+floor(EPS_NEIGH/2.0); p++) {
				  if ((p >= 0) && (p < _target->GetX()) && (q >= 0) && (q < _target->GetY()) && (r >= 0) && (r < _target->GetZ())) {
				    eps_target += sqrt((targetGradient(p, q, r, 0) * targetGradient(p, q, r, 0)) + (targetGradient(p, q, r, 1) * targetGradient(p, q, r, 1)) + (targetGradient(p, q, r, 2) * targetGradient(p, q, r, 2)));
				    n++;
				  }
				}
			  }
			}
		  }
		  if (n > 0) eps_target /= double(n);

		  norm = sqrt((eps_target * eps_target) + (targetGradient(i, j, k, 0) * targetGradient(i, j, k, 0)) + (targetGradient(i, j, k, 1) * targetGradient(i, j, k, 1)) + (targetGradient(i, j, k, 2) * targetGradient(i, j, k, 2)));
		  _normalisedGradientTarget(i, j, k, 0) = targetGradient(i, j, k, 0) / norm;
   		  _normalisedGradientTarget(i, j, k, 1) = targetGradient(i, j, k, 1) / norm;
   		  _normalisedGradientTarget(i, j, k, 2) = targetGradient(i, j, k, 2) / norm;
        }
   	  }
   	}
  }

  cout << "done" << endl;

  attr = _source->GetImageAttributes();
  attr._t = 3;
  _normalisedGradientSource.Initialize(attr);

  cout << "Computing normalised gradient source ... "; cout.flush();

  if (_SimilarityMeasure == NGS) {

	//Compute epsilon for source image using the whole image
	eps_source = 0;
	n = 0;
	for (k = 0; k < _source->GetZ(); k++) {
	  for (j = 0; j < _source->GetY(); j++) {
	    for (i = 0; i < _source->GetX(); i++) {
	      if (_source->GetZ() == 1) {
	  	    eps_source += sqrt((_sourceGradient(i, j, 0, 0) * _sourceGradient(i, j, 0, 0)) + (_sourceGradient(i, j, 0, 1) * _sourceGradient(i, j, 0, 1)) + (_sourceGradient(i, j, 0, 2) * _sourceGradient(i, j, 0, 2)));
	  	  } else {
	  	    eps_source += sqrt((_sourceGradient(i, j, k, 0) * _sourceGradient(i, j, k, 0)) + (_sourceGradient(i, j, k, 1) * _sourceGradient(i, j, k, 1)) + (_sourceGradient(i, j, k, 2) * _sourceGradient(i, j, k, 2)));
	  	  }
	      n++;
	  	}
	  }
	}
	if (n > 0) eps_source /= double(n);

	if (eps_source == 0) {
	  cerr << "Source image has no structure!" << endl;
	  exit(1);
	}

	for (k = 0; k < _source->GetZ(); k++) {
	  for (j = 0; j < _source->GetY(); j++) {
		for (i = 0; i < _source->GetX(); i++) {
		  norm = sqrt((eps_source * eps_source) + (_sourceGradient(i, j, k, 0) * _sourceGradient(i, j, k, 0)) + (_sourceGradient(i, j, k, 1) * _sourceGradient(i, j, k, 1)) + (_sourceGradient(i, j, k, 2) * _sourceGradient(i, j, k, 2)));
		  _normalisedGradientSource(i, j, k, 0) = _sourceGradient(i, j, k, 0) / norm;
		  _normalisedGradientSource(i, j, k, 1) = _sourceGradient(i, j, k, 1) / norm;
		  _normalisedGradientSource(i, j, k, 2) = _sourceGradient(i, j, k, 2) / norm;
		}
	  }
	}

  } else {

    for (k = 0; k < _source->GetZ(); k++) {
      for (j = 0; j < _source->GetY(); j++) {
        for (i = 0; i < _source->GetX(); i++) {
		  eps_source = 1e-4;
		  n = 0;

		  //Compute epsilon for source image using neighbourhood patches
	      if (_source->GetZ() == 1) {
			for (q = j-floor(EPS_NEIGH/2.0); q <= j+floor(EPS_NEIGH/2.0); q++) {
			  for (p = i-floor(EPS_NEIGH/2.0); p <= i+floor(EPS_NEIGH/2.0); p++) {
			    if ((p >= 0) && (p < _source->GetX()) && (q >= 0) && (q < _source->GetY())) {
			      eps_source += sqrt((_sourceGradient(p, q, 0, 0) * _sourceGradient(p, q, 0, 0)) + (_sourceGradient(p, q, 0, 1) * _sourceGradient(p, q, 0, 1)) + (_sourceGradient(p, q, 0, 2) * _sourceGradient(p, q, 0, 2)));
			      n++;
			    }
			  }
	        }
		  } else {
			for (r = k-floor(EPS_NEIGH/2.0); r <= k+floor(EPS_NEIGH/2.0); r++) {
			  for (q = j-floor(EPS_NEIGH/2.0); q <= j+floor(EPS_NEIGH/2.0); q++) {
			    for (p = i-floor(EPS_NEIGH/2.0); p <= i+floor(EPS_NEIGH/2.0); p++) {
				  if ((p >= 0) && (p < _source->GetX()) && (q >= 0) && (q < _source->GetY()) && (r >= 0) && (r < _source->GetZ())) {
				    eps_source += sqrt((_sourceGradient(p, q, r, 0) * _sourceGradient(p, q, r, 0)) + (_sourceGradient(p, q, r, 1) * _sourceGradient(p, q, r, 1)) + (_sourceGradient(p, q, r, 2) * _sourceGradient(p, q, r, 2)));
				    n++;
				  }
			    }
			  }
		    }
		  }
	      if (n > 0) eps_source /= double(n);

		  norm = sqrt((eps_source * eps_source) + (_sourceGradient(i, j, k, 0) * _sourceGradient(i, j, k, 0)) + (_sourceGradient(i, j, k, 1) * _sourceGradient(i, j, k, 1)) + (_sourceGradient(i, j, k, 2) * _sourceGradient(i, j, k, 2)));
		  _normalisedGradientSource(i, j, k, 0) = _sourceGradient(i, j, k, 0) / norm;
		  _normalisedGradientSource(i, j, k, 1) = _sourceGradient(i, j, k, 1) / norm;
		  _normalisedGradientSource(i, j, k, 2) = _sourceGradient(i, j, k, 2) / norm;
        }
      }
    }
  }
  cout << "done" << endl;

  if(_DebugFlag == true) {
	cout << "Writing debug files to disk ... "; cout.flush();

    char buffer[500];

    sprintf(buffer, "/homes/sp2010/biomedic/debug/target_%d.nii.gz", _CurrentLevel);
    _normalisedGradientTarget.Write(buffer);

    sprintf(buffer, "/homes/sp2010/biomedic/debug/source_%d.nii.gz", _CurrentLevel);
    _normalisedGradientSource.Write(buffer);

    sprintf(buffer, "/homes/sp2010/biomedic/debug/targetGrad_%d.nii.gz", _CurrentLevel);
    targetGradient.Write(buffer);

    sprintf(buffer, "/homes/sp2010/biomedic/debug/sourceGrad_%d.nii.gz", _CurrentLevel);
    _sourceGradient.Write(buffer);

    sprintf(buffer, "/homes/sp2010/biomedic/debug/targetREAL_%d.nii.gz", _CurrentLevel);
    _target->Write(buffer);

    sprintf(buffer, "/homes/sp2010/biomedic/debug/sourceREAL_%d.nii.gz", _CurrentLevel);
    _source->Write(buffer);

    sprintf(buffer, "/homes/sp2010/biomedic/debug/distanceMask_%d.nii.gz", _CurrentLevel);
    _distanceMask.Write(buffer);

    cout << "done" << endl;
  }

  cout << "Computing source second order derivatives ... "; cout.flush();

  //Compute "source" gradients
  tmp = _normalisedGradientSource.GetFrame(0);
  gradient.SetInput (&tmp);
  gradient.SetOutput(&_normalisedGradientSourceGradient[0]);
  gradient.SetPadding(MIN_GREY);
  gradient.Run();

  tmp = _normalisedGradientSource.GetFrame(1);
  gradient.SetInput (&tmp);
  gradient.SetOutput(&_normalisedGradientSourceGradient[1]);
  gradient.SetPadding(MIN_GREY);
  gradient.Run();

  tmp = _normalisedGradientSource.GetFrame(2);
  gradient.SetInput (&tmp);
  gradient.SetOutput(&_normalisedGradientSourceGradient[2]);
  gradient.SetPadding(MIN_GREY);
  gradient.Run();

  cout << "done" << endl;

  // Set up gradient of "source" image (transformed)
  attr = _target->GetImageAttributes();
  attr._t = 3;
  _transformedNormalisedGradientSourceGradient[0].Initialize(attr);
  _transformedNormalisedGradientSourceGradient[1].Initialize(attr);
  _transformedNormalisedGradientSourceGradient[2].Initialize(attr);
}


void irtkImageGradientFreeFormRegistration2::UpdateSource()
{
  double *ptrX;
  double *ptrY;
  double *ptrZ;
  double *vector1;
  double *vector2;
  double *result1;
  double *result2;
  double *result3;
  double x, y, z, t1, t2, u1, u2, v1, v2;
  int a, b, c, i, j, k, offset1, offset2, offset3, offset4, offset5, offset6, offset7, offset8;

  IRTK_START_TIMING();

  // Print debugging information
  this->Debug("irktImageGradientFreeFormRegistration2::UpdateSource()");

  // Generate transformed tmp image
  _transformedNormalisedGradientSource = _normalisedGradientTarget;

  // Calculate offsets for fast pixel access
  offset1 = 0;
  offset2 = 1;
  offset3 = this->_source->GetX();
  offset4 = this->_source->GetX()+1;
  offset5 = this->_source->GetX()*this->_source->GetY();
  offset6 = this->_source->GetX()*this->_source->GetY()+1;
  offset7 = this->_source->GetX()*this->_source->GetY()+this->_source->GetX();
  offset8 = this->_source->GetX()*this->_source->GetY()+this->_source->GetX()+1;

  vector1 = new double[3];
  vector2 = new double[3];
  result1 = new double[3];
  result2 = new double[3];
  result3 = new double[3];

  double *ptr2disp = _displacementLUT;
  double *ptr2latt = _latticeCoordLUT;
  if ((_target->GetZ() == 1) && (_source->GetZ() == 1)) {
    for (j = 0; j < _target->GetY(); j++) {
      for (i = 0; i < _target->GetX(); i++) {
        if (_distanceMask.Get(i, j, 0) == 0) {
          x = ptr2latt[0];
		  y = ptr2latt[1];
		  z = 0;
		  _affd->FFD2D(x, y);
		  x += ptr2disp[0];
		  y += ptr2disp[1];
		  z += ptr2disp[2];
		  _source->WorldToImage(x, y, z);

          // Check whether transformed point is inside volume
          if ((x > 0) && (x < _source->GetX()-1) &&
              (y > 0) && (y < _source->GetY()-1)) {

            if (_InterpolationMode == Interpolation_Linear) {
              // Calculated integer coordinates
              a  = int(x);
              b  = int(y);

              // Calculated fractional coordinates
              t1 = x - a;
              u1 = y - b;
              t2 = 1 - t1;
              u2 = 1 - u1;

              // Spherical linear interpolation in source image
              ptrX = _normalisedGradientSource.GetPointerToVoxels(a, b, 0, 0);
              ptrY = _normalisedGradientSource.GetPointerToVoxels(a, b, 0, 1);

              // ************************************************************ //

              vector1[0] = ptrX[offset2];
              vector1[1] = ptrY[offset2];
              vector1[2] = 0;

              vector2[0] = ptrX[offset4];
              vector2[1] = ptrY[offset4];
              vector2[2] = 0;

              this->Slerp(result1, vector1, vector2, u2, u1);

              // ************************************************************ //

              vector1[0] = ptrX[offset1];
			  vector1[1] = ptrY[offset1];
			  vector1[2] = 0;

			  vector2[0] = ptrX[offset3];
			  vector2[1] = ptrY[offset3];
			  vector2[2] = 0;

			  this->Slerp(result2, vector1, vector2, u2, u1);

			  // ************************************************************ //

			  this->Slerp(result1, result1, result2, t1, t2);

              _transformedNormalisedGradientSource(i, j, 0, 0) = result1[0];
              _transformedNormalisedGradientSource(i, j, 0, 1) = result1[1];
            } else if (_InterpolationMode == Interpolation_NN) {
              // Calculated NN integer coordinates
			  a  = round(x);
			  b  = round(y);

			  // NN interpolation in source image
			  ptrX = _normalisedGradientSource.GetPointerToVoxels(a, b, 0, 0);
			  ptrY = _normalisedGradientSource.GetPointerToVoxels(a, b, 0, 1);
              _transformedNormalisedGradientSource(i, j, 0, 0) = ptrX[0];
              _transformedNormalisedGradientSource(i, j, 0, 1) = ptrY[0];
            } else {
              cerr << "irktImageGradientFreeFormRegistration2::UpdateSource(): Interpolation not implemented" << endl;
              exit(1);
            }
          } else {
        	_transformedNormalisedGradientSource(i, j, 0, 0) = -2;
        	_transformedNormalisedGradientSource(i, j, 0, 1) = -2;
          }
        } else {
          _transformedNormalisedGradientSource(i, j, 0, 0) = -2;
          _transformedNormalisedGradientSource(i, j, 0, 1) = -2;
        }
        ptr2disp += 3;
        ptr2latt += 3;
      }
    }
  } else {
    for (k = 0; k < _target->GetZ(); k++) {
      for (j = 0; j < _target->GetY(); j++) {
        for (i = 0; i < _target->GetX(); i++) {
          if (_distanceMask.Get(i, j, k) == 0) {
        	x = ptr2latt[0];
			y = ptr2latt[1];
			z = ptr2latt[2];
			_affd->FFD3D(x, y, z);
			x += ptr2disp[0];
			y += ptr2disp[1];
			z += ptr2disp[2];
			_source->WorldToImage(x, y, z);

            // Check whether transformed point is inside volume
            if ((x > 0) && (x < _source->GetX()-1) &&
                (y > 0) && (y < _source->GetY()-1) &&
                (z > 0) && (z < _source->GetZ()-1)) {

			  if (_InterpolationMode == Interpolation_Linear) {
                // Calculated integer coordinates
                a  = int(x);
                b  = int(y);
                c  = int(z);

                // Calculated fractional coordinates
                t1 = x - a;
                u1 = y - b;
                v1 = z - c;
                t2 = 1 - t1;
                u2 = 1 - u1;
                v2 = 1 - v1;

                // Spherical linear interpolation in source image
			    ptrX = _normalisedGradientSource.GetPointerToVoxels(a, b, c, 0);
			    ptrY = _normalisedGradientSource.GetPointerToVoxels(a, b, c, 1);
			    ptrZ = _normalisedGradientSource.GetPointerToVoxels(a, b, c, 2);

			    // ************************************************************ //

			    vector1[0] = ptrX[offset2];
			    vector1[1] = ptrY[offset2];
			    vector1[2] = ptrZ[offset2];

			    vector2[0] = ptrX[offset6];
			    vector2[1] = ptrY[offset6];
			    vector2[2] = ptrZ[offset6];

			    this->Slerp(result1, vector1, vector2, v2, v1);

			    // ************************************************************ //

			    vector1[0] = ptrX[offset4];
			    vector1[1] = ptrY[offset4];
			    vector1[2] = ptrZ[offset4];

			    vector2[0] = ptrX[offset8];
			    vector2[1] = ptrY[offset8];
			    vector2[2] = ptrZ[offset8];

			    this->Slerp(result2, vector1, vector2, v2, v1);

			    // ************************************************************ //

			    this->Slerp(result1, result1, result2, u2, u1);

			    // ************************************************************ //

			    vector1[0] = ptrX[offset1];
				vector1[1] = ptrY[offset1];
				vector1[2] = ptrZ[offset1];

				vector2[0] = ptrX[offset5];
				vector2[1] = ptrY[offset5];
				vector2[2] = ptrZ[offset5];

				this->Slerp(result2, vector1, vector2, v2, v1);

				// ************************************************************ //

				vector1[0] = ptrX[offset3];
				vector1[1] = ptrY[offset3];
				vector1[2] = ptrZ[offset3];

				vector2[0] = ptrX[offset7];
				vector2[1] = ptrY[offset7];
				vector2[2] = ptrZ[offset7];

				this->Slerp(result3, vector1, vector2, v2, v1);

				// ************************************************************ //

				this->Slerp(result2, result2, result3, u2, u1);

				// ************************************************************ //

				this->Slerp(result1, result1, result2, t1, t2);

				_transformedNormalisedGradientSource(i, j, k, 0) = result1[0];
				_transformedNormalisedGradientSource(i, j, k, 1) = result1[1];
				_transformedNormalisedGradientSource(i, j, k, 2) = result1[2];
			  } else if (_InterpolationMode == Interpolation_NN) {
				// Calculated NN integer coordinates
				a  = round(x);
				b  = round(y);
				c  = round(z);

				// NN interpolation in source image
				ptrX = _normalisedGradientSource.GetPointerToVoxels(a, b, c, 0);
				ptrY = _normalisedGradientSource.GetPointerToVoxels(a, b, c, 1);
				ptrZ = _normalisedGradientSource.GetPointerToVoxels(a, b, c, 2);
				_transformedNormalisedGradientSource(i, j, k, 0) = ptrX[0];
				_transformedNormalisedGradientSource(i, j, k, 1) = ptrY[0];
				_transformedNormalisedGradientSource(i, j, k, 2) = ptrZ[0];
			  } else {
				cerr << "irktImageGradientFreeFormRegistration2::UpdateSource(): Interpolation not implemented" << endl;
				exit(1);
			  }
            } else {
              _transformedNormalisedGradientSource(i, j, k, 0) = -2;
              _transformedNormalisedGradientSource(i, j, k, 1) = -2;
              _transformedNormalisedGradientSource(i, j, k, 2) = -2;
            }
          } else {
        	_transformedNormalisedGradientSource(i, j, k, 0) = -2;
        	_transformedNormalisedGradientSource(i, j, k, 1) = -2;
        	_transformedNormalisedGradientSource(i, j, k, 2) = -2;
          }
          ptr2disp += 3;
		  ptr2latt += 3;
        }
      }
    }
  }

  delete vector1;
  delete vector2;
  delete result1;
  delete result2;
  delete result3;

  IRTK_END_TIMING();
}

void irtkImageGradientFreeFormRegistration2::UpdateSourceAndGradient()
{
  double *ptrX;
  double *ptrY;
  double *ptrZ;
  double *ptrXX;
  double *ptrXY;
  double *ptrXZ;
  double *ptrYX;
  double *ptrYY;
  double *ptrYZ;
  double *ptrZX;
  double *ptrZY;
  double *ptrZZ;
  double *vector1;
  double *vector2;
  double *result1;
  double *result2;
  double *result3;
  double *vector1_grad;
  double *vector2_grad;
  double *result1_gradX;
  double *result2_gradX;
  double *result3_gradX;
  double *result1_gradY;
  double *result2_gradY;
  double *result3_gradY;
  double *result1_gradZ;
  double *result2_gradZ;
  double *result3_gradZ;
  double x, y, z, t1, t2, u1, u2, v1, v2;
  int a, b, c, i, j, k, offset1, offset2, offset3, offset4, offset5, offset6, offset7, offset8;

  IRTK_START_TIMING();

  // Print debugging information
  this->Debug("irktImageGradientFreeFormRegistration2::UpdateSourceAndGradient()");

  // Generate transformed tmp image
  _transformedNormalisedGradientSource = _normalisedGradientTarget;

  // Calculate offsets for fast pixel access
  offset1 = 0;
  offset2 = 1;
  offset3 = this->_source->GetX();
  offset4 = this->_source->GetX()+1;
  offset5 = this->_source->GetX()*this->_source->GetY();
  offset6 = this->_source->GetX()*this->_source->GetY()+1;
  offset7 = this->_source->GetX()*this->_source->GetY()+this->_source->GetX();
  offset8 = this->_source->GetX()*this->_source->GetY()+this->_source->GetX()+1;

  vector1 = new double[3];
  vector2 = new double[3];
  result1 = new double[3];
  result2 = new double[3];
  result3 = new double[3];
  vector1_grad = new double[3];
  vector2_grad = new double[3];
  result1_gradX = new double[3];
  result2_gradX = new double[3];
  result3_gradX = new double[3];
  result1_gradY = new double[3];
  result2_gradY = new double[3];
  result3_gradY = new double[3];
  result1_gradZ = new double[3];
  result2_gradZ = new double[3];
  result3_gradZ = new double[3];

  double *ptr2disp = _displacementLUT;
  double *ptr2latt = _latticeCoordLUT;
  if ((_target->GetZ() == 1) && (_source->GetZ() == 1)) {
    for (j = 0; j < _target->GetY(); j++) {
      for (i = 0; i < _target->GetX(); i++) {
        if (_distanceMask.Get(i, j, 0) == 0) {
          x = ptr2latt[0];
		  y = ptr2latt[1];
		  z = 0;
		  _affd->FFD2D(x, y);
		  x += ptr2disp[0];
		  y += ptr2disp[1];
		  z += ptr2disp[2];
		  _source->WorldToImage(x, y, z);

          // Check whether transformed point is inside volume
          if ((x > 0) && (x < _source->GetX()-1) &&
              (y > 0) && (y < _source->GetY()-1)) {

			if (_InterpolationMode == Interpolation_Linear) {
              // Calculated integer coordinates
              a  = int(x);
              b  = int(y);

              // Calculated fractional coordinates
              t1 = x - a;
              u1 = y - b;
              t2 = 1 - t1;
              u2 = 1 - u1;

              // Spherical linear interpolation in source and gradient images
			  ptrX = _normalisedGradientSource.GetPointerToVoxels(a, b, 0, 0);
			  ptrY = _normalisedGradientSource.GetPointerToVoxels(a, b, 0, 1);
			  ptrXX = _normalisedGradientSourceGradient[0].GetPointerToVoxels(a, b, 0, 0);
			  ptrXY = _normalisedGradientSourceGradient[0].GetPointerToVoxels(a, b, 0, 1);
			  ptrYX = _normalisedGradientSourceGradient[1].GetPointerToVoxels(a, b, 0, 0);
			  ptrYY = _normalisedGradientSourceGradient[1].GetPointerToVoxels(a, b, 0, 1);

			  // ************************************************************ //

			  vector1[0] = ptrX[offset2];
			  vector1[1] = ptrY[offset2];
			  vector1[2] = 0;

			  vector2[0] = ptrX[offset4];
			  vector2[1] = ptrY[offset4];
			  vector2[2] = 0;

			  this->Slerp(result1, vector1, vector2, u2, u1);

			  vector1_grad[0] = ptrXX[offset2];
			  vector1_grad[1] = ptrYX[offset2];
			  vector1_grad[2] = 0;

			  vector2_grad[0] = ptrXX[offset4];
			  vector2_grad[1] = ptrYX[offset4];
			  vector2_grad[2] = 0;

			  this->SlerpGradient(result1_gradX, vector1_grad, vector2_grad, vector1, vector2, u2, u1);

			  vector1_grad[0] = ptrXY[offset2];
			  vector1_grad[1] = ptrYY[offset2];
			  vector1_grad[2] = 0;

			  vector2_grad[0] = ptrXY[offset4];
			  vector2_grad[1] = ptrYY[offset4];
			  vector2_grad[2] = 0;

			  this->SlerpGradient(result1_gradY, vector1_grad, vector2_grad, vector1, vector2, u2, u1);

			  // ************************************************************ //

			  vector1[0] = ptrX[offset1];
			  vector1[1] = ptrY[offset1];
			  vector1[2] = 0;

			  vector2[0] = ptrX[offset3];
			  vector2[1] = ptrY[offset3];
			  vector2[2] = 0;

			  this->Slerp(result2, vector1, vector2, u2, u1);

			  vector1_grad[0] = ptrXX[offset1];
			  vector1_grad[1] = ptrYX[offset1];
			  vector1_grad[2] = 0;

			  vector2_grad[0] = ptrXX[offset3];
			  vector2_grad[1] = ptrYX[offset3];
			  vector2_grad[2] = 0;

			  this->SlerpGradient(result2_gradX, vector1_grad, vector2_grad, vector1, vector2, u2, u1);

			  vector1_grad[0] = ptrXY[offset1];
			  vector1_grad[1] = ptrYY[offset1];
			  vector1_grad[2] = 0;

			  vector2_grad[0] = ptrXY[offset3];
			  vector2_grad[1] = ptrYY[offset3];
			  vector2_grad[2] = 0;

			  this->SlerpGradient(result2_gradY, vector1_grad, vector2_grad, vector1, vector2, u2, u1);

			  // ************************************************************ //

			  // In the mixture, slerp gradients first so we have the values of result1 and result2
			  this->SlerpGradient(result1_gradX, result1_gradX, result2_gradX, result1, result2, t1, t2);
			  this->SlerpGradient(result1_gradY, result1_gradY, result2_gradY, result1, result2, t1, t2);
			  this->Slerp(result1, result1, result2, t1, t2);

			  _transformedNormalisedGradientSource(i, j, 0, 0) = result1[0];
			  _transformedNormalisedGradientSource(i, j, 0, 1) = result1[1];

			  _transformedNormalisedGradientSourceGradient[0](i, j, 0, 0) = result1_gradX[0];
			  _transformedNormalisedGradientSourceGradient[0](i, j, 0, 1) = result1_gradY[0];

			  _transformedNormalisedGradientSourceGradient[1](i, j, 0, 0) = result1_gradX[1];
			  _transformedNormalisedGradientSourceGradient[1](i, j, 0, 1) = result1_gradY[1];

			} else if (_InterpolationMode == Interpolation_NN) {
			  // Calculated NN integer coordinates
			  a  = round(x);
			  b  = round(y);

			  // NN interpolation in source image
			  ptrX = _normalisedGradientSource.GetPointerToVoxels(a, b, 0, 0);
			  ptrY = _normalisedGradientSource.GetPointerToVoxels(a, b, 0, 1);
			  _transformedNormalisedGradientSource(i, j, 0, 0) = ptrX[0];
			  _transformedNormalisedGradientSource(i, j, 0, 1) = ptrY[0];

			  // NN interpolation in gradient image
			  ptrXX = _normalisedGradientSourceGradient[0].GetPointerToVoxels(a, b, 0, 0);
			  _transformedNormalisedGradientSourceGradient[0](i, j, 0, 0) = ptrXX[0];
			  ptrXY = _normalisedGradientSourceGradient[0].GetPointerToVoxels(a, b, 0, 1);
			  _transformedNormalisedGradientSourceGradient[0](i, j, 0, 1) = ptrXY[0];
			  ptrYX = _normalisedGradientSourceGradient[1].GetPointerToVoxels(a, b, 0, 0);
			  _transformedNormalisedGradientSourceGradient[1](i, j, 0, 0) = ptrYX[0];
			  ptrYY = _normalisedGradientSourceGradient[1].GetPointerToVoxels(a, b, 0, 1);
			  _transformedNormalisedGradientSourceGradient[1](i, j, 0, 1) = ptrYY[0];
			} else {
			  cerr << "irktImageGradientFreeFormRegistration2::UpdateSource(): Interpolation not implemented" << endl;
			  exit(1);
			}
          } else {
        	_transformedNormalisedGradientSource(i, j, 0, 0) = -2;
        	_transformedNormalisedGradientSource(i, j, 0, 1) = -2;
        	_transformedNormalisedGradientSourceGradient[0](i, j, 0, 0) = 0;
        	_transformedNormalisedGradientSourceGradient[0](i, j, 0, 1) = 0;
        	_transformedNormalisedGradientSourceGradient[1](i, j, 0, 0) = 0;
        	_transformedNormalisedGradientSourceGradient[1](i, j, 0, 1) = 0;
          }
        } else {
          _transformedNormalisedGradientSource(i, j, 0, 0) = -2;
		  _transformedNormalisedGradientSource(i, j, 0, 1) = -2;
		  _transformedNormalisedGradientSourceGradient[0](i, j, 0, 0) = 0;
		  _transformedNormalisedGradientSourceGradient[0](i, j, 0, 1) = 0;
		  _transformedNormalisedGradientSourceGradient[1](i, j, 0, 0) = 0;
		  _transformedNormalisedGradientSourceGradient[1](i, j, 0, 1) = 0;
        }
        ptr2disp += 3;
		ptr2latt += 3;
      }
    }
  } else {
    for (k = 0; k < _target->GetZ(); k++) {
      for (j = 0; j < _target->GetY(); j++) {
        for (i = 0; i < _target->GetX(); i++) {
          if (_distanceMask.Get(i, j, k) == 0) {
        	x = ptr2latt[0];
			y = ptr2latt[1];
			z = ptr2latt[2];
			_affd->FFD3D(x, y, z);
			x += ptr2disp[0];
			y += ptr2disp[1];
			z += ptr2disp[2];
			_source->WorldToImage(x, y, z);

            // Check whether transformed point is inside volume
            if ((x > 0) && (x < _source->GetX()-1) &&
                (y > 0) && (y < _source->GetY()-1) &&
                (z > 0) && (z < _source->GetZ()-1)) {

			  if (_InterpolationMode == Interpolation_Linear) {
                // Calculated integer coordinates
                a  = int(x);
                b  = int(y);
                c  = int(z);

                // Calculated fractional coordinates
                t1 = x - a;
                u1 = y - b;
                v1 = z - c;
                t2 = 1 - t1;
                u2 = 1 - u1;
                v2 = 1 - v1;

                // Spherical linear interpolation in source image
			    ptrX = _normalisedGradientSource.GetPointerToVoxels(a, b, c, 0);
			    ptrY = _normalisedGradientSource.GetPointerToVoxels(a, b, c, 1);
			    ptrZ = _normalisedGradientSource.GetPointerToVoxels(a, b, c, 2);
			    ptrXX = _normalisedGradientSourceGradient[0].GetPointerToVoxels(a, b, c, 0);
				ptrXY = _normalisedGradientSourceGradient[0].GetPointerToVoxels(a, b, c, 1);
				ptrXZ = _normalisedGradientSourceGradient[0].GetPointerToVoxels(a, b, c, 2);
				ptrYX = _normalisedGradientSourceGradient[1].GetPointerToVoxels(a, b, c, 0);
				ptrYY = _normalisedGradientSourceGradient[1].GetPointerToVoxels(a, b, c, 1);
				ptrYZ = _normalisedGradientSourceGradient[1].GetPointerToVoxels(a, b, c, 2);
				ptrZX = _normalisedGradientSourceGradient[2].GetPointerToVoxels(a, b, c, 0);
				ptrZY = _normalisedGradientSourceGradient[2].GetPointerToVoxels(a, b, c, 1);
				ptrZZ = _normalisedGradientSourceGradient[2].GetPointerToVoxels(a, b, c, 2);

			    // ************************************************************ //

				vector1[0] = ptrX[offset2];
				vector1[1] = ptrY[offset2];
				vector1[2] = ptrZ[offset2];

				vector2[0] = ptrX[offset6];
				vector2[1] = ptrY[offset6];
				vector2[2] = ptrZ[offset6];

				this->Slerp(result1, vector1, vector2, v2, v1);

				vector1_grad[0] = ptrXX[offset2];
				vector1_grad[1] = ptrYX[offset2];
				vector1_grad[2] = ptrZX[offset2];

				vector2_grad[0] = ptrXX[offset6];
				vector2_grad[1] = ptrYX[offset6];
				vector2_grad[2] = ptrZX[offset6];

				this->SlerpGradient(result1_gradX, vector1_grad, vector2_grad, vector1, vector2, v2, v1);

				vector1_grad[0] = ptrXY[offset2];
				vector1_grad[1] = ptrYY[offset2];
				vector1_grad[2] = ptrZY[offset2];

				vector2_grad[0] = ptrXY[offset6];
				vector2_grad[1] = ptrYY[offset6];
				vector2_grad[2] = ptrZY[offset6];

				this->SlerpGradient(result1_gradY, vector1_grad, vector2_grad, vector1, vector2, v2, v1);

				vector1_grad[0] = ptrXZ[offset2];
				vector1_grad[1] = ptrYZ[offset2];
				vector1_grad[2] = ptrZZ[offset2];

				vector2_grad[0] = ptrXZ[offset6];
				vector2_grad[1] = ptrYZ[offset6];
				vector2_grad[2] = ptrZZ[offset6];

				this->SlerpGradient(result1_gradZ, vector1_grad, vector2_grad, vector1, vector2, v2, v1);

				// ************************************************************ //

				vector1[0] = ptrX[offset4];
				vector1[1] = ptrY[offset4];
				vector1[2] = ptrZ[offset4];

				vector2[0] = ptrX[offset8];
				vector2[1] = ptrY[offset8];
				vector2[2] = ptrZ[offset8];

				this->Slerp(result2, vector1, vector2, v2, v1);

				vector1_grad[0] = ptrXX[offset4];
				vector1_grad[1] = ptrYX[offset4];
				vector1_grad[2] = ptrZX[offset4];

				vector2_grad[0] = ptrXX[offset8];
				vector2_grad[1] = ptrYX[offset8];
				vector2_grad[2] = ptrZX[offset8];

				this->SlerpGradient(result2_gradX, vector1_grad, vector2_grad, vector1, vector2, v2, v1);

				vector1_grad[0] = ptrXY[offset4];
				vector1_grad[1] = ptrYY[offset4];
				vector1_grad[2] = ptrZY[offset4];

				vector2_grad[0] = ptrXY[offset8];
				vector2_grad[1] = ptrYY[offset8];
				vector2_grad[2] = ptrZY[offset8];

				this->SlerpGradient(result2_gradY, vector1_grad, vector2_grad, vector1, vector2, v2, v1);

				vector1_grad[0] = ptrXZ[offset4];
				vector1_grad[1] = ptrYZ[offset4];
				vector1_grad[2] = ptrZZ[offset4];

				vector2_grad[0] = ptrXZ[offset8];
				vector2_grad[1] = ptrYZ[offset8];
				vector2_grad[2] = ptrZZ[offset8];

				this->SlerpGradient(result2_gradZ, vector1_grad, vector2_grad, vector1, vector2, v2, v1);

				// ************************************************************ //

				// In the mixture, slerp gradients first so we have the values of result1 and result2
			    this->SlerpGradient(result1_gradX, result1_gradX, result2_gradX, result1, result2, u2, u1);
			    this->SlerpGradient(result1_gradY, result1_gradY, result2_gradY, result1, result2, u2, u1);
			    this->SlerpGradient(result1_gradZ, result1_gradZ, result2_gradZ, result1, result2, u2, u1);
				this->Slerp(result1, result1, result2, u2, u1);

				// ************************************************************ //

				vector1[0] = ptrX[offset1];
				vector1[1] = ptrY[offset1];
				vector1[2] = ptrZ[offset1];

				vector2[0] = ptrX[offset5];
				vector2[1] = ptrY[offset5];
				vector2[2] = ptrZ[offset5];

				this->Slerp(result2, vector1, vector2, v2, v1);

				vector1_grad[0] = ptrXX[offset1];
				vector1_grad[1] = ptrYX[offset1];
				vector1_grad[2] = ptrZX[offset1];

				vector2_grad[0] = ptrXX[offset5];
				vector2_grad[1] = ptrYX[offset5];
				vector2_grad[2] = ptrZX[offset5];

				this->SlerpGradient(result2_gradX, vector1_grad, vector2_grad, vector1, vector2, v2, v1);

				vector1_grad[0] = ptrXY[offset1];
				vector1_grad[1] = ptrYY[offset1];
				vector1_grad[2] = ptrZY[offset1];

				vector2_grad[0] = ptrXY[offset5];
				vector2_grad[1] = ptrYY[offset5];
				vector2_grad[2] = ptrZY[offset5];

				this->SlerpGradient(result2_gradY, vector1_grad, vector2_grad, vector1, vector2, v2, v1);

				vector1_grad[0] = ptrXZ[offset1];
				vector1_grad[1] = ptrYZ[offset1];
				vector1_grad[2] = ptrZZ[offset1];

				vector2_grad[0] = ptrXZ[offset5];
				vector2_grad[1] = ptrYZ[offset5];
				vector2_grad[2] = ptrZZ[offset5];

				this->SlerpGradient(result2_gradZ, vector1_grad, vector2_grad, vector1, vector2, v2, v1);

				// ************************************************************ //

				vector1[0] = ptrX[offset3];
				vector1[1] = ptrY[offset3];
				vector1[2] = ptrZ[offset3];

				vector2[0] = ptrX[offset7];
				vector2[1] = ptrY[offset7];
				vector2[2] = ptrZ[offset7];

				this->Slerp(result3, vector1, vector2, v2, v1);

				vector1_grad[0] = ptrXX[offset3];
				vector1_grad[1] = ptrYX[offset3];
				vector1_grad[2] = ptrZX[offset3];

				vector2_grad[0] = ptrXX[offset7];
				vector2_grad[1] = ptrYX[offset7];
				vector2_grad[2] = ptrZX[offset7];

				this->SlerpGradient(result3_gradX, vector1_grad, vector2_grad, vector1, vector2, v2, v1);

				vector1_grad[0] = ptrXY[offset3];
				vector1_grad[1] = ptrYY[offset3];
				vector1_grad[2] = ptrZY[offset3];

				vector2_grad[0] = ptrXY[offset7];
				vector2_grad[1] = ptrYY[offset7];
				vector2_grad[2] = ptrZY[offset7];

				this->SlerpGradient(result3_gradY, vector1_grad, vector2_grad, vector1, vector2, v2, v1);

				vector1_grad[0] = ptrXZ[offset3];
				vector1_grad[1] = ptrYZ[offset3];
				vector1_grad[2] = ptrZZ[offset3];

				vector2_grad[0] = ptrXZ[offset7];
				vector2_grad[1] = ptrYZ[offset7];
				vector2_grad[2] = ptrZZ[offset7];

				this->SlerpGradient(result3_gradZ, vector1_grad, vector2_grad, vector1, vector2, v2, v1);

				// ************************************************************ //

				// In the mixture, slerp gradients first so we have the values of result2 and result3
				this->SlerpGradient(result2_gradX, result2_gradX, result3_gradX, result2, result3, u2, u1);
				this->SlerpGradient(result2_gradY, result2_gradY, result3_gradY, result2, result3, u2, u1);
				this->SlerpGradient(result2_gradZ, result2_gradZ, result3_gradZ, result2, result3, u2, u1);
				this->Slerp(result2, result2, result3, u2, u1);

				// ************************************************************ //

				// In the mixture, slerp gradients first so we have the values of result1 and result2
				this->SlerpGradient(result1_gradX, result1_gradX, result2_gradX, result1, result2, t1, t2);
				this->SlerpGradient(result1_gradY, result1_gradY, result2_gradY, result1, result2, t1, t2);
				this->SlerpGradient(result1_gradZ, result1_gradZ, result2_gradZ, result1, result2, t1, t2);
				this->Slerp(result1, result1, result2, t1, t2);

				_transformedNormalisedGradientSource(i, j, k, 0) = result1[0];
				_transformedNormalisedGradientSource(i, j, k, 1) = result1[1];
				_transformedNormalisedGradientSource(i, j, k, 2) = result1[2];

				_transformedNormalisedGradientSourceGradient[0](i, j, k, 0) = result1_gradX[0];
				_transformedNormalisedGradientSourceGradient[0](i, j, k, 1) = result1_gradY[0];
				_transformedNormalisedGradientSourceGradient[0](i, j, k, 2) = result1_gradZ[0];

				_transformedNormalisedGradientSourceGradient[1](i, j, k, 0) = result1_gradX[1];
				_transformedNormalisedGradientSourceGradient[1](i, j, k, 1) = result1_gradY[1];
				_transformedNormalisedGradientSourceGradient[1](i, j, k, 2) = result1_gradZ[1];

				_transformedNormalisedGradientSourceGradient[2](i, j, k, 0) = result2_gradX[2];
				_transformedNormalisedGradientSourceGradient[2](i, j, k, 1) = result2_gradY[2];
				_transformedNormalisedGradientSourceGradient[2](i, j, k, 2) = result2_gradZ[2];

			  } else if (_InterpolationMode == Interpolation_NN) {
				// Calculated NN integer coordinates
				a  = round(x);
				b  = round(y);
				c  = round(z);

				// NN interpolation in source image
				ptrX = _normalisedGradientSource.GetPointerToVoxels(a, b, c, 0);
				ptrY = _normalisedGradientSource.GetPointerToVoxels(a, b, c, 1);
				ptrZ = _normalisedGradientSource.GetPointerToVoxels(a, b, c, 2);
				_transformedNormalisedGradientSource(i, j, k, 0) = ptrX[0];
				_transformedNormalisedGradientSource(i, j, k, 1) = ptrY[0];
				_transformedNormalisedGradientSource(i, j, k, 2) = ptrZ[0];

				ptrXX = _normalisedGradientSourceGradient[0].GetPointerToVoxels(a, b, c, 0);
				_transformedNormalisedGradientSourceGradient[0](i, j, k, 0) = ptrXX[0];
				ptrXY = _normalisedGradientSourceGradient[0].GetPointerToVoxels(a, b, c, 1);
				_transformedNormalisedGradientSourceGradient[0](i, j, k, 1) = ptrXY[0];
				ptrXZ = _normalisedGradientSourceGradient[0].GetPointerToVoxels(a, b, c, 2);
				_transformedNormalisedGradientSourceGradient[0](i, j, k, 2) = ptrXZ[0];
				ptrYX = _normalisedGradientSourceGradient[1].GetPointerToVoxels(a, b, c, 0);
				_transformedNormalisedGradientSourceGradient[1](i, j, k, 0) = ptrYX[0];
				ptrYY = _normalisedGradientSourceGradient[1].GetPointerToVoxels(a, b, c, 1);
				_transformedNormalisedGradientSourceGradient[1](i, j, k, 1) = ptrYY[0];
				ptrYZ = _normalisedGradientSourceGradient[1].GetPointerToVoxels(a, b, c, 2);
				_transformedNormalisedGradientSourceGradient[1](i, j, k, 2) = ptrYZ[0];
				ptrZX = _normalisedGradientSourceGradient[2].GetPointerToVoxels(a, b, c, 0);
				_transformedNormalisedGradientSourceGradient[2](i, j, k, 0) = ptrZX[0];
				ptrZY = _normalisedGradientSourceGradient[2].GetPointerToVoxels(a, b, c, 1);
				_transformedNormalisedGradientSourceGradient[2](i, j, k, 1) = ptrZY[0];
				ptrZZ = _normalisedGradientSourceGradient[2].GetPointerToVoxels(a, b, c, 2);
				_transformedNormalisedGradientSourceGradient[2](i, j, k, 2) = ptrZZ[0];
			  } else {
				cerr << "irktImageGradientFreeFormRegistration2::UpdateSource(): Interpolation not implemented" << endl;
				exit(1);
			  }
            } else {
			  _transformedNormalisedGradientSource(i, j, k, 0) = -2;
			  _transformedNormalisedGradientSource(i, j, k, 1) = -2;
			  _transformedNormalisedGradientSource(i, j, k, 2) = -2;
			  _transformedNormalisedGradientSourceGradient[0](i, j, k, 0) = 0;
			  _transformedNormalisedGradientSourceGradient[0](i, j, k, 1) = 0;
			  _transformedNormalisedGradientSourceGradient[0](i, j, k, 2) = 0;
			  _transformedNormalisedGradientSourceGradient[1](i, j, k, 0) = 0;
			  _transformedNormalisedGradientSourceGradient[1](i, j, k, 1) = 0;
			  _transformedNormalisedGradientSourceGradient[1](i, j, k, 2) = 0;
			  _transformedNormalisedGradientSourceGradient[2](i, j, k, 0) = 0;
			  _transformedNormalisedGradientSourceGradient[2](i, j, k, 1) = 0;
			  _transformedNormalisedGradientSourceGradient[2](i, j, k, 2) = 0;
		    }
		  } else {
		    _transformedNormalisedGradientSource(i, j, k, 0) = -2;
		    _transformedNormalisedGradientSource(i, j, k, 1) = -2;
		    _transformedNormalisedGradientSource(i, j, k, 2) = -2;
		    _transformedNormalisedGradientSourceGradient[0](i, j, k, 0) = 0;
			_transformedNormalisedGradientSourceGradient[0](i, j, k, 1) = 0;
			_transformedNormalisedGradientSourceGradient[0](i, j, k, 2) = 0;
			_transformedNormalisedGradientSourceGradient[1](i, j, k, 0) = 0;
			_transformedNormalisedGradientSourceGradient[1](i, j, k, 1) = 0;
			_transformedNormalisedGradientSourceGradient[1](i, j, k, 2) = 0;
			_transformedNormalisedGradientSourceGradient[2](i, j, k, 0) = 0;
			_transformedNormalisedGradientSourceGradient[2](i, j, k, 1) = 0;
			_transformedNormalisedGradientSourceGradient[2](i, j, k, 2) = 0;
		  }
          ptr2disp += 3;
		  ptr2latt += 3;
        }
      }
    }
  }

  delete vector1;
  delete vector2;
  delete result1;
  delete result2;
  delete result3;
  delete vector1_grad;
  delete vector2_grad;
  delete result1_gradX;
  delete result2_gradX;
  delete result3_gradX;
  delete result1_gradY;
  delete result2_gradY;
  delete result3_gradY;
  delete result1_gradZ;
  delete result2_gradZ;
  delete result3_gradZ;

  IRTK_END_TIMING();
}

void irtkImageGradientFreeFormRegistration2::EvaluateGradient2D(double *gradient)
{
    double basis, pos[3], offset;
    int i, j, i1, i2, j1, j2, k1, k2, x, y, index, index2;

    // Initialize gradient to zero
    for (i = 0; i < _affd->NumberOfDOFs(); i++) {
        gradient[i] = 0;
    }

    // Loop over control points
    for (y = 0; y < _affd->GetY(); y++) {
        for (x = 0; x < _affd->GetX(); x++) {

            // Compute DoFs corresponding to the control point
            index  = _affd->LatticeToIndex(x, y, 0);
            index2 = index+_affd->GetX()*_affd->GetY()*_affd->GetZ();

            // Check if any DoF corresponding to the control point is active
            if ((_affd->irtkTransformation::GetStatus(index) == _Active) || (_affd->irtkTransformation::GetStatus(index2) == _Active)) {

                // If so, calculate bounding box of control point in image coordinates
                _affd->BoundingBoxImage(_target, index, i1, j1, k1, i2, j2, k2, 1.0);

                // Loop over all voxels in the target (reference) volume
                for (j = j1; j <= j2; j++) {
                    for (i = i1; i <= i2; i++) {

                        // Check whether reference point is valid
                    	if ((_distanceMask(i, j, 0) == 0) && (_transformedNormalisedGradientSource(i, j, 0, 0) >= -1)) {

                            // Convert position from voxel coordinates to world coordinates
                            pos[0] = i;
                            pos[1] = j;
                            pos[2] = 0;
                            _target->ImageToWorld(pos[0], pos[1], pos[2]);

                            // Convert world coordinates into lattice coordinates
                            _affd->WorldToLattice(pos[0], pos[1], pos[2]);

                            // Compute B-spline tensor product at pos
                            basis = _affd->B(pos[0] - x) * _affd->B(pos[1] - y);

                            // Convert voxel-based gradient into gradient with respect to parameters (chain rule)
                            //
                            // NOTE: This currently assumes that the control points displacements are aligned with the world coordinate displacements
                            //
                            gradient[index]  += basis * _similarityGradient(i, j, 0, 0);
                            gradient[index2] += basis * _similarityGradient(i, j, 0, 1);
                        }
                    }
                }
            }
        }
    }
}

void irtkImageGradientFreeFormRegistration2::EvaluateGradient3D(double *gradient)
{
    double basis, *ptr;
    int i, j, k, i1, i2, j1, j2, k1, k2, x, y, z, index, index2, index3;

    // Initialize gradient to zero
    for (i = 0; i < _affd->NumberOfDOFs(); i++) {
        gradient[i] = 0;
    }

    // Loop over control points
    for (z = 0; z < _affd->GetZ(); z++) {
        for (y = 0; y < _affd->GetY(); y++) {
            for (x = 0; x < _affd->GetX(); x++) {

                // Compute DoFs corresponding to the control point
                index  = _affd->LatticeToIndex(x, y, z);
                index2 = index+_affd->GetX()*_affd->GetY()*_affd->GetZ();
                index3 = index+2*_affd->GetX()*_affd->GetY()*_affd->GetZ();

                // Check if any DoF corresponding to the control point is active
                if ((_affd->irtkTransformation::GetStatus(index) == _Active) || (_affd->irtkTransformation::GetStatus(index2) == _Active) || (_affd->irtkTransformation::GetStatus(index3) == _Active)) {

                    // If so, calculate bounding box of control point in image coordinates
                    _affd->BoundingBoxImage(_target, index, i1, j1, k1, i2, j2, k2, 1);

                    // Loop over all voxels in the target (reference) volume
                    //
                    // NOTE: This currently assumes that the control point lattice is aligned with the target image
                    //
                    for (k = k1; k <= k2; k++) {
                        for (j = j1; j <= j2; j++) {
                            ptr = &(_latticeCoordLUT[3 * (k * (_target->GetX()*_target->GetY()) + j * _target->GetX() + i1)]);
                            for (i = i1; i <= i2; i++) {
                                // Check whether reference point is valid
                            	if ((_distanceMask(i, j, k) == 0) && (_transformedNormalisedGradientSource(i, j, k, 0) >= -1)) {
                                    // Compute B-spline tensor product at current position
                                    basis = _affd->B(ptr[0] - x) * _affd->B(ptr[1] - y) * _affd->B(ptr[2] - z);

                                    // Convert voxel-based gradient into gradient with respect to parameters (chain rule)
                                    //
                                    // NOTE: This currently assumes that the control points displacements are aligned with the world coordinate displacements
                                    //
                                    gradient[index]  += basis * _similarityGradient(i, j, k, 0);
                                    gradient[index2] += basis * _similarityGradient(i, j, k, 1);
                                    gradient[index3] += basis * _similarityGradient(i, j, k, 2);
                                }
                                ptr += 3;
                            }
                        }
                    }
                }
            }
        }
    }
}

double irtkImageGradientFreeFormRegistration2::Evaluate()
{
  double tmp, similarity;

  // Print debugging information
  this->Debug("irtkImageGradientFreeFormRegistration2::Evaluate");

  // Evaluate similarity
  switch (_SimilarityMeasure) {
    case NGD:
   	  similarity = this->EvaluateNGD();
      break;
    case NGP:
   	  similarity = this->EvaluateNGP();
      break;
    case NGS:
      similarity = this->EvaluateNGS();
      break;
    default:
   	  similarity = 0;
      cerr << this->NameOfClass() << "::Evaluate: No such metric implemented" << endl;
      exit(1);
  }
  cout << "Similarity = " << similarity << "\t";

  // Add penalty for smoothness
  if (this->_Lambda1 > 0) {
    tmp = this->_Lambda1*this->SmoothnessPenalty();
    cout << "Bending = " << tmp << "\t";
    similarity += tmp;
  }
  // Add penalty for volume preservation
  if (this->_Lambda2 > 0) {
    tmp = this->_Lambda2*this->VolumePreservationPenalty();
    cout << "Volume = " << tmp;
    similarity += tmp;
  }
  if ((this->_Lambda1 > 0) || (this->_Lambda2 > 0)) cout << endl;

  //Return similarity measure + penalty terms
  return similarity;
}

double irtkImageGradientFreeFormRegistration2::EvaluateGradient(double *gradient)
{
  double norm, max_length, mx, my, mz;
  int i, j, k, x, y, z, index, index2, index3;
  static double *g = NULL, *h = NULL, gg, dgg, gamma;

  IRTK_START_TIMING();

  // Compute gradient with respect to displacements
  // Allocate memory for metric
  switch (_SimilarityMeasure) {
    case NGD:
      this->EvaluateGradientNGD();
      break;
    case NGP:
      this->EvaluateGradientNGP();
      break;
    case NGS:
	  this->EvaluateGradientNGS();
	  break;
    default:
      cerr << this->NameOfClass() << "::Evaluate: No such metric implemented" << endl;
      exit(1);
  }

  // Extract matrix for reorientation of gradient
  irtkMatrix m = _source->GetImageToWorldMatrix();

  // Reorient gradient
  for (k = 0; k < _target->GetZ(); k++) {
    for (j = 0; j < _target->GetY(); j++) {
      for (i = 0; i < _target->GetX(); i++) {
    	if ((_distanceMask.Get(i, j, k) == 0) && (_transformedNormalisedGradientSource(i, j, k, 0) >= -1)) {
          mx = m(0, 0) * _similarityGradient(i, j, k, 0) + m(0, 1) * _similarityGradient(i, j, k, 1) + m(0, 2) * _similarityGradient(i, j, k, 2);
          my = m(1, 0) * _similarityGradient(i, j, k, 0) + m(1, 1) * _similarityGradient(i, j, k, 1) + m(1, 2) * _similarityGradient(i, j, k, 2);
          mz = m(2, 0) * _similarityGradient(i, j, k, 0) + m(2, 1) * _similarityGradient(i, j, k, 1) + m(2, 2) * _similarityGradient(i, j, k, 2);
          _similarityGradient(i, j, k, 0) = mx;
          _similarityGradient(i, j, k, 1) = my;
          _similarityGradient(i, j, k, 2) = mz;
        }
      }
    }
  }

  if (_affd->GetZ() == 1) {
    this->EvaluateGradient2D(gradient);
  } else {
    this->EvaluateGradient3D(gradient);
  }

  // Update gradient
  if (_CurrentIteration == 0) {
    // First iteration, so let's initialize
    if (g != NULL) delete []g;
    g = new double [_affd->NumberOfDOFs()];
    if (h != NULL) delete []h;
    h = new double [_affd->NumberOfDOFs()];
    for (i = 0; i < _affd->NumberOfDOFs(); i++) {
      g[i] = -gradient[i];
      h[i] = g[i];
    }
  } else {
    // Update gradient direction to be conjugate
    gg = 0;
    dgg = 0;
    for (i = 0; i < _affd->NumberOfDOFs(); i++) {
      gg  += g[i]*h[i];
      dgg += (gradient[i]+g[i])*gradient[i];
    }
    gamma = dgg/gg;
    for (i = 0; i < _affd->NumberOfDOFs(); i++) {
      g[i] = -gradient[i];
      h[i] = g[i] + gamma*h[i];
      gradient[i] = -h[i];
    }
  }

  if (this->_Lambda1 > 0) {
    this->SmoothnessPenaltyGradient(gradient);
  }

  if (this->_Lambda2 > 0) {
    this->VolumePreservationPenaltyGradient(gradient);
  }

  // Calculate maximum of gradient vector
  max_length = 0;
  for (z = 0; z < _affd->GetZ(); z++) {
    for (y = 0; y < _affd->GetY(); y++) {
      for (x = 0; x < _affd->GetX(); x++) {
        index  = _affd->LatticeToIndex(x, y, z);
        index2 = index+_affd->GetX()*_affd->GetY()*_affd->GetZ();
        index3 = index+2*_affd->GetX()*_affd->GetY()*_affd->GetZ();
        norm = sqrt(gradient[index] * gradient[index] + gradient[index2] * gradient[index2] + gradient[index3] * gradient[index3]);
        if (norm > max_length) max_length = norm;
      }
    }
  }

  // Deal with active and passive control points
  for (i = 0; i < _affd->NumberOfDOFs(); i++) {
	if (_affd->irtkTransformation::GetStatus(i) == _Passive) {
	  gradient[i] = 0;
	  _affd->Put(i,0);
	}
  }

  IRTK_END_TIMING();

  return max_length;
}

double irtkImageGradientFreeFormRegistration2::EvaluateNGD()
{
  int i, j, k, n;
  double sum;

  // Print debugging information
  this->Debug("irtkImageGradientFreeFormRegistration2::EvaluateNGD");

  // Initialize metric
  n = 0;
  sum = 0;

  // Compute metric
  for (k = 0; k < _normalisedGradientTarget.GetZ(); k++) {
    for (j = 0; j < _normalisedGradientTarget.GetY(); j++) {
      for (i = 0; i < _normalisedGradientTarget.GetX(); i++) {
    	if ((_distanceMask(i, j, k) == 0) && (_transformedNormalisedGradientSource(i, j, k, 0) >= -1)) {
          sum += (_normalisedGradientTarget(i, j, k, 0) - _transformedNormalisedGradientSource(i, j, k, 0)) * (_normalisedGradientTarget(i, j, k, 0) - _transformedNormalisedGradientSource(i, j, k, 0)) +
        		 (_normalisedGradientTarget(i, j, k, 1) - _transformedNormalisedGradientSource(i, j, k, 1)) * (_normalisedGradientTarget(i, j, k, 1) - _transformedNormalisedGradientSource(i, j, k, 1)) +
        		 (_normalisedGradientTarget(i, j, k, 2) - _transformedNormalisedGradientSource(i, j, k, 2)) * (_normalisedGradientTarget(i, j, k, 2) - _transformedNormalisedGradientSource(i, j, k, 2));
          n++;
        }
      }
    }
  }

  // Return similarity measure
  if (n > 0) {
    return -sum / double(n);
  } else {
    cerr << "irtkImageGradientFreeFormRegistration2::EvaluateNGD: No samples available" << endl;
    return 0;
  }
}

void irtkImageGradientFreeFormRegistration2::EvaluateGradientNGD()
{
  int i, j, k;

  // Print debugging information
  this->Debug("irtkImageGradientFreeFormRegistration2::EvaluateGradientNGD");

  // Compute gradient
  for (k = 0; k < _normalisedGradientTarget.GetZ(); k++) {
    for (j = 0; j < _normalisedGradientTarget.GetY(); j++) {
      for (i = 0; i < _normalisedGradientTarget.GetX(); i++) {
    	if ((_distanceMask(i, j, k) == 0) && (_transformedNormalisedGradientSource(i, j, k, 0) >= -1)) {
          _similarityGradient(i, j, k, 0) = 2.0 * ((_normalisedGradientTarget(i, j, k, 0) - _transformedNormalisedGradientSource(i, j, k, 0)) * _transformedNormalisedGradientSourceGradient[0](i, j, k, 0) +
        		                                   (_normalisedGradientTarget(i, j, k, 1) - _transformedNormalisedGradientSource(i, j, k, 1)) * _transformedNormalisedGradientSourceGradient[1](i, j, k, 0) +
        		                                   (_normalisedGradientTarget(i, j, k, 2) - _transformedNormalisedGradientSource(i, j, k, 2)) * _transformedNormalisedGradientSourceGradient[2](i, j, k, 0));

          _similarityGradient(i, j, k, 1) = 2.0 * ((_normalisedGradientTarget(i, j, k, 0) - _transformedNormalisedGradientSource(i, j, k, 0)) * _transformedNormalisedGradientSourceGradient[0](i, j, k, 1) +
												   (_normalisedGradientTarget(i, j, k, 1) - _transformedNormalisedGradientSource(i, j, k, 1)) * _transformedNormalisedGradientSourceGradient[1](i, j, k, 1) +
												   (_normalisedGradientTarget(i, j, k, 2) - _transformedNormalisedGradientSource(i, j, k, 2)) * _transformedNormalisedGradientSourceGradient[2](i, j, k, 1));

          _similarityGradient(i, j, k, 2) = 2.0 * ((_normalisedGradientTarget(i, j, k, 0) - _transformedNormalisedGradientSource(i, j, k, 0)) * _transformedNormalisedGradientSourceGradient[0](i, j, k, 2) +
												   (_normalisedGradientTarget(i, j, k, 1) - _transformedNormalisedGradientSource(i, j, k, 1)) * _transformedNormalisedGradientSourceGradient[1](i, j, k, 2) +
												   (_normalisedGradientTarget(i, j, k, 2) - _transformedNormalisedGradientSource(i, j, k, 2)) * _transformedNormalisedGradientSourceGradient[2](i, j, k, 2));
        }
      }
    }
  }
}

double irtkImageGradientFreeFormRegistration2::EvaluateNGP()
{
  int i, j, k, n;
  double sum;

  // Print debugging information
  this->Debug("irtkImageGradientFreeFormRegistration2::EvaluateNGP");

  // Initialize metric
  n = 0;
  sum = 0;

  // Compute metric
  for (k = 0; k < _normalisedGradientTarget.GetZ(); k++) {
    for (j = 0; j < _normalisedGradientTarget.GetY(); j++) {
      for (i = 0; i < _normalisedGradientTarget.GetX(); i++) {
        if ((_distanceMask(i, j, k) == 0) && (_transformedNormalisedGradientSource(i, j, k, 0) >= -1)) {
          sum += (_normalisedGradientTarget(i, j, k, 0) * _transformedNormalisedGradientSource(i, j, k, 0)) +
         		 (_normalisedGradientTarget(i, j, k, 1) * _transformedNormalisedGradientSource(i, j, k, 1)) +
          		 (_normalisedGradientTarget(i, j, k, 2) * _transformedNormalisedGradientSource(i, j, k, 2));
          n++;
        }
      }
    }
  }

  // Return similarity measure
  if (n > 0) {
	return sum / double(n);
  } else {
    cerr << "irtkImageGradientFreeFormRegistration2::EvaluateNGP: No samples available" << endl;
    return 0;
  }
}

void irtkImageGradientFreeFormRegistration2::EvaluateGradientNGP()
{
  int i, j, k;

  // Print debugging information
  this->Debug("irtkImageGradientFreeFormRegistration2::EvaluateGradientNGP");

  // Compute gradient
  for (k = 0; k < _normalisedGradientTarget.GetZ(); k++) {
    for (j = 0; j < _normalisedGradientTarget.GetY(); j++) {
      for (i = 0; i < _normalisedGradientTarget.GetX(); i++) {
    	if ((_distanceMask(i, j, k) == 0) && (_transformedNormalisedGradientSource(i, j, k, 0) >= -1)) {
          _similarityGradient(i, j, k, 0) = (_normalisedGradientTarget(i, j, k, 0) * _transformedNormalisedGradientSourceGradient[0](i, j, k, 0)) +
        		                            (_normalisedGradientTarget(i, j, k, 1) * _transformedNormalisedGradientSourceGradient[1](i, j, k, 0)) +
        		                            (_normalisedGradientTarget(i, j, k, 2) * _transformedNormalisedGradientSourceGradient[2](i, j, k, 0));

		  _similarityGradient(i, j, k, 1) = (_normalisedGradientTarget(i, j, k, 0) * _transformedNormalisedGradientSourceGradient[0](i, j, k, 1)) +
											(_normalisedGradientTarget(i, j, k, 1) * _transformedNormalisedGradientSourceGradient[1](i, j, k, 1)) +
											(_normalisedGradientTarget(i, j, k, 2) * _transformedNormalisedGradientSourceGradient[2](i, j, k, 1));

		  _similarityGradient(i, j, k, 2) = (_normalisedGradientTarget(i, j, k, 0) * _transformedNormalisedGradientSourceGradient[0](i, j, k, 2)) +
											(_normalisedGradientTarget(i, j, k, 1) * _transformedNormalisedGradientSourceGradient[1](i, j, k, 2)) +
											(_normalisedGradientTarget(i, j, k, 2) * _transformedNormalisedGradientSourceGradient[2](i, j, k, 2));
        }
      }
    }
  }
}

double irtkImageGradientFreeFormRegistration2::EvaluateNGS()
{
  int i, j, k, n;
  double sum;

  // Print debugging information
  this->Debug("irtkImageGradientFreeFormRegistration2::EvaluateNGS");

  // Initialize metric
  n = 0;
  sum = 0;

  // Compute metric
  for (k = 0; k < _normalisedGradientTarget.GetZ(); k++) {
    for (j = 0; j < _normalisedGradientTarget.GetY(); j++) {
      for (i = 0; i < _normalisedGradientTarget.GetX(); i++) {
    	if ((_distanceMask(i, j, k) == 0) && (_transformedNormalisedGradientSource(i, j, k, 0) >= -1)) {
          sum += (_normalisedGradientTarget(i, j, k, 0) * _transformedNormalisedGradientSource(i, j, k, 0)) +
           		 (_normalisedGradientTarget(i, j, k, 1) * _transformedNormalisedGradientSource(i, j, k, 1)) +
            	 (_normalisedGradientTarget(i, j, k, 2) * _transformedNormalisedGradientSource(i, j, k, 2))
            	 *
            	 (_normalisedGradientTarget(i, j, k, 0) * _transformedNormalisedGradientSource(i, j, k, 0)) +
            	 (_normalisedGradientTarget(i, j, k, 1) * _transformedNormalisedGradientSource(i, j, k, 1)) +
            	 (_normalisedGradientTarget(i, j, k, 2) * _transformedNormalisedGradientSource(i, j, k, 2));
          n++;
        }
      }
    }
  }

  // Return similarity measure
  if (n > 0) {
	return sum / double(n);
  } else {
    cerr << "irtkImageGradientFreeFormRegistration2::EvaluateNGS: No samples available" << endl;
    return 0;
  }
}

void irtkImageGradientFreeFormRegistration2::EvaluateGradientNGS()
{
  int i, j, k;

  // Print debugging information
  this->Debug("irtkImageGradientFreeFormRegistration2::EvaluateGradientNGS");

  // Compute gradient
  for (k = 0; k < _normalisedGradientTarget.GetZ(); k++) {
    for (j = 0; j < _normalisedGradientTarget.GetY(); j++) {
      for (i = 0; i < _normalisedGradientTarget.GetX(); i++) {
        if ((_distanceMask(i, j, k) == 0) && (_transformedNormalisedGradientSource(i, j, k, 0) >= -1)) {
          _similarityGradient(i, j, k, 0) = 2.0 * ((_normalisedGradientTarget(i, j, k, 0) * _transformedNormalisedGradientSource(i, j, k, 0)) +
												   (_normalisedGradientTarget(i, j, k, 1) * _transformedNormalisedGradientSource(i, j, k, 1)) +
												   (_normalisedGradientTarget(i, j, k, 2) * _transformedNormalisedGradientSource(i, j, k, 2)))
												  *
												  ((_normalisedGradientTarget(i, j, k, 0) * _transformedNormalisedGradientSourceGradient[0](i, j, k, 0)) +
	        		                               (_normalisedGradientTarget(i, j, k, 1) * _transformedNormalisedGradientSourceGradient[1](i, j, k, 0)) +
	        		                               (_normalisedGradientTarget(i, j, k, 2) * _transformedNormalisedGradientSourceGradient[2](i, j, k, 0)));

		  _similarityGradient(i, j, k, 1) = 2.0 * ((_normalisedGradientTarget(i, j, k, 0) * _transformedNormalisedGradientSource(i, j, k, 0)) +
												   (_normalisedGradientTarget(i, j, k, 1) * _transformedNormalisedGradientSource(i, j, k, 1)) +
												   (_normalisedGradientTarget(i, j, k, 2) * _transformedNormalisedGradientSource(i, j, k, 2)))
												  *
												  ((_normalisedGradientTarget(i, j, k, 0) * _transformedNormalisedGradientSourceGradient[0](i, j, k, 1)) +
												   (_normalisedGradientTarget(i, j, k, 1) * _transformedNormalisedGradientSourceGradient[1](i, j, k, 1)) +
												   (_normalisedGradientTarget(i, j, k, 2) * _transformedNormalisedGradientSourceGradient[2](i, j, k, 1)));

		  _similarityGradient(i, j, k, 2) = 2.0 * ((_normalisedGradientTarget(i, j, k, 0) * _transformedNormalisedGradientSource(i, j, k, 0)) +
												   (_normalisedGradientTarget(i, j, k, 1) * _transformedNormalisedGradientSource(i, j, k, 1)) +
												   (_normalisedGradientTarget(i, j, k, 2) * _transformedNormalisedGradientSource(i, j, k, 2)))
												  *
												  ((_normalisedGradientTarget(i, j, k, 0) * _transformedNormalisedGradientSourceGradient[0](i, j, k, 2)) +
												   (_normalisedGradientTarget(i, j, k, 1) * _transformedNormalisedGradientSourceGradient[1](i, j, k, 2)) +
												   (_normalisedGradientTarget(i, j, k, 2) * _transformedNormalisedGradientSourceGradient[2](i, j, k, 2)));
        }
      }
    }
  }
}

