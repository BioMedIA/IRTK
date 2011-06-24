/*=========================================================================

  Library   : Image Registration Toolkit (IRTK)
  Module    : $Id$
  Copyright : Imperial College, Department of Computing
              Visual Information Processing (VIP), 2008 onwards
  Date      : $Date$
  Version   : $Revision$
  Changes   : $Author$

=========================================================================*/

#include <irtkCardiac.h>

#include <irtkGradientImage.h>

#include <irtkGradientDescentConstrainedOptimizer.h>

#include <irtkResamplingWithPadding.h>

#include <irtkGaussianBlurring.h>

#include <irtkSegmentationFunction.h>

#undef HAS_TBB

// Used as temporary memory for transformed intensities
extern irtkGreyImage **_mtmpImage;

irtkGreyImage **_utmpImage;

irtkRealImage *_mweight;

// The original target and source images
extern irtkGreyImage **tmp_mtarget, **tmp_msource;

irtkGreyImage **tmp_mutarget,**tmp_musource;

#include <irtkMultiThreadedImageFreeFormRegistration.h>

irtkCardiac3DImageFreeFormRegistration::irtkCardiac3DImageFreeFormRegistration()
{
  // Print debugging information
  this->Debug("irtkCardiac3DImageFreeFormRegistration::irtkCardiac3DImageFreeFormRegistration");

  // Set inputs
  _utarget = NULL;
  _usource = NULL;

  // Set metric
  _umetric = NULL;

  // Allocate interpolation object
  _uinterpolator = NULL;
}

irtkCardiac3DImageFreeFormRegistration::~irtkCardiac3DImageFreeFormRegistration()
{

}

void irtkCardiac3DImageFreeFormRegistration::Initialize()
{
  int i, j, k, t;
  
  // Print debugging information
  this->Debug("irtkCardiac3DImageFreeFormRegistration::Initialize");

  // Initialize base class
  this->irtkMultipleImageFreeFormRegistration::Initialize();

  // Check that the t-dimensions are the same for both images.
  // Do not currently support different t-dimensions.
  for (i = 0; i < _numberOfuImages; i++) {
    if (_utarget[i]->GetT() != _usource[i]->GetT()) {
      cerr << this->NameOfClass() << "::Initialize() : Images no. " << i << " have different t-dimensions." << endl;
      exit(1);
    }
  }

  _usource_x1 = new double[_numberOfuImages];
  _usource_y1 = new double[_numberOfuImages];
  _usource_z1 = new double[_numberOfuImages];
  _usource_x2 = new double[_numberOfuImages];
  _usource_y2 = new double[_numberOfuImages];
  _usource_z2 = new double[_numberOfuImages];
  _uinterpolator = new irtkInterpolateImageFunction *[_numberOfuImages];

  // Initialize pointers
  _utmpImage         = NULL;
  _uaffdLookupTable  = NULL;
  _umffdLookupTable  = NULL;
  _mweight           = NULL;
  _myoprob           = NULL;
  _omega             = NULL;
  _comega            = NULL;

  _mweight = new irtkRealImage(_target[0]->GetImageAttributes());
  _myoprob = new irtkRealImage(_target[0]->GetImageAttributes());
  _omega = new irtkRealImage(_target[0]->GetImageAttributes());

  //this->ExtendThreshold(_threshold,14);
  //evaluate weight
  if (_target[0]->GetX() != _threshold->GetX() ||
	  _target[0]->GetY() != _threshold->GetY() ||
	  _target[0]->GetZ() != _threshold->GetZ()){
		  cerr << "untagged first target's dimension must equal to segmentation's dimension! sorry for the trouble please use transformation tools to achieve this" <<endl;
		  exit(1);
  }
  this->EvaluateWeight(_mweight,_target[0],_threshold);
  _mweight->Write("weight.nii");
  irtkRealImage tmp(_threshold->GetImageAttributes());
  for (t = 0; t < _threshold->GetT(); t++) {
	  for (k = 0; k < _threshold->GetZ(); k++) {
		  for (j = 0; j < _threshold->GetY(); j++) {
			  for (i = 0; i < _threshold->GetX(); i++) {
				  if(_threshold->GetAsDouble(i,j,k,t) > 2){
					  tmp.PutAsDouble(i,j,k,t,1 - _mweight->GetAsDouble(i,j,k,t));
				  }
			  }
		  }
	  }
  }
  tmp.Write("tweight.nii");
  irtkGaussian guassian;
  double denom;
  //evaluate myocardium probability
  if (this->_Lambda2 > 0) {
	  this->EvaluateMyoProb1(_myoprob,_threshold,guassian,denom);
  }
  //extend threshold
  ExtendThreshold (_threshold, 3);
  _threshold->Write("thresholdafterdilation.nii.gz");
  //evaluate myocardium probability
  if (this->_Lambda2 > 0) {
	  this->EvaluateMyoProb2(_myoprob,_threshold,guassian,denom);
  }
  //resample threshold and weight using isotropic shape based interpolation
  irtkInterpolateImageFunction *interpolator = new irtkShapeBasedInterpolateImageFunction;
  irtkInterpolateImageFunction *interpolatorw = new irtkLinearInterpolateImageFunction;
  irtkInterpolateImageFunction *interpolatorp = new irtkLinearInterpolateImageFunction;
  irtkInterpolateImageFunction *interpolatoro = new irtkLinearInterpolateImageFunction;
  //padding
  for (t = 0; t < _threshold->GetT(); t++) {
	  for (k = 0; k < _threshold->GetZ(); k++) {
		  for (j = 0; j < _threshold->GetY(); j++) {
			  for (i = 0; i < _threshold->GetX(); i++) {
				  if(_threshold->GetAsDouble(i,j,k,t) == 3){
					  _threshold->PutAsDouble(i,j,k,t,2);
				  }else if(_threshold->GetAsDouble(i,j,k,t) == 2){
					  _threshold->PutAsDouble(i,j,k,t,3);
				  }
				  _mweight->PutAsDouble(i,j,k,t,_mweight->GetAsDouble(i,j,k,t)*1000.0);
			  }
		  }
	  }
  }
  //resample to isotropic use shape based interpolation
  double xsize, ysize, zsize, size;
  _threshold->GetPixelSize(&xsize, &ysize, &zsize);
  size = xsize;
  size = (size < ysize) ? size : ysize;
  size = (size < zsize) ? size : zsize;

  irtkResampling<irtkGreyPixel> resampling(size,size,size);
  resampling.SetInput (_threshold);
  resampling.SetOutput(_threshold);
  resampling.SetInterpolator(interpolator);
  resampling.Run();
  irtkResampling<irtkRealPixel> resamplingw(size,size,size);
  resamplingw.SetInput (_mweight);
  resamplingw.SetOutput(_mweight);
  resamplingw.SetInterpolator(interpolatorw);
  resamplingw.Run();

  if (this->_Lambda2 > 0) {
	  irtkResampling<irtkRealPixel> resamplingp(size,size,size);
	  resamplingw.SetInput (_myoprob);
	  resamplingw.SetOutput(_myoprob);
	  resamplingw.SetInterpolator(interpolatorp);
	  resamplingw.Run();
	  _myoprob->Write("myoprob.nii.gz");

	  irtkResampling<irtkRealPixel> resamplingo(size,size,size);
	  resamplingw.SetInput (_omega);
	  resamplingw.SetOutput(_omega);
	  resamplingw.SetInterpolator(interpolatoro);
	  resamplingw.Run();
  }

  //padding
  for (t = 0; t < _threshold->GetT(); t++) {
		for (k = 0; k < _threshold->GetZ(); k++) {
			for (j = 0; j < _threshold->GetY(); j++) {
				for (i = 0; i < _threshold->GetX(); i++) {
					if(_threshold->GetAsDouble(i,j,k,t) == 2){
						_threshold->PutAsDouble(i,j,k,t,3);
					}else if(_threshold->GetAsDouble(i,j,k,t) == 3){
						_threshold->PutAsDouble(i,j,k,t,2);
					}
					_mweight->PutAsDouble(i,j,k,t,_mweight->GetAsDouble(i,j,k,t)/1000.0);
				}
			}
		}
	}
  delete interpolator;
  delete interpolatorw;
  delete interpolatorp;
  delete interpolatoro;
}

void irtkCardiac3DImageFreeFormRegistration::Initialize(int level)
{
	int i, j, k, l, t, n;
	double x, y, z;
	float *ptr;
	double dx, dy, dz, temp;
	irtkGreyPixel target_min, target_max, target_nbins;
	irtkGreyPixel source_min, source_max, source_nbins;

  // Print debugging information
  this->Debug("irtkCardiac3DImageFreeFormRegistration::Initialize(int)");

  // Initialize base class
  this->irtkMultipleImageFreeFormRegistration::Initialize(level);

  
  // Allocate memory for temporary images
  tmp_mutarget = new irtkGreyImage*[_numberOfuImages];
  tmp_musource = new irtkGreyImage*[_numberOfuImages];

  // Allocate memory for metric
  _umetric = new irtkSimilarityMetric*[_numberOfuImages];

  for (n = 0; n < _numberOfuImages; n++) {

      target_max = MIN_GREY;
      target_min = MAX_GREY;
      source_max = MIN_GREY;
      source_min = MAX_GREY;

      // Copy source and target to temp space
      tmp_mutarget[n] = new irtkGreyImage(*_utarget[n]);
      tmp_musource[n] = new irtkGreyImage(*_usource[n]);

      // Swap source and target with temp space copies
      swap(tmp_mutarget[n], _utarget[n]);
      swap(tmp_musource[n], _usource[n]);

      // Blur images if necessary
      if (_TargetBlurring[level] > 0) {
          cout << "Blurring utarget ... "; cout.flush();
          irtkGaussianBlurringWithPadding<irtkGreyPixel> blurring(_uTargetBlurring[level], _TargetPadding);
          blurring.SetInput (_utarget[n]);
          blurring.SetOutput(_utarget[n]);
          blurring.Run();
          cout << "done" << endl;
      }

      if (_SourceBlurring[level] > 0) {
          cout << "Blurring usource ... "; cout.flush();
          irtkGaussianBlurring<irtkGreyPixel> blurring(_uSourceBlurring[level]);
          blurring.SetInput (_usource[n]);
          blurring.SetOutput(_usource[n]);
          blurring.Run();
          cout << "done" << endl;
      }

      _utarget[n]->GetPixelSize(&dx, &dy, &dz);
      temp = fabs(_uTargetResolution[0][0]-dx) + fabs(_uTargetResolution[0][1]-dy) + fabs(_uTargetResolution[0][2]-dz);

      if (level > 0 || temp > 0.000001) {
          cout << "Resampling utarget ... "; cout.flush();
          // Create resampling filter
          irtkResamplingWithPadding<irtkGreyPixel> resample(_uTargetResolution[level][0],
              _uTargetResolution[level][1],
              _uTargetResolution[level][2],
              _TargetPadding);
          resample.SetInput (_utarget[n]);
          resample.SetOutput(_utarget[n]);
          resample.Run();
          cout << "done" << endl;
      }

      _usource[n]->GetPixelSize(&dx, &dy, &dz);
      temp = fabs(_uSourceResolution[0][0]-dx) + fabs(_uSourceResolution[0][1]-dy) + fabs(_uSourceResolution[0][2]-dz);

      if (level > 0 || temp > 0.000001) {
          cout << "Resampling usource ... "; cout.flush();
          // Create resampling filter
          irtkResamplingWithPadding<irtkGreyPixel> resample(_uSourceResolution[level][0],
              _uSourceResolution[level][1],
              _uSourceResolution[level][2], MIN_GREY);

          resample.SetInput (_usource[n]);
          resample.SetOutput(_usource[n]);
          resample.Run();
          cout << "done" << endl;
      }

      // Find out the min and max values in target image, ignoring padding
      for (t = 0; t < _utarget[n]->GetT(); t++) {
          for (k = 0; k < _utarget[n]->GetZ(); k++) {
              for (j = 0; j < _utarget[n]->GetY(); j++) {
                  for (i = 0; i < _utarget[n]->GetX(); i++) {
                      if (_utarget[n]->Get(i, j, k, t) > _TargetPadding) {
                          if (_utarget[n]->Get(i, j, k, t) > target_max)
                              target_max = _utarget[n]->Get(i, j, k, t);
                          if (_utarget[n]->Get(i, j, k, t) < target_min)
                              target_min = _utarget[n]->Get(i, j, k, t);
                      } else {
                          _utarget[n]->Put(i, j, k, t, _TargetPadding);
                      }
                  }
              }
          }
      }

      // Find out the min and max values in source image, ignoring padding
      for (t = 0; t < _usource[n]->GetT(); t++) {
          for (k = 0; k < _usource[n]->GetZ(); k++) {
              for (j = 0; j < _usource[n]->GetY(); j++) {
                  for (i = 0; i < _usource[n]->GetX(); i++) {
                      if (_usource[n]->Get(i, j, k, t) > source_max)
                          source_max = _usource[n]->Get(i, j, k, t);
                      if (_usource[n]->Get(i, j, k, t) < source_min)
                          source_min = _usource[n]->Get(i, j, k, t);
                  }
              }
          }
      }
      // Check whether dynamic range of data is not to large
      if (target_max - target_min > MAX_GREY) {
          cerr << this->NameOfClass()
              << "::Initialize: Dynamic range of target is too large" << endl;
          exit(1);
      } else {
          for (t = 0; t < _utarget[n]->GetT(); t++) {
              for (k = 0; k < _utarget[n]->GetZ(); k++) {
                  for (j = 0; j < _utarget[n]->GetY(); j++) {
                      for (i = 0; i < _utarget[n]->GetX(); i++) {
                          if (_utarget[n]->Get(i, j, k, t) > _TargetPadding) {
                              _utarget[n]->Put(i, j, k, t, _utarget[n]->Get(i, j, k, t) - target_min);
                          } else {
                              _utarget[n]->Put(i, j, k, t, _TargetPadding - 1);
                          }
                      }
                  }
              }
          }
      }

      if ((_SimilarityMeasure == SSD) || (_SimilarityMeasure == CC) ||
          (_SimilarityMeasure == LC)  || (_SimilarityMeasure == K) || (_SimilarityMeasure == ML)) {
              if (source_max - target_min > MAX_GREY) {
                  cerr << this->NameOfClass()
                      << "::Initialize: Dynamic range of source is too large" << endl;
                  exit(1);
              } else {
                  for (t = 0; t < _usource[n]->GetT(); t++) {
                      for (k = 0; k < _usource[n]->GetZ(); k++) {
                          for (j = 0; j < _usource[n]->GetY(); j++) {
                              for (i = 0; i < _usource[n]->GetX(); i++) {
                                  _usource[n]->Put(i, j, k, t, _usource[n]->Get(i, j, k, t) - target_min);
                              }
                          }
                      }
                  }
              }
      } else {
          if (source_max - source_min > MAX_GREY) {
              cerr << this->NameOfClass()
                  << "::Initialize: Dynamic range of source is too large" << endl;
              exit(1);
          } else {
              for (t = 0; t < _usource[n]->GetT(); t++) {
                  for (k = 0; k < _usource[n]->GetZ(); k++) {
                      for (j = 0; j < _usource[n]->GetY(); j++) {
                          for (i = 0; i < _usource[n]->GetX(); i++) {
                              _usource[n]->Put(i, j, k, t, _usource[n]->Get(i, j, k, t) - source_min);
                          }
                      }
                  }
              }
          }
      }

      // Pad target image if necessary
      irtkPadding(*_utarget[n], _TargetPadding);

      switch (_SimilarityMeasure) {
      case SSD:
          _umetric[n] = new irtkSSDSimilarityMetric;
          break;
      case CC:
          // Rescale images by an integer factor if necessary
          _umetric[n] = new irtkCrossCorrelationSimilarityMetric;
          break;
      case JE:
          // Rescale images by an integer factor if necessary
          target_nbins = irtkCalculateNumberOfBins(_utarget[n], _NumberOfBins,
              target_min, target_max);
          source_nbins = irtkCalculateNumberOfBins(_usource[n], _NumberOfBins,
              source_min, source_max);
          _umetric[n] = new irtkJointEntropySimilarityMetric(target_nbins, source_nbins);
          break;
      case MI:
          // Rescale images by an integer factor if necessary
          target_nbins = irtkCalculateNumberOfBins(_utarget[n], _NumberOfBins,
              target_min, target_max);
          source_nbins = irtkCalculateNumberOfBins(_usource[n], _NumberOfBins,
              source_min, source_max);
          _umetric[n] = new irtkMutualInformationSimilarityMetric(target_nbins, source_nbins);
          break;
      case NMI:
          // Rescale images by an integer factor if necessary
          target_nbins = irtkCalculateNumberOfBins(_utarget[n], _NumberOfBins,
              target_min, target_max);
          source_nbins = irtkCalculateNumberOfBins(_usource[n], _NumberOfBins,
              source_min, source_max);
          _umetric[n] = new irtkNormalisedMutualInformationSimilarityMetric(target_nbins, source_nbins);
          break;
      case CR_XY:
          // Rescale images by an integer factor if necessary
          target_nbins = irtkCalculateNumberOfBins(_utarget[n], _NumberOfBins,
              target_min, target_max);
          source_nbins = irtkCalculateNumberOfBins(_usource[n], _NumberOfBins,
              source_min, source_max);
          _umetric[n] = new irtkCorrelationRatioXYSimilarityMetric(target_nbins, source_nbins);
          break;
      case CR_YX:
          // Rescale images by an integer factor if necessary
          target_nbins = irtkCalculateNumberOfBins(_utarget[n], _NumberOfBins,
              target_min, target_max);
          source_nbins = irtkCalculateNumberOfBins(_usource[n], _NumberOfBins,
              source_min, source_max);
          _umetric[n] = new irtkCorrelationRatioYXSimilarityMetric(target_nbins, source_nbins);
          break;
      case LC:
          _umetric[n] = new irtkLabelConsistencySimilarityMetric;
          break;
      case K:
          // Rescale images by an integer factor if necessary
          target_nbins = irtkCalculateNumberOfBins(_utarget[n], _NumberOfBins,
              target_min, target_max);
          source_nbins = irtkCalculateNumberOfBins(_usource[n], _NumberOfBins,
              source_min, source_max);
          _umetric[n] = new irtkKappaSimilarityMetric(target_nbins, source_nbins);
          break;
      case ML:
          // Rescale images by an integer factor if necessary
          _umetric[n] = new irtkMLSimilarityMetric(classification);
          if (_umetric[n]==NULL) {
              cerr<<"Please, do not forget to set the ML metric!!!"<<endl;
          }
          break;
      default:
          cerr<<"Can not recognize similarity metric type!!!"<<endl;
          exit(1);
      }
      // Print some debugging information
      cout << "uTarget image (reference) no. " << n << endl;
      _utarget[n]->Print();
      cout << "Range is from " << target_min << " to " << target_max << endl;

      cout << "uSource image (transform) no. " << n << endl;
      _usource[n]->Print();
      cout << "Range is from " << source_min << " to " << source_max << endl;
  }

  for (n = 0; n < _numberOfuImages; n++) {
    // Setup the interpolator
    _uinterpolator[n] = irtkInterpolateImageFunction::New(_InterpolationMode, _usource[n]);

    // Setup interpolation for the source image
    _uinterpolator[n]->SetInput(_usource[n]);
    _uinterpolator[n]->Initialize();

    // Calculate the source image domain in which we can interpolate
    _uinterpolator[n]->Inside(_usource_x1[n], _usource_y1[n], _usource_z1[n],
                             _usource_x2[n], _usource_y2[n], _usource_z2[n]);
  }

  // Allocate memory for metric
  _utmpMetricA = new irtkSimilarityMetric*[_numberOfuImages];
  _utmpMetricB = new irtkSimilarityMetric*[_numberOfuImages];
  // Allocate memory for temp image
  _utmpImage = new irtkGreyImage*[_numberOfuImages];

  // Allocate memory for lookup tables
  _uaffdLookupTable = new float*[_numberOfuImages];
  _umffdLookupTable = new float*[_numberOfuImages];

  for (l = 0; l < _numberOfuImages; l++) {
	    // Allocate memory for metric
  _utmpMetricA[l] = irtkSimilarityMetric::New(_umetric[l]);
  _utmpMetricB[l] = irtkSimilarityMetric::New(_umetric[l]);
    _utmpImage[l] = new irtkGreyImage(_utarget[l]->GetX(),
                                     _utarget[l]->GetY(),
                                     _utarget[l]->GetZ(),
                                     _utarget[l]->GetT());

    n = _utarget[l]->GetNumberOfVoxels() * 3 / _utarget[l]->GetT();

    // Allocate memory for lookup table for single-level FFD
    _uaffdLookupTable[l]  = new float[n];

    // Allocate memory for lookup table for multi-level FFD
    _umffdLookupTable[l]  = new float[n];

    // Initialize lookup table for multi-level FFD (this is done only once)
    ptr = _umffdLookupTable[l];
    for (k = 0; k < _utarget[l]->GetZ(); k++) {
      for (j = 0; j < _utarget[l]->GetY(); j++) {
        for (i = 0; i < _utarget[l]->GetX(); i++) {
          x = i;
          y = j;
          z = k;
          _utarget[l]->ImageToWorld(x, y, z);
          _mffd->Transform(x, y, z);
          ptr[0] = x;
          ptr[1] = y;
          ptr[2] = z;
          ptr += 3;
        }
      }
    }
  }

  // Allocate memory for comega
  _comega = new double[_affd->NumberOfDOFs()];
  for (i = 0; i < _affd->NumberOfDOFs(); i ++)
	  _comega[i] = 0;

  // Padding of FFD
  //irtkPadding(*_threshold, 2 , _affd);
  //irtkPadding(tmp_mutarget, this->_TargetPadding, _affd, _numberOfuImages);
}

void irtkCardiac3DImageFreeFormRegistration::Finalize()
{
  // Print debugging information
  this->Debug("irtkCardiac3DImageFreeFormRegistration::Finalize");

  delete []_usource_x1;
  delete []_usource_y1;
  delete []_usource_z1;
  delete []_usource_x2;
  delete []_usource_y2;
  delete []_usource_z2;
  delete []_uinterpolator;

  delete _mweight;
  delete _omega;
  delete _myoprob;

  // Finalize base class
  this->irtkMultipleImageFreeFormRegistration::Finalize();
}

void irtkCardiac3DImageFreeFormRegistration::Finalize(int level)
{
  // Print debugging information
  this->Debug("irtkCardiac3DImageFreeFormRegistration::Finalize(int)");
  // Finalize base class
  this->irtkMultipleImageFreeFormRegistration::Finalize(level);

  int n;
  // Swap source and target back with temp space copies (see Initialize)
  swap(tmp_mutarget, _utarget);
  swap(tmp_musource, _usource);

#ifdef HAS_TBB
  irtkSimilarityMetric *umetric;
  while (queue.size() > 0) {
    queue.pop(umetric);
    delete umetric;
  }
#endif
  for (n = 0; n < _numberOfuImages; n++) {
	  delete tmp_mutarget[n];
	  delete tmp_musource[n];
	  delete _utmpImage[n];
	  delete _uaffdLookupTable[n];
	  delete _umffdLookupTable[n];
      delete _uinterpolator[n];
	  delete _utmpMetricA[n];
	  delete _utmpMetricB[n];
  }
  delete []tmp_mutarget;
  delete []tmp_musource;

  delete []_utmpImage;
  delete []_utmpMetricA;
  delete []_utmpMetricB;
  delete []_uaffdLookupTable;
  delete []_umffdLookupTable;
  delete []_comega;
}

double irtkCardiac3DImageFreeFormRegistration::VolumePreservationPenalty()
{
	int i, j, k, index;
	double x, y, z, penalty, jacobian;
	irtkMatrix jac,tmp_jac;

	penalty = 0;
	for (k = 0; k < _affd->GetZ(); k++) {
		for (j = 0; j < _affd->GetY(); j++) {
			for (i = 0; i < _affd->GetX(); i++) {
				x = i;
				y = j;
				z = k;
				index = _affd->LatticeToIndex(x,y,z);
				if(_comega[index] > 0){
					_affd->LatticeToWorld(x, y, z);
					_affd->Jacobian(tmp_jac,x,y,z);
					// Calculate jacobian
					jac.Initialize(3, 3);
					_mffd->LocalJacobian(jac, x, y, z);

					// Subtract identity matrix
					tmp_jac(0, 0) = tmp_jac(0, 0) - 1;
					tmp_jac(1, 1) = tmp_jac(1, 1) - 1;
					tmp_jac(2, 2) = tmp_jac(2, 2) - 1;

					// Add jacobian
					jac += tmp_jac;
					// Determinant of Jacobian of deformation derivatives
					jacobian = (jac(0, 0)*jac(1, 1)*jac(2, 2) + jac(0, 1)*jac(1, 2)*jac(2, 0) +
						jac(0, 2)*jac(1, 0)*jac(2, 1) - jac(0, 2)*jac(1, 1)*jac(2, 0) -
						jac(0, 0)*jac(1, 2)*jac(2, 1) - jac(0, 1)*jac(1, 0)*jac(2, 2));
					if(jacobian < 0.0000001) jacobian = 0.0000001;
					//jacobian = _affd->irtkTransformation::Jacobian(x, y, z);
					//if (jacobian < 0.001)
					//jacobian = 0.001;
					// Torsten Rohlfing et al. MICCAI'01 (w/o scaling correction):
					penalty += _comega[index]*fabs(log(jacobian));
				}
			}
		}
	}

	// Normalize sum by number of DOFs
	return -penalty;
}

double irtkCardiac3DImageFreeFormRegistration::VolumePreservationPenalty(int index)
{	
    int i, j, k, i1, j1, k1, i2, j2, k2, count;
    double x, y, z, jacobian, penalty;

    _affd->IndexToLattice(index, i, j, k);
    penalty = 0;
    count = 0;
    k1 = (k-1)>0?(k-1):0;
    j1 = (j-1)>0?(j-1):0;
    i1 = (i-1)>0?(i-1):0;
    k2 = (k+2) < _affd->GetZ()? (k+2) : _affd->GetZ();
    j2 = (j+2) < _affd->GetY()? (j+2) : _affd->GetY();
    i2 = (i+2) < _affd->GetX()? (i+2) : _affd->GetX();
    for (k = k1; k < k2; k++) {
        for (j = j1; j < j2; j++) {
            for (i = i1; i < i2; i++) {
                x = i;
                y = j;
                z = k;
                index = _affd->LatticeToIndex(i,j,k);
                if(_comega[index] > 0){
                    x = i;
                    y = j;
                    z = k;
                    _affd->LatticeToWorld(x, y, z);
                    irtkMatrix jac,tmp_jac;
                    _affd->Jacobian(tmp_jac,x,y,z);
                    // Calculate jacobian
                    jac.Initialize(3, 3);
                    _mffd->LocalJacobian(jac, x, y, z);

                    // Subtract identity matrix
                    tmp_jac(0, 0) = tmp_jac(0, 0) - 1;
                    tmp_jac(1, 1) = tmp_jac(1, 1) - 1;
                    tmp_jac(2, 2) = tmp_jac(2, 2) - 1;

                    // Add jacobian
                    jac += tmp_jac;
                    // Determinant of Jacobian of deformation derivatives
                    jacobian = (jac(0, 0)*jac(1, 1)*jac(2, 2) + jac(0, 1)*jac(1, 2)*jac(2, 0) +
                        jac(0, 2)*jac(1, 0)*jac(2, 1) - jac(0, 2)*jac(1, 1)*jac(2, 0) -
                        jac(0, 0)*jac(1, 2)*jac(2, 1) - jac(0, 1)*jac(1, 0)*jac(2, 2));
                    if(jacobian < 0.0000001) jacobian = 0.0000001;
                    // Normalize sum by number of weights
                    penalty += _comega[index]*fabs(log(jacobian));
                    count ++;
                }
            }
        }
    }
    if(count > 0)
        return -penalty/count;
    else
        return 0;

	
}

void irtkCardiac3DImageFreeFormRegistration::UpdateLUT()
{
	this->irtkMultipleImageFreeFormRegistration::UpdateLUT();
	int i, j, k, n;
	double x, y, z;
	float *ptr2mffd;
	float *ptr2affd;

	// Print debugging information
	this->Debug("irtkCardiac3DImageFreeFormRegistration::UpdateLUT");

	for (n = 0; n < _numberOfuImages; n++) {
		ptr2affd = _uaffdLookupTable[n];
		ptr2mffd = _umffdLookupTable[n];
		for (k = 0; k < _utarget[n]->GetZ(); k++) {
			for (j = 0; j < _utarget[n]->GetY(); j++) {
				for (i = 0; i < _utarget[n]->GetX(); i++) {
					x = i;
					y = j;
					z = k;
					_utarget[n]->ImageToWorld(x, y, z);
					_affd->LocalDisplacement(x, y, z);
					ptr2affd[0] = x + ptr2mffd[0];
					ptr2affd[1] = y + ptr2mffd[1];
					ptr2affd[2] = z + ptr2mffd[2];
					ptr2mffd += 3;
					ptr2affd += 3;
				}
			}
		}
	}
}

double irtkCardiac3DImageFreeFormRegistration::Evaluate()
{
#ifndef HAS_TBB
  // Image coordinates
  int i, j, k, t, n;
  // World coordinates
  double x, y, z, wx, wy, wz, l1, l2;
  // Pointer to reference data
  irtkGreyPixel *ptr2target,*ptr2utarget;
  irtkGreyPixel *ptr2tmp,*ptr2utmp,threshold;
  irtkRealPixel weight;
  float *ptr;
  double *sweight = new double[_numberOfImages];
  double *tweight = new double[_numberOfuImages];
#endif

  // Print debugging information
  this->Debug("irtkCardiac3DImageFreeFormRegistration::Evaluate");

#ifdef HAS_TBB
  irtkMultiThreadedImageFreeFormRegistrationEvaluate evaluate(this);
  parallel_reduce(blocked_range<int>(0, _target->GetZ(), 1), evaluate);
#else
  l1 = 0; l2 = 0;
  for (n = 0; n < _numberOfImages; n++) {
	  // Initialize metric
  _metric[n]->Reset();
  sweight[n] = 0;
    // Loop over all voxels in the target (reference) volume
    ptr2target = _target[n]->GetPointerToVoxels();
    ptr2tmp    = _mtmpImage[n]->GetPointerToVoxels();
    for (t = 0; t < _target[n]->GetT(); t++) {
      ptr        = _mffdLookupTable[n];
      for (k = 0; k < _target[n]->GetZ(); k++) {
        for (j = 0; j < _target[n]->GetY(); j++) {
          for (i = 0; i < _target[n]->GetX(); i++) {
            // Check whether reference point is valid
            if (*ptr2target > _TargetPadding) {
              x = i;
              y = j;
              z = k;
              _target[n]->ImageToWorld(x, y, z);
			  wx = x; wy = y; wz = z;
			  _mweight->WorldToImage(wx,wy,wz);
              if(wx > 0 && wy > 0 && wz > 0
                  && wx < _mweight->GetX() - 1
                  && wy < _mweight->GetY() - 1
                  && wz < _mweight->GetZ() - 1 ){
                      weight = _mweight->GetAsDouble(round(wx),round(wy),round(wz));
              }else{
                  weight = 0;
              }
              _affd->LocalDisplacement(x, y, z);
              x += ptr[0];
              y += ptr[1];
              z += ptr[2];
              _source[n]->WorldToImage(x, y, z);
              // Check whether transformed point is inside volume
              if ((x > _source_x1[n]) && (x < _source_x2[n]) &&
                  (y > _source_y1[n]) && (y < _source_y2[n]) &&
                  ((_target[n]->GetZ() == 1 && round(z) == 0)
				  ||( _target[n]->GetZ() != 1 
				  && z > _source_z1[n] && z < _source_z2[n]))) {
					  // Add sample to metric
					  if ( weight > 0 ) {
						  *ptr2tmp =  round(_interpolator[n]->EvaluateInside(x, y, z, t));
						  _metric[n]->Add(*ptr2target, *ptr2tmp, weight);
						  sweight[n] += weight;
						  l1 += weight;
					  }
              } else {
                *ptr2tmp = -1;
              }
            }
            // Increment pointers to next voxel
            ptr2tmp++;
            ptr2target++;
            ptr += 3;
          }
        }
      }
    }
  }

  for (n = 0; n < _numberOfuImages; n++) {
	    _umetric[n]->Reset();
	    tweight[n] = 0;
    // Loop over all voxels in the target (reference) volume
    ptr2utarget = _utarget[n]->GetPointerToVoxels();
    ptr2utmp    = _utmpImage[n]->GetPointerToVoxels();
    for (t = 0; t < _utarget[n]->GetT(); t++) {
      ptr        = _umffdLookupTable[n];
      for (k = 0; k < _utarget[n]->GetZ(); k++) {
        for (j = 0; j < _utarget[n]->GetY(); j++) {
          for (i = 0; i < _utarget[n]->GetX(); i++) {
            // Check whether reference point is valid
            if (*ptr2utarget > _TargetPadding) {
              x = i; y = j; z = k;
              _utarget[n]->ImageToWorld(x, y, z);
			  wx = x; wy = y; wz = z;
			  _threshold->WorldToImage(wx,wy,wz);
              if(wx > 0 && wy > 0 && wz > 0
                  && wx < _threshold->GetX() - 1
                  && wy < _threshold->GetY() - 1
                  && wz < _threshold->GetZ() - 1 ){
                      threshold = _threshold->GetAsDouble(round(wx),round(wy),round(wz));
              }else{
                  threshold = 0;
              }
              wx = x; wy = y; wz = z;
              _mweight->WorldToImage(wx,wy,wz);
              if(wx > 0 && wy > 0 && wz > 0
                  && wx < _mweight->GetX() - 1
                  && wy < _mweight->GetY() - 1
                  && wz < _mweight->GetZ() - 1 ){
                  weight = _mweight->GetAsDouble(round(wx),round(wy),round(wz));
              }else{
                  weight = 0;
              }
              _affd->LocalDisplacement(x, y, z);
              x += ptr[0];
              y += ptr[1];
              z += ptr[2];
              _usource[n]->WorldToImage(x, y, z);
              // Check whether transformed point is inside volume
              if ((x > _usource_x1[n]) && (x < _usource_x2[n]) &&
                  (y > _usource_y1[n]) && (y < _usource_y2[n]) &&
                  ((_utarget[n]->GetZ() == 1 && round(z) == 0)
				  ||( _utarget[n]->GetZ() != 1 
				  && z > _usource_z1[n] && z < _usource_z2[n]))) {
					  // Add sample to metric
					  if (threshold > 2) {
						  *ptr2utmp =  round(_uinterpolator[n]->EvaluateInside(x, y, z, t));
						  _umetric[n]->Add(*ptr2utarget, *ptr2utmp,1.0 - weight);
						  tweight[n] += 1.0 - weight;
						  l2 += 1.0 - weight;
					  }
              } else {
                *ptr2utmp = -1;
              }
            }
            // Increment pointers to next voxel
            ptr2utmp++;
            ptr2utarget++;
            ptr += 3;
          }
        }
      }
    }
  }
#endif

  // Evaluate similarity measure
  double similarity = combine_mysimilarity(combine_mysimilarity(_metric,sweight,_numberOfImages)
	  ,combine_mysimilarity(_umetric,tweight,_numberOfuImages),l1,l2);

  // Add penalty for landmark regulation
  if (this->_Lregu > 0) {
	  similarity += this->_Lregu * this->LandMarkPenalty(-1);
  }
  // Add penalty for smoothness
  if (this->_Lambda1 > 0) {
    similarity += this->_Lambda1*this->SmoothnessPenalty();
  }
  // Add penalty for volume preservation
  if (this->_Lambda2 > 0) {
    similarity += this->_Lambda2*this->VolumePreservationPenalty();
  }
  // Add penalty for topology preservation
  if (this->_Lambda3 > 0) {
    similarity += this->_Lambda3*this->TopologyPreservationPenalty();
  }

  delete []sweight;
  delete []tweight;

  // Return similarity measure + penalty terms
  return similarity;
}

double irtkCardiac3DImageFreeFormRegistration::EvaluateDerivative(int index, double step)
{
  float *ptr;
  irtkPoint p1, p2;
  double bi, bj, bk, dx, dy, dz, p[3],x,y,z,wx,wy,wz,l1a,l1b,l2a,l2b;
  int i, j, k, i1, i2, j1, j2, k1, k2, dim, t, n,min,max;
  irtkGreyPixel *ptr2target, *ptr2tmp,threshold;
  irtkRealPixel weight;
  irtkSimilarityMetric **tmpMetricA, **tmpMetricB;
  irtkGreyPixel *ptr2utarget, *ptr2utmp;
  irtkSimilarityMetric **utmpMetricA, **utmpMetricB;
  double *sweight = new double[_numberOfImages];
  double *tweight = new double[_numberOfuImages];

  // Print debugging information
  this->Debug("irtkCardiac3DImageFreeFormRegistration::EvaluateDerivative(int, double)");

#ifdef HAS_TBB
  // Create similarity metric if necessary
  if (queue.pop_if_present(tmpMetricA) == false) {
    tmpMetricA = irtkSimilarityMetric::New(_metric);
  }
  // Create similarity metric if necessary
  if (queue.pop_if_present(tmpMetricB) == false) {
    tmpMetricB = irtkSimilarityMetric::New(_metric);
  }
    // Create similarity metric if necessary
  if (queue.pop_if_present(utmpMetricA) == false) {
    utmpMetricA = irtkSimilarityMetric::New(_umetric);
  }
  // Create similarity metric if necessary
  if (queue.pop_if_present(utmpMetricB) == false) {
    utmpMetricB = irtkSimilarityMetric::New(_umetric);
  }
#else
  tmpMetricA = _tmpMetricA;
  tmpMetricB = _tmpMetricB;
  utmpMetricA = _utmpMetricA;
  utmpMetricB = _utmpMetricB;
#endif

  // Calculate whether this DOF corresponds to x, y or z-displacement
  dim = int(index / (_affd->GetX()*_affd->GetY()*_affd->GetZ()));
  l1a=l2a=l1b=l2b=0;

  for (n = 0; n < _numberOfImages; n++) {
	  // Initialize metrics for forward and backward derivative steps
	  tmpMetricA[n]->Reset(_metric[n]);
	  tmpMetricB[n]->Reset(_metric[n]);
      sweight[n] = 0;
	  

    // Calculate bounding box of control point in world coordinates
    _affd->BoundingBox(index, p1, p2);
    _target[0]->WorldToImage(p1);
    _target[0]->WorldToImage(p2);

    // Calculate bounding box of control point in image coordinates
    _affd->MultiBoundingBox(_target[n], index, i1, j1, k1, i2, j2, k2, 1.0 / _SpeedupFactor);

    // Calculate incremental changes in lattice coordinates when looping
    // over target
    dx = (FFDLOOKUPTABLESIZE-1)/(p2._x-p1._x);
    dy = (FFDLOOKUPTABLESIZE-1)/(p2._y-p1._y);
    dz = (FFDLOOKUPTABLESIZE-1)/(p2._z-p1._z);

	min = round((FFDLOOKUPTABLESIZE-1)*(0.5 - 0.5/_SpeedupFactor));
	max = round((FFDLOOKUPTABLESIZE-1)*(0.5 + 0.5/_SpeedupFactor));

    // Loop over all voxels in the target (reference) volume
	for (t = 0; t < _target[n]->GetT(); t++) {
		for (k = k1; k <= k2; k++) {
			for (j = j1; j <= j2; j++) {
				ptr2target = _target[n]->GetPointerToVoxels(i1, j, k, t);
				ptr        = &(_affdLookupTable[n][3*_target[n]->VoxelToIndex(i1, j, k)]);
				ptr2tmp  = _mtmpImage[n]->GetPointerToVoxels(i1, j, k, t);
				for (i = i1; i <= i2; i++) {
					x = i; y = j; z = k;
					_target[n]->ImageToWorld(x,y,z);
					_target[0]->WorldToImage(x,y,z);
					if(round(dz*(z-p1._z))<=max && round(dz*(z-p1._z))>=min){
						bk = step * _localLookupTable[round(dz*(z-p1._z))];
						if(round(dy*(y-p1._y))<=max && round(dy*(y-p1._y))>=min){
							bj = bk * _localLookupTable[round(dy*(y-p1._y))];
							// Check whether reference point is valid
							if (*ptr2target > _TargetPadding && round(dx*(x-p1._x))<=max && round(dx*(x-p1._x))>=min) {
									bi = bj * _localLookupTable[round(dx*(x-p1._x))];
									// Delete old samples from both metrics
									wx = i;	wy = j;	wz = k;
									// Convert transformed point to image coordinates
									_target[n]->ImageToWorld(wx,wy,wz);
									_mweight->WorldToImage(wx,wy,wz);
									if(wx > 0 && wy > 0 && wz > 0
										&& wx < _mweight->GetX() - 1
										&& wy < _mweight->GetY() - 1
										&& wz < _mweight->GetZ() - 1 ){
											weight = _mweight->GetAsDouble(round(wx),round(wy),round(wz));
									}else{
										weight = 0;
									}
									if (*ptr2tmp != -1 && *ptr2target > _TargetPadding && weight > 0) {
										tmpMetricA[n]->Delete(*ptr2target, *ptr2tmp, weight);
										tmpMetricB[n]->Delete(*ptr2target, *ptr2tmp, weight);
									}

									p[0] = ptr[0];
									p[1] = ptr[1];
									p[2] = ptr[2];
									p[dim] += bi;

									_source[n]->WorldToImage(p[0], p[1], p[2]);

									// Check whether transformed point is inside volume
									if ((p[0] > _source_x1[n]) && (p[0] < _source_x2[n]) &&
										(p[1] > _source_y1[n]) && (p[1] < _source_y2[n]) &&
										((_target[n]->GetZ() == 1 && round(p[2]) == 0)
										||( _target[n]->GetZ() != 1 
										&& p[2] > _source_z1[n] && p[2] < _source_z2[n]))) {
											if (*ptr2target > _TargetPadding && weight > 0) {
											// Add sample to metric
											tmpMetricA[n]->Add(*ptr2target, round(_interpolator[n]->EvaluateInside(p[0], p[1], p[2], t)),weight);
											l1a += weight;
											sweight[n] += weight;
											}
									}

									p[0] = ptr[0];
									p[1] = ptr[1];
									p[2] = ptr[2];
									p[dim] -= bi;

									_source[n]->WorldToImage(p[0], p[1], p[2]);

									// Check whether transformed point is inside volume
									if ((p[0] > _source_x1[n]) && (p[0] < _source_x2[n]) &&
										(p[1] > _source_y1[n]) && (p[1] < _source_y2[n]) &&
										((_target[n]->GetZ() == 1 && round(p[2]) == 0)
										||( _target[n]->GetZ() != 1 
										&& p[2] > _source_z1[n] && p[2] < _source_z2[n]))) {
											if (*ptr2target > _TargetPadding && weight > 0) {
												// Add sample to metric
												tmpMetricB[n]->Add(*ptr2target, round(_interpolator[n]->EvaluateInside(p[0], p[1], p[2], t)),weight);
												l1b += weight;
												sweight[n] += weight;
											}
									}
							}
						}
						// Increment pointers to next voxel
						ptr2target++;
						ptr2tmp++;
						ptr += 3;
					}
				}
			}
		}
	}
  }

  for (n = 0; n < _numberOfuImages; n++) {
	  utmpMetricA[n]->Reset(_umetric[n]);
	  utmpMetricB[n]->Reset(_umetric[n]);
	  tweight[n] = 0;

    // Calculate bounding box of control point in world coordinates
    _affd->BoundingBox(index, p1, p2);
    _target[0]->WorldToImage(p1);
    _target[0]->WorldToImage(p2);

    // Calculate bounding box of control point in image coordinates
    _affd->MultiBoundingBox(_utarget[n], index, i1, j1, k1, i2, j2, k2, 1.0 / _SpeedupFactor);

    // Calculate incremental changes in lattice coordinates when looping
    // over target
    dx = (FFDLOOKUPTABLESIZE-1)/(p2._x-p1._x);
    dy = (FFDLOOKUPTABLESIZE-1)/(p2._y-p1._y);
    dz = (FFDLOOKUPTABLESIZE-1)/(p2._z-p1._z);

	min = round((FFDLOOKUPTABLESIZE-1)*(0.5 - 0.5/_SpeedupFactor));
	max = round((FFDLOOKUPTABLESIZE-1)*(0.5 + 0.5/_SpeedupFactor));

    // Loop over all voxels in the target (reference) volume
	for (t = 0; t < _utarget[n]->GetT(); t++) {
		for (k = k1; k <= k2; k++) {
			for (j = j1; j <= j2; j++) {
				ptr2utarget = _utarget[n]->GetPointerToVoxels(i1, j, k, t);
				ptr        = &(_uaffdLookupTable[n][3*_utarget[n]->VoxelToIndex(i1, j, k)]);
				ptr2utmp  = _utmpImage[n]->GetPointerToVoxels(i1, j, k, t);
				for (i = i1; i <= i2; i++) {
					x = i; y = j; z = k;
					_utarget[n]->ImageToWorld(x,y,z);
					_target[0]->WorldToImage(x,y,z);
					if(round(dz*(z-p1._z))<=max && round(dz*(z-p1._z))>=min){
						bk = step * _localLookupTable[round(dz*(z-p1._z))];
						if(round(dy*(y-p1._y))<=max && round(dy*(y-p1._y))>=min){
							bj = bk * _localLookupTable[round(dy*(y-p1._y))];
							// Check whether reference point is valid
							if (*ptr2utarget > _TargetPadding && round(dx*(x-p1._x))<=max && round(dx*(x-p1._x))>=min) {
									bi = bj * _localLookupTable[round(dx*(x-p1._x))];
									// Delete old samples from both metrics
									wx = i; wy = j;	wz = k;
									// Convert transformed point to image coordinates
									_utarget[n]->ImageToWorld(wx,wy,wz);
									_threshold->WorldToImage(wx,wy,wz);
									if(wx > 0 && wy > 0 && wz > 0
										&& wx < _threshold->GetX() - 1
										&& wy < _threshold->GetY() - 1
										&& wz < _threshold->GetZ() - 1 ){
											threshold = _threshold->GetAsDouble(round(wx),round(wy),round(wz));
									}else{
										threshold = 0;
									}
                                    wx = i; wy = j;	wz = k;
                                    // Convert transformed point to image coordinates
                                    _utarget[n]->ImageToWorld(wx,wy,wz);
                                    _mweight->WorldToImage(wx,wy,wz);
                                    if(wx > 0 && wy > 0 && wz > 0
                                        && wx < _mweight->GetX() - 1
                                        && wy < _mweight->GetY() - 1
                                        && wz < _mweight->GetZ() - 1 ){
                                            weight = _mweight->GetAsDouble(round(wx),round(wy),round(wz));
                                    }else{
                                        weight = 0;
                                    }
									if (*ptr2utmp != -1 && *ptr2utarget > _TargetPadding && (threshold > 2)) {
										utmpMetricA[n]->Delete(*ptr2utarget, *ptr2utmp, 1.0 - weight);
										utmpMetricB[n]->Delete(*ptr2utarget, *ptr2utmp, 1.0 - weight);
									}

									p[0] = ptr[0];
									p[1] = ptr[1];
									p[2] = ptr[2];
									p[dim] += bi;

									_usource[n]->WorldToImage(p[0], p[1], p[2]);

									// Check whether transformed point is inside volume
									if ((p[0] > _usource_x1[n]) && (p[0] < _usource_x2[n]) &&
										(p[1] > _usource_y1[n]) && (p[1] < _usource_y2[n]) &&
										((_utarget[n]->GetZ() == 1 && round(p[2]) == 0)
										||( _utarget[n]->GetZ() != 1 
										&& p[2] > _usource_z1[n] && p[2] < _usource_z2[n]))) {
											if (*ptr2utarget > _TargetPadding  && (threshold > 2)) {
											// Add sample to metric
											utmpMetricA[n]->Add(*ptr2utarget, round(_uinterpolator[n]->EvaluateInside(p[0], p[1], p[2], t)),1.0 - weight);
											l2a += 1.0 - weight;
											tweight[n] += 1.0 - weight;
											}
									}

									p[0] = ptr[0];
									p[1] = ptr[1];
									p[2] = ptr[2];
									p[dim] -= bi;

									_usource[n]->WorldToImage(p[0], p[1], p[2]);

									// Check whether transformed point is inside volume
									if ((p[0] > _usource_x1[n]) && (p[0] < _usource_x2[n]) &&
										(p[1] > _usource_y1[n]) && (p[1] < _usource_y2[n]) &&
										((_utarget[n]->GetZ() == 1 && round(p[2]) == 0)
										||( _utarget[n]->GetZ() != 1 
										&& p[2] > _usource_z1[n] && p[2] < _usource_z2[n]))) {

											if (*ptr2utarget > _TargetPadding  && (threshold > 2)) {
											// Add sample to metric
											utmpMetricB[n]->Add(*ptr2utarget, round(_uinterpolator[n]->EvaluateInside(p[0], p[1], p[2], t)),1.0 - weight);
											l2b += 1.0 - weight;
											tweight[n] += 1.0 - weight;
											}
									}
							}
						}
					}
					// Increment pointers to next voxel
					ptr2utarget++;
					ptr2utmp++;
					ptr += 3;
				}
			}
		}
	}
  }

  // Save value of DOF for which we calculate the derivative
  double dof = _affd->Get(index);

  // Evaluate similarity measure
  double similarityA =  combine_mysimilarity(combine_mysimilarity(tmpMetricA,sweight,_numberOfImages)
	  ,combine_mysimilarity(utmpMetricA,tweight,_numberOfuImages),l1a+l1b,l2a+l2b);

  // Add penalties
  _affd->Put(index, dof + step);

   // Add penalty for landmark regulation
  if (this->_Lregu > 0) {
	  similarityA += this->_Lregu * this->LandMarkPenalty(index);
  }
  // Smoothness
  if (this->_Lambda1 > 0) {
    similarityA += this->_Lambda1*this->SmoothnessPenalty(index);
  }
  // Volume preservation
  if (this->_Lambda2 > 0) {
    similarityA += this->_Lambda2*this->VolumePreservationPenalty(index);
  }
  // Topology preservation
  if (this->_Lambda3 > 0) {
    similarityA += this->_Lambda3*this->TopologyPreservationPenalty(index);
  }

  // Evaluate similarity measure
  double similarityB = combine_mysimilarity(combine_mysimilarity(tmpMetricB,sweight,_numberOfImages)
	  ,combine_mysimilarity(utmpMetricB,tweight,_numberOfuImages),l1a+l1b,l2a+l2b);

  // Add penalties
  _affd->Put(index, dof - step);

  // Add penalty for landmark regulation
  if (this->_Lregu > 0) {
	  similarityB += this->_Lregu * this->LandMarkPenalty(index);
  }
  // Smoothness
  if (this->_Lambda1 > 0) {
    similarityB += this->_Lambda1*this->SmoothnessPenalty(index);
  }
  // Volume preservation
  if (this->_Lambda2 > 0) {
    similarityB += this->_Lambda2*this->VolumePreservationPenalty(index);
  }
  // Topology preservation
  if (this->_Lambda3 > 0) {
    similarityB += this->_Lambda3*this->TopologyPreservationPenalty(index);
  }

  // Restore value of DOF for which we calculate the derivative
  _affd->Put(index, dof);

#ifdef HAS_TBB
  queue.push(tmpMetricA);
  queue.push(tmpMetricB);
#endif

  delete []sweight;
  delete []tweight;

  return similarityA - similarityB;
}

double irtkCardiac3DImageFreeFormRegistration::EvaluateGradient(float step, float *dx)
{
  int i;
  double norm;

  // Update lookup table
  this->UpdateLUT();

#ifdef HAS_TBB
  parallel_for(blocked_range<int>(0, _affd->NumberOfDOFs(), 1), irtkMultiThreadedImageFreeFormRegistrationEvaluateGradient(this, dx, step));
#else
  for (i = 0; i < _affd->NumberOfDOFs(); i++) {
    if (_affd->irtkTransformation::GetStatus(i) == _Active) {
      dx[i] = this->EvaluateDerivative(i, step);
    } else {
      dx[i] = 0;
    }
  }
#endif

  // Calculate norm of vector
  norm = 0;
  for (i = 0; i < _affd->NumberOfDOFs(); i++) {
    norm += dx[i] * dx[i];
  }

  // Normalize vector
  norm = sqrt(norm);
  if (norm > 0) {
    for (i = 0; i < _affd->NumberOfDOFs(); i++) {
      dx[i] /= norm;
    }
  } else {
    for (i = 0; i < _affd->NumberOfDOFs(); i++) {
      dx[i] = 0;
    }
  }

  return norm;
}

void irtkCardiac3DImageFreeFormRegistration::GuessParameter()
{
  int i;
  double xsize, ysize, zsize;

  this->irtkMultipleImageFreeFormRegistration::GuessParameter();

  _utarget[0]->GetPixelSize(&xsize, &ysize, &zsize);

  // Default target parameters
  _uTargetBlurring[0]      = GuessResolution(xsize, ysize, zsize) / 2.0;
  _uTargetResolution[0][0] = GuessResolution(xsize, ysize, zsize);
  _uTargetResolution[0][1] = GuessResolution(xsize, ysize, zsize);
  _uTargetResolution[0][2] = GuessResolution(xsize, ysize, zsize);

  for (i = 1; i < _NumberOfLevels; i++) {
    _uTargetBlurring[i]      = _uTargetBlurring[i-1] * 2;
    _uTargetResolution[i][0] = _uTargetResolution[i-1][0] * 2;
    _uTargetResolution[i][1] = _uTargetResolution[i-1][1] * 2;
    _uTargetResolution[i][2] = _uTargetResolution[i-1][2] * 2;
  }

  // Read source pixel size
  _usource[0]->GetPixelSize(&xsize, &ysize, &zsize);

  // Default source parameters
  _uSourceBlurring[0]      = GuessResolution(xsize, ysize, zsize) / 2.0;
  _uSourceResolution[0][0] = GuessResolution(xsize, ysize, zsize);
  _uSourceResolution[0][1] = GuessResolution(xsize, ysize, zsize);
  _uSourceResolution[0][2] = GuessResolution(xsize, ysize, zsize);

  for (i = 1; i < _NumberOfLevels; i++) {
    _uSourceBlurring[i]      = _uSourceBlurring[i-1] * 2;
    _uSourceResolution[i][0] = _uSourceResolution[i-1][0] * 2;
    _uSourceResolution[i][1] = _uSourceResolution[i-1][1] * 2;
    _uSourceResolution[i][2] = _uSourceResolution[i-1][2] * 2;
  }
}

bool irtkCardiac3DImageFreeFormRegistration::Read(char *buffer1, char *buffer2, int &level)
{
  int ok = false;
  int i,n;
  double dx,dy,dz;

  // uTarget blurring
  if (strstr(buffer1, "uTarget blurring (in mm)") != NULL) {
    if (level == -1) {
      for (i = 0; i < MAX_NO_RESOLUTIONS; i++) {
        this->_uTargetBlurring[i] = pow(2.0, double(i)) * atof(buffer2);
      }
    } else {
      this->_uTargetBlurring[level] = atof(buffer2);
    }
    ok = true;
  }
  // uTarget resolution
  if (strstr(buffer1, "uTarget resolution (in mm)") != NULL) {
    _utarget[0]->GetPixelSize(&dx, &dy, &dz);
    if (level == -1) {
      n = sscanf(buffer2, "%lf %lf %lf", &(this->_uTargetResolution[0][0]),  &(this->_uTargetResolution[0][1]),  &(this->_uTargetResolution[0][2]));
      if (n == 1) {
        this->_uTargetResolution[0][1] = this->_uTargetResolution[0][0];
        this->_uTargetResolution[0][2] = this->_uTargetResolution[0][0];
      }
      if (this->_uTargetResolution[0][0] == 0) this->_uTargetResolution[0][0] = dx;
      if (this->_uTargetResolution[0][1] == 0) this->_uTargetResolution[0][1] = dy;
      if (this->_uTargetResolution[0][2] == 0) this->_uTargetResolution[0][2] = dz;
      for (i = 0; i < MAX_NO_RESOLUTIONS; i++) {
        this->_uTargetResolution[i][0] = pow(2.0, double(i)) * this->_uTargetResolution[0][0];
        this->_uTargetResolution[i][1] = pow(2.0, double(i)) * this->_uTargetResolution[0][1];
        this->_uTargetResolution[i][2] = pow(2.0, double(i)) * this->_uTargetResolution[0][2];
      }
    } else {
      n = sscanf(buffer2, "%lf %lf %lf", &(this->_uTargetResolution[level][0]),  &(this->_uTargetResolution[level][1]),  &(this->_uTargetResolution[level][2]));
      if (n == 1) {
        this->_uTargetResolution[level][1] = this->_uTargetResolution[level][0];
        this->_uTargetResolution[level][2] = this->_uTargetResolution[level][0];
      }
      if (this->_uTargetResolution[level][0] == 0) this->_uTargetResolution[level][0] = dx;
      if (this->_uTargetResolution[level][1] == 0) this->_uTargetResolution[level][1] = dy;
      if (this->_uTargetResolution[level][2] == 0) this->_uTargetResolution[level][2] = dz;
    }
    ok = true;
  }
  // uSource blurring
  if (strstr(buffer1, "uSource blurring (in mm)") != NULL) {
    if (level == -1) {
      for (i = 0; i < MAX_NO_RESOLUTIONS; i++) {
        this->_uSourceBlurring[i] = pow(2.0, double(i)) * atof(buffer2);
      }
    } else {
      this->_uSourceBlurring[level] = atof(buffer2);
    }
    ok = true;
  }
  // uSource resolution
  if (strstr(buffer1, "uSource resolution (in mm)") != NULL) {
    _usource[0]->GetPixelSize(&dx, &dy, &dz);
    if (level == -1) {
      n = sscanf(buffer2, "%lf %lf %lf", &(this->_uSourceResolution[0][0]),  &(this->_uSourceResolution[0][1]),  &(this->_uSourceResolution[0][2]));
      if (n == 1) {
        this->_uSourceResolution[0][1] = this->_uSourceResolution[0][0];
        this->_uSourceResolution[0][2] = this->_uSourceResolution[0][0];
      }
      if (this->_uSourceResolution[0][0] == 0) this->_uSourceResolution[0][0] = dx;
      if (this->_uSourceResolution[0][1] == 0) this->_uSourceResolution[0][1] = dy;
      if (this->_uSourceResolution[0][2] == 0) this->_uSourceResolution[0][2] = dz;
      for (i = 0; i < MAX_NO_RESOLUTIONS; i++) {
        this->_uSourceResolution[i][0] = pow(2.0, double(i)) * this->_uSourceResolution[0][0];
        this->_uSourceResolution[i][1] = pow(2.0, double(i)) * this->_uSourceResolution[0][1];
        this->_uSourceResolution[i][2] = pow(2.0, double(i)) * this->_uSourceResolution[0][2];
      }
    } else {
      n = sscanf(buffer2, "%lf %lf %lf", &(this->_uSourceResolution[level][0]),  &(this->_uSourceResolution[level][1]),  &(this->_uSourceResolution[level][2]));
      if (n == 1) {
        this->_uSourceResolution[level][1] = this->_uSourceResolution[level][0];
        this->_uSourceResolution[level][2] = this->_uSourceResolution[level][0];
      }
      if (this->_uSourceResolution[level][0] == 0) this->_uSourceResolution[level][0] = dx;
      if (this->_uSourceResolution[level][1] == 0) this->_uSourceResolution[level][1] = dy;
      if (this->_uSourceResolution[level][2] == 0) this->_uSourceResolution[level][2] = dz;
    }
    ok = true;
  }

  if (ok == false) {
    return this->irtkMultipleImageFreeFormRegistration::Read(buffer1, buffer2, level);
  } else {
    return ok;
  }
}

void irtkCardiac3DImageFreeFormRegistration::Write(ostream &to)
{
	int i;
	to << "\n#\n# Cardiac 3d Registration parameters\n#\n\n";
	for (i = 0; i < this->_NumberOfLevels; i++) {
		to << "\n#\n# Cardiac 3d Registration parameters for resolution level " << i+1 << "\n#\n\n";
		to << "Resolution level                  = " << i+1 << endl;
		to << "uTarget blurring (in mm)           = " << this->_uTargetBlurring[i] << endl;
		to << "uTarget resolution (in mm)         = " << this->_uTargetResolution[i][0] << " " << this->_uTargetResolution[i][1] << " " << this->_uTargetResolution[i][2] << endl;
		to << "uSource blurring (in mm)           = " << this->_uSourceBlurring[i] << endl;
		to << "uSource resolution (in mm)         = " << this->_uSourceResolution[i][0] << " " << this->_uSourceResolution[i][1] << " " << this->_uSourceResolution[i][2] << endl;
	}

	this->irtkMultipleImageFreeFormRegistration::Write(to);
}

double irtkCardiac3DImageFreeFormRegistration::WeightFunction (double edge,double edgemax,double threshold){

	double weight = 0;
	weight = fabs(edge*threshold)/edgemax;
	if(weight > 1)
		weight = 1;
	if(weight < 0)
		weight = 0;
	return weight;
}

void irtkCardiac3DImageFreeFormRegistration::EvaluateWeight (irtkRealImage *weight, irtkGreyImage *target, irtkGreyImage *threshold)
{
	int i,j,k,t;
	double xsize,ysize,zsize;
	irtkRealPixel *ptr2weight,*ptr2tmp;
	irtkGreyPixel *ptr2threshold,*ptr2target,*ptr2edge,*edgemax,*edgemin;

	//Innitialize buffers
	irtkGreyImage *tmpedge = new irtkGreyImage(target->GetX(),
		target->GetY(),
		target->GetZ(),
		target->GetT());
	irtkRealImage *tmpthresholdedge = new irtkRealImage(threshold->GetX(),
		threshold->GetY(),
		threshold->GetZ(),
		threshold->GetT());
	ptr2edge = tmpedge->GetPointerToVoxels();
	ptr2tmp = tmpthresholdedge->GetPointerToVoxels();
	for (t = 0; t < target->GetT(); t++) {
		for (k = 0; k < target->GetZ(); k++) {
			for (j = 0; j < target->GetY(); j++) {
				for (i = 0; i < target->GetX(); i++) {
					*ptr2edge = 0;
				}
			}
		}
	}
	for (t = 0; t < threshold->GetT(); t++) {
		for (k = 0; k < threshold->GetZ(); k++) {
			for (j = 0; j < threshold->GetY(); j++) {
				for (i = 0; i < threshold->GetX(); i++) {
					*ptr2tmp = 0;
				}
			}
		}
	}
	irtkRealImage *tmpthreshold = new irtkRealImage(*threshold);

	//evaluate edges
	ptr2threshold = threshold->GetPointerToVoxels();
	ptr2tmp = tmpthreshold->GetPointerToVoxels();
	for (t = 0; t < threshold->GetT(); t++) {
		for (k = 0; k < threshold->GetZ(); k++) {
			for (j = 0; j < threshold->GetY(); j++) {
				for (i = 0; i < threshold->GetX(); i++) {
					if (*ptr2threshold == 3)
						*ptr2tmp = 8;
					else
						*ptr2tmp = 0;
					ptr2threshold++;
					ptr2tmp++;

				}
			}
		}
	}

	edgemax = new irtkGreyPixel();
	edgemin = new irtkGreyPixel();

	// Compute gradient
	irtkGradientImage<irtkRealPixel> tmpgradient;
	tmpgradient.SetPadding(-1);
	tmpgradient.SetInput(tmpthreshold);
	tmpgradient.SetOutput(tmpthresholdedge);
	tmpgradient.Run();

	//tmpthresholdedge->Write("1.nii");

	cout << "Blurring Edge probalility ... "; cout.flush();
	tmpthresholdedge->GetPixelSize(&xsize, &ysize, &zsize);
	irtkGaussianBlurringWithPadding<irtkRealPixel> blurring(sqrt(xsize*xsize+ysize*ysize+zsize*zsize), -1);
	blurring.SetInput (tmpthresholdedge);
	blurring.SetOutput(tmpthresholdedge);
	blurring.Run();
	cout << "done" << endl;

	//tmpthresholdedge->Write("2.nii");
	
	irtkGradientImage<irtkGreyPixel> edgegradient;
	edgegradient.SetPadding(-1);
	edgegradient.SetInput(target);
	edgegradient.SetOutput(tmpedge);
	edgegradient.Run();

	tmpedge->GetMinMax(edgemin,edgemax);

	//tmpedge->Write("3.nii");

	//caculate the weight
	ptr2weight = weight->GetPointerToVoxels();
	ptr2target = target->GetPointerToVoxels();
	ptr2edge = tmpedge->GetPointerToVoxels();
	ptr2tmp = tmpthresholdedge->GetPointerToVoxels();
	for (t = 0; t < target->GetT(); t++) {
		for (k = 0; k < target->GetZ(); k++) {
			for (j = 0; j < target->GetY(); j++) {
				for (i = 0; i < target->GetX(); i++) {
					if (*ptr2target >= 0 && *ptr2edge > 0 && *ptr2tmp > 0)
						*ptr2weight = WeightFunction(*ptr2edge,*edgemax,*ptr2tmp);
					else
						*ptr2weight = 0;

					ptr2weight++;
					ptr2target++;
					ptr2edge++;
					ptr2tmp++;
				}
			}
		}
	}
	blurring.SetSigma(sqrt(xsize*xsize+ysize*ysize));
	blurring.SetInput (weight);
	blurring.SetOutput(weight);
	blurring.Run();
	delete tmpedge;
	delete tmpthresholdedge;
	delete tmpthreshold;
	delete edgemax;
	delete edgemin;
}

void irtkCardiac3DImageFreeFormRegistration::EvaluateMyoProb1 (irtkRealImage *weight, irtkGreyImage *threshold, irtkGaussian &gaussian, double &denom)
{
	int i,j,k,t;
	double mean,sigma;
	mean = 0; sigma = 0; denom = 0;
	for (t = 0; t < threshold->GetT(); t++) {
		for (k = 0; k < threshold->GetZ(); k++) {
			for (j = 0; j < threshold->GetY(); j++) {
				for (i = 0; i < threshold->GetX(); i++) {
					if(threshold->GetAsDouble(i,j,k,t) == 3){
						mean += _target[0]->GetAsDouble(i,j,k,t);
						denom++;
					}
				}
			}
		}
	}
	mean = mean/denom;
	for (t = 0; t < threshold->GetT(); t++) {
		for (k = 0; k < threshold->GetZ(); k++) {
			for (j = 0; j < threshold->GetY(); j++) {
				for (i = 0; i < threshold->GetX(); i++) {
					if(threshold->GetAsDouble(i,j,k,t) == 3){
						sigma += pow((_target[0]->GetAsDouble(i,j,k,t) - mean),2);
					}
				}
			}
		}
	}
	sigma = sigma/denom;

	gaussian.Initialise(mean,sigma);

	denom = gaussian.Evaluate(mean);
}

void irtkCardiac3DImageFreeFormRegistration::EvaluateMyoProb2 (irtkRealImage *weight, irtkGreyImage *threshold, irtkGaussian &gaussian, double &denom)
{
	int i,j,k,t;

	for (t = 0; t < threshold->GetT(); t++) {
		for (k = 0; k < threshold->GetZ(); k++) {
			for (j = 0; j < threshold->GetY(); j++) {
				for (i = 0; i < threshold->GetX(); i++) {
					if(threshold->GetAsDouble(i,j,k,t) == 3){
						weight->PutAsDouble(i,j,k,t,gaussian.Evaluate(_target[0]->GetAsDouble(i,j,k,t))/denom);
					}else{
						weight->PutAsDouble(i,j,k,t,0);
					}
				}
			}
		}
	}

}

void irtkCardiac3DImageFreeFormRegistration::EvaluateOmega (){
	int i;
	for(i=0;i<_affd->NumberOfDOFs();i++)
		this->EvaluateOmega(i);
}

void irtkCardiac3DImageFreeFormRegistration::EvaluateOmega (int index){
	irtkPoint p1, p2;
	int i, j, k, i1, i2, j1, j2, k1, k2, dim, t, min,max;
	double bi, bj, bk, dx, dy, dz, x,y,z;
	irtkRealPixel *ptr;

	dim = int(index / (_affd->GetX()*_affd->GetY()*_affd->GetZ()));
	_comega[index] = 0;


	// Calculate bounding box of control point in world coordinates
	_affd->BoundingBox(index, p1, p2);
	_target[0]->WorldToImage(p1);
	_target[0]->WorldToImage(p2);
	// Calculate bounding box of control point in image coordinates
	_affd->MultiBoundingBox(_threshold, index, i1, j1, k1, i2, j2, k2, 1.0 / _SpeedupFactor);

	// Calculate incremental changes in lattice coordinates when looping
	// over target
	dx = (FFDLOOKUPTABLESIZE-1)/(p2._x-p1._x);
	dy = (FFDLOOKUPTABLESIZE-1)/(p2._y-p1._y);
	dz = (FFDLOOKUPTABLESIZE-1)/(p2._z-p1._z);

	min = round((FFDLOOKUPTABLESIZE-1)*(0.5 - 0.5/_SpeedupFactor));
	max = round((FFDLOOKUPTABLESIZE-1)*(0.5 + 0.5/_SpeedupFactor));

	// Loop over all voxels in the target (reference) volume
	for (t = 0; t < _omega->GetT(); t++) {
		for (k = k1; k <= k2; k++) {
			for (j = j1; j <= j2; j++) {
				ptr = _omega->GetPointerToVoxels(i1, j, k, t);
				for (i = i1; i <= i2; i++) {
					if(*ptr > 0){
						x = i; y = j; z = k;
						_omega->ImageToWorld(x,y,z);
						_target[0]->WorldToImage(x,y,z);
						if(round(dz*(z-p1._z))<=max && round(dz*(z-p1._z))>=min){
							bk = _localLookupTable[round(dz*(z-p1._z))];
							if(round(dy*(y-p1._y))<=max && round(dy*(y-p1._y))>=min){
								bj = bk * _localLookupTable[round(dy*(y-p1._y))];
								// Check whether reference point is valid
								if (*ptr > 0 && round(dx*(x-p1._x))<=max && round(dx*(x-p1._x))>=min) {
									bi = bj * _localLookupTable[round(dx*(x-p1._x))];
									_comega[index] += bi*(*ptr);
								}
							}
						}
					}

					// Increment pointers to next voxel
					ptr ++;
			}
		}
	}
}
}

void irtkCardiac3DImageFreeFormRegistration::UpdateOmega ()
{
	int i, j, k;
	double x, y, z, jacobian, count;
	count = 0;
	irtkMatrix jac,tmp_jac;
	/*for (k = 0; k < _omega->GetZ(); k++) {
		for (j = 0; j < _omega->GetY(); j++) {
			for (i = 0; i < _omega->GetX(); i++) {
				if(_myoprob->GetAsDouble(i,j,k) > 0.01){
					_omega->PutAsDouble(i,j,k,_myoprob->GetAsDouble(i,j,k));
				}else{
					_omega->PutAsDouble(i,j,k,0.01);
				}
				count += _omega->GetAsDouble(i,j,k);
			}
		}
	}*/
	for (k = 0; k < _omega->GetZ(); k++) {
		for (j = 0; j < _omega->GetY(); j++) {
			for (i = 0; i < _omega->GetX(); i++) {
				x = i;
				y = j;
				z = k;
				_omega->ImageToWorld(x, y, z);
				// myocardium adaptive jacobian weight
				_affd->Jacobian(tmp_jac,x,y,z);
				// Calculate jacobian
				jac.Initialize(3, 3);
				_mffd->LocalJacobian(jac, x, y, z);

				// Subtract identity matrix
				tmp_jac(0, 0) = tmp_jac(0, 0) - 1;
				tmp_jac(1, 1) = tmp_jac(1, 1) - 1;
				tmp_jac(2, 2) = tmp_jac(2, 2) - 1;

				// Add jacobian
				jac += tmp_jac;
				// Determinant of Jacobian of deformation derivatives
				jacobian = (jac(0, 0)*jac(1, 1)*jac(2, 2) + jac(0, 1)*jac(1, 2)*jac(2, 0) +
					jac(0, 2)*jac(1, 0)*jac(2, 1) - jac(0, 2)*jac(1, 1)*jac(2, 0) -
					jac(0, 0)*jac(1, 2)*jac(2, 1) - jac(0, 1)*jac(1, 0)*jac(2, 2));

				if(jacobian < 0.0000001) 
					jacobian = 0.0000001;

				if(_myoprob->GetAsDouble(i,j,k) > 0){
					_omega->PutAsDouble(i,j,k,
						fabs(log(jacobian))
						/(fabs(log(_myoprob->GetAsDouble(i,j,k)))+0.3)
						);
				}else{
					_omega->PutAsDouble(i,j,k,fabs(log(jacobian))/100.0);
				}
				if(_omega->GetAsDouble(i,j,k)<0.001)
					_omega->PutAsDouble(i,j,k,0);
				count += _omega->GetAsDouble(i,j,k);
			}
		}
	}	
	//normalize
	if(count>0){
		for (k = 0; k < _omega->GetZ(); k++) {
			for (j = 0; j < _omega->GetY(); j++) {
				for (i = 0; i < _omega->GetX(); i++) {
					_omega->PutAsDouble(i,j,k,_omega->GetAsDouble(i,j,k)/count);
				}
			}
		}
	}
}

void irtkCardiac3DImageFreeFormRegistration::ExtendThreshold (irtkGreyImage *threshold, int n)
{
	int i,j,k,t;
	int a,b,c;
	double x,y,z;
	irtkGreyImage tmp(threshold->GetImageAttributes());
		for (k = 0; k < _utarget[0]->GetZ(); k++) {
			for (j = 0; j < _utarget[0]->GetY(); j++) {
				for (i = 0; i < _utarget[0]->GetX(); i++) {
					x = i; y = j; z = k;
					_utarget[0]->ImageToWorld(x,y,z);
					threshold->WorldToImage(x,y,z);
					if(round(x) > 0 && round(x) < threshold->GetX()
						&& round(y) > 0 && round(y) < threshold->GetY()
						&& round(z) > 0 && round(z) < threshold->GetZ()
						&& threshold->GetAsDouble(round(x),round(y),round(z)) == 3){
						for(a=i-n;a<i+n+1;a++){
							for(b=j-n;b<j+n+1;b++){
								for(c=k;c<k+2;c++){
									x = a; y = b; z = c;
									_utarget[0]->ImageToWorld(x,y,z);
									threshold->WorldToImage(x,y,z);
									if(round(x) > 0 && round(x) < threshold->GetX()
										&& round(y) > 0 && round(y) < threshold->GetY()
										&& round(z) > 0 && round(z) < threshold->GetZ()
										&& threshold->GetAsDouble(round(x),round(y),round(z)) < 3
										){
										tmp.PutAsDouble(round(x),round(y),round(z),3);
									}
								}
							}
						}

					}
				}
			}
		}
	for (t = 0; t < threshold->GetT(); t++) {
		for (k = 0; k < threshold->GetZ(); k++) {
			for (j = 0; j < threshold->GetY(); j++) {
				for (i = 0; i < threshold->GetX(); i++) {
					if(tmp.GetAsDouble(i,j,k,t) == 3){
						threshold->PutAsDouble(i,j,k,t,3);
					}
				}
			}
		}
	}
}

void irtkCardiac3DImageFreeFormRegistration::Run()
{
  int i, j, level;
  char buffer[256];
  double step, epsilon, delta, maxChange;

  // Print debugging information
  this->Debug("irtkMultipleImageRegistration::Run");

  if (_source == NULL) {
    cerr << "Registration::Run: Filter has no source input" << endl;
    exit(1);
  }

  if (_target == NULL) {
    cerr << "Registration::Run: Filter has no target input" << endl;
    exit(1);
  }

  if (_transformation == NULL) {
    cerr << "irtkMultipleImageRegistration::Run: Filter has no transformation output" << endl;
    exit(1);
  }

  // Do the initial set up for all levels
  this->Initialize();

  // Loop over levels
  for (level = _NumberOfLevels-1; level >= 0; level--) {


    // Initial step size
    step = _LengthOfSteps[level];

    // Print resolution level
    cout << "Resolution level no. " << level+1 << " (step sizes ";
    cout << step << " to " << step / pow(2.0, static_cast<double>(_NumberOfSteps[level]-1)) << ")\n";

    // Initial Delta
    delta = _Delta[level];
    cout << "Delta values : " << delta << " to ";
    cout << delta / pow(2.0, static_cast<double>(_NumberOfSteps[level]-1)) << "\n";

#ifdef HISTORY
    history->Clear();
#endif

    // Initialize for this level
    this->Initialize(level);

#ifdef HAS_TBB
    task_scheduler_init init(tbb_no_threads);

    tick_count t_start = tick_count::now();
#endif

    // Run the registration filter at this resolution
    for (i = 0; i < _NumberOfSteps[level]; i++) {
		//Update Omega
		//&& (level == 0 || level < _NumberOfLevels-1)
		if (this->_Lambda2 > 0) {
			this->UpdateOmega();
			this->EvaluateOmega();
			_omega->Write("omega.nii.gz");
		}
      for (j = 0; j < _NumberOfIterations[level]; j++) {

        cout << "Iteration = " << j + 1 << " (out of " << _NumberOfIterations[level];
        cout << "), step size = " << step << endl;

        // Optimize at lowest level of resolution
        _optimizer->SetStepSize(step);
        _optimizer->SetEpsilon(_Epsilon);
        _optimizer->Run(epsilon, maxChange);

        // Check whether we made any improvement or not
        if (epsilon > _Epsilon && maxChange > delta) {
          sprintf(buffer, "log_%.3d_%.3d_%.3d.dof", level, i+1, j+1);
          if (_DebugFlag == true) _transformation->Write(buffer);
          this->Print();
        } else {
          sprintf(buffer, "log_%.3d_%.3d_%.3d.dof", level, i+1, j+1);
          if (_DebugFlag == true) _transformation->Write(buffer);
          this->Print();
          break;
        }
      }
      step = step / 2;
      delta = delta / 2.0;
    }

#ifdef HAS_TBB

    tick_count t_end = tick_count::now();
    if (tbb_debug) cout << this->NameOfClass() << " = " << (t_end - t_start).seconds() << " secs." << endl;
    init.terminate();

#endif
    // Do the final cleaning up for this level
    this->Finalize(level);

#ifdef HISTORY
    history->Print();
#endif

  }

  // Do the final cleaning up for all levels
  this->Finalize();
}
