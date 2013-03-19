/*=========================================================================

  Library   : Image Registration Toolkit (IRTK)
  Module    : $Id: irtkImageRegistration.cc 686 2012-10-15 14:43:52Z ws207 $
  Copyright : Imperial College, Department of Computing
              Visual Information Processing (VIP), 2008 onwards
  Date      : $Date: 2012-10-15 15:43:52 +0100 (Mon, 15 Oct 2012) $
  Version   : $Revision: 686 $
  Changes   : $Author: ws207 $

=========================================================================*/

#include <irtkRegistration.h>
#include <irtkImageRegistration.h>
#include <irtkImageRegistrationWithPadding.h>
#include <irtkGradientDescentConstrainedOptimizer.h>


//irtkGreyImage *tmp_target, *tmp_source;

irtkImageRegistrationWithPadding::irtkImageRegistrationWithPadding() : irtkImageRegistration()
{
  _SourcePadding   = MIN_GREY;
}


void irtkImageRegistrationWithPadding::Initialize(int level)
{
  int i, j, k, t;
  double dx, dy, dz, temp;
  irtkGreyPixel target_min, target_max, target_nbins;
  irtkGreyPixel source_min, source_max, source_nbins;

  // Copy source and target to temp space
  tmp_target = new irtkGreyImage(*_target);
  tmp_source = new irtkGreyImage(*_source);

  // Swap source and target with temp space copies
  swap(tmp_target, _target);
  swap(tmp_source, _source);

  // Blur images if necessary
  if (_TargetBlurring[level] > 0) {
    cout << "Blurring target ... "; cout.flush();
    irtkGaussianBlurringWithPadding<irtkGreyPixel> blurring(_TargetBlurring[level], _TargetPadding);
    blurring.SetInput (_target);
    blurring.SetOutput(_target);
    blurring.Run();
    cout << "done" << endl;
  }

  if (_SourceBlurring[level] > 0) {
    cout << "Blurring source ... "; cout.flush();
    irtkGaussianBlurringWithPadding<irtkGreyPixel> blurring(_SourceBlurring[level],_SourcePadding);
    blurring.SetInput (_source);
    blurring.SetOutput(_source);
    blurring.Run();
    cout << "done" << endl;
  }

  _target->GetPixelSize(&dx, &dy, &dz);
  temp = fabs(_TargetResolution[0][0]-dx) + fabs(_TargetResolution[0][1]-dy) + fabs(_TargetResolution[0][2]-dz);

  if (level > 0 || temp > 0.000001) {
    cout << "Resampling target ... "; cout.flush();
    // Create resampling filter
    irtkResamplingWithPadding<irtkGreyPixel> resample(_TargetResolution[level][0],
        _TargetResolution[level][1],
        _TargetResolution[level][2],
        _TargetPadding);
    resample.SetInput (_target);
    resample.SetOutput(_target);
    resample.Run();
    cout << "done" << endl;
  }

  _source->GetPixelSize(&dx, &dy, &dz);
  temp = fabs(_SourceResolution[0][0]-dx) + fabs(_SourceResolution[0][1]-dy) + fabs(_SourceResolution[0][2]-dz);

  if (level > 0 || temp > 0.000001) {
    cout << "Resampling source ... "; cout.flush();
    // Create resampling filter
    irtkResamplingWithPadding<irtkGreyPixel> resample(_SourceResolution[level][0],
        _SourceResolution[level][1],
        _SourceResolution[level][2], _SourcePadding);

    resample.SetInput (_source);
    resample.SetOutput(_source);
    resample.Run();
    cout << "done" << endl;
  }

  // Find out the min and max values in target image, ignoring padding
  target_max = MIN_GREY;
  target_min = MAX_GREY;
  for (t = 0; t < _target->GetT(); t++) {
    for (k = 0; k < _target->GetZ(); k++) {
      for (j = 0; j < _target->GetY(); j++) {
        for (i = 0; i < _target->GetX(); i++) {
          if (_target->Get(i, j, k, t) > _TargetPadding) {
            if (_target->Get(i, j, k, t) > target_max)
              target_max = _target->Get(i, j, k, t);
            if (_target->Get(i, j, k, t) < target_min)
              target_min = _target->Get(i, j, k, t);
          } else {
            _target->Put(i, j, k, t, _TargetPadding);
          }
        }
      }
    }
  }

  // Find out the min and max values in source image, ignoring padding
  source_max = MIN_GREY;
  source_min = MAX_GREY;
  for (t = 0; t < _source->GetT(); t++) {
    for (k = 0; k < _source->GetZ(); k++) {
      for (j = 0; j < _source->GetY(); j++) {
        for (i = 0; i < _source->GetX(); i++) {
          if (_source->Get(i, j, k, t) > _SourcePadding){
            if (_source->Get(i, j, k, t) > source_max)
              source_max = _source->Get(i, j, k, t);
            if (_source->Get(i, j, k, t) < source_min)
              source_min = _source->Get(i, j, k, t);
	  } else {
	    _source->Put(i, j, k, t, _SourcePadding);
	  }
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
    for (t = 0; t < _target->GetT(); t++) {
      for (k = 0; k < _target->GetZ(); k++) {
        for (j = 0; j < _target->GetY(); j++) {
          for (i = 0; i < _target->GetX(); i++) {
            if (_target->Get(i, j, k, t) > _TargetPadding) {
              _target->Put(i, j, k, t, _target->Get(i, j, k, t) - target_min);
            } else {
              _target->Put(i, j, k, t, -1);
            }
          }
        }
      }
    }
  }

    if (source_max - source_min > MAX_GREY) {
      cerr << this->NameOfClass()
           << "::Initialize: Dynamic range of source is too large" << endl;
      exit(1);
    } else {
      for (t = 0; t < _source->GetT(); t++) {
        for (k = 0; k < _source->GetZ(); k++) {
          for (j = 0; j < _source->GetY(); j++) {
            for (i = 0; i < _source->GetX(); i++) {
              if (_source->Get(i, j, k, t) > _SourcePadding) {
                _source->Put(i, j, k, t, _source->Get(i, j, k, t) - source_min);
	      } else {
		_source->Put(i, j, k, t, -1);  
	      }
            }
          }
        }
      }
    }

/*if ((_SimilarityMeasure == SSD) || (_SimilarityMeasure == CC) ||
      (_SimilarityMeasure == LC)  || (_SimilarityMeasure == K) || (_SimilarityMeasure == ML)) {
    if (source_max - target_min > MAX_GREY) {
      cerr << this->NameOfClass()
           << "::Initialize: Dynamic range of source is too large" << endl;
      exit(1);
    } else {
      for (t = 0; t < _source->GetT(); t++) {
        for (k = 0; k < _source->GetZ(); k++) {
          for (j = 0; j < _source->GetY(); j++) {
            for (i = 0; i < _source->GetX(); i++) {
              _source->Put(i, j, k, t, _source->Get(i, j, k, t) - source_min);
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
      for (t = 0; t < _source->GetT(); t++) {
        for (k = 0; k < _source->GetZ(); k++) {
          for (j = 0; j < _source->GetY(); j++) {
            for (i = 0; i < _source->GetX(); i++) {
              _source->Put(i, j, k, t, _source->Get(i, j, k, t) - source_min);
            }
          }
        }
      }
    }
  }
  */

  // Pad target image if necessary
  irtkPadding(*_target, _TargetPadding);

  // Allocate memory for metric
  switch (_SimilarityMeasure) {
  case SSD:
    _metric = new irtkSSDSimilarityMetric;
    break;
  case CC:
    // Rescale images by an integer factor if necessary
    _metric = new irtkCrossCorrelationSimilarityMetric;
    break;
  case JE:
    // Rescale images by an integer factor if necessary
    target_nbins = irtkCalculateNumberOfBins(_target, _NumberOfBins,
                   target_min, target_max);
    source_nbins = irtkCalculateNumberOfBins(_source, _NumberOfBins,
                   source_min, source_max);
    _metric = new irtkJointEntropySimilarityMetric(target_nbins, source_nbins);
    break;
  case MI:
    // Rescale images by an integer factor if necessary
    target_nbins = irtkCalculateNumberOfBins(_target, _NumberOfBins,
                   target_min, target_max);
    source_nbins = irtkCalculateNumberOfBins(_source, _NumberOfBins,
                   source_min, source_max);
    _metric = new irtkMutualInformationSimilarityMetric(target_nbins, source_nbins);
    break;
  case NMI:
    // Rescale images by an integer factor if necessary
    target_nbins = irtkCalculateNumberOfBins(_target, _NumberOfBins,
                   target_min, target_max);
    source_nbins = irtkCalculateNumberOfBins(_source, _NumberOfBins,
                   source_min, source_max);
    _metric = new irtkNormalisedMutualInformationSimilarityMetric(target_nbins, source_nbins);
    break;
  case CR_XY:
    // Rescale images by an integer factor if necessary
    target_nbins = irtkCalculateNumberOfBins(_target, _NumberOfBins,
                   target_min, target_max);
    source_nbins = irtkCalculateNumberOfBins(_source, _NumberOfBins,
                   source_min, source_max);
    _metric = new irtkCorrelationRatioXYSimilarityMetric(target_nbins, source_nbins);
    break;
  case CR_YX:
    // Rescale images by an integer factor if necessary
    target_nbins = irtkCalculateNumberOfBins(_target, _NumberOfBins,
                   target_min, target_max);
    source_nbins = irtkCalculateNumberOfBins(_source, _NumberOfBins,
                   source_min, source_max);
    _metric = new irtkCorrelationRatioYXSimilarityMetric(target_nbins, source_nbins);
    break;
  case LC:
    _metric = new irtkLabelConsistencySimilarityMetric;
    break;
  case K:
    // Rescale images by an integer factor if necessary
    target_nbins = irtkCalculateNumberOfBins(_target, _NumberOfBins,
                   target_min, target_max);
    source_nbins = irtkCalculateNumberOfBins(_source, _NumberOfBins,
                   source_min, source_max);
    _metric = new irtkKappaSimilarityMetric(target_nbins, source_nbins);
    break;
  case ML:
    // Rescale images by an integer factor if necessary
    _metric = new irtkMLSimilarityMetric(classification);
    if (_metric==NULL) {
      cerr<<"Please, do not forget to set the ML metric!!!"<<endl;
    }
    break;
  }

  // Setup the interpolator - currently only linear supported
  _interpolator = irtkInterpolateImageFunction::New(Interpolation_Linear, _source);
  //_interpolator = irtkInterpolateImageFunction::New(_InterpolationMode, _source);

  // Setup interpolation for the source image
  _interpolator->SetInput(_source);
  _interpolator->Initialize();

  // Calculate the source image domain in which we can interpolate
  _interpolator->Inside(_source_x1, _source_y1, _source_z1,
                        _source_x2, _source_y2, _source_z2);

  // Setup the optimizer
  switch (_OptimizationMethod) {
  case DownhillDescent:
    _optimizer = new irtkDownhillDescentOptimizer;
    break;
  case GradientDescent:
    _optimizer = new irtkGradientDescentOptimizer;
    break;
  case GradientDescentConstrained:
    _optimizer = new irtkGradientDescentConstrainedOptimizer;
    break;
  case SteepestGradientDescent:
    _optimizer = new irtkSteepestGradientDescentOptimizer;
    break;
  case ConjugateGradientDescent:
    _optimizer = new irtkConjugateGradientDescentOptimizer;
    break;
  default:
    cerr << "Unkown optimizer" << endl;
    exit(1);
  }
  _optimizer->SetTransformation(_transformation);
  _optimizer->SetRegistration(this);

  // Print some debugging information
  cout << "Target image (reference)" << endl;
  _target->Print();
  cout << "Range is from " << target_min << " to " << target_max << endl;

  cout << "Source image (transform)" << endl;
  _source->Print();
  cout << "Range is from " << source_min << " to " << source_max << endl;

  // Print initial transformation
  cout << "Initial transformation for level = " << level+1 << endl;;
  _transformation->Print();

}
