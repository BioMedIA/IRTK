/*=========================================================================

  Library   : Image Registration Toolkit (IRTK)
  Module    : $Id$
  Copyright : Imperial College, Department of Computing
              Visual Information Processing (VIP), 2008 onwards
  Date      : $Date$
  Version   : $Revision$
  Changes   : $Author$

=========================================================================*/

#include <irtkTransformation.h>

#ifdef HAS_TBB

#include <tbb/tick_count.h>

class irtkMultiThreadedImageTransformation
{

  /// Time frame to transform
  int _toutput, _tinput;

  /// Pointer to image transformation class
  irtkImageTransformation *_imagetransformation;

public:

  irtkMultiThreadedImageTransformation(irtkImageTransformation *imagetransformation,  int toutput, int tinput) {
    _toutput = toutput;
    _tinput  = tinput;
    _imagetransformation = imagetransformation;
  }

  void operator()(const blocked_range<int> &r) const {
    int i, j, k;
    double x, y, z, time;

    time = _imagetransformation->_output->ImageToTime(_toutput);

    for (k = r.begin(); k != r.end(); k++) {

      for (j = 0; j < _imagetransformation->_output->GetY(); j++) {
        for (i = 0; i < _imagetransformation->_output->GetX(); i++) {
          if (_imagetransformation->_output->GetAsDouble(i, j, k, _toutput) > _imagetransformation->_TargetPaddingValue) {
            x = i;
            y = j;
            z = k;
            // Transform point into world coordinates
            _imagetransformation->_output->ImageToWorld(x, y, z);
            // Transform point
            if (_imagetransformation->_Invert == False) {
              _imagetransformation->_transformation->Transform(x, y, z, time);
            } else {
              _imagetransformation->_transformation->Inverse(x, y, z, time);
            }
            // Transform point into image coordinates
            _imagetransformation->_input->WorldToImage(x, y, z);
            // Check whether transformed point is in FOV of input
            if ((x > -0.5) && (x < _imagetransformation->_input->GetX()-0.5) &&
                (y > -0.5) && (y < _imagetransformation->_input->GetY()-0.5) &&
                (z > -0.5) && (z < _imagetransformation->_input->GetZ()-0.5)) {
            	_imagetransformation->_output->PutAsDouble(i, j, k, _toutput, _imagetransformation->_ScaleFactor * _imagetransformation->_interpolator->Evaluate(x, y, z, _tinput) + _imagetransformation->_Offset);
            } else {
              // Fill with padding value
            	_imagetransformation->_output->PutAsDouble(i, j, k, _toutput, _imagetransformation->_SourcePaddingValue);
            }
          } else {
            // Fill with padding value
          	_imagetransformation->_output->PutAsDouble(i, j, k, _toutput, _imagetransformation->_SourcePaddingValue);
          }
        }
      }
    }
  }
};

#endif

irtkImageTransformation::irtkImageTransformation()
{
  // Set input and output
  _input  = NULL;
  _output = NULL;

  // Set transformation
  _transformation = NULL;

  // Set interpolator
  _interpolator = NULL;

  // Set padding value
  _TargetPaddingValue = voxel_limits<double>::min();

  // Set padding value
  _SourcePaddingValue = 0;

  // Scale factor
  _ScaleFactor = 1;
  
  // Offset
  _Offset = 0;
  
  // Set invert mode
  _Invert = False;
}

irtkImageTransformation::~irtkImageTransformation()
{
  // Set input and output
  _input  = NULL;
  _output = NULL;

  // Set transformation
  _transformation = NULL;

  // Set interpolator
  _interpolator = NULL;

  // Set padding value
  _TargetPaddingValue = -std::numeric_limits<double>::max();

  // Set padding value
  _SourcePaddingValue = 0;

  // Set invert mode
  _Invert = False;
}

irtkImageTransformation *irtkImageTransformation::New(irtkTransformation *transformation)
{
  irtkImageTransformation *imagetransformation;

  if (strcmp(transformation->NameOfClass(), "irtkHomogeneousTransformation") == 0) {
    imagetransformation = new irtkImageHomogeneousTransformation;
    imagetransformation->SetTransformation(transformation);
    return imagetransformation;
  }
  if (strcmp(transformation->NameOfClass(), "irtkRigidTransformation") == 0) {
    imagetransformation = new irtkImageHomogeneousTransformation;
    imagetransformation->SetTransformation(transformation);
    return imagetransformation;
  }
  if (strcmp(transformation->NameOfClass(), "irtkAffineTransformation") == 0) {
    imagetransformation = new irtkImageHomogeneousTransformation;
    imagetransformation->SetTransformation(transformation);
    return imagetransformation;
  }
  imagetransformation = new irtkImageTransformation;
  imagetransformation->SetTransformation(transformation);
  return imagetransformation;
}

void irtkImageTransformation::SetTransformation(irtkTransformation *transformation)
{
  if (transformation != NULL) {
    _transformation = transformation;
  } else {
    cerr << "irtkImageTransformation::SetInput: Transformation is NULL\n";
    exit(1);
  }
}

void irtkImageTransformation::SetInput(irtkImage *image)
{
  if (image != NULL) {
    _input = image;
  } else {
    cerr << "irtkImageTransformation::SetInput: Input is NULL\n";
    exit(1);
  }
}

void irtkImageTransformation::SetInput(irtkImage *image, irtkTransformation *transformation)
{
  if ((image != NULL) && (transformation != NULL)) {
    _input = image;
    _transformation = transformation;
  } else {
    cerr << "irtkImageTransformation::SetInput: Input is NULL\n";
    exit(1);
  }
}

void irtkImageTransformation::SetOutput(irtkImage *image)
{
  if (image != NULL) {
    _output = image;
  } else {
    cerr << "irtkImageTransformation::SetOutput: Output is NULL\n";
    exit(1);
  }
}

void irtkImageTransformation::Run()
{
  int i, j, k, l;

#ifdef HAS_TBB
  double t;
#else
  double x, y, z, t;
#endif

  // Check inputs and outputs
  if (_input == NULL) {
    cerr << "irtkImageTransformation::Run: Filter has no input" << endl;
    exit(1);
  }

  if (_output == NULL) {
    cerr << "irtkImageTransformation::Run: Filter has no output" << endl;
    exit(1);
  }

  if (_transformation == NULL) {
    cerr << "irtkImageTransformation::Run: Filter has no transformation" << endl;
    exit(1);
  }

  if (_interpolator == NULL) {
    cerr << "irtkImageTransformation::Run: Filter has no interpolator" << endl;
    exit(1);
  }

  if (_input->IsEmpty() == True) {
    cerr << "irtkImageTransformation::Run: Input is empty" << endl;
    exit(1);
  }

  if (_input == _output) {
    cerr << "irtkImageTransformation::Run: Input equals output" << endl;
    exit(1);
  }

  // Setup interpolation
  _interpolator->SetInput(_input);
  _interpolator->Initialize();

#ifdef HAS_TBB
  task_scheduler_init init(tbb_no_threads);

  tick_count t_start = tick_count::now();
#endif

  // Calculate transformation
  for (l = 0; l < _output->GetT(); l++) {
    t = round(this->_input->TimeToImage(this->_output->ImageToTime(l)));

    if ((t >= 0) || (t < this->_input->GetT())) {

#ifdef HAS_TBB
      parallel_for(blocked_range<int>(0, _output->GetZ(), 1), irtkMultiThreadedImageTransformation(this, l, t));
#else

      double time = this->_output->ImageToTime(l);

      for (k = 0; k < _output->GetZ(); k++) {
        for (j = 0; j < _output->GetY(); j++) {
          for (i = 0; i < _output->GetX(); i++) {
            if (this->_output->GetAsDouble(i, j, k, l) > _TargetPaddingValue) {
              x = i;
              y = j;
              z = k;
              // Transform point into world coordinates
              _output->ImageToWorld(x, y, z);
              // Transform point
              if (_Invert == False) {
                _transformation->Transform(x, y, z, time);
              } else {
                _transformation->Inverse(x, y, z, time);
              }
              // Transform point into image coordinates
              _input->WorldToImage(x, y, z);
              // Check whether transformed point is in FOV of input
              if ((x > -0.5) && (x < _input->GetX()-0.5) &&
                  (y > -0.5) && (y < _input->GetY()-0.5) &&
                  (z > -0.5) && (z < _input->GetZ()-0.5)) {
              	this->_output->PutAsDouble(i, j, k, l, _ScaleFactor * _interpolator->Evaluate(x, y, z, t) + _Offset);
              } else {
                // Fill with padding value
              	this->_output->PutAsDouble(i, j, k, l, _SourcePaddingValue);
              }
            } else {
              // Fill with padding value
            	this->_output->PutAsDouble(i, j, k, l, _SourcePaddingValue);
            }
          }
        }
      }

#endif

    } else {
      for (k = 0; k < this->_output->GetZ(); k++) {
        for (j = 0; j < this->_output->GetY(); j++) {
          for (i = 0; i < this->_output->GetX(); i++) {
          	this->_output->PutAsDouble(i, j, k, l, this->_SourcePaddingValue);
          }
        }
      }
    }
  }


#ifdef HAS_TBB

  tick_count t_end = tick_count::now();
  if (tbb_debug) cout << "irtkImageTransformation = " << (t_end - t_start).seconds() << " secs." << endl;
  init.terminate();

#endif

}


