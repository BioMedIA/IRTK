/*=========================================================================

  Library   : Image Registration Toolkit (IRTK)
  Module    : $Id$
  Copyright : Imperial College, Department of Computing
              Visual Information Processing (VIP), 2009 onwards
  Date      : $Date$
  Version   : $Revision$
  Changes   : $Author$

=========================================================================*/

#include <irtkRegistration.h>

#ifdef HAS_VTK

void irtkModelRigidRegistration::GuessParameter()
{
	int i;
	double xsize, ysize, zsize;

	if (_image == NULL) {
		cerr << "irtkModelRigidRegistration::GuessParameter: Image not found" << endl;
		exit(1);
	}

	// Default parameters for registration
	_NumberOfLevels = 1;

	// Default parameters for optimization
	_OptimizationMethod = GradientDescent;
	_Epsilon = 0.0001;

	// Read target pixel size
	_image->GetPixelSize(&xsize, &ysize, &zsize);

	// Default target parameters
	_ImageBlurring[0] = GuessResolution(xsize, ysize, zsize) / 2.0;
	_ImageResolution[0][0] = GuessResolution(xsize, ysize, zsize);
	_ImageResolution[0][1] = GuessResolution(xsize, ysize, zsize);
	_ImageResolution[0][2] = GuessResolution(xsize, ysize, zsize);

	for (i = 1; i < _NumberOfLevels; i++) {
		_ImageBlurring[i] = _ImageBlurring[i-1] * 2;
		_ImageResolution[i][0] = _ImageResolution[i-1][0] * 2;
		_ImageResolution[i][1] = _ImageResolution[i-1][1] * 2;
		_ImageResolution[i][2] = _ImageResolution[i-1][2] * 2;
	}

	// Remaining parameters
	for (i = 0; i < _NumberOfLevels; i++) {
		_NumberOfIterations[i] = 20;
		_NumberOfSteps[i] = 4;
		_LengthOfSteps[i] = 2 * pow(2.0, i);
	}

}

void irtkModelRigidRegistration::Initialize()
{
	// Call base class
	this->irtkModelRegistration::Initialize();

	// Invert rigid transformation (to be backwards compatible)
	((irtkRigidTransformation *)_transformation)->Invert();
	((irtkRigidTransformation *)_transformation)->UpdateParameter();
}

void irtkModelRigidRegistration::Finalize()
{
	// Call base class
	this->irtkModelRegistration::Finalize();

	// Invert rigid transformation (to be backwards compatible)
	((irtkRigidTransformation *)_transformation)->Invert();
	((irtkRigidTransformation *)_transformation)->UpdateParameter();
}

double irtkModelRigidRegistration::Evaluate()
{
	int i;
	double point[3], normal[3], profile[MAX_PROFILE], similarity;

	// Print debugging information
	this->Debug("irtkModelRigidRegistration::Evaluate");

	// Get intensity profile
	_model->GetPointData()->SetActiveScalars("IntensityProfile");
	vtkDataArray *profiles = _model->GetPointData()->GetScalars();

	// Get normals
	vtkDataArray *normals = _model->GetPointData()->GetNormals();

	// Initialize similarity
	similarity = 0;
	for (i = 0; i < _model->GetNumberOfPoints(); i++) {
		_model->GetPoints()->GetPoint (i, point);
		_transformation->Transform(point[0], point[1], point[2]);
    // Get normal
    if (normals != NULL) {
      normals->GetTuple(i, normal);
      if (profiles != NULL) {
      	profiles->GetTuple(i, profile);
        // Compute metric
        similarity += _metric->Evaluate(point, normal, profile);
      } else {
        // Compute metric
        similarity += _metric->Evaluate(point, normal);
      }
    } else {
      similarity += _metric->Evaluate(point);
    }
	}

	return similarity / _model->GetNumberOfPoints();
}

#endif
