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

#include <vtkPolyDataNormals.h>

// The original image
extern irtkGreyImage *tmp_image;

irtkModelFreeFormRegistration::irtkModelFreeFormRegistration() {
	// Print debugging information
	this->Debug("irtkModelFreeFormRegistration::irtkModelFreeFormRegistration");

	// Default optimization
	_OptimizationMethod = GradientDescent;

	// Default parameters for non-rigid registration
	_Lambda1 = 0;
	_Lambda2 = 0;
	_Lambda3 = 0;
	_DX = 20;
	_DY = 20;
	_DZ = 20;
	_Subdivision = True;
}

void irtkModelFreeFormRegistration::GuessParameter() {
	int i;
	double xsize, ysize, zsize, spacing;

	if ((_image == NULL) || (_model == NULL)) {
		cerr
				<< "irtkModelFreeFormRegistration::GuessParameter: Image and/or model not found"
				<< endl;
		exit(1);
	}

	// Default parameters for registration
	_NumberOfLevels = 1;

	// Default parameters for optimization
	_OptimizationMethod = GradientDescent;
	_Epsilon = 0.0001;

	// Read target pixel size
	_image->GetPixelSize(&xsize, &ysize, &zsize);

	// Use xsize as spacing
	spacing = xsize;

	// Default target parameters
	_ImageBlurring[0] = GuessResolution(xsize, ysize, zsize) / 2.0;
	_ImageResolution[0][0] = GuessResolution(xsize, ysize, zsize);
	_ImageResolution[0][1] = GuessResolution(xsize, ysize, zsize);
	_ImageResolution[0][2] = GuessResolution(xsize, ysize, zsize);

	for (i = 1; i < _NumberOfLevels; i++) {
		_ImageBlurring[i] = _ImageBlurring[i - 1] * 2;
		_ImageResolution[i][0] = _ImageResolution[i - 1][0] * 2;
		_ImageResolution[i][1] = _ImageResolution[i - 1][1] * 2;
		_ImageResolution[i][2] = _ImageResolution[i - 1][2] * 2;
	}

	// Default parameters for non-rigid registration
	_Lambda1 = 0;
	_Lambda2 = 0;
	_Lambda3 = 0;
	_DX = _image->GetX() * spacing / 10.0;
	_DY = _image->GetX() * spacing / 10.0;
	_DZ = _image->GetX() * spacing / 10.0;
	_Subdivision = True;

	// Remaining parameters
	for (i = 0; i < _NumberOfLevels; i++) {
		_NumberOfIterations[i] = 10;
		_NumberOfSteps[i] = 4;
		_LengthOfSteps[i] = _DX / 8.0 * pow(2.0, i);
	}
}

void irtkModelFreeFormRegistration::Initialize() {
	// Print debugging information
	this->Debug("irtkModelFreeFormRegistration::Initialize");

	// Initialize base class
	this->irtkModelRegistration::Initialize();

	// Pointer to multi-level FFD
	_mffd = (irtkMultiLevelFreeFormTransformation *) _transformation;

	// Create FFD
	if (_mffd->NumberOfLevels() == 0) {
		_affd = new irtkBSplineFreeFormTransformation(*_image, this->_DX,
				this->_DY, this->_DZ);
	} else {
		_affd
				= (irtkBSplineFreeFormTransformation *) _mffd->PopLocalTransformation();
	}
}

void irtkModelFreeFormRegistration::Initialize(int level) {
	int a, b, c, i, j, l, m, n;
	double point[3];

	// Print debugging information
	this->Debug("irtkModelFreeFormRegistration::Initialize(int)");

	// Initialize base class
	this->irtkModelRegistration::Initialize(level);

	// Tell optimizer which transformation to optimize
	_optimizer->SetTransformation(_affd);

	// Create list of indices for each control point
	_id = new list<int> [_affd->GetX() * _affd->GetY() * _affd->GetZ()];

	for (i = 0; i < _model->GetNumberOfPoints(); i++) {
		_model->GetPoints()->GetPoint(i, point);
		_affd->WorldToLattice(point[0], point[1], point[2]);
		l = int(floor(point[0]));
		m = int(floor(point[1]));
		n = int(floor(point[2]));
		for (a = -1; a < 3; a++) {
			for (b = -1; b < 3; b++) {
				for (c = -1; c < 3; c++) {
					if ((l + a >= 0) && (l + a < _affd->GetX()) && (m + b >= 0) && (m + b
							< _affd->GetY()) && (n + c >= 0) && (n + c < _affd->GetZ())) {
						j = _affd->LatticeToIndex(l + a, m + b, n + c);
						_id[j].push_back(i);
					}
				}
			}
		}
	}

	// Create list of similarities
	_similarity = new double[_model->GetNumberOfPoints()];
}

void irtkModelFreeFormRegistration::Finalize()
{
	// Print debugging information
	this->Debug("irtkModelFreeFormRegistration::Finalize");

	// Push local transformation back on transformation stack
	_mffd->PushLocalTransformation(_affd);

	// Finalize base class
	this->irtkModelRegistration::Finalize();
}

void irtkModelFreeFormRegistration::Finalize(int level)
{
	// Print debugging information
	this->Debug("irtkModelFreeFormRegistration::Finalize(int)");

	// Finalize base class
	this->irtkModelRegistration::Finalize(level);

	// Check if we are not at the lowest level of resolution
	if (level != 0) {
		if (this->_Subdivision == True) {
			_affd->Subdivide();
		} else {
			// Push local transformation back on transformation stack
			_mffd->PushLocalTransformation(_affd);

			// Create new FFD
			_affd = new irtkBSplineFreeFormTransformation(*_image, this->_DX / pow(
					2.0, this->_NumberOfLevels - level), this->_DY / pow(2.0,
					this->_NumberOfLevels - level), this->_DZ / pow(2.0,
					this->_NumberOfLevels - level));
		}
	}

	// Delete similarities
	delete _similarity;

	// Delete list of ids
	delete []_id;
}

void irtkModelFreeFormRegistration::UpdateLUT() {
	/*
	 int i;
	 double point[3], normal[3]

	 for (i = 0; i < _model->GetNumberOfPoints(); i++) {
	 _model->GetPoints()->GetPoint (i, point);
	 _transform->Transform(m, point[0], point[1], point[2]);
	 }
	 */

	vtkPolyDataNormals *filter = vtkPolyDataNormals::New();
	filter->SetInput(_model);
	filter->SplittingOff();
	filter->Update();
	_model->GetPointData()->SetNormals(
			filter->GetOutput()->GetPointData()->GetNormals());
	//  filter->Delete();
}

double irtkModelFreeFormRegistration::SmoothnessPenalty() {
	int i, j, k;
	double x, y, z, penalty;

	penalty = 0;
	for (k = 0; k < _affd->GetZ(); k++) {
		for (j = 0; j < _affd->GetY(); j++) {
			for (i = 0; i < _affd->GetX(); i++) {
				x = i;
				y = j;
				z = k;
				_affd->LatticeToWorld(x, y, z);
				penalty += _affd->Bending(x, y, z);
			}
		}
	}
	return -penalty / _affd->NumberOfDOFs();
}

double irtkModelFreeFormRegistration::SmoothnessPenalty(int index) {
	int i, j, k;
	double x, y, z;

	_affd->IndexToLattice(index, i, j, k);
	x = i;
	y = j;
	z = k;
	_affd->LatticeToWorld(x, y, z);
	return -_affd->Bending(x, y, z);
}

double irtkModelFreeFormRegistration::VolumePreservationPenalty() {
	int i, j, k;
	double x, y, z, penalty;

	penalty = 0;
	for (k = 0; k < _affd->GetZ(); k++) {
		for (j = 0; j < _affd->GetY(); j++) {
			for (i = 0; i < _affd->GetZ(); i++) {
				x = i;
				y = j;
				z = k;
				_affd->LatticeToWorld(x, y, z);
				// Torsten Rohlfing et al. MICCAI'01 (w/o scaling correction):
				penalty += fabs(log(_affd->irtkTransformation::Jacobian(x, y, z)));
			}
		}
	}

	// Normalize sum by number of DOFs
	return penalty / (double) _affd->NumberOfDOFs();
}

double irtkModelFreeFormRegistration::VolumePreservationPenalty(int index) {
	int i, j, k;
	double x, y, z;

	_affd->IndexToLattice(index, i, j, k);
	x = i;
	y = j;
	z = k;
	_affd->LatticeToWorld(x, y, z);
	return fabs(log(_affd->irtkTransformation::Jacobian(x, y, z)));
}

double irtkModelFreeFormRegistration::TopologyPreservationPenalty() {
	int i, j, k;
	double x, y, z, jac, penalty;

	penalty = 0;
	for (k = 0; k < _affd->GetZ() - 1; k++) {
		for (j = 0; j < _affd->GetY() - 1; j++) {
			for (i = 0; i < _affd->GetZ() - 1; i++) {
				x = i + 0.5;
				y = j + 0.5;
				z = k + 0.5;
				_affd->LatticeToWorld(x, y, z);
				jac = _affd->irtkTransformation::Jacobian(x, y, z);
				if (jac < 0.3) {
					penalty += 10 * jac * jac + 0.1 / (jac * jac) - 2.0;
				}
			}
		}
	}
	return -penalty;
}

double irtkModelFreeFormRegistration::TopologyPreservationPenalty(int index) {
	int i, j, k, l, m, n;
	double x, y, z, jac, penalty;

	penalty = 0;
	for (l = 0; l <= 1; l++) {
		for (m = 0; m <= 1; m++) {
			for (n = 0; n <= 1; n++) {
				_affd->IndexToLattice(index, i, j, k);
				x = i + l - 0.5;
				y = j + m - 0.5;
				z = k + n - 0.5;
				_affd->LatticeToWorld(x, y, z);
				jac = _affd->irtkTransformation::Jacobian(x, y, z);
				if (jac < 0.3) {
					penalty += 10 * jac * jac + 0.1 / (jac * jac) - 2.0;
				}
			}
		}
	}
	return -penalty;
}

double irtkModelFreeFormRegistration::Evaluate() {
	int i;
	double point[3], normal[3], profile[MAX_PROFILE], tmp, similarity;

	// Print debugging information
	this->Debug("irtkModelFreeFormRegistration::Evaluate");

	// Update normals
	this->UpdateLUT();

	// Get intensity profile
	_model->GetPointData()->SetActiveScalars("IntensityProfile");
	vtkDataArray *profiles = _model->GetPointData()->GetScalars();

	// Get normals
	vtkDataArray *normals = _model->GetPointData()->GetNormals();

	// Initialize similarity
	similarity = 0;
	for (i = 0; i < _model->GetNumberOfPoints(); i++) {
		// Get transformed point
		_model->GetPoints()->GetPoint(i, point);
		_affd->Transform(point[0], point[1], point[2]);
		// Get normal
		if (normals != NULL) {
			normals->GetTuple(i, normal);
			if (profiles != NULL) {
				profiles->GetTuple(i, profile);
				// Compute metric
				tmp = _metric->Evaluate(point, normal, profile);
			} else {
				// Compute metric
				tmp = _metric->Evaluate(point, normal);
			}
		} else {
			tmp = _metric->Evaluate(point);
		}
		// Add metric
		similarity += tmp;
		// Store metric
		_similarity[i] = tmp;
	}

	// Add penalty for smoothness
	if (this->_Lambda1 > 0) {
		similarity += this->_Lambda1 * this->SmoothnessPenalty();
	}
	// Add penalty for volume preservation
	if (this->_Lambda2 > 0) {
		similarity += this->_Lambda2 * this->VolumePreservationPenalty();
	}
	// Add penalty for topology preservation
	if (this->_Lambda3 > 0) {
		similarity += this->_Lambda3 * this->TopologyPreservationPenalty();
	}
	// Return similarity measure + penalty terms
	return similarity;
}

double irtkModelFreeFormRegistration::EvaluateDerivative(int index, double step) {
	int j;
	list<int>::iterator i;
	double point[3], normal[3], profile[MAX_PROFILE], similarityA, similarityB;

	// Print debugging information
	this->Debug("irtkModelFreeFormRegistration::EvaluateDerivative(int, double)");

	// Get intensity profile
	_model->GetPointData()->SetActiveScalars("IntensityProfile");
	vtkDataArray *profiles = _model->GetPointData()->GetScalars();

	// Get normals
	vtkDataArray *normals = _model->GetPointData()->GetNormals();

	if (index < _affd->GetX() * _affd->GetY() * _affd->GetZ()) {
		j = index;
	} else {
		if (index < 2 * _affd->GetX() * _affd->GetY() * _affd->GetZ()) {
			j = index - _affd->GetX() * _affd->GetY() * _affd->GetZ();
		} else {
			j = index - 2 * _affd->GetX() * _affd->GetY() * _affd->GetZ();
		}
	}

	// Save value of DOF for which we calculate the derivative
	double dof = _affd->Get(index);

	// Add penalties
	_affd->Put(index, dof + step);

	// Initialize similarity
	similarityA = 0;
	for (i = _id[j].begin(); i != _id[j].end(); ++i) {
		_model->GetPoints()->GetPoint(*i, point);
		_affd->Transform(point[0], point[1], point[2]);
		// Get normal
		if (normals != NULL) {
			normals->GetTuple(*i, normal);
			if (profiles != NULL) {
				profiles->GetTuple(*i, profile);
				// Compute metric
				similarityA -= _similarity[*i];
				similarityA += _metric->Evaluate(point, normal, profile);
			} else {
				// Compute metric
				similarityA -= _similarity[*i];
				similarityA += _metric->Evaluate(point, normal);
			}
		} else {
			similarityA -= _similarity[*i];
			similarityA += _metric->Evaluate(point);
		}
	}

	// Smoothness
	if (this->_Lambda1 > 0) {
		similarityA += this->_Lambda1 * this->SmoothnessPenalty(index);
	}
	// Volume preservation
	if (this->_Lambda2 > 0) {
		similarityA += this->_Lambda2 * this->VolumePreservationPenalty(index);
	}
	// Topology preservation
	if (this->_Lambda3 > 0) {
		similarityA += this->_Lambda3 * this->TopologyPreservationPenalty(index);
	}

	// Add penalties
	_affd->Put(index, dof - step);

	// Initialize similarity
	similarityB = 0;
	for (i = _id[j].begin(); i != _id[j].end(); ++i) {
		_model->GetPoints()->GetPoint(*i, point);
		_affd->Transform(point[0], point[1], point[2]);
		// Get normal
		if (normals != NULL) {
			normals->GetTuple(*i, normal);
			if (profiles != NULL) {
				profiles->GetTuple(*i, profile);
				// Compute metric
				similarityB -= _similarity[*i];
				similarityB += _metric->Evaluate(point, normal, profile);
			} else {
				// Compute metric
				similarityB -= _similarity[*i];
				similarityB += _metric->Evaluate(point, normal);
			}
		} else {
			similarityB -= _similarity[*i];
			similarityB += _metric->Evaluate(point);
		}
	}

	// Smoothness
	if (this->_Lambda1 > 0) {
		similarityB += this->_Lambda1 * this->SmoothnessPenalty(index);
	}
	// Volume preservation
	if (this->_Lambda2 > 0) {
		similarityB += this->_Lambda2 * this->VolumePreservationPenalty(index);
	}
	// Topology preservation
	if (this->_Lambda3 > 0) {
		similarityB += this->_Lambda3 * this->TopologyPreservationPenalty(index);
	}

	// Restore value of DOF for which we calculate the derivative
	_affd->Put(index, dof);

	return similarityA - similarityB;
}

double irtkModelFreeFormRegistration::EvaluateGradient(float step, float *dx) {
	int i;
	double norm;

	// Update lookup table
	this->UpdateLUT();

	for (i = 0; i < _affd->NumberOfDOFs(); i++) {
		if (_affd->irtkTransformation::GetStatus(i) == _Active) {
			dx[i] = this->EvaluateDerivative(i, step);
		} else {
			dx[i] = 0;
		}
	}

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

Bool irtkModelFreeFormRegistration::Read(char *buffer1, char *buffer2,
		int &level) {
	int ok = False;

	if ((strstr(buffer1, "Lambda ") != NULL) || (strstr(buffer1, "Lambda1")
			!= NULL)) {
		this->_Lambda1 = atof(buffer2);
		cout << "Lambda 1 is ... " << this->_Lambda1 << endl;
		ok = True;
	}
	if (strstr(buffer1, "Lambda2") != NULL) {
		this->_Lambda2 = atof(buffer2);
		cout << "Lambda 2 is ... " << this->_Lambda2 << endl;
		ok = True;
	}
	if (strstr(buffer1, "Lambda3") != NULL) {
		this->_Lambda3 = atof(buffer2);
		cout << "Lambda 3 is ... " << this->_Lambda3 << endl;
		ok = True;
	}
	if (strstr(buffer1, "Control point spacing in X") != NULL) {
		this->_DX = atof(buffer2);
		cout << "Control point spacing in X is ... " << this->_DX << endl;
		ok = True;
	}
	if (strstr(buffer1, "Control point spacing in Y") != NULL) {
		this->_DY = atof(buffer2);
		cout << "Control point spacing in Y is ... " << this->_DY << endl;
		ok = True;
	}
	if (strstr(buffer1, "Control point spacing in Z") != NULL) {
		this->_DZ = atof(buffer2);
		cout << "Control point spacing in Z is ... " << this->_DZ << endl;
		ok = True;
	}
	if (strstr(buffer1, "Subdivision") != NULL) {
		if ((strcmp(buffer2, "False") == 0) || (strcmp(buffer2, "No") == 0)) {
			this->_Subdivision = False;
			cout << "Subdivision is ... false" << endl;
		} else {
			if ((strcmp(buffer2, "True") == 0) || (strcmp(buffer2, "Yes") == 0)) {
				this->_Subdivision = True;
				cout << "Subdivision is ... true" << endl;
			} else {
				cerr << "Can't read boolean value = " << buffer2 << endl;
				exit(1);
			}
		}
		ok = True;
	}

	if (ok == False) {
		return this->irtkModelRegistration::Read(buffer1, buffer2, level);
	} else {
		return ok;
	}
}

void irtkModelFreeFormRegistration::Write(ostream &to) {
	to << "\n#\n# Non-rigid registration parameters\n#\n\n";
	to << "Lambda1                           = " << this->_Lambda1 << endl;
	to << "Lambda2                           = " << this->_Lambda2 << endl;
	to << "Lambda3                           = " << this->_Lambda3 << endl;
	to << "Control point spacing in X        = " << this->_DX << endl;
	to << "Control point spacing in Y        = " << this->_DY << endl;
	to << "Control point spacing in Z        = " << this->_DZ << endl;
	if (_Subdivision == True) {
		to << "Subdivision                       = True" << endl;
	} else {
		to << "Subdivision                       = False" << endl;
	}

	this->irtkModelRegistration::Write(to);
}

#endif
