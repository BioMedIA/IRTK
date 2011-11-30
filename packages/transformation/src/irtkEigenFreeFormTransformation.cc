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

irtkEigenFreeFormTransformation::irtkEigenFreeFormTransformation() : irtkBSplineFreeFormTransformation()
{
  int i, n;

  // Initialize element labels
  n = this->NumberOfElements();
  _label = new int [n];
  for (i = 0; i < n; i++) {
    _label[i] = 0;
  }

  // Initialize fields
  _EigenValues  = new irtkVector;
  _EigenVectors = new irtkMatrix;
  _ShapeVector  = new irtkVector;
}

irtkEigenFreeFormTransformation::irtkEigenFreeFormTransformation(const irtkBSplineFreeFormTransformation &ffd, const irtkVector &EigenValues, const irtkMatrix &EigenVectors, const irtkVector &ShapeVector) : irtkBSplineFreeFormTransformation(ffd)
{
  int i, n;

  // Initialize element labels
  n = this->NumberOfElements();
  _label = new int [n];
  for (i = 0; i < n; i++) {
    _label[i] = 0;
  }

  // Initialize fields
  _EigenValues  = new irtkVector(EigenValues);
  _EigenVectors = new irtkMatrix(EigenVectors);
  _ShapeVector  = new irtkVector(ShapeVector);
}

irtkEigenFreeFormTransformation::irtkEigenFreeFormTransformation(const irtkEigenFreeFormTransformation &ffd) : irtkBSplineFreeFormTransformation(ffd)
{
  int i, n;

  // Copy labels
  n = ffd.NumberOfElements();
  _label = new int [n];
  for (i = 0; i < n; i++) {
    _label[i] = ffd._label[i];
  }

  // Copy fields
  _EigenValues  = ffd._EigenValues;
  _EigenVectors = ffd._EigenVectors;
  _ShapeVector  = ffd._ShapeVector;
}

irtkEigenFreeFormTransformation::~irtkEigenFreeFormTransformation()
{
  // Be good
  if (_label        != NULL) delete [] _label;
  if (_EigenValues  != NULL) delete _EigenValues;
  if (_EigenVectors != NULL) delete _EigenVectors;
  if (_ShapeVector  != NULL) delete _ShapeVector;
  _label        = NULL;
  _EigenValues  = NULL;
  _EigenVectors = NULL;
  _ShapeVector  = NULL;
}

void irtkEigenFreeFormTransformation::Put(int index, double shape)
{
  int i, j, k, n;
  double x, y, z, sqrt_lambda, shape_diff;

  if (index >= this->NumberOfDOFs()) {
    cerr << "irtkEigenFreeFormTransformation::Put: No such dof" << endl;
    exit(1);
  }

  // Only update control points for nonrigid modes
  if (this->GetStatus(index) == _Active) {

    // Compute lambda
    sqrt_lambda = sqrt(_EigenValues->Get(index));

    // Compute shape difference weight
    shape_diff = (shape - _ShapeVector->Get(index)) * sqrt_lambda;

    // Update control points
    n = 0;
    for (k = 0; k < _z; k++) {
      for (j = 0; j < _y; j++) {
        for (i = 0; i < _x; i++) {
          this->irtkBSplineFreeFormTransformation::Get(i, j, k, x, y, z);
          // Add difference between new and old shape parameter
          this->irtkBSplineFreeFormTransformation::Put(i, j, k,  x + shape_diff * _EigenVectors->Get(3*n, index), y + shape_diff * _EigenVectors->Get(3*n+1, index), z + shape_diff * _EigenVectors->Get(3*n+2, index));
          n++;
        }
      }
    }

    // Change shape vector
    _ShapeVector->Put(index, shape);
  }
}

//
// Overloaded I/O methods
//

istream& irtkEigenFreeFormTransformation::Import(istream& is)
{
  char buffer[255];
  int i, n;

  // Read keyword
  is >> buffer;

  // Call parent I/O
  ((irtkBSplineFreeFormTransformation *)this)->Import(is);
  //
  // Read element labels
  //

  // Read label keyword
  is >> buffer;
  if (strcmp(buffer, "LABEL:") != 0) {
    cerr << "istream& operator>> irtkEigenFreeFormTransformation: ";
    cerr << "Expecting LABEL: attribute " << endl;
    exit(1);
  }

  // Skip rest of line
  is.get(buffer, 255);
  is.clear();
  is.seekg(1, ios::cur);

  // Get number of elements
  n = this->NumberOfElements();

  // Allocate temporary memory
  unsigned char *label = new unsigned char[n];

  // Read binary data for element label
  is.read((char *)&(label[0]), n*sizeof(unsigned char));

  // Delete old labels
  if (_label != NULL) delete []_label;

  // Allocate new labels
  _label = new int[n];

  // Convert label
  for (i = 0; i < n; i++) {
    _label[i] = (int)label[i];
  }

  // Free temporary memory
  delete [] label;

  // Skip rest of line
  is.get(buffer, 255);
  is.clear();
  is.seekg(1, ios::cur);

  //
  // Read in EigenSystem
  //

  is >> *(_EigenValues);
  is >> *(_EigenVectors);
  is >> *(_ShapeVector);

  return is;
}

ostream& irtkEigenFreeFormTransformation::Export(ostream& os)
{
  int i, n;

  // Write keyword
  os << "EFFD:" << endl;

  // Call parent I/O
  ((irtkBSplineFreeFormTransformation *)this)->Export(os);

  //
  // Write out element labels
  //

  // Get number of elements
  n = NumberOfElements();

  // Allocate
  unsigned char *label = new unsigned char [n*sizeof(unsigned char)];

  // Convert
  for (i = 0; i < n; i++) {
    label[i] = (unsigned char)_label[i];
  }

  // Write binary data for element label
  os << "LABEL:" << endl;
  os.write((char *)&(label[0]), n*sizeof(unsigned char));

  // Write end of line
  os << endl;

  // Free temporary memory
  delete [] label;

  //
  // Write out EigenSystem
  //

  os << *(_EigenValues)  << endl;
  os << *(_EigenVectors) << endl;
  os << *(_ShapeVector)  << endl;

  return os;
}

int irtkEigenFreeFormTransformation::CheckHeader(char *name)
{
  char buffer[255];

  ifstream from(name);

  if (!from) {
    cerr << "irtkEigenFreeFormTransformation"
         << "::CheckHeader: Can't open file " << name << "\n";
    exit(1);
  }

  // Read and check keyword
  from >> buffer;
  if (strcmp(buffer, "EFFD:") == 0) {
    return true;
  } else {
    return false;
  }
}

void irtkEigenFreeFormTransformation::Print()
{
  int i, n;

  // Print keyword
  cout << "EFFD:" << endl;

  // Print FFD
  this->irtkBSplineFreeFormTransformation::Print();

  // Print labels
  n = this->NumberOfElements();
  cout << "LABEL:\n";
  for (i = 0; i < n; i++) {
    cout << _label[i] << " ";
  }
  cout << endl;

  // Print EigenSystem
  _EigenValues->Print();
  _EigenVectors->Print();
  _ShapeVector->Print();
}

irtkCifstream& irtkEigenFreeFormTransformation::Import(irtkCifstream& from)
{
  double *data;
  int i, j, k, index;
  unsigned int magic_no, trans_type;

  // Read magic no. for transformations
  from.ReadAsUInt(&magic_no, 1);
  if (magic_no != IRTKTRANSFORMATION_MAGIC) {
    cerr << "irtkEigenFreeFormTransformation::Import: Not a vaild transformation file" << endl;
    exit(1);
  }

  // Read transformation type
  from.ReadAsUInt(&trans_type, 1);
  if (trans_type != IRTKTRANSFORMATION_EIGEN_FFD) {
    cerr << "irtkEigenFreeFormTransformation::Import: Not a vaild Eigen FFD transformation" << endl;
    exit(1);
  }

  // Free memory if necessary
  _data = this->Deallocate(_data, _x, _y, _z);
  delete []_status;

  // Read no of control points
  from.ReadAsInt(&_x, 1);
  from.ReadAsInt(&_y, 1);
  from.ReadAsInt(&_z, 1);

  // Read orientation of bounding box
  from.ReadAsDouble(_xaxis, 3);
  from.ReadAsDouble(_yaxis, 3);

  // Read spacing of bounding box
  from.ReadAsDouble(&_dx, 1);
  from.ReadAsDouble(&_dy, 1);
  from.ReadAsDouble(&_dz, 1);

  // Read spacing of bounding box
  from.ReadAsDouble(&_origin._x, 1);
  from.ReadAsDouble(&_origin._y, 1);
  from.ReadAsDouble(&_origin._z, 1);

  // Initialize control points
  _data = this->Allocate(_data, _x, _y, _z);

  // Initialize control points
  for (i = -2; i < _x+2; i++) {
    for (j = -2; j < _y+2; j++) {
      for (k = -2; k < _z+2; k++) {
        _data[k][j][i]._x = 0;
        _data[k][j][i]._y = 0;
        _data[k][j][i]._z = 0;
      }
    }
  }

  // Allocate temporary memory
  data = new double[3*_x*_y*_z];

  // Read control point data
  from.ReadAsDouble(data, 3*_x*_y*_z);

  // Convert data
  index = 0;
  for (i = 0; i < _x; i++) {
    for (j = 0; j < _y; j++) {
      for (k = 0; k < _z; k++) {
        _data[k][j][i]._x = data[index];
        _data[k][j][i]._y = data[index+1];
        _data[k][j][i]._z = data[index+2];
        index += 3;
      }
    }
  }

  // Free temporary memory
  delete []data;

  // Initialize memory for control point status
  _status = new _Status[3*_x*_y*_z];

  // Read control point status
  from.ReadAsInt((int *)_status, 3*_x*_y*_z);

  // Update transformation matrix
  this->UpdateMatrix();

  return from;
}

irtkCofstream& irtkEigenFreeFormTransformation::Export(irtkCofstream& to)
{
  double *data;
  int i, j, k, index;
  unsigned int magic_no, trans_type;

  // Write magic no. for transformations
  magic_no = IRTKTRANSFORMATION_MAGIC;
  to.WriteAsUInt(&magic_no, 1);

  // Write transformation type
  trans_type = IRTKTRANSFORMATION_EIGEN_FFD;
  to.WriteAsUInt(&trans_type, 1);

  // Write no of control points
  to.WriteAsInt(&_x, 1);
  to.WriteAsInt(&_y, 1);
  to.WriteAsInt(&_z, 1);

  // Write orientation of bounding box
  to.WriteAsDouble(_xaxis, 3);
  to.WriteAsDouble(_yaxis, 3);

  // Write spacing of bounding box
  to.WriteAsDouble(&_dx, 1);
  to.WriteAsDouble(&_dy, 1);
  to.WriteAsDouble(&_dz, 1);

  // Write spacing of bounding box
  to.WriteAsDouble(&_origin._x, 1);
  to.WriteAsDouble(&_origin._y, 1);
  to.WriteAsDouble(&_origin._z, 1);

  // Allocate temporary memory
  data = new double[3*_x*_y*_z];

  // Convert data
  index = 0;
  for (i = 0; i < _x; i++) {
    for (j = 0; j < _y; j++) {
      for (k = 0; k < _z; k++) {
        _data[k][j][i]._x = data[index];
        _data[k][j][i]._y = data[index+1];
        _data[k][j][i]._z = data[index+2];
        index += 3;
      }
    }
  }

  // Write control point data
  to.WriteAsDouble(data, 3*_x*_y*_z);

  // Free temporary memory
  delete []data;

  // Write control point status
  to.WriteAsInt((int *)_status, 3*_x*_y*_z);

  return to;
}
