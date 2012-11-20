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

#include <newt2.h>

irtkMultiLevelFreeFormTransformation::irtkMultiLevelFreeFormTransformation()
{
  int i;

  // Initialize local transformation
  for (i = 0; i < MAX_TRANS; i++) {
    _localTransformation[i] = NULL;
  }

  _NumberOfLevels = 0;
}

irtkMultiLevelFreeFormTransformation::irtkMultiLevelFreeFormTransformation(const irtkMultiLevelFreeFormTransformation &transformation) : irtkAffineTransformation(transformation)
{
  int i;

  // Initialize local transformation
  for (i = 0; i < transformation._NumberOfLevels; i++) {
    if (strcmp(transformation._localTransformation[i]->NameOfClass(),
               "irtkBSplineFreeFormTransformation3D") == 0) {
      irtkBSplineFreeFormTransformation *ffd = dynamic_cast<irtkBSplineFreeFormTransformation *>(transformation._localTransformation[i]);
      _localTransformation[i] = new irtkBSplineFreeFormTransformation(*ffd);
    }
    if (strcmp(transformation._localTransformation[i]->NameOfClass(),
               "irtkLinearFreeFormTransformation") == 0) {
      irtkLinearFreeFormTransformation *ffd = dynamic_cast<irtkLinearFreeFormTransformation *>(transformation._localTransformation[i]);
      _localTransformation[i] = new irtkLinearFreeFormTransformation(*ffd);
    }
    if (strcmp(transformation._localTransformation[i]->NameOfClass(),
               "irtkBSplineFreeFormTransformation4D") == 0) {
      irtkBSplineFreeFormTransformation4D *ffd = dynamic_cast<irtkBSplineFreeFormTransformation4D *>(transformation._localTransformation[i]);
      _localTransformation[i] = new irtkBSplineFreeFormTransformation4D(*ffd);
    }
    if (strcmp(transformation._localTransformation[i]->NameOfClass(),
			   "irtkBSplineFreeFormTransformationPeriodic") == 0) {
	  irtkBSplineFreeFormTransformationPeriodic *ffd = dynamic_cast<irtkBSplineFreeFormTransformationPeriodic *>(transformation._localTransformation[i]);
	  _localTransformation[i] = new irtkBSplineFreeFormTransformationPeriodic(*ffd);
	}
  }

  _NumberOfLevels = transformation._NumberOfLevels;
}

irtkMultiLevelFreeFormTransformation::irtkMultiLevelFreeFormTransformation(const irtkRigidTransformation &transformation) : irtkAffineTransformation(transformation)
{
  int i;

  // Initialize local transformation
  for (i = 0; i < MAX_TRANS; i++) {
    _localTransformation[i] = NULL;
  }

  _NumberOfLevels = 0;
}

irtkMultiLevelFreeFormTransformation::irtkMultiLevelFreeFormTransformation(const irtkAffineTransformation &transformation) : irtkAffineTransformation(transformation)
{
  int i;

  // Initialize local transformation
  for (i = 0; i < MAX_TRANS; i++) {
    _localTransformation[i] = NULL;
  }

  _NumberOfLevels = 0;
}

irtkMultiLevelFreeFormTransformation::~irtkMultiLevelFreeFormTransformation()
{
  int i;

  // Delete local transformations
  for (i = 0; i < _NumberOfLevels; i++) {
    delete _localTransformation[i];
  }

  _NumberOfLevels = 0;
}

void irtkMultiLevelFreeFormTransformation::Transform(double &x, double &y, double &z, double t)
{
  int i;
  double u, v, w, dx, dy, dz;

  // Initialize displacement
  dx = 0;
  dy = 0;
  dz = 0;

  // Compute local transformation
  for (i = 0; i < _NumberOfLevels; i++) {

    // Reset
    u = x;
    v = y;
    w = z;

    // Displacement
    _localTransformation[i]->LocalDisplacement(u, v, w, t);

    // Add to displacement
    dx += u;
    dy += v;
    dz += w;
  }

  // Compute global transformation
  this->irtkAffineTransformation::Transform(x, y, z, t);

  // Compute sum
  x += dx;
  y += dy;
  z += dz;
}

void irtkMultiLevelFreeFormTransformation::Displacement(double& x, double& y, double& z, double t)
{
  // Global displacement.
  double globalX, globalY, globalZ;

  globalX = x;
  globalY = y;
  globalZ = z;

  this->GlobalDisplacement(globalX, globalY, globalZ, t);

  // Local displacement.
  double localX, localY, localZ;

  localX = x;
  localY = y;
  localZ = z;

  this->LocalDisplacement(localX, localY, localZ, t);

  x = globalX + localX;
  y = globalY + localY;
  z = globalZ + localZ;
}

void irtkMultiLevelFreeFormTransformation::Transform(int n, double &x, double &y, double &z, double t)
{
  int i;
  double u, v, w, dx, dy, dz;

  if (n > _NumberOfLevels) {
    cerr << "irtkMultiLevelFreeFormTransformation::Transform: No such "
    << "transformation" << endl;
    exit(1);
  }

  // Initialize displacement
  dx = 0;
  dy = 0;
  dz = 0;

  // Compute local transformation
  for (i = 0; i < n; i++) {

    // Reset
    u = x;
    v = y;
    w = z;

    // Displacement
    _localTransformation[i]->LocalDisplacement(u, v, w, t);

    // Add to displacement
    dx += u;
    dy += v;
    dz += w;
  }

  // Compute global transformation
  this->irtkAffineTransformation::Transform(x, y, z, t);

  // Compute sum
  x += dx;
  y += dy;
  z += dz;
}

void irtkMultiLevelFreeFormTransformation::GlobalTransform(double &x, double &y, double &z, double t)
{
  // Compute global transformation
  this->irtkAffineTransformation::Transform(x, y, z, t);
}

void irtkMultiLevelFreeFormTransformation::GlobalDisplacement(double &x, double &y, double &z, double t)
{
  double u, v, w;

  u = x;
  v = y;
  w = z;

  // Compute global transformation
  this->irtkAffineTransformation::Transform(u, v, w, t);

  x = u - x;
  y = v - y;
  z = w - z;
}

void irtkMultiLevelFreeFormTransformation::LocalTransform(double &x, double &y, double &z, double t)
{
  int i;
  double u, v, w, dx, dy, dz;

  // Initialize displacement
  dx = 0;
  dy = 0;
  dz = 0;

  // Compute local transformation
  for (i = 0; i < _NumberOfLevels; i++) {

    // Reset
    u = x;
    v = y;
    w = z;

    // Displacement
    _localTransformation[i]->LocalDisplacement(u, v, w, t);

    // Add to displacement
    dx += u;
    dy += v;
    dz += w;
  }

  // Compute sum
  x += dx;
  y += dy;
  z += dz;
}

void irtkMultiLevelFreeFormTransformation::LocalTransform(int n, double &x, double &y, double &z, double t)
{
  int i;
  double u, v, w, dx, dy, dz;

  if (n > _NumberOfLevels) {
    cerr << "irtkMultiLevelFreeFormTransformation::LocalTransform: No such "
    << "transformation" << endl;
    exit(1);
  }

  // Initialize displacement
  dx = 0;
  dy = 0;
  dz = 0;

  // Compute local transformation
  for (i = 0; i < n; i++) {

    // Reset
    u = x;
    v = y;
    w = z;

    // Displacement
    _localTransformation[i]->LocalDisplacement(u, v, w, t);

    // Add to displacement
    dx += u;
    dy += v;
    dz += w;
  }

  // Compute sum
  x += dx;
  y += dy;
  z += dz;
}

void irtkMultiLevelFreeFormTransformation::LocalDisplacement(int n, double &x, double &y, double &z, double t)
{
  int i;
  double u, v, w, dx, dy, dz;

  if (n > _NumberOfLevels) {
    cerr << "irtkMultiLevelFreeFormTransformation::LocalDisplacement: No such "
    << "transformation" << endl;
    exit(1);
  }

  // Initialize displacement
  dx = 0;
  dy = 0;
  dz = 0;

  // Compute local transformation
  for (i = 0; i < n; i++) {

    // Reset
    u = x;
    v = y;
    w = z;

    // Displacement
    _localTransformation[i]->LocalDisplacement(u, v, w, t);

    // Add to displacement
    dx += u;
    dy += v;
    dz += w;
  }

  // Compute sum
  x = dx;
  y = dy;
  z = dz;
}

void irtkMultiLevelFreeFormTransformation::LocalDisplacement(double &x, double &y, double &z, double t)
{
  int i;
  double u, v, w, dx, dy, dz;

  // Initialize displacement
  dx = 0;
  dy = 0;
  dz = 0;

  // Compute local transformation
  for (i = 0; i < _NumberOfLevels; i++) {

    // Reset
    u = x;
    v = y;
    w = z;

    // Displacement
    _localTransformation[i]->LocalDisplacement(u, v, w, t);

    // Add to displacement
    dx += u;
    dy += v;
    dz += w;
  }

  // Compute sum
  x = dx;
  y = dy;
  z = dz;
}

void irtkMultiLevelFreeFormTransformation::InterpolateGlobalDisplacement(irtkBSplineFreeFormTransformation3D *f)
{
	int i, j, k, count, noOfCPs;
	double x, y, z;

  noOfCPs = f->NumberOfDOFs() / 3;

  double *dx = new double[noOfCPs];
  double *dy = new double[noOfCPs];
  double *dz = new double[noOfCPs];

  count = 0;
  for (k = 0; k < f->GetZ(); k++){
  	for (j = 0; j < f->GetY(); j++){
  		for (i = 0; i < f->GetX(); i++){

  			// Reset.
  			f->Put(i, j, k, 0.0f, 0.0f, 0.0f);

  			x = i;
  			y = j;
  			z = k;

  			f->LatticeToWorld(x, y, z);
  		  this->GlobalDisplacement(x,y,z);

  		  dx[count] = x;
  		  dy[count] = y;
  		  dz[count] = z;
  		  count++;
  		}
  	}
  }

  // Compute B-spline coefficients that provide an interpolation
  // of the global affine component.
  f->Interpolate(dx, dy, dz);

  delete [] dx;
  delete [] dy;
  delete [] dz;

  // Some reporting.
  if (f->GetX() < 4 || f->GetY() < 4 || f->GetZ() < 4){
  	cerr << "irtkMultiLevelFreeFormTransformation::InterpolateGlobalDisplacement : ";
  	cerr << "Very small lattice for interpolation. Result likely to be inaccurate." << endl;
  	return;
  }

  double totalVol, effectiveVol;

  totalVol = f->GetX() * f->GetXSpacing() + f->GetY() * f->GetYSpacing() + f->GetZ() * f->GetZSpacing();
  effectiveVol = (f->GetX()-4) * f->GetXSpacing() + (f->GetY()-4) * f->GetYSpacing() + (f->GetZ()-4) * f->GetZSpacing();

  cout << "irtkMultiLevelFreeFormTransformation::InterpolateGlobalDisplacement : ";
  cout << "Accurate interpolation of affine transformation over ";
  printf("% .1f %% of lattice volume\n", 100.0 * effectiveVol / totalVol);

}

void irtkMultiLevelFreeFormTransformation::MergeGlobalIntoLocalDisplacement()
{
	int i;

  irtkBSplineFreeFormTransformation3D *ffd =
  		dynamic_cast<irtkBSplineFreeFormTransformation3D *> (this->GetLocalTransformation(0));

  if (ffd == NULL){
  	cerr << "irtkFluidFreeFormTransformation::MergeGlobalIntoLocalDisplacement: ";
  	cerr << "FFD should be of type irtkBSplineFreeFormTransformation3D" << endl;
  	exit(1);
  }

  // Get a copy for making a FFD interpolation of the global affine component.
  irtkBSplineFreeFormTransformation3D *ffdCopy = new irtkBSplineFreeFormTransformation3D(*ffd);

  this->InterpolateGlobalDisplacement(ffdCopy);

  // Add the calculated coefficents to the coefficients of the
  // first FFD in the given MFFD.
  for (i = 0; i < ffd->NumberOfDOFs(); i++){
  	ffd->Put(i, ffd->Get(i) + ffdCopy->Get(i));
  }

  // Reset matrix previously used for global transform to identity
  this->Reset();

  // Clean up.
  delete ffdCopy;
}

void irtkMultiLevelFreeFormTransformation::Jacobian(irtkMatrix &jac, double x, double y, double z, double t)
{
  int i;
  irtkMatrix tmp_jac;

  // Compute global jacobian
  this->GlobalJacobian(jac, x, y, z, t);

  // Compute local jacobian
  for (i = 0; i < _NumberOfLevels; i++) {

    // Calculate jacobian
    _localTransformation[i]->Jacobian(tmp_jac, x, y, z, t);

    // Subtract identity matrix
    tmp_jac(0, 0) = tmp_jac(0, 0) - 1;
    tmp_jac(1, 1) = tmp_jac(1, 1) - 1;
    tmp_jac(2, 2) = tmp_jac(2, 2) - 1;

    // Add jacobian
    jac += tmp_jac;
  }
}

void irtkMultiLevelFreeFormTransformation::LocalJacobian(irtkMatrix &jac, double x, double y, double z, double t)
{
  int i;
  irtkMatrix tmp_jac(3, 3);

  // Initialize to identity
  jac.Ident();

  // Compute local jacobian
  for (i = 0; i < _NumberOfLevels; i++) {

    // Calculate jacobian
    _localTransformation[i]->Jacobian(tmp_jac, x, y, z, t);


    // Subtract identity matrix
    tmp_jac(0, 0) = tmp_jac(0, 0) - 1;
    tmp_jac(1, 1) = tmp_jac(1, 1) - 1;
    tmp_jac(2, 2) = tmp_jac(2, 2) - 1;

    // Add jacobian
    jac += tmp_jac;
  }
}

double irtkMultiLevelFreeFormTransformation::Bending(double x, double y, double z)
{
  int i;
  double val = 0.0;
  
  for (i = 0; i < _NumberOfLevels; i++){
    val += _localTransformation[i]->Bending(x, y, z);
  }
  return val;
}

void irtkMultiLevelFreeFormTransformation::ApproximateAsNew(double *x1, double *y1, double *z1, double *x2, double *y2, double *z2, int no)
{
    int level, i;
    double error, rms_error, max_error;

    // Loop over all levels
    rms_error = 0;
    for (level = 0; level < _NumberOfLevels; level++) {

        // Approximate current level
        if (dynamic_cast<irtkBSplineFreeFormTransformation *>(_localTransformation[level]) != NULL) {
            cout << "Approximating FFD at level " << level+1 << endl;
            ((irtkBSplineFreeFormTransformation *)_localTransformation[level])->ApproximateAsNew(x1, y1, z1, x2, y2, z2, no);
        } else {
            cerr << "Only supported for irtkBSplineFreeFormTransformation" << endl;
            return;
        }

        // Calculate error
        rms_error = 0;
        max_error = 0;

        // Loop over points
        for (i = 0;  i < no; i++) {

            // Calculate error
            error = sqrt(x2[i]*x2[i] + y2[i]*y2[i] + z2[i]*z2[i]);

            // Calculate maximum error
            rms_error += error;
            if (error > max_error) {
                max_error = error;
            }
        }

        rms_error /= (double)no;
        cout << "Mean error: " << rms_error << endl;
        cout << "Max. error: " << max_error << endl;
    }
}

double irtkMultiLevelFreeFormTransformation::Approximate(double *x1, double *y1, double *z1, double *x2, double *y2, double *z2, int no)
{
  int level, i;
  double error, rms_error, max_error;

  // Loop over all levels
  rms_error = 0;
  for (level = 0; level < _NumberOfLevels; level++) {

    // Approximate current level
    if (dynamic_cast<irtkFreeFormTransformation3D *>(_localTransformation[level]) != NULL) {
      cout << "Approximating FFD at level " << level+1 << endl;
      ((irtkFreeFormTransformation3D *)_localTransformation[level])->Approximate(x1, y1, z1, x2, y2, z2, no);
    } else {
      cerr << "Not supported for 4D transformations" << endl;
      return 0;
    }

    // Calculate error
    rms_error = 0;
    max_error = 0;

    // Loop over points
    for (i = 0;  i < no; i++) {

      // Calculate error
      error = sqrt(x2[i]*x2[i] + y2[i]*y2[i] + z2[i]*z2[i]);

      // Calculate maximum error
      rms_error += error;
      if (error > max_error) {
        max_error = error;
      }
    }

    rms_error /= (double)no;
    cout << "Mean error: " << rms_error << endl;
    cout << "Max. error: " << max_error << endl;
  }

  // Return error
  return rms_error;
}

int irtkMultiLevelFreeFormTransformation::CheckHeader(char *name)
{
  char buffer[255];

  ifstream from(name);

  if (!from) {
    cerr << "irtkMultiLevelFreeFormTransformation::CheckHeader: Can't open file "
    << name << "\n";
    exit(1);
  }

  // Read keyword
  from >> buffer;
  if ((strcmp(buffer, "MFFD:") != 0) && (strcmp(buffer, "TC:") != 0)) {
    return false;
  }

  return true;
}

istream& irtkMultiLevelFreeFormTransformation::Import(istream &from)
{
  int i;
  char buffer[255];

  // Read keyword
  from >> buffer;

  // Delete old local transformations
  for (i = 0; i < _NumberOfLevels; i++) {
    if (_localTransformation[i] != NULL) delete _localTransformation[i];
  }

  // Check whether this is a transformation in MFFD file format
  if ((strcmp(buffer, "MFFD:") == 0) || (strcmp(buffer, "TC:") == 0)) {

    // Read number of levels
    from >> _NumberOfLevels;

    // Read global transformation
    this->irtkAffineTransformation::Import(from);

    // Read local transformations
    for (i = 0; i < _NumberOfLevels; i++) {

      streampos pos = from.tellg();

      // Read local transformation keyword
      from >> buffer;
      if ((strcmp(buffer, "FFD_BINARY:") == 0) ||
          (strcmp(buffer, "AFFD_BINARY:") == 0)) {

        // Go back and read out local transformation
        from.seekg(pos, ios::beg);
        _localTransformation[i] = new irtkBSplineFreeFormTransformation;
        ((irtkBSplineFreeFormTransformation *)_localTransformation[i])->Import(from);

      } else {
        if ((strcmp(buffer, "LinearFFD_BINARY:") == 0) ||
            (strcmp(buffer, "LinearFFD:") == 0)) {

          // Go back and read out local transformation
          from.seekg(pos, ios::beg);

          _localTransformation[i] = new irtkLinearFreeFormTransformation;
          ((irtkLinearFreeFormTransformation *)_localTransformation[i])->Import(from);

        } else {
          if ((strcmp(buffer, "EFFD:") == 0)) {

            // Go back and read out local transformation
            from.seekg(pos, ios::beg);
            _localTransformation[i] = new irtkEigenFreeFormTransformation;
            ((irtkEigenFreeFormTransformation *)_localTransformation[i])->Import(from);
          } else {
            cerr << "irtkMultiLevelFreeFormTransformation::Import: Unknown file format "
            << buffer << endl;
            exit(1);
          }
        }
      }
    }
  } else {
    cerr << "irtkMultiLevelFreeFormTransformation::Import: Unknown file format " << endl;
    exit(1);
  }

  return from;
}

ostream& irtkMultiLevelFreeFormTransformation::Export(ostream& to)
{
  int i;

  // Write keyword and number of levels
  to << "MFFD: " << _NumberOfLevels << endl;

  // Write global transformation
  this->irtkAffineTransformation::Export(to);

  // Write each local transformation
  for (i = 0; i < _NumberOfLevels; i++) {
    if (strcmp(_localTransformation[i]->NameOfClass(), "irtkBSplineFreeFormTransformation3D") == 0) {
      ((irtkBSplineFreeFormTransformation *)_localTransformation[i])->Export(to);
    } else {
      if (strcmp(_localTransformation[i]->NameOfClass(), "irtkLinearFreeFormTransformation") == 0) {
        ((irtkLinearFreeFormTransformation *)_localTransformation[i])->Export(to);
      } else {
        ((irtkEigenFreeFormTransformation *)_localTransformation[i])->Export(to);
      }
    }
  }

  return to;
}

irtkCifstream& irtkMultiLevelFreeFormTransformation::Read(irtkCifstream& from)
{
  int i, offset;
  unsigned int magic_no, trans_type;

  // Delete old local transformations
  for (i = 0; i < _NumberOfLevels; i++) {
    if (_localTransformation[i] != NULL) delete _localTransformation[i];
  }

  // Read magic no. for transformations
  from.ReadAsUInt(&magic_no, 1);
  if (magic_no != IRTKTRANSFORMATION_MAGIC) {
    cerr << "irtkMultiLevelFreeFormTransformation::Read: Not a valid transformation file" << endl;
    exit(1);
  }

  // Read transformation type
  from.ReadAsUInt(&trans_type, 1);
  if (trans_type != IRTKTRANSFORMATION_MFFD) {
    cerr << "irtkMultiLevelFreeFormTransformation::Read: Not a valid multilevel FFD transformation" << endl;
    exit(1);
  }

  // Read number of local transformations
  from.ReadAsInt(&_NumberOfLevels, 1);

  // Read global transformation
  this->irtkAffineTransformation::Read(from);

  // Read local transformations
  for (i = 0; i < _NumberOfLevels; i++) {

    // Remember current file location
    offset = from.Tell();

    // Read magic no. for transformations
    from.ReadAsUInt(&magic_no, 1);
    if (magic_no != IRTKTRANSFORMATION_MAGIC) {
      cerr << "irtkMultiLevelFreeFormTransformation::Read: Not a valid transformation found at = " << offset << endl;
      exit(1);
    }

    // Read transformation type
    from.ReadAsUInt(&trans_type, 1);
    switch (trans_type) {
      case IRTKTRANSFORMATION_BSPLINE_FFD:
      case IRTKTRANSFORMATION_BSPLINE_FFD_EXT1:
        from.Seek(offset);
        _localTransformation[i] = new irtkBSplineFreeFormTransformation;
        ((irtkBSplineFreeFormTransformation *)_localTransformation[i])->Read(from);
        break;
      case IRTKTRANSFORMATION_BSPLINE_FFD_4D:
        from.Seek(offset);
        _localTransformation[i] = new irtkBSplineFreeFormTransformation4D;
        ((irtkBSplineFreeFormTransformation4D *)_localTransformation[i])->Read(from);
        break;
      case IRTKTRANSFORMATION_EIGEN_FFD:
        from.Seek(offset);
        _localTransformation[i] = new irtkEigenFreeFormTransformation;
        ((irtkEigenFreeFormTransformation *)_localTransformation[i])->Read(from);
        break;
      case IRTKTRANSFORMATION_LINEAR_FFD:
      case IRTKTRANSFORMATION_LINEAR_FFD_EXT1:
        from.Seek(offset);
        _localTransformation[i] = new irtkLinearFreeFormTransformation;
        ((irtkLinearFreeFormTransformation *)_localTransformation[i])->Read(from);
        break;
      case IRTKTRANSFORMATION_PERIODIC:
    	from.Seek(offset);
    	_localTransformation[i] = new irtkBSplineFreeFormTransformationPeriodic;
    	((irtkBSplineFreeFormTransformationPeriodic *)_localTransformation[i])->Read(from);
    	break;
      default:
        cerr << "irtkMultiLevelFreeFormTransformation::Read: No a valid transformation type at = " << offset << endl;
        exit(1);

    }
  }

  return from;
}

irtkCofstream& irtkMultiLevelFreeFormTransformation::Write(irtkCofstream& to)
{
  int i;
  unsigned int magic_no, trans_type;

  // Write magic no. for transformations
  magic_no = IRTKTRANSFORMATION_MAGIC;
  to.WriteAsUInt(&magic_no, 1);

  // Write transformation type
  trans_type = IRTKTRANSFORMATION_MFFD;
  to.WriteAsUInt(&trans_type, 1);

  // Write transformation type
  to.WriteAsInt(&_NumberOfLevels, 1);

  // Write global transformation
  this->irtkAffineTransformation::Write(to);

  // Write local transformations
  for (i = 0; i < _NumberOfLevels; i++) {
    _localTransformation[i]->Write(to);
  }

  return to;
}

bool irtkMultiLevelFreeFormTransformation::IsIdentity()
{
  int i;

  if (this->irtkAffineTransformation::IsIdentity() != true) return false;
  for (i = 0; i < _NumberOfLevels; i++) {
    if (_localTransformation[i]->IsIdentity() != true) return false;
  }
  return true;
}

double irtkMultiLevelFreeFormTransformation::Inverse(double &x, double &y, double &z, double t, double tolerance)
{
  int check;
  double error;

  // Initialize global variables
  irtkTransformationPointer  = this;
  x_invert = x;
  y_invert = y;
  z_invert = z;

  // Pointer to B-spline wrapper
  void (*Newton_function)(int, float [], float []) = irtkTransformationEvaluate;

  // Calculate initial estimate using affine transformation
  this->irtkHomogeneousTransformation::Inverse(x, y, z, t);

  // Inverse
  float invert[3], f_invert[3];
  invert[0] = x;
  invert[1] = y;
  invert[2] = z;

  // Numerically approximate the inverse transformation
  newt2(invert-1, 3, &check, Newton_function);

  // Calculate error
  irtkTransformationEvaluate(3, invert-1, f_invert-1);
  error = sqrt(f_invert[0]*f_invert[0]+f_invert[1]*f_invert[1]+f_invert[2]*f_invert[2]);
  if (error > tolerance) {
    cout << "irtkMultiLevelFreeFormTransformation::Inverse: RMS error = " << error << "\n";
  }

  // Set output to solution
  x = invert[0];
  y = invert[1];
  z = invert[2];

  return error;
}

void irtkMultiLevelFreeFormTransformation::Print()
{
  int i;

  cout << "irtkMultiLevelFreeFormTransformation" << endl;

  // Print global transformations
  this->irtkAffineTransformation::Print();

  // Print local transformations
  for (i = 0; i < _NumberOfLevels; i++) {
    _localTransformation[i]->Print();
  }

}

void irtkMultiLevelFreeFormTransformation::CombineLocalTransformation()
{
  int i;
  irtkFreeFormTransformation *first, *second;

  // Loop over local transformations
  first  = NULL;
  second = NULL;
  while(_NumberOfLevels>1) {
    first  = this->PopLocalTransformation();
    second = this->PopLocalTransformation();
    if(first->NumberOfDOFs() == second->NumberOfDOFs()) {
      for(i = 0; i < first->NumberOfDOFs(); i++) {
        second->Put(i, first->Get(i) + second->Get(i));
      }
    } else {
      cerr << "irtkMultiLevelFreeFormTransformation::CombineLocalTransformation: Only implemented for transformations with equal no. of DOFs" <<endl;
      exit(1);
    }
    this->PushLocalTransformation(second);
    // Delete transformation
    delete first;
  }
}

