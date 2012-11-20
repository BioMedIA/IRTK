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

#include <irtkFluidFreeFormTransformation.h>

#ifdef HAS_SUBDIVISION
#include <irtkLatticeFreeFormTransformation.h>
#include <irtkMultiFrameLatticeFreeFormTransformation.h>
#endif


irtkTransformation::irtkTransformation(): _status(NULL)
{
}

irtkTransformation::~irtkTransformation()
{
  // Free memory for DOF status if necessary
  if (_status != NULL) {
    delete []_status;
    _status = NULL;
  }
}


irtkTransformation *irtkTransformation::New(irtkTransformation *transformation)
{
  if (dynamic_cast <irtkRigidTransformation *>(transformation) != NULL) return
      new irtkRigidTransformation(*(dynamic_cast <irtkRigidTransformation *>(transformation)));

  cerr << "irtkTransformation::New: Cannot determine transformation type" << endl;
  return NULL;
}

irtkTransformation *irtkTransformation::New(char *name)
{
  unsigned int magic_no, trans_type;
  irtkTransformation *transformation = NULL;

  // File format is not in old form, so try new file format
  irtkCifstream from;
  from.Open(name);

  // Read magic no. for transformations
  from.ReadAsUInt(&magic_no, 1);
  if (magic_no != IRTKTRANSFORMATION_MAGIC) {

    // Transformation is in old-style format
    cerr << "\nWARNING:\n" << endl;
    cerr << "DOF file is in old format. Trying to convert and import DOF file. Note that this conversion will ";
    cerr << "only work for DOF files produced by voxel-based registration algorithms (rreg, areg, nreg, etc). ";
    cerr << "This conversion will not work for DOF files produced by surface registration algorithms (srreg, ";
    cerr << "sarg, snreg, etc). To convert these DOF files, please use the dofimport tool. \n\n";

    if (irtkRigidTransformation::CheckHeader(name) == true) {
      transformation = new irtkRigidTransformation;
      transformation->Import(name);
      return transformation;
    }

    if (irtkAffineTransformation::CheckHeader(name) == true) {
      transformation = new irtkAffineTransformation;
      transformation->Import(name);
      return transformation;
    }

    if (irtkBSplineFreeFormTransformation::CheckHeader(name) == true) {
      transformation = new irtkBSplineFreeFormTransformation;
      transformation->Import(name);
      return transformation;
    }

    if (irtkLinearFreeFormTransformation::CheckHeader(name) == true) {
      transformation = new irtkLinearFreeFormTransformation;
      transformation->Import(name);
      return transformation;
    }

    if (irtkEigenFreeFormTransformation::CheckHeader(name) == true) {
      transformation = new irtkEigenFreeFormTransformation;
      transformation->Import(name);
      return transformation;
    }

    if (irtkMultiLevelFreeFormTransformation::CheckHeader(name) == true) {
      transformation = new irtkMultiLevelFreeFormTransformation;
      transformation->Import(name);
      return transformation;
    }

#ifdef HAS_SUBDIVISION
    if (irtkLatticeFreeFormTransformation::CheckHeader(name) == true) {
      transformation = new irtkLatticeFreeFormTransformation;
      transformation->Import(name);
      return transformation;
    }

    if (irtkMultiFrameLatticeFreeFormTransformation::CheckHeader(name) == true) {
      transformation = new irtkMultiFrameLatticeFreeFormTransformation;
      transformation->Import(name);
      return transformation;
    }
#endif

    cerr << "irtkTransformation::New: Unknown transformation file format" << endl;
    exit(1);
  }

  // Read transformation type
  from.ReadAsUInt(&trans_type, 1);
  switch (trans_type) {
  case IRTKTRANSFORMATION_HOMOGENEOUS:
    transformation = new irtkHomogeneousTransformation;
    transformation->Read(name);
    return transformation;
  case IRTKTRANSFORMATION_RIGID:
    transformation = new irtkRigidTransformation;
    transformation->Read(name);
    return transformation;
  case IRTKTRANSFORMATION_AFFINE:
    transformation = new irtkAffineTransformation;
    transformation->Read(name);
    return transformation;
  case IRTKTRANSFORMATION_BSPLINE_FFD:
    transformation = new irtkBSplineFreeFormTransformation;
    transformation->Read(name);
    return transformation;
  case IRTKTRANSFORMATION_BSPLINE_FFD_4D:
    transformation = new irtkBSplineFreeFormTransformation4D;
    transformation->Read(name);
    return transformation;
  case IRTKTRANSFORMATION_LINEAR_FFD:
  case IRTKTRANSFORMATION_LINEAR_FFD_EXT1:
     transformation = new irtkLinearFreeFormTransformation;
    transformation->Read(name);
    return transformation;
  case IRTKTRANSFORMATION_MFFD:
    transformation = new irtkMultiLevelFreeFormTransformation;
    transformation->Read(name);
    return transformation;
  case IRTKTRANSFORMATION_FLUID:
    transformation = new irtkFluidFreeFormTransformation;
    transformation->Read(name);
    return transformation;
  case IRTKTRANSFORMATION_EIGEN_FFD:
    transformation = new irtkEigenFreeFormTransformation;
    transformation->Read(name);
    return transformation;
  case IRTKTRANSFORMATION_PERIODIC:
	transformation = new irtkBSplineFreeFormTransformationPeriodic;
	transformation->Read(name);
	return transformation;
#ifdef HAS_SUBDIVISION
  case IRTKTRANSFORMATION_LATTICE_FFD:
    transformation = new irtkLatticeFreeFormTransformation;
    transformation->Read(name);
    return transformation;
  case IRTKTRANSFORMATION_MULTI_FRAME_LATTICE_FFD:
    transformation = new irtkMultiFrameLatticeFreeFormTransformation;
    transformation->Read(name);
    return transformation;
#endif
  default:
    cerr << "irtkTransformation::New: Unknown transformation file format: " << trans_type << endl;
    exit(1);
  }
}

void irtkTransformation::Displacement(irtkGenericImage<double> &image)
{
	int i, j, k;
	double x, y, z;

  // Calculate displacement field
  for (k = 0; k < image.GetZ(); k++) {
    for (j = 0; j < image.GetY(); j++) {
      for (i = 0; i < image.GetX(); i++) {
        x = i;
        y = j;
        z = k;
        // Transform point into world coordinates
        image.ImageToWorld(x, y, z);
        // Calculate displacement
        this->Displacement(x, y, z);
        // Store displacement
        image(i, j, k, 0) = x;
        image(i, j, k, 1) = y;
        image(i, j, k, 2) = z;
       }
    }
  }
}

void irtkTransformation::Read(char *name)
{
  irtkCifstream from;
  unsigned int magic_no;

  // Open file
  from.Open(name);

  // Read magic no. for transformations
  from.ReadAsUInt(&magic_no, 1);
  if (magic_no != IRTKTRANSFORMATION_MAGIC) {
    // Transformation is in old-style format
    cerr << "\nWARNING:\n" << endl;
    cerr << "DOF file is in old format. Trying to convert and import DOF file. Note that this conversion will ";
    cerr << "only work for DOF files produced by voxel-based registration algorithms (rreg, areg, nreg, etc). ";
    cerr << "This conversion will not work for DOF files produced by surface registration algorithms (srreg, ";
    cerr << "sarg, snreg, etc). To convert these DOF files, please use the dofimport tool. \n\n";

    // Try to import file
    this->Import(name);
  } else {

    // Rewind file
    from.Seek(0);

    // Read file
    this->Read(from);
  }
}

void irtkTransformation::Write(char *name)
{
  irtkCofstream to;

  // Open file
  to.Open(name);

  // Export file
  this->Write(to);

  to.Close();
}

void irtkTransformation::Draw()
{
}



