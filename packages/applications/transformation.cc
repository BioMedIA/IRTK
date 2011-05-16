/*=========================================================================

Library   : Image Registration Toolkit (IRTK)
Module    : $Id$
Copyright : Imperial College, Department of Computing
Visual Information Processing (VIP), 2008 onwards
Date      : $Date$
Version   : $Revision$
Changes   : $Author$

=========================================================================*/

#include <irtkImage.h>
#include <irtkFileToImage.h>

#include <irtkRegionFilter.h>

#include <irtkTransformation.h>

#include <irtkResampling.h>

// Default filenames
char *input_name = NULL, *output_name = NULL, *dof_name  = NULL;

std::string dataTypeName(int dataType)
{
	std::string typeName;

	switch (dataType){
  case IRTK_VOXEL_CHAR:
    typeName = "char";
    break;
  case IRTK_VOXEL_UNSIGNED_CHAR:
    typeName = "unsigned char";
    break;
  case IRTK_VOXEL_SHORT:
    typeName = "short";
    break;
  case IRTK_VOXEL_UNSIGNED_SHORT:
    typeName = "unsigned short";
    break;
  case IRTK_VOXEL_INT:
    typeName = "int";
    break;
  case IRTK_VOXEL_UNSIGNED_INT:
    typeName = "unsigned int";
    break;
  case IRTK_VOXEL_FLOAT:
    typeName = "float";
    break;
  case IRTK_VOXEL_DOUBLE:
    typeName = "double";
    break;
  case IRTK_VOXEL_RGB:
    typeName = "RGB";
    break;
  default:
    typeName = "unknown";
	}

	return typeName;
}


void usage()
{
	cerr << "Usage: transformation [source] [output] <options>\n" << endl;
	cerr << "where <options> is one or more of the following:\n" << endl;
	cerr << "<-dofin file>      Transformation" << endl;
	cerr << "<-target file>     Target image" << endl;
	cerr << "<-2d>              Project transformation in 2D" << endl;
	cerr << "<-Rx1 value>       Region of interest" << endl;
	cerr << "<-Ry1 value>       Region of interest" << endl;
	cerr << "<-Rz1 value>       Region of interest" << endl;
	cerr << "<-Rx2 value>       Region of interest" << endl;
	cerr << "<-Ry2 value>       Region of interest" << endl;
	cerr << "<-Rz2 value>       Region of interest" << endl;
	cerr << "<-Tx1 value>       Region of interest in target image" << endl;
	cerr << "<-Ty1 value>       Region of interest in target image" << endl;
	cerr << "<-Tz1 value>       Region of interest in target image" << endl;
	cerr << "<-Tx2 value>       Region of interest in target image" << endl;
	cerr << "<-Ty2 value>       Region of interest in target image" << endl;
	cerr << "<-Tz2 value>       Region of interest in target image" << endl;
	cerr << "<-Sx1 value>       Region of interest in source image" << endl;
	cerr << "<-Sy1 value>       Region of interest in source image" << endl;
	cerr << "<-Sz1 value>       Region of interest in source image" << endl;
	cerr << "<-Sx2 value>       Region of interest in source image" << endl;
	cerr << "<-Sy2 value>       Region of interest in source image" << endl;
	cerr << "<-Sz2 value>       Region of interest in source image" << endl;
	cerr << "<-Tp value>        Target padding value" << endl;
	cerr << "<-Sp value>        Source padding value" << endl;
	cerr << "<-invert>          Invert transformation" << endl;
	cerr << "<-nn>              Nearst Neighbor interpolation" << endl;
	cerr << "<-linear>          Linear interpolation" << endl;
	cerr << "<-bspline>         B-spline interpolation" << endl;
	cerr << "<-cspline>         Cubic spline interpolation" << endl;
	cerr << "<-sbased>          Shape based interpolation" << endl;
	cerr << "<-sinc>            Sinc interpolation" << endl;
	exit(1);
}

int main(int argc, char **argv)
{
	int ok, invert, twod;
	int target_padding, source_padding;
	int target_x1, target_y1, target_z1, target_x2, target_y2, target_z2;
	int source_x1, source_y1, source_z1, source_x2, source_y2, source_z2;
	irtkTransformation *transformation = NULL;
	irtkImage *source = NULL;
	irtkImage *target = NULL;
	irtkImageFunction *interpolator = NULL;

	int targetType = -1;
	int sourceType = -1;

	// Check command line
	if (argc < 3) {
		usage();
	}

	// Parse image
	input_name  = argv[1];
	argc--;
	argv++;
	output_name = argv[1];
	argc--;
	argv++;

	// Read image
	cout << "Reading image ... "; cout.flush();
	source = irtkImage::New(input_name);
	irtkFileToImage *inputReader = irtkFileToImage::New(input_name);
	sourceType = inputReader->GetDataType();
	cout << "done" << endl;

	// Fix ROI
	target_x1 = -1;
	target_y1 = -1;
	target_z1 = -1;
	target_x2 = -1;
	target_y2 = -1;
	target_z2 = -1;
	source_x1 = -1;
	source_y1 = -1;
	source_z1 = -1;
	source_x2 = -1;
	source_y2 = -1;
	source_z2 = -1;

  // Other options
  invert = false;
  twod = false;
  source_padding = 0;
  target_padding = MIN_GREY;

	while (argc > 1) {
		ok = false;
		if ((ok == false) && (strcmp(argv[1], "-dofin") == 0)) {
			argc--;
			argv++;
			dof_name = argv[1];
			argc--;
			argv++;
			ok = true;
		}
		if ((ok == false) && (strcmp(argv[1], "-target") == 0)) {
			argc--;
			argv++;
			target = irtkImage::New(argv[1]);
			irtkFileToImage *targetReader = irtkFileToImage::New(argv[1]);
			targetType = targetReader->GetDataType();
			argc--;
			argv++;
			ok = true;
		}
		if ((ok == false) && (strcmp(argv[1], "-Rx1") == 0)) {
			argc--;
			argv++;
			target_x1 = source_x1 = atoi(argv[1]);
			argc--;
			argv++;
			ok = true;
		}
		if ((ok == false) && (strcmp(argv[1], "-Rx2") == 0)) {
			argc--;
			argv++;
			target_x2 = source_x2 = atoi(argv[1]);
			argc--;
			argv++;
			ok = true;
		}
		if ((ok == false) && (strcmp(argv[1], "-Ry1") == 0)) {
			argc--;
			argv++;
			target_y1 = source_y1 = atoi(argv[1]);
			argc--;
			argv++;
			ok = true;
		}
		if ((ok == false) && (strcmp(argv[1], "-Ry2") == 0)) {
			argc--;
			argv++;
			target_y2 = source_y2 = atoi(argv[1]);
			argc--;
			argv++;
			ok = true;
		}
		if ((ok == false) && (strcmp(argv[1], "-Rz1") == 0)) {
			argc--;
			argv++;
			target_z1 = source_z1 = atoi(argv[1]);
			argc--;
			argv++;
			ok = true;
		}
		if ((ok == false) && (strcmp(argv[1], "-Rz2") == 0)) {
			argc--;
			argv++;
			target_z2 = source_z2 = atoi(argv[1]);
			argc--;
			argv++;
			ok = true;
		}
		if ((ok == false) && (strcmp(argv[1], "-Tx1") == 0)) {
			argc--;
			argv++;
			target_x1 = atoi(argv[1]);
			argc--;
			argv++;
			ok = true;
		}
		if ((ok == false) && (strcmp(argv[1], "-Tx2") == 0)) {
			argc--;
			argv++;
			target_x2 = atoi(argv[1]);
			argc--;
			argv++;
			ok = true;
		}
		if ((ok == false) && (strcmp(argv[1], "-Ty1") == 0)) {
			argc--;
			argv++;
			target_y1 = atoi(argv[1]);
			argc--;
			argv++;
			ok = true;
		}
		if ((ok == false) && (strcmp(argv[1], "-Ty2") == 0)) {
			argc--;
			argv++;
			target_y2 = atoi(argv[1]);
			argc--;
			argv++;
			ok = true;
		}
		if ((ok == false) && (strcmp(argv[1], "-Tz1") == 0)) {
			argc--;
			argv++;
			target_z1 = atoi(argv[1]);
			argc--;
			argv++;
			ok = true;
		}
		if ((ok == false) && (strcmp(argv[1], "-Tz2") == 0)) {
			argc--;
			argv++;
			target_z2 = atoi(argv[1]);
			argc--;
			argv++;
			ok = true;
		}
		if ((ok == false) && (strcmp(argv[1], "-Sx1") == 0)) {
			argc--;
			argv++;
			source_x1 = atoi(argv[1]);
			argc--;
			argv++;
			ok = true;
		}
		if ((ok == false) && (strcmp(argv[1], "-Sx2") == 0)) {
			argc--;
			argv++;
			source_x2 = atoi(argv[1]);
			argc--;
			argv++;
			ok = true;
		}
		if ((ok == false) && (strcmp(argv[1], "-Sy1") == 0)) {
			argc--;
			argv++;
			source_y1 = atoi(argv[1]);
			argc--;
			argv++;
			ok = true;
		}
		if ((ok == false) && (strcmp(argv[1], "-Sy2") == 0)) {
			argc--;
			argv++;
			source_y2 = atoi(argv[1]);
			argc--;
			argv++;
			ok = true;
		}
		if ((ok == false) && (strcmp(argv[1], "-Sz1") == 0)) {
			argc--;
			argv++;
			source_z1 = atoi(argv[1]);
			argc--;
			argv++;
			ok = true;
		}
		if ((ok == false) && (strcmp(argv[1], "-Sz2") == 0)) {
			argc--;
			argv++;
			source_z2 = atoi(argv[1]);
			argc--;
			argv++;
			ok = true;
		}
		if ((ok == false) && (strcmp(argv[1], "-Tp") == 0)) {
			argc--;
			argv++;
			target_padding = atoi(argv[1]);
			argc--;
			argv++;
			ok = true;
		}
		if ((ok == false) && (strcmp(argv[1], "-Sp") == 0)) {
			argc--;
			argv++;
			source_padding = atoi(argv[1]);
			argc--;
			argv++;
			ok = true;
		}
		if ((ok == false) && (strcmp(argv[1], "-invert") == 0)) {
			argc--;
			argv++;
			invert = true;
			ok = true;
		}
		if ((ok == false) && (strcmp(argv[1], "-2d") == 0)) {
			argc--;
			argv++;
			twod = true;
			ok = true;
		}
		if ((ok == false) && (strcmp(argv[1], "-nn") == 0)) {
			argc--;
			argv++;
			interpolator = new irtkNearestNeighborInterpolateImageFunction;
			ok = true;
		}
		if ((ok == false) && (strcmp(argv[1], "-linear") == 0)) {
			argc--;
			argv++;
			interpolator = new irtkLinearInterpolateImageFunction;
			ok = true;
		}
		if ((ok == false) && (strcmp(argv[1], "-bspline") == 0)) {
			argc--;
			argv++;
			interpolator = new irtkBSplineInterpolateImageFunction;
			ok = true;
		}
		if ((ok == false) && (strcmp(argv[1], "-cspline") == 0)) {
			argc--;
			argv++;
			interpolator = new irtkCSplineInterpolateImageFunction;
			ok = true;
		}
		if ((ok == false) && (strcmp(argv[1], "-sinc") == 0)) {
			argc--;
			argv++;
			interpolator = new irtkSincInterpolateImageFunction;
			ok = true;
		}
		if ((ok == false) && (strcmp(argv[1], "-sbased") == 0)) {
			argc--;
			argv++;
			interpolator = new irtkShapeBasedInterpolateImageFunction;
			ok = true;
		}
		if (ok == false) {
			cerr << "Can not parse argument " << argv[1] << endl;
			usage();
		}
	}

	// Create default interpolator
	if (interpolator == NULL) {
		interpolator = new irtkNearestNeighborInterpolateImageFunction;
	}

	// If there is no target image use copy of source image as target image
	if (target == NULL) {
		target = irtkImage::New(source);
	}else{
		irtkImageAttributes atr = target->GetImageAttributes();
		atr._t = source->GetT();
		target->Initialize(atr);

		if (sourceType != targetType){
			cerr << "\nWarning! source and target image have different data types:" << endl;
			cerr << "Converting data from " << dataTypeName(sourceType) << " to " << dataTypeName(targetType) << ".\n" << endl;
		}
	}

	// Compute region of interest for target image
	if (target_x1 == -1) target_x1 = 0;
	if (target_y1 == -1) target_y1 = 0;
	if (target_z1 == -1) target_z1 = 0;
	if (target_x2 == -1) target_x2 = target->GetX();
	if (target_y2 == -1) target_y2 = target->GetY();
	if (target_z2 == -1) target_z2 = target->GetZ();

	// Compute region of interest for source image
	if (source_x1 == -1) source_x1 = 0;
	if (source_y1 == -1) source_y1 = 0;
	if (source_z1 == -1) source_z1 = 0;
	if (source_x2 == -1) source_x2 = source->GetX();
	if (source_y2 == -1) source_y2 = source->GetY();
	if (source_z2 == -1) source_z2 = source->GetZ();

	// If there is an region of interest for the target image, use it
	if ((target_x1 != 0) || (target_x2 != target->GetX()) ||
		(target_y1 != 0) || (target_y2 != target->GetY()) ||
		(target_z1 != 0) || (target_z2 != target->GetZ())) {
			irtkRegionFilter *region = new irtkRegionFilter;
			region->SetInput (target);
			region->SetOutput(target);
			region->PutRegion(target_x1, target_y1, target_z1, 0, target_x2, target_y2, target_z2, 1);
			region->Run();
	}

	// If there is an region of interest for the source image, use it
	if ((source_x1 != 0) || (source_x2 != source->GetX()) ||
		(source_y1 != 0) || (source_y2 != source->GetY()) ||
		(source_z1 != 0) || (source_z2 != source->GetZ())) {
			irtkRegionFilter *region = new irtkRegionFilter;
			region->SetInput (source);
			region->SetOutput(source);
			region->PutRegion(source_x1, source_y1, source_z1, 0, source_x2, source_y2, source_z2, 1);
			region->Run();
	}

	if (dof_name != NULL) {
		// Read transformation
		transformation = irtkTransformation::New(dof_name);
	} else {
		// Create identity transformation
		transformation = new irtkRigidTransformation;
	}


	// Create image transformation
	irtkImageTransformation *imagetransformation =
		new irtkImageTransformation;

	imagetransformation->SetInput (source, transformation);
	imagetransformation->SetOutput(target);
	imagetransformation->PutTargetPaddingValue(target_padding);
	imagetransformation->PutSourcePaddingValue(source_padding);
	imagetransformation->PutInterpolator(interpolator);

  // Do inverse transformation if necessary
  if (invert == true) {
    imagetransformation->InvertOn();
  } else {
    imagetransformation->InvertOff();
  }

  // Do inverse transformation if necessary
  if (twod == true) {
    imagetransformation->TwoDOn();
  } else {
    imagetransformation->TwoDOff();
  }

	// Transform image
	imagetransformation->Run();

	// Write the final transformation estimate
	target->Write(output_name);
}
