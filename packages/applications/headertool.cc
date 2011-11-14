/*=========================================================================

Library   : Image Registration Toolkit (IRTK)
Module    : $Id: headertool.cc 235 2010-10-18 09:25:20Z dr $
Copyright : Imperial College, Department of Computing
Visual Information Processing (VIP), 2008 onwards
Date      : $Date: 2010-10-18 10:25:20 +0100 (一, 18 十月 2010) $
Version   : $Revision: 235 $
Changes   : $Author: dr $

=========================================================================*/

#include <irtkImage.h>

#include <irtkFileToImage.h>

#include <irtkTransformation.h>

char *input_name = NULL, *output_name = NULL, *dof_name  = NULL, *output_header = NULL;

void usage()
{
	cerr << "Usage: headertool [in] [out] <options>\n" << endl;
	cerr << "where <options> can be one or more of the following:\n";
	cerr << "<-size        dx dy dz>                      \t Voxel size   (in mm)\n";
	cerr << "<-tsize       dt>                            \t Voxel size   (in ms)\n";
	cerr << "<-orientation x1 x2 x3   y1 y2 y3  z1 z2 z3> \t Image orientation\n";
	cerr << "<-origin      x  y  z>                       \t Image origin (in mm)\n";
	cerr << "<-torigin     t>                             \t Image origin (in ms)\n";
	cerr << "<-outputheader matrixfilename>               \t Output header in matrix\n";
	cerr << "<-target      image>                         \t Copy header from target image\n\n";
	cerr << "<-dofin       transformation>                \t Apply transformation to axis, spacing and origin information in the" << endl;
	cerr << "	                                            \t header. Transformation may be rigid or affine and with no shearing.\n";
	exit(1);
}

int main(int argc, char **argv)
{
	int i, ok;
	double xsize, ysize, zsize, tsize, xaxis[3], yaxis[3], zaxis[3], origin[4];
	irtkTransformation *transformation = NULL;

	if (argc < 3) {
		usage();
	}

	// Parse filenames
	input_name  = argv[1];
	argc--;
	argv++;
	output_name = argv[1];
	argc--;
	argv++;

	// Read image
	irtkFileToImage *reader = irtkFileToImage::New(input_name);
	irtkBaseImage *image = reader->GetOutput();

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
		if ((ok == false) && (strcmp(argv[1], "-outputheader") == 0)) {
			argc--;
			argv++;
			output_header = argv[1];
			argc--;
			argv++;
			ok = true;
		}
		if ((ok == false) && (strcmp(argv[1], "-target") == 0)) {
			argc--;
			argv++;
			irtkGreyImage target(argv[1]);
			argc--;
			argv++;
			target.GetPixelSize(&xsize, &ysize, &zsize);
			image->PutPixelSize(xsize, ysize, zsize);
			target.GetOrientation(xaxis, yaxis, zaxis);
			image->PutOrientation(xaxis, yaxis, zaxis);
			target.GetOrigin(origin[0], origin[1], origin[2]);
			image->PutOrigin(origin[0], origin[1], origin[2]);
			ok = true;
		}
		if ((ok == false) && (strcmp(argv[1], "-size") == 0)) {
			argc--;
			argv++;
			xsize = atof(argv[1]);
			argc--;
			argv++;
			ysize = atof(argv[1]);
			argc--;
			argv++;
			zsize = atof(argv[1]);
			argc--;
			argv++;
			ok = true;

			// Set voxel size
			image->PutPixelSize(xsize, ysize, zsize);
		}
		if ((ok == false) && (strcmp(argv[1], "-tsize") == 0)) {
			argc--;
			argv++;
			tsize = atof(argv[1]);
			argc--;
			argv++;
			ok = true;

			// Set voxel size
			image->GetPixelSize(&xsize, &ysize, &zsize);
			image->PutPixelSize(xsize, ysize, zsize, tsize);
		}
		if ((ok == false) && (strcmp(argv[1], "-orientation") == 0)) {
			argc--;
			argv++;
			xaxis[0] = atof(argv[1]);
			argc--;
			argv++;
			xaxis[1] = atof(argv[1]);
			argc--;
			argv++;
			xaxis[2] = atof(argv[1]);
			argc--;
			argv++;
			yaxis[0] = atof(argv[1]);
			argc--;
			argv++;
			yaxis[1] = atof(argv[1]);
			argc--;
			argv++;
			yaxis[2] = atof(argv[1]);
			argc--;
			argv++;
			zaxis[0] = atof(argv[1]);
			argc--;
			argv++;
			zaxis[1] = atof(argv[1]);
			argc--;
			argv++;
			zaxis[2] = atof(argv[1]);
			argc--;
			argv++;
			ok = true;

			// Set orientation (third argument now required)
			image->PutOrientation(xaxis, yaxis, zaxis);
		}
		if ((ok == false) && (strcmp(argv[1], "-origin") == 0)) {
			argc--;
			argv++;
			origin[0] = atof(argv[1]);
			argc--;
			argv++;
			origin[1] = atof(argv[1]);
			argc--;
			argv++;
			origin[2] = atof(argv[1]);
			argc--;
			argv++;
			ok = true;

			// Set origin
			image->PutOrigin(origin[0], origin[1], origin[2]);
		}
		if ((ok == false) && (strcmp(argv[1], "-torigin") == 0)) {
			argc--;
			argv++;
			origin[3] = atof(argv[1]);
			argc--;
			argv++;
			ok = true;

			// Set origin
			image->GetOrigin(origin[0], origin[1], origin[2]);
			image->PutOrigin(origin[0], origin[1], origin[2], origin[3]);
		}

		if (ok == false) {
			cout << "Can't parse argument: " << argv[1] << endl;
			usage();
		}
	}

	// Transform header
	if (dof_name != NULL) {

		// Read transformation
		transformation = irtkTransformation::New(dof_name);

		if (strcmp(transformation->NameOfClass(), "irtkAffineTransformation") == 0) {
			irtkAffineTransformation *affineTransf = dynamic_cast<irtkAffineTransformation*> (transformation);

			if (fabs(affineTransf->GetShearXY()) +
					fabs(affineTransf->GetShearXZ()) +
					fabs(affineTransf->GetShearYZ()) > 0.001){
				cerr << "Affine transformation provided contains shearing : Cannot be used to modify header" << endl;
				exit(1);
			}
		} else if (strcmp(transformation->NameOfClass(), "irtkRigidTransformation") != 0) {
			cerr<<"header tool: Can only modify header with a rigid or a no-shear affine transformation"<<endl;
			exit(1);
		}

		irtkImageAttributes attr = image->GetImageAttributes();

		// Origin:

		transformation->Transform(attr._xorigin, attr._yorigin, attr._zorigin);

		// Grid spacings:

		irtkVector v(3);
		irtkVector u(3);

		// Zero vector.
		u.Initialize(3);
		transformation->Transform(u(0), u(1), u(2));

		for (i = 0; i < 3; i++){
			v(i) = attr._xaxis[i];
		}
		transformation->Transform(v(0), v(1), v(2));
		v = v - u;
		attr._dx = attr._dx * v.Norm();

		for (i = 0; i < 3; i++){
			v(i) = attr._yaxis[i];
		}
		transformation->Transform(v(0), v(1), v(2));
		v = v - u;
		attr._dy = attr._dy * v.Norm();

		for (i = 0; i < 3; i++){
			v(i) = attr._zaxis[i];
		}
		transformation->Transform(v(0), v(1), v(2));
		v = v - u;
		attr._dz = attr._dz * v.Norm();

		// Axes:

		// Isolate rotation part of transformation.
		irtkRigidTransformation rotation;
		for (i = 3; i < 6; i++){
			rotation.Put(i, transformation->Get(i));
		}

		rotation.Transform(attr._xaxis[0], attr._xaxis[1], attr._xaxis[2]);
		rotation.Transform(attr._yaxis[0], attr._yaxis[1], attr._yaxis[2]);
		rotation.Transform(attr._zaxis[0], attr._zaxis[1], attr._zaxis[2]);

		// Grid size:

		// Remains the same so no need to do anything.

		// Update image attributes
		image->PutImageAttributes(attr);
	}

	image->Write(output_name);

	if(output_header != NULL){
		irtkMatrix header(4,4);
		header = image->GetImageToWorldMatrix();
		header.Write(output_header);
	}

	return 0;
}

