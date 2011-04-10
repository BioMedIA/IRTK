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
	cerr << "<-dofin       dof>                           \t Change header from transformation\n";
	cerr << "<-outputheader matrixfilename>               \t Output header in matrix\n";
	cerr << "<-target      image>                         \t Copy header from target image\n\n";
	exit(1);
}

int main(int argc, char **argv)
{
	int i, ok;
	double xsize, ysize, zsize, tsize, xaxis[3], yaxis[3], zaxis[3], origin[4];
	double sx, sy, sz, tansxy, tansxz, tansyz;
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
		irtkMatrix itmat(4, 4);
		irtkMatrix tmp1(4, 4);
		irtkMatrix tmp2(4, 4);
		irtkMatrix transfMat;

		itmat = image->GetImageToWorldMatrix();

		if (strcmp(transformation->NameOfClass(), "irtkRigidTransformation") == 0) {
			irtkRigidTransformation *rigidTransf = dynamic_cast<irtkRigidTransformation*> (transformation);

			transfMat = rigidTransf->GetMatrix();
			transfMat = transfMat * itmat;
			//Scale
			irtkMatrix copy(4, 4);
			copy = transfMat;

			// Get scale manipulating the columns of the upper left 3x3
			// sub-matrix.
			irtkVector col_0, col_1, col_2;
			col_0.Initialize(3);
			col_1.Initialize(3);
			col_2.Initialize(3);
			for (i = 0; i < 3; ++i) {
				col_0(i) = copy(i, 0);
				col_1(i) = copy(i, 1);
				col_2(i) = copy(i, 2);
			}

			// Compute X scale factor and normalize first col.
			sx = col_0.Norm();
			col_0 /= sx;

			// Actually, tansxy and col_1 are still to large by a factor of sy.
			// Now, compute Y scale and normalize 2nd col and rescale tansxy.
			sy = col_1.Norm();
			col_1  /= sy;

			// Actually, tansxz, tansyz and col_2 are still to large by a factor of
			// sz.  Next, get Z scale, normalize 3rd col and scale tansxz and tansyz.
			sz = col_2.Norm();
			col_2  /= sz;

			// At this point, the columns are orthonormal.  Check for a coordinate
			// system flip.  If the determinant is -1, then negate the matrix and the
			// scaling factors.
			/*irtkVector col_1_x_col_2;
			col_1_x_col_2.Initialize(3);
			col_1_x_col_2 = col_1.CrossProduct(col_2);

			if (col_0.ScalarProduct(col_1_x_col_2) < 0) {
				sx *= -1;
				sy *= -1;
				sz *= -1;
				col_0 *= -1;
				col_1 *= -1;
				col_2 *= -1;
			}*/

			xaxis[0] = col_0(0); xaxis[1] = col_0(1); xaxis[2] = col_0(2);
			yaxis[0] = col_1(0); yaxis[1] = col_1(1); yaxis[2] = col_1(2);
			zaxis[0] = col_2(0); zaxis[1] = col_2(1); zaxis[2] = col_2(2);

			//Rotation
			image->PutOrientation(xaxis, yaxis, zaxis);
			//Translation
			image->GetOrigin(origin[0], origin[1], origin[2]);
			rigidTransf->Transform(origin[0], origin[1], origin[2]);
			image->PutOrigin(origin[0], origin[1], origin[2]);
		} else if (strcmp(transformation->NameOfClass(), "irtkAffineTransformation") == 0) {
			irtkAffineTransformation *affineTransf = dynamic_cast<irtkAffineTransformation*> (transformation);

			transfMat = affineTransf->GetMatrix();
			transfMat = transfMat * itmat;
			//Scale
			irtkMatrix copy(4, 4);
			copy = transfMat;

			// Get scale and shear by manipulating the columns of the upper left 3x3
			// sub-matrix.
			irtkVector col_0, col_1, col_2;
			col_0.Initialize(3);
			col_1.Initialize(3);
			col_2.Initialize(3);
			for (i = 0; i < 3; ++i) {
				col_0(i) = copy(i, 0);
				col_1(i) = copy(i, 1);
				col_2(i) = copy(i, 2);
			}

			// Compute X scale factor and normalize first col.
			sx = col_0.Norm();
			col_0 /= sx;

			// Compute XY shear factor and make 2nd col orthogonal to 1st.
			tansxy = col_0.ScalarProduct(col_1);
			col_1 = col_1 - col_0 * tansxy;

			// Actually, tansxy and col_1 are still to large by a factor of sy.
			// Now, compute Y scale and normalize 2nd col and rescale tansxy.
			sy = col_1.Norm();
			col_1  /= sy;
			tansxy /= sy;

			// Compute XZ and YZ shears, orthogonalize 3rd col
			tansxz = col_0.ScalarProduct(col_2);
			col_2 = col_2 - col_0 * tansxz;

			tansyz = col_1.ScalarProduct(col_2);
			col_2 = col_2 - col_1 * tansyz;

			// Actually, tansxz, tansyz and col_2 are still to large by a factor of
			// sz.  Next, get Z scale, normalize 3rd col and scale tansxz and tansyz.
			sz = col_2.Norm();
			col_2  /= sz;
			tansxz /= sz;
			tansyz /= sz;

			xaxis[0] = col_0(0); xaxis[1] = col_0(1); xaxis[2] = col_0(2);
			yaxis[0] = col_1(0); yaxis[1] = col_1(1); yaxis[2] = col_1(2);
			zaxis[0] = col_2(0); zaxis[1] = col_2(1); zaxis[2] = col_2(2);

			//Scale
			image->GetPixelSize(&xsize, &ysize, &zsize);
			image->PutPixelSize(xsize*sx, ysize*sy, zsize*sz);
			//Rotation
			image->PutOrientation(xaxis, yaxis, zaxis);
			//Translation
			image->GetOrigin(origin[0], origin[1], origin[2]);
			affineTransf->Transform(origin[0], origin[1], origin[2]);
			image->PutOrigin(origin[0], origin[1], origin[2]);
		} else {
			cerr<<"header tool can't parse more than affine transformation"<<endl;
			exit(1);
		}



	}

	image->Write(output_name);

	if(output_header != NULL){
		irtkMatrix header(4,4);
		header = image->GetImageToWorldMatrix();
		//irtkImageAttributes attr = image->GetImageAttributes();
		//header.Ident();
		//header(0, 0) = attr._xaxis[0];
		//header(1, 0) = attr._xaxis[1];
		//header(2, 0) = attr._xaxis[2];
		//header(0, 1) = attr._yaxis[0];
		//header(1, 1) = attr._yaxis[1];
		//header(2, 1) = attr._yaxis[2];
		//header(0, 2) = attr._zaxis[0];
		//header(1, 2) = attr._zaxis[1];
		//header(2, 2) = attr._zaxis[2];
		header.Write(output_header);
	}

	return 0;
}

