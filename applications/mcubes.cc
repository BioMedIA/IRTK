/*=========================================================================

  Library   : Image Registration Toolkit (IRTK)
  Module    : $Id$
  Copyright : Imperial College, Department of Computing
              Visual Information Processing (VIP), 2008 onwards
  Date      : $Date$
  Version   : $Revision$
  Changes   : $Author$

=========================================================================*/

#ifdef HAS_VTK

#include <irtkImage.h>

#include <irtkGaussianBlurring.h>
#include <irtkResampling.h>

#include <vtkPointData.h>
#include <vtkPolyData.h>
#include <vtkPolyDataWriter.h>
#include <vtkDecimatePro.h>
#include <vtkSmoothPolyDataFilter.h>
#include <vtkMarchingCubes.h>
#include <vtkPolyDataNormals.h>

char *input_name, *output_name;

void usage()
{

	cerr << "  Usage: mcubes [image] [polydata] [threshold] <options>\n" << endl;

	cerr << "  Run a marching cubes filter on a an image at a given threshold" << endl;
	cerr << "  in order to generate a VTK surface mesh of the boundary.\n" << endl;

	cerr << "  Where <options> are one or more of the following:" << endl;
	cerr << "  <-blur> sigma         Blur input image with kernel size sigma before running filter." << endl;
	cerr << "  <-size> x y z         Resample image to the given voxel size before running." << endl;
	cerr << "  <-isotropic>          Resample image to an isotropic voxel size (minimum of input" << endl;
	cerr << "                        voxel x, y and z dimensions)" << endl;
	cerr << "  <-decimate>           Decimate the vertices in the output surface." << endl;
    cerr << "  <-smooth> iterations  Apply a number of iterations of (Laplacian) smoothing to " << endl;
    cerr << "                        the resulting surface." << endl;
	cerr << "  <-normals> on|off     Choose whether to generate normals (default) or not." << endl;
	cerr << "  <-gradients> on|off   Choose whether to generate gradients or not (default)." << endl;
	cerr << "  <-ascii>              Write output mesh in ASCII format (default:  binary)." << endl;
	cerr << "  <-sepSurfs>           Genarate surfaces for two contours." << endl;
	cerr << "  <-ignoreGeom>         Ignore the geometry of the given image (orientation etc)." << endl;
	cerr << "  <-close>              Put zeros around the image to generate a closed surface(s)."<<endl;
	cerr << "  <-close_x>            Close in x direction only."<<endl;
	cerr << "  <-close_y>            Close in y direction only."<<endl;
	cerr << "  <-close_z>            Close in z direction only."<<endl;
	cerr << " " << endl;
	exit(1);
}

int main(int argc, char **argv)
{
	int ok, i, j, k, bASCII = false, ignoreGeometry, iterations;
	int cx,cy,cz;
	float threshold;
	double xaxis[3], yaxis[3], zaxis[3], point[3];
	irtkPoint origin;

	vtkDecimatePro *decimate = NULL;
	vtkSmoothPolyDataFilter *smooth = NULL;

	ignoreGeometry = 0;
	cx = 0;
	cy = 0;
	cz = 0;
	iterations = 10;

	if (argc < 4) {
		usage();
	}

	// Parse parameters
	input_name = argv[1];
	argv++;
	argc--;
	output_name = argv[1];
	argv++;
	argc--;
	threshold = atof(argv[1]);
	argv++;
	argc--;

	// Read image
	irtkRealImage image;
	image.Read(input_name);

	// Set up marching cubes filter
	vtkMarchingCubes *mcubes = vtkMarchingCubes::New();
	mcubes->SetValue(0, threshold);
	mcubes->ComputeNormalsOn();
	mcubes->ComputeGradientsOff();

	// Parse remaining arguments
	while (argc > 1) {
		ok = false;
		if ((!ok) && (strcmp(argv[1], "-decimate") == 0)) {
			argc--;
			argv++;
			decimate = vtkDecimatePro::New();
			ok = true;
		}
		if ((!ok) && (strcmp(argv[1], "-smooth") == 0)) {
			argc--;
			argv++;
			smooth = vtkSmoothPolyDataFilter::New();
			iterations = atoi(argv[1]);
			argc--;
			argv++;
			ok = true;
		}
		if ((!ok) && (strcmp(argv[1], "-ignoreGeom") == 0)) {
			argc--;
			argv++;
			ignoreGeometry = 1;
			ok = true;
		}
		if ((!ok) && (strcmp(argv[1], "-close") == 0)) {
			argc--;
			argv++;
			cx = 1;
			cy = 1;
			cz = 1;
			ok = true;
		}
		if ((!ok) && (strcmp(argv[1], "-close_x") == 0)) {
			argc--;
			argv++;
			cx = 1;
			ok = true;
		}
		if ((!ok) && (strcmp(argv[1], "-close_y") == 0)) {
			argc--;
			argv++;
			cy = 1;
			ok = true;
		}
		if ((!ok) && (strcmp(argv[1], "-close_z") == 0)) {
			argc--;
			argv++;
			cz = 1;
			ok = true;
		}
		if ((!ok) && (strcmp(argv[1], "-gradients") == 0) && (strcmp(argv[2], "on") == 0)) {
			argc--;
			argv++;
			argc--;
			argv++;
			mcubes->ComputeGradientsOn();
			ok = true;
		}
		if ((!ok) && (strcmp(argv[1], "-gradients") == 0) && (strcmp(argv[2], "off") == 0)) {
			argc--;
			argv++;
			argc--;
			argv++;
			mcubes->ComputeGradientsOff();
			ok = true;
		}
		if ((!ok) && (strcmp(argv[1], "-normals") == 0) && (strcmp(argv[2], "on") == 0)) {
			argc--;
			argv++;
			argc--;
			argv++;
			mcubes->ComputeNormalsOn();
			ok = true;
		}
		if ((!ok) && (strcmp(argv[1], "-normals") == 0) && (strcmp(argv[2], "off") == 0)) {
			argc--;
			argv++;
			argc--;
			argv++;
			mcubes->ComputeNormalsOff();
			ok = true;
		}
		if ((!ok) && (strcmp(argv[1], "-blur") == 0)) {
			argc--;
			argv++;
			// Blur image
			cerr << "Blurring image with sigma = " << atof(argv[1]) << endl;
			irtkGaussianBlurring<irtkRealPixel> gaussianBlurring(atof(argv[1]));
			gaussianBlurring.SetInput (&image);
			gaussianBlurring.SetOutput(&image);
			gaussianBlurring.Run();
			argc--;
			argv++;
			ok = true;
		}
		if ((!ok) && (strcmp(argv[1], "-isotropic") == 0)) {
			argc--;
			argv++;
			// Resample image to isotropic voxels (smalles voxel dimension)
			double xsize, ysize, zsize, size;
			image.GetPixelSize(&xsize, &ysize, &zsize);
			size = xsize;
			size = (size < ysize) ? size : ysize;
			size = (size < zsize) ? size : zsize;
			cerr << "Resampling image to isotropic voxel size (in mm): "
					<< size << endl;
			irtkResampling<irtkRealPixel> resampling(size, size, size);
			irtkLinearInterpolateImageFunction interpolator;
			resampling.SetInput (&image);
			resampling.SetOutput(&image);
			resampling.SetInterpolator(&interpolator);
			resampling.Run();
			ok = true;
		}
		if ((!ok) && (strcmp(argv[1], "-size") == 0)) {
			argc--;
			argv++;
			// Resample image
			cerr << "Resampling image to voxel size (in mm): "
					<< atof(argv[1]) << "," << atof(argv[2]) << "," << atof(argv[3]) << endl;
			irtkResampling<irtkRealPixel> resampling(atof(argv[1]), atof(argv[2]), atof(argv[3]));
			irtkLinearInterpolateImageFunction interpolator;
			resampling.SetInput (&image);
			resampling.SetOutput(&image);
			resampling.SetInterpolator(&interpolator);
			resampling.Run();
			argc -= 3;
			argv += 3;
			ok = true;
		}
		if ((!ok) && (strcmp(argv[1], "-ascii") == 0)) {
			argc--;
			argv++;
			bASCII = true;
			ok = true;
		}
		if ((!ok) && (strcmp(argv[1], "-sepSurfs") == 0)) {
			mcubes->SetNumberOfContours(2);
			mcubes->SetValue(1, threshold+0.1);
		}
		if (!ok) {
			cerr << "Cannot parse argument " << argv[1] << endl;
			usage();
		}
	}

	// Convert image to VTK, taking the differences in vtk/irtk image geometry into account
	vtkStructuredPoints *vtkimage = vtkStructuredPoints::New();

	if(ignoreGeometry == 1){
		irtkImageAttributes tmpatr;
		image.PutOrientation(tmpatr._xaxis,tmpatr._yaxis,tmpatr._zaxis);
		image.PutOrigin(tmpatr._xorigin,tmpatr._yorigin,tmpatr._zorigin);
	}

	irtkRealImage dummy;
	irtkImageAttributes attr = image.GetImageAttributes();
	attr._x += 2*cx;
	attr._y += 2*cy;
	attr._z += 2*cz;

	dummy.Initialize(attr);

	for(i=0; i<dummy.GetX();i++)
		for(j=0; j<dummy.GetY();j++)
			for(k=0; k<dummy.GetZ();k++){
				dummy.Put(i,j,k,threshold-10.0);
			}

	for(i=0; i<image.GetX();i++)
		for(j=0; j<image.GetY();j++)
			for(k=0; k<image.GetZ();k++)
			{
				dummy.Put(i+cx,j+cy,k+cz,image.Get(i,j,k));
			}

	xaxis[0] = 1; xaxis[1] = 0; xaxis[2] = 0;
	yaxis[0] = 0; yaxis[1] = 1; yaxis[2] = 0;
	zaxis[0] = 0; zaxis[1] = 0; zaxis[2] = 1;

	dummy.PutPixelSize(1, 1, 1);
	dummy.PutOrientation(xaxis, yaxis, zaxis);
	dummy.ImageToVTK(vtkimage);

	// Set as input to MC
	mcubes->SetInput(vtkimage);

	// Specify output
	vtkPolyData *output = NULL;

	// Let's go to work
	if (decimate != NULL) {
		cout << "Decimating ... \n";
		decimate->SetInputConnection(mcubes->GetOutputPort());
		if (smooth != NULL) {
			cout << "Smoothing ... \n";
			smooth->SetInputConnection(decimate->GetOutputPort());
			smooth->SetNumberOfIterations(iterations);
			output = smooth->GetOutput();
		} else {
			output = decimate->GetOutput();
		}
	} else if (smooth != NULL) {
		cout << "Smoothing ... \n";
		smooth->SetNumberOfIterations(iterations);
		smooth->SetInputConnection(mcubes->GetOutputPort());
		output = smooth->GetOutput();
	} else {
		output = mcubes->GetOutput();
	}

	output->Update();

	// Now transform between vtk and image coordinate systems
	for (i = 0; i < output->GetNumberOfPoints(); i++) {
		output->GetPoint(i, point);

		dummy.WorldToImage(point[0], point[1], point[2]);

		point[0] = point[0] - cx;
		point[1] = point[1] - cy;
		point[2] = point[2] - cz;

		image.ImageToWorld(point[0], point[1], point[2]);

		output->GetPoints()->SetPoint(i, point);
	}

	output->Modified();

	// Recalculate normals, if applicable
	if (output->GetPointData()->GetNormals() != NULL) {
		vtkPolyDataNormals *filter = vtkPolyDataNormals::New();
		filter->SetInput(output);
		filter->Update(); // absolutely necessary!
		output->GetPointData()->SetNormals(filter->GetOutput()->GetPointData()->GetNormals());
		filter->Delete();
	}

	// Write result
	vtkPolyDataWriter *writer = vtkPolyDataWriter::New();
	writer->SetInput(output);
	writer->SetFileName(output_name);

	if (!bASCII) {
		writer->SetFileTypeToBinary();
	}

	writer->Write();

	// Be good
	if (smooth != NULL) {
		smooth->Delete();
	}
	if (decimate != NULL) {
		decimate->Delete();
	}
	vtkimage->Delete();
	mcubes->Delete();
	writer->Delete();

	return 0;
}

#else

#include <irtkImage.h>

int main( int argc, char *argv[] )
{
	cerr << argv[0] << " this program needs to be compiled with vtk enabled." << endl;
	return 0;
}

#endif

