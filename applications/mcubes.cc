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

// vtk includes
#include <vtkPointData.h>
#include <vtkPolyData.h>
#include <vtkPolyDataWriter.h>
//#include <vtkDecimate.h>
#include <vtkDecimatePro.h>
#include <vtkSmoothPolyDataFilter.h>
#include <vtkMarchingCubes.h>
#include <vtkPolyDataNormals.h>

char *input_name, *output_name;

void usage()
{

  cerr << "Usage: mcubes [image] [polydata] [threshold] <-decimate> <-smooth> <-normals on|off> <-gradients on|off> <-blur sigma> <-isotropic> <-size x y z> <-ascii>\n";

  exit(1);
}

int main(int argc, char **argv)
{
  int ok, i, bASCII = False;
  float threshold;
  double xaxis[3], yaxis[3], zaxis[3], point[3];
  irtkPoint origin;

  vtkDecimatePro *decimate = NULL;
  vtkSmoothPolyDataFilter *smooth = NULL;

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
  irtkGreyImage image;
  image.Read(input_name);

  // Set up marching cubes filter
  vtkMarchingCubes *mcubes = vtkMarchingCubes::New();
  mcubes->SetValue(0, threshold);
  mcubes->ComputeNormalsOn();
  mcubes->ComputeGradientsOff();

  // Parse remaining arguments
  while (argc > 1) {
    ok = False;
    if ((!ok) && (strcmp(argv[1], "-decimate") == 0)) {
      argc--;
      argv++;
      decimate = vtkDecimatePro::New();
      ok = True;
    }
    if ((!ok) && (strcmp(argv[1], "-smooth") == 0)) {
      argc--;
      argv++;
      smooth = vtkSmoothPolyDataFilter::New();
      ok = True;
    }
    if ((!ok) && (strcmp(argv[1], "-gradients") == 0) && (strcmp(argv[2], "on") == 0)) {
      argc--;
      argv++;
      argc--;
      argv++;
      mcubes->ComputeGradientsOn();
      ok = True;
    }
    if ((!ok) && (strcmp(argv[1], "-gradients") == 0) && (strcmp(argv[2], "off") == 0)) {
      argc--;
      argv++;
      argc--;
      argv++;
      mcubes->ComputeGradientsOff();
      ok = True;
    }
    if ((!ok) && (strcmp(argv[1], "-normals") == 0) && (strcmp(argv[2], "on") == 0)) {
      argc--;
      argv++;
      argc--;
      argv++;
      mcubes->ComputeNormalsOn();
      ok = True;
    }
    if ((!ok) && (strcmp(argv[1], "-normals") == 0) && (strcmp(argv[2], "off") == 0)) {
      argc--;
      argv++;
      argc--;
      argv++;
      mcubes->ComputeNormalsOff();
      ok = True;
    }
    if ((!ok) && (strcmp(argv[1], "-blur") == 0)) {
      argc--;
      argv++;
      // Blur image
      cerr << "Blurring image with sigma = " << atof(argv[1]) << endl;
      irtkGaussianBlurring<irtkGreyPixel> gaussianBlurring(atof(argv[1]));
      gaussianBlurring.SetInput (&image);
      gaussianBlurring.SetOutput(&image);
      gaussianBlurring.Run();
      argc--;
      argv++;
      ok = True;
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
      irtkResampling<irtkGreyPixel> resampling(size, size, size);
      irtkLinearInterpolateImageFunction interpolator;
      resampling.SetInput (&image);
      resampling.SetOutput(&image);
      resampling.SetInterpolator(&interpolator);
      resampling.Run();
      ok = True;
    }
    if ((!ok) && (strcmp(argv[1], "-size") == 0)) {
      argc--;
      argv++;
      // Resample image
      cerr << "Resampling image to voxel size (in mm): "
           << atof(argv[1]) << "," << atof(argv[2]) << "," << atof(argv[3]) << endl;
      irtkResampling<irtkGreyPixel> resampling(atof(argv[1]), atof(argv[2]), atof(argv[3]));
      irtkLinearInterpolateImageFunction interpolator;
      resampling.SetInput (&image);
      resampling.SetOutput(&image);
      resampling.SetInterpolator(&interpolator);
      resampling.Run();
      argc -= 3;
      argv += 3;
      ok = True;
    }
    if ((!ok) && (strcmp(argv[1], "-ascii") == 0)) {
      argc--;
      argv++;
      bASCII = True;
      ok = True;
    }
    if (!ok) {
      cerr << "Cannot parse argument " << argv[1] << endl;
      usage();
    }
  }

  // Convert image to VTK, taking the differences in vtk/irtk image geometry into account
  vtkStructuredPoints *vtkimage = vtkStructuredPoints::New();
  irtkGreyImage dummy(image);
  xaxis[0] = 1; xaxis[1] = 0; xaxis[2] = 0;
  yaxis[0] = 0; yaxis[1] = 1; yaxis[2] = 0;
  zaxis[0] = 0; zaxis[1] = 0; zaxis[2] = 1;
  dummy.PutPixelSize(1, 1, 1);
  dummy.PutOrigin(0, 0, 0);
  dummy.PutOrientation(xaxis, yaxis, zaxis);
  dummy.ImageToVTK(vtkimage);

  // Set as input to MC
  mcubes->SetInput(vtkimage);

  // Specify output
  vtkPolyData *output = NULL;

  // Let's go to work
  if (decimate != NULL) {
    cout << "Decimating ... \n";
    decimate->SetInput(mcubes->GetOutput());
    if (smooth != NULL) {
      cout << "Smoothing ... \n";
      smooth->SetInput(decimate->GetOutput());
      output = smooth->GetOutput();
    } else {
      output = decimate->GetOutput();
    }
  } else if (smooth != NULL) {
    cout << "Smoothing ... \n";
    smooth->SetInput(mcubes->GetOutput());
    output = smooth->GetOutput();
  } else {
    output = mcubes->GetOutput();
  }
  output->Update(); // just in case

  // Now transform between vtk and image coordinate systems
  for (i = 0; i < output->GetNumberOfPoints(); i++) {
    output->GetPoint(i, point);
    dummy.WorldToImage(point[0], point[1], point[2]);
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
    filter->Delete(); // be good
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

