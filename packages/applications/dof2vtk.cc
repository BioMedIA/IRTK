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

#include <vtkFloatArray.h>
#include <vtkPointData.h>
#include <vtkStructuredGrid.h>
#include <vtkStructuredGridWriter.h>

#include <irtkImage.h>
#include <irtkTransformation.h>
#include <irtkResampling.h>

// Default filenames
char *image_name = NULL, *input_name = NULL, *output_name = NULL;

void usage()
{
  cerr << "Usage: dof2vtk [dof] [vtk]                 Load dof from file output dof in world coordinate in vtk format\n" << endl;
  cerr << "<-image image>                             Load image from file output dof in the reference image coordinate\n"										  ;
  cerr << "<-mask value>                              Mask intensity\n";
  cerr << "<-rmatr>                                   Remove Orientation and Origion info\n";
  cerr << "<-isotropic>							      Resample the image according to the smallest resolution Linear interpolator\n";
  cerr << "<-deformation>							  Use deformation of the material point other then deformation on control points\n";
  exit(1);
}

int main(int argc, char **argv)
{
  int x, y, z;
  double p1[3], p2[3];
  int mask,rmatr,isotropic,ok,deformation;

  if (argc < 3) {
    usage();
  }
 
  input_name  = argv[1];
  argc--;
  argv++;
  output_name = argv[1];
  argc--;
  argv++;

  mask = -1;
  rmatr = 0;
  isotropic = 0;
  deformation = 0;

  while (argc > 1) {
	  ok = false;
	  if (strcmp(argv[1], "-rmatr") == 0) {
		  argc--;
		  argv++;
		  rmatr = 1;
		  ok = true;
	  }else if (strcmp(argv[1], "-mask") == 0) {
		  argc--;
		  argv++;
		  mask = atoi(argv[1]);
		  argc--;
		  argv++;
		  ok = true;
	  }else if (strcmp(argv[1], "-image") == 0) {
		  argc--;
		  argv++;
		  image_name = argv[1];
		  argc--;
		  argv++;
		  ok = true;
	  }else if (strcmp(argv[1], "-isotropic") == 0) {
      argc--;
      argv++;
	  isotropic = 1;
	  ok = true;	 
      }
	  else if (strcmp(argv[1], "-deformation") == 0) {
      argc--;
      argv++;
	  deformation = 1;
	  ok = true;	 
      }
	  else if (!ok) {
		  cerr << "Invalid option : " << argv[1] << endl;
		  exit(1);
	  }
  }

  // Read transformation
  irtkTransformation *transform = irtkTransformation::New(input_name);
  // Create initial multi-level free-form deformation
  irtkMultiLevelFreeFormTransformation *mffd = NULL;
  irtkBSplineFreeFormTransformation *affd = NULL;
  irtkGreyImage image,timage;

  // Set up vtk points
  vtkPoints *points = vtkPoints::New();

  // Set up vtk vectors
  vtkFloatArray *vectors = vtkFloatArray::New();
  vectors->SetNumberOfComponents(3);

  if(image_name){
	  // Read image
	  image.Read(image_name);
	  timage.Read(image_name);

  // Remove Attributes
	  if(rmatr == 1){
		  irtkImageAttributes tmpatr;
		  image.PutOrientation(tmpatr._xaxis,tmpatr._yaxis,tmpatr._zaxis);
		  image.PutOrigin(tmpatr._xorigin,tmpatr._yorigin,tmpatr._zorigin);
      }

	  if(isotropic == 1){
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
	  }

	  // Initialize point structure with transformed point positions
	  for (z = 0; z < image.GetZ(); z++) {
		  for (y = 0; y < image.GetY(); y++) {
			  for (x = 0; x < image.GetX(); x++) {
				  p1[0] = x;
				  p1[1] = y;
				  p1[2] = z;
				  timage.ImageToWorld(p1[0], p1[1], p1[2]);
				  p2[0] = p1[0];
				  p2[1] = p1[1];
				  p2[2] = p1[2];
				  transform->Transform(p2[0], p2[1], p2[2]);
				  if(rmatr == 1){
					  p1[0] = x;
					  p1[1] = y;
					  p1[2] = z;
					  image.ImageToWorld(p1[0], p1[1], p1[2]);
					  timage.WorldToImage(p2[0],p2[1],p2[2]);
					  image.ImageToWorld(p2[0],p2[1],p2[2]);
				  }
				  p2[0] -= p1[0];
				  p2[1] -= p1[1];
				  p2[2] -= p1[2];
				  points->InsertNextPoint(p1);
				  if(image.GetAsDouble(x,y,z) > mask || mask == -1){
					  vectors->InsertNextTuple(p2);
				  }else{
					  p2[0] = 0;
					  p2[1] = 0;
					  p2[2] = 0;
					  vectors->InsertNextTuple(p2);
				  }
			  }
		  }
	  }
  }else{
	  if (strcmp(transform->NameOfClass(), "irtkMultiLevelFreeFormTransformation") == 0) {
		  mffd = new irtkMultiLevelFreeFormTransformation(*((irtkMultiLevelFreeFormTransformation *)transform));
		  affd = (irtkBSplineFreeFormTransformation *)mffd->PopLocalTransformation();
	  } else {
		  cerr << "Input transformation is not of type multi-level free form deformation" << endl;
		  exit(1);
	  }
	  // Initialize point structure with transformed point positions
	  for (z = 0; z < affd->GetZ(); z++) {
		  for (y = 0; y < affd->GetY(); y++) {
			  for (x = 0; x < affd->GetX(); x++) {
				  p1[0] = x;
				  p1[1] = y;
				  p1[2] = z;
				  affd->LatticeToWorld(p1[0],p1[1],p1[2]);
				  if(deformation == 1){
					  p2[0] = p1[0];
					  p2[1] = p1[1];
					  p2[2] = p1[2];
					  transform->Transform(p2[0], p2[1], p2[2]);
					  p2[0] -= p1[0];
					  p2[1] -= p1[1];
					  p2[2] -= p1[2];					  				  
				  }else{
					  affd->Get(x,y,z,p2[0],p2[1],p2[2]);
				  }
				  points->InsertNextPoint(p1);
				  vectors->InsertNextTuple(p2);
			  }
		  }
	  }
  }
  // Allocate objects for vtkStructuredGrid format
  vtkStructuredGrid *grid = vtkStructuredGrid::New();
  vtkStructuredGridWriter *writer = vtkStructuredGridWriter::New();

  // Set structured grid
  if(image_name)
	  grid->SetDimensions(image.GetX(), image.GetY(), image.GetZ());
  else
	  grid->SetDimensions(affd->GetX(), affd->GetY(), affd->GetZ());

  grid->SetPoints(points);
  grid->GetPointData()->SetVectors(vectors);

  // Write structured grid
  writer->SetInput(grid);
  writer->SetFileName(output_name);
  writer->SetFileTypeToBinary();
  writer->SetVectorsName("deformation vectors");
  writer->Update();

}

#else

#include <irtkImage.h>

int main( int argc, char *argv[] )
{
  cerr << argv[0] << " needs to be compiled with the VTK library " << endl;
}
#endif
