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

// Default filenames
char *image_name = NULL, *input_name = NULL, *output_name = NULL;

void usage()
{
  cerr << "Usage: dof2vtk [image] [dof] [vtk]\n" << endl;
  exit(1);
}

int main(int argc, char **argv)
{
  int x, y, z;
  double p1[3], p2[3];

  if (argc < 3) {
    usage();
  }

  // Parse file names
  image_name  = argv[1];
  argc--;
  argv++;
  input_name  = argv[1];
  argc--;
  argv++;
  output_name = argv[1];
  argc--;
  argv++;

  // Read image
  irtkGreyImage image(image_name);

  // Read transformation
  irtkTransformation *transform = irtkTransformation::New(input_name);

  // Set up vtk points
  vtkPoints *points = vtkPoints::New();

  // Set up vtk vectors
  vtkFloatArray *vectors = vtkFloatArray::New();
  vectors->SetNumberOfComponents(3);

  // Initialize point structure with transformed point positions
  for (z = 0; z < image.GetZ(); z++) {
    for (y = 0; y < image.GetY(); y++) {
      for (x = 0; x < image.GetX(); x++) {
        p1[0] = x;
        p1[1] = y;
        p1[2] = z;
        image.ImageToWorld(p1[0], p1[1], p1[2]);
        points->InsertNextPoint(p1);
        p2[0] = p1[0];
        p2[1] = p1[1];
        p2[2] = p1[2];
        transform->Transform(p2[0], p2[1], p2[2]);
        p2[0] -= p1[0];
        p2[1] -= p1[1];
        p2[2] -= p1[2];
        vectors->InsertNextTuple(p2);
      }
    }
  }

  // Allocate objects for vtkStructuredGrid format
  vtkStructuredGrid *grid = vtkStructuredGrid::New();
  vtkStructuredGridWriter *writer = vtkStructuredGridWriter::New();

  // Set structured grid
  grid->SetDimensions(image.GetX(), image.GetY(), image.GetZ());
  grid->SetPoints(points);
  grid->GetPointData()->SetVectors(vectors);

  // Write structured grid
  writer->SetInput(grid);
  writer->SetFileName(output_name);
  writer->SetFileTypeToBinary();
  writer->SetVectorsName("vectors");
  writer->Update();

}

#else

#include <irtkImage.h>

int main( int argc, char *argv[] )
{
  cerr << argv[0] << " needs to be compiled with the VTK library " << endl;
}
#endif
