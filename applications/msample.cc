/*=========================================================================

  Library   : Image Registration Toolkit (IRTK)
  Module    : $Id$
  Copyright : Imperial College, Department of Computing
              Visual Information Processing (VIP), 2009 onwards
  Date      : $Date$
  Version   : $Revision$
  Changes   : $Author$

=========================================================================*/

#include <irtkImage.h>
#include <irtkImageFunction.h>

#ifdef HAS_VTK

#include <vtkPointData.h>
#include <vtkPolyData.h>
#include <vtkPolyDataReader.h>
#include <vtkPolyDataWriter.h>
#include <vtkDoubleArray.h>

char *input_name = NULL, *output_name = NULL, *image_name = NULL;

int main(int argc, char **argv)
{
  int i, j, n;
  double *profile;
  double x, y, z, ds, point[3], normal[3];

  if (argc != 4) {
    cerr << "Usage: msample [input] [image] [output] -n [number of steps (default 10)] -ds [step size (default 1)]" << endl;
    exit(1);
  }

  // Parse filenames
  input_name = argv[1];
  argv++;
  argc--;
  image_name = argv[1];
  argv++;
  argc--;
  output_name = argv[1];
  argv++;
  argc--;

  n  = 10;
  ds = 1;
  
  // Allocate memory for intensity profile
  profile = new double[2*n+1];

  // Read model
  vtkPolyDataReader *reader = vtkPolyDataReader::New();
  reader->SetFileName(input_name);
  reader->Modified();
  reader->Update();
  vtkPolyData *model = vtkPolyData::New();
  model = reader->GetOutput();
  model->Update();

  // Extract normals
  if (model->GetPointData()->GetNormals() == NULL) {
    cerr << "Model has no normals" << endl;
    exit(1);
  }

  // Read image
  irtkGreyImage image;
  image.Read(image_name);
  irtkLinearInterpolateImageFunction interpolator;
  interpolator.SetInput(&image);
  interpolator.Initialize();
  
  // Allocate memory
  vtkDoubleArray *array = vtkDoubleArray::New();
  array->SetNumberOfTuples(model->GetNumberOfPoints());
  array->SetNumberOfComponents(2*n+1);
  array->SetName("IntensityProfile");
  
  for (i = 0; i < model->GetNumberOfPoints(); i++) {
    model->GetPoints()->GetPoint (i, point);
    model->GetPointData()->GetNormals()->GetTuple(i, normal);
    image.WorldToImage(point[0], point[1], point[2]);
    for (j = 0; j < 2*n+1; j++) {
      x = point[0] + (j - n) * ds * normal[0];
      y = point[1] + (j - n) * ds * normal[1];
      z = point[2] + (j - n) * ds * normal[2];
      profile[j] = interpolator.Evaluate(x, y, z);
    }
    array->InsertTupleValue(i, profile);
  }
  model->GetPointData()->SetScalars(array);
  
  vtkPolyDataWriter *writer = vtkPolyDataWriter::New();
  writer->SetInput(model);
  writer->SetFileName(output_name);
  writer->Write();
}

#else

#include <irtkImage.h>

int main( int argc, char *argv[] )
{
  cerr << argv[0] << " this program needs to be compiled with vtk enabled." << endl;
  return 0;
}

#endif
