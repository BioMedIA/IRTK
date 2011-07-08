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
#include <irtkImageFunction.h>
#include <irtkTransformation.h>
#include <vtkPointLocator.h>

// Default filenames
char *cmr1_name = NULL, *cmr2_name = NULL, *input_name = NULL, *output_name = NULL;

void usage()
{
  cerr << "Usage: cardiacsurfacecompare [dof] [vtk] [cmr1] [cmr2] \n"<<             
	      "Compare cmr1 with cmr2 using dof input and output in vtk format\n" << endl;
  exit(1);
}

int main(int argc, char **argv)
{
  int i, x, y, z;
  double p1[3], p2[3];
  vtkIdType j;

  if (argc < 5) {
    usage();
  }

  input_name  = argv[1];
  argc--;
  argv++;
  output_name = argv[1];
  argc--;
  argv++;
  cmr1_name = argv[1];
  argc--;
  argv++;
  cmr2_name = argv[1];
  argc--;
  argv++;

  // Read transformation
  irtkTransformation *transform = irtkTransformation::New(input_name);

  // Read model
  vtkPolyDataReader *reader = vtkPolyDataReader::New();
  reader->SetFileName(cmr1_name);
  reader->Modified();
  reader->Update();
  vtkPolyData *model = vtkPolyData::New();
  model = reader->GetOutput();
  model->Update();

  vtkPolyDataReader *reader_count = vtkPolyDataReader::New();
  reader_count->SetFileName(cmr2_name);
  reader_count->Modified();
  reader_count->Update();
  vtkPolyData *model_count = vtkPolyData::New();
  model_count = reader_count->GetOutput();
  model_count->Update();

  vtkPointLocator *locator = vtkPointLocator::New(); 
  locator->SetDataSet(model_count); // data represents the surface 
  locator->BuildLocator(); 

  for (i = 0; i < model->GetNumberOfPoints(); i++) {
      model->GetPoints()->GetPoint(i,p1);
      transform->Transform(p1[0],p1[1],p1[2]);
      j = locator->FindClosestPoint(p1);
      double difference = *(model->GetPointData()->GetScalars()->GetTuple(i))-*(model_count->GetPointData()->GetScalars()->GetTuple(j));
      model->GetPointData()->GetScalars()->SetTuple(i,&difference);
  }

  vtkPolyDataWriter *writer = vtkPolyDataWriter::New();
  writer->SetInput(model);
  writer->SetFileName(output_name);
  writer->Write();
  writer->Delete();
  reader->Delete();
  reader_count->Delete();
  delete transform;
}

#else

#include <irtkImage.h>

int main( int argc, char *argv[] )
{
  cerr << argv[0] << " needs to be compiled with the VTK library " << endl;
}
#endif
