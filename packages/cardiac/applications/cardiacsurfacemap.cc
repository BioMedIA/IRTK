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
#include <vtkCellLocator.h>

// Default filenames
char *target_name = NULL, *source_name = NULL, *dof_name = NULL, *output_name = NULL;

void usage()
{
  cerr << "Usage: cardiacsurfacecompare [target] [source] [dof] [output] \n"<<             
	      "Map target to source using points in target and location in source with the given dof in vtk format\n" << endl;
  exit(1);
}

int main(int argc, char **argv)
{
  int i, x, y, z;
  double p1[3], p2[3], distance;
  vtkIdType k;
  int j;

  if (argc < 5) {
    usage();
  }

  target_name = argv[1];
  argc--;
  argv++;
  source_name = argv[1];
  argc--;
  argv++;
  dof_name  = argv[1];
  argc--;
  argv++;
  output_name = argv[1];
  argc--;
  argv++;

  // Read transformation
  irtkTransformation *transform = irtkTransformation::New(dof_name);

  // Read model
  vtkPolyDataReader *reader = vtkPolyDataReader::New();
  reader->SetFileName(target_name);
  reader->Modified();
  reader->Update();
  vtkPolyData *model = vtkPolyData::New();
  model = reader->GetOutput();
  model->Update();

  vtkPolyDataReader *reader_count = vtkPolyDataReader::New();
  reader_count->SetFileName(source_name);
  reader_count->Modified();
  reader_count->Update();
  vtkPolyData *model_count = vtkPolyData::New();
  model_count = reader_count->GetOutput();
  model_count->Update();

  vtkCellLocator *locator = vtkCellLocator::New(); 
  locator->SetDataSet(model_count); // data represents the surface 
  locator->BuildLocator(); 

  for (i = 0; i < model->GetNumberOfPoints(); i++) {
      model->GetPoints()->GetPoint(i,p1);
      transform->Transform(p1[0],p1[1],p1[2]);
	  locator->FindClosestPoint(p1,p2,k,j,distance);
	  model->GetPoints()->SetPoint(i,p2);
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
