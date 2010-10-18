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

#include <irtkTransformation.h>

#include <vtkPolyData.h>
#include <vtkPointData.h>
#include <vtkPolyDataReader.h>
#include <vtkPolyDataWriter.h>
#include <vtkPolyDataNormals.h>

// Default filenames
char *input_name = NULL, *output_name = NULL, *dof_name  = NULL;

void usage()
{
  cerr << "Usage: transformation [source] [output] <options>\n" << endl;
  cerr << "where <options> is one or more of the following:\n" << endl;
  cerr << "<-dofin file>      Transformation" << endl;
  exit(1);
}

int main(int argc, char **argv)
{
  int i, ok;
  irtkTransformation *transformation = NULL;

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

  // Read surface
  vtkPolyDataReader *surface_reader = vtkPolyDataReader::New();
  surface_reader->SetFileName(input_name);
  surface_reader->Modified();
  surface_reader->Update();
  vtkPolyData *surface = surface_reader->GetOutput();

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
    if (ok == false) {
      cerr << "Can not parse argument " << argv[1] << endl;
      usage();
    }
  }

  if (dof_name != NULL) {
    // Read transformation
    transformation = irtkTransformation::New(dof_name);
  } else {
    // Create identity transformation
    transformation = new irtkRigidTransformation;
  }

  if (surface->GetPoints()->GetData()->GetNumberOfComponents() == 3) {
    double coord[3];
    for (i = 0; i < surface->GetNumberOfPoints(); i++) {
      (surface->GetPoints())->GetPoint(i, coord);
      transformation->Transform(coord[0], coord[1], coord[2]);
      surface->GetPoints()->SetPoint(i, coord);
    }
    surface->Modified();
  } else if (surface->GetPoints()->GetData()->GetNumberOfComponents() == 4) {
    double coord[4];
    for (i = 0; i < surface->GetNumberOfPoints(); i++) {
      (surface->GetPoints())->GetPoint(i, coord);
      transformation->Transform(coord[0], coord[1], coord[2]);
      surface->GetPoints()->SetPoint(i, coord);
    }
    surface->Modified();
  }

  // Recalculate normals, if applicable
  if (surface->GetPointData()->GetNormals() != NULL) {
    vtkPolyDataNormals *filter = vtkPolyDataNormals::New();
    filter->SetInput(surface);
    filter->Update(); // absolutely necessary!
    surface->GetPointData()->SetNormals(filter->GetOutput()->GetPointData()->GetNormals());
    filter->Delete(); // be good
  }

  vtkPolyDataWriter *writer = vtkPolyDataWriter::New();
  writer->SetInput(surface);
  writer->SetFileName(output_name);
  writer->Write();
}

#else

#include <irtkImage.h>

int main( int argc, char *argv[] )
{
  cerr << argv[0] << " needs to be compiled with the VTK library " << endl;
}
#endif
