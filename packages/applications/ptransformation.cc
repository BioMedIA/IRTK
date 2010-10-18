/*=========================================================================

  Library   : Image Registration Toolkit (IRTK)
  Module    : $Id$
  Copyright : Imperial College, Department of Computing
              Visual Information Processing (VIP), 2008 onwards
  Date      : $Date$
  Version   : $Revision$
  Changes   : $Author$

=========================================================================*/

#include <irtkImage.h>

#include <irtkTransformation.h>

#ifdef HAS_VTK

#include <vtkPolyDataWriter.h>
#include <vtkPolyDataReader.h>
#include <vtkPoints.h>
#include <vtkPolyData.h>

// Default filenames
char *input_name = NULL, *output_name = NULL;
char *dof_name  = NULL;

void usage()
{
  cerr << "Usage: ptransformation [source] [output] <options>\n"
       << endl;
  cerr << "<-dofin file>      Transformation" << endl;
  cerr << "<-invert>          Invert transformation" << endl;
  exit(1);
}

int main(int argc, char **argv)
{
  int ok, invert;
  irtkPointSet output, input;
  irtkTransformation *transformation;

  // Check command line
  if (argc < 3) {
    usage();
  }


  // Parse input and output point lists
  input_name = argv[1];
  argc--;
  argv++;
  output_name = argv[1];
  argc--;
  argv++;

  // Default parameters
  invert = false;

  // Parse arguments
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
    if ((ok == false) && (strcmp(argv[1], "-invert") == 0)) {
      argc--;
      argv++;
      invert = true;
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

  // Read point lists
  vtkPolyDataReader *reader = vtkPolyDataReader::New();
  reader->SetFileName(input_name);
  reader->Update();

  vtkPolyData* surface = reader->GetOutput();
  vtkPoints*   points  = surface->GetPoints();

  for (int i=0; i < points->GetNumberOfPoints(); i++) {
    double p[3];
    points->GetPoint(i,p);
    if (invert == false) {
      transformation->Transform(p[0],p[1],p[2]);
    } else {
      transformation->Inverse(p[0],p[1],p[2]);
    }
    points->SetPoint(i,p);
  }

  // Write the final set
  vtkPolyDataWriter   *writer = vtkPolyDataWriter::New();
  writer->SetFileName(output_name);
  writer->SetInput(surface);
  writer->Write();
}

#else

int main(int argc, char **argv)
{
  cerr << "ptransformation: this program needs to be compiled with vtk enabled.\n";
  return 0;
}
#endif
