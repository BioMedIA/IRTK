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

// vtk includes
#include <vtkPoints.h>
#include <vtkStructuredGrid.h>
#include <vtkPointData.h>
#include <vtkPolyData.h>
#include <vtkPolyDataReader.h>
#include <vtkStructuredGridReader.h>
#include <vtkUnstructuredGrid.h>
#include <vtkUnstructuredGridReader.h>

vtkPointSet *data;
vtkPoints  *points;
vtkDataArray *vectors;

char *dofin_name, *dofout_name, *vtk_name;

vtkPointSet *read(char *file)
{
  int i;
  char buffer[256];
  vtkPointSet *pset;

  ifstream is(file);
  if (!is) {
    cerr << "Can't open file " << file << endl;
    exit(1);
  }

  for (i = 0; i < 3; i++) {
    is.getline(buffer, 256);
  }
  is >> buffer >> buffer;

  if (strcmp(buffer, "POLYDATA") == 0) {

    // Read vtkPolyData object
    vtkPolyDataReader *reader = vtkPolyDataReader::New();
    reader->SetFileName(file);
    reader->Update();
    pset = reader->GetOutput();
    pset->Register(pset);
    reader->Delete();
  } else {
    if (strcmp(buffer, "UNSTRUCTURED_GRID") == 0) {
      // Read vtkUnstructuredGrid object
      vtkUnstructuredGridReader *reader = vtkUnstructuredGridReader::New();
      reader->SetFileName(file);
      reader->Update();
      pset = reader->GetOutput();
      pset->Register(pset);
      reader->Delete();
    } else {
      if (strcmp(buffer, "STRUCTURED_GRID") == 0) {
        // Read vtkStructuredGrid object
        vtkStructuredGridReader *reader = vtkStructuredGridReader::New();
        reader->SetFileName(file);
        reader->Update();
        pset = reader->GetOutput();
        pset->Register(pset);
        reader->Delete();
      } else {
        cerr << "Unknown VTK data type" << endl;
        exit(1);
      }
    }
  }
  return pset;
}

void usage()
{
  cerr << "Usage: approximate [vtk] [dofin] [dofout] <options>\n";
  cerr << "where <options> is one or more of the following:\n" << endl;
  cerr << "<-invert>       Approximate inverse" << endl;
  exit(1);
}

int main(int argc, char **argv)
{
  int i, invert, ok;
  double p[3], v[3];

  // Check command line
  if (argc < 4) {
    usage();
  }

  vtk_name    = argv[1];
  argv++;
  argc--;
  dofin_name  = argv[1];
  argv++;
  argc--;
  dofout_name = argv[1];
  argv++;
  argc--;

  // Parse remaining parameters
  invert = false;
  while (argc > 1) {
    ok = false;
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

  // Read VTK data
  data = read(vtk_name);

  // Read transformation
  irtkMultiLevelFreeFormTransformation *mffd =
    new irtkMultiLevelFreeFormTransformation;
  mffd->irtkTransformation::Read(dofin_name);

  // Get points and vectors
  points  = data->GetPoints();
  vectors = data->GetPointData()->GetVectors();
  if (vectors == NULL) {
    cerr << "approximate: displacement vectors missing" << endl;
    exit(1);
  }

  // Allocate some memory
  double *x1 = new double[points->GetNumberOfPoints()];
  double *y1 = new double[points->GetNumberOfPoints()];
  double *z1 = new double[points->GetNumberOfPoints()];
  double *x2 = new double[points->GetNumberOfPoints()];
  double *y2 = new double[points->GetNumberOfPoints()];
  double *z2 = new double[points->GetNumberOfPoints()];

  // Say what we are doing
  if (invert == true) {
    cout << "Approximating inverse transformation" << endl;
  } else {
    cout << "Approximating transformation" << endl;
    // Invert global transformation
    mffd->Invert();
  }

  // Convert points
  for (i = 0; i < points->GetNumberOfPoints(); i++) {
    points->GetPoint(i, p);
    vectors->GetTuple(i, v);

    // Take global transformation into account
    v[0] += p[0];
    v[1] += p[1];
    v[2] += p[2];
    mffd->GlobalTransform(p[0], p[1], p[2]);
    v[0] = v[0] - p[0];
    v[1] = v[1] - p[1];
    v[2] = v[2] - p[2];
    points->GetPoint(i, p);

    if (invert == true) {
      x1[i] = p[0] + v[0];
      y1[i] = p[1] + v[1];
      z1[i] = p[2] + v[2];
      x2[i] = -v[0];
      y2[i] = -v[1];
      z2[i] = -v[2];
    } else {
      x1[i] = p[0];
      y1[i] = p[1];
      z1[i] = p[2];
      x2[i] = v[0];
      y2[i] = v[1];
      z2[i] = v[2];
    }
  }

  // Invert global transformation back
  mffd->Invert();

  // Approximate transformations
  mffd->Approximate(x1, y1, z1, x2, y2, z2, points->GetNumberOfPoints());

  // Write transformation collection
  mffd->irtkTransformation::Write(dofout_name);
}

#else

int main(int argc, char **argv)
{
  cerr << "approximate: this program needs to be compiled with vtk enabled.\n";
  return 0;
}
#endif
