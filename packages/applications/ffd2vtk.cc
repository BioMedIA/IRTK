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

// vtk includes
#include <vtkPoints.h>
#include <vtkPointData.h>
#include <vtkDataArray.h>
#include <vtkFloatArray.h>
#include <vtkStructuredGrid.h>
#include <vtkStructuredGridWriter.h>

// Default filenames
char *input_name = NULL, *output_name = NULL;

void usage()
{
  cerr << "Usage: ffd2vtk [dof] [vtk] <-displacement>\n" << endl;
  exit(1);
}

int main(int argc, char **argv)
{
  double p[3], v[3];
  int nx, ny, nz, x, y, z, ok, displacement;

  if (argc < 3) {
    usage();
  }

  // Parse file names
  input_name  = argv[1];
  argc--;
  argv++;
  output_name = argv[1];
  argc--;
  argv++;

  // Read multi-level FFD
  irtkMultiLevelFreeFormTransformation *mffd =  new irtkMultiLevelFreeFormTransformation;
  mffd->irtkTransformation::Read(input_name);

  if (mffd->NumberOfLevels() > 1) {
    cerr << "Can't handle FFDs with more than one level" << endl;
    exit(1);
  }

  // Default parameters
  displacement = false;

  // Parse remaining arguments
  while (argc > 1) {
    ok = false;
    if ((ok == false) && (strcmp(argv[1], "-displacement") == 0)) {
      argc--;
      argv++;
      displacement = true;
      ok = true;
    }
    if (ok == false) {
      cerr << "Can't parse argument " << argv[1] << endl;
      usage();
    }
  }

  // Extract FFD and get lattice dimensions
  irtkFreeFormTransformation3D *ffd = dynamic_cast<irtkFreeFormTransformation3D *>(mffd->GetLocalTransformation(0));

  if (ffd == NULL) {
    cerr << "Free-form transformation is not 3D" << endl;
    exit(1);
  }

  nx = ffd->GetX();
  ny = ffd->GetY();
  nz = ffd->GetZ();

  // Set up vtk points
  vtkPoints *points = vtkPoints::New();

  // Set up vtk vectors
  vtkFloatArray *vectors = vtkFloatArray::New();
  vectors->SetNumberOfComponents(3);

  // Remind user what we are doing
  if (displacement == true) {
    cerr << "Converting displacements to VTK" << endl;
  } else {
    cerr << "Converting control points to VTK" << endl;
  }
  // Ignore any translations
  irtkMatrix m = mffd->GetMatrix();

  // Initialize point structure with transformed point positions
  for (z = 0; z < nz; z++) {
    for (y = 0; y < ny; y++) {
      for (x = 0; x < nx; x++) {
        p[0] = x;
        p[1] = y;
        p[2] = z;
        ffd->LatticeToWorld(p[0], p[1], p[2]);
        points->InsertNextPoint(p);
        if (displacement == true) {
          ffd->LocalDisplacement(p[0], p[1], p[2]);
        } else {
          ffd->Get(x, y, z, p[0], p[1], p[2]);
        }
        v[0] = m(0, 0)*p[0] + m(0, 1)*p[1] + m(0, 2)*p[2];
        v[1] = m(1, 0)*p[0] + m(1, 1)*p[1] + m(1, 2)*p[2];
        v[2] = m(2, 0)*p[0] + m(2, 1)*p[1] + m(2, 2)*p[2];
        vectors->InsertNextTuple(v);
      }
    }
  }

  // Allocate objects for vtkStructuredGrid format
  vtkStructuredGrid *grid = vtkStructuredGrid::New();
  vtkStructuredGridWriter *writer = vtkStructuredGridWriter::New();

  // Set structured grid
  grid->SetDimensions(nx, ny, nz);
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
