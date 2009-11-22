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
#include <vtkStructuredGrid.h>
#include <vtkStructuredGridReader.h>


// Default filenames
char *input_name = NULL, *output_name = NULL;

void usage()
{
  cerr << "Usage: vtk2ffd [vtk] [dof] <-interpolate>\n" << endl;
  exit(1);
}

int main(int argc, char **argv)
{
  int i, n, nx, ny, nz, ok, interpolate;
  double x1, y1, z1, x2, y2, z2, dx, dy, dz, xaxis[3], yaxis[3], zaxis[3], v[3];


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

  // Default parameters
  interpolate = False;

  // Parse remaining arguments
  while (argc > 1) {
    ok = False;
    if ((ok == False) && (strcmp(argv[1], "-interpolate") == 0)) {
      argc--;
      argv++;
      interpolate = True;
      ok = True;
    }
    if (ok == False) {
      cerr << "Can't parse argument " << argv[1] << endl;
      usage();
    }
  }

  vtkStructuredGrid *grid = vtkStructuredGrid::New();
  vtkStructuredGridReader *reader = vtkStructuredGridReader::New();
  reader->SetFileName(input_name);
  reader->Update();
  grid = reader->GetOutput();

  // Extract FFD and get lattice dimensions
  nx = grid->GetDimensions()[0];
  ny = grid->GetDimensions()[1];
  nz = grid->GetDimensions()[2];
  n  = nx*ny*nz;
  x1 = grid->GetPoint(0)[0];
  y1 = grid->GetPoint(0)[1];
  z1 = grid->GetPoint(0)[2];
  x2 = grid->GetPoint(n-1)[0];
  y2 = grid->GetPoint(n-1)[1];
  z2 = grid->GetPoint(n-1)[2];

  // Calculate control point spacing
  dx = sqrt(pow(grid->GetPoint(0)[0] - grid->GetPoint(1)[0], 2.0) +
            pow(grid->GetPoint(0)[1] - grid->GetPoint(1)[1], 2.0) +
            pow(grid->GetPoint(0)[2] - grid->GetPoint(1)[2], 2.0));
  dy = sqrt(pow(grid->GetPoint(0)[0] - grid->GetPoint(nx)[0], 2.0) +
            pow(grid->GetPoint(0)[1] - grid->GetPoint(nx)[1], 2.0) +
            pow(grid->GetPoint(0)[2] - grid->GetPoint(nx)[2], 2.0));
  dz = sqrt(pow(grid->GetPoint(0)[0] - grid->GetPoint(nx*ny)[0], 2.0) +
            pow(grid->GetPoint(0)[1] - grid->GetPoint(nx*ny)[1], 2.0) +
            pow(grid->GetPoint(0)[2] - grid->GetPoint(nx*ny)[2], 2.0));

  // Calculate x-axis
  xaxis[0] = (grid->GetPoint(1)[0] - grid->GetPoint(0)[0]) / dx;
  xaxis[1] = (grid->GetPoint(1)[1] - grid->GetPoint(0)[1]) / dx;
  xaxis[2] = (grid->GetPoint(1)[2] - grid->GetPoint(0)[2]) / dx;

  // Calculate z-axis
  yaxis[0] = (grid->GetPoint(nx)[0] - grid->GetPoint(0)[0]) / dy;
  yaxis[1] = (grid->GetPoint(nx)[1] - grid->GetPoint(0)[1]) / dy;
  yaxis[2] = (grid->GetPoint(nx)[2] - grid->GetPoint(0)[2]) / dy;

  // Calculate z-axis
  zaxis[0] = (grid->GetPoint(nx*ny)[0] - grid->GetPoint(0)[0]) / dz;
  zaxis[1] = (grid->GetPoint(nx*ny)[1] - grid->GetPoint(0)[1]) / dz;
  zaxis[2] = (grid->GetPoint(nx*ny)[2] - grid->GetPoint(0)[2]) / dz;

  irtkBSplineFreeFormTransformation *affd = new irtkBSplineFreeFormTransformation(x1, y1, z1, x2, y2, z2, dx, dy, dz, xaxis, yaxis, zaxis);

  irtkMultiLevelFreeFormTransformation *mffd = new irtkMultiLevelFreeFormTransformation;

  if (interpolate == False) {
    // vectors contain control point values
    i = 0;
    vtkDataArray* vectors;
    vectors = grid->GetPointData()->GetVectors();
    for (int z = 0; z < nz; z++) {
      for (int y = 0; y < ny; y++) {
        for (int x = 0; x < nx; x++) {
          affd->Put(x, y, z,
                    vectors->GetTuple(i)[0],
                    vectors->GetTuple(i)[1],
                    vectors->GetTuple(i)[2]);
          i++;
        }
      }
    }
  } else {
    // Allocate memory for displacement values
    double *xdata = new double[n];
    double *ydata = new double[n];
    double *zdata = new double[n];

    // Get pointer to points and vectors
    vtkDataArray* vectors = grid->GetPointData()->GetVectors();

    // Copy displacement values
    for (i = 0; i < n; i++) {
      vectors->GetTuple(i, v);
      xdata[i] = v[0];
      ydata[i] = v[1];
      zdata[i] = v[2];
    }
    affd->Interpolate(xdata, ydata, zdata);
  }

  mffd->PushLocalTransformation(affd);
  mffd->irtkTransformation::Write(output_name);
}

#else

#include <irtkImage.h>

int main( int argc, char *argv[] )
{
  cerr << argv[0] << " needs to be compiled with the VTK library " << endl;
}
#endif
