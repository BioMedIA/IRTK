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

int main(int argc, char **argv)
{
  irtkMatrix jac;
  int i, j, k, nx, ny, nz;
  double x, y, z, v1[3];
  float p[3], v2[3];

  // Read first transformation
  irtkMultiLevelFreeFormTransformation *mffd1 = new irtkMultiLevelFreeFormTransformation;
  mffd1->irtkTransformation::Read(argv[1]);

  // Read second transformation
  irtkMultiLevelFreeFormTransformation *mffd2 = new irtkMultiLevelFreeFormTransformation;
  mffd2->irtkTransformation::Read(argv[2]);

  // Read third transformation
  irtkMultiLevelFreeFormTransformation *mffd3 = new irtkMultiLevelFreeFormTransformation;
  mffd3->irtkTransformation::Read(argv[3]);

  // Extract FFD and get lattice dimensions
  irtkFreeFormTransformation3D *affd = dynamic_cast<irtkFreeFormTransformation3D *>(mffd3->GetLocalTransformation(0));

  if (affd == NULL) {
    cerr << "Free-form transformation is not 3D" << endl;
    exit(1);
  }
  nx = affd->GetX();
  ny = affd->GetY();
  nz = affd->GetZ();

  // Set up vtk
  vtkPoints *points = vtkPoints::New();
  vtkDataArray *vectors = vtkFloatArray::New();

  // Ignore any translations
  irtkMatrix m = mffd2->GetMatrix();

  // Loop over points in target
  for (k = 0; k < nz; k++) {
    for (j = 0; j < ny; j++) {
      for (i = 0; i < nx; i++) {
        x = i;
        y = j;
        z = k;

        // Transform point from lattice coordinates to target coordinates
        affd->LatticeToWorld(x, y, z);

        p[0] = x;
        p[1] = y;
        p[2] = z;
        points->InsertNextPoint(p);

        // Calculate total jacobian at the point
        double val = mffd1->irtkTransformation::Jacobian(x, y, z);

        if ((val > 0.1) && (val < 5)) {

          mffd1->Jacobian(jac, x, y, z);
          jac.Invert();

          // Transform target world coordinates to source coordinates
          mffd1->Transform(x, y, z);

          // Store world coords
          v1[0] = x;
          v1[1] = y;
          v1[2] = z;

          // Tranform point
          mffd2->LocalDisplacement(v1[0], v1[1], v1[2]);

          v2[0] = m(0, 0)*v1[0] + m(0, 1)*v1[1] + m(0, 2)*v1[2];
          v2[1] = m(1, 0)*v1[0] + m(1, 1)*v1[1] + m(1, 2)*v1[2];
          v2[2] = m(2, 0)*v1[0] + m(2, 1)*v1[1] + m(2, 2)*v1[2];
          v1[0] = v2[0];
          v1[1] = v2[1];
          v1[2] = v2[2];

          // Transform deformation vector
          v2[0]=jac(0,0)*v1[0]+jac(0,1)*v1[1]+jac(0,2)*v1[2];
          v2[1]=jac(1,0)*v1[0]+jac(1,1)*v1[1]+jac(1,2)*v1[2];
          v2[2]=jac(2,0)*v1[0]+jac(2,1)*v1[1]+jac(2,2)*v1[2];

        } else {
          v2[0] = 0;
          v2[1] = 0;
          v2[2] = 0;
        }
        vectors->InsertNextTuple(v2);
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
  writer->SetFileName(argv[4]);
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


