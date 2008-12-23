/*=========================================================================

  Library   : Image Registration Toolkit (IRTK)
  Module    : $Id$
  Copyright : Imperial College, Department of Computing
              Visual Information Processing (VIP), 2008 onwards
  Date      : $Date$
  Version   : $Revision$
  Changes   : $Author$

=========================================================================*/

#if (defined HAS_VTK && defined HAS_CONTRIB)

#include <stdio.h>
#include <string.h>
#include <irtkEigenAnalysis.h>
#include <irtkGeometry.h>

// Minimum norm of an eigenvector
#define MIN_NORM 0.01

int iNoOfLandmarks;

// vtk includes
#include <vtkPointSet.h>
#include <vtkPoints.h>
#include <vtkPointData.h>
#include <vtkStructuredGrid.h>
#include <vtkPolyData.h>
#include <vtkPolyDataReader.h>
#include <vtkPolyDataWriter.h>
#include <vtkStructuredGridReader.h>
#include <vtkStructuredGridWriter.h>
#include <vtkUnstructuredGrid.h>
#include <vtkUnstructuredGridReader.h>
#include <vtkUnstructuredGridWriter.h>

vtkPointSet *data;
vtkPoints  *points;
vtkDataArray *vectors;

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

void looEigenAnalysis(irtkMatrix& D, int leave, irtkMatrix& Mean, irtkMatrix& Eigenvector, int iNoOfDatasets)
{
  int i, j, k, ii, jj;

  // Form matrix T = (1/iNoOfDatasets) D^T * D which has dimensions
  // iNoOfDatasets-1 x iNoOfDatasets-1
  irtkMatrix T;
  T.Initialize(iNoOfDatasets-1,iNoOfDatasets-1);
  Mean.Initialize(3*iNoOfLandmarks,1);
  Mean *= 0.0;

  // Compute mean
  for (i = 0; i<iNoOfDatasets; i++) {
    if (i == leave)
      continue;
    for (k = 0; k<3*iNoOfLandmarks; k++)
      Mean(k,0) += D(k,i);
  }
  Mean /= iNoOfDatasets-1;

  // Form matrix T = (1/iNoOfDatasets) D^T * D which has dimensions
  // iNoOfDatasets-1 x iNoOfDatasets-1
  for (i = 0; i < iNoOfDatasets; i++) {
    if (i==leave)
      continue;
    ii = (i > leave) ? (i-1) : i;
    for (j = 0; j < iNoOfDatasets; j++) {
      if (j==leave)
        continue;
      jj = (j > leave) ? (j-1) : j;
      T(ii,jj) = 0.0;
      for (k = 0; k < 3*iNoOfLandmarks; k++)
        T(ii,jj) += (D(k,i)-Mean(k,0))*(D(k,j)-Mean(k,0));
    }
  }
  T /= iNoOfDatasets-1;

  // Compute the eigen decomposition
  irtkEigenAnalysis ea(iNoOfDatasets-1);
  for (i = 0; i < iNoOfDatasets-1; i++)
    for (j = 0; j < iNoOfDatasets-1; j++)
      ea.Matrix(i,j) = T(i,j);
  ea.DecrSortEigenStuff();

  // Convert the eigenvectors of T to true eigenvectors of the
  // covariance matrix C = D*D^T.
  Eigenvector.Initialize(3*iNoOfLandmarks,iNoOfDatasets-1);
  Eigenvector *= 0.0;
  for (i = 0; i< iNoOfDatasets-1 ; i++) {
    for (k = 0; k <  3*iNoOfLandmarks; k++) {
      for (j = 0; j < iNoOfDatasets; j++) {
        if (j == leave)
          continue;
        jj = (j > leave) ? (j-1) : j;
        Eigenvector(k,i) += D(k,j)*ea.Eigenvector(jj,i);
      }
    }
  }

  float fTotalVar   = 0;
  float fCummulated = 0;
  for (i = 0; i < iNoOfDatasets-1; i++) {
    fTotalVar += ea.Eigenvalue(i);
  }
  for (i = 0; i < iNoOfDatasets-1; i++) {
    fCummulated += 100 * ea.Eigenvalue(i) / fTotalVar;
  }

  // Normalize eigenvectors
  for (i = 0; i < iNoOfDatasets-1; i++) {
    float fNorm = 0.0;
    for (j = 0; j < 3*iNoOfLandmarks; j++) {
      fNorm += Eigenvector(j,i) * Eigenvector(j,i);
    }
    fNorm = sqrt(fNorm);
    if (100 * ea.Eigenvalue(i) / fTotalVar > MIN_NORM) {
      for (j = 0; j < 3*iNoOfLandmarks; j++) {
        Eigenvector(j,i) /= fNorm;
      }
    } else {
      for (j = 0; j < 3*iNoOfLandmarks; j++) {
        Eigenvector(j,i) = 0;
      }
    }
  }
}

float looReconstruct(irtkMatrix& D, irtkMatrix& Mean, irtkMatrix& Fi,
                     int leave, int modes)
{
  int i, j;
  float error;
  irtkMatrix b, dX, dXhat;

  // Displacement wrt mean: dX = X - Mean
  dX.Initialize(Mean.Rows(), 1);
  dXhat.Initialize(Mean.Rows(), 1);
  b.Initialize(modes, 1);
  for (j = 0; j < 3*iNoOfLandmarks; j++) {
    dX(j,0) = D(j, leave) - Mean(j, 0);
  }

  // Optimal (LSE) parameters: b = Fi^T dX
  for (i = 0; i < modes; i++) {
    b(i, 0) = 0.0;
    for (j = 0; j <  3*iNoOfLandmarks ; j++) {
      b(i, 0) += Fi(j, i) * dX(j, 0);
    }
  }
  // Reconstructed displacements: dXhat = Fi b
  for (j = 0; j < 3*iNoOfLandmarks; j++) {
    dXhat(j, 0) = 0.0;
    for (i = 0; i < modes; i++) {
      dXhat(j, 0) += Fi(j, i) * b(i, 0);
    }
  }

  /*
  // Optimal (LSE) parameters: b = Fi^T dX
  Fi.Transpose();
  b = Fi * dX;
  for (i = modes; i < b.Rows(); i++){
    b(i, 0) = 0;
  }

  // Reconstructed displacements: dXhat = Fi b
  Fi.Transpose();
  dXhat = Fi * b;
  */

  // Reconstruction error: Error = |dX-dXhat|
  error = 0;
  for (j = 0; j < 3*iNoOfLandmarks; j++) {
    error += (dX(j,0)-dXhat(j,0))*(dX(j,0)-dXhat(j,0));
  }
  return sqrt(error/iNoOfLandmarks);
}

int main(int argc, char *argv[])
{
  int i, j;
  int iNoOfDatasets;
  irtkMatrix D, Ev, ev;
  irtkVector MeanShape;

  if ((argc < 3)) {
    cout << " usage: pdm-loo landmarks1 ... landmarks_N \n";
    return 1;
  }

  iNoOfDatasets  = argc-1;
  iNoOfLandmarks = 0;
  cerr << " There are " << iNoOfDatasets << " data sets for training"
       << endl;

  // Read all landmarks in matrix D
  for (i = 0; i < iNoOfDatasets; i++) {
    cerr << " Including landmarks in " << argv[i+1] << endl;
    data    = read(argv[i+1]);
    points  = data->GetPoints();
    vectors = data->GetPointData()->GetVectors();

    if (i == 0) {
      iNoOfLandmarks = points->GetNumberOfPoints();
      D.Initialize(3*iNoOfLandmarks,iNoOfDatasets);
      cerr << " There are " << iNoOfDatasets
           << " datasets with " << points->GetNumberOfPoints()
           << " landmarks." << endl;
    } else {
      if ((points->GetNumberOfPoints()) != iNoOfLandmarks) {
        cerr << "Datasets must contain the same number of landmarks" << endl;
        exit(1);
      }
    }

    for (j = 0; j < iNoOfLandmarks; j++) {
      double p[3], v[3];
      if (vectors == NULL) {
        points->GetPoint(j,p);
        D(3*j,i)   = p[0];
        D(3*j+1,i) = p[1];
        D(3*j+2,i) = p[2];
      } else {
        vectors->GetTuple(j,v);
        D(3*j,i)   = v[0];
        D(3*j+1,i) = v[1];
        D(3*j+2,i) = v[2];
      }
    }
  }

  // Leave one out
  irtkMatrix S, Fi, Mean;

  for (int modes=0; modes<iNoOfDatasets; modes++) {
    float fErrorMean = 0;
    float fErrorSD   = 0;
    for (int leave = 0; leave<iNoOfDatasets; leave++) {
      float Error;
      //      looEigenAnalysis(D, iNoOfDatasets, Mean, Fi, iNoOfDatasets + 1);
      //      Error = looReconstruct(D, Mean, Fi, leave, modes);
      looEigenAnalysis(D, leave, Mean, Fi, iNoOfDatasets);
      Error = looReconstruct(D, Mean, Fi, leave, modes);
      //      cerr << " Error loo " << leave+1 << " = " << Error << endl;
      fErrorMean += Error;
      fErrorSD   += Error*Error;
    }
    fErrorMean /= iNoOfDatasets;
    fErrorSD    = sqrt(fErrorSD/iNoOfDatasets - fErrorMean * fErrorMean);
    cerr << " Avg reconst error for " << modes << " modes = " << fErrorMean
         << " with SD = " << fErrorSD << endl;
  }

}

#else

#include <irtkImage.h>

int main( int argc, char *argv[] )
{
  cerr << argv[0] << " needs to be compiled with the contrib and VTK library " << endl;
}
#endif

