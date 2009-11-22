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
#include <irtkImage.h>
#include <irtkEigenAnalysis.h>
#include <irtkTransformation.h>

// Minimum norm of an eigenvector
#define MIN_NORM 0.01

// vtk includes
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

void write(char *file, vtkPointSet *pset)
{
  if (pset->IsA("vtkStructuredGrid")) {
    vtkStructuredGridWriter *writer = vtkStructuredGridWriter::New();
    writer->SetInput((vtkStructuredGrid *)pset);
    writer->SetFileName(file);
    writer->SetFileTypeToBinary();
    writer->Update();
  } else {
    if (pset->IsA("vtkUnstructuredGrid")) {
      vtkUnstructuredGridWriter *writer = vtkUnstructuredGridWriter::New();
      writer->SetInput((vtkUnstructuredGrid *)pset);
      writer->SetFileName(file);
      writer->SetFileTypeToBinary();
      writer->Update();
    } else {
      if (pset->IsA("vtkPolyData")) {
        vtkPolyDataWriter *writer = vtkPolyDataWriter::New();
        writer->SetInput((vtkPolyData *)pset);
        writer->SetFileName(file);
        writer->SetFileTypeToBinary();
        writer->Update();
      } else {
        cerr << "Unknown VTK data type" << endl;
        exit(1);
      }
    }
  }
}

int main(int argc, char **argv)
{
  int i, j, k, iNoOfLandmarks, iNoOfDatasets;
  irtkMatrix M, T, Eigenvector;
  irtkVector MeanShape, Eigenvalues;

  if ((argc < 5)) {
    cout << "usage: pcanalysis landmarks1 ... landmarks_N landmarks_mean "
         << "eigen_values eigen_vectors\n";
    return 1;
  }

  iNoOfDatasets  = argc-4;
  iNoOfLandmarks = 0;
  cerr << " There are " << iNoOfDatasets << " data sets for training"
       << endl;

  // Read all landmarks in matrix M
  for (i = 0; i < iNoOfDatasets; i++) {
    cerr << " Including landmarks in " << argv[i+1] << endl;
    data    = read(argv[i+1]);
    points  = data->GetPoints();
    vectors = data->GetPointData()->GetVectors();

    if (i == 0) {
      iNoOfLandmarks = points->GetNumberOfPoints();
      M.Initialize(3*iNoOfLandmarks,iNoOfDatasets);
      MeanShape.Initialize(3*iNoOfLandmarks);
      for (j = 0; j < 3*iNoOfLandmarks; j++) {
        MeanShape(j) = 0.0;
      }
      cerr << " There are " << iNoOfDatasets
           << " datasets with "
           << points->GetNumberOfPoints()
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
        M(3*j,i)   = p[0];
        M(3*j+1,i) = p[1];
        M(3*j+2,i) = p[2];
        MeanShape(3*j)   += p[0];
        MeanShape(3*j+1) += p[1];
        MeanShape(3*j+2) += p[2];
      } else {
        vectors->GetTuple(j,v);
        M(3*j,i)   = v[0];
        M(3*j+1,i) = v[1];
        M(3*j+2,i) = v[2];
        MeanShape(3*j)   += v[0];
        MeanShape(3*j+1) += v[1];
        MeanShape(3*j+2) += v[2];
      }
    }
  }


  MeanShape = MeanShape*(1.0/iNoOfDatasets);

  // Substract the mean
  for (i = 0; i < iNoOfDatasets; i++) {
    for (j = 0; j < iNoOfLandmarks; j++) {
      M(3*j,i)    -= MeanShape(3*j);
      M(3*j+1,i)  -= MeanShape(3*j+1);
      M(3*j+2,i)  -= MeanShape(3*j+2);
    }
  }

  // Form matrix T = (1/iNoOfDatasets) M^T * M which has dimensions
  // iNoOfDatasets x iNoOfDatasets
  T.Initialize(iNoOfDatasets,iNoOfDatasets);
  for (i = 0; i < iNoOfDatasets; i++) {
    for (j = 0; j < iNoOfDatasets; j++) {
      T(i, j) = 0.0;
      for (k = 0; k < 3 * iNoOfLandmarks; k++)
        T(i, j) += M(k, i)*M(k, j);
    }
  }
  T /= iNoOfDatasets;

  // Compute the eigen decomposition
  irtkEigenAnalysis ea(iNoOfDatasets);
  for (i = 0; i < iNoOfDatasets; i++) {
    for (int j=0; j < iNoOfDatasets; j++) {
      ea.Matrix(i,j) = T(i,j);
    }
  }
  ea.DecrSortEigenStuff();

  // Convert the eigenvectors of T to true eigenvectors of the
  // covariance matrix C = M*M^T.
  Eigenvector.Initialize(3*iNoOfLandmarks,iNoOfDatasets);
  Eigenvector *= 0.0;
  for (i = 0; i < iNoOfDatasets ; i++) {
    for (j = 0; j < 3*iNoOfLandmarks; j++) {
      for (k = 0; k < iNoOfDatasets; k++) {
        Eigenvector(j,i) += M(j,k)*ea.Eigenvector(k,i);
      }
    }
  }

  Eigenvalues.Initialize(iNoOfDatasets);
  for (i = 0; i < iNoOfDatasets ; i++) {
    Eigenvalues(i) = ea.Eigenvalue(i);
  }

  float fTotalVar   = 0;
  float fCummulated = 0;
  for (i = 0; i < iNoOfDatasets; i++) {
    fTotalVar += ea.Eigenvalue(i);
  }
  for (i = 0; i < iNoOfDatasets; i++) {
    cout<< "Mode:" << i << " Eigenvalue:"<<ea.Eigenvalue(i)<<endl;

    fCummulated += 100 * ea.Eigenvalue(i) / fTotalVar;
    cout << " Mode [" << i << "] explains "
         << 100 * ea.Eigenvalue(i) / fTotalVar
         << " % ("
         << fCummulated
         << " %) of shape variance." << endl;
  }

  // Normalize eigenvectors
  for (i = 0; i < iNoOfDatasets ; i++) {
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
      Eigenvalues(i) = 0;
    }
  }

  cerr << " Writing eigenvalues to " << argv[argc-2] << endl;
  Eigenvalues.Write(argv[argc-2]);

  cerr << " Writing eigenshapes to " << argv[argc-1] << endl;
  Eigenvector.Write(argv[argc-1]);

  // Convert mean back to VTK
  for (j = 0; j < iNoOfLandmarks; j++) {
    double p[3], v[3];
    if (vectors == NULL) {
      p[0] = MeanShape(3*j);
      p[1] = MeanShape(3*j+1);
      p[2] = MeanShape(3*j+2);
      points->SetPoint(j,p);
    } else {
      v[0] = MeanShape(3*j);
      v[1] = MeanShape(3*j+1);
      v[2] = MeanShape(3*j+2);
      vectors->SetTuple(j,v);
    }
  }

  cerr << " Writing mean shape to " << argv[argc-3] << endl;
  write(argv[argc-3], data);
}

#else

#include <irtkImage.h>

int main( int argc, char *argv[] )
{
  cerr << argv[0] << " needs to be compiled with the contrib and VTK library"
       << endl;
}
#endif
