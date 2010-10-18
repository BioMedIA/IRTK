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

#include <irtkImage.h>

#include <irtkTransformation.h>

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

irtkMatrix Ev;
irtkVector ev;

vtkPointSet *data;
vtkPoints  *points;
vtkDataArray *vectors;

char *dofin_name = NULL, *dofout_name = NULL;
char *eigenvector_name = NULL, *eigenvalue_name = NULL;

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

void usage()
{
  cerr << "Usage: pcmodes [input] [output] [eigenvalue] [eigenvectors] "
       << "<options>\n" << endl;
  cerr << "where <options> is one or more of the following:\n" << endl;
  cerr << "<-mode [i]>         Show i-th mode of variation" << endl;
  cerr << "<-no [i]>           Show i sampels of mode of variation" << endl;
  cerr << "<-min_sigma [i]>    Minimum sigma" << endl;
  cerr << "<-max_sigma [i]>    Maximum sigma" << endl;
  cerr << "<-with-mean>        With mean transformation" << endl;
  cerr << "<-without-mean>     Without mean transformation" << endl;
  exit(1);
}

int main(int argc, char **argv)
{
  int i, j, mode, mean, no, ok;
  double *x, *y, *z, sigma, min_sigma, max_sigma;

  if (argc < 5) {
    usage();
  }

  dofin_name = argv[1];
  argv++;
  argc--;

  dofout_name = argv[1];
  argv++;
  argc--;

  eigenvalue_name = argv[1];
  argv++;
  argc--;

  eigenvector_name = argv[1];
  argv++;
  argc--;

  // Read mean
  cout << "Reading mean from " << dofin_name << endl;
  data    = read(dofin_name);
  points  = data->GetPoints();
  vectors = data->GetPointData()->GetVectors();

  // Read eigenvalues
  cout << "Reading eigenvalues from " << eigenvalue_name << endl;
  ev.Read(eigenvalue_name);

  // Read eigenvectors
  cout << "Reading eigenvectors from " << eigenvector_name << endl;
  Ev.Read(eigenvector_name);

  // Default parameters
  min_sigma = -1;
  max_sigma = +1;
  mode      =  0;
  no        = 20;
  mean      = true;

  // Parse remaining parameters
  while (argc > 1) {
    ok = false;
    if ((ok == false) && (strcmp(argv[1], "-with-mean") == 0)) {
      argc--;
      argv++;
      mean = true;
      ok = true;
    }
    if ((ok == false) && (strcmp(argv[1], "-without-mean") == 0)) {
      argc--;
      argv++;
      mean = false;
      ok = true;
    }
    if ((ok == false) && (strcmp(argv[1], "-min_sigma") == 0)) {
      argc--;
      argv++;
      min_sigma = atof(argv[1]);
      argc--;
      argv++;
      ok = true;
    }
    if ((ok == false) && (strcmp(argv[1], "-max_sigma") == 0)) {
      argc--;
      argv++;
      max_sigma = atof(argv[1]);
      argc--;
      argv++;
      ok = true;
    }
    if ((ok == false) && (strcmp(argv[1], "-mode") == 0)) {
      argc--;
      argv++;
      mode = atoi(argv[1]);
      argc--;
      argv++;
      ok = true;
    }
    if ((ok == false) && (strcmp(argv[1], "-no") == 0)) {
      argc--;
      argv++;
      no = atoi(argv[1]);
      argc--;
      argv++;
      ok = true;
    }
    if (ok == false) {
      cerr << "Can not parse argument " << argv[1] << endl;
      usage();
    }
  }

  if (vectors == NULL) {
    cerr << "Can't find vectors" << endl;

  }

  x = new double[points->GetNumberOfPoints()];
  y = new double[points->GetNumberOfPoints()];
  z = new double[points->GetNumberOfPoints()];
  for (j = 0; j < points->GetNumberOfPoints(); j++) {
    double v[3];
    if (vectors != NULL)
      vectors -> GetTuple(j, v);
    else
      points  -> GetPoint(j, v);

    x[j] = v[0];
    y[j] = v[1];
    z[j] = v[2];
  }




  for (i = 0; i < no+1; i++) {
    // Calculate sigma
    sigma = i * (max_sigma - min_sigma) / double(no) + min_sigma;

    if (mean == true) {
      // Say what we are calculating
      cerr << "Plotting mode = " << mode << " (sigma = " << sigma << ") "
           << "with mean" << endl;
      cerr << "For i: "<< i <<" sigma*sqrt(lamda) is equal:"<<sigma*sqrt(ev(mode));
      for (j = 0; j < points->GetNumberOfPoints(); j++) {
        double v[3];
        v[0] = x[j] + sigma*sqrt(ev(mode))*Ev(3*j,mode);
        v[1] = y[j] + sigma*sqrt(ev(mode))*Ev(3*j+1,mode);
        v[2] = z[j] + sigma*sqrt(ev(mode))*Ev(3*j+2,mode);
        if (vectors!= NULL)
          vectors->SetTuple(j, v);
        else
          points->SetPoint(j,v);
      }
    } else {
      // Say what we are calculating
      cerr << "Plotting mode = " << mode << " (sigma = " << sigma << ") "
           << "without mean" << endl;
      cerr << "For i: "<< i <<" sigma*sqrt(lamda) is equal:"<<sigma*sqrt(ev(mode));
      for (j = 0; j < points->GetNumberOfPoints(); j++) {
        double v[3];
        v[0] = sigma*sqrt(ev(mode))*Ev(3*j,mode);
        v[1] = sigma*sqrt(ev(mode))*Ev(3*j+1,mode);
        v[2] = sigma*sqrt(ev(mode))*Ev(3*j+2,mode);
        if (vectors!=NULL)
          vectors->SetTuple(j, v);
        else
          points->SetPoint(j,v);
      }
    }

    // Write file
    char buffer[255];
    sprintf (buffer, "%s_mode=%d_%.2d_sigma=%+.1f.vtk",
             dofout_name, mode, i, sigma);
    cerr << "Saving to = " << buffer << endl;
    write(buffer, data);
  }
}

#else

#include <irtkImage.h>

int main( int argc, char *argv[] )
{
  cerr << argv[0] << " needs to be compiled with the contrib and vtk library " << endl;
}
#endif
