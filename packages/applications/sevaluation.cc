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

#include <vtkPolyData.h>
#include <vtkPolyDataReader.h>
#include <vtkCleanPolyData.h>
#include <vtkPolyDataWriter.h>

#include <irtkTransformation.h>

#include <irtkSurfaceRegistration.h>

char *_target_name = NULL, *_source_name = NULL;
char *dofin_name = NULL, *resultout_name = NULL;

void usage()
{
  cerr << "Usage: sevaluation [target] [source] <options> \n" << endl;
  cerr << "where <options> is one or more of the following:\n" << endl;
  cerr << "<-locator>          Locator: 0 = cell locator, 1 = point locator, 2 = kd-tree locator (default = 1)" << endl;
  cerr << "<-symmetric>        Use symmetric distance (default OFF)" << endl;
  cerr << "<-ignoreedges>      Ignores edges in ICP (default OFF)" << endl;
  cerr << "<-RMS>			   Use rms distance instead of mean distance" << endl;
  cerr << "<-output file>      Result output file" << endl;
  cerr << "<-scalar>           Put the distance as scalar back to target" << endl;
  exit(1);
}

int main(int argc, char **argv)
{
  int i, n, id, locatorType, rms, scalar, ok;
  bool ignoreEdges, symmetricDistance;
  double error, rerror, source_point[3], target_point[3];
  irtkLocator *target_locator, *source_locator;
  irtkTransformation *transformation = NULL;
  vtkDoubleArray *array = NULL;

  if (argc < 3) {
    usage();
  }

  // Default parameters
  locatorType = 1;
  ok = 0;
  rms = 0;
  scalar = 0;
  ignoreEdges = false;
  symmetricDistance = false;

  // Parse filenames
  _target_name = argv[1];
  argv++;
  argc--;
  _source_name = argv[1];
  argv++;
  argc--;

  // Parse remaining parameters
  while (argc > 1) {
    ok = false;
    if ((ok == false) && (strcmp(argv[1], "-locator") == 0)) {
      argc--;
      argv++;
      locatorType = atoi(argv[1]);
      argc--;
      argv++;
      ok = true;
    }
    if ((ok == false) && (strcmp(argv[1], "-ignoreedges") == 0)) {
      argc--;
      argv++;
      ignoreEdges = true;
      ok = true;
    }
    if ((ok == false) && (strcmp(argv[1], "-dofin") == 0)) {
      argc--;
      argv++;
      dofin_name = argv[1];
      argc--;
      argv++;
      ok = true;
    }
    if ((ok == false) && (strcmp(argv[1], "-symmetric") == 0)) {
      argc--;
      argv++;
      symmetricDistance = true;
      ok = true;
    }
	if ((ok == false) && (strcmp(argv[1], "-RMS") == 0)){
		argc--;
		argv++;
		rms = 1;
		ok = true;
	}
	if ((ok == false) && (strcmp(argv[1], "-output") == 0)){
		argc--;
		argv++;
		resultout_name = argv[1];
		argc--;
		argv++;
		ok = true;
	}
	if ((ok == false) && (strcmp(argv[1], "-scalar") == 0)){
		argc--;
		argv++;
		scalar = true;
		ok = true;
	}
    if (ok == false) {
      cerr << "Can not parse argument " << argv[1] << endl;
      usage();
    }
  }

  // Target pipeline
  cout << "Reading target ... " << _target_name << endl;
  vtkPolyDataReader *target_reader = vtkPolyDataReader::New();
  target_reader->SetFileName(_target_name);
  target_reader->Modified();
  target_reader->Update();
  vtkPolyData *target = vtkPolyData::New();
  target = target_reader->GetOutput();
  target->Update();

  // Source pipeline
  cout << "Reading source ... " << _source_name << endl;
  vtkPolyDataReader *source_reader = vtkPolyDataReader::New();
  source_reader->SetFileName(_source_name);
  source_reader->Modified();
  source_reader->Update();
  vtkPolyData *source = vtkPolyData::New();
  source = source_reader->GetOutput();
  source->Update();

  if (ignoreEdges == true) {
    MarkBoundary(target);
    MarkBoundary(source);
  }

  // Check if to do symmetric registration
  if (symmetricDistance == true) {
    // Create target locator
    target_locator = new irtkLocator;
    target_locator->SelectLocatorType(locatorType);
    target_locator->SetDataSet(target);
    // Create source locator
    source_locator = new irtkLocator;
    source_locator->SelectLocatorType(locatorType);
    source_locator->SetDataSet(source);
  } else {
    target_locator = NULL;
    // Create source locator
    source_locator = new irtkLocator;
    source_locator->SelectLocatorType(locatorType);
    source_locator->SetDataSet(source);
  }

  if (ignoreEdges && (target->GetPointData ()->GetScalars () == NULL)) {
    cout << "irtkSurfaceRegistration::Optimize: Asked to ignore edges but no edge outline has been defined" << endl;
    exit (1);
  }

  if(scalar){
	  array = vtkDoubleArray::New();
	  array->SetNumberOfTuples(target->GetNumberOfPoints());
	  array->SetNumberOfComponents(1);
	  array->SetName("DistanceProfile");
  }

  n = 0;
  error = 0;
  for (i = 0; i < target->GetNumberOfPoints(); i++) {
    target->GetPoints()->GetPoint (i, target_point);
    source_point[0] = target_point[0];
    source_point[1] = target_point[1];
    source_point[2] = target_point[2];
    id = source_locator->FindClosestPoint (source_point);
    if (ignoreEdges) {
      if ((id > 0) && (id < source->GetNumberOfPoints()-1)) {
		  rerror = pow(target_point[0] -source_point[0], 2) 
			  + pow(target_point[1] - source_point[1], 2) 
			  +  pow(target_point[2] - source_point[2], 2);
		  if(!rms){
			  rerror = sqrt(rerror);
		  }
		  error += rerror;
		  if(scalar){
			  array->InsertTupleValue(i, &rerror);
		  }
        n++;
      }
    } else {
		rerror = pow(target_point[0] -source_point[0], 2) 
			+ pow(target_point[1] - source_point[1], 2) 
			+  pow(target_point[2] - source_point[2], 2);
		if(!rms){
			rerror = sqrt(rerror);
		}
		error += rerror;
		if(scalar){
			array->InsertTupleValue(i, &rerror);
		}
		n++;
    }
  }
  if (symmetricDistance == true) {
    for (i = 0; i < source->GetNumberOfPoints(); i++) {
      source->GetPoints()->GetPoint (i, source_point);
      target_point[0] = source_point[0];
      target_point[1] = source_point[1];
      target_point[2] = source_point[2];
      id = target_locator->FindClosestPoint (target_point);
      if (ignoreEdges) {
        if ((id > 0) && (id < target->GetNumberOfPoints()-1)) {
			rerror = pow(target_point[0] -source_point[0], 2) 
				+ pow(target_point[1] - source_point[1], 2) 
				+  pow(target_point[2] - source_point[2], 2);
			if(!rms){
				rerror = sqrt(rerror);
			}
			error += rerror;
			n++;
        }
      } else {
		  rerror = pow(target_point[0] -source_point[0], 2) 
			  + pow(target_point[1] - source_point[1], 2) 
			  +  pow(target_point[2] - source_point[2], 2);
		  if(!rms){
			  rerror = sqrt(rerror);
		  }
		  error += rerror;
		  n++;
      }
    }
  }

  if(scalar){
	  target->GetPointData()->SetScalars(array);
	  vtkPolyDataWriter *writer = vtkPolyDataWriter::New();
	  writer->SetInput(target);
	  writer->SetFileName(_target_name);
	  writer->Write();
	  writer->Delete();
	  array->Delete();
  }
  double meanerror,stderror;
  stderror = 0;

  // output to cout
  if(rms == 0){
	  meanerror = error/n;
	  cout << "Mean distance is " << meanerror << " mm" << endl;
	  double tstd = 0;
	   for (i = 0; i < target->GetNumberOfPoints(); i++) {
	  target->GetPoints()->GetPoint (i, target_point);
	  source_point[0] = target_point[0];
	  source_point[1] = target_point[1];
	  source_point[2] = target_point[2];
	  id = source_locator->FindClosestPoint (source_point);
		  
			  tstd = pow(target_point[0] -source_point[0], 2) 
				  + pow(target_point[1] - source_point[1], 2) 
				  +  pow(target_point[2] - source_point[2], 2);
			  tstd = sqrt(tstd);
			  tstd = pow(tstd - meanerror,2);
			  stderror += tstd;
	  }
	  stderror=sqrt(stderror/double(target->GetNumberOfPoints()));
	  cout<<"standard deviation is " << stderror  << " mm" << endl;
  }else{
	  meanerror = sqrt(error/n);
	  cout << "RMS distance is " << meanerror << " mm" << endl;
	  stderror = 0;
  }

  if(resultout_name){
	  cerr << "Writing Results: " << resultout_name << endl;
	  ofstream fout(resultout_name,ios::app);	  
	  fout << meanerror << " " << stderror << endl;
	  fout.close();
  }

  //Final clean up
  if (symmetricDistance == true)
	  delete target_locator;
  delete source_locator;
  target_reader->Delete();
  source_reader->Delete();
}

#else

#include <irtkImage.h>

int main( int argc, char *argv[] )
{
  cerr << argv[0] << " needs to be compiled with the VTK library " << endl;
}
#endif
