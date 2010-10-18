#ifdef HAS_VTK

#include <irtkRegistration.h>

// Default filenames
char *source_name = NULL, *target_name = NULL;
char *resultout_name = NULL;

void usage()
{
  cerr << "Usage: surfacedistance [target] [source] <options>\n" << endl;
  cerr << "Evaluate the surfacedistance between target and source\n" << endl;
  cerr << "<-RMS>				Use rms distance instead of mean distance" << endl;
  cerr << "<-output file>       Result output file" << endl;
  exit(1);
}

int main(int argc, char **argv)
{
  int ok, i, j, rms;
  double error,rerror,errorx,errory,errorz, dx, dy, dz;
  vtkPolyData *target,*source;
  vtkPolyDataReader *source_reader,*target_reader;
  double point1[3],point2[3];

  // Check command line
  if (argc < 3){
    usage();
  }

  rms = 0;

  // Parse source and target point lists
  target_name = argv[1];
  argc--;
  argv++;
  source_name = argv[1];
  argc--;
  argv++;

  // Read target
  target_reader = vtkPolyDataReader::New();
  target_reader->SetFileName(target_name);
  target_reader->Modified();
  target_reader->Update();
  target = vtkPolyData::New();
  target = target_reader->GetOutput();
  target->Update();
  // Read source
  source_reader = vtkPolyDataReader::New();
  source_reader->SetFileName(source_name);
  source_reader->Modified();
  source_reader->Update();
  source = vtkPolyData::New();
  source = source_reader->GetOutput();
  source->Update();

  // Parse remaining parameters
  while (argc > 1){
    ok = false;
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
    if (ok == false){
      cerr << "Can not parse argument " << argv[1] << endl;
      usage();
    }
  } 

  // Build a locator 
  vtkPointLocator *pointLocator = vtkPointLocator::New();
  pointLocator->SetDataSet(source_reader->GetOutput());
  pointLocator->BuildLocator();

  error = 0;
  for (i = 0; i < target->GetNumberOfPoints(); i++){
	  target->GetPoints()->GetPoint(i,point1);
	  vtkIdType ptId;
	  ptId = pointLocator->FindClosestPoint(point1);
	  source->GetPoints()->GetPoint(ptId,point2);
	  rerror = pow(point1[0] - point2[0], 2) 
		  + pow(point1[1] - point2[1], 2) +  pow(point1[2] - point2[2], 2);
	  if(!rms){
		  rerror = sqrt(rerror);
	  }
	  error += rerror;
  }
  
  // output to cout
  if(rms == 0){
      cout << "Mean distance is " << error/double(target->GetNumberOfPoints()) << " mm" << endl;
  }else{
	  cout << "RMS distance is " << sqrt(error/double(target->GetNumberOfPoints())) << " mm" << endl;
  }

  if(resultout_name){
	  cerr << "Writing Results: " << resultout_name << endl;

	  ofstream fout(resultout_name,ios::app);
	  if(rms == 0){
		  fout << error/double(target->GetNumberOfPoints()) << endl;
	  }else{
		  fout << sqrt(error/double(target->GetNumberOfPoints())) << endl;
	  }
	  fout.close();
  }

  //Final clean up
  pointLocator->Delete();
  target_reader->Delete();
  source_reader->Delete();
}

#else

#include <irtkImage.h>

int main( int argc, char *argv[] ){
  cerr << argv[0] << " needs to be compiled with the VTK library " << endl;
}
#endif

