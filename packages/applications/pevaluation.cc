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

char *_target_name = NULL, *_source_name = NULL;
//char *dofin_name = NULL, *dofout_name = NULL;
char *resultout_name = NULL;

void usage()
{
  cerr << "Usage: pevaluation [target] [source]" << endl;
  cerr << "<-output file>       Result output file" << endl;
  exit(1);
}

int main(int argc, char **argv)
{
  int i, ok;
  double error, source_point[3], target_point[3];

  if (argc < 3) {
    usage();
  }

  // Parse filenames
  _target_name = argv[1];
  argv++;
  argc--;
  _source_name = argv[1];
  argv++;
  argc--;

  // Parse remaining parameters
  while (argc > 1){
    ok = false;
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

  if (target->GetNumberOfPoints() != source->GetNumberOfPoints()) {
    cerr << "Landmarks should have equal number of points" << endl;
    exit(1);
  }

  error = 0;
  for (i = 0; i < target->GetNumberOfPoints(); i++) {
    target->GetPoints()->GetPoint (i, target_point);
    source->GetPoints()->GetPoint (i, source_point);
    error += sqrt((target_point[0] - source_point[0]) * (target_point[0] - source_point[0]) +
                  (target_point[1] - source_point[1]) * (target_point[1] - source_point[1]) +
                  (target_point[2] - source_point[2]) * (target_point[2] - source_point[2]));
    cout << "Error for landmark " << i + 1 << " = " <<
         sqrt((target_point[0] - source_point[0]) * (target_point[0] - source_point[0]) +
              (target_point[1] - source_point[1]) * (target_point[1] - source_point[1]) +
              (target_point[2] - source_point[2]) * (target_point[2] - source_point[2])) << endl;

  }

  error /= (double)target->GetNumberOfPoints();
  cout << "RMS = " << error << " mm" << endl;

  double tstd = 0, stderror = 0;
  for (i = 0; i < target->GetNumberOfPoints(); i++){
	  target->GetPoints()->GetPoint (i, target_point);
	  source->GetPoints()->GetPoint (i, source_point);
	  tstd = sqrt((target_point[0] - source_point[0]) * (target_point[0] - source_point[0]) +
		  (target_point[1] - source_point[1]) * (target_point[1] - source_point[1]) +
		  (target_point[2] - source_point[2]) * (target_point[2] - source_point[2]));
	  tstd = pow(tstd - error,2);
	  stderror += tstd;
  }
  stderror=sqrt(stderror/double(target->GetNumberOfPoints()));
  cout<<"standard deviation is " << stderror  << " mm" << endl;

  if(resultout_name){
	  cerr << "Writing Results: " << resultout_name << endl;
	  ofstream fout(resultout_name,ios::app);	  
	  fout << error << " " << stderror << endl;
	  fout.close();
  }
}

#else

#include <irtkImage.h>

int main( int argc, char *argv[] )
{
  cerr << argv[0] << " needs to be compiled with the VTK library " << endl;
}
#endif
