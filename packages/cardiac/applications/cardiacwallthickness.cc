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

char *_target_name = NULL, *_source_name = NULL, *_myo_name = NULL;

void usage()
{
	cerr << "Usage: cardiacwallthickness [endo.vtk] [epi.vtk] \n" << endl;
	cerr << "-myocardium                 [myo.vtk]" << endl;
	exit(1);
}

int main(int argc, char **argv)
{
	int i, n, id;
	double rerror, source_point[3], target_point[3], myo_point[3];
	irtkLocator *source_locator;
	irtkLocator *target_locator;
	vtkDoubleArray *array_source = NULL;
	vtkDoubleArray *array_target = NULL;
	vtkDoubleArray *array_myo = NULL;

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

	while (argc > 1){
		int ok = false;
		if ((ok == false) && (strcmp(argv[1], "-myocardium") == 0)){
			argv++;
			argc--;
			ok = true;
			_myo_name = argv[1];
			argv++;
			argc--;
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

	// Myocardial pipeline
	vtkPolyData *myo = NULL;
	vtkPolyDataReader *myo_reader = NULL;
	if(_myo_name != NULL){
		cout << "Reading myocardium ... " << _myo_name << endl;
		myo_reader = vtkPolyDataReader::New();
		myo_reader->SetFileName(_myo_name);
		myo_reader->Modified();
		myo_reader->Update();
		myo = vtkPolyData::New();
		myo = myo_reader->GetOutput();
		myo->Update();
	}

	// Create source locator
	source_locator = new irtkLocator;
	source_locator->SelectLocatorType(1);
	source_locator->SetDataSet(source);

	// Create target locator
	target_locator = new irtkLocator;
	target_locator->SelectLocatorType(1);
	target_locator->SetDataSet(target);

	array_target = vtkDoubleArray::New();
	array_target->SetNumberOfTuples(target->GetNumberOfPoints());
	array_target->SetNumberOfComponents(1);
	array_target->SetName("WallThickness");

	array_source = vtkDoubleArray::New();
	array_source->SetNumberOfTuples(source->GetNumberOfPoints());
	array_source->SetNumberOfComponents(1);
	array_source->SetName("WallThickness");

	n = 0;
	for (i = 0; i < target->GetNumberOfPoints(); i++) {
		target->GetPoints()->GetPoint (i, target_point);
		source_point[0] = target_point[0];
		source_point[1] = target_point[1];
		source_point[2] = target_point[2];
		id = source_locator->FindClosestPoint (source_point);

		rerror = pow(target_point[0] -source_point[0], 2) 
			+ pow(target_point[1] - source_point[1], 2) 
			+  pow(target_point[2] - source_point[2], 2);
		rerror = sqrt(rerror);
		array_target->InsertTupleValue(i, &rerror);
		n++;
	}

	target->GetPointData()->SetScalars(array_target);

	vtkPolyDataWriter *writer = vtkPolyDataWriter::New();
	writer->SetInput(target);
	writer->SetFileName(_target_name);
	writer->Write();
	writer->Delete();

	n = 0;
	for (i = 0; i < source->GetNumberOfPoints(); i++) {
		source->GetPoints()->GetPoint (i, target_point);
		source_point[0] = target_point[0];
		source_point[1] = target_point[1];
		source_point[2] = target_point[2];
		id = target_locator->FindClosestPoint (source_point);

		rerror = pow(target_point[0] -source_point[0], 2) 
			+ pow(target_point[1] - source_point[1], 2) 
			+  pow(target_point[2] - source_point[2], 2);
		rerror = sqrt(rerror);
		array_source->InsertTupleValue(i, &rerror);
		n++;
	}

	source->GetPointData()->SetScalars(array_source);

	vtkPolyDataWriter *writersource = vtkPolyDataWriter::New();
	writersource->SetInput(source);
	writersource->SetFileName(_source_name);
	writersource->Write();
	writersource->Delete();

	//find points for myo
	if(myo != NULL){

		array_myo = vtkDoubleArray::New();
		array_myo->SetNumberOfTuples(myo->GetNumberOfPoints());
		array_myo->SetNumberOfComponents(1);
		array_myo->SetName("WallThickness");

		n = 0;
		for (i = 0; i < myo->GetNumberOfPoints(); i++) {
			myo->GetPoints()->GetPoint (i, target_point);
			source_point[0] = target_point[0];
			source_point[1] = target_point[1];
			source_point[2] = target_point[2];
			id = target_locator->FindClosestPoint (source_point);

			rerror = pow(target_point[0] -source_point[0], 2) 
				+ pow(target_point[1] - source_point[1], 2) 
				+  pow(target_point[2] - source_point[2], 2);
			rerror = sqrt(rerror);

			source_point[0] = target_point[0];
			source_point[1] = target_point[1];
			source_point[2] = target_point[2];
			int id2 = source_locator->FindClosestPoint (source_point);

			double rerror2;
			rerror2 = pow(target_point[0] -source_point[0], 2) 
				+ pow(target_point[1] - source_point[1], 2) 
				+  pow(target_point[2] - source_point[2], 2);
			rerror2 = sqrt(rerror2);

			if(rerror2 < rerror){
				rerror = *array_source->GetTuple(id2);
			}else{
				rerror = *array_target->GetTuple(id);
			}

			array_myo->InsertTupleValue(i, &rerror);
			n++;
		}

		myo->GetPointData()->SetScalars(array_myo);

		vtkPolyDataWriter *writermyo = vtkPolyDataWriter::New();
		writermyo->SetInput(myo);
		writermyo->SetFileName(_myo_name);
		writermyo->Write();
		writermyo->Delete();
	}

	array_target->Delete();
	array_source->Delete();
	if(array_myo != NULL)
		array_myo->Delete();

	delete source_locator;
	delete target_locator;
	target_reader->Delete();
	source_reader->Delete();
	if(myo_reader != NULL){
		myo_reader->Delete();
	}

}

#else

#include <irtkImage.h>

int main( int argc, char *argv[] )
{
	cerr << argv[0] << " needs to be compiled with the VTK library " << endl;
}
#endif
