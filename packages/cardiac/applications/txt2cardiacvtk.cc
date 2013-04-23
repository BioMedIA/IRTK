/*=========================================================================

Library   : Image Registration Toolkit (IRTK)
Module    : $Id: vtk2txt.cc 745 2012-12-14 12:17:33Z ws207 $
Copyright : Imperial College, Department of Computing
Visual Information Processing (VIP), 2009 onwards
Date      : $Date: 2012-12-14 12:17:33 +0000 (Fri, 14 Dec 2012) $
Version   : $Revision: 745 $
Changes   : $Author: ws207 $

=========================================================================*/

#include <irtkImage.h>
#include <irtkImageFunction.h>

#ifdef HAS_VTK



char *input_name = NULL, *reference_name = NULL, *output_name = NULL;

void usage(){
	cerr << "Usage: txt2cardiacvtk [input txt] [reference vtk] [output vtk] " << endl;
	cerr << "Options: -scalar if set presume the txt is 4D (3D + scalar)." << endl;
	cerr << "Options: -spatial if set use reference vtk's position." << endl;
	exit(1);
}

int main(int argc, char **argv)
{
	int i,ok,scalar = false,spatial = false;
	double point[3];

	if (argc < 3) {
		usage();
	}

	// Parse filenames
	input_name = argv[1];
	argv++;
	argc--;
	reference_name = argv[1];
	argv++;
	argc--;
	output_name = argv[1];
	argv++;
	argc--;

	// Parse remaining parameters
	while (argc > 1) {
		ok = false;
		if ((ok == false) && (strcmp(argv[1], "-scalar") == 0)){
			argc--;
			argv++;
			scalar = true;
			ok = true;
		}
		if ((ok == false) && (strcmp(argv[1], "-spatial") == 0)){
			argc--;
			argv++;
			spatial = true;
			ok = true;
		}
		if (ok == false) {
			cerr << "Can not parse argument " << argv[1] << endl;
			usage();
		}
	}

	// Read model
	vtkPolyDataReader *reader = vtkPolyDataReader::New();
	reader->SetFileName(reference_name);
	reader->Modified();
	reader->Update();
	vtkPolyData *model = vtkPolyData::New();
	model = reader->GetOutput();
	model->Update();

	// Read scalar if there's scalar
	vtkDoubleArray *array = NULL;
	if(scalar == true){
		if(model->GetPointData()->HasArray("DistanceProfile")){
			array = (vtkDoubleArray *)model->GetPointData()->GetArray("DistanceProfile");
		}else if(model->GetPointData()->HasArray("WallThickness")){
			array = (vtkDoubleArray *)model->GetPointData()->GetArray("WallThickness");
		}else {
			array = vtkDoubleArray::New();
			array->SetNumberOfTuples(model->GetNumberOfPoints());
			array->SetNumberOfComponents(1);
			array->SetName("WallThickness");
		}
	}

	ifstream input;
	input.open(input_name);

	double scalarv;

	for (i = 0; i < model->GetNumberOfPoints(); i++) {

		input >> point[0];
		input >> point[1];
		input >> point[2];

		if(spatial == false){
			model->GetPoints()->SetPoint (i, point);
		}
		// if there is scalar output scalar
		if(array != NULL){
			input  >> scalarv;

			array->SetTuple(i,&scalarv);
		}

	}

	if(array != NULL){
		model->GetPointData()->SetScalars(array);
	}

	vtkPolyDataWriter *writer = vtkPolyDataWriter::New();
	writer->SetInput(model);
	writer->SetFileName(output_name);
	writer->Write();
}

#else

#include <irtkImage.h>

int main( int argc, char *argv[] )
{
	cerr << argv[0] << " this program needs to be compiled with vtk enabled." << endl;
	return 0;
}

#endif
