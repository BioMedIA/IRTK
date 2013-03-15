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
	cerr << "Usage: txt2cardiacvtk [input txt] [reference vtk] [output txt] " << endl;
	exit(1);
}

int main(int argc, char **argv)
{
	int i;
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
	if(model->GetPointData()->HasArray("DistanceProfile")){
		array = (vtkDoubleArray *)model->GetPointData()->GetArray("DistanceProfile");
	}else if(model->GetPointData()->HasArray("WallThickness")){
		array = (vtkDoubleArray *)model->GetPointData()->GetArray("WallThickness");
	}

	ifstream input;
	input.open(input_name);
	
	double scalar;

	for (i = 0; i < model->GetNumberOfPoints(); i++) {
		input >> point[0];
		input >> point[1];
		input >> point[2];

		model->GetPoints()->SetPoint (i, point);
		// if there is scalar output scalar
		if(array != NULL){
			input  >> scalar;

			array->SetTuple(i,&scalar);
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
