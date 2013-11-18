/*=========================================================================

Library   : Image Registration Toolkit (IRTK)
Module    : $Id$
Copyright : Imperial College, Department of Computing
Visual Information Processing (VIP), 2009 onwards
Date      : $Date$
Version   : $Revision$
Changes   : $Author$

=========================================================================*/

#include <irtkImage.h>
#include <irtkImageFunction.h>
#include <irtkVTKFunctions.h>

#ifdef HAS_VTK

char *input_name = NULL, *output_name = NULL, *array_name = "WallThickness", *connection_name = NULL;

void usage(){
	cerr << "Usage: vtk2txt [input] [output] " << endl;
	cerr << "Options:" << endl;
	cerr << "  -arrayname <name>      Array name of the input ply format." << endl;
	cerr << "  -connection <output>   Output connection to a txt file." << endl;
	exit(1);
}

int main(int argc, char **argv)
{
	int i, ok;
	double point[3];

	if (argc < 3) {
		usage();
	}

	// Parse filenames
	input_name = argv[1];
	argv++;
	argc--;
	output_name = argv[1];
	argv++;
	argc--;

	while (argc > 1) {
		ok = false;
		if ((ok == false) && (strcmp(argv[1], "-arrayname") == 0)) {
			argc--;
			argv++;
			ok = true;
			array_name = argv[1];
			argc--;
			argv++;
		}
		if ((ok == false) && (strcmp(argv[1], "-connection") == 0)) {
			argc--;
			argv++;
			ok = true;
			connection_name = argv[1];
			argc--;
			argv++;
		}
		if (ok == false) {
			cerr << "Can not parse argument " << argv[1] << endl;
			usage();
		}
	}

	// Read model
	vtkPolyDataReader *reader = vtkPolyDataReader::New();
	reader->SetFileName(input_name);
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
	}else if(model->GetPointData()->HasArray(array_name)){
		array = (vtkDoubleArray *)model->GetPointData()->GetArray(array_name);
	}

	ofstream fout;
	fout.open(output_name);

	for (i = 0; i < model->GetNumberOfPoints(); i++) {
		model->GetPoints()->GetPoint (i, point);

		fout << point[0] << " " << point[1] << " " << point[2];
		// if there is scalar output scalar
		if(array != NULL){
			fout  << " " << *array->GetTuple(i);
		}
		fout << endl;
	}

	fout.close();
	fout.clear();

	if(connection_name != NULL){
		
		fout.open(connection_name);

		for (i = 0; i < model->GetNumberOfPoints(); i++) {
			fout << i;

			vtkIdList *dconnectlist = vtkIdList::New();
			//Find neighbor
			GetConnectedVertices(model,i,dconnectlist);
			//Get number of neighbor
			int n = dconnectlist->GetNumberOfIds();
			
			for (int j = 0; j < n; j++){
				fout << " " << dconnectlist->GetId(j);
			}
			dconnectlist->Delete();
			
			fout << endl;
		}

		fout.close();
	}

}

#else

#include <irtkImage.h>

int main( int argc, char *argv[] )
{
	cerr << argv[0] << " this program needs to be compiled with vtk enabled." << endl;
	return 0;
}

#endif
