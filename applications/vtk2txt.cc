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

#ifdef HAS_VTK



char *input_name = NULL, *output_name = NULL;

void usage(){
	cerr << "Usage: vtk2txt [input] [output] " << endl;
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
	output_name = argv[1];
	argv++;
	argc--;

	// Read model
	vtkPolyDataReader *reader = vtkPolyDataReader::New();
	reader->SetFileName(input_name);
	reader->Modified();
	reader->Update();
	vtkPolyData *model = vtkPolyData::New();
	model = reader->GetOutput();
	model->Update();

    ofstream fout(output_name,ios::app);

    for (i = 0; i < model->GetNumberOfPoints(); i++) {
        model->GetPoints()->GetPoint (i, point);

        fout << point[0] << " " << point[1] << " " << point[2] << endl;
    }

    fout.close();
		
}

#else

#include <irtkImage.h>

int main( int argc, char *argv[] )
{
	cerr << argv[0] << " this program needs to be compiled with vtk enabled." << endl;
	return 0;
}

#endif
