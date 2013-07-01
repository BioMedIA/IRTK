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
#include <vtkPLYReader.h>

#ifdef HAS_VTK



char *input_name = NULL, *output_name = NULL;

void usage(){
	cerr << "Usage: ply2vtk [input] [output] " << endl;
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
	vtkPLYReader *reader = vtkPLYReader::New();
	reader->SetFileName(input_name);
	reader->Modified();
	reader->Update();
	vtkPolyData *model = vtkPolyData::New();
	model = reader->GetOutput();
	model->Update();

	vtkPolyDataWriter *writer = vtkPolyDataWriter::New();
	writer->SetInput(model);
	writer->SetFileTypeToASCII();
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
