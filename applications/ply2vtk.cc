/*=========================================================================

Library   : Image Registration Toolkit (IRTK)
Module    : $Id$
Copyright : Imperial College, Department of Computing
Visual Information Processing (VIP), 2009 onwards
Date      : $Date$
Version   : $Revision$
Changes   : $Author$

Copyright (c) 1999-2014 and onwards, Imperial College London
All rights reserved.
See LICENSE for details

=========================================================================*/

#include <irtkImage.h>
#include <irtkImageFunction.h>

#ifdef HAS_VTK

#include <vtkPLYReader.h>


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
#if VTK_MAJOR_VERSION < 6
	model->Update();
#endif


	vtkPolyDataWriter *writer = vtkPolyDataWriter::New();
#if VTK_MAJOR_VERSION >= 6
  	writer->SetInputData(model);
#else //VTK_MAJOR_VERSION >= 6
  	writer->SetInput(model);
#endif //VTK_MAJOR_VERSION >= 6

	writer->SetFileTypeToASCII();
	writer->SetFileName(output_name);
	writer->Write();

}

#else

int main( int argc, char *argv[] )
{
	cerr << argv[0] << " this program needs to be compiled with vtk enabled." << endl;
	return 0;
}

#endif
