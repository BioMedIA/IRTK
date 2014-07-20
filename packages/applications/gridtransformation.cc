/*=========================================================================

  Library   : Image Registration Toolkit (IRTK)
  Module    : $Id$
  Copyright : Imperial College, Department of Computing
              Visual Information Processing (VIP), 2008 onwards
  Date      : $Date$
  Version   : $Revision$
  Changes   : $Author$

Copyright (c) IXICO LIMITED
All rights reserved.
See COPYRIGHT for details

=========================================================================*/

#include <irtkImage.h>

#include <irtkTransformation.h>

#ifdef HAS_VTK

#include <vtkUnstructuredGridWriter.h>
#include <vtkUnstructuredGridReader.h>
#include <vtkPoints.h>
#include <vtkUnstructuredGrid.h>

// Default filenames
char *input_name = NULL, *output_name = NULL;
char *dof_name  = NULL, *irtkoutput_name = NULL;

void usage()
{
  cerr << "Usage: gridtransformation [input] [output] <options>\n"
       << endl;
  cerr << "<-dofin file>                Transformation" << endl;
  cerr << "<-invert>                    Invert transformation" << endl;
  cerr << "<-source image>              reference source image" << endl;
  cerr << "<-target image>              reference target image" << endl;
  cerr << "<-partial number>            Transform 0-number points copy rest"<<endl;
  cerr << "<-time> [frame time]         Transformation is 4D use time"    << endl;
  cerr << "<-irtk filename>             output the result in irtk format"<<endl;
  exit(1);
}

int main(int argc, char **argv)
{
  int ok, invert,sourceon,targeton, pnumber;
  float time;
  irtkPointSet output, input;
  irtkTransformation *transformation;
  irtkGreyImage target,source;
  irtkPointSet irtkpointset;

  // Check command line
  if (argc < 3) {
    usage();
  }
  pnumber = 0;
  time = -1;


  // Parse input and output point lists
  input_name = argv[1];
  argc--;
  argv++;
  output_name = argv[1];
  argc--;
  argv++;

  // Default parameters
  invert = false;
  targeton = false;
  sourceon = false;

  // Parse arguments
  while (argc > 1) {
    ok = false;
    if ((ok == false) && (strcmp(argv[1], "-dofin") == 0)) {
      argc--;
      argv++;
      dof_name = argv[1];
      argc--;
      argv++;
      ok = true;
    }
    if ((ok == false) && (strcmp(argv[1], "-invert") == 0)) {
      argc--;
      argv++;
      invert = true;
      ok = true;
    }
	if ((ok == false) && (strcmp(argv[1], "-source") == 0)) {
      argc--;
      argv++;
      ok = true;
	  sourceon = true;
	  source.Read(argv[1]);
	  argc--;
      argv++;
    }
	if ((ok == false) && (strcmp(argv[1], "-target") == 0)) {
		argc--;
		argv++;
		targeton = true;
		target.Read(argv[1]);
		ok = true;
		argc--;
		argv++;
    }
	if ((ok == false) && (strcmp(argv[1], "-partial") == 0)) {
		argc--;
		argv++;
		pnumber = atoi(argv[1]);
		ok = true;
		argc--;
		argv++;
    }
    if ((ok == false) && (strcmp(argv[1], "-irtk") == 0)) {
        argc--;
        argv++;
        irtkoutput_name = argv[1];
        ok = true;
        argc--;
        argv++;
    }
	if ((ok == false) && (strcmp(argv[1], "-time") == 0)) {
		argc--;
		argv++;
		time = atof(argv[1]);
		argc--;
		argv++;
		cout << "evaluation for time: " << time << endl;
		ok = true;
	}
    if (ok == false) {
      cerr << "Can not parse argument " << argv[1] << endl;
      usage();
    }
  }

  if (dof_name != NULL) {
    // Read transformation
    transformation = irtkTransformation::New(dof_name);
  } else {
    // Create identity transformation
    transformation = new irtkRigidTransformation;
  }

  // Read point lists
  vtkUnstructuredGridReader *reader = vtkUnstructuredGridReader::New();
  reader->SetFileName(input_name);
  reader->Update();

  vtkUnstructuredGrid* surface = reader->GetOutput();
  vtkPoints*   points  = surface->GetPoints();

  for (int i=0; i < points->GetNumberOfPoints(); i++) {
    double p[3];
    points->GetPoint(i,p);
	if(sourceon&&targeton){
		target.WorldToImage(p[0],p[1],p[2]);
		//p[1] = target.GetY() - p[1];
		source.ImageToWorld(p[0],p[1],p[2]);
	}
	if(i<pnumber || pnumber == 0){
		if (invert == false) {
			if(time >= 0){
				//TFFD with time
				transformation->Transform(p[0],p[1],p[2],time);
			}else{
				transformation->Transform(p[0],p[1],p[2]);
			}
		} else {
			transformation->Inverse(p[0],p[1],p[2]);
		}
	}
    points->SetPoint(i,p);
    irtkpointset.Add(p);
  }

  // Write the final set
  vtkUnstructuredGridWriter   *writer = vtkUnstructuredGridWriter::New();
  writer->SetFileName(output_name);
#if VTK_MAJOR_VERSION >= 6
  writer->SetInputData(surface);
#else //VTK_MAJOR_VERSION >= 6
  writer->SetInput(surface);
#endif //VTK_MAJOR_VERSION >= 6

  writer->Write();

  if(irtkoutput_name){
      irtkpointset.WriteVTK(irtkoutput_name);
  }
}

#else
#include <irtkCommon.h>
int main(int argc, char **argv)
{
  cerr << "ptransformation: this program needs to be compiled with vtk enabled.\n";
  return 0;
}
#endif
