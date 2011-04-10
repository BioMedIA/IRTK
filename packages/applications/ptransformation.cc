/*=========================================================================

  Library   : Image Registration Toolkit (IRTK)
  Module    : $Id$
  Copyright : Imperial College, Department of Computing
              Visual Information Processing (VIP), 2008 onwards
  Date      : $Date$
  Version   : $Revision$
  Changes   : $Author$

=========================================================================*/

#include <irtkImage.h>

#include <irtkTransformation.h>

#ifdef HAS_VTK

#include <vtkPolyDataWriter.h>
#include <vtkPolyDataReader.h>
#include <vtkPoints.h>
#include <vtkPolyData.h>

// Default filenames
char *input_name = NULL, *output_name = NULL;
char *dof_name  = NULL;

void usage()
{
  cerr << "Usage: ptransformation [input] [output] <options>\n"
       << endl;
  cerr << "<-dofin file>      Transformation" << endl;
  cerr << "<-invert>          Invert transformation" << endl;
  cerr << "<-source image>    reference source image" << endl;
  cerr << "<-target image>    reference target image" << endl;
  cerr << "<-partial number>  Transform 0-number points copy rest"<<endl;
  exit(1);
}

int main(int argc, char **argv)
{
  int ok, invert,sourceon,targeton, pnumber;
  irtkPointSet output, input;
  irtkTransformation *transformation;
  irtkGreyImage target,source;

  // Check command line
  if (argc < 3) {
    usage();
  }
  pnumber = 0;


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
  vtkPolyDataReader *reader = vtkPolyDataReader::New();
  reader->SetFileName(input_name);
  reader->Update();

  vtkPolyData* surface = reader->GetOutput();
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
			transformation->Transform(p[0],p[1],p[2]);
		} else {
			transformation->Inverse(p[0],p[1],p[2]);
		}
	}
    points->SetPoint(i,p);
  }

  // Write the final set
  vtkPolyDataWriter   *writer = vtkPolyDataWriter::New();
  writer->SetFileName(output_name);
  writer->SetInput(surface);
  writer->Write();
}

#else

int main(int argc, char **argv)
{
  cerr << "ptransformation: this program needs to be compiled with vtk enabled.\n";
  return 0;
}
#endif
