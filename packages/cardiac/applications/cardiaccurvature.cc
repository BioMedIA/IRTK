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

#include <irtkImage.h>
#include <irtkVTKFunctions.h>
#include <vtkCurvatures.h>
#include <vtkDecimatePro.h>
#include <vtkSmoothPolyDataFilter.h>

char *input_name = NULL, *output_name = NULL;

void usage()
{
  cerr << "Usage: cardiaccurvature [input] [output] \n" << endl;
  cerr << "evaluate the curvature of endocardial surface \n" << endl;
  cerr << "  <-decimate>           Decimate the vertices in the output surface." << endl;
  cerr << "  <-smooth> iterations  Apply a number of iterations of (Laplacian) smoothing to " << endl;
  cerr << "                        the resulting surface." << endl;
  exit(1);
}

int main(int argc, char **argv)
{
    int iterations,ok;
    vtkDecimatePro *decimate = NULL;
    vtkSmoothPolyDataFilter *smooth = NULL;

    if (argc < 3) {
        usage();
    }
    input_name  = argv[1];
    argc--;
    argv++;
    output_name = argv[1];
    argc--;
    argv++;

    iterations = 0;

    while (argc > 1) {
        ok = false;
        if ((!ok) && (strcmp(argv[1], "-decimate") == 0)) {
            argc--;
            argv++;
            decimate = vtkDecimatePro::New();
            ok = true;
        }
        if ((!ok) && (strcmp(argv[1], "-smooth") == 0)) {
            argc--;
            argv++;
            smooth = vtkSmoothPolyDataFilter::New();
            iterations = atoi(argv[1]);
            argc--;
            argv++;
            ok = true;
        }
        if (!ok) {
            cerr << "Cannot parse argument " << argv[1] << endl;
            usage();
        }
    }

    vtkPolyDataReader *reader = vtkPolyDataReader::New();
    reader->SetFileName(input_name);
    reader->Modified();
    reader->Update();

    vtkCurvatures *curvature = vtkCurvatures::New();
    curvature->SetCurvatureTypeToMean();

    // Let's go to work
    if (decimate != NULL) {
        cout << "Decimating ... \n";
        decimate->SetInputConnection(reader->GetOutputPort());
        if (smooth != NULL) {
            cout << "Smoothing ... \n";
            smooth->SetInputConnection(decimate->GetOutputPort());
            smooth->SetNumberOfIterations(iterations);
            curvature->SetInputConnection(smooth->GetOutputPort());
        } else {
            curvature->SetInputConnection(decimate->GetOutputPort());
        }
    } else if (smooth != NULL) {
        cout << "Smoothing ... \n";
        smooth->SetNumberOfIterations(iterations);
        smooth->SetInputConnection(reader->GetOutputPort());
        curvature->SetInputConnection(smooth->GetOutputPort());
    } else {
        curvature->SetInputConnection(reader->GetOutputPort());
    }

    vtkPolyDataWriter *writer = vtkPolyDataWriter::New();
    writer->SetInputConnection(curvature->GetOutputPort());
    writer->SetFileName(output_name);
    writer->Write();
    writer->Delete();
    reader->Delete();
    curvature->Delete();
}

#else

#include <irtkImage.h>

int main( int argc, char *argv[] )
{
  cerr << argv[0] << " needs to be compiled with the VTK library " << endl;
}
#endif
