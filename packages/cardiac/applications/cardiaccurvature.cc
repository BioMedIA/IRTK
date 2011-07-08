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
#include <vtkMath.h>

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
    int i,j,n,iterations,ok;
    double mean,std;
    vtkDecimatePro *decimate = NULL;

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
    curvature->SetCurvatureTypeToGaussian();

    // Let's go to work
    if (decimate != NULL) {
        cout << "Decimating ... \n";
        decimate->SetInputConnection(reader->GetOutputPort());
        curvature->SetInputConnection(decimate->GetOutputPort());      
    } else {
        curvature->SetInputConnection(reader->GetOutputPort());
    }

    curvature->Update();
    // Regulate
    vtkPolyData *model = curvature->GetOutput();

    mean = 0; std = 0;
    for(i = 0; i < model->GetNumberOfPoints(); i++){
        mean += *(model->GetPointData()->GetScalars()->GetTuple(i));
    }
    mean = mean/model->GetNumberOfPoints();
    for(i = 0; i < model->GetNumberOfPoints(); i++){
        std += pow((*(model->GetPointData()->GetScalars()->GetTuple(i)) - mean),2);
    }
    std = sqrt(std/model->GetNumberOfPoints());
    double max,min;
    //max = mean+2*std;
    //min = mean-2*std;
    max = 0.01;
    min = -0.01;
    for(i = 0; i < model->GetNumberOfPoints(); i++){
        if(*(model->GetPointData()->GetScalars()->GetTuple(i))>max)
            model->GetPointData()->GetScalars()->SetTuple(i,&max);
        if(*(model->GetPointData()->GetScalars()->GetTuple(i))<min)
            model->GetPointData()->GetScalars()->SetTuple(i,&min);

    }

    // smooth
    if(iterations > 0){
        double p1[3],p2[3],weight;
        vtkDoubleArray *array = (vtkDoubleArray*)model->GetPointData()->GetScalars();
        vtkDoubleArray *narray = vtkDoubleArray::New();
        while(iterations > 0){
            for (i = 0; i < model->GetNumberOfPoints(); i++) {
                model->GetPoints()->GetPoint(i,p1);
                vtkIdList *list = vtkIdList::New();
                //Find neighbor
                GetConnectedVertices(model,i,list);
                //Get number of neighbor
                n = list->GetNumberOfIds();
                double value, *tmp;
                weight = 4.0;
                value = (*(array->GetTuple(i)))*4.0;
                for (j = 0; j < n; j++){
                    model->GetPoints()->GetPoint(list->GetId(j),p2);
                    double distance;
                    distance = vtkMath::Distance2BetweenPoints(p1,p2);
                    tmp = array->GetTuple(list->GetId(j));
                    if(distance > 0){
                        if(distance < 0.25)
                            distance = 0.25;
                        value += (*tmp)/(distance);
                        weight += 1.0/distance;
                    }
                }
                value = value / weight;
                narray->InsertTupleValue(i, &value);
                list->Delete();
            }
            iterations --;
        }
        array->DeepCopy(narray);
        narray->Delete();
    }

    vtkPolyDataWriter *writer = vtkPolyDataWriter::New();
    writer->SetInput(model);
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
