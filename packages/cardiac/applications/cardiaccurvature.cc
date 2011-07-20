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
#include <vtkTetra.h>
#include <vtkCellCenters.h>
#include <vtkCellLocator.h>
#include <vtkGenericCell.h>

char *input_name = NULL, *output_name = NULL, *mask_name = NULL;

void usage()
{
  cerr << "Usage: cardiaccurvature [input] [output] \n" << endl;
  cerr << "evaluate the curvature of endocardial surface \n" << endl;
  cerr << "  <-mask>   vtk value  Apply a mask to the curvature, if maskvalue < value do not calculate curvature " << endl;
  cerr << "  <-smooth> iterations  Apply a number of iterations of (Laplacian) smoothing to " << endl;
  cerr << "                        the resulting surface." << endl;
  exit(1);
}

int main(int argc, char **argv)
{
    int i,j,n,m,iterations,ok;
    double p1[3],p2[3],p3[3],p4[3],c[3],weight,distance,maskvalue,norm[3][3];
    vtkIdType k,k1,k2,k3;

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
    maskvalue = 0;

    while (argc > 1) {
        ok = false;
        if ((!ok) && (strcmp(argv[1], "-smooth") == 0)) {
            argc--;
            argv++;
            iterations = atoi(argv[1]);
            argc--;
            argv++;
            ok = true;
        }
        if ((!ok) && (strcmp(argv[1], "-mask") == 0)) {
            argc--;
            argv++;
            mask_name = argv[1];
            argc--;
            argv++;
            maskvalue = atof(argv[1]);
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
    reader->Update();

    vtkPolyData *output_model = vtkPolyData::New();

    vtkPolyDataReader *maskreader = vtkPolyDataReader::New();
    vtkPolyData *mask = NULL;
    if(mask_name != NULL){
        maskreader->SetFileName(mask_name);
        maskreader->Update();
        mask = maskreader->GetOutput();
    }     

    // Let's go to work
    cout << "evaluating the radial curvature ... \n";
    // Regulate
    output_model->DeepCopy(reader->GetOutput());

    // smooth
    double mean,max,min;

    vtkDecimatePro *decimate = vtkDecimatePro::New(); 
    decimate->SetInputConnection(reader->GetOutputPort());
    decimate->SetTargetReduction(0.5);
    decimate->Update();

    vtkPolyData *dmodel = decimate->GetOutput();

    vtkCellCenters *ccenter = vtkCellCenters::New(); 
    ccenter->SetInputConnection(decimate->GetOutputPort());
    ccenter->VertexCellsOn();
    ccenter->Update();

    vtkPolyData *cmodel = ccenter->GetOutput();

    // circumsphere replace angle curvature with radial curvature
    vtkTetra *tetra = vtkTetra::New();
    // move data to output
    vtkCellLocator *locator = vtkCellLocator::New(); 
    vtkDoubleArray *array = vtkDoubleArray::New();
    locator->SetDataSet(output_model); // data represents the surface 
    locator->BuildLocator(); 

    //mask
    vtkPointLocator *masklocator = vtkPointLocator::New(); 
    if(mask != NULL){
        masklocator->SetDataSet(mask);
        masklocator->BuildLocator();
    }

    mean = 8;
    for (i = 0; i < dmodel->GetNumberOfCells(); i++) {
        // from the center of the cell, stretch the cell and anchor on the surface!
        // anchored points and the center of the cell is used to generate a circumsphere
        // the radius of the circumsphere is the radial curvature!
        ok = 1;
        cmodel->GetPoints()->GetPoint(i,p1);
        if(mask != NULL){
            k = masklocator->FindClosestPoint(p1);
            if(*(mask->GetPointData()->GetScalars()->GetTuple(k)) < maskvalue){
                ok = 0;
            }
        }
        vtkIdList *list;
        //Find neighbor
        list = dmodel->GetCell(i)->GetPointIds();
        //Get number of neighbor
        n = list->GetNumberOfIds();
        if(n == 3){
            if(ok == 1){
                dmodel->GetPoints()->GetPoint(list->GetId(0),p2);
                dmodel->GetPoints()->GetPoint(list->GetId(1),p3);
                dmodel->GetPoints()->GetPoint(list->GetId(2),p4);
                //normalize p
                norm[0][0] = p2[0] - p1[0];
                norm[0][1] = p2[1] - p1[1];
                norm[0][2] = p2[2] - p1[2];
                norm[1][0] = p3[0] - p1[0];
                norm[1][1] = p3[1] - p1[1];
                norm[1][2] = p3[2] - p1[2];
                norm[2][0] = p4[0] - p1[0];
                norm[2][1] = p4[1] - p1[1];
                norm[2][2] = p4[2] - p1[2];
                for(j = 0; j < 3; j++){
                    distance = sqrt(pow(norm[j][0],2)+pow(norm[j][1],2)+pow(norm[j][2],2));
                    if(distance > 0){
                        for(k = 0; k < 3; k++){
                            norm[j][k] = norm[j][k] / distance;
                        }
                    }
                }
                //reevalute p with 5 mm distance
                p2[0] = p1[0] + norm[0][0] * mean;
                p2[1] = p1[1] + norm[0][1] * mean;
                p2[2] = p1[2] + norm[0][2] * mean;
                p3[0] = p1[0] + norm[1][0] * mean;
                p3[1] = p1[1] + norm[1][1] * mean;
                p3[2] = p1[2] + norm[1][2] * mean;
                p4[0] = p1[0] + norm[2][0] * mean;
                p4[1] = p1[1] + norm[2][1] * mean;
                p4[2] = p1[2] + norm[2][2] * mean;
                //find cloest point on surface
                c[0] = p2[0];
                c[1] = p2[1];
                c[2] = p2[2];
                locator->FindClosestPoint(c,p2,k1,j,distance);
                if(mask != NULL){
                    k = masklocator->FindClosestPoint(p2);
                    if(*(mask->GetPointData()->GetScalars()->GetTuple(k)) < maskvalue){
                        //anchor point on masked region do not use it.
                        p2[0] = c[0];
                        p2[1] = c[1];
                        p2[2] = c[2];
                        k1 = -1;
                    }
                }
                c[0] = p3[0];
                c[1] = p3[1];
                c[2] = p3[2];
                locator->FindClosestPoint(c,p3,k2,j,distance);
                if(mask != NULL){
                    k = masklocator->FindClosestPoint(p3);
                    if(*(mask->GetPointData()->GetScalars()->GetTuple(k)) < maskvalue){
                        //anchor point on masked region do not use it.
                        p3[0] = c[0];
                        p3[1] = c[1];
                        p3[2] = c[2];
                        k2 = -1;
                    }
                }
                c[0] = p4[0];
                c[1] = p4[1];
                c[2] = p4[2];
                locator->FindClosestPoint(c,p4,k3,j,distance);
                if(mask != NULL){
                    k = masklocator->FindClosestPoint(p4);
                    if(*(mask->GetPointData()->GetScalars()->GetTuple(k)) < maskvalue){
                        //anchor point on masked region do not use it.
                        p4[0] = c[0];
                        p4[1] = c[1];
                        p4[2] = c[2];
                        k3 = -1;
                    }
                }
                //evalute curvature
                if(k1 == k2 && k2 == k3){
                    weight = 0;
                }else{
                    weight = 1.0/tetra->Circumsphere(p1,p2,p3,p4,c);
                }
            }else{
                weight = 0;
            }
            array->InsertTupleValue(i,&weight);
        }else{
            cerr << "not triangle cell!" <<endl;
        }
    }
    cmodel->GetPointData()->SetScalars(array);
    tetra->Delete();

    max = 0.002;
    min = -0.002;
    for(i = 0; i < cmodel->GetNumberOfPoints(); i++){
        if(*(cmodel->GetPointData()->GetScalars()->GetTuple(i))>max)
            cmodel->GetPointData()->GetScalars()->SetTuple(i,&max);
        if(*(cmodel->GetPointData()->GetScalars()->GetTuple(i))<min)
            cmodel->GetPointData()->GetScalars()->SetTuple(i,&min);
    }
    //cmodel to dmodel
    vtkPointLocator *plocator = vtkPointLocator::New(); 
    vtkDoubleArray *narray;
    narray = vtkDoubleArray::New();
    plocator->SetDataSet(cmodel);
    plocator->BuildLocator();
    for (i = 0; i < dmodel->GetNumberOfPoints(); i++) {
        dmodel->GetPoints()->GetPoint(i,p1);
        k = plocator->FindClosestPoint(p1);
        double value;
        value = *(array->GetTuple(k));
        narray->InsertNextTupleValue(&value);
    }
    dmodel->GetPointData()->SetScalars(narray);
    double value, *tmp;
    // smooth 
    m = 8;
    while(m > 0){
        array = (vtkDoubleArray*)dmodel->GetPointData()->GetScalars();
        narray = NULL;
        narray = vtkDoubleArray::New();
        for (i = 0; i < dmodel->GetNumberOfPoints(); i++) {
            dmodel->GetPoints()->GetPoint(i,p1);
            vtkIdList *connectlist = vtkIdList::New();
            //Find neighbor
            GetConnectedVertices(dmodel,i,connectlist);
            //Get number of neighbor
            n = connectlist->GetNumberOfIds();
            weight = 1.0;
            value = (*(array->GetTuple(i)))*1.0;
            for (j = 0; j < n; j++){
                dmodel->GetPoints()->GetPoint(connectlist->GetId(j),p2);
                distance = vtkMath::Distance2BetweenPoints(p1,p2);
                tmp = array->GetTuple(connectlist->GetId(j));
                if(distance < 1)
                    distance = 1;
                value += (*tmp)/(distance);
                weight += 1.0/distance;
            }
            value = value / weight;
            narray->InsertTupleValue(i, &value);
            connectlist->Delete();
        }
        m --;
        array->DeepCopy(narray);
        narray->Delete();
    }

    //dmodel to outputmodel
    vtkPointLocator *pdlocator = vtkPointLocator::New(); 
    pdlocator->SetDataSet(dmodel);
    pdlocator->BuildLocator();

    narray = NULL;
    narray = vtkDoubleArray::New();
    for (i = 0; i < output_model->GetNumberOfPoints(); i++) {
        output_model->GetPoints()->GetPoint(i,p1);
        k = pdlocator->FindClosestPoint(p1);
        vtkIdList *dconnectlist = vtkIdList::New();
        //Find neighbor
        GetConnectedVertices(dmodel,k,dconnectlist);
        //Get number of neighbor
        n = dconnectlist->GetNumberOfIds();
        dmodel->GetPoints()->GetPoint(k,p2);
        distance = vtkMath::Distance2BetweenPoints(p1,p2);
        if(distance < 1) distance = 1;
        weight = 1.0/distance;
        value = (*(dmodel->GetPointData()->GetScalars()->GetTuple(k)))*weight;
        for (j = 0; j < n; j++){
            dmodel->GetPoints()->GetPoint(dconnectlist->GetId(j),p2);            
            distance = vtkMath::Distance2BetweenPoints(p1,p2);
            tmp = dmodel->GetPointData()->GetScalars()->GetTuple(dconnectlist->GetId(j));
            if(distance < 1) distance = 1;
            value += (*tmp)/(distance);
            weight += 1.0/distance;
        }
        value = value * 10000.0 / weight + 1.0;
        narray->InsertNextTupleValue(&value);
        dconnectlist->Delete();
    }
    output_model->GetPointData()->SetScalars(narray);

    //final smooth if necessary
    while(iterations > 0){
        array = (vtkDoubleArray*)output_model->GetPointData()->GetScalars();
        narray = NULL;
        narray = vtkDoubleArray::New();
        for (i = 0; i < output_model->GetNumberOfPoints(); i++) {
            output_model->GetPoints()->GetPoint(i,p1);
            vtkIdList *connectlist = vtkIdList::New();
            //Find neighbor
            GetConnectedVertices(output_model,i,connectlist);
            //Get number of neighbor
            n = connectlist->GetNumberOfIds();
            weight = 1.0;
            value = (*(array->GetTuple(i)))*1.0;
            for (j = 0; j < n; j++){
                output_model->GetPoints()->GetPoint(connectlist->GetId(j),p2);
                distance = vtkMath::Distance2BetweenPoints(p1,p2);
                tmp = array->GetTuple(connectlist->GetId(j));
                if(distance < 1)
                    distance = 1;
                value += (*tmp)/(distance);
                weight += 1.0/distance;
            }
            value = value / weight;
            narray->InsertTupleValue(i, &value);
            connectlist->Delete();
        }
        iterations --;
        array->DeepCopy(narray);
        narray->Delete();
    }

    cout << "... done! \n";

    vtkPolyDataWriter *writer = vtkPolyDataWriter::New();
    writer->SetInput(output_model);
    writer->SetFileName(output_name);
    writer->Write();
    writer->Delete();
    locator->Delete();
    plocator->Delete();
    pdlocator->Delete();
    masklocator->Delete();
    reader->Delete();
    ccenter->Delete();
    decimate->Delete();
    maskreader->Delete();
}

#else

#include <irtkImage.h>

int main( int argc, char *argv[] )
{
  cerr << argv[0] << " needs to be compiled with the VTK library " << endl;
}
#endif
