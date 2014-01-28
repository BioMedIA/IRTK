/*=========================================================================

  Library   : Image Registration Toolkit (IRTK)
  Module    : $Id$
  Copyright : Imperial College, Department of Computing
              Visual Information Processing (VIP), 2008 onwards
  Date      : $Date$
  Version   : $Revision$
  Changes   : $Author$

Copyright (c) 1999-2014 and onwards, Imperial College London
All rights reserved.
See LICENSE for details

=========================================================================*/

// Image types and other standard things
#include <iostream>
#include <cstdlib> 
#include "vtkUnstructuredGrid.h"
#include <vtkUnstructuredGridReader.h>
#include <vtkUnstructuredGridWriter.h>
#include <vtkSmartPointer.h>
#include <vtkFloatArray.h>
#include <vtkPointData.h>
#include <vtkCellData.h>
#include <vtkIdList.h>
#include <irtkImage.h>
#include <irtkTransformation.h>

char* refrencename = NULL, *groundtruthname = NULL, *dofname = NULL, *outputmeshname = NULL;
char* longaxisname = NULL, *outputname = NULL, *ahameshname = NULL, *referenceimagename = NULL;

void usage()
{
    cerr << "Usage : miccai2012evaluation referencemesh.vtk groundTruth.vtk outputMesh.vtk ahamesh.vtk referenceimage.vtk" << std::endl;
    cerr << "-dofin [inputDof.dof.gz]           transformation between reference and groundtruth" << endl;
    cerr << "-strain                            evaluate the strain"    << endl;
    cerr << "-radialstrain [axis landmark name] evaluate the radial strain"    << endl;
    cerr << "-outputname   [output file name]   output the result to a txt file"    << endl;
    cerr << "-time         [frame time]         transformation is 4D use time"    << endl;
    cerr << "-invert                            invert the transformation" << endl;
    exit(1);
}

int main( int argc, char * argv[] )
{

    int ok,invert,strain;
    double error,x,y,z,time,aha;
    double p1[3],p2[3],radial[3];
    double myDisplacement[3];
    double trueDisplacement[3];
    double errorDisplacement[3];
    irtkGreyImage image;

    if (argc < 5)
    {
        usage();
        exit(EXIT_FAILURE);
    }
    const unsigned int MAX_PATH_LENGTH = 1024;
    time = -1;
    invert = false;
    strain = false;

    refrencename  = argv[1];
    argc--;
    argv++;
    groundtruthname = argv[1];
    argc--;
    argv++;
    outputmeshname = argv[1];
    argc--;
    argv++;
    ahameshname = argv[1];
    argc--;
    argv++;
    referenceimagename = argv[1];
    argc--;
    argv++;

    while (argc > 1) {
        ok = false;
        if ((ok == false) && (strcmp(argv[1], "-radialstrain") == 0)) {
            argc--;
            argv++;
            longaxisname = argv[1];
            argc--;
            argv++;
            ok = true;
        }
        if ((ok == false) && (strcmp(argv[1], "-strain") == 0)) {
            argc--;
            argv++;
            strain = true;
            ok = true;
        }
        if ((ok == false) && (strcmp(argv[1], "-dofin") == 0)) {
            argc--;
            argv++;
            dofname = argv[1];
            argc--;
            argv++;
            ok = true;
        }
        if ((ok == false) && (strcmp(argv[1], "-outputname") == 0)) {
            argc--;
            argv++;
            outputname = argv[1];
            argc--;
            argv++;
            ok = true;
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
        if ((ok == false) && (strcmp(argv[1], "-invert") == 0)) {
            argc--;
            argv++;
            invert = true;
            ok = true;
        }
        if (ok == false) {
            cerr << "Can not parse argument " << argv[1] << endl;
            usage();
        }
    }

    image.Read(referenceimagename);

    // reading mesh with regions
    vtkSmartPointer<vtkUnstructuredGridReader> readerAHA = vtkSmartPointer<vtkUnstructuredGridReader>::New();
    readerAHA->SetFileName(ahameshname);
    readerAHA->Update();

    //// Compute aha region for every point
    vtkSmartPointer<vtkIdList> idList = vtkSmartPointer<vtkIdList>::New();
    vtkSmartPointer<vtkFloatArray> ahaArray = vtkSmartPointer<vtkFloatArray>::New();
    ahaArray->SetName("aha");
    ahaArray->SetNumberOfComponents(1);

    for (int p=0; p<readerAHA->GetOutput()->GetPoints()->GetNumberOfPoints(); p++)
    {
        // Get aha region associated to this point.
        // Since aha regions are specified per cell, we first query for the cells and thake the first one
        idList->Reset();
        readerAHA->GetOutput()->GetPointCells(p, idList);
        float ahaRegion = (float)readerAHA->GetOutput()->GetCellData()->GetArray("Zones")->GetTuple1(idList->GetId(0));
        ahaArray->InsertNextTuple1(ahaRegion);
    }

    // Loading the first mesh of the list (this assumes you have computed displacements to the first frame)
    vtkSmartPointer<vtkUnstructuredGridReader> readerFirst = vtkSmartPointer<vtkUnstructuredGridReader>::New();
    readerFirst->SetFileName(refrencename);
    readerFirst->Update();
    vtkPoints* pointsFirst = readerFirst->GetOutput()->GetPoints();

    // Read unstructured grid in the file list
    vtkSmartPointer<vtkUnstructuredGridReader> reader = vtkSmartPointer<vtkUnstructuredGridReader>::New();
    reader->SetFileName(groundtruthname);
    reader->Update();
    vtkPoints* points = reader->GetOutput()->GetPoints();

    // Output mesh allocation (including array for storing error)
    vtkSmartPointer<vtkUnstructuredGrid> ug = vtkSmartPointer<vtkUnstructuredGrid>::New();
    ug->DeepCopy(reader->GetOutput());
    vtkSmartPointer<vtkFloatArray> errorArray = vtkSmartPointer<vtkFloatArray>::New();
    errorArray->SetName("error");
    errorArray->SetNumberOfComponents(3);
    // Set up strain value
    vtkSmartPointer<vtkFloatArray> scalarvectors = vtkSmartPointer<vtkFloatArray>::New();
    scalarvectors->SetNumberOfComponents(1);
    //ug->GetPointData()->AddArray(ahaArray);

    // Read dof
    irtkTransformation *transformation;
    if (dofname != NULL) {
        // Read transformation
        transformation = irtkTransformation::New(dofname);
    } else {
        // Create identity transformation
        transformation = new irtkRigidTransformation;
    }

    // Read long axis landmark
    irtkPointSet landmarks;
    if(longaxisname){
        landmarks.ReadVTK(longaxisname);
        if(landmarks.Size() == 2){
            p1[0] = landmarks(0)._x;
            p1[1] = landmarks(0)._y;
            p1[2] = landmarks(0)._z;

            // convert p2 to normal
            p2[0] = landmarks(1)._x - landmarks(0)._x;
            p2[1] = landmarks(1)._y - landmarks(0)._y;
            p2[2] = landmarks(1)._z - landmarks(0)._z;
            double distance = sqrt(pow(p2[0],2)+pow(p2[1],2)+pow(p2[2],2));
            if(distance > 0){
                for(int j = 0; j < 3; j++){
                    p2[j] = p2[j] / distance;
                }
            }
        }else{
            cerr << "Longitudinal axis needs to be defined by 2 points, apex and mid valve" << endl;
            exit(1);
        }
    }

    ofstream fout;

    if(outputname){
        fout.open(outputname,ios::app);
    }

    // Compute error for each point of the mesh
    for (int p=0; p<points->GetNumberOfPoints(); p++)
    {

        error = 0;
        myDisplacement[0] = 0;
        myDisplacement[1] = 0;
        myDisplacement[2] = 0;

        x = pointsFirst->GetPoint(p)[0];
        y = pointsFirst->GetPoint(p)[1];
        z = pointsFirst->GetPoint(p)[2];

        image.WorldToImage(x,y,z);

        if(image.GetAsDouble(x,y,z) > 1000){

			if(invert == true){
				x = points->GetPoint(p)[0];
				y = points->GetPoint(p)[1];
				z = points->GetPoint(p)[2];
			}else{
				x = pointsFirst->GetPoint(p)[0];
				y = pointsFirst->GetPoint(p)[1];
				z = pointsFirst->GetPoint(p)[2];
			}

			if(time >= 0){
				//TFFD with time
				transformation->LocalDisplacement(x,y,z,time);
			}else{
				//FFD without time  
				transformation->LocalDisplacement(x,y,z);
			}
			myDisplacement[0] = x;
			myDisplacement[1] = y;
			myDisplacement[2] = z;

			for (int d=0; d<3; d++)
			{
				// Groundtruth displacement
				if(invert == true){
					trueDisplacement[d] = pointsFirst->GetPoint(p)[d] - points->GetPoint(p)[d];
				}else{
					trueDisplacement[d] = points->GetPoint(p)[d] - pointsFirst->GetPoint(p)[d];
				}

				// Error made on displacement
				errorDisplacement[d] = myDisplacement[d] - trueDisplacement[d];

				//calculate error
				error += errorDisplacement[d]*errorDisplacement[d];
			}

			error = sqrt(error);
		}

		if(outputname){
			fout << error << " ";
		}

        aha = *(ahaArray->GetTuple(p));

        if(aha < 1 || aha > 17){
            error = 1000;
            myDisplacement[0] = 1000;
            myDisplacement[1] = 1000;
            myDisplacement[2] = 1000;
        }

        // Insert error in the array
        errorArray->InsertNextTuple(errorDisplacement);
        scalarvectors->InsertNextTuple(&error);
    }

    if(outputname){;
        fout.close();
    }

    ug->GetPointData()->SetVectors(errorArray);
    ug->GetPointData()->SetScalars(scalarvectors);

    // Write output mesh
    vtkSmartPointer<vtkUnstructuredGridWriter> writer =vtkSmartPointer<vtkUnstructuredGridWriter>::New();
    writer->SetFileName(outputmeshname);
    writer->SetInput(ug);
    writer->Update();

    if(outputname){
        ofstream fout(outputname,ios::app);
        fout << endl;
        fout.close();
    }

    exit(EXIT_SUCCESS);

