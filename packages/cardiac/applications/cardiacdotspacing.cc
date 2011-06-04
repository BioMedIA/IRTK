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
#include <irtkRegistration.h>

#ifdef HAS_VTK

#include <vtkPointData.h>
#include <vtkPolyData.h>
#include <vtkPolyDataReader.h>
#include <vtkPolyDataWriter.h>
#include <vtkDoubleArray.h>
#include <vtkDecimatePro.h>

#endif

char *input_name = NULL, *output_name = NULL, *dofin_name = NULL, *out_name = NULL;

/// example: cardiacdotspacing 1 surface\endo1.vtk surface\streching.vtk -dofin transformation\lnreg.dof -outsurface surface\epi1.vtk
void usage()
{
	cerr << "Usage: cardiacdotspacing [mode] [input pointset/mesh] [output txt/mesh]\n" << endl;
	cerr << "mode 0 input pointset output txt mode 1 input mesh output mesh needs dofin\n" << endl;
	cerr << "pointset must have a size of 16 representing 16 segmentation sectors\n" << endl;
	cerr << "-dofin [transformation]       mesh dotspacing option, input transformation\n" << endl;
	cerr << "-outsurface [name]       mesh dotspacing option, input transformation\n" << endl;
	exit(1);
}

int main(int argc, char **argv)
{
	if( argc < 4 ) usage();

	int mode = atoi(argv[1]);
	argv++;
	argc--;

	if(mode == 0){

		int i;
		irtkPointSet input;
		char* output;
		double bep[17],dx,dy,dz,d;

		input.ReadVTK(argv[1]);
		if(input.Size() != 16) usage();
		argv++;
		argc--;
		output = argv[1];
		argv++;
		argc--;

		bep[16] = 0;

		for(i=0;i<16;i++){
			if(i==0||i==6){
				dx = input(i)._x - input(i+1)._x;
				dy = input(i)._y - input(i+1)._y;
				dz = input(i)._z - input(i+1)._z;
				d = sqrt(pow(dx,2)+pow(dy,2)+pow(dz,2));
				dx = input(i)._x - input(i+5)._x;
				dy = input(i)._y - input(i+5)._y;
				dz = input(i)._z - input(i+5)._z;
				d += sqrt(pow(dx,2)+pow(dy,2)+pow(dz,2));
				bep[i] = d/2;
			}else if(i==5||i==11){
				dx = input(i)._x - input(i-5)._x;
				dy = input(i)._y - input(i-5)._y;
				dz = input(i)._z - input(i-5)._z;
				d = sqrt(pow(dx,2)+pow(dy,2)+pow(dz,2));
				dx = input(i)._x - input(i-1)._x;
				dy = input(i)._y - input(i-1)._y;
				dz = input(i)._z - input(i-1)._z;
				d += sqrt(pow(dx,2)+pow(dy,2)+pow(dz,2));
				bep[i] = d/2;
			}else if(i==12){
				dx = input(i)._x - input(i+1)._x;
				dy = input(i)._y - input(i+1)._y;
				dz = input(i)._z - input(i+1)._z;
				d = sqrt(pow(dx,2)+pow(dy,2)+pow(dz,2));
				dx = input(i)._x - input(i+3)._x;
				dy = input(i)._y - input(i+3)._y;
				dz = input(i)._z - input(i+3)._z;
				d += sqrt(pow(dx,2)+pow(dy,2)+pow(dz,2));
				bep[i] = d/2;
			}else if(i==15){
				dx = input(i)._x - input(i-3)._x;
				dy = input(i)._y - input(i-3)._y;
				dz = input(i)._z - input(i-3)._z;
				d = sqrt(pow(dx,2)+pow(dy,2)+pow(dz,2));
				dx = input(i)._x - input(i-1)._x;
				dy = input(i)._y - input(i-1)._y;
				dz = input(i)._z - input(i-1)._z;
				d += sqrt(pow(dx,2)+pow(dy,2)+pow(dz,2));
				bep[i] = d/2;
			}else{
				dx = input(i)._x - input(i+1)._x;
				dy = input(i)._y - input(i+1)._y;
				dz = input(i)._z - input(i+1)._z;
				d = sqrt(pow(dx,2)+pow(dy,2)+pow(dz,2));
				dx = input(i)._x - input(i-1)._x;
				dy = input(i)._y - input(i-1)._y;
				dz = input(i)._z - input(i-1)._z;
				d += sqrt(pow(dx,2)+pow(dy,2)+pow(dz,2));
				bep[i] = d/2;
			}
		}
		remove(output);

		ofstream fout(output,ios::app);
		for(i=0; i<17; i++)
			fout << bep[i] <<" ";
		fout << endl;
		fout.close();
	}else if(mode == 1){
#ifdef HAS_VTK
		int i, j, ok;
		double n, ds[2], point1[3], point2[3],tpoint1[3],tpoint2[3];
		double *stretch = new double[1];

		// Parse filenames
		input_name = argv[1];
		argv++;
		argc--;
		output_name = argv[1];
		argv++;
		argc--;

		while (argc > 1) {
			ok = false;
			if ((ok == false) && (strcmp(argv[1], "-dofin") == 0)) {
				argc--;
				argv++;
				dofin_name = argv[1];
				argc--;
				argv++;
				ok = true;
			}
			if ((ok == false) && (strcmp(argv[1], "-outsurface") == 0)) {
				argc--;
				argv++;
				out_name = argv[1];
				argc--;
				argv++;
				ok = true;
			}
			if (ok == false) {
				cerr << "Can't parse argument " << argv[1] << endl;
				usage();
			}
		}

		// Create initial multi-level free-form deformation
		irtkTransformation *transform = NULL;
		if (dofin_name != NULL){
			transform = irtkTransformation::New(dofin_name);
		} else {
			cerr << "mode 1 must have transformation" << endl;
            exit(1);
		}

		// Read model
		vtkPolyDataReader *reader = vtkPolyDataReader::New();
		reader->SetFileName(input_name);
		reader->Modified();
		reader->Update();
		vtkPolyData *model = vtkPolyData::New();
		vtkPolyData *valuemodel = vtkPolyData::New();
		vtkPolyData *decimatemodel = vtkPolyData::New();
		
		model = reader->GetOutput();
		model->Update();
		valuemodel->CopyStructure(model);
		valuemodel->Update();

		if(out_name != NULL){
			// load out surface
			vtkPolyDataReader *reader_out = vtkPolyDataReader::New();
			reader_out->SetFileName(out_name);
			reader_out->Modified();
			reader_out->Update();
			vtkPolyData *model_out = vtkPolyData::New();

			model_out = reader->GetOutput();
			model_out->Update();

			//valuemodel = (model + model_out)/2
			for (i = 0; i < model->GetNumberOfPoints(); i++) {
				int closest,j;
				double distance = 1000000;
				double point[3],point2[3],normal[3];
				model->GetPoints()->GetPoint (i, point);
                closest = i;
				for(j = 0; j < model_out->GetNumberOfPoints(); j++){
					model_out->GetPoints()->GetPoint (j, point2);
					if(distance > sqrt(pow(point[0] - point2[0],2)+
						pow(point[1] - point2[1],2)+pow(point[2] - point2[2],2))){
							distance = sqrt(pow(point[0] - point2[0],2)+
								pow(point[1] - point2[1],2)+pow(point[2] - point2[2],2));
							closest = j;
					}
				}
				model_out->GetPoints()->GetPoint (closest, point2);
				normal[0] = point2[0] + point[0];
				normal[1] = point2[1] + point[1];
				normal[2] = point2[2] + point[2];
				for(j = 0; j < 3; j++){
					normal[j] = normal[j] / 2;
				}
				valuemodel->GetPoints()->SetPoint(i,normal);
			}
				reader_out->Delete();			
		}

		// decimate model
		vtkDecimatePro *decimate = vtkDecimatePro::New();
		cout << "Decimating ... \n";
		decimate->SetInput(valuemodel);
		decimatemodel = decimate->GetOutput();
		decimatemodel->Update();

		// Allocate memory
		cout << "Evaluating Point Spacing ... \n";
		vtkDoubleArray *decimatearray = vtkDoubleArray::New();
		decimatearray->SetNumberOfTuples(decimatemodel->GetNumberOfPoints());
		decimatearray->SetNumberOfComponents(1);
		decimatearray->SetName("stretch precentage");
		vtkDoubleArray *ndecimatearray = vtkDoubleArray::New();
		ndecimatearray->SetNumberOfTuples(decimatemodel->GetNumberOfPoints());
		ndecimatearray->SetNumberOfComponents(1);
		ndecimatearray->SetName("stretch precentage");

		for (i = 0; i < decimatemodel->GetNumberOfPoints(); i++) {
			vtkIdList *list = vtkIdList::New();
			//Find neighbor
			GetConnectedVertices(decimatemodel,i,list);
			//Get number of neighbor
			n = list->GetNumberOfIds();
			//Get current vertices
			decimatemodel->GetPoints()->GetPoint (i, point1);
			//transform current vertices
			tpoint1[0] = point1[0]; tpoint1[1] = point1[1]; tpoint1[2] = point1[2];
			transform->Transform(tpoint1[0],tpoint1[1],tpoint1[2]);
			//Initialize ds
			ds[0] = 0; ds[1] = 0;
			for (j = 0; j < n; j++){
			decimatemodel->GetPoints()->GetPoint (list->GetId(j), point2);
			ds[0] += sqrt(pow(point2[2] - point1[2], 2)
				+ pow(point2[1] - point1[1], 2) + pow(point2[0] - point1[0], 2))/n;
			//transform neighbor vertices
			tpoint2[0] = point2[0]; tpoint2[1] = point2[1]; tpoint2[2] = point2[2];
			transform->Transform(tpoint2[0],tpoint2[1],tpoint2[2]);
			//evaluate new distance
			ds[1] += sqrt(pow(tpoint2[2] - tpoint1[2], 2)
				+ pow(tpoint2[1] - tpoint1[1], 2) + pow(tpoint2[0] - tpoint1[0], 2))/n;
			}
			//stretch measure relative change of the vertices spacing
			*stretch = (ds[1] - ds[0]) / ds[0];
			decimatearray->InsertTupleValue(i, stretch);
			list->Delete();
		}

		// normalize		
		for (i = 0; i < decimatemodel->GetNumberOfPoints(); i++) {
			vtkIdList *list = vtkIdList::New();
			//Find neighbor
			GetConnectedVertices(decimatemodel,i,list);
			//Get number of neighbor
			n = list->GetNumberOfIds();
			double value = 0, *tmp;
			for (j = 0; j < n; j++){
				tmp = decimatearray->GetTuple(list->GetId(j));
				value += *tmp;
			}
			//stretch measure relative change of the vertices spacing
			*stretch = value / n;
			ndecimatearray->InsertTupleValue(i, stretch);
			list->Delete();
		}

		vtkDoubleArray *array = vtkDoubleArray::New();
		array->SetNumberOfTuples(decimatemodel->GetNumberOfPoints());
		array->SetNumberOfComponents(1);
		array->SetName("stretch precentage");

		for (i = 0; i < model->GetNumberOfPoints(); i++) {
			vtkIdList *list = vtkIdList::New();
			//Get current vertices
			valuemodel->GetPoints()->GetPoint (i, point1);
			//Find closest decimate point
			int closest,j;
			double distance = 1000000;
            closest = i;
			for(j = 0; j < decimatemodel->GetNumberOfPoints(); j++){
				decimatemodel->GetPoints()->GetPoint (j, point2);
				if(distance > sqrt(pow(point1[0] - point2[0],2)+
					pow(point1[1] - point2[1],2)+pow(point1[2] - point2[2],2))){
					distance = sqrt(pow(point1[0] - point2[0],2)+
					pow(point1[1] - point2[1],2)+pow(point1[2] - point2[2],2));
					closest = j;
				}
			}

			//Find neighbor
			GetConnectedVertices(decimatemodel,closest,list);
			//Get number of neighbor
			n = list->GetNumberOfIds();
			//Initialize ds
			ds[0] = 0; ds[1] = 0;
			for (j = 0; j < n; j++){
				ndecimatearray->GetTuple(list->GetId(j),stretch);
				decimatemodel->GetPoints()->GetPoint (list->GetId(j), point2);
				ds[0] += 1/(sqrt(pow(point2[2] - point1[2], 2)
					+ pow(point2[1] - point1[1], 2) + pow(point2[0] - point1[0], 2))+0.0001);
				ds[1] += *stretch /(sqrt(pow(point2[2] - point1[2], 2)
					+ pow(point2[1] - point1[1], 2) + pow(point2[0] - point1[0], 2))+0.0001);
			}
			//stretch measure relative change of the vertices spacing
			*stretch = ds[1] / ds[0];

            //ndecimatearray->GetTuple(closest,stretch);
			array->InsertTupleValue(i, stretch);
			list->Delete();
		}

		model->GetPointData()->SetScalars(array);

		vtkPolyDataWriter *writer = vtkPolyDataWriter::New();
		writer->SetInput(model);
		writer->SetFileName(output_name);
		writer->Write();

		delete []stretch;
		decimate->Delete();
		array->Delete();
		decimatearray->Delete();
		ndecimatearray->Delete();
        delete transform;
#else

cerr << argv[0] << " this program needs to be compiled with vtk enabled." << endl;
return 0;

#endif
	}

}

