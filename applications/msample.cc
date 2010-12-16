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



char *input_name = NULL, *output_name = NULL, *image_name = NULL, *count_name = NULL;

void GetConnectedVertices(vtkSmartPointer<vtkPolyData> mesh, int seed, vtkSmartPointer<vtkIdList> connectedVertices)
{

	//get all cells that vertex 'seed' is a part of
	vtkSmartPointer<vtkIdList> cellIdList =
		vtkSmartPointer<vtkIdList>::New();
	mesh->GetPointCells(seed, cellIdList);

	//cout << "There are " << cellIdList->GetNumberOfIds() << " cells that use point " << seed << endl;
		//loop through all the cells that use the seed point
		for(vtkIdType i = 0; i < cellIdList->GetNumberOfIds(); i++)
		{

			vtkCell* cell = mesh->GetCell(cellIdList->GetId(i));
			//cout << "The cell has " << cell->GetNumberOfEdges() << " edges." << endl;

			//if the cell doesn't have any edges, it is a line
			if(cell->GetNumberOfEdges() <= 0)
			{
				continue;
			}

			for(vtkIdType e = 0; e < cell->GetNumberOfEdges(); e++)
			{
				vtkCell* edge = cell->GetEdge(e);

				vtkIdList* pointIdList = edge->GetPointIds();
				//cout << "This cell uses " << pointIdList->GetNumberOfIds() <<" points" << endl;
				/*
				for(vtkIdType p = 0; p < pointIdList->GetNumberOfIds(); p++)
				{
				cout << "Edge " << i << " uses point " << pointIdList->GetId(p) << endl;
				}
				*/
				if(pointIdList->GetId(0) == seed || pointIdList->GetId(1) == seed)
				{
					if(pointIdList->GetId(0) == seed)
					{
						connectedVertices->InsertUniqueId(pointIdList->GetId(1));
					}
					else
					{
						connectedVertices->InsertUniqueId(pointIdList->GetId(0));
					}
				}
			}
		}
		//cout << "There are " << connectedVertices->GetNumberOfIds() << " points connected to point " << seed << endl;
}

void usage(){
	cerr << "Usage: msample [input] [image] [output] " << endl;
	cerr << "-n [number of steps (default 10)]" << endl;
	cerr << "-ds [step size (default 1)]" << endl;
	cerr << "-scalar use one scalar instead of array on the norm"<<endl;
	cerr << "-abs take the abs value of the scalar"<<endl;
	cerr << "-count [another mesh] sample along between two mesh, evaluate percentage > 0"<<endl;
	cerr << "-blur blur the result or not"<<endl;
	cerr << "-fill [value] filling the holes on the surface with n iterations"<<endl;
	exit(1);
}

int main(int argc, char **argv)
{
	int i, j, k, l, n, m, son, abson, blur, fill, ok;
	double *profile, count;
	double x, y, z, ds, point[3], normal[3];

	if (argc < 4) {
		usage();
	}

	// Parse filenames
	input_name = argv[1];
	argv++;
	argc--;
	image_name = argv[1];
	argv++;
	argc--;
	output_name = argv[1];
	argv++;
	argc--;

	n  = 10;
	ds = 1;
	son = 0;
	abson = 0;
	count = 0;
	blur = 0;
	fill = 0;

	while (argc > 1) {
		ok = false;
		if ((ok == false) && (strcmp(argv[1], "-n") == 0)) {
			argc--;
			argv++;
			n = atoi(argv[1]);
			argc--;
			argv++;
			ok = true;
		}
		if ((ok == false) && (strcmp(argv[1], "-count") == 0)) {
			argc--;
			argv++;
			count_name = argv[1];
			argc--;
			argv++;
			ok = true;
		}
		if ((ok == false) && (strcmp(argv[1], "-ds") == 0)) {
			argc--;
			argv++;
			ds = atof(argv[1]);
			argc--;
			argv++;
			ok = true;
		}
		if ((ok == false) && (strcmp(argv[1], "-scalar") == 0)) {
			argc--;
			argv++;
			son = 1;
			ok = true;
		}
		if ((ok == false) && (strcmp(argv[1], "-abs") == 0)) {
			argc--;
			argv++;
			abson = 1;
			ok = true;
		}
		if ((ok == false) && (strcmp(argv[1], "-blur") == 0)) {
			argc--;
			argv++;
			blur = 1;
			ok = true;
		}
		if ((ok == false) && (strcmp(argv[1], "-fill") == 0)) {
			argc--;
			argv++;
			fill = atoi(argv[1]);
			argc--;
			argv++;
			ok = true;
		}
		if (ok == false) {
			cerr << "Can't parse argument " << argv[1] << endl;
			usage();
		}
	}

	// Allocate memory for intensity profile
	if(son)
		profile = new double[1];
	else
		profile = new double[2*n+1];

	// Read model
	vtkPolyDataReader *reader = vtkPolyDataReader::New();
	reader->SetFileName(input_name);
	reader->Modified();
	reader->Update();
	vtkPolyData *model = vtkPolyData::New();
	model = reader->GetOutput();
	model->Update();

	// Extract normals
	if (model->GetPointData()->GetNormals() == NULL) {
		cerr << "Model has no normals" << endl;
		exit(1);
	}

	// Read image
	irtkGreyImage image;
	image.Read(image_name);
	irtkNearestNeighborInterpolateImageFunction interpolator;
	interpolator.SetInput(&image);
	interpolator.Initialize();

	// Allocate memory
	if(count_name){
		vtkPolyDataReader *reader_count = vtkPolyDataReader::New();
		reader_count->SetFileName(count_name);
		reader_count->Modified();
		reader_count->Update();
		vtkPolyData *model_count = vtkPolyData::New();
		model_count = reader_count->GetOutput();
		model_count->Update();
		vtkDoubleArray *array = vtkDoubleArray::New();
		array->SetNumberOfTuples(model->GetNumberOfPoints());
		array->SetNumberOfComponents(1);
		array->SetName("IntensityCountProfile");
		vtkDoubleArray *narray = vtkDoubleArray::New();
		narray->SetNumberOfTuples(model->GetNumberOfPoints());
		narray->SetNumberOfComponents(1);
		narray->SetName("IntensityCountProfile");
		int count2;

		// Build a locator 
		vtkPointLocator *pointLocator = vtkPointLocator::New();
		pointLocator->SetDataSet(reader_count->GetOutput());
		pointLocator->BuildLocator();

		for (i = 0; i < model->GetNumberOfPoints(); i++) {
			double distance,point2[3];
			model->GetPoints()->GetPoint (i, point);
			vtkIdType ptId;
			ptId = pointLocator->FindClosestPoint(point);
			model_count->GetPoints()->GetPoint(ptId,point2);
			normal[0] = point2[0] - point[0];
			normal[1] = point2[1] - point[1];
			normal[2] = point2[2] - point[2];
			distance = sqrt(pow(normal[0],2)+pow(normal[1],2)+pow(normal[2],2));
			if(distance > 0){
				for(j = 0; j < 3; j++){
					normal[j] = normal[j] / distance;
				}
				ds = distance / n;
				count = 0; count2 = 0;
				if(!son){
					for (j = 0; j < n; j++) {
						x = point[0] + (j + 1) * ds * normal[0];
						y = point[1] + (j + 1) * ds * normal[1];
						z = point[2] + (j + 1) * ds * normal[2];
						image.WorldToImage(x,y,z);
						if(interpolator.Evaluate(x, y, z) >= 0){
							count += interpolator.Evaluate(x, y, z);
							count2 ++;
						}
					}
					if(count2 > 0){
						count = count / count2;
					}
				}else{
					x = point[0] + n / 2 * ds * normal[0];
					y = point[1] + n / 2 * ds * normal[1];
					z = point[2] + n / 2 * ds * normal[2];
					image.WorldToImage(x,y,z);
					if(interpolator.Evaluate(x, y, z) >= 0){
						count = interpolator.Evaluate(x, y, z);
					}
				}
			}			
			array->InsertTupleValue(i, &count);
		}
		for (i = 0; i < model->GetNumberOfPoints(); i++) {
			narray->InsertTupleValue(i, array->GetTuple(i));
		}
		// blur
		if(blur == 1){
			for (i = 0; i < model->GetNumberOfPoints(); i++) {
				vtkIdList *list = vtkIdList::New();
				//Find neighbor
				GetConnectedVertices(model,i,list);
				//Get number of neighbor
				n = list->GetNumberOfIds();
				double value = 0, *tmp;
				for (j = 0; j < n; j++){
					tmp = array->GetTuple(list->GetId(j));
					value += *tmp;
				}
				//stretch measure relative change of the vertices spacing
				count = value / n;
				narray->InsertTupleValue(i, &count);
				list->Delete();
			}
		}
		if(fill > 0){
			//find zero ids, and replace them using neighbor's value
			vtkIdList *zerolist = vtkIdList::New();
			for (i = 0; i < model->GetNumberOfPoints(); i++) {
				if(*array->GetTuple(i) <= 0)
					zerolist->InsertUniqueId(i);
			}
			for(l = 0; l < fill; l++){
				m = zerolist->GetNumberOfIds();
				for (k = 0; k < m; k++){
					vtkIdList *list = vtkIdList::New();
					//Find neighbor
					GetConnectedVertices(model,zerolist->GetId(k),list);
					//Get number of neighbor
					n = list->GetNumberOfIds();
					double value = 0, *tmp;
					for (j = 0; j < n; j++){
						tmp = array->GetTuple(list->GetId(j));
						if(*tmp > 0){
							value = *tmp;
						}
					}
					//stretch measure relative change of the vertices spacing
					if (value > 0){
						narray->InsertTupleValue(zerolist->GetId(k), &value);
						zerolist->DeleteId(zerolist->GetId(k));
					}
					list->Delete();
				}
				for (i = 0; i < model->GetNumberOfPoints(); i++) {
					array->InsertTupleValue(i, narray->GetTuple(i));
				}
			}
			zerolist->Delete();

			//find nonzerolist replace them using zero so only holes with radius smaller then the number of iterations are filled
			vtkIdList *nonzerolist = vtkIdList::New();
			for (i = 0; i < model->GetNumberOfPoints(); i++) {
				if(*array->GetTuple(i) > 0)
					nonzerolist->InsertUniqueId(i);
			}
			for(l = 0; l < fill; l++){
				m = nonzerolist->GetNumberOfIds();
				for (k = 0; k < m; k++){
					vtkIdList *list = vtkIdList::New();
					//Find neighbor
					GetConnectedVertices(model,nonzerolist->GetId(k),list);
					//Get number of neighbor
					n = list->GetNumberOfIds();
					double value = 1, *tmp;
					for (j = 0; j < n; j++){
						tmp = array->GetTuple(list->GetId(j));
						if(*tmp <= 0){
							value = *tmp;
						}
					}
					//stretch measure relative change of the vertices spacing
					if (value <= 0){
						narray->InsertTupleValue(nonzerolist->GetId(k), &value);
						nonzerolist->DeleteId(nonzerolist->GetId(k));
					}
					list->Delete();
				}
				for (i = 0; i < model->GetNumberOfPoints(); i++) {
					array->InsertTupleValue(i, narray->GetTuple(i));
				}
			}
			nonzerolist->Delete();

		}
		model->GetPointData()->SetScalars(narray);
		array->Delete();
		narray->Delete();
		reader_count->Delete();
	}else{
		vtkDoubleArray *array = vtkDoubleArray::New();
		array->SetNumberOfTuples(model->GetNumberOfPoints());
		array->SetNumberOfComponents(2*n+1);
		array->SetName("IntensityProfile");

		if(!son){
			for (i = 0; i < model->GetNumberOfPoints(); i++) {
				model->GetPoints()->GetPoint (i, point);
				model->GetPointData()->GetNormals()->GetTuple(i, normal);
				for (j = 0; j < 2*n+1; j++) {
					x = point[0] + (j - n) * ds * normal[0];
					y = point[1] + (j - n) * ds * normal[1];
					z = point[2] + (j - n) * ds * normal[2];
					image.WorldToImage(x,y,z);
					if(abson)
						profile[j] = abs(interpolator.Evaluate(x, y, z));
					else
						profile[j] = interpolator.Evaluate(x, y, z);
				}
				array->InsertTupleValue(i, profile);
			}
			model->GetPointData()->SetScalars(array);
			array->Delete();
		}else{
			vtkDoubleArray *array = vtkDoubleArray::New();
			array->SetNumberOfTuples(model->GetNumberOfPoints());
			array->SetNumberOfComponents(1);
			array->SetName("IntensityValue");

			for (i = 0; i < model->GetNumberOfPoints(); i++) {
				model->GetPoints()->GetPoint (i, point);
				image.WorldToImage(point[0], point[1], point[2]);
				if(abson)
					*profile = abs(interpolator.Evaluate(point[0], point[1], point[2]));
				else
					*profile = interpolator.Evaluate(point[0], point[1], point[2]);
				array->InsertTupleValue(i, profile);
			}
			model->GetPointData()->SetScalars(array);
			array->Delete();
		}
	}

	vtkPolyDataWriter *writer = vtkPolyDataWriter::New();
	writer->SetInput(model);
	writer->SetFileName(output_name);
	writer->Write();
	reader->Delete();

	delete []profile;
}

#else

#include <irtkImage.h>

int main( int argc, char *argv[] )
{
	cerr << argv[0] << " this program needs to be compiled with vtk enabled." << endl;
	return 0;
}

#endif
