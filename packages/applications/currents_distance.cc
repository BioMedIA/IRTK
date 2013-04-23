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
#include <irtkTransformation.h>

#include <vtkFloatArray.h>
#include <vtkTriangle.h>
#include <vtkCellArray.h>
#include <vtkPointData.h>
#include <vtkPolyData.h>
#include <vtkPolyDataReader.h>
#include <vtkPolyDataWriter.h>
#include <vtkPointLocator.h>

#include <vtkTriangleFilter.h>
#include <vtkMath.h>

char *input_nameA = NULL;
char *input_nameB = NULL;
char *maskName = NULL;
char *dofinA_name = NULL;
char *resultout_name = NULL;

void usage()
{
	cerr << " currents_distance [inputA] [inputB]" << endl;
	cerr << " " << endl;
	cerr << " Estimate the 'currents' type distance between a pair of surfaces." << endl;
	cerr << " See Glaunes, IPMI 2005 and Linh Ha, MICCAI 2010. Implemented by Paul Aljabar" << endl;
	cerr << " " << endl;
	cerr << " Options:" << endl;
	cerr << " " << endl;
	cerr << " -mask name          : Name of scalars mask. Restrict distance summation to faces with positive values" << endl;
	cerr << "                       of the named mask. Scalars with the same name must be present in both surfaces" << endl;
	cerr << " -sigma val          : Kernel width in mm." << endl;
	cerr << " -kernelFraction val : Kernel width as a fraction of the mean of the estimated radii of each of the surfaces" << endl;
	cerr << " -verbose            : Give a bit more output." << endl;
	cerr << " -dofinA  transf     : Apply the given transformation to surface A before processing." << endl;
	cerr << " -output file        : Result output file" << endl;
	cerr << " " << endl;
	exit(1);
}

int main(int argc, char **argv)
{
	int i, j, jj, ind;
	bool ok, useMask;
	double searchRadius;
	double pt1[3], pt2[3], pt3[3], centroid[3];
	double e1[3], e2[3], e3[3];
	double n1[3], n2[3];
	double volA, volB, rA, rB, sigmaKer, kernelFraction;
	int nearestPtCount;
	int noOfFacesA, noOfFacesB;
	double boundsA[6], boundsB[6];
	double gaussConst1, gaussConst2;
	double val;

	bool verbose = false;

	// Size of kernel as a fraction of the 'radius' of the surfaces
	kernelFraction = 0.05;

	sigmaKer = -1.0;


	if (argc < 2){
		usage();
	}

	input_nameA = argv[1];
	argc--;
	argv++;
	input_nameB = argv[1];
	argc--;
	argv++;


	// Parse remaining arguments
	while (argc > 1){
		ok = false;
		if ((!ok) && (strcmp(argv[1], "-kernelFraction") == 0)) {
			argc--;
			argv++;
			kernelFraction = atof(argv[1]);
			argc--;
			argv++;
			ok = true;
		}
		if ((!ok) && (strcmp(argv[1], "-sigma") == 0)) {
			argc--;
			argv++;
			sigmaKer = atof(argv[1]);
			argc--;
			argv++;
			ok = true;
		}
		if ((!ok) && (strcmp(argv[1], "-mask") == 0)) {
			argc--;
			argv++;
			maskName = argv[1];
			argc--;
			argv++;
			ok = true;
		}
		if ((!ok) && (strcmp(argv[1], "-verbose") == 0)) {
			argc--;
			argv++;
			verbose = true;
			ok = true;
		}
		if ((ok == false) && (strcmp(argv[1], "-output") == 0)){
			argc--;
			argv++;
			resultout_name = argv[1];
			argc--;
			argv++;
			ok = true;
		}
		if ((!ok) && (strcmp(argv[1], "-dofinA") == 0)) {
			argc--;
			argv++;
			dofinA_name = argv[1];
			argc--;
			argv++;
			ok = true;
		}
		if (!ok){
			cerr << "Cannot parse argument " << argv[1] << endl;
			usage();
		}
	}

	vtkObject::GlobalWarningDisplayOff();

#ifdef  WIN32
	SetErrorMode(SEM_FAILCRITICALERRORS | SEM_NOGPFAULTERRORBOX);
#endif //  WIN32

	// Read surface
	vtkPolyDataReader *surface_readerA = vtkPolyDataReader::New();
	surface_readerA->SetFileName(input_nameA);
	surface_readerA->Modified();
	surface_readerA->Update();

	vtkPolyDataReader *surface_readerB = vtkPolyDataReader::New();
	surface_readerB->SetFileName(input_nameB);
	surface_readerB->Modified();
	surface_readerB->Update();


	vtkPolyData *surfaceA;
	surfaceA = surface_readerA->GetOutput();
	surfaceA->Update();
	surfaceA->BuildCells();
	surfaceA->BuildLinks();

	vtkPolyData *surfaceB;
	surfaceB = surface_readerB->GetOutput();
	surfaceB->Update();
	surfaceB->BuildCells();
	surfaceB->BuildLinks();

	if (dofinA_name != NULL){
		// Apply the transformation to the points of surface A.
		irtkTransformation *transformation = irtkTransformation::New(dofinA_name);
		for (i = 0; i < surfaceA->GetNumberOfPoints(); i++) {
			(surfaceA->GetPoints())->GetPoint(i, pt1);
			transformation->Transform(pt1[0], pt1[1], pt1[2]);
			surfaceA->GetPoints()->SetPoint(i, pt1);
		}
		surfaceA->Modified();
	}

	// Derive a new polydata set with points equal to the
	// centers of the faces of the input surfaces and vectors
	// associated with each point that are the normals of the faces.

	vtkTriangleFilter *triFilterA = vtkTriangleFilter::New();
	triFilterA->SetInput(surfaceA);
	triFilterA->Update();
	vtkPolyData *trianglesA = vtkPolyData::New();
	trianglesA = triFilterA->GetOutput();
	trianglesA->BuildCells();
	trianglesA->BuildLinks();

	vtkTriangleFilter *triFilterB = vtkTriangleFilter::New();
	triFilterB->SetInput(surfaceB);
	triFilterB->Update();
	vtkPolyData *trianglesB = vtkPolyData::New();
	trianglesB = triFilterB->GetOutput();
	trianglesB->BuildCells();
	trianglesB->BuildLinks();





	vtkFloatArray *maskA = vtkFloatArray::New();
	vtkFloatArray *maskB = vtkFloatArray::New();

	useMask = false;

	if (maskName != NULL){
		maskA = (vtkFloatArray*) trianglesA->GetPointData()->GetArray(maskName, ind);
		if (ind == -1 || maskA == NULL){
			cerr << "Scalars unavailable with name " << maskName << " in surface " << input_nameA << endl;
			exit(1);
		}
		maskB = (vtkFloatArray*) trianglesB->GetPointData()->GetArray(maskName, ind);
		if (ind == -1 || maskB == NULL){
			cerr << "Scalars unavailable with name " << maskName << " in surface " << input_nameB << endl;
			exit(1);
		}

		useMask =true;
	}



	trianglesA->GetBounds(boundsA);
	trianglesB->GetBounds(boundsB);

	if (boundsA[0] > boundsB[1] ||
		boundsB[0] > boundsA[1] ||
		boundsA[2] > boundsB[3] ||
		boundsB[2] > boundsA[3] ||
		boundsA[4] > boundsB[5] ||
		boundsB[4] > boundsA[5]){
			cout << "Surfaces' bounding boxes do not intersect." << endl;
			exit(1);
	}



	vtkCellArray *facesA = trianglesA->GetPolys();
	noOfFacesA = facesA->GetNumberOfCells();

	vtkCellArray *facesB = trianglesB->GetPolys();
	noOfFacesB = facesB->GetNumberOfCells();


	vtkPoints *centresA = vtkPoints::New();
	centresA->SetNumberOfPoints(noOfFacesA);

	vtkPoints *centresB = vtkPoints::New();
	centresB->SetNumberOfPoints(noOfFacesB);


	vtkFloatArray *normalsA = vtkFloatArray::New();
	normalsA->SetNumberOfComponents(3);
	normalsA->SetNumberOfTuples(noOfFacesA);

	vtkFloatArray *normalsB = vtkFloatArray::New();
	normalsB->SetNumberOfComponents(3);
	normalsB->SetNumberOfTuples(noOfFacesB);


	vtkFloatArray *maskOutA = vtkFloatArray::New();
	maskOutA->SetNumberOfComponents(1);
	maskOutA->SetNumberOfTuples(noOfFacesA);

	vtkFloatArray *maskOutB = vtkFloatArray::New();
	maskOutB->SetNumberOfComponents(1);
	maskOutB->SetNumberOfTuples(noOfFacesB);




	vtkIdType nptsForFace = 0;
	vtkIdType *ptIdsForFace;
	vtkIdType ptID1, ptID2, ptID3;


	facesA->InitTraversal();

	for (i = 0; i < noOfFacesA; ++i){

		facesA->GetNextCell(nptsForFace, ptIdsForFace);

		if (nptsForFace != 3){
			cerr << "Error: faces must have three points each." << endl;
			exit(0);
		}

		ptID1 = ptIdsForFace[0];
		ptID2 = ptIdsForFace[1];
		ptID3 = ptIdsForFace[2];

		trianglesA->GetPoint(ptID1, pt1);
		trianglesA->GetPoint(ptID2, pt2);
		trianglesA->GetPoint(ptID3, pt3);

		for (j = 0; j < 3; ++j){
			e1[j] = pt2[j] - pt1[j];
			e2[j] = pt3[j] - pt1[j];
			centroid[j] = (pt1[j] + pt2[j] + pt3[j]) / 3.0;
		}

		centresA->SetPoint(i, centroid);

		vtkMath::Cross(e1, e2, e3);
		normalsA->SetTuple3(i, e3[0], e3[1], e3[2]);

		val = 0;
		if (useMask){
			if( maskA->GetTuple1(ptID1) > 0 )
				val++;
			if( maskA->GetTuple1(ptID2) > 0 )
				val++;
			if( maskA->GetTuple1(ptID3) > 0 )
				val++;
		}

		maskOutA->SetTuple1(i, (val > 1) ? 1 : 0);

	}


	facesB->InitTraversal();

	// Note where this loop starts and ends.
	for (i = 0; i < noOfFacesB; ++i){

		facesB->GetNextCell(nptsForFace, ptIdsForFace);

		if (nptsForFace != 3){
			cerr << "Error: faces must have three points each." << endl;
			exit(0);
		}

		ptID1 = ptIdsForFace[0];
		ptID2 = ptIdsForFace[1];
		ptID3 = ptIdsForFace[2];

		trianglesB->GetPoint(ptID1, pt1);
		trianglesB->GetPoint(ptID2, pt2);
		trianglesB->GetPoint(ptID3, pt3);

		for(j = 0; j < 3; ++j){
			e1[j] = pt2[j] - pt1[j];
			e2[j] = pt3[j] - pt1[j];
			centroid[j] = (pt1[j] + pt2[j] + pt3[j]) / 3.0;
		}

		centresB->SetPoint(i, centroid);

		// Cross product is skew symmetric for the second surface, we do e2 x e1 while we did
		// e1 x e2 for the first surface - i.e.  the sense is negated.
		vtkMath::Cross(e2, e1, e3);
		normalsB->SetTuple3(i, e3[0], e3[1], e3[2]);

		val = 0;
		if (useMask){
			if( maskB->GetTuple1(ptID1) > 0 )
				val++;
			if( maskB->GetTuple1(ptID2) > 0 )
				val++;
			if( maskB->GetTuple1(ptID3) > 0 )
				val++;
		}
		maskOutB->SetTuple1(i, (val > 1) ? 1 : 0);
	}


	normalsA->SetName("faceNormals");
	maskOutA->SetName("maskOut");

	normalsB->SetName("faceNormals");
	maskOutB->SetName("maskOut");


	vtkPolyData *currentA = vtkPolyData::New();
	currentA->SetPoints(centresA);
	currentA->GetPointData()->AddArray(normalsA);
	currentA->GetPointData()->AddArray(maskOutA);
	currentA->Update();

	vtkPolyData *currentB = vtkPolyData::New();
	currentB->SetPoints(centresB);
	currentB->GetPointData()->AddArray(normalsB);
	currentB->GetPointData()->AddArray(maskOutB);
	currentB->Update();


	vtkPointLocator *point_locatorA = vtkPointLocator::New();
	point_locatorA->SetNumberOfPointsPerBucket(5);
	point_locatorA->SetDataSet(currentA);
	point_locatorA->BuildLocator();

	vtkPointLocator *point_locatorB = vtkPointLocator::New();
	point_locatorB->SetNumberOfPointsPerBucket(5);
	point_locatorB->SetDataSet(currentB);
	point_locatorB->BuildLocator();


	if (sigmaKer < 0){
		// Need to establish a suitable radius for the kernel on centre
		// to centre distances.

		volA = (boundsA[1] - boundsA[0]) *
			(boundsA[3] - boundsA[2]) *
			(boundsA[5] - boundsA[4]);

		volB = (boundsB[1] - boundsB[0]) *
			(boundsB[3] - boundsB[2]) *
			(boundsB[5] - boundsB[4]);

		// Assume points distributed in a roughly spherical arrangement.
		rA = pow(volA * 3.0 / 4.0 / M_PI, (1.0/3.0));
		rB = pow(volB * 3.0 / 4.0 / M_PI, (1.0/3.0));


		sigmaKer = 0.5 * (rA + rB) * kernelFraction;

		if (verbose){
			cout << "Estimated radii  : " << rA << " and " << rB << endl;
			cout << "Sigma for kernel : " << sigmaKer << endl;
		}
	} else {
		if (verbose){
			cout << "Sigma for kernel given : " << sigmaKer << endl;
		}
	}

	// Constants for kernel
	gaussConst1 = 1 / sigmaKer / sqrt(2 * M_PI);
	gaussConst2 = -1.0 / sigmaKer / sigmaKer;

	searchRadius = 2.5 * sigmaKer;


	vtkIdList *nearestPtIDs = vtkIdList::New();

	double pt2ptDistSq, totalDist;
	double aDotB, bDotA, aDotA, bDotB;

	// The currents distance that we seek to measure and return.
	totalDist = 0.0;


	// inter surface distance:
	for (i = 0; i < noOfFacesA; ++i){

		if (useMask && maskOutA->GetTuple1(i) <= 0){
			continue;
		}

		currentA->GetPoint(i, pt1);

		point_locatorB->FindPointsWithinRadius(searchRadius, pt1, nearestPtIDs);

		nearestPtCount = nearestPtIDs->GetNumberOfIds();

		// Current normal.
		normalsA->GetTuple(i, n1);

		for (j = 0; j < nearestPtCount; ++j){
			jj = nearestPtIDs->GetId(j);

			currentB->GetPoint(jj, pt2);
			normalsB->GetTuple(jj, n2);

			pt2ptDistSq = vtkMath::Distance2BetweenPoints(pt1, pt2);

			val = gaussConst1 * exp(gaussConst2 * pt2ptDistSq);
			val *= vtkMath::Dot(n1, n2);

			totalDist += val;

		}
	}

	aDotB = totalDist;
	totalDist = 0.0;

	for (i = 0; i < noOfFacesB; ++i){

		if (useMask && maskOutB->GetTuple1(i) <= 0){
			continue;
		}

		currentB->GetPoint(i, pt1);

		point_locatorA->FindPointsWithinRadius(searchRadius, pt1, nearestPtIDs);

		nearestPtCount = nearestPtIDs->GetNumberOfIds();

		// Current normal.
		normalsB->GetTuple(i, n1);

		for (j = 0; j < nearestPtCount; ++j){
			jj = nearestPtIDs->GetId(j);

			currentA->GetPoint(jj, pt2);
			normalsA->GetTuple(jj, n2);

			pt2ptDistSq = vtkMath::Distance2BetweenPoints(pt1, pt2);

			val = gaussConst1 * exp(gaussConst2 * pt2ptDistSq);
			val *= vtkMath::Dot(n1, n2);

			totalDist += val;

		}
	}

	bDotA = totalDist;
	totalDist = 0.0;

	// Within dists for checking : to drop

	for (i = 0; i < noOfFacesA; ++i){

		if (useMask && maskOutA->GetTuple1(i) <= 0){
			continue;
		}

		currentA->GetPoint(i, pt1);

		point_locatorA->FindPointsWithinRadius(searchRadius, pt1, nearestPtIDs);

		nearestPtCount = nearestPtIDs->GetNumberOfIds();

		// Current normal.
		normalsA->GetTuple(i, n1);

		for (j = 0; j < nearestPtCount; ++j){
			jj = nearestPtIDs->GetId(j);

			currentA->GetPoint(jj, pt2);
			normalsA->GetTuple(jj, n2);

			pt2ptDistSq = vtkMath::Distance2BetweenPoints(pt1, pt2);

			val = gaussConst1 * exp(gaussConst2 * pt2ptDistSq);
			val *= vtkMath::Dot(n1, n2);

			totalDist += val;

		}
	}

	aDotA = totalDist;
	totalDist = 0.0;

	for (i = 0; i < noOfFacesB; ++i){

		if (useMask && maskOutB->GetTuple1(i) <= 0){
			continue;
		}

		currentB->GetPoint(i, pt1);

		point_locatorB->FindPointsWithinRadius(searchRadius, pt1, nearestPtIDs);

		nearestPtCount = nearestPtIDs->GetNumberOfIds();

		// Current normal.
		normalsB->GetTuple(i, n1);

		for (j = 0; j < nearestPtCount; ++j){
			jj = nearestPtIDs->GetId(j);

			currentB->GetPoint(jj, pt2);
			normalsB->GetTuple(jj, n2);

			pt2ptDistSq = vtkMath::Distance2BetweenPoints(pt1, pt2);

			val = gaussConst1 * exp(gaussConst2 * pt2ptDistSq);
			val *= vtkMath::Dot(n1, n2);

			totalDist += val;

		}
	}

	bDotB = totalDist;

	// Various flavours of measures:
	cout << aDotA << " " << bDotB << " " << aDotB << " " << bDotA << endl;
	cout << aDotA - aDotB - bDotA + bDotB << endl;
	cout << aDotA + aDotB + bDotA + bDotB << endl;
	cout << sqrt(fabs(aDotB + bDotA)) << endl;

	if(resultout_name){
		cerr << "Writing Results: " << resultout_name << endl;
		ofstream fout(resultout_name,ios::app);	  
		fout << sqrt(fabs(aDotB + bDotA)) << endl;
		fout.close();
	}

	return 0;
}

#else

#include <irtkImage.h>

int main( int argc, char *argv[] ){
	cerr << argv[0] << " needs to be compiled with the VTK library " << endl;
}
#endif
