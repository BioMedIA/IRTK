#ifdef HAS_VTK
#define Normal 0
#define Radial 1
#define Longitudinal 2
#define Circumferential 3
#include <irtkImage.h>
#include <irtkImageFunction.h>
#include <irtkTransformation.h>

// Default filenames
char *endo_name = NULL, *epi_name = NULL, *input_name = NULL, *output_name = NULL, *axis_name = NULL;

void usage()
{
  cerr << "Usage: cardiacdof2vtk [dof] [vtk] [endo] [epi] \n"<<             
	      "Load dof from file output dof in vtk format between endo and epi cardial surfaces\n" << endl;
  cerr << "<-axis [landmarks]>       Define longitudinal axis based on two landmarks output in radial circumferential and longitudinal coorndinate\n" << endl;										  ;
  cerr << "<-longitudinal>           Define output in longitudinal direction\n" <<endl;
  cerr << "<-circumferential>        Define output in circumferential direction\n" <<endl;
  cerr << "<-radial>                 Define output in radial direction\n" <<endl;
  cerr << "<-strain>				 Output Strain" <<endl;
  exit(1);
}

int main(int argc, char **argv)
{
  int i, j, x, y, z;
  double p1[3], p2[3],normal[3],axis[3],circle[3];
  double ds,distance;
  int ok,mode,strain;

  if (argc < 5) {
    usage();
  }
  mode = Normal;
  strain = false;

  input_name  = argv[1];
  argc--;
  argv++;
  output_name = argv[1];
  argc--;
  argv++;
  endo_name = argv[1];
  argc--;
  argv++;
  epi_name = argv[1];
  argc--;
  argv++;

  while (argc > 1) {
	  ok = false;
	  if (strcmp(argv[1], "-strain") == 0) {
		  argc--;
		  argv++;
		  strain = true;
		  ok = true;	 
	  } else if (strcmp(argv[1], "-radial") == 0) {
		  argc--;
		  argv++;
		  mode = Radial;
		  ok = true;	 
	  } else if (strcmp(argv[1], "-longitudinal") == 0) {
		  argc--;
		  argv++;
		  mode = Longitudinal;
		  ok = true;	 
	  } else if (strcmp(argv[1], "-circumferential") == 0) {
		  argc--;
		  argv++;
		  mode = Circumferential;
		  ok = true;	 
	  } else if (strcmp(argv[1], "-axis") == 0) {
		  argc--;
		  argv++;
		  axis_name = argv[1];
		  argc--;
		  argv++;
		  ok = true;	 
      } else if (!ok) {
		  cerr << "Invalid option : " << argv[1] << endl;
		  exit(1);
	  }
  }

  // Read transformation
  irtkTransformation *transform = irtkTransformation::New(input_name);

  // Set up vtk points
  vtkPoints *points = vtkPoints::New();

  // Set up vtk vectors
  vtkDoubleArray *vectors = vtkDoubleArray::New();
  vectors->SetNumberOfComponents(3);
  // Set up strain value
  vtkDoubleArray *strainvectors = vtkDoubleArray::New();
  strainvectors->SetNumberOfComponents(1);
  // Set up strain tensor
  vtkDoubleArray *straintensor = vtkDoubleArray::New();
  straintensor->SetNumberOfComponents(9);

  // Read model
  vtkPolyDataReader *reader = vtkPolyDataReader::New();
  reader->SetFileName(endo_name);
  reader->Modified();
  reader->Update();
  vtkPolyData *model = vtkPolyData::New();
  model = reader->GetOutput();
  model->Update();

  vtkPolyDataReader *reader_count = vtkPolyDataReader::New();
  reader_count->SetFileName(epi_name);
  reader_count->Modified();
  reader_count->Update();
  vtkPolyData *model_count = vtkPolyData::New();
  model_count = reader_count->GetOutput();
  model_count->Update();

  irtkPointSet axislandmarks;
  if(axis_name != NULL){
	  axislandmarks.ReadVTK(axis_name);
	  if(axislandmarks.Size() == 2){
		  axis[0] = axislandmarks(1)._x - axislandmarks(0)._x;
		  axis[1] = axislandmarks(1)._y - axislandmarks(0)._y;
		  axis[2] = axislandmarks(1)._z - axislandmarks(0)._z;
		  distance = sqrt(pow(axis[0],2)+pow(axis[1],2)+pow(axis[2],2));
		  if(distance > 0){
			  for(j = 0; j < 3; j++){
				  axis[j] = axis[j] / distance;
			  }
		  }
	  }else{
		  cerr << "Longitudinal axis needs to be defined by 2 points, apex and mid valve" << endl;
		  exit(1);
	  }
  }else{
	  if(mode == Longitudinal || mode == Circumferential){
		cerr << "Longitudinal and Circumferential motion calculation needs definition of longitudinal axis" <<endl;
		exit(1);
	  }
  }

  // Build a locator 
  vtkPointLocator *pointLocator = vtkPointLocator::New();
  pointLocator->SetDataSet(reader_count->GetOutput());
  pointLocator->BuildLocator();

  for (i = 0; i < model->GetNumberOfPoints(); i++) {
	  model->GetPoints()->GetPoint (i, p1);
	  vtkIdType ptId;
	  ptId = pointLocator->FindClosestPoint(p1);
	  model_count->GetPoints()->GetPoint(ptId,p2);
	  normal[0] = p2[0] - p1[0];
	  normal[1] = p2[1] - p1[1];
	  normal[2] = p2[2] - p1[2];
	  distance = sqrt(pow(normal[0],2)+pow(normal[1],2)+pow(normal[2],2));
	  if(distance > 0){
		  for(j = 0; j < 3; j++){
			  normal[j] = normal[j] / distance;
		  }
		  ds = distance / 2;
		  x = p1[0] + ds * normal[0];
		  y = p1[1] + ds * normal[1];
		  z = p1[2] + ds * normal[2];
		  p1[0] = x;
		  p1[1] = y;
		  p1[2] = z;
		  p2[0] = p1[0];
		  p2[1] = p1[1];
		  p2[2] = p1[2];
		  transform->Transform(p2[0], p2[1], p2[2]);
		  p2[0] -= p1[0];
		  p2[1] -= p1[1];
		  p2[2] -= p1[2];
		  points->InsertNextPoint(p1);

		  if(mode == Radial){
			  ds = p2[0]*normal[0] + p2[1]*normal[1] + p2[2]*normal[2];
			  p2[0] = ds * normal[0];
			  p2[1] = ds * normal[1];
			  p2[2] = ds * normal[2];
		  }else if(mode == Longitudinal){
			  ds = p2[0]*axis[0] + p2[1]*axis[1] + p2[2]*axis[2];
			  p2[0] = ds * axis[0];
			  p2[1] = ds * axis[1];
			  p2[2] = ds * axis[2];
		  }else if(mode == Circumferential){
			  //(a2b3 ? a3b2, a3b1 ? a1b3, a1b2 ? a2b1)
			  circle[0] = axis[1]*normal[2] - axis[2]*normal[1];
			  circle[1] = axis[2]*normal[0] - axis[0]*normal[2];
			  circle[2] = axis[0]*normal[1] - axis[1]*normal[0];
			  ds = p2[0]*circle[0] + p2[1]*circle[1] + p2[2]*circle[2];
			  p2[0] = ds * circle[0];
			  p2[1] = ds * circle[1];
			  p2[2] = ds * circle[2];
		  }

		  vectors->InsertNextTuple(p2);

		  if(strain){
			 double result_of_strain,tensor_of_strain[9];
			 irtkMatrix jac,strainm;
			 jac.Initialize(3, 3);
			 strainm.Initialize(3,3);
			 transform->LocalJacobian(jac,p1[0],p1[1],p1[2]);
			 strainm = jac;
			 strainm.Transpose();
			 strainm = strainm*jac;

			 for(x=0;x<9;x++){
				 tensor_of_strain[x] = strainm(x/3,x%3);
			 }
			 strainm(0,0) = strainm(0,0) - 1;
			 strainm(1,1) = strainm(1,1) - 1;
			 strainm(2,2) = strainm(2,2) - 1;

			 irtkMatrix spt,sp;
			 spt.Initialize(1,3);
			 sp.Initialize(3,1);
			 if(mode == Radial){
				 sp(0,0) = normal[0]; spt(0,0) = normal[0];
				 sp(1,0) = normal[1]; spt(0,1) = normal[1];
				 sp(2,0) = normal[2]; spt(0,2) = normal[2];
				 strainm = spt*strainm*sp;
			 }else if(mode == Longitudinal){
				 sp(0,0) = axis[0]; spt(0,0) = axis[0];
				 sp(1,0) = axis[1]; spt(0,1) = axis[1];
				 sp(2,0) = axis[2]; spt(0,2) = axis[2];
				 strainm = spt*strainm*sp;
			 }else if(mode == Circumferential){
				 sp(0,0) = circle[0]; spt(0,0) = circle[0];
				 sp(1,0) = circle[1]; spt(0,1) = circle[1];
				 sp(2,0) = circle[2]; spt(0,2) = circle[2];
				 strainm = spt*strainm*sp;
			 }
			 result_of_strain = strainm.Det();
			 strainvectors->InsertNextTuple(&result_of_strain);
			 straintensor->InsertNextTuple(tensor_of_strain);
		  }
	  }			
  }
  vtkPolyData *output = vtkPolyData::New();
  output->SetPoints(points);
  output->GetPointData()->SetVectors(vectors);
  if(strain){
	  output->GetPointData()->SetScalars(strainvectors);
	  output->GetPointData()->SetTensors(straintensor);
  }
  vtkPolyDataWriter *writer = vtkPolyDataWriter::New();
  writer->SetInput(output);
  writer->SetFileName(output_name);
  writer->Write();
  writer->Delete();
  reader->Delete();
  reader_count->Delete();
  delete transform;
}

#else

#include <irtkImage.h>

int main( int argc, char *argv[] )
{
  cerr << argv[0] << " needs to be compiled with the VTK library " << endl;
}
#endif
