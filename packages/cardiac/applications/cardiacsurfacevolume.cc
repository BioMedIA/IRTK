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
#include <irtkRegistration.h>

#ifdef HAS_VTK

#include <vtkPointData.h>
#include <vtkPolyData.h>
#include <vtkPolyDataReader.h>
#include <vtkPolyDataWriter.h>
#include <vtkDoubleArray.h>
#include <vtkPointLocator.h>
#include <vtkIdList.h>
#include <vtkLine.h>

char *input_name = NULL, *output_name = NULL, *image_name = NULL, *count_name = NULL;

void usage(){
	cerr << "Usage: cardiacsurfacevolume [input surface] [reference image] [outputfile] " << endl;
	cerr << "17 AHA segmentation a center line will be generated from the 17 apex and 0 basal, reference image is the manual segmentation, blood pool > 0" << endl;
	exit(1);
}

double vote(irtkGreyImage& image, int x, int y, int z){
	int pool[17],i,j,k,allcount,countzero,countnon,maxcountnon;
	for(i = 0; i < 17; i++){
		pool[i] = 0;
	}
	allcount = 0;
	for(i = x-1; i < x+2; i++){
		for(j = y-1; j < y+2; j++){
			for(k = z-1; k < z+2; k++){
				if((i < image.GetX()) && (i >= 0) && (j < image.GetY()) 
					&& (j >= 0) && (k < image.GetZ()) && (k >= 0)){
						if(image.GetAsDouble(i,j,k) >= 0 && image.GetAsDouble(i,j,k) < 17){
							pool[round(image.GetAsDouble(i,j,k))]++;
							allcount++;
						}
				}
			}
		}
	}
	countzero = pool[0];
	countnon = 0;
	maxcountnon = 0;
	k = 0;
	for(i = 1; i < 17; i++){
		countnon += pool[i];
		if(pool[i] > maxcountnon){
			maxcountnon = pool[i];
			k = i;
		}
	}
	if(countzero >= allcount-2)
		return 0;
	if(maxcountnon >= allcount-2)
		return k;
	if(image.GetAsDouble(x,y,z) != 0){
		return image.GetAsDouble(x,y,z);
	}
	if(countzero < countnon){
		return k;
	}
    return 0;
}

int main(int argc, char **argv)
{
	int i, j;
	double *profile,weight1,weight2,distance,x,y,z,**surfacepoints,**insertionpoints;
	double point[3],point1[3],point2[3],pinsertion[3],normal[3],volume[16],iregion[6];
	irtkGreyPixel *ptr,*tmpptr;
	vtkLine *linefunction;

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

	// Allocate memory for intensity profile
	profile = new double[1];

	// Read model
	vtkPolyDataReader *reader = vtkPolyDataReader::New();
	reader->SetFileName(input_name);
	reader->Modified();
	reader->Update();
	vtkPolyData *model = vtkPolyData::New();
	model = reader->GetOutput();
	model->Update();

	//Initialize pionts
	surfacepoints = new double*[model->GetNumberOfPoints()];
	insertionpoints = new double*[model->GetNumberOfPoints()];
	for (i = 0; i < model->GetNumberOfPoints(); i++) {
		surfacepoints[i] = new double[3];
		insertionpoints[i] = new double[3];
	}

	// Read image
	irtkGreyImage image;
	image.Read(image_name);
    irtkGreyImage segmentation;
    segmentation.Read(image_name);
	image.Initialize(image.GetImageAttributes());

	// Define point1 point2
	// Evaluate number of weights
	weight1 = 0; weight2 = 0;
	vtkDoubleArray *narray = vtkDoubleArray::New();
	narray->SetNumberOfTuples(model->GetNumberOfPoints());
	narray->SetNumberOfComponents(1);
	narray->DeepCopy(model->GetPointData()->GetScalars());
	for (i = 0; i < model->GetNumberOfPoints(); i++) {
		narray->GetTuple(i,profile);
		if(*profile == 0){
			weight1 ++;
		}
		if(*profile == 17){
			weight2 ++;
		}
	}
	point1[0] = 0; point1[1] = 0; point1[2] = 0;
	point2[0] = 0; point2[1] = 0; point2[2] = 0;
	for (i = 0; i < model->GetNumberOfPoints(); i++) {
		narray->GetTuple(i,profile);
		if(*profile == 0){
			model->GetPoint(i,point);
			for(j=0;j<3;j++){
				point1[j] += point[j]/weight1;
			}
		}
		if(*profile == 17){
			model->GetPoint(i,point);
			for(j=0;j<3;j++){
				point2[j] += point[j]/weight2;
			}
		}
	}

	// for all points insertion line point1 and point2
	linefunction = vtkLine::New();
	// Evaluate step size
	weight2 = min(min(image.GetXSize(),image.GetYSize()),image.GetZSize());
	irtkImageAttributes attr = image.GetImageAttributes();
	// start loop
	iregion[0] = attr._x; iregion[1] = attr._y; iregion[2] = attr._z;
	iregion[3] = 0; iregion[4] = 0; iregion[5] = 0;

	for (i = 0; i < model->GetNumberOfPoints(); i++) {
		//get point
		model->GetPoint(i,point);
		x = point[0];
		y = point[1];
		z = point[2];
		image.WorldToImage(x,y,z);
		if(round(x) < iregion[0] && round(x) >= 0){
			iregion[0] = round(x);
		}
		if(round(y) < iregion[1] && round(y) >= 0){
			iregion[1] = round(y);
		}
		if(round(z) < iregion[2] && round(z) >= 0){
			iregion[2] = round(z);
		}
		if(round(x) > iregion[3] && round(x) < attr._x){
			iregion[3] = round(x);
		}
		if(round(y) > iregion[4] && round(y) < attr._y){
			iregion[4] = round(y);
		}
		if(round(z) > iregion[5] && round(z) < attr._z){
			iregion[5] = round(z);
		}
		//Get insertion point to the center
		distance = linefunction->DistanceToLine(point,point1,point2,weight1,pinsertion);
		for(j = 0; j < 3; j++){
			surfacepoints[i][j] = point[j];
			insertionpoints[i][j] = pinsertion[j];
		}			
	}

	for(x = iregion[0]; x < iregion[3] + 1; x++){
		for(y = iregion[1]; y < iregion[4] + 1; y++){
			for(z = iregion[2]; z < iregion[5] + 1; z++){
				double maxdistance = 100000;
				point[0] = x; point[1] = y; point[2] = z;
				image.ImageToWorld(point[0],point[1],point[2]);
				if(segmentation.GetAsDouble(x,y,z) > 0){
                    for (i = 0; i < model->GetNumberOfPoints(); i++) {
                        for(j = 0; j < 3; j++){
                            point1[j] = surfacepoints[i][j];
                            point2[j] = insertionpoints[i][j];
                        }
                        distance = linefunction->DistanceToLine(point,point1,point2,weight1,pinsertion);
                        for(j = 0; j < 3; j++){
                            normal[j] = pinsertion[j] - point[j];
                        }
                        weight2 = pow(normal[0],2)+pow(normal[1],2)+pow(normal[2],2);
                        if(abs(weight2 - distance) < 0.1 && distance < maxdistance){
                            maxdistance = distance;
                            narray->GetTuple(i,profile);
                            if(*profile != 0 && *profile != 17){
                                image.PutAsDouble(x,y,z,*profile);
                            }
                        }
                    }
				}
			}
		}
	}

	delete []profile;
	linefunction->Delete();

	//finally fill the volume with neighbor vote
	irtkGreyImage tmpimage;
	for(i=0 ; i < 5; i++){
		tmpimage.Initialize(image.GetImageAttributes());
		ptr = image.GetPointerToVoxels();
		tmpptr = tmpimage.GetPointerToVoxels();
		for(j=0; j<image.GetNumberOfVoxels(); j++){
			*tmpptr = *ptr;
			ptr++;
			tmpptr++;
		}
		for(x = iregion[0]; x < iregion[3] + 1; x++){
			for(y = iregion[1]; y < iregion[4] + 1; y++){
				for(z = iregion[2]; z < iregion[5] + 1; z++){
					image.PutAsDouble(x,y,z,vote(tmpimage,x,y,z));
				}
			}
		}
	}


	// write volume out 
	char buffer[255];
	sprintf(buffer, "%s.nii", input_name);
	image.Write(buffer);

	// evaluate volume based on the projected image.
	for(i=1; i<17; i++){
		volume[i-1] = 0;
		ptr = image.GetPointerToVoxels();
		for(j=0; j<image.GetNumberOfVoxels(); j++){
			if(*ptr == i){
				volume[i-1] += image.GetXSize()*image.GetYSize()*image.GetZSize()/1000.0;
			}
			ptr++;
		}
	}

	ofstream fileout(output_name,ios::app);
	for(i=0; i<16; i++){
		fileout << volume[i] <<" ";
	}
	fileout << " 1 " << endl;
	fileout.close();

	for (i = 0; i < model->GetNumberOfPoints(); i++) {
		delete []surfacepoints[i];
		delete []insertionpoints[i];
	}
	delete []surfacepoints;
	delete []insertionpoints;
}

#else

#include <irtkImage.h>

int main( int argc, char *argv[] )
{
	cerr << argv[0] << " this program needs to be compiled with vtk enabled." << endl;
	return 0;
}

#endif
