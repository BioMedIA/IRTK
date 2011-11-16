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
#include <irtkImageFunction.h>
#include <irtkResampling.h>
#include <irtkResamplingWithPadding.h>

void usage()
{
	cerr << "Usage: combineimages [reference] [input 1 ... input n] [output]\n" << endl;
	exit(1);
}

int main(int argc, char **argv)
{
	int i, j, k, t, time, m;
	double x,y,z;

	// Determine how many volumes we have
	t = argc-3;

	if (t < 1) usage();

	cout << "Combine images from " << t << " images" << endl;

	irtkGreyImage* reference = new irtkGreyImage(argv[1]);

	double xsize, ysize, zsize, size;
	reference->GetPixelSize(&xsize, &ysize, &zsize);
	size = xsize;
	size = (size < ysize) ? size : ysize;
	size = (size < zsize) ? size : zsize;
	cerr << "Resampling reference to isotropic voxel size (in mm): "
		<< size << endl;
	irtkResamplingWithPadding<irtkGreyPixel> resampling(size, size, size,-1);
	resampling.SetInput ((irtkGreyImage*)(reference));
	resampling.SetOutput((irtkGreyImage*)(reference));
	resampling.Run();
	cerr << "done.."<<endl;

	irtkGreyImage* input = new irtkGreyImage[t];  

	int *input_x = new int[t];
	int *input_y = new int[t];
	int *input_z = new int[t];

	//add (1/2)
	irtkInterpolateImageFunction **interpolator = new irtkInterpolateImageFunction *[t];

	irtkRealImage* output = new irtkRealImage(reference->GetImageAttributes());

	// Read remaining images
	for (i = 0; i < t; i++) {

		cout << "Reading " << argv[i+2] << endl;
		input[i].Read(argv[i+2]);
		interpolator[i] = irtkInterpolateImageFunction::New(Interpolation_Linear, &input[i]);
		interpolator[i]->SetInput(&input[i]);
		interpolator[i]->Initialize();
		input_x[i] = input[i].GetX() - 1;
		input_y[i] = input[i].GetY() - 1;
		input_z[i] = input[i].GetZ() - 1;
	}


	cout << "Inserting volumes into sequence" << endl;
	// Setup the interpolator
	// Compute domain on which the linear interpolation is defined
	irtkRealPixel* ptr2tmp = NULL;

	for (time = 0; time < reference->GetT(); time++) {
		for (k = 0; k < reference->GetZ(); k++) {
			for (j = 0; j < reference->GetY(); j++) {
				for (i = 0; i < reference->GetX(); i++) {
					ptr2tmp = output->GetPointerToVoxels(i, j, k, time);
					*ptr2tmp = 1;
					x = i;
					y = j;
					z = k;
					reference->ImageToWorld(x, y, z);
					int tmp = 0;
					for (m = 0; m < t; m++) {
						input[m].WorldToImage(x, y, z);
						// Check whether transformed point is inside volume
						if ((x > 0) && (x < input_x[m]) &&
							(y > 0) && (y < input_y[m]) &&
							(z > 0) && (z < input_z[m])) {
								// Add sample to metric			  
								*ptr2tmp =  double(*ptr2tmp) * double(interpolator[m]->EvaluateInside(x, y, z, time)) / 100.0;
								tmp ++;
						}
						if (input_z[m] == 0 && round(z) == 0){
							if ((x > 0) && (x < input_x[m]) &&
								(y > 0) && (y < input_y[m])) {
									// Add sample to metric		
									*ptr2tmp =  double(*ptr2tmp) * double(interpolator[m]->EvaluateInside(x, y, 0, time)) / 100.0;
									tmp ++;
							}
						}
						input[m].ImageToWorld(x, y, z);
					}
					if(tmp != 0){
						*ptr2tmp = pow(double(*ptr2tmp) , 1.0 / double(tmp))*100.0;
					}else{
						*ptr2tmp = -1;
					}
				}
			}
		}
	}

	// Write image
	cout << "Writing output to " << argv[t+2] << endl;
	output->Write(argv[t+2]);

	delete reference;
	delete output;
	delete []input;
	delete []input_x;
	delete []input_y;
	delete []input_z;
	delete []interpolator;
}

