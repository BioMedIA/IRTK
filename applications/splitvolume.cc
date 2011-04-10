/*=========================================================================

Library   : Image Registration Toolkit (IRTK)
Module    : $Id: makevolume.cc 2 2008-12-23 12:40:14Z dr $
Copyright : Imperial College, Department of Computing
Visual Information Processing (VIP), 2008 onwards
Date      : $Date: 2008-12-23 12:40:14 +0000 (‰∫? 23 ÂçÅ‰∫åÊú?2008) $
Version   : $Revision: 2 $
Changes   : $Author: dr $

=========================================================================*/

#include <irtkImage.h>

#ifdef USE_VXL
// need to include vxl here
#else
#include <nr.h>
#endif

void usage()
{
	cerr << "Usage: splitvolume [input] [out] <options>\n" << endl;
	cerr << "Options : -sub sub the volume with average intensity for late images, to mark out the infarct" <<endl;
	cerr << "Options : -ref [reference] use reference's orientation and origin" <<endl;
	cerr << "Options : -sequence split the volume into independent time frames" <<endl;
	cerr << "Options : -slice split the volume into independent slices (each slice can have a number of time frames)" <<endl;
	exit(1);
}

int main(int argc, char **argv)
{
	int t,x,y,z,ok;
	int subaverage = 0;
	int sequenceonly = 0;
	int sliceonly = 0;
	irtkGreyImage *ref = NULL;
	// Check command line
	if (argc < 3) {
		usage();
	}
	char *output = NULL;
	irtkGreyImage *image = new irtkGreyImage(argv[1]);
	argc--;
	argv++;
	output = argv[1];
	argc--;
	argv++;

	while (argc > 1) {
		ok = false;
		if ((ok == false) && (strcmp(argv[1], "-ref") == 0)) {
			argc--;
			argv++;
			ref = new irtkGreyImage(argv[1]);
			argc--;
			argv++;
			ok = true;
		}
		if ((ok == false) && (strcmp(argv[1], "-sub") == 0)) {
			argc--;
			argv++;
			subaverage = 1;
			ok = true;
		}
		if ((ok == false) && (strcmp(argv[1], "-sequence") == 0)) {
			argc--;
			argv++;
			sequenceonly = 1;
			ok = true;
		}
		if ((ok == false) && (strcmp(argv[1], "-slice") == 0)) {
			argc--;
			argv++;
			sliceonly = 1;
			ok = true;
		}
		if (ok == false) {
			cerr << "Can not parse argument " << argv[1] << endl;
			usage();
		}
	}  

	if(sliceonly){
		irtkImageAttributes attr;
			irtkImageAttributes attrt;
		for (z = 0; z < image->GetZ(); z++) {
					if(ref != NULL){
						attr = ref->GetImageAttributes();
						attrt = image->GetImageAttributes();
						image->PutOrientation(attr._xaxis,attr._yaxis,attr._zaxis);
						image->PutOrigin(attr._xorigin,attr._yorigin,attr._zorigin);
						attr._dx = attrt._dx; attr._dy = attrt._dy;
						attr._dz = attrt._dz; attr._x = attrt._x;
						attr._y = attrt._y; attr._z = 1;
					}else{
						attr = image->GetImageAttributes();
						attr._z = 1;
					}
					image->WorldToImage(attr._xorigin,attr._yorigin,attr._zorigin);
					attr._zorigin = attr._zorigin - image->GetZ() / 2.0 + 0.5 + z;
					image->ImageToWorld(attr._xorigin,attr._yorigin,attr._zorigin);
					irtkGreyImage *target = new irtkGreyImage(attr);
					target->PutOrigin(attr._xorigin,attr._yorigin,attr._zorigin);
					for (t = 0; t < image->GetT(); t++) {
						for (y = 0; y < target->GetY(); y++) {
							for (x = 0; x < target->GetX(); x++) {
								target->Put(x, y, 0, t, image->Get(x, y, z, t));
							}
						}
					}
					if ( subaverage){
						int average = target->GetAverage();
						for (t = 0; t < target->GetT(); t++) {
							for (y = 0; y < target->GetY(); y++) {
								for (x = 0; x < target->GetX(); x++) {
									int temp = target->Get(x, y, 0, 0) - average;
									if (temp < 0) temp = 0;
									target->Put(x, y, 0, 0, temp);
								}
							}
						}
					}
					char buffer[255];
					sprintf(buffer, "%s%.4d.nii", output, z);
					target->Write(buffer);
					delete target;
				}
	}else{
		for (t = 0; t < image->GetT(); t++) {
			// Combine images
			irtkImageAttributes attr;
			irtkImageAttributes attrt;
			if(sequenceonly){
				if(ref != NULL){
					attr = ref->GetImageAttributes();
					attr._t = 1; attr._dt = 1;
					attrt = image->GetImageAttributes();
					attr._dx = attrt._dx; attr._dy = attrt._dy;
					attr._dz = attrt._dz; attr._x = attrt._x;
					attr._y = attrt._y; attr._z = attrt._z;
				}else{
					attr = image->GetImageAttributes();
					attr._t = 1; attr._dt = 1;
				}
				irtkGreyImage *target = new irtkGreyImage(attr);

				for (z = 0; z < target->GetZ(); z++) {
					for (y = 0; y < target->GetY(); y++) {
						for (x = 0; x < target->GetX(); x++) {
							target->Put(x, y, z, 0, image->Get(x, y, z, t));
						}
					}
				}
				if ( subaverage){
					int average = target->GetAverage();
					for (z = 0; z < target->GetZ(); z++) {
						for (y = 0; y < target->GetY(); y++) {
							for (x = 0; x < target->GetX(); x++) {
								int temp = target->Get(x, y, z, 0) - average;
								if (temp < 0) temp = 0;
								target->Put(x, y, z, 0, temp);
							}
						}
					}
				}

				char buffer[255];
				sprintf(buffer, "%s%.2d.nii", output, t);
				target->Write(buffer);
				delete target;
			}else{
				for (z = 0; z < image->GetZ(); z++) {
					if(ref != NULL){
						attr = ref->GetImageAttributes();
						attr._t = 1; attr._dt = 1;
						attrt = image->GetImageAttributes();
						image->PutOrientation(attr._xaxis,attr._yaxis,attr._zaxis);
						image->PutOrigin(attr._xorigin,attr._yorigin,attr._zorigin);
						attr._dx = attrt._dx; attr._dy = attrt._dy;
						attr._dz = attrt._dz; attr._x = attrt._x;
						attr._y = attrt._y; attr._z = 1;
					}else{
						attr = image->GetImageAttributes();
						attr._t = 1; attr._dt = 1;
						attr._z = 1;
					}
					image->WorldToImage(attr._xorigin,attr._yorigin,attr._zorigin);
					attr._zorigin = attr._zorigin - image->GetZ() / 2.0 + 0.5 + z;
					image->ImageToWorld(attr._xorigin,attr._yorigin,attr._zorigin);
					irtkGreyImage *target = new irtkGreyImage(attr);
					target->PutOrigin(attr._xorigin,attr._yorigin,attr._zorigin);
					for (y = 0; y < target->GetY(); y++) {
						for (x = 0; x < target->GetX(); x++) {
							target->Put(x, y, 0, 0, image->Get(x, y, z, t));
						}
					}
					if ( subaverage){
						int average = target->GetAverage();
						for (y = 0; y < target->GetY(); y++) {
							for (x = 0; x < target->GetX(); x++) {
								int temp = target->Get(x, y, 0, 0) - average;
								if (temp < 0) temp = 0;
								target->Put(x, y, 0, 0, temp);
							}
						}
					}

					char buffer[255];
					sprintf(buffer, "%s%.2d%.2d.nii", output, t, z);
					target->Write(buffer);
					delete target;
				}
			}
		}
	}
	delete image;
}

