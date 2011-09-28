/*=========================================================================

Library   : Image Registration Toolkit (IRTK)
Module    : $Id$
Copyright : Imperial College, Department of Computing
Visual Information Processing (VIP), 2008 onwards
Date      : $Date$
Version   : $Revision$
Changes   : $Author$

=========================================================================*/

#include <irtkTagFunction.h>

#include <irtkImageToFile.h>

void usage()
{
	cerr << "Usage: tracktag [inputsequence] [last landmark] [output prefix include directory] [ds(searching space)] [track 0 or center 1 or center after track 2]\n" << endl;
	exit(1);
}

int main(int argc, char **argv)
{
	
	int i, j, x, y, z, t, toggle;
	char buffer[255];
	if(argc < 5){
		usage();
		exit(1);
	}

	// Determine how many volumes we have
	irtkGreyImage inputsequence;

	cout << "Reading image: " << argv[1] << endl;
	inputsequence.Read(argv[1]);

	irtkImageAttributes inputSeqAttributes=inputsequence.GetImageAttributes();

	t = inputSeqAttributes._t;

	inputSeqAttributes._t = 1;

	irtkGreyImage* input = new irtkGreyImage[t];

	// Read last landmark
	irtkPointSet last;
	irtkPointSet *rest = new irtkPointSet[t-1];
	cout << "Reading landmark: " << argv[2] << endl;
	last.ReadVTK(argv[2]);

	char* filename;
	filename = argv[3];

	int ds;
	ds = atoi(argv[4]);

	if(argc == 6)
		toggle = atoi(argv[5]);
	else
		toggle = 0;

	if (toggle == 1){
		for (i=0;i<t-1;i++){
			sprintf(buffer, "%s%.2d.vtk",filename,i+1);
			rest[i].ReadVTK(buffer);
		}
	}

	cout << "Tracking tags" << endl;
	// initialize last image.
	input[t-1].Initialize(inputSeqAttributes);
	for (z = 0; z < inputSeqAttributes._z; z++) {
		for (y = 0; y < inputSeqAttributes._y; y++) {
			for (x = 0; x < inputSeqAttributes._x; x++) {
				input[t-1](x, y, z) = inputsequence(x,y,z,t-1);
			}
		}
	}
	if(toggle != 1){
		irtkTagFunction *tracktag = new irtkTagFunction();
		for (i = t-2; i >= 0; i--) {
			cout << "Tracking Volume: " << i+1 << " ..." << endl;
			input[i].Initialize(inputSeqAttributes);
			for (z = 0; z < inputSeqAttributes._z; z++) {
				for (y = 0; y < inputSeqAttributes._y; y++) {
					for (x = 0; x < inputSeqAttributes._x; x++) {
						input[i](x, y, z) = inputsequence(x,y,z,i);
					}
				}
			}			

			if(i == t - 2){
				tracktag->SetPointSet(last);
			}else{
				tracktag->SetPointSet(rest[i + 1]);
			}

			tracktag->SetInput(NULL,NULL,&input[i+1],&input[i],NULL);
			
			tracktag->Track(toggle,ds);

			rest[i-1] = tracktag->GetOutput();
			
			cout << "Output Landmark " << i << " ..." << endl;
			sprintf(buffer, "%s%.2d.vtk",filename,i);
			rest[i].WriteVTK(buffer);
		}
		delete tracktag;
	}else{
		for (i = 1; i < t; i++) {
			cout << "Centering Volume: " << i+1 << " ..." << endl;
			input[i].Initialize(inputSeqAttributes);
			for (z = 0; z < inputSeqAttributes._z; z++) {
				for (y = 0; y < inputSeqAttributes._y; y++) {
					for (x = 0; x < inputSeqAttributes._x; x++) {
						input[i](x, y, z) = inputsequence(x,y,z,i);
					}
				}
			}
			for (j = 0; j < last.Size(); j++){
				input[i].GravityCenter(rest[i-1](j),ds);
			}
			cout << "Output Landmark " << i << " ..." << endl;
			char buffer[255];
			sprintf(buffer, "%s%.2d.vtk",filename,i);
			rest[i-1].WriteVTK(buffer);
		}
	}

	delete[] input;
	delete[] rest;
}

