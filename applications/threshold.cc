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

void usage()
{
	cerr << "Usage: threshold [input] [output] [threshold]" << endl;
	exit(1);
}

template<typename VoxelType> void threshold(irtkGenericImage<VoxelType> &image, double t)
{
	int x, y, z;

	for (z = 0; z < image.GetZ(); z++) {
		for (y = 0; y < image.GetY(); y++) {
			for (x = 0; x < image.GetX(); x++) {
				if (image(x, y, z) <= t) {
					image(x, y, z) = 0;
				} else {
					image(x, y, z) = 1;
				}
			}
		}
	}

}

int main(int argc, char **argv)
{
	double t;

	if (argc != 4) {
		usage();
	}

	// Read image
	irtkImage *image = irtkImage::New(argv[1]);

	// Read threshold value
	t = atof(argv[3]);

	// Apply threshold
	switch (image->GetScalarType()) {
		case IRTK_VOXEL_CHAR: {
			threshold(*(dynamic_cast<irtkGenericImage<char> *>(image)), t);
			break;
		}
		case IRTK_VOXEL_UNSIGNED_CHAR: {
			threshold(*(dynamic_cast<irtkGenericImage<unsigned char> *>(image)), t);
			break;
		}
		case IRTK_VOXEL_SHORT: {
			threshold(*(dynamic_cast<irtkGenericImage<short> *>(image)), t);
			break;
		}
		case IRTK_VOXEL_UNSIGNED_SHORT: {
			threshold(*(dynamic_cast<irtkGenericImage<unsigned short> *>(image)), t);
			break;
		}
		case IRTK_VOXEL_FLOAT: {
			threshold(*(dynamic_cast<irtkGenericImage<float> *>(image)), t);
			break;
		}
		case IRTK_VOXEL_DOUBLE: {
			threshold(*(dynamic_cast<irtkGenericImage<double> *>(image)), t);
			break;
		}
		default:
			cerr << "Voxel type not supported" << endl;
			exit(1);
	}

	//threshold(image, t);

	// Write image
	image->Write(argv[2]);

	return 0;
}
