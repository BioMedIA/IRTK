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

char *input_name = NULL, *seg_name = NULL;
char *out_name = NULL;

void usage()
{
	cerr << "Usage: cardiaccrop [input] [segmentation] [output]" << endl;
	cerr << "-boarder [size] reserved size in x y direction"    << endl;
    cerr << "this app crops the cardiac image using the given segmentation"    << endl;
	exit(1);
}

int main( int argc, char** argv )
{
	int ok,i,j,k,minx,maxx,miny,maxy,size;
	// Check command line
	if (argc < 4) {
		usage();
	}

     size = 10;

	// Parse source and target images
	input_name = argv[1];
	argc--;
	argv++;
     seg_name = argv[1];
     argc--;
     argv++;
	out_name = argv[1];
	argc--;
	argv++;

	// Parse remaining parameters
	while (argc > 1) {
		ok = false;
		if ((ok == false) && (strcmp(argv[1], "-boarder") == 0)) {
			argc--;
			argv++;
			ok = true;
			size = atoi(argv[1]);
			argc--;
			argv++;
		}
		if (ok == false) {
			cerr << "Can not parse argument " << argv[1] << endl;
			usage();
		}
	}

	irtkGreyImage input,segmentation,output;

    input.Read(input_name);
    segmentation.Read(seg_name);
    maxx = 0; minx = input.GetX();
    maxy = 0; miny = input.GetY();

    for ( k = 0; k < input.GetZ(); k++){
        for ( j = 0; j< input.GetY(); j++){
            for ( i = 0; i<input.GetX(); i++){
                if(segmentation.GetAsDouble(i,j,k) > 0){
                    if(i > maxx)
                        maxx = i;
                    if(i < minx)
                        minx = i;
                    if(j > maxy)
                        maxy = j;
                    if(j < miny)
                        miny = j;
                }
            }
        }
    }

    minx -= size;
    if(minx < 0) 
        minx = 0;
    maxx += size;
    if(maxx > input.GetX() - 1) 
        minx = input.GetX() - 1;
    miny -= size;
    if(miny < 0) 
        miny = 0;
    maxy -= size;
    if(maxy > input.GetY() - 1) 
        miny = input.GetY() - 1;

    output = input.GetRegion(minx,miny,0,maxx+1,maxy+1,input.GetZ());

	output.Write(out_name);
}
