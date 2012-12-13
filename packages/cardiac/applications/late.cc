/*=========================================================================

  Library   : Image Registration Toolkit (IRTK)
  Module    : $Id$
  Copyright : Imperial College, Department of Computing
              Visual Information Processing (VIP), 2008 onwards
  Date      : $Date$
  Version   : $Revision$
  Changes   : $Author$

=========================================================================*/
#include <irtkBep.h>
#ifdef HAS_VTK
char *late_name = NULL, *seg_name = NULL;
char *out_name = NULL;

void usage()
{
	cerr << "Usage: late [lateimage] [segmentation output]" << endl;
	cerr << "-threshold [thresholdname] myocardium mask for segmentation"    << endl;
	exit(1);
}

int main( int argc, char** argv )
{
	int ok,i,j,k;
	irtkBep cf;
	// Check command line
	if (argc < 3) {
		usage();
	}

	// Parse source and target images
	late_name = argv[1];
	argc--;
	argv++;
	out_name = argv[1];
	argc--;
	argv++;

	// Read target image
	cout << "Reading target ... "; cout.flush();
	irtkRealImage **late = new irtkRealImage *[1];
	late[0] = new irtkRealImage(late_name);
	cout << "done" << endl;
	irtkRealImage output(late[0]->GetImageAttributes());

	// Parse remaining parameters
	while (argc > 1) {
		ok = false;
		if ((ok == false) && (strcmp(argv[1], "-threshold") == 0)) {
			argc--;
			argv++;
			ok = true;
			seg_name = argv[1];
			argc--;
			argv++;
		}
		if (ok == false) {
			cerr << "Can not parse argument " << argv[1] << endl;
			usage();
		}
	}

	irtkGreyImage segmentation;
	if(seg_name != NULL){
		segmentation.Read(seg_name);
		for ( k = 0; k < late[0]->GetZ(); k++){
			for ( j = 0; j< late[0]->GetY(); j++){
				for ( i = 0; i<late[0]->GetX(); i++){
					if(segmentation.GetAsDouble(i,j,k) > 0){
					}else
						late[0]->PutAsDouble(i,j,k,-1);
				}
			}
		}
	}
	cf.EvaluateInfarction(&output, late);
	for ( k = 0; k < output.GetZ(); k++){
		for ( j = 0; j< output.GetY(); j++){
			for ( i = 0; i< output.GetX(); i++){
				if(segmentation.GetAsDouble(i,j,k) > 0 || late[0]->GetAsDouble(i,j,k) > 0){
				}else
					output.PutAsDouble(i,j,k,0);
			}
		}
	}
	output.Write(out_name);
}
#else
int main( int argc, char** argv ){
  cerr << "No VTK" <<endl;
  exit(1);
}
#endif
