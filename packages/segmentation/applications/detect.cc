/*=========================================================================
 
  Library   : Image Registration Toolkit (IRTK)
  Module    : $Id: convert.cc 107 2009-12-01 22:03:57Z dr $
  Copyright : Imperial College, Department of Computing
              Visual Information Processing (VIP), 2008 onwards
  Date      : $Date: 2009-12-01 22:03:57 +0000 (‰∫? 01 ÂçÅ‰∫åÊú?2009) $
  Version   : $Revision: 107 $
  Changes   : $Author: dr $
 
=========================================================================*/

#include <irtkSegmentationFunction.h>

char *input_name = NULL, *output_name = NULL;

void usage()
{
  cerr << "Usage: detect [in] [out] <options>\n\n";
  cerr << "where <options can be one or more of the following:\n";
  cerr << "<-scale value>						  Scale output region by value\n";
  cerr << "<-marginal value>					  Detect Object in Low(1) Med(3) High(2) intensity range default 0 detect in all intensity range\n";
  cerr << "<-cascade name>			              Detect object using given cascade haar feature classifier\n";
  exit(1);
}

int main(int argc, char **argv)
{
  double scale;
  int ok,image_type,marginalvalue;
  irtkSegmentationFunction cf;
#ifdef HAS_OPENCV
  CvHaarClassifierCascade *_classifier = NULL;
#endif

  if (argc < 3) {
    usage();
  }

  // Parse filenames
  input_name  = argv[1];
  argc--;
  argv++;
  output_name = argv[1];
  argc--;
  argv++;

  // Parse remaining options
  scale = 1;
  marginalvalue = 0;
  image_type = IRTK_VOXEL_SHORT;

  while (argc > 1) {
    ok = false;
	if (strcmp(argv[1], "-char") == 0) {
      image_type = IRTK_VOXEL_CHAR;
      argc--;
      argv++;
      ok = true;
    } else if (strcmp(argv[1], "-uchar") == 0) {
      image_type = IRTK_VOXEL_UNSIGNED_CHAR;
      argc--;
      argv++;
      ok = true;
    } else if (strcmp(argv[1], "-short") == 0) {
      image_type = IRTK_VOXEL_SHORT;
      argc--;
      argv++;
      ok = true;
    } else if (strcmp(argv[1], "-ushort") == 0) {
      image_type = IRTK_VOXEL_UNSIGNED_SHORT;
      argc--;
      argv++;
      ok = true;
    } else if (strcmp(argv[1], "-float") == 0) {
      image_type = IRTK_VOXEL_FLOAT;
      argc--;
      argv++;
      ok = true;
    } else if (strcmp(argv[1], "-double") == 0) {
      image_type = IRTK_VOXEL_DOUBLE;
      argc--;
      argv++;
      ok = true;
    } else if((strcmp(argv[1], "-scale") == 0)){
		argc--;
		argv++;
		scale = atof(argv[1]);
		argc--;
		argv++;
	} else if((strcmp(argv[1], "-marginal") == 0)){
		argc--;
		argv++;
		marginalvalue = atoi(argv[1]);
		argc--;
		argv++;
	}
#ifdef HAS_OPENCV
	else if ((strcmp(argv[1], "-cascade") == 0)) {
		argc--;
		argv++;
		_classifier = (CvHaarClassifierCascade*)cvLoad( argv[1], 0, 0, 0 );
		cout << "Cascade file position is ... " << argv[1] << endl;
		if( !_classifier )
		{
			cerr<<"ERROR: Could not load classifier: "<<argv[1]<<endl;
		}
		argc--;
		argv++;
		ok = true;
	}
#endif
	else if (!ok) {
      cerr << "Invalid option : " << argv[1] << endl;
      exit(1);
    }
  }
  
  // Read image
  irtkGenericImage<float> image(input_name);

  irtkAffineTransformation interestregion;

  // Detect Object
#ifdef HAS_OPENCV
  if(_classifier){
	  double x,y,z;
	  irtkImageAttributes atr = image.GetImageAttributes();
	  atr._t = 1;
	  irtkRealImage threshold(atr);
	  irtkGreyImage timage(atr);
	  for(k = 0; k < atr._z; k++){
		  for(j = 0; j < atr._y; j++){
			  for(i = 0; i < atr._x; i++){
				  timage.PutAsDouble(i,j,k,image.GetAsDouble(i,j,k));				
			  }
		  }
	  }
	  if(marginalvalue != 0){
		  cf.EvaluateThreshold(&threshold, &timage,0,3);
		  for(k = 0; k < atr._z; k++){
			  for(j = 0; j < atr._y; j++){
				  for(i = 0; i < atr._x; i++){
					  if(threshold.GetAsDouble(i,j,k) == marginalvalue)
						  threshold.PutAsDouble(i,j,k,1);	
					  else
						  threshold.PutAsDouble(i,j,k,0);
				  }
			  }
		  }
		  interestregion = cf.DetectObject(&threshold,&timage,_classifier,timage.GetXSize(),scale,60);
	  }else{
		  interestregion = cf.DetectObject(&timage,_classifier,timage.GetXSize(),scale,60);
	  }
	  irtkPoint ip1,ip2;
	  irtkGreyImage *interest = new irtkGreyImage(timage.GetImageAttributes());
	  cf.GenerateBox(interestregion, interest, ip1, ip2, &timage);
	  image = image.GetRegion(ip1._x,ip1._y,0,ip2._x,ip2._y,image.GetZ());
	  delete interest;
  }
#endif

  // Convert image
  switch (image_type) {
    case IRTK_VOXEL_CHAR: {
        irtkGenericImage<char> output = image;
        output.Write(output_name);
      }
      break;
    case IRTK_VOXEL_UNSIGNED_CHAR: {
        irtkGenericImage<unsigned char> output = image;
        output.Write(output_name);
      }
      break;
    case IRTK_VOXEL_SHORT: {
        irtkGenericImage<short> output = image;
        output.Write(output_name);
      }
      break;
    case IRTK_VOXEL_UNSIGNED_SHORT: {
        irtkGenericImage<unsigned short> output = image;
        output.Write(output_name);
      }
      break;
    case IRTK_VOXEL_FLOAT: {
        irtkGenericImage<float> output = image;
        output.Write(output_name);
        break;
      }
    case IRTK_VOXEL_DOUBLE: {
        irtkGenericImage<double> output = image;
        output.Write(output_name);
        break;
      }
    default:
      cerr << "Unknown voxel type for output format" << endl;
      exit(1);
  }
}
