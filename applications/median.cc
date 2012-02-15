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
#include <irtkFileToImage.h>

#include <irtkNonLocalMedianFilter.h>

char *input_name = NULL, *output_name = NULL;

void usage()
{
  cerr << "Usage: median [in] [out] [sigma] <options>" << endl;
  cerr << "where <options> are one or more of the following:\n";
  cerr << "\t<-short>           Set data type of output to short integers." << endl;
  cerr << "\t<-float>           Set data type of output to floating point." << endl;
  exit(1);
}

int main(int argc, char **argv)
{
  int ok;
  double sigma;

  if (argc < 4) {
    usage();
  }

  // Parse parameters
  input_name  = argv[1];
  argc--;
  argv++;
  output_name = argv[1];
  argc--;
  argv++;
  sigma = atof(argv[1]);
  argc--;
  argv++;

  // Determine input image data type.
  irtkFileToImage *reader = irtkFileToImage::New(input_name);
  irtkBaseImage *memorycleaner = reader->GetOutput();
  int dataType = reader->GetDataType();

  // Default

  while (argc > 1) {
    ok = false;
    if ((ok == false) && (strcmp(argv[1], "-short") == 0)) {
      argc--;
      argv++;
      // Possible override of input image data type.
      dataType = IRTK_VOXEL_INT;
      ok = true;
    }
    if ((ok == false) && (strcmp(argv[1], "-float") == 0)) {
      argc--;
      argv++;
      // Possible override of input image data type.
      dataType = IRTK_VOXEL_FLOAT;
      ok = true;
    }
    if (ok == false) {
      cerr << "Unknown option: " << argv[1] << endl;
      usage();
    }
  }

  // Blur image
  if (dataType >= IRTK_VOXEL_CHAR && dataType <= IRTK_VOXEL_UNSIGNED_INT)
  {
  	// Read input
  	irtkGreyImage input;
  	input.Read(input_name);

    irtkNonLocalMedianFilter<irtkGreyPixel> median(sigma);
    median.SetInput (&input);
    median.SetOutput(&input);
    median.Run();

  	input.Write(output_name);
  }
  else if (dataType >= IRTK_VOXEL_FLOAT && dataType <= IRTK_VOXEL_DOUBLE)
  {
  	// Read input
  	irtkRealImage input;
  	input.Read(input_name);

    irtkNonLocalMedianFilter<irtkRealPixel> median(sigma);
    median.SetInput (&input);
    median.SetOutput(&input);
    median.Run();
  	// Write image
  	input.Write(output_name);
  }
  else
  {
  	cout << "median: Unknown data type, quitting." << endl;
  	exit(1);
  }

  //clean memory
  delete memorycleaner;
  delete reader;

  return 0;
}
