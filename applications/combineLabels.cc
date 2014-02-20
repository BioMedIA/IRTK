/*=========================================================================

  Library   : Image Registration Toolkit (IRTK)
  Module    : $Id$
  Copyright : Imperial College, Department of Computing
              Visual Information Processing (VIP), 2008 onwards
  Date      : $Date$
  Version   : $Revision$
  Changes   : $Author$

=========================================================================*/

// super-UGLY, but avoid wrecking IRTK build process
#include "src/combineLabels_core.cc"
#include <stdexcept>

char *  output_name = NULL;
char ** input_names = NULL;
char *  mask_name   = NULL;

void usage()
{
  cerr << " " << endl;
  cerr << " Usage: combineLabels [output] [N] [input1] .. [inputN] <-options>" << endl;
  cerr << " " << endl;
  cerr << " Combine [N] label images into a consensus labelling using the vote rule." << endl;
  cerr << " Labels for voxels where the vote is tied for are decided randomly." << endl;
  cerr << " " << endl;
  cerr << " The images [input1] .. [inputN] are assumed to contain labellings on the" << endl;
  cerr << " same voxel grid.  For each voxel, the modal label from the images is used " << endl;
  cerr << " to generate a label in the output image." << endl;
  cerr << " " << endl;
  cerr << " Options:" << endl;
  cerr << " -u [filename]     Write out a volume showing the unanimous voxels." << endl;
  cerr << " -pad [value]      Padding value, default = -1." << endl;
  cerr << " -seed [value]     Setting the same value over several runs will generate the same results (sets the seed of the pseudorandom number generator), default = current time (non-reproducible behaviour)." << endl;
  cerr << " " << endl;
  exit(1);
}

int main(int argc, char **argv)
{
  int numberOfClassifiers, i, v, ok, voxels, contendedVoxelCount, contendedVoxelIndex, equivocalCount;
  irtkGreyPixel *pIn, *pIn_0, *pOut;
  irtkGreyImage input, output;
  irtkGreyImage *input_0;
  irtkByteImage mask;
  irtkBytePixel *pMask;
  int writeMask = false;
  int pad = -1;
  int inSeed = time(NULL);

  // Check command line
  if (argc < 4){
    usage();
  }

  output_name = argv[1];
  argc--;
  argv++;

  numberOfClassifiers = atoi(argv[1]);
  argc--;
  argv++;

  input_names = new char *[numberOfClassifiers];
  for (i = 0; i < numberOfClassifiers; i++){
    input_names[i] = argv[1];
    argc--;
    argv++;
  }

  // TODO change to getopt
  while (argc > 1){
    ok = false;
    if ((ok == false) && (strcmp(argv[1], "-u") == 0)){
      argc--;      argv++;
      mask_name = argv[1];
      writeMask = true;
      argc--;      argv++;
      ok = true;
    }
    if ((ok == false) && (strcmp(argv[1], "-pad") == 0)){
      argc--;      argv++;
      pad = atoi(argv[1]);
      argc--;      argv++;
      ok = true;
    }
    if ((ok == false) && (strcmp(argv[1], "-seed") == 0)){
      argc--;      argv++;
      inSeed = atoi(argv[1]);
      argc--;      argv++;
      ok = true;
    }

    if (ok == false){
      cerr << "Can not parse argument " << argv[1] << endl;
      usage();
    }
  }

  // Mask has a value of 1 for all uncontended voxels.
  mask.Read(input_names[0]);
  voxels = mask.GetNumberOfVoxels();

  // Reset mask so that all voxels are marked as unanimous.
  pMask = mask.GetPointerToVoxels();
  for (i = 0; i < voxels; ++i){
    *pMask = 1;
    ++pMask;
  }

  input_0 = new irtkGreyImage(input_names[0]);

  for (i = 1; i < numberOfClassifiers; ++i){

    input.Read(input_names[i]);

    pIn   = input.GetPointerToVoxels();
    pIn_0 = input_0->GetPointerToVoxels();
    pMask = mask.GetPointerToVoxels();

      for (v = 0; v < voxels; ++v){

        if (*pIn != *pIn_0 && (*pIn > pad || *pIn_0 > pad) ){
          *pMask = 0;
        }
        ++pIn;
        ++pIn_0;
        ++pMask;
      }
  }

  // Do we need to write out the `unanimask'?
  if (writeMask == true && mask_name != NULL){
    mask.Write(mask_name);
  }


  // free some memory.
  delete input_0;

  contendedVoxelCount = 0;
  pMask = mask.GetPointerToVoxels();

  for (i = 0; i < voxels; ++i){
    if (*pMask == 0){
      ++contendedVoxelCount;
    }
    ++pMask;
  }

  cout << "Number of contended voxels = " << contendedVoxelCount << " = ";
  cout << 100.0 * contendedVoxelCount / ((double) voxels) << "%" << endl;

  countMap *counts = new countMap[contendedVoxelCount];

  for (i = 0; i < numberOfClassifiers; ++i){

    input.Read(input_names[i]);

    pIn = input.GetPointerToVoxels();
    pMask = mask.GetPointerToVoxels();

    contendedVoxelIndex = 0;
    for (v = 0; v < voxels; ++v){

      //Is the current voxel contended?
      if (*pMask == 0){
        if (*pIn > pad){
          ++counts[contendedVoxelIndex][*pIn];
        }
        ++contendedVoxelIndex;
      }
      ++pIn;
      ++pMask;
    }
  }

  // Initialise output image using first input image, all inputs assumed
  // same size etc..
  output.Read(input_names[0]);
  pOut  = output.GetPointerToVoxels();
  pMask = mask.GetPointerToVoxels();

  contendedVoxelIndex = 0;
  for (v = 0; v < voxels; ++v){

    //Is the current voxel contentded?
    if (*pMask == 0){
      *pOut = getMostPopular(counts[contendedVoxelIndex]);
      ++contendedVoxelIndex;
    }

    ++pOut;
    ++pMask;
  }

  // Get ready for random stuff.
  // check argument -seed to get a reproducible output
  // good quality PRNG
  boost::mt19937 rng;   
  rng.seed(inSeed);

  contendedVoxelIndex = 0;
  equivocalCount = 0;

  pOut  = output.GetPointerToVoxels();
  pMask = mask.GetPointerToVoxels();
  for (v = 0; v < voxels; ++v){

    //Is the current voxel contended?
    if (*pMask == 0){
      if (isEquivocal(counts[contendedVoxelIndex])){
        ++equivocalCount;
	try {
	  *pOut = decideOnTie(counts[contendedVoxelIndex], rng);
	}
	catch (const std::invalid_argument &e) { std::cerr << e.what() << std::endl; }
      }
      ++contendedVoxelIndex;
    }
    ++pOut;
    ++pMask;
  }

  cout << "Number of equivocal voxels = " << equivocalCount << " = ";
  cout << 100.0 * equivocalCount / ((double) voxels) << "%" << endl;

  output.Write(output_name);

  delete [] counts;
}
