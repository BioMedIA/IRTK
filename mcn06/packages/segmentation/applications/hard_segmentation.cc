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

#include <irtkProbabilisticAtlas.h>


char *output_name;

void usage()
{
  cerr << "Usage: hard_segmentation [n] [atlas 1 ... atlas n] [output]" << endl;
  exit(1);
}

int main(int argc, char **argv)
{
  int i, n;

  if (argc < 4) {
    usage();
    exit(1);
  }

  // Number of tissues
  n = atoi(argv[1]);
  argc--;
  argv++;

  // Probabilistic atlas
  irtkRealImage **atlas = new irtkRealImage*[n];


  // Read atlas for each tissue
  for (i = 0; i < n; i++) {
    atlas[i] = new irtkRealImage;
    atlas[i]->Read(argv[1]);
    cerr << "Image " << i <<" = " << argv[1] <<endl;
    argc--;
    argv++;
  }

  // File name for output
  output_name = argv[1];
  argc--;
  argv++;


  irtkProbabilisticAtlas p_atlas;
  p_atlas.AddProbabilityMaps(n, atlas);
  p_atlas.ComputeHardSegmentation();
  p_atlas.WriteHardSegmentation(output_name);

}

