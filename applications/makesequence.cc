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
  cerr << "Usage: makesequence [input 1 ... input n] [output] <options>\n" << endl;
  exit(1);
}

int main(int argc, char **argv)
{
  int i, x, y, z, t;

  // Determine how many volumes we have
  t = argc-2;

  if (t < 1) usage();

  cout << "Making sequence from " << t << " volumes" << endl;

  irtkGreyImage* input = new irtkGreyImage[t];

  // Read first image
  cout << "Reading " << argv[1] << endl;
  input[0].Read(argv[1]);

  // Read remaining images
  for (i = 1; i < t; i++) {

    cout << "Reading " << argv[i+1] << endl;
    input[i].Read(argv[i+1]);

    if (!(input[0].GetImageAttributes() == input[i].GetImageAttributes())) {
      cerr << "Mismatch of volume geometry" << endl;
      exit(1);
    }
  }

  irtkGreyImage output(input[0].GetImageAttributes());

  cout << "Inserting volumes into sequence" << endl;
  for (i = 0; i < t; i++) {
    cout << "Volume " << i+1 << " ..." << endl;
    for (z = 0; z < output.GetZ(); z++) {
      for (y = 0; y < output.GetY(); y++) {
        for (x = 0; x < output.GetX(); x++) {
          output(x, y, z, i) = input[i](x, y, z);
        }
      }
    }
  }

  // Write image
  cout << "Writing sequence to " << argv[t+1] << endl;
  output.Write(argv[t+1]);

  delete[] input;
}

