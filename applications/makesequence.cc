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
  cerr << endl;
  cerr << "Options:" << endl;
  cerr << "  -dt <float>      Temporal resolution of sequence." << endl;
  cerr << "  -torigin <float> Temporal origin of sequence." << endl;
  exit(1);
}

int main(int argc, char **argv)
{
  int i, x, y, z, t;
  irtkImageAttributes ipt0_at;     //add (1/2)
  irtkImageAttributes iptI_at;     //add (1/2)

  // Parse options
  double dt      = 1;
  double torigin = 0;

  for (i = argc - 2; i > 0; i -= 2) {
    if (argv[i][0] != '-') break;
    if (strcmp(argv[i], "-dt") == 0) {
      dt = atof(argv[i+1]);
    } else if (strcmp(argv[i], "-torigin") == 0) {
      torigin = atof(argv[i+1]);
    } else break;
    argc -= 2;
  }

  // Determine how many volumes we have
  t = argc-2;

  if (t < 1) usage();

  cout << "Making sequence from " << t << " volumes" << endl;

  irtkImage ** input = new irtkImage*[t];

  // Read first image
  cout << "Reading " << argv[1] << endl;
  input[0] = irtkImage::New(argv[1]);
  ipt0_at = input[0]->GetImageAttributes();

  if (ipt0_at._t > 1) {
    cerr << "Input images may not be four-dimensional" << endl;
    exit(1);
  }

  // Read remaining images
  for (i = 1; i < t; i++) {

    cout << "Reading " << argv[i+1] << endl;
    input[i] = irtkImage::New(argv[i+1]);

    iptI_at = input[i]->GetImageAttributes();
    if (iptI_at._t > 1) {
      cerr << "Input images may not be four-dimensional" << endl;
      exit(1);
    }
    if ((ipt0_at._x!=iptI_at._x)||(ipt0_at._y!=iptI_at._y)||(ipt0_at._z!=iptI_at._z)) {
      cerr << "Mismatch of volume geometry" << endl;
      exit(1);
    }
  }

  irtkImage          *output              = NULL;
  irtkImageAttributes OutputSeqAttributes = ipt0_at;
  OutputSeqAttributes._t       = t;
  OutputSeqAttributes._dt      = dt;
  OutputSeqAttributes._torigin = torigin;

  // Convert image
  switch (input[0]->GetScalarType()) {
    case IRTK_VOXEL_CHAR: {
        output = new irtkGenericImage<char> (OutputSeqAttributes);
      }
      break;
    case IRTK_VOXEL_UNSIGNED_CHAR: {
        output = new irtkGenericImage<unsigned char> (OutputSeqAttributes);
      }
      break;
    case IRTK_VOXEL_SHORT: {
        output = new irtkGenericImage<short> (OutputSeqAttributes);
      }
      break;
    case IRTK_VOXEL_UNSIGNED_SHORT: {
        output = new irtkGenericImage<unsigned short> (OutputSeqAttributes);
      }
      break;
    case IRTK_VOXEL_FLOAT: {
        output = new irtkGenericImage<float> (OutputSeqAttributes);
        break;
      }
    case IRTK_VOXEL_DOUBLE: {
        output = new irtkGenericImage<double> (OutputSeqAttributes);
        break;
      }
    default:
      cerr << "Unknown voxel type for output format" << endl;
      exit(1);
  }

  cout << "Inserting volumes into sequence" << endl;
  for (i = 0; i < t; i++) {
    cout << "Volume " << i+1 << " ..." << endl;
    for (z = 0; z < output->GetZ(); z++) {
      for (y = 0; y < output->GetY(); y++) {
        for (x = 0; x < output->GetX(); x++) {
          output->PutAsDouble(x, y, z, i, input[i]->GetAsDouble(x, y, z));
        }
      }
    }
  }

  // Write image
  cout << "Writing sequence to " << argv[t+1] << endl;
  output->Write(argv[t+1]);

  delete[] input;
}

