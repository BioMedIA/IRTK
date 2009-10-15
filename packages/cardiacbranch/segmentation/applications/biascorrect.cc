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
#include <irtkMultiChannelImage.h>

#include <irtkBiasField.h>

char *output;
bool logtransformed =false;

void usage()
{
  cerr << "Usage: biascorrect [image] [bias] [output] <-padding> <-logtransformed>" << endl;
  exit(1);
}

int main(int argc, char **argv)
{
  double x, y, z, bias;
  int i, j, k, ok, padding;

  if (argc < 4) {
    usage();
    exit(1);
  }

  // Input image
  irtkGreyImage image;
  image.Read(argv[1]);
  argc--;
  argv++;

  // Create bias field
  irtkBSplineBiasField *biasfield = new irtkBSplineBiasField;
  biasfield->Read(argv[1]);
  argc--;
  argv++;

  // Output file name
  output = argv[1];
  argc--;
  argv++;

  // Default parameters
  padding    = MIN_GREY;

  // Parse remaining parameters
  while (argc > 1) {
    ok = False;
    if ((ok == False) && (strcmp(argv[1], "-padding") == 0)) {
      argc--;
      argv++;
      padding = atoi(argv[1]);
      argc--;
      argv++;
      ok = True;
    }
    if ((ok == False) && (strcmp(argv[1], "-logtransformed") == 0)) {
      argc--;
      argv++;
      logtransformed = true;
      ok = True;
    }
    if (ok == False) {
      cerr << "Can not parse argument " << argv[1] << endl;
      usage();
    }
  }

  irtkMultiChannelImage mch;
  mch.SetPadding(padding);
  mch.AddImage(image);
  if (!logtransformed) mch.Log(0);



  for (k = 0; k < (mch.GetImage(0)).GetZ(); k++) {
    for (j = 0; j < (mch.GetImage(0)).GetY(); j++) {
      for (i = 0; i < (mch.GetImage(0)).GetX(); i++) {
        x = i;
        y = j;
        z = k;
        (mch.GetImage(0)).ImageToWorld(x, y, z);
        bias = biasfield->Bias(x, y, z);
        if ((mch.GetImage(0))(i, j, k) != padding) {
          (mch.GetImage(0))(i, j, k) = (mch.GetImage(0))(i, j, k) - bias;
        }
      }
    }
  }

  mch.Exp(0);
  (mch.GetImage(0)).Write(output);
}


