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

#ifdef HAS_CONTRIB

#include <irtkEuclideanDistanceTransform.h>

char *input_name = NULL, *output_name = NULL;

void usage()
{
  cerr << "Usage: dmap [in] [out] <-3D/-2D>" << endl;
  exit(1);
}

int main(int argc, char **argv)
{
  int x, y, z, ok;
  irtkRealImage input, inputA, inputB, outputA, outputB;

  if (argc < 3) {
    usage();
  }

  input_name  = argv[1];
  argc--;
  argv++;
  output_name = argv[1];
  argc--;
  argv++;

  // Default mode
  irtkEuclideanDistanceTransform<irtkRealPixel> *edt = NULL;
  while (argc > 1) {
    ok = False;
    if ((ok == False) && (strcmp(argv[1], "-2D") == 0)) {
      argc--;
      argv++;
      edt = new irtkEuclideanDistanceTransform<irtkRealPixel>
      (irtkEuclideanDistanceTransform<irtkRealPixel>::irtkDistanceTransform2D);
      ok = True;
    }
    if ((ok == False) && (strcmp(argv[1], "-3D") == 0)) {
      argc--;
      argv++;
      edt = new irtkEuclideanDistanceTransform<irtkRealPixel>
      (irtkEuclideanDistanceTransform<irtkRealPixel>::irtkDistanceTransform3D);
      ok = True;
    }
    if (ok == False) {
      cerr << "Can't parse argument " << argv[1] << endl;
      usage();
    }
  }

  if (edt == NULL) {
    edt = new irtkEuclideanDistanceTransform<irtkRealPixel>
    (irtkEuclideanDistanceTransform<irtkRealPixel>::irtkDistanceTransform3D);
  }

  // Read input
  input.Read(input_name);

  // Threshold image
  inputA = input;
  inputB = input;
  for (z = 0; z < input.GetZ(); z++) {
    for (y = 0; y < input.GetY(); y++) {
      for (x = 0; x < input.GetX(); x++) {
        if (input(x, y, z) > 0.5) {
          inputA(x, y, z) = 1;
          inputB(x, y, z) = 0;
        } else {
          inputA(x, y, z) = 0;
          inputB(x, y, z) = 1;
        }
      }
    }
  }

  // Calculate EDT
  cout << "Doing outside DT" << endl;
  edt->SetInput (& inputA);
  edt->SetOutput(&outputA);
  edt->Run();
  cout << "Doing inside DT" << endl;
  edt->SetInput (& inputB);
  edt->SetOutput(&outputB);
  edt->Run();

  for (z = 0; z < input.GetZ(); z++) {
    for (y = 0; y < input.GetY(); y++) {
      for (x = 0; x < input.GetX(); x++) {
        outputA(x, y, z)  = sqrt(outputA(x, y, z)) - sqrt(outputB(x, y, z));
      }
    }
  }

  // Write image
  outputA.Write(output_name);

  return 0;
}

#else

int main(int argc, char **argv)
{
  cerr << "Needs to be compiled with HAS_CONTRIB" << endl;
  return 0;
}

#endif
