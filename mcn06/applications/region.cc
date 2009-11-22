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

char *input_name = NULL, *output_name = NULL;

void usage()
{
  cerr << "Usage: region [in] [out] <-Rx1 x1> <-Ry1 y1> <-Rz1 z1> <-Rt1 zt1> "
       << "<-Rx2 x2> <-Ry2 y2> <-Rz2 z2> <-Rt2 t2> " << endl;
  exit(1);
}

int main(int argc, char **argv)
{
  Bool ok;
  irtkGreyImage in, out;
  int x1, x2, y1, t1, y2, z1, z2, t2;

  if (argc < 3) {
    usage();
  }

  input_name  = argv[1];
  argc--;
  argv++;
  output_name = argv[1];
  argc--;
  argv++;

  // Read input
  in.Read(input_name);

  // Default roi
  x1 = 0;
  y1 = 0;
  z1 = 0;
  t1 = 0;
  x2 = in.GetX();
  y2 = in.GetY();
  z2 = in.GetZ();
  t2 = in.GetT();

  while (argc > 1) {
    ok = False;
    if ((ok == False) && (strcmp(argv[1], "-Rx1") == 0)) {
      argc--;
      argv++;
      x1 = atoi(argv[1]);
      argc--;
      argv++;
      ok = True;
    }
    if ((ok == False) && (strcmp(argv[1], "-Rx2") == 0)) {
      argc--;
      argv++;
      x2 = atoi(argv[1]);
      argc--;
      argv++;
      ok = True;
    }
    if ((ok == False) && (strcmp(argv[1], "-Ry1") == 0)) {
      argc--;
      argv++;
      y1 = atoi(argv[1]);
      argc--;
      argv++;
      ok = True;
    }
    if ((ok == False) && (strcmp(argv[1], "-Ry2") == 0)) {
      argc--;
      argv++;
      y2 = atoi(argv[1]);
      argc--;
      argv++;
      ok = True;
    }
    if ((ok == False) && (strcmp(argv[1], "-Rz1") == 0)) {
      argc--;
      argv++;
      z1 = atoi(argv[1]);
      argc--;
      argv++;
      ok = True;
    }
    if ((ok == False) && (strcmp(argv[1], "-Rz2") == 0)) {
      argc--;
      argv++;
      z2 = atoi(argv[1]);
      argc--;
      argv++;
      ok = True;
    }
    if ((ok == False) && (strcmp(argv[1], "-Rt1") == 0)) {
      argc--;
      argv++;
      t1 = atoi(argv[1]);
      argc--;
      argv++;
      ok = True;
    }
    if ((ok == False) && (strcmp(argv[1], "-Rt2") == 0)) {
      argc--;
      argv++;
      t2 = atoi(argv[1]);
      argc--;
      argv++;
      ok = True;
    }
  }

  // Get region
  out = in.GetRegion(x1, y1, z1, t1, x2, y2, z2, t2);

  // Write region
  out.Write(output_name);

  return 0;
}
