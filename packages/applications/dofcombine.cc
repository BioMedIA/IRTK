/*=========================================================================

  Library   : Image Registration Toolkit (IRTK)
  Module    : $Id$
  Copyright : Imperial College, Department of Computing
              Visual Information Processing (VIP), 2008 onwards
  Date      : $Date$
  Version   : $Revision$
  Changes   : $Author$

Copyright (c) IXICO LIMITED
All rights reserved.
See COPYRIGHT for details

=========================================================================*/

#include <irtkRegistration.h>

char *dof1_name = NULL;
char *dof2_name = NULL;
char *dof3_name = NULL;

void usage()
{
  cerr << "Usage: dofcombine [input doffile 1] [input doffile 2] [output doffile] [-invert1] [-invert2]" << endl;
  exit(1);
}

int main(int argc, char **argv)
{
  int ok, invert1, invert2;
  irtkMatrix m1, m2, m3;
  irtkAffineTransformation transformation1, transformation2, transformation3;

  // Check arguments
  if (argc < 4) usage();

  // Parse arguments
  dof1_name = argv[1];
  argc--;
  argv++;
  dof2_name = argv[1];
  argc--;
  argv++;
  dof3_name = argv[1];
  argc--;
  argv++;

  // Default
  invert1 = false;
  invert2 = false;

  while (argc > 1) {
    ok = false;
    if ((ok == false) && (strcmp(argv[1], "-invert1") == 0)) {
      argc--;
      argv++;
      invert1 = true;
      ok = true;
    }
    if ((ok == false) && (strcmp(argv[1], "-invert2") == 0)) {
      argc--;
      argv++;
      invert2 = true;
      ok = true;
    }
    if (ok == false) {
      cerr << "Can not parse argument " << argv[1] << endl;
      usage();
    }
  }

  // Read transformations
  transformation1.irtkTransformation::Read(dof1_name);
  m1 = transformation1.GetMatrix();
  transformation2.irtkTransformation::Read(dof2_name);
  m2 = transformation2.GetMatrix();

  // Invert if necessary
  if (invert1 == true) m1.Invert();
  if (invert2 == true) m2.Invert();

  // Combine
  m3 = m1 * m2;
  transformation3.PutMatrix(m3);
  transformation3.UpdateParameter();

  // Write transformation
  transformation3.irtkTransformation::Write(dof3_name);
  transformation3.Print();
}
