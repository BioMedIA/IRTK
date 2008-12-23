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

#include <irtkTransformation.h>

void usage()
{
  cerr << "Usage: dofimport [old dof file] [new dof file]" << endl;
  exit(1);
}

int main(int argc, char **argv)
{
  irtkTransformation *transformation;

  if (argc != 3) {
    usage();
  }

  // Read transformation
  transformation = irtkTransformation::New(argv[1]);

  // Write transformation
  transformation->Write(argv[2]);
}
