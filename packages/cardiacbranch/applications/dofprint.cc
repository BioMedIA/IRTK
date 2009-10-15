/*=========================================================================

  Library   : Image Registration Toolkit (IRTK)
  Module    : $Id$
  Copyright : Imperial College, Department of Computing
              Visual Information Processing (VIP), 2008 onwards
  Date      : $Date$
  Version   : $Revision$
  Changes   : $Author$

=========================================================================*/

#include <irtkRegistration.h>

int main(int argc, char **argv)
{
  irtkTransformation *transformation;

  // Check arguments
  if (argc != 2) {
    cerr << "Usage: dofprint [doffile]" << endl;
    exit(1);
  }

  transformation = irtkTransformation::New(argv[1]);
  transformation->Print();
}
