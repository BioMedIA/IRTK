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

char *dofin_name  = NULL;
char *dofout_name = NULL;

void usage()
{
  cerr << "Usage: dofinvert [dofin] [dofout]" << endl;
  exit(1);
}

int main(int argc, char **argv)
{
  // Create transformation
  irtkAffineTransformation *transformation = new irtkAffineTransformation;

  // Check arguments
  if (argc != 3) usage();

  // Parse arguments
  dofin_name = argv[1];
  argc--;
  argv++;
  dofout_name = argv[1];
  argc--;
  argv++;

  // Read transform
  transformation->irtkTransformation::Read(dofin_name);
  transformation->Print();

  // Invert transformation
  cout << "Inverting transformation ..." << endl;
  transformation->Invert();
  transformation->UpdateParameter();

  // Write transform
  transformation->irtkTransformation::Write(dofout_name);
  transformation->Print();
}
