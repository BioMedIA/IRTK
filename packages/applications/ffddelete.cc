/*=========================================================================

  Library   : Image Registration Toolkit (IRTK)
  Module    : $Id$
  Copyright : Imperial College, Department of Computing
              Visual Information Processing (VIP), 2009 onwards
  Date      : $Date$
  Version   : $Revision$
  Changes   : $Author$

=========================================================================*/

#include <irtkImage.h>

#include <irtkTransformation.h>

char *dofin_name  = NULL;
char *dofout_name = NULL;

irtkFreeFormTransformation *ffd;
irtkMultiLevelFreeFormTransformation *mffd;

void usage()
{
  cerr << "Usage: ffddelete [dofin] [dofout]" << endl;
  cout << endl;
  exit(1);
}

int main(int argc, char **argv)
{
  // Check command line
  if (argc < 3) {
    usage();
  }

  dofin_name  = argv[1];
  argc--;
  argv++;
  dofout_name = argv[1];
  argc--;
  argv++;

  // Read transformation
  cout << "Reading transformation ... "; cout.flush();
  mffd = new irtkMultiLevelFreeFormTransformation;
  mffd->irtkTransformation::Read(dofin_name);
  cout << "done" << endl;

  // Default level is the last level
  if (mffd->NumberOfLevels() > 0){
  	mffd->PopLocalTransformation();
  } else {
  	cerr << "ffddelete: Transformation has 0 ffds" << endl;
		exit(1);
  }

  // Write transformation
  mffd->irtkTransformation::Write(dofout_name);
}
