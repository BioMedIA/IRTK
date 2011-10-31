/*=========================================================================

  Library   : Image Registration Toolkit (IRTK)
  Module    : $Id$
  Copyright : Imperial College, Department of Computing
              Visual Information Processing (VIP), 2011 onwards
  Date      : $Date$
  Version   : $Revision$
  Changes   : $Author$

=========================================================================*/

#include <irtkImage.h>

#include <irtkTransformation.h>

void usage()
{
  cerr << "Usage: ffdcompose [ffd1_in] [ffd2_in] [ffd_out]\n" << endl;
  cerr << "ffdcompose computes the composition T of two FFDs T_1 and T_2 such that T = T_1 o T_2\n" << endl;
  exit(1);
}

int main(int argc, char **argv)
{
	int i;

	if (argc != 4){
		usage();
	}

	cout << "Reading T_1 = " << argv[1] << endl;
	irtkTransformation *t1 = irtkTransformation::New(argv[1]);
	argv++;
	argc--;

	cout << "Reading T_2 = " << argv[1] << endl;
  irtkTransformation *t2 = irtkTransformation::New(argv[1]);
	argv++;
	argc--;

	// Print out warning
	cout << "Note the current implementation will ignore the global (affine) transformation component" << endl;

	// Convert first transformation
  irtkMultiLevelFreeFormTransformation *mffd1 = dynamic_cast<irtkMultiLevelFreeFormTransformation *>(t1);
  if (mffd1 == NULL){
  	cerr << "Transformation T_1 is not of type irtkMultiLevelFreeFormTransformation or irtkFluidFreeFormTransformation" << endl;
  	exit(1);
  }

  // Check if first transformation is either a irtkFluidFreeFormTransformation or has only one level
  if ((mffd1->NumberOfLevels() > 1) && (strcmp(mffd1->NameOfClass(), "irtkMultiLevelFreeFormTransformation") == 0)){
  	cerr << "Transformation T_1 has more than one level and is a irtkMultiLevelFreeFormTransformation (should be irtkFluidFreeFormTransformation)" << endl;
  	exit(1);
  }

	// Convert second transformation
  irtkMultiLevelFreeFormTransformation *mffd2 = dynamic_cast<irtkMultiLevelFreeFormTransformation *>(t2);
  if (mffd2 == NULL){
  	cerr << "Transformation T_2 is not of type irtkMultiLevelFreeFormTransformation or irtkFluidFreeFormTransformation" << endl;
  	exit(1);
  }

  // Check if second transformation has only one level
  if (mffd2->NumberOfLevels() > 1){
  	cerr << "Transformation T_2 has more than one level (which is not allowed)" << endl;
  	exit(1);
  }

	// Create fluid free-form deformation
	irtkFluidFreeFormTransformation *t = new irtkFluidFreeFormTransformation();
	for (i = 0; i < mffd1->NumberOfLevels(); i++) t->PushLocalTransformation(mffd1->GetLocalTransformation(i));
	t->PushLocalTransformation(mffd2->PopLocalTransformation());

	// Write local transformation
	cout << "Writing T = T_1 o T_2 = " << argv[1] << endl;
	t->irtkTransformation::Write(argv[1]);
}

