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
  cerr << "Usage: ffdcompose [T1] [T2] [mffd_out]\n" << endl;
  cerr << "ffdcompose computes the composition T of two MFFDs T_1 and T_2 such that" << endl;
  cerr << "T(x) = T_2 o T_1(x) \n" << endl;
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

	// Convert second transformation
  irtkMultiLevelFreeFormTransformation *mffd2 = dynamic_cast<irtkMultiLevelFreeFormTransformation *>(t2);

  if (mffd2 == NULL){
  	cerr << "Transformation T_2 is not of type irtkMultiLevelFreeFormTransformation or irtkFluidFreeFormTransformation" << endl;
  	exit(1);
  }

  // If a transformation has more than one level, it must be a Fluid Free Form Transformation.
  if ((mffd1->NumberOfLevels() > 1) && (strcmp(mffd1->NameOfClass(), "irtkMultiLevelFreeFormTransformation") == 0)){
  	cerr << "Transformation T_1 has more than one level and is a irtkMultiLevelFreeFormTransformation " << endl;
  	cerr << "(should be irtkFluidFreeFormTransformation)" << endl;
  	exit(1);
  }

  if ((mffd2->NumberOfLevels() > 1) && (strcmp(mffd2->NameOfClass(), "irtkMultiLevelFreeFormTransformation") == 0)){
  	cerr << "Transformation T_2 has more than one level and is a irtkMultiLevelFreeFormTransformation " << endl;
  	cerr << "(should be irtkFluidFreeFormTransformation)" << endl;
  	exit(1);
  }

  // Details of the first FFD in T2.
  irtkBSplineFreeFormTransformation3D *ffd =
  		dynamic_cast<irtkBSplineFreeFormTransformation3D *> (mffd2->GetLocalTransformation(0));
  if (ffd == NULL){
  	cerr << "Input FFDs should be of type irtkBSplineFreeFormTransformation3D" << endl;
  	exit(1);
  }

	// Create fluid free-form deformation
	irtkFluidFreeFormTransformation *t = new irtkFluidFreeFormTransformation();

	mffd2->MergeGlobalIntoLocalDisplacement();

	// Push all of the FFDs in T1
	for (i = 0; i < mffd1->NumberOfLevels(); i++){
		t->PushLocalTransformation(mffd1->GetLocalTransformation(i));
	}

	// Push all of the FFDs in T2
	for (i = 0; i < mffd2->NumberOfLevels(); i++){
		t->PushLocalTransformation(mffd2->GetLocalTransformation(i));
	}

	// Global matrix from T1.
	t->PutMatrix(mffd1->GetMatrix());

	// Write local transformation
	cout << "Writing T = T_2 o T_1 to " << argv[1] << endl;
	t->irtkTransformation::Write(argv[1]);

}

