/*=========================================================================

  Library   : Image Registration Toolkit (IRTK)
  Module    : $Id$
  Copyright : Imperial College, Department of Computing
              Visual Information Processing (VIP), 2008 onwards
  Date      : $Date$
  Version   : $Revision$
  Changes   : $Author$

  =========================================================================*/

#include <irtkTransformation.h>

#define MAX_DOFS 10000
#define EPSILON 0.001

void usage()
{
  cerr << "Usage: dofaverage [dofout] [dofin1...N] <options>\n" << endl;
  cerr << "Average a set of input affine transformations.\n" << endl;
  cerr << endl;

  cerr << "<options> can be one or more of the followings:" << endl;
  cerr << "<-identity weight>         Add an identity transformation before taking the average" << endl;
  cerr << "<-prefix directory>        Directory as the prefix for the file names" << endl;
  cerr << endl;
  exit(1);
}

int main(int argc, char **argv)
{
  // Check command line
  if(argc < 3){
    usage();
  }

  // Parse arguments
  char *dofout = argv[1];
  argc--;
  argv++;

  char **dofin_name = new char *[MAX_DOFS];
  double *dofin_weight = new double[MAX_DOFS];
  double total_weight = 0;

  int dofCount = 0;
  while((argc > 1) && (argv[1][0] != '-')) {
    dofin_name[dofCount] = argv[1];
    dofin_weight[dofCount] = 1;
    total_weight += 1;
    dofCount++;
    argc--;
    argv++;
  }

  if(dofCount > 0){
    cout << "Read " << dofCount << " transformations from command line." << endl;
  }

  // Parse options
  bool ok;
  bool identity = false;
  double identity_weight = 0;
  char *prefix_name = NULL;
  
  while(argc > 1){
    ok = false;
    if((ok == false) && (strcmp(argv[1], "-identity") == 0)){
      identity = true;
      argc--;
      argv++;
      identity_weight = atof(argv[1]);
      argc--;
      argv++;
      ok = true;
    }
    if((ok == false) && (strcmp(argv[1], "-prefix") == 0)){
      argc--;
      argv++;
      prefix_name = argv[1];
      argc--;
      argv++;
      ok = true;
    }
    if(ok == false){
      cerr << "Can not parse argument " << argv[1] << endl;
      usage();
    }
  }

  // Process the input filenames
  for(int i=0; i<dofCount; i++){
    char buffer[255];
    if(prefix_name != NULL){
      sprintf(buffer, "%s/%s", prefix_name, dofin_name[i]);
    }else{
      sprintf(buffer, "%s", dofin_name[i]);
    }
    dofin_name[i] = strdup(buffer);
  }

  // The affine transformation matrices
  irtkMatrix *globMats;
  irtkMatrix globMatAv(4, 4);
  globMatAv.Ident();

  // Allocate memory for input transformations and one extra for the identity assumed from the target to itself
  if(identity == true){
    globMats = new irtkMatrix[dofCount + 1];
  }else{
    globMats = new irtkMatrix[dofCount];
  }

  // Read the affine transformations
  irtkAffineTransformation *transform = NULL;

  for(int i=0; i<dofCount; i++) {
    cout << "Reading " << dofin_name[i] << " for affine averaging" << endl;
    transform = dynamic_cast<irtkAffineTransformation *>(irtkTransformation::New(dofin_name[i]));
    globMats[i].Initialize(4, 4);
    globMats[i] = transform->GetMatrix();
    delete transform;
  }

  // Compute the average transformation
  // Please refer to the following papers for more information
  // [1] Marc Alexa. Linear combination of transformations.
  // [2] Maria Kuklisova-Murgasova et al. A dynamic 4D probabilistic atlas of the developing brain.
  if(identity == true){
    // Put the identity matrix at the end
    globMats[dofCount].Initialize(4, 4);
    globMats[dofCount].Ident();

    dofin_weight[dofCount] = identity_weight;
    total_weight += identity_weight;

    globMatAv = FrechetMean(globMats, dofin_weight, dofCount+1, 20);
  }else{
    globMatAv = FrechetMean(globMats, dofin_weight, dofCount, 20);
  }

  // Write the average transform
  transform = new irtkAffineTransformation;
  transform->PutMatrix(globMatAv);
  transform->irtkTransformation::Write(dofout);

  // Free memory
  delete []dofin_name;
  dofin_name = NULL;
  delete []dofin_weight;
  dofin_weight = NULL;
  delete []globMats;
  globMats = NULL;
  delete transform;
  transform = NULL;

  return 0;
}
