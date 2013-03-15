/*=========================================================================

  Library   : Image Registration Toolkit (IRTK)
  Module    : $Id$
  Copyright : Imperial College, Department of Computing
              Visual Information Processing (VIP), 2008 onwards
  Date      : $Date$
  Version   : $Revision$
  Changes   : $Author$

=========================================================================*/

#ifdef HAS_VTK

#include <irtkRegistration.h>

// Default filenames
char *source_name = NULL, *target_name = NULL, *dofout_name = NULL;

void usage()
{
  cerr << "Usage: prreg [target] [source]\n" << endl;
  cerr << "<-dofout    file>    Final transformation estimate" << endl;
  cerr << "<-xy_only>           Allow only xy translation" << endl;
  exit(1);
}

int main(int argc, char **argv)
{
  int ok, i, xyonly, transonly;
  double error;
  irtkPointSet target, source;

  // Check command line
  if (argc < 3) {
    usage();
  }

  xyonly = false;
  transonly = false;

  // Create registration filter
  irtkPointRigidRegistration *registration = new irtkPointRigidRegistration;

  // Create transformation
  irtkRigidTransformation *transformation = new irtkRigidTransformation;

  // Parse source and target point lists
  target_name = argv[1];
  argc--;
  argv++;
  source_name = argv[1];
  argc--;
  argv++;

  // Read point list
  target.ReadVTK(target_name);

  // Read point list
  source.ReadVTK(source_name);

  // Parse remaining parameters
  while (argc > 1) {
    ok = false;
    if ((ok == false) && (strcmp(argv[1], "-dofout") == 0)) {
      argc--;
      argv++;
      dofout_name = argv[1];
      argc--;
      argv++;
      ok = true;
    }
	if ((ok == false) && (strcmp(argv[1], "-xy_only") == 0)) {
      argc--;
      argv++;
      xyonly = true;
      ok = true;
    }
	if ((ok == false) && (strcmp(argv[1], "-trans_only") == 0)) {
      argc--;
      argv++;
      transonly = true;
      ok = true;
    }
    if (ok == false) {
      cerr << "Can not parse argument " << argv[1] << endl;
      usage();
    }
  }

  // Set input and output for the registration filter
  irtkPointSet tmp1;
  irtkPointSet tmp2;

  double x1,y1,z1,x2,y2,z2;

  x1 = 0; y1 = 0; z1 = 0;
  x2 = 0; y2 = 0; z2 = 0;

  for(i = 0; i < target.Size(); i++){
	  irtkPoint p = target(i);
	  x1 += p._x;
	  y1 += p._y;
	  z1 += p._z;
  }

  x1 = x1 / target.Size();
  y1 = y1 / target.Size();
  z1 = z1 / target.Size();

  for(i = 0; i < source.Size(); i++){
	  irtkPoint p = source(i);
	  x2 += p._x;
	  y2 += p._y;
	  z2 += p._z;
  }

  x2 = x2 / source.Size();
  y2 = y2 / source.Size();
  z2 = z2 / source.Size();

  for(i = 0; i < target.Size(); i++){
	  irtkPoint p = target(i);
	  if(transonly == true){
		  p._x = x1;
		  p._y = y1;
		  p._z = z1;
	  }
	  if(xyonly == true){
		  p._z = 0;
	  }
	  tmp1.Add(p);
  }
  for(i = 0; i < source.Size(); i++){
	  irtkPoint p = source(i);
	  if(transonly == true){
		  p._x = x2;
		  p._y = y2;
		  p._z = z2;
	  }
	  if(xyonly == true){
		  p._z = 0;
	  }
	  tmp2.Add(p);
  }

  registration->SetInput(&tmp1, &tmp2);
  registration->SetOutput((irtkRigidTransformation *)transformation);

  // Run registration filter
  registration->Run();

  // Calculate residual error
  transformation->irtkTransformation::Transform(target);

  error = 0;
  for (i = 0; i < target.Size(); i++) {
    irtkPoint p1 = target(i);
    irtkPoint p2 = source(i);
    error += sqrt(pow(double(p1._x - p2._x), 2.0) +
                  pow(double(p1._y - p2._y), 2.0) +
                  pow(double(p1._z - p2._z), 2.0));
  }
  cout << "RMS is " << error/target.Size() << " mm" << endl;

  // Write the final transformation estimate
  if (dofout_name != NULL) transformation->irtkTransformation::Write(dofout_name);
}

#else

#include <irtkImage.h>

int main( int argc, char *argv[] )
{
  cerr << argv[0] << " needs to be compiled with the VTK library " << endl;
}
#endif
