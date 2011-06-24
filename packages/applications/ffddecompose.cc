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
#include <irtkImageFunction.h>
#include <irtkTransformation.h>
#include <irtkRegistration.h>

// Default filenames
char *dofOut_name = NULL;
char *dofIn_name  = NULL;
char *dofBase_name = NULL;
char *image_name = NULL;

void usage()
{
  cerr << " Usage: ffddecompose [dofIn] [dofBase] [image] [dofOut] <options>\n" << endl;
  cerr << " T_in = T_base * T_output dofOut's dimension is based on the input image" << endl;
  cerr << " " << endl;
  exit(1);
}

int main(int argc, char **argv)
{
  int i, j, k, xdim, ydim, zdim, numberOfPs, count;
  double dx, dy, dz, x, y, z, xb, yb, zb, xi, yi, zi, xStore, yStore, zStore;

  // Check command line
  if (argc < 5) {
    usage();
  }

  // Parse dofOut_name
  dofIn_name = argv[1];
  argc--;
  argv++;
  dofBase_name = argv[1];
  argc--;
  argv++;
  image_name = argv[1];
  argc--;
  argv++;
  dofOut_name = argv[1];
  argc--;
  argv++;

  irtkTransformation *t_in = irtkTransformation::New(dofIn_name);
  irtkTransformation *t_base = irtkTransformation::New(dofBase_name);

  irtkMultiLevelFreeFormTransformation *mffd_out = new irtkMultiLevelFreeFormTransformation;

  irtkImage *image = irtkImage::New(image_name);

  dx = image->GetXSize()*32.0;
  dy = image->GetYSize()*32.0;
  dz = image->GetZSize()*32.0;

  // Extract FFD and get lattice dimensions
  irtkBSplineFreeFormTransformation *affd_out = new irtkBSplineFreeFormTransformation(*image, dx, dy, dz);

  xdim = image->GetX();
  ydim = image->GetY();
  zdim = image->GetZ();

  numberOfPs = xdim * ydim * zdim;

  // Space to store the control point displacements.
  double *xdata = new double[numberOfPs];
  double *ydata = new double[numberOfPs];
  double *zdata = new double[numberOfPs];
  double *xpos = new double[numberOfPs];
  double *ypos = new double[numberOfPs];
  double *zpos = new double[numberOfPs];
  count = 0;

  // Loop for each control point in the target
  for (k = 0; k < zdim; k++) {
    for (j = 0; j < ydim; j++) {
      for (i = 0; i < xdim; i++) {

        x = i;
        y = j;
        z = k;

        // Transform point from lattice coordinates to target coordinates
        image->ImageToWorld(x, y, z);

        xStore = x;
        yStore = y;
        zStore = z;

        // Find base Transformation and new position
        t_base->Transform(x,y,z);
        // base transformation
        xb = x - xStore;
        yb = y - yStore;
        zb = z - zStore;
        // new position
        xpos[count] = x;
        ypos[count] = y;
        zpos[count] = z;

        // Find input transformation
        x = xStore;
        y = yStore;
        z = zStore;

        // Transform point
        t_in->Transform(x,y,z);
        // input transformation
        xi = x - xStore;
        yi = y - yStore;
        zi = z - zStore;


        // Displacement of the control point.
        xdata[count] = xi - xb;
        ydata[count] = yi - yb;
        zdata[count] = zi - zb;
        count ++;

      }
    }
  }

 
  // Interpolate the ffd and write dof
  for(i = 0; i < 4; i++){
      affd_out->Approximate(xpos,ypos,zpos,xdata, ydata, zdata, numberOfPs);
      affd_out->Subdivide();
  }
  mffd_out->PushLocalTransformation(affd_out);
  mffd_out->irtkTransformation::Write(dofOut_name);

  delete [] xdata;
  delete [] ydata;
  delete [] zdata;
  delete [] xpos;
  delete [] ypos;
  delete [] zpos;
  delete mffd_out;
  delete t_in;
  delete t_base;
}

