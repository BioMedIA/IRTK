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

#ifdef USE_VXL
// need to include vxl here
#else
#include <nr.h>
#endif

void usage()
{
  cerr << "Usage: makevolume [input 1 ... input n] [out] <options>\n" << endl;
  exit(1);
}

int main(int argc, char **argv)
{
  int i, j, k, l, x1, y1, x2, y2, z, t1, t2;
  double xsize1, ysize1, zsize1, tsize1, xsize2, ysize2, zsize2, tsize2, zsize, a, b, c;
  double xaxis1[3], yaxis1[3], zaxis1[3], xaxis2[3], yaxis2[3], zaxis2[3];
  irtkPoint origin1, origin2;

  // Determine how many slices we have
  z = argc-2;

  if (z < 1) usage();

  cout << "Making volume from " << z << " slices" << endl;

  long unsigned* index = new unsigned long[z];
  float* distance = new float[z];
  irtkGreyImage* input = new irtkGreyImage[z];

  // Read first image
  cout << "Reading " << argv[1] << endl;
  input[0].Read(argv[1]);
  input[0].GetPixelSize(&xsize1, &ysize1, &zsize1, &tsize1);
  input[0].GetOrientation(xaxis1, yaxis1, zaxis1);
  origin1 = input[0].GetOrigin();
  x1 = input[0].GetX();
  y1 = input[0].GetY();
  t1 = input[0].GetT();

  // Calculate distance
  distance[0] = 0;
  cout << "Distance is " << distance[0] << endl;

  // Read remaining images
  for (i = 1; i < z; i++) {

    cout << "Reading " << argv[i+1] << endl;
    input[i].Read(argv[i+1]);
    input[i].GetPixelSize(&xsize2, &ysize2, &zsize2, &tsize2);
    input[i].GetOrientation(xaxis2, yaxis2, zaxis2);
    x2 = input[i].GetX();
    y2 = input[i].GetY();
    t2 = input[i].GetT();
    origin2 = input[i].GetOrigin();

    // Calculate distance (was: zaxis)
    distance[i] = (origin1._x - origin2._x) * -zaxis1[0] +
                  (origin1._y - origin2._y) * -zaxis1[1] +
                  (origin1._z - origin2._z) * -zaxis1[2];
    cout << "Distance is " << distance[i] << endl;

    if (i > 1) {
      if ((x1 != x2) || (y1 != y2) || (t1 != t2)) {
        cerr << "Image dimensions are different" << endl;
        exit(1);
      }
      if ((xsize1 != xsize2) || (ysize1 != ysize2) || (tsize1 != tsize2)) {
        cerr << "Voxel sizes are different" << endl;
        exit(1);
      }
      a = xaxis1[0] * xaxis2[0] + xaxis1[1] * xaxis2[1] + xaxis1[2] * xaxis2[2] - 1;
      b = yaxis1[0] * yaxis2[0] + yaxis1[1] * yaxis2[1] + yaxis1[2] * yaxis2[2] - 1;
      c = zaxis1[0] * zaxis2[0] + zaxis1[1] * zaxis2[1] + zaxis1[2] * zaxis2[2] - 1;

      // Added zaxis check (jas)
      if ((fabs(a) > 0.001) || (fabs(b) > 0.001) || (fabs(c) > 0.001)) {
        cerr << "Image orientations are different " << a << " " << b << " " << c << endl;
        exit(1);
      }
    }
  }
#ifdef USE_VXL
  cerr << "Not implemented in the VXL library." << endl;
#else
  if (fabs(distance[0] - distance[1]) > 0) {
    indexx(z, distance-1, index-1);
  } else {
    for (i = 0; i < z; i++) {
      index[i] = i+1;
    }
  }
#endif

  if ((z > 2) && (fabs(distance[index[1]-1] - distance[index[2]-1]) > 0)) {
    // True slice thickness is distance no. 2
    zsize = fabs(distance[index[1]-1] - distance[index[2]-1]);
  } else {
    zsize = zsize1;
  }

  irtkGreyImage output(x1, y1, z, t1);
  output.PutPixelSize(xsize1, ysize1, zsize, tsize1);
  output.PutOrientation(xaxis1, yaxis1, zaxis1); // Added zaxis1 - jas
  output.PutOrigin(input[index[0]-1].GetOrigin() + (input[index[z-1]-1].GetOrigin() - input[index[0]-1].GetOrigin()) / 2.0);

  cout << "Inserting slices into volume" << endl;
  for (k = 0; k < z; k++) {
    cout << "Slice " << index[k]-1 << " ..." << endl;
    for (l = 0; l < t1; l++) {
      for (j = 0; j < y1; j++) {
        for (i = 0; i < x1; i++) {
          output(i, j, k, l) = input[index[k]-1].Get(i, j, 0, l);
        }
      }
    }
  }

  // Write image
  cout << "Writing volume to " << argv[z+1] << endl;
  output.Write(argv[z+1]);

  delete[] index;
  delete[] distance;
  delete[] input;
}

