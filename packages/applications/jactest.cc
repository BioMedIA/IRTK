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

#include <nr.h>
#include <nrutil.h>

#define DELTA 0.001

int main(int, char **)
{
  int i, j, k;
  long seed;
  double tmp, jac1[3], jac2[3];
  irtkPoint p, p1, p2;
  irtkPointSet pset;

  // Create transformation
  irtkAffineTransformation trans;

  // Generate 100 different points
  for (i = 0; i < 100; i++) {
    pset.Add(irtkPoint((ran2(&seed) - 0.5) * 100, (ran2(&seed) - 0.5) * 100, (ran2(&seed) - 0.5) * 100));
  }

  // Loop over 100 random transformations
  for (i = 0; i < 100; i++) {

    // Generate random transformation parameters
    trans.Put(0, (ran2(&seed) - 0.5) * 100);
    trans.Put(1, (ran2(&seed) - 0.5) * 100);
    trans.Put(2, (ran2(&seed) - 0.5) * 100);
    trans.Put(3, ran2(&seed) * 360 - 180);
    trans.Put(4, ran2(&seed) * 360 - 180);
    trans.Put(5, ran2(&seed) * 360 - 180);
    trans.Put(6, (ran2(&seed) + 0.5) * 100);
    trans.Put(7, (ran2(&seed) + 0.5) * 100);
    trans.Put(8, (ran2(&seed) + 0.5) * 100);
    trans.Put(9,  ran2(&seed) * 90 - 45);
    trans.Put(10, ran2(&seed) * 90 - 45);
    trans.Put(11, ran2(&seed) * 90 - 45);

    // Loop over dofs
    for (j = 0; j < trans.NumberOfDOFs(); j++) {

      // Loop over random points
      for (k = 0; k < 100; k++) {
        p = pset(k);

        // Compute Jacobian analytically
        trans.JacobianDOFs(jac1, j, p._x, p._y, p._z);

        // Compute Jacobian numerically
        tmp = trans.Get(j);
        trans.Put(j, tmp + DELTA);
        p1 = p;
        trans.Transform(p1._x, p1._y, p1._z);
        trans.Put(j, tmp - DELTA);
        p2 = p;
        trans.Transform(p2._x, p2._y, p2._z);
        trans.Put(j, tmp);
        jac2[0] = (p1._x - p2._x) / (2.0 * DELTA);
        jac2[1] = (p1._y - p2._y) / (2.0 * DELTA);
        jac2[2] = (p1._z - p2._z) / (2.0 * DELTA);

        // Compute error
        double a, b, c;
        a = sqrt((jac1[0] - jac2[0])*(jac1[0] - jac2[0]));
        b = sqrt((jac1[1] - jac2[1])*(jac1[1] - jac2[1]));
        c = sqrt((jac1[2] - jac2[2])*(jac1[2] - jac2[2]));
        if ((a > 0.0000001) || (b > 0.0000001) || (c > 0.0000001)) {
          cout << i << " " << j << " " << a << " " << b << " " << c << endl;
        }
      }
    }
  }
}
