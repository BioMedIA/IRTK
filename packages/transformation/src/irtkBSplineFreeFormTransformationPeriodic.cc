/*=========================================================================

  Library   : Image Registration Toolkit (IRTK)
  Module    : $Id: irtkNonLocalMedianFilter.h 552 2012-02-15 15:18:09Z ws207 $
  Copyright : Imperial College, Department of Computing
              Visual Information Processing (VIP), 2008 onwards
  Date      : $Date: 2012-02-15 15:18:09 +0000 (Wed, 15 Feb 2012) $
  Version   : $Revision: 552 $
  Changes   : $Author: ws207 $

=========================================================================*/

#include <irtkTransformation.h>

#define LUTSIZE (double)(FFDLOOKUPTABLESIZE-1)

double irtkBSplineFreeFormTransformationPeriodic::LookupTable   [FFDLOOKUPTABLESIZE][4];

double irtkBSplineFreeFormTransformationPeriodic::LookupTable_I [FFDLOOKUPTABLESIZE][4];

double irtkBSplineFreeFormTransformationPeriodic::LookupTable_II[FFDLOOKUPTABLESIZE][4];

irtkBSplineFreeFormTransformationPeriodic::irtkBSplineFreeFormTransformationPeriodic()
{
    int i, x, y, z;

    _periodic = true;

    // Initialize control point domain
    _origin._x = 0.5;
    _origin._y = 0.5;
    _origin._z = 0.5;

    // Initialize control point domain (time)
    _tMin = 0;
    _tMax = 1;

    // Initialize x-axis
    _xaxis[0] = 1;
    _xaxis[1] = 0;
    _xaxis[2] = 0;

    // Initialize y-axis
    _yaxis[0] = 0;
    _yaxis[1] = 1;
    _yaxis[2] = 0;

    // Update z-axis
    _zaxis[0] = _xaxis[1]*_yaxis[2] - _xaxis[2]*_yaxis[1];
    _zaxis[1] = _xaxis[2]*_yaxis[0] - _xaxis[0]*_yaxis[2];
    _zaxis[2] = _xaxis[0]*_yaxis[1] - _xaxis[1]*_yaxis[0];

    // Initialize control point dimensions
    _x = 2;
    _y = 2;
    _z = 2;
    _t = 2;

    // Initialize control point spacing
    _dx = 1;
    _dy = 1;
    _dz = 1;
    _dt = 1;

    // Initialize transformation matrix
    _matL2W = irtkMatrix(4, 4);
    _matW2L = irtkMatrix(4, 4);

    // Update transformation matrix
    this->UpdateMatrix();

    // Intialize memory for control point values
    _xdata = this->Allocate(_xdata, _x, _y, _z, _t);
    _ydata = this->Allocate(_ydata, _x, _y, _z, _t);
    _zdata = this->Allocate(_zdata, _x, _y, _z, _t);

    // Initialize memory for control point status
    _status = new _Status[3*_x*_y*_z*_t];
    for (i = 0; i < 3*_x*_y*_z*_t; i++) {
        _status[i] = _Active;
    }

    // set control points for t==_t-1 to passive
    for (z = 0; z < _z; z++) {
        for (y = 0; y < _y; y++) {
            for (x = 0; x < _x; x++) {
                this->PutStatusCP(x, y, z, _t-1, _Passive, _Passive, _Passive);
            }
        }
    }

    // Initialize lookup table
    for (i = 0; i < FFDLOOKUPTABLESIZE; i++) {
        this->LookupTable[i][0]   = this->B0(i/LUTSIZE);
        this->LookupTable[i][1]   = this->B1(i/LUTSIZE);
        this->LookupTable[i][2]   = this->B2(i/LUTSIZE);
        this->LookupTable[i][3]   = this->B3(i/LUTSIZE);
        this->LookupTable_I[i][0] = this->B0_I(i/LUTSIZE);
        this->LookupTable_I[i][1] = this->B1_I(i/LUTSIZE);
        this->LookupTable_I[i][2] = this->B2_I(i/LUTSIZE);
        this->LookupTable_I[i][3] = this->B3_I(i/LUTSIZE);
        this->LookupTable_II[i][0] = this->B0_II(i/LUTSIZE);
        this->LookupTable_II[i][1] = this->B1_II(i/LUTSIZE);
        this->LookupTable_II[i][2] = this->B2_II(i/LUTSIZE);
        this->LookupTable_II[i][3] = this->B3_II(i/LUTSIZE);
    }
}

irtkBSplineFreeFormTransformationPeriodic::irtkBSplineFreeFormTransformationPeriodic(irtkBaseImage &image, double dx, double dy, double dz, double dt)
{
    int i, x ,y, z;
    double x1, y1, z1, x2, y2, z2;

    _periodic = true;

    // Figure out FOV
    x1 = 0;
    y1 = 0;
    z1 = 0;
    x2 = image.GetX()-1;
    y2 = image.GetY()-1;
    z2 = image.GetZ()-1;
    image.ImageToWorld(x1, y1, z1);
    image.ImageToWorld(x2, y2, z2);

    // Initialize control point domain
    _origin._x = (x2 + x1) / 2.0;
    _origin._y = (y2 + y1) / 2.0;
    _origin._z = (z2 + z1) / 2.0;

    // Initialize control point domain
    _tMin = 0;
    _tMax = 1;

    // Initialize x-axis and y-axis
    image.GetOrientation(_xaxis, _yaxis, _zaxis);

    double a = x1 * _xaxis[0] + y1 * _xaxis[1] + z1 * _xaxis[2];
    double b = x1 * _yaxis[0] + y1 * _yaxis[1] + z1 * _yaxis[2];
    double c = x1 * _zaxis[0] + y1 * _zaxis[1] + z1 * _zaxis[2];
    x1 = a;
    y1 = b;
    z1 = c;
    a = x2 * _xaxis[0] + y2 * _xaxis[1] + z2 * _xaxis[2];
    b = x2 * _yaxis[0] + y2 * _yaxis[1] + z2 * _yaxis[2];
    c = x2 * _zaxis[0] + y2 * _zaxis[1] + z2 * _zaxis[2];
    x2 = a;
    y2 = b;
    z2 = c;

    // Initialize control point dimensions
    _x = round((x2 - x1) / dx) + 1;
    _y = round((y2 - y1) / dy) + 1;
    _z = round((z2 - z1) / dz) + 1;
    _t = round((_tMax - _tMin) / dt) + 1;

    // Initialize control point spacing
    _dx = (x2 - x1) / (_x - 1);
    _dy = (y2 - y1) / (_y - 1);
    _dz = (z2 - z1) / (_z - 1);
    _dt = (_tMax - _tMin) / (_t - 1);

    // Initialize transformation matrix
    _matL2W = irtkMatrix(4, 4);
    _matW2L = irtkMatrix(4, 4);

    // Update transformation matrix
    this->UpdateMatrix();

    // Intialize memory for control point values
    _xdata = this->Allocate(_xdata, _x, _y, _z, _t);
    _ydata = this->Allocate(_ydata, _x, _y, _z, _t);
    _zdata = this->Allocate(_zdata, _x, _y, _z, _t);

    // Initialize memory for control point status
    _status = new _Status[3*_x*_y*_z*_t];
    for (i = 0; i < 3*_x*_y*_z*_t; i++) {
        _status[i] = _Active;
    }

    // set control points for t==_t-1 to passive
    for (z = 0; z < _z; z++) {
        for (y = 0; y < _y; y++) {
            for (x = 0; x < _x; x++) {
                this->PutStatusCP(x, y, z, _t-1, _Passive, _Passive, _Passive);
            }
        }
    }

    // Initialize lookup table
    for (i = 0; i < FFDLOOKUPTABLESIZE; i++) {
        this->LookupTable[i][0]   = this->B0(i/LUTSIZE);
        this->LookupTable[i][1]   = this->B1(i/LUTSIZE);
        this->LookupTable[i][2]   = this->B2(i/LUTSIZE);
        this->LookupTable[i][3]   = this->B3(i/LUTSIZE);
        this->LookupTable_I[i][0] = this->B0_I(i/LUTSIZE);
        this->LookupTable_I[i][1] = this->B1_I(i/LUTSIZE);
        this->LookupTable_I[i][2] = this->B2_I(i/LUTSIZE);
        this->LookupTable_I[i][3] = this->B3_I(i/LUTSIZE);
        this->LookupTable_II[i][0] = this->B0_II(i/LUTSIZE);
        this->LookupTable_II[i][1] = this->B1_II(i/LUTSIZE);
        this->LookupTable_II[i][2] = this->B2_II(i/LUTSIZE);
        this->LookupTable_II[i][3] = this->B3_II(i/LUTSIZE);
    }
}

irtkBSplineFreeFormTransformationPeriodic::irtkBSplineFreeFormTransformationPeriodic(double x1, double y1, double z1, double t1,
    double x2, double y2, double z2, double t2,
    double dx, double dy, double dz, double dt,
    double* xaxis, double* yaxis, double* zaxis)
{
    int i, x, y, z;

    _periodic = true;

    // Initialize control point domain
    _origin._x = (x2 + x1) / 2.0;
    _origin._y = (y2 + y1) / 2.0;
    _origin._z = (z2 + z1) / 2.0;

    _tMin = t1;
    _tMax = t2;

    _xaxis[0] = xaxis[0];
    _xaxis[1] = xaxis[1];
    _xaxis[2] = xaxis[2];
    _yaxis[0] = yaxis[0];
    _yaxis[1] = yaxis[1];
    _yaxis[2] = yaxis[2];
    _zaxis[0] = zaxis[0];
    _zaxis[1] = zaxis[1];
    _zaxis[2] = zaxis[2];

    double a = x1 * _xaxis[0] + y1 * _xaxis[1] + z1 * _xaxis[2];
    double b = x1 * _yaxis[0] + y1 * _yaxis[1] + z1 * _yaxis[2];
    double c = x1 * _zaxis[0] + y1 * _zaxis[1] + z1 * _zaxis[2];
    x1 = a;
    y1 = b;
    z1 = c;
    a = x2 * _xaxis[0] + y2 * _xaxis[1] + z2 * _xaxis[2];
    b = x2 * _yaxis[0] + y2 * _yaxis[1] + z2 * _yaxis[2];
    c = x2 * _zaxis[0] + y2 * _zaxis[1] + z2 * _zaxis[2];
    x2 = a;
    y2 = b;
    z2 = c;

    // Initialize control point dimensions
    _x = round((x2 - x1) / dx) + 1;
    _y = round((y2 - y1) / dy) + 1;
    _z = round((z2 - z1) / dz) + 1;
    _t = round((_tMax - _tMin) / dt) + 1;

    // Initialize control point spacing
    _dx = (x2 - x1) / (_x - 1);
    _dy = (y2 - y1) / (_y - 1);
    _dz = (z2 - z1) / (_z - 1);
    _dt = (_tMax - _tMin) / (_t - 1);

    // Initialize transformation matrix
    _matL2W = irtkMatrix(4, 4);
    _matW2L = irtkMatrix(4, 4);

    // Update transformation matrix
    this->UpdateMatrix();

    // Intialize memory for control point values
    _xdata = this->Allocate(_xdata, _x, _y, _z, _t);
    _ydata = this->Allocate(_ydata, _x, _y, _z, _t);
    _zdata = this->Allocate(_zdata, _x, _y, _z, _t);

    // Initialize memory for control point status
    _status = new _Status[3*_x*_y*_z*_t];
    for (i = 0; i < 3*_x*_y*_z*_t; i++) {
        _status[i] = _Active;
    }

    // set control points for t==_t-1 to passive
    for (z = 0; z < _z; z++) {
        for (y = 0; y < _y; y++) {
            for (x = 0; x < _x; x++) {
                this->PutStatusCP(x, y, z, _t-1, _Passive, _Passive, _Passive);
            }
        }
    }

    // Initialize lookup table
    for (i = 0; i < FFDLOOKUPTABLESIZE; i++) {
        this->LookupTable[i][0]   = this->B0(i/LUTSIZE);
        this->LookupTable[i][1]   = this->B1(i/LUTSIZE);
        this->LookupTable[i][2]   = this->B2(i/LUTSIZE);
        this->LookupTable[i][3]   = this->B3(i/LUTSIZE);
        this->LookupTable_I[i][0] = this->B0_I(i/LUTSIZE);
        this->LookupTable_I[i][1] = this->B1_I(i/LUTSIZE);
        this->LookupTable_I[i][2] = this->B2_I(i/LUTSIZE);
        this->LookupTable_I[i][3] = this->B3_I(i/LUTSIZE);
        this->LookupTable_II[i][0] = this->B0_II(i/LUTSIZE);
        this->LookupTable_II[i][1] = this->B1_II(i/LUTSIZE);
        this->LookupTable_II[i][2] = this->B2_II(i/LUTSIZE);
        this->LookupTable_II[i][3] = this->B3_II(i/LUTSIZE);
    }
}

irtkBSplineFreeFormTransformationPeriodic::irtkBSplineFreeFormTransformationPeriodic(const irtkBSplineFreeFormTransformationPeriodic &ffd) : irtkFreeFormTransformation4D(ffd)
{
    int i, j, k, l;

    _periodic = true;

    // Initialize origin
    _origin = ffd._origin;

    // Initialize control point dimensions
    _x = ffd._x;
    _y = ffd._y;
    _z = ffd._z;
    _t = ffd._t;

    // Intialize control point spacing
    _dx = ffd._dx;
    _dy = ffd._dy;
    _dz = ffd._dz;
    _dt = ffd._dt;

    // Intialize time domain
    _tMin = ffd._tMin;
    _tMax = ffd._tMax;

    // Initialize x-axis
    _xaxis[0] = ffd._xaxis[0];
    _xaxis[1] = ffd._xaxis[1];
    _xaxis[2] = ffd._xaxis[2];

    // Initialize y-axis
    _yaxis[0] = ffd._yaxis[0];
    _yaxis[1] = ffd._yaxis[1];
    _yaxis[2] = ffd._yaxis[2];

    // Initialize z-axis
    _zaxis[0] = ffd._zaxis[0];
    _zaxis[1] = ffd._zaxis[1];
    _zaxis[2] = ffd._zaxis[2];

    // Initialize transformation matrix
    _matL2W = irtkMatrix(4, 4);
    _matW2L = irtkMatrix(4, 4);

    // Update transformation matrix
    this->UpdateMatrix();

    // Initialize memory for control point values
    _xdata = this->Allocate(_xdata, _x, _y, _z, _t);
    _ydata = this->Allocate(_ydata, _x, _y, _z, _t);
    _zdata = this->Allocate(_zdata, _x, _y, _z, _t);
    for (i = -4; i < _x+4; i++) {
        for (j = -4; j < _y+4; j++) {
            for (k = -4; k < _z+4; k++) {
                for (l = -4; l < _t+4; l++) {
                    _xdata[l][k][j][i] = ffd._xdata[l][k][j][i];
                    _ydata[l][k][j][i] = ffd._ydata[l][k][j][i];
                    _zdata[l][k][j][i] = ffd._zdata[l][k][j][i];
                }
            }
        }
    }

    // Initialize memory for control point status
    _status = new _Status[3*_x*_y*_z*_t];
    for (i = 0; i < 3*_x*_y*_z*_t; i++) {
        _status[i] = ffd._status[i];
    }

    // Initialize lookup table
    for (i = 0; i < FFDLOOKUPTABLESIZE; i++) {
        this->LookupTable[i][0]   = this->B0(i/LUTSIZE);
        this->LookupTable[i][1]   = this->B1(i/LUTSIZE);
        this->LookupTable[i][2]   = this->B2(i/LUTSIZE);
        this->LookupTable[i][3]   = this->B3(i/LUTSIZE);
        this->LookupTable_I[i][0] = this->B0_I(i/LUTSIZE);
        this->LookupTable_I[i][1] = this->B1_I(i/LUTSIZE);
        this->LookupTable_I[i][2] = this->B2_I(i/LUTSIZE);
        this->LookupTable_I[i][3] = this->B3_I(i/LUTSIZE);
        this->LookupTable_II[i][0] = this->B0_II(i/LUTSIZE);
        this->LookupTable_II[i][1] = this->B1_II(i/LUTSIZE);
        this->LookupTable_II[i][2] = this->B2_II(i/LUTSIZE);
        this->LookupTable_II[i][3] = this->B3_II(i/LUTSIZE);
    }
}

irtkBSplineFreeFormTransformationPeriodic::~irtkBSplineFreeFormTransformationPeriodic()
{
    // Free memory for control points if necessary
    if (_xdata != NULL) _xdata = this->Deallocate(_xdata, _x, _y, _z, _t);
    if (_ydata != NULL) _ydata = this->Deallocate(_ydata, _x, _y, _z, _t);
    if (_zdata != NULL) _zdata = this->Deallocate(_zdata, _x, _y, _z, _t);

    _x = 0;
    _y = 0;
    _z = 0;
}

void irtkBSplineFreeFormTransformationPeriodic::FFD1(double &x, double &y, double &z, double time) const
{
    double s, t, u, v, B_I, B_J, B_K, B_L;
    int a, b, c, d, i, j, k, l, S, T, U, V;
    int tp; // time point used for periodicity
    double timep;

    //  // Check if there is some work to do
    //  if ((x < -2) || (y < -2) || (z < -2) || (time < -2) || (x > _x+1) || (y > _y+1) || (z > _z+1) || (time > _t+1)) {
    //    x = 0;
    //    y = 0;
    //    z = 0;
    //    return;
    //  }
    // Check if there is some work to do
    if ((x < -2) || (y < -2) || (z < -2) || (x > _x+1) || (y > _y+1) || (z > _z+1)) {
        x = 0;
        y = 0;
        z = 0;
        return;
    }
    // periodic model:
    timep = time;
    if (_periodic) {
        while (timep < 0)
            timep += double(_t-1);
        while (timep >= _t-1)
            timep -= double(_t-1);
    } else {
        if ((time < -2) || (time > _t+1)) {
            x = 0;
            y = 0;
            z = 0;
            return;
        }
    }

    // Now calculate the real stuff
    a = (int)floor(x)-1;
    b = (int)floor(y)-1;
    c = (int)floor(z)-1;
    d = (int)floor(timep)-1;

    s = x-(a+1);
    t = y-(b+1);
    u = z-(c+1);
    v = timep-(d+1);

    S = round(LUTSIZE*s);
    T = round(LUTSIZE*t);
    U = round(LUTSIZE*u);
    V = round(LUTSIZE*v);

    // Initialize displacement
    x = 0;
    y = 0;
    z = 0;

    for (l = 0; l < 4; l++) {
        B_L = this->LookupTable[V][l];
        for (k = 0; k < 4; k++) {
            B_K = this->LookupTable[U][k] * B_L;
            for (j = 0; j < 4; j++) {
                B_J = this->LookupTable[T][j] * B_K;
                for (i = 0; i < 4; i++) {
                    B_I = this->LookupTable[S][i] * B_J;

                    // periodic model:
                    tp = ( (d+l)<0 ) 	   ? (d+l+_t-1) : (d+l);
                    tp = ( (d+l)>=_t-1 ) ? (d+l-_t+1) : (d+l);

                    x += B_I * _xdata[tp][c+k][b+j][a+i];
                    y += B_I * _ydata[tp][c+k][b+j][a+i];
                    z += B_I * _zdata[tp][c+k][b+j][a+i];
                }
            }
        }
    }
}

void irtkBSplineFreeFormTransformationPeriodic::FFD2(double &x, double &y, double &z, double time) const
{
    double s, t, u, v, B_I, B_J, B_K, B_L;
    int a, b, c, d, i, j, k, l, S, T, U, V;
    int tp; // time point used for periodicity

    // Now calculate the real stuff
    a = (int)floor(x)-1;
    b = (int)floor(y)-1;
    c = (int)floor(z)-1;
    d = (int)floor(time)-1;

    s = x-(a+1);
    t = y-(b+1);
    u = z-(c+1);
    v = time-(d+1);

    S = round(LUTSIZE*s);
    T = round(LUTSIZE*t);
    U = round(LUTSIZE*u);
    V = round(LUTSIZE*v);

    // Initialize displacement
    x = 0;
    y = 0;
    z = 0;

    for (l = 0; l < 4; l++) {
        B_L = this->LookupTable[V][l];
        for (k = 0; k < 4; k++) {
            B_K = this->LookupTable[U][k] * B_L;
            for (j = 0; j < 4; j++) {
                B_J = this->LookupTable[T][j] * B_K;
                for (i = 0; i < 4; i++) {
                    B_I = this->LookupTable[S][i] * B_J;

                    // periodic model:
                    tp = ( (d+l)<0 ) 	   ? (d+l+_t-1) : (d+l);
                    tp = ( (d+l)>=_t-1 ) ? (d+l-_t+1) : (d+l);

                    x += B_I * _xdata[tp][c+k][b+j][a+i];
                    y += B_I * _ydata[tp][c+k][b+j][a+i];
                    z += B_I * _zdata[tp][c+k][b+j][a+i];
                }
            }
        }
    }
}

double irtkBSplineFreeFormTransformationPeriodic::Approximate(double *x1, double *y1, double *z1, double *t1, double *x2, double *y2, double *z2, int no)
{
    int a, b, c, d, i, j, k, l, I, J, K, L, S, T, U, V, index;
    double s, t, u, v, x, y, z, B_I, B_J, B_K, B_L, basis, basis2, error, phi, norm, time;
    int tp; // time point used for periodicity

    // Allocate memory
    double ****dx = NULL;
    double ****dy = NULL;
    double ****dz = NULL;
    double ****ds = NULL;
    dx = Allocate(dx, _x, _y, _z, _t);
    dy = Allocate(dy, _x, _y, _z, _t);
    dz = Allocate(dz, _x, _y, _z, _t);
    ds = Allocate(ds, _x, _y, _z, _t);

    // Subtract displacements which are approximated by current control points
    for (index = 0; index < no; index++) {
        x = x1[index];
        y = y1[index];
        z = z1[index];
        this->LocalDisplacement(x, y, z, t1[index]);
        x2[index] -= x;
        y2[index] -= y;
        z2[index] -= z;
    }

    // Initialize data structures
    for (l = -2; l < _t+2; l++) {
        for (k = -2; k < _z+2; k++) {
            for (j = -2; j < _y+2; j++) {
                for (i = -2; i < _x+2; i++) {
                    dx[l][k][j][i] = 0;
                    dy[l][k][j][i] = 0;
                    dz[l][k][j][i] = 0;
                    ds[l][k][j][i] = 0;
                }
            }
        }
    }

    // Initial loop: Calculate change of control points
    for (index = 0; index < no; index++) {
        x = x1[index];
        y = y1[index];
        z = z1[index];
        time = t1[index];
        this->WorldToLattice(x, y, z);
        time = this->TimeToLattice(time);

        // periodic model:
        if (_periodic) {
            while (time < 0)
                time += double(_t-1);
            while (time >= _t-1)
                time -= double(_t-1);
        }

        a = (int)floor(x);
        b = (int)floor(y);
        c = (int)floor(z);
        d = (int)floor(time);
        s = x-a;
        t = y-b;
        u = z-c;
        v = time-d;
        S = round(LUTSIZE*s);
        T = round(LUTSIZE*t);
        U = round(LUTSIZE*u);
        V = round(LUTSIZE*v);
        norm = 0;
        for (l = 0; l < 4; l++) {
            B_L = this->LookupTable[V][l];
            for (k = 0; k < 4; k++) {
                B_K = B_L * this->LookupTable[U][k];
                for (j = 0; j < 4; j++) {
                    B_J = B_K * this->LookupTable[T][j];
                    for (i = 0; i < 4; i++) {
                        B_I = B_J * this->LookupTable[S][i];
                        norm += B_I * B_I;
                    }
                }
            }
        }
        for (l = 0; l < 4; l++) {
            B_L = this->LookupTable[V][l];
            L = l + d - 1;

            // periodic model:
            L = ( (L)<0 )     ? (L+_t-1) : (L);
            L = ( (L)>=_t-1 ) ? (L-_t+1) : (L);

            if ((L >= -2) && (L < _t+2)) {
                for (k = 0; k < 4; k++) {
                    B_K = B_L * this->LookupTable[U][k];
                    K = k + c - 1;
                    if ((K >= 0) && (K < _z+2)) {
                        for (j = 0; j < 4; j++) {
                            B_J = B_K * this->LookupTable[T][j];
                            J = j + b - 1;
                            if ((J >= -2) && (J < _y+2)) {
                                for (i = 0; i < 4; i++) {
                                    B_I = B_J * this->LookupTable[S][i];
                                    I = i + a - 1;
                                    if ((I >= -2) && (I < _x+2)) {
                                        basis = B_I / norm;
                                        basis2 = B_I * B_I;
                                        phi = x2[index] * basis;
                                        dx[L][K][J][I] += basis2 * phi;
                                        phi = y2[index] * basis;
                                        dy[L][K][J][I] += basis2 * phi;
                                        phi = z2[index] * basis;
                                        dz[L][K][J][I] += basis2 * phi;
                                        ds[L][K][J][I] += basis2;
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
    }

    // Add displacements which are approximated by current control points
    for (index = 0; index < no; index++) {
        x = x1[index];
        y = y1[index];
        z = z1[index];
        this->LocalDisplacement(x, y, z, t1[index]);
        x2[index] += x;
        y2[index] += y;
        z2[index] += z;
    }

    // Final loop: Calculate new control points
    for (l = -2; l < _t+2; l++) {
        for (k = -2; k < _z+2; k++) {
            for (j = -2; j < _y+2; j++) {
                for (i = -2; i < _x+2; i++) {
                    if (ds[l][k][j][i] > 0) {
                        _xdata[l][k][j][i] += dx[l][k][j][i] / ds[l][k][j][i];
                        _ydata[l][k][j][i] += dy[l][k][j][i] / ds[l][k][j][i];
                        _zdata[l][k][j][i] += dz[l][k][j][i] / ds[l][k][j][i];
                    }
                }
            }
        }
    }

    // Calculate residual error
    error = 0;
    double max = 0;
    for (index = 0; index < no; index++) {
        x = x1[index];
        y = y1[index];
        z = z1[index];
        this->LocalDisplacement(x, y, z, t1[index]);
        x2[index] -= x;
        y2[index] -= y;
        z2[index] -= z;
        // Calculate error
        error += sqrt(x2[index]*x2[index]+y2[index]*y2[index]+z2[index]*z2[index]);
        if (sqrt(x2[index]*x2[index]+y2[index]*y2[index]+z2[index]*z2[index]) > max) max = sqrt(x2[index]*x2[index]+y2[index]*y2[index]+z2[index]*z2[index]);
    }
    error = error / (double)no;
    cout << max << endl;

    // Deallocate memory
    Deallocate(dx, _x, _y, _z, _t);
    Deallocate(dy, _x, _y, _z, _t);
    Deallocate(dz, _x, _y, _z, _t);
    Deallocate(ds, _x, _y, _z, _t);

    // Return error
    return error;
}

void irtkBSplineFreeFormTransformationPeriodic::Jacobian(irtkMatrix &jac, double x, double y, double z, double t)
{
    this->LocalJacobian(jac, x, y, z, t);
}

void irtkBSplineFreeFormTransformationPeriodic::GlobalJacobian(irtkMatrix &jac, double, double, double, double)
{
    // Jacobian matrix is 3 x 3
    jac.Initialize(3, 3);

    // Set matrix to identity
    jac(0, 0) = 1;
    jac(1, 1) = 1;
    jac(2, 2) = 1;
}

void irtkBSplineFreeFormTransformationPeriodic::LocalJacobian(irtkMatrix &jac, double x, double y, double z, double t)
{
    int i, j, k, l;
    int floor_x, floor_y, floor_z, floor_t;
    int I, J, K, L;
    int IND_X, IND_Y, IND_Z, IND_T;

    double frac_x, frac_y, frac_z, frac_t;
    double coeff;
    double B_L, B_K, B_J, B_I, B_K_I, B_J_I, B_I_I;

    // The transformation maps (x, y, z, t) to (Tx, Ty, Tz)
    // Find the partial derivatives of the transformation Tx, Ty and Tz w.r.t x, y and z
    //     dTz/dz dTy/dz dTx/dz dTz/dy dTy/dy dTx/dy dTz/dx dTy/dx dTx/dx
    double z_k=0, y_k=0, x_k=0, z_j=0, y_j=0, x_j=0, z_i=0, y_i=0, x_i=0;

    // Partial derivatives with respect to time are not provided by this function.
    // I.e. none of dTx/dt , dTy/dt or dTz/dt is provided (they should be zero, right?).

    // Convert to lattice coordinates
    this->WorldToLattice(x, y, z);
    t = this->TimeToLattice(t);

    // Compute derivatives
    floor_x = (int)floor(x);
    floor_y = (int)floor(y);
    floor_z = (int)floor(z);
    floor_t = (int)floor(t);

    frac_x = x-floor_x;
    frac_y = y-floor_y;
    frac_z = z-floor_z;
    frac_t = t-floor_t;

    IND_X = round(LUTSIZE*frac_x);
    IND_Y = round(LUTSIZE*frac_y);
    IND_Z = round(LUTSIZE*frac_z);
    IND_T = round(LUTSIZE*frac_t);

    for (l = 0; l < 4; l++){
        L = l + floor_t - 1;

        // periodic model:
        L = ( (L)<0 )     ? (L+_t-1) : (L);
        L = ( (L)>=_t-1 ) ? (L-_t+1) : (L);


        if ((L >= 0) && (L < _t)){
            B_L   = this->LookupTable[IND_T][l];
            // We are only returning the first three columns of the
            // Jacobian so do not need the following commented out bit
            // (the full Jacobian is a 4x3 matrix and we return a 3x3 one)

            // B_L_I = this->LookupTable_I[IND_T][l];

            for (k = 0; k < 4; k++) {
                K = k + floor_z - 1;
                if ((K >= 0) && (K < _z)) {
                    B_K   = this->LookupTable[IND_Z][k];
                    B_K_I = this->LookupTable_I[IND_Z][k];
                    for (j = 0; j < 4; j++) {
                        J = j + floor_y - 1;
                        if ((J >= 0) && (J < _y)) {
                            B_J   = this->LookupTable[IND_Y][j];
                            B_J_I = this->LookupTable_I[IND_Y][j];
                            for (i = 0; i < 4; i++) {
                                I = i + floor_x - 1;
                                if ((I >= 0) && (I < _x)) {
                                    B_I   = this->LookupTable[IND_X][i];
                                    B_I_I = this->LookupTable_I[IND_X][i];
                                    coeff = B_I_I * B_J * B_K * B_L;
                                    x_i += _xdata[L][K][J][I] * coeff;
                                    y_i += _ydata[L][K][J][I] * coeff;
                                    z_i += _zdata[L][K][J][I] * coeff;
                                    coeff = B_I * B_J_I * B_K * B_L;
                                    x_j += _xdata[L][K][J][I] * coeff;
                                    y_j += _ydata[L][K][J][I] * coeff;
                                    z_j += _zdata[L][K][J][I] * coeff;
                                    coeff = B_I * B_J * B_K_I * B_L;
                                    x_k += _xdata[L][K][J][I] * coeff;
                                    y_k += _ydata[L][K][J][I] * coeff;
                                    z_k += _zdata[L][K][J][I] * coeff;
                                    // We are only returning the first three columns of the
                                    // Jacobian so do not need the following
                                    //  coeff = B_I * B_J * B_K * B_L_I;
                                    // 	x_l += _xdata[L][K][J][I] * coeff;
                                    //  y_l += _ydata[L][K][J][I] * coeff;
                                    //  z_l += _zdata[L][K][J][I] * coeff;
                                }
                            }
                        }
                    }
                }
            }

        }
    }

    // First 3 columns of the Jacobian matrix form a 3 x 3 sub-matrix
    jac.Initialize(3, 3);

    // Get deformation derivatives
    jac(0, 0) = x_i;
    jac(1, 0) = y_i;
    jac(2, 0) = z_i;
    jac(0, 1) = x_j;
    jac(1, 1) = y_j;
    jac(2, 1) = z_j;
    jac(0, 2) = x_k;
    jac(1, 2) = y_k;
    jac(2, 2) = z_k;

    // Convert derivatives to world coordinates
    jac = jac * _matW2L(0, 0, 3, 3);
    jac(0, 0) += 1;
    jac(1, 1) += 1;
    jac(2, 2) += 1;
}

double irtkBSplineFreeFormTransformationPeriodic::Bending(double, double, double, double)
{
    cerr << "irtkBSplineFreeFormTransformationPeriodic::Bending: Not implemented yet" << endl;
    exit(1);
}

double irtkBSplineFreeFormTransformationPeriodic::Bending2D(int i, int j, int t)
{
    int I, J;
    double B_J, B_I, B_J_I, B_I_I, B_J_II, B_I_II, v;

    // Derivatives
    double y_jj=0, x_jj=0, y_ii=0, x_ii=0, y_ij=0, x_ij=0;

    // Values of the B-spline basis functions and its derivative (assuming that we compute the bending energy only at the control point location c = i, j, k)
    double b[3]    = {1.0/6.0, 2.0/3.0, 1.0/6.0};
    double b_i[3]  = {-0.5, 0, 0.5};
    double b_ii[3] = {1.0, -2.0, 1.0};

    for (J = j-1; J < j+2; J++) {
        // Get B-spline basis
        B_J    = b   [J-(j-1)];
        B_J_I  = b_i [J-(j-1)];
        B_J_II = b_ii[J-(j-1)];
        for (I = i-1; I < i+2; I++) {
            // Get B-spline basis
            B_I    = b   [I-(i-1)];
            B_I_I  = b_i [I-(i-1)];
            B_I_II = b_ii[I-(i-1)];

            v = B_I_II * B_J;
            y_ii += _ydata[t][0][J][I] * v;
            x_ii += _xdata[t][0][J][I] * v;

            v = B_I * B_J_II;
            y_jj += _ydata[t][0][J][I] * v;
            x_jj += _xdata[t][0][J][I] * v;

            v = B_I_I * B_J_I;
            y_ij += _ydata[t][0][J][I] * v;
            x_ij += _xdata[t][0][J][I] * v;
        }
    }

    // Compute bending
    return (x_ii*x_ii + x_jj*x_jj + y_ii*y_ii + y_jj*y_jj + 2*(x_ij*x_ij + y_ij*y_ij));
}

double irtkBSplineFreeFormTransformationPeriodic::Bending3D(int i, int j, int k, int t)
{
    int I, J, K;
    double B_K, B_J, B_I, B_K_I, B_J_I, B_I_I, B_K_II, B_J_II, B_I_II, v;

    // Derivatives
    double z_kk=0, y_kk=0, x_kk=0, z_jj=0, y_jj=0, x_jj=0, z_ii=0, y_ii=0, x_ii=0;
    double z_ij=0, y_ij=0, x_ij=0, z_ik=0, y_ik=0, x_ik=0, z_jk=0, y_jk=0, x_jk=0;

    // Values of the B-spline basis functions and its derivative (assuming that we compute the bending energy only at the control point location c = i, j, k)
    double b[3]    = {1.0/6.0, 2.0/3.0, 1.0/6.0};
    double b_i[3]  = {-0.5, 0, 0.5};
    double b_ii[3] = {1.0, -2.0, 1.0};

    for (K = k-1; K < k+2; K++) {
        // Get B-spline basis
        B_K    = b   [K-(k-1)];
        B_K_I  = b_i [K-(k-1)];
        B_K_II = b_ii[K-(k-1)];
        for (J = j-1; J < j+2; J++) {
            // Get B-spline basis
            B_J    = b   [J-(j-1)];
            B_J_I  = b_i [J-(j-1)];
            B_J_II = b_ii[J-(j-1)];
            for (I = i-1; I < i+2; I++) {
                // Get B-spline basis
                B_I    = b   [I-(i-1)];
                B_I_I  = b_i [I-(i-1)];
                B_I_II = b_ii[I-(i-1)];

                v = B_I * B_J * B_K_II;
                z_kk += _zdata[t][K][J][I] * v;
                y_kk += _ydata[t][K][J][I] * v;
                x_kk += _xdata[t][K][J][I] * v;

                v = B_I * B_J_II * B_K;
                z_jj += _zdata[t][K][J][I] * v;
                y_jj += _ydata[t][K][J][I] * v;
                x_jj += _xdata[t][K][J][I] * v;

                v = B_I_II * B_J * B_K;
                z_ii += _zdata[t][K][J][I] * v;
                y_ii += _ydata[t][K][J][I] * v;
                x_ii += _xdata[t][K][J][I] * v;

                v = B_I_I * B_J_I * B_K;
                z_ij += _zdata[t][K][J][I] * v;
                y_ij += _ydata[t][K][J][I] * v;
                x_ij += _xdata[t][K][J][I] * v;

                v = B_I_I * B_J * B_K_I;
                z_ik += _zdata[t][K][J][I] * v;
                y_ik += _ydata[t][K][J][I] * v;
                x_ik += _xdata[t][K][J][I] * v;

                v = B_I * B_J_I * B_K_I;
                z_jk += _zdata[t][K][J][I] * v;
                y_jk += _ydata[t][K][J][I] * v;
                x_jk += _xdata[t][K][J][I] * v;
            }
        }
    }

    // Compute bending
    return (x_ii*x_ii + x_jj*x_jj + x_kk*x_kk +
        y_ii*y_ii + y_jj*y_jj + y_kk*y_kk +
        z_ii*z_ii + z_jj*z_jj + z_kk*z_kk +
        2*(x_ij*x_ij + y_ij*y_ij + z_ij*z_ij +
        x_ik*x_ik + y_ik*y_ik + z_ik*z_ik +
        x_jk*x_jk + y_jk*y_jk + z_jk*z_jk));
}

double irtkBSplineFreeFormTransformationPeriodic::Bending()
{
    int i, j, k, t;
    double bending;

    bending = 0;
    if (_z == 1) {
        for (t = 0; t < _t; t++){
            for (j = 0; j < _y; j++) {
                for (i = 0; i < _x; i++) {
                    bending += this->Bending2D(i, j, t);
                }
            }
        }
    } else {
        for (t = 0; t < _t; t++){
            for (k = 0; k < _z; k++) {
                for (j = 0; j < _y; j++) {
                    for (i = 0; i < _x; i++) {
                        bending += this->Bending3D(i, j, k, t);
                    }
                }
            }
        }
    }
    return bending;
}

void irtkBSplineFreeFormTransformationPeriodic::BendingGradient2D(double *gradient)
{
    int I, J, index, i, j, k, t, n = _x*_y*_z*_t;
    double B_J, B_I, B_J_I, B_I_I, B_J_II, B_I_II, v, tmp[3];

    // Derivatives
    double ****xd_jj = NULL;
    xd_jj = this->Allocate(xd_jj, _x, _y, _z, _t);
    double ****yd_jj = NULL;
    yd_jj = this->Allocate(yd_jj, _x, _y, _z, _t);
    double ****zd_jj = NULL;
    zd_jj = this->Allocate(zd_jj, _x, _y, _z, _t);
    double ****xd_ii = NULL;
    xd_ii = this->Allocate(xd_ii, _x, _y, _z, _t);
    double ****yd_ii = NULL;
    yd_ii = this->Allocate(yd_ii, _x, _y, _z, _t);
    double ****zd_ii = NULL;
    zd_ii = this->Allocate(zd_ii, _x, _y, _z, _t);
    double ****xd_ij = NULL;
    xd_ij = this->Allocate(xd_ij, _x, _y, _z, _t);
    double ****yd_ij = NULL;
    yd_ij = this->Allocate(yd_ij, _x, _y, _z, _t);
    double ****zd_ij = NULL;
    zd_ij = this->Allocate(zd_ij, _x, _y, _z, _t);

    // Values of the B-spline basis functions and its derivative (assuming that we compute the bending energy only at the control point location i, j, k)
    double b[3]    = {1.0/6.0, 2.0/3.0, 1.0/6.0};
    double b_i[3]  = {-0.5, 0, 0.5};
    double b_ii[3] = {1.0, -2.0, 1.0};

    for (t = 0; t < _t; t++){
        for (k = 0; k < _z; k++) {
            for (j = 0; j < _y; j++) {
                for (i = 0; i < _x; i++) {
                    for (J = j-1; J < j+2; J++) {
                        // Get B-spline basis
                        B_J    = b   [J-(j-1)];
                        B_J_I  = b_i [J-(j-1)];
                        B_J_II = b_ii[J-(j-1)];
                        for (I = i-1; I < i+2; I++) {
                            // Get B-spline basis
                            B_I    = b   [I-(i-1)];
                            B_I_I  = b_i [I-(i-1)];
                            B_I_II = b_ii[I-(i-1)];

                            v = B_I * B_J_II;
                            yd_jj[t][k][j][i] += 2 * _ydata[t][k][J][I] * v;
                            xd_jj[t][k][j][i] += 2 * _xdata[t][k][J][I] * v;

                            v = B_I_II * B_J;
                            yd_ii[t][k][j][i] += 2 * _ydata[t][k][J][I] * v;
                            xd_ii[t][k][j][i] += 2 * _xdata[t][k][J][I] * v;

                            v = B_I_I * B_J_I;
                            yd_ij[t][k][j][i] += 4 * _ydata[t][k][J][I] * v;
                            xd_ij[t][k][j][i] += 4 * _xdata[t][k][J][I] * v;

                        }
                    }
                }
            }
        }
    }


    for (t = 0; t < _t; t++){
        for (k = 0; k < _z; k++) {
            for (j = 0; j < _y; j++) {
                for (i = 0; i < _x; i++) {
                    // Initialize tmp variables
                    tmp[0] = 0;
                    tmp[1] = 0;
                    tmp[2] = 0;
                    for (J = j-1; J < j+2; J++) {
                        // Get B-spline basis
                        B_J    = b   [J-(j-1)];
                        B_J_I  = b_i [J-(j-1)];
                        B_J_II = b_ii[J-(j-1)];
                        for (I = i-1; I < i+2; I++) {
                            // Get B-spline basis
                            B_I    = b   [I-(i-1)];
                            B_I_I  = b_i [I-(i-1)];
                            B_I_II = b_ii[I-(i-1)];

                            v = B_I * B_J_II;
                            tmp[1] += yd_jj[t][k][J][I] * v;
                            tmp[0] += xd_jj[t][k][J][I] * v;

                            v = B_I_II * B_J;
                            tmp[1] += yd_ii[t][k][J][I] * v;
                            tmp[0] += xd_ii[t][k][J][I] * v;

                            v = B_I_I * B_J_I;
                            tmp[1] += yd_ij[t][k][J][I] * v;
                            tmp[0] += xd_ij[t][k][J][I] * v;

                        }
                    }
                    index = this->LatticeToIndex(i, j, k, t);
                    gradient[index]     += -tmp[0];
                    gradient[index+n]   += -tmp[1];
                    gradient[index+2*n] += 0;
                }
            }
        }
    }

    this->Deallocate(xd_jj, _x, _y, _z,_t);
    this->Deallocate(yd_jj, _x, _y, _z,_t);
    this->Deallocate(zd_jj, _x, _y, _z,_t);
    this->Deallocate(xd_ii, _x, _y, _z,_t);
    this->Deallocate(yd_ii, _x, _y, _z,_t);
    this->Deallocate(zd_ii, _x, _y, _z,_t);
    this->Deallocate(xd_ij, _x, _y, _z,_t);
    this->Deallocate(yd_ij, _x, _y, _z,_t);
    this->Deallocate(zd_ij, _x, _y, _z,_t);
}

void irtkBSplineFreeFormTransformationPeriodic::BendingGradient3D(double *gradient)
{
    int I, J, K, index, i, j, k, t, n = _x*_y*_z*_t;
    double B_K, B_J, B_I, B_K_I, B_J_I, B_I_I, B_K_II, B_J_II, B_I_II, v, tmp[3];

    // Derivatives
    double ****xd_kk = NULL;
    xd_kk = this->Allocate(xd_kk, _x, _y, _z, _t);
    double ****xd_jj = NULL;
    xd_jj = this->Allocate(xd_jj, _x, _y, _z, _t);
    double ****xd_ii = NULL;
    xd_ii = this->Allocate(xd_ii, _x, _y, _z, _t);
    double ****xd_ij = NULL;
    xd_ij = this->Allocate(xd_ij, _x, _y, _z, _t);
    double ****xd_ik = NULL;
    xd_ik = this->Allocate(xd_ik, _x, _y, _z, _t);
    double ****xd_jk = NULL;
    xd_jk = this->Allocate(xd_jk, _x, _y, _z, _t);

    double ****yd_kk = NULL;
    yd_kk = this->Allocate(yd_kk, _x, _y, _z, _t);
    double ****yd_jj = NULL;
    yd_jj = this->Allocate(yd_jj, _x, _y, _z, _t);
    double ****yd_ii = NULL;
    yd_ii = this->Allocate(yd_ii, _x, _y, _z, _t);
    double ****yd_ij = NULL;
    yd_ij = this->Allocate(yd_ij, _x, _y, _z, _t);
    double ****yd_ik = NULL;
    yd_ik = this->Allocate(yd_ik, _x, _y, _z, _t);
    double ****yd_jk = NULL;
    yd_jk = this->Allocate(yd_jk, _x, _y, _z, _t);

    double ****zd_kk = NULL;
    zd_kk = this->Allocate(zd_kk, _x, _y, _z, _t);
    double ****zd_jj = NULL;
    zd_jj = this->Allocate(zd_jj, _x, _y, _z, _t);
    double ****zd_ii = NULL;
    zd_ii = this->Allocate(zd_ii, _x, _y, _z, _t);
    double ****zd_ij = NULL;
    zd_ij = this->Allocate(zd_ij, _x, _y, _z, _t);
    double ****zd_ik = NULL;
    zd_ik = this->Allocate(zd_ik, _x, _y, _z, _t);
    double ****zd_jk = NULL;
    zd_jk = this->Allocate(zd_jk, _x, _y, _z, _t);

    // Values of the B-spline basis functions and its derivative (assuming that we compute the bending energy only at the control point location i, j, k)
    double b[3]    = {1.0/6.0, 2.0/3.0, 1.0/6.0};
    double b_i[3]  = {-0.5, 0, 0.5};
    double b_ii[3] = {1.0, -2.0, 1.0};

    for (t = 0; t < _t; t++){
        for (k = 0; k < _z; k++) {
            for (j = 0; j < _y; j++) {
                for (i = 0; i < _x; i++) {
                    for (K = k-1; K < k+2; K++) {
                        // Get B-spline basis
                        B_K    = b   [K-(k-1)];
                        B_K_I  = b_i [K-(k-1)];
                        B_K_II = b_ii[K-(k-1)];
                        for (J = j-1; J < j+2; J++) {
                            // Get B-spline basis
                            B_J    = b   [J-(j-1)];
                            B_J_I  = b_i [J-(j-1)];
                            B_J_II = b_ii[J-(j-1)];
                            for (I = i-1; I < i+2; I++) {
                                // Get B-spline basis
                                B_I    = b   [I-(i-1)];
                                B_I_I  = b_i [I-(i-1)];
                                B_I_II = b_ii[I-(i-1)];

                                v = B_I * B_J * B_K_II;
                                zd_kk[t][k][j][i] += 2 * _zdata[t][K][J][I] * v;
                                yd_kk[t][k][j][i] += 2 * _ydata[t][K][J][I] * v;
                                xd_kk[t][k][j][i] += 2 * _xdata[t][K][J][I] * v;

                                v = B_I * B_J_II * B_K;
                                zd_jj[t][k][j][i] += 2 * _zdata[t][K][J][I] * v;
                                yd_jj[t][k][j][i] += 2 * _ydata[t][K][J][I] * v;
                                xd_jj[t][k][j][i] += 2 * _xdata[t][K][J][I] * v;

                                v = B_I_II * B_J * B_K;
                                zd_ii[t][k][j][i] += 2 * _zdata[t][K][J][I] * v;
                                yd_ii[t][k][j][i] += 2 * _ydata[t][K][J][I] * v;
                                xd_ii[t][k][j][i] += 2 * _xdata[t][K][J][I] * v;

                                v = B_I_I * B_J_I * B_K;
                                zd_ij[t][k][j][i] += 4 * _zdata[t][K][J][I] * v;
                                yd_ij[t][k][j][i] += 4 * _ydata[t][K][J][I] * v;
                                xd_ij[t][k][j][i] += 4 * _xdata[t][K][J][I] * v;

                                v = B_I_I * B_J * B_K_I;
                                zd_ik[t][k][j][i] += 4 * _zdata[t][K][J][I] * v;
                                yd_ik[t][k][j][i] += 4 * _ydata[t][K][J][I] * v;
                                xd_ik[t][k][j][i] += 4 * _xdata[t][K][J][I] * v;

                                v = B_I * B_J_I * B_K_I;
                                zd_jk[t][k][j][i] += 4 * _zdata[t][K][J][I] * v;
                                yd_jk[t][k][j][i] += 4 * _ydata[t][K][J][I] * v;
                                xd_jk[t][k][j][i] += 4 * _xdata[t][K][J][I] * v;
                            }
                        }
                    }
                }
            }
        }
    }

    for(t = 0; t < _t; t++){
        for (k = 0; k < _z; k++) {
            for (j = 0; j < _y; j++) {
                for (i = 0; i < _x; i++) {
                    // Initialize tmp variables
                    tmp[0] = 0;
                    tmp[1] = 0;
                    tmp[2] = 0;
                    for (K = k-1; K < k+2; K++) {
                        // Get B-spline basis
                        B_K    = b   [K-(k-1)];
                        B_K_I  = b_i [K-(k-1)];
                        B_K_II = b_ii[K-(k-1)];
                        for (J = j-1; J < j+2; J++) {
                            // Get B-spline basis
                            B_J    = b   [J-(j-1)];
                            B_J_I  = b_i [J-(j-1)];
                            B_J_II = b_ii[J-(j-1)];
                            for (I = i-1; I < i+2; I++) {
                                // Get B-spline basis
                                B_I    = b   [I-(i-1)];
                                B_I_I  = b_i [I-(i-1)];
                                B_I_II = b_ii[I-(i-1)];

                                v = B_I * B_J * B_K_II;
                                tmp[2] += zd_kk[t][K][J][I] * v;
                                tmp[1] += yd_kk[t][K][J][I] * v;
                                tmp[0] += xd_kk[t][K][J][I] * v;

                                v = B_I * B_J_II * B_K;
                                tmp[2] += zd_jj[t][K][J][I] * v;
                                tmp[1] += yd_jj[t][K][J][I] * v;
                                tmp[0] += xd_jj[t][K][J][I] * v;

                                v = B_I_II * B_J * B_K;
                                tmp[2] += zd_ii[t][K][J][I] * v;
                                tmp[1] += yd_ii[t][K][J][I] * v;
                                tmp[0] += xd_ii[t][K][J][I] * v;

                                v = B_I_I * B_J_I * B_K;
                                tmp[2] += zd_ij[t][K][J][I] * v;
                                tmp[1] += yd_ij[t][K][J][I] * v;
                                tmp[0] += xd_ij[t][K][J][I] * v;

                                v = B_I_I * B_J * B_K_I;
                                tmp[2] += zd_ik[t][K][J][I] * v;
                                tmp[1] += yd_ik[t][K][J][I] * v;
                                tmp[0] += xd_ik[t][K][J][I] * v;

                                v = B_I * B_J_I * B_K_I;
                                tmp[2] += zd_jk[t][K][J][I] * v;
                                tmp[1] += yd_jk[t][K][J][I] * v;
                                tmp[0] += xd_jk[t][K][J][I] * v;
                            }
                        }
                    }
                    index = this->LatticeToIndex(i, j, k, t);
                    gradient[index]     += -tmp[0];
                    gradient[index+n]   += -tmp[1];
                    gradient[index+2*n] += -tmp[2];
                }
            }
        }
    }

    this->Deallocate(xd_kk, _x, _y, _z, _t);
    this->Deallocate(xd_jj, _x, _y, _z, _t);
    this->Deallocate(xd_ii, _x, _y, _z, _t);
    this->Deallocate(xd_ij, _x, _y, _z, _t);
    this->Deallocate(xd_ik, _x, _y, _z, _t);
    this->Deallocate(xd_jk, _x, _y, _z, _t);

    this->Deallocate(yd_kk, _x, _y, _z, _t);
    this->Deallocate(yd_jj, _x, _y, _z, _t);
    this->Deallocate(yd_ii, _x, _y, _z, _t);
    this->Deallocate(yd_ij, _x, _y, _z, _t);
    this->Deallocate(yd_ik, _x, _y, _z, _t);
    this->Deallocate(yd_jk, _x, _y, _z, _t);

    this->Deallocate(zd_kk, _x, _y, _z, _t);
    this->Deallocate(zd_jj, _x, _y, _z, _t);
    this->Deallocate(zd_ii, _x, _y, _z, _t);
    this->Deallocate(zd_ij, _x, _y, _z, _t);
    this->Deallocate(zd_ik, _x, _y, _z, _t);
    this->Deallocate(zd_jk, _x, _y, _z, _t);
}

void irtkBSplineFreeFormTransformationPeriodic::BendingGradient(double *gradient)
{
    if (_z == 1) {
        this->BendingGradient2D(gradient);
    } else {
        this->BendingGradient3D(gradient);
    }
}

int irtkBSplineFreeFormTransformationPeriodic::CheckHeader(char *)
{
    return false;
}

void irtkBSplineFreeFormTransformationPeriodic::Interpolate(double *, double *, double *)
{
    cerr << "irtkBSplineFreeFormTransformationPeriodic::Interpolate: Not implemented yet" << endl;
    exit(1);
}

void irtkBSplineFreeFormTransformationPeriodic::Subdivide()
{
    int noTime = 1;
    if (noTime){
        this->Subdivide3D();
    }
}

void irtkBSplineFreeFormTransformationPeriodic::Subdivide2D()
{
    int i, j, k, l, i1, j1, k1, i2, j2, k2;

    // Weights for subdivision
    double w[2][3];
    w[1][0] = 0;
    w[1][1] = 1.0/2.0;
    w[1][2] = 1.0/2.0;
    w[0][0] = 1.0/8.0;
    w[0][1] = 6.0/8.0;
    w[0][2] = 1.0/8.0;

    // Allocate memory for new control points
    double ****x = NULL;
    double ****y = NULL;
    double ****z = NULL;
    x = this->Allocate(x, 2*_x-1, 2*_y-1, _z, _t);
    y = this->Allocate(y, 2*_x-1, 2*_y-1, _z, _t);
    z = this->Allocate(z, 2*_x-1, 2*_y-1, _z, _t);

    for (i = 0; i < _x; i++) {
        for (j = 0; j < _y; j++) {
            for (k = 0; k < _z; k++) {
                for (l = 0; l < _t; l++) {
                    for (i1 = 0; i1 < 2; i1++) {
                        for (j1 = 0; j1 < 2; j1++) {
                            x[l][k][2*j+j1][2*i+i1] = 0;
                            y[l][k][2*j+j1][2*i+i1] = 0;
                            z[l][k][2*j+j1][2*i+i1] = 0;
                            for (i2 = 0; i2 < 3; i2++) {
                                for (j2 = 0; j2 < 3; j2++) {
                                    x[l][k][2*j+j1][2*i+i1] += w[i1][i2] * w[j1][j2] * _xdata[l][k][j+j2-1][i+i2-1];
                                    y[l][k][2*j+j1][2*i+i1] += w[i1][i2] * w[j1][j2] * _ydata[l][k][j+j2-1][i+i2-1];
                                    z[l][k][2*j+j1][2*i+i1] += w[i1][i2] * w[j1][j2] * _zdata[l][k][j+j2-1][i+i2-1];
                                }
                            }
                        }
                    }
                }
            }
        }
    }

    // Deallocate points
    this->Deallocate(_xdata, _x, _y, _z, _t);
    this->Deallocate(_ydata, _x, _y, _z, _t);
    this->Deallocate(_zdata, _x, _y, _z, _t);
    delete []_status;

    // Update pointers to control points
    _xdata = x;
    _ydata = y;
    _zdata = z;

    // Increase number of control points
    _x = 2*_x - 1;
    _y = 2*_y - 1;
    _z = _z;
    _t = _t;

    // Recalculate control point spacing
    _dx /= 2.0;
    _dy /= 2.0;
    _dz  = _dz;
    _dt  = _dt;

    // Update transformation matrix
    this->UpdateMatrix();

    // Initialize memory for control point status
    _status = new _Status[3*_x*_y*_z*_t];
    for (i = 0; i < 3*_x*_y*_z*_t; i++) {
        _status[i] = _Active;
    }
}

void irtkBSplineFreeFormTransformationPeriodic::Subdivide3D()
{
    int i, j, k, l, i1, j1, k1, i2, j2, k2;

    // Weights for subdivision
    double w[2][3];
    w[1][0] = 0;
    w[1][1] = 1.0/2.0;
    w[1][2] = 1.0/2.0;
    w[0][0] = 1.0/8.0;
    w[0][1] = 6.0/8.0;
    w[0][2] = 1.0/8.0;

    // Allocate memory for new control points
    double ****x = NULL;
    double ****y = NULL;
    double ****z = NULL;
    x = this->Allocate(x, 2*_x-1, 2*_y-1, 2*_z-1, _t);
    y = this->Allocate(y, 2*_x-1, 2*_y-1, 2*_z-1, _t);
    z = this->Allocate(z, 2*_x-1, 2*_y-1, 2*_z-1, _t);

    for (i = 0; i < _x; i++) {
        for (j = 0; j < _y; j++) {
            for (k = 0; k < _z; k++) {
                for (l = 0; l < _t; l++) {
                    for (i1 = 0; i1 < 2; i1++) {
                        for (j1 = 0; j1 < 2; j1++) {
                            for (k1 = 0; k1 < 2; k1++) {
                                x[l][2*k+k1][2*j+j1][2*i+i1] = 0;
                                y[l][2*k+k1][2*j+j1][2*i+i1] = 0;
                                z[l][2*k+k1][2*j+j1][2*i+i1] = 0;
                                for (i2 = 0; i2 < 3; i2++) {
                                    for (j2 = 0; j2 < 3; j2++) {
                                        for (k2 = 0; k2 < 3; k2++) {
                                            x[l][2*k+k1][2*j+j1][2*i+i1] += w[i1][i2] * w[j1][j2] * w[k1][k2] * _xdata[l][k+k2-1][j+j2-1][i+i2-1];
                                            y[l][2*k+k1][2*j+j1][2*i+i1] += w[i1][i2] * w[j1][j2] * w[k1][k2] * _ydata[l][k+k2-1][j+j2-1][i+i2-1];
                                            z[l][2*k+k1][2*j+j1][2*i+i1] += w[i1][i2] * w[j1][j2] * w[k1][k2] * _zdata[l][k+k2-1][j+j2-1][i+i2-1];
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
    }

    // Deallocate points
    this->Deallocate(_xdata, _x, _y, _z, _t);
    this->Deallocate(_ydata, _x, _y, _z, _t);
    this->Deallocate(_zdata, _x, _y, _z, _t);
    delete []_status;

    // Update pointers to control points
    _xdata = x;
    _ydata = y;
    _zdata = z;

    // Increase number of control points
    _x = 2*_x - 1;
    _y = 2*_y - 1;
    _z = 2*_z - 1;
    _t = _t;

    // Recalculate control point spacing
    _dx /= 2.0;
    _dy /= 2.0;
    _dz /= 2.0;
    _dt  = _dt;

    // Update transformation matrix
    this->UpdateMatrix();

    // Initialize memory for control point status
    _status = new _Status[3*_x*_y*_z*_t];
    for (i = 0; i < 3*_x*_y*_z*_t; i++) {
        _status[i] = _Active;
    }
}

void irtkBSplineFreeFormTransformationPeriodic::Subdivide4D()
{
    int i, j, k, l, i1, j1, k1, l1, i2, j2, k2, l2;
    int ip, jp, kp, lp;

    // Weights for subdivision
    double w[2][3];
    w[1][0] = 0;
    w[1][1] = 1.0/2.0;
    w[1][2] = 1.0/2.0;
    w[0][0] = 1.0/8.0;
    w[0][1] = 6.0/8.0;
    w[0][2] = 1.0/8.0;

    // Allocate memory for new control points
    double ****x = NULL;
    double ****y = NULL;
    double ****z = NULL;
    x = this->Allocate(x, 2*_x-1, 2*_y-1, 2*_z-1, 2*_t-1);
    y = this->Allocate(y, 2*_x-1, 2*_y-1, 2*_z-1, 2*_t-1);
    z = this->Allocate(z, 2*_x-1, 2*_y-1, 2*_z-1, 2*_t-1);

    for (i = 0; i < _x; i++) {
        for (j = 0; j < _y; j++) {
            for (k = 0; k < _z; k++) {
                for (l = 0; l < _t; l++) {
                    for (i1 = 0; i1 < 2; i1++) {
                        for (j1 = 0; j1 < 2; j1++) {
                            for (k1 = 0; k1 < 2; k1++) {
                                for (l1 = 0; l1 < 2; l1++) {
                                    x[2*l+l1][2*k+k1][2*j+j1][2*i+i1] = 0;
                                    y[2*l+l1][2*k+k1][2*j+j1][2*i+i1] = 0;
                                    z[2*l+l1][2*k+k1][2*j+j1][2*i+i1] = 0;
                                    for (i2 = 0; i2 < 3; i2++) {
                                        for (j2 = 0; j2 < 3; j2++) {
                                            for (k2 = 0; k2 < 3; k2++) {
                                                for (l2 = 0; l2 < 3; l2++) {
                                                    ip = i+i2-1;
                                                    jp = j+j2-1;
                                                    kp = k+k2-1;
                                                    // periodic model (time)
                                                    lp = l+l2-1;
                                                    if (_periodic) {
                                                        if (lp < 0)
                                                            lp += _t-1;
                                                        if (lp >= _t-1)
                                                            lp -= _t-1;
                                                    }

                                                    x[2*l+l1][2*k+k1][2*j+j1][2*i+i1] += w[i1][i2] * w[j1][j2] * w[k1][k2] * w[l1][l2] * _xdata[lp][kp][jp][ip];
                                                    y[2*l+l1][2*k+k1][2*j+j1][2*i+i1] += w[i1][i2] * w[j1][j2] * w[k1][k2] * w[l1][l2] * _ydata[lp][kp][jp][ip];
                                                    z[2*l+l1][2*k+k1][2*j+j1][2*i+i1] += w[i1][i2] * w[j1][j2] * w[k1][k2] * w[l1][l2] * _zdata[lp][kp][jp][ip];
                                                }
                                            }
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
    }

    // Deallocate points
    this->Deallocate(_xdata, _x, _y, _z, _t);
    this->Deallocate(_ydata, _x, _y, _z, _t);
    this->Deallocate(_zdata, _x, _y, _z, _t);
    delete []_status;

    // Update pointers to control points
    _xdata = x;
    _ydata = y;
    _zdata = z;

    // Increase number of control points
    _x = 2*_x - 1;
    _y = 2*_y - 1;
    _z = 2*_z - 1;
    _t = 2*_t - 1;

    // Recalculate control point spacing
    _dx /= 2.0;
    _dy /= 2.0;
    _dz /= 2.0;
    _dt /= 2.0;

    // Update transformation matrix
    this->UpdateMatrix();

    // Initialize memory for control point status
    _status = new _Status[3*_x*_y*_z*_t];
    for (i = 0; i < 3*_x*_y*_z*_t; i++) {
        _status[i] = _Active;
    }
}

void irtkBSplineFreeFormTransformationPeriodic::BoundingBoxCP(int index, irtkPoint &p1, irtkPoint &p2, double fraction) const
{
    int i, j, k, l;

    if (index >= _x*_y*_z*_t) {
        index -= _x*_y*_z*_t;
        if (index >= _x*_y*_z*_t) {
            index -= _x*_y*_z*_t;
        }
    }
    l = index/(_z*_y*_x);
    k = index%(_z*_y*_x)/(_y*_x);
    j = index%(_z*_y*_x)%(_y*_x)/_x;
    i = index%(_z*_y*_x)%(_y*_x)%_x;
    p1 = irtkPoint(i-2*fraction, j-2*fraction, k-2*fraction);
    p2 = irtkPoint(i+2*fraction, j+2*fraction, k+2*fraction);
    this->LatticeToWorld(p1);
    this->LatticeToWorld(p2);
}

void irtkBSplineFreeFormTransformationPeriodic::BoundingBoxCP(int index, double &x1, double &y1, double &z1, double &t1, double &x2, double &y2, double &z2, double &t2, double fraction) const
{
    if (index >= _x*_y*_z*_t) {
        index -= _x*_y*_z*_t;
        if (index >= _x*_y*_z*_t) {
            index -= _x*_y*_z*_t;
        }
    }
    x1 = index%(_z*_y*_x)%(_y*_x)%_x-2*fraction;
    y1 = index%(_z*_y*_x)%(_y*_x)/_x-2*fraction;
    z1 = index%(_z*_y*_x)/(_y*_x)-2*fraction;
    t1 = index/(_z*_y*_x)-2*fraction;					///////////////////////////////////////////////////
    x2 = index%(_z*_y*_x)%(_y*_x)%_x+2*fraction;
    y2 = index%(_z*_y*_x)%(_y*_x)/_x+2*fraction;
    z2 = index%(_z*_y*_x)/(_y*_x)+2*fraction;
    t2 = index/(_z*_y*_x)+2*fraction;					///////////////////////////////////////////////////
    this->LatticeToWorld(x1, y1, z1);
    this->LatticeToWorld(x2, y2, z2);
    //  t1 = this->LatticeToTime(t1);
    //  t2 = this->LatticeToTime(t2);
}

void irtkBSplineFreeFormTransformationPeriodic::BoundingBoxImage(irtkGreyImage *image, int index, int &i1, int &j1, int &k1, int &i2, int &j2, int &k2, int &t1, int &t2, double fraction) const
{
    cerr << "irtkBSplineFreeFormTransformationPeriodic::BoundingBoxImage: Not applicable for this Transformation" << endl;
    exit(1);
}

void irtkBSplineFreeFormTransformationPeriodic::BoundingBoxImage(irtkGreyImage *image, int index, int &i1, int &j1, int &k1, int &i2, int &j2, int &k2, double fraction) const
{
    cerr << "irtkBSplineFreeFormTransformationPeriodic::BoundingBoxImage: Not applicable for this Transformation" << endl;
    exit(1);
}

void irtkBSplineFreeFormTransformationPeriodic::BoundingBoxImage(irtkGreyImage *image, int index, int &i1, int &j1, int &k1, int &i2, int &j2, int &k2, double &t1, double &t2, double fraction) const
{
    double x1, y1, z1, x2, y2, z2;

    // Calculate bounding box in world coordinates
    this->BoundingBoxCP(index, x1, y1, z1, t1, x2, y2, z2, t2, fraction);

    // Transform world coordinates to image coordinates
    image->WorldToImage(x1, y1, z1);
    image->WorldToImage(x2, y2, z2);
    //  t1 = image->TimeToImage(t1);
    //  t2 = image->TimeToImage(t1);
    //  cout<<"[BondingBox] x: "<<x1<<" - "<<x2<<" ; y: "<<y1<<" - "<<y2<<" ; z: "<<z1<<" - "<<z2<<endl;
    // Calculate bounding box in image coordinates
    i1 = (x1 < 0) ? 0 : int(x1)+1;
    j1 = (y1 < 0) ? 0 : int(y1)+1;
    k1 = (z1 < 0) ? 0 : int(z1)+1;
    //  l1 = (t1 < 0) ? 0 : int(t1)+1;
    i2 = (x2 < 0) ? -1 : int(x2);
    j2 = (y2 < 0) ? -1 : int(y2);
    k2 = (z2 < 0) ? -1 : int(z2);
    i2 = (i2 >= image->GetX()) ? image->GetX()-1 : i2;
    j2 = (j2 >= image->GetY()) ? image->GetY()-1 : j2;
    k2 = (k2 >= image->GetZ()) ? image->GetZ()-1 : k2;
    //  i2 = (int(x2) >= image->GetX()) ? image->GetX()-1 : int(x2);
    //  j2 = (int(y2) >= image->GetY()) ? image->GetY()-1 : int(y2);
    //  k2 = (int(z2) >= image->GetZ()) ? image->GetZ()-1 : int(z2);
    //  l2 = (int(t2) >= image->GetT()) ? image->GetT()-1 : int(t2);
    //  cout<<"[BondingBox] i: "<<i1<<" - "<<i2<<" ; j: "<<j1<<" - "<<j2<<" ; k: "<<k1<<" - "<<k2<<endl;
}


void irtkBSplineFreeFormTransformationPeriodic::Print()
{
    // Print keyword and no. of DOFs
    cout << "BSplineFFD4D: " << this->NumberOfDOFs() << endl;

    // Print no. of control points
    cout << "Control points: " << _x << " x " << _y << " x " << _z << " x " << _t << endl;
    cout << "Spacing: " << _dx << " x " << _dy << " x " << _dz << " x " << _dt << endl;
    cout << "Time: " << _tMin << " to " << _tMax << endl;
    cout << "Origin: " << _origin._x << " " << _origin._y << " " << _origin._z << endl;
    // Print x-axis
    cout << "X-axis is " << _xaxis[0] << " " << _xaxis[1] << " " << _xaxis[2] << endl;
    // Print x-axis
    cout << "Y-axis is " << _yaxis[0] << " " << _yaxis[1] << " " << _yaxis[2] << endl;
    // Print x-axis
    cout << "Z-axis is " << _zaxis[0] << " " << _zaxis[1] << " " << _zaxis[2] << endl;
}

irtkCifstream& irtkBSplineFreeFormTransformationPeriodic::Read(irtkCifstream& from)
{
    double *data;
    int i, j, k, l, index;
    unsigned int magic_no, trans_type;

    // Read magic no. for transformations
    from.ReadAsUInt(&magic_no, 1);
    if (magic_no != IRTKTRANSFORMATION_MAGIC) {
        cerr << "irtkBSplineFreeFormTransformationPeriodic::Read: Not a vaild transformation file" << endl;
        exit(1);
    }

    // Read transformation type
    from.ReadAsUInt(&trans_type, 1);
    if (trans_type != IRTKTRANSFORMATION_PERIODIC) {
        cerr << "irtkBSplineFreeFormTransformationPeriodic::Read: Not a vaild B-Spline FFD transformation" << endl;
        exit(1);
    }

    // Free memory if necessary
    _xdata = Deallocate(_xdata, _x, _y, _z, _t);
    _ydata = Deallocate(_ydata, _x, _y, _z, _t);
    _zdata = Deallocate(_zdata, _x, _y, _z, _t);
    delete []_status;

    // Read no of control points
    from.ReadAsInt(&_x, 1);
    from.ReadAsInt(&_y, 1);
    from.ReadAsInt(&_z, 1);
    from.ReadAsInt(&_t, 1);

    // Read orientation of bounding box
    from.ReadAsDouble(_xaxis, 3);
    from.ReadAsDouble(_yaxis, 3);
    from.ReadAsDouble(_zaxis, 3);

    // Read spacing of bounding box
    from.ReadAsDouble(&_dx, 1);
    from.ReadAsDouble(&_dy, 1);
    from.ReadAsDouble(&_dz, 1);
    from.ReadAsDouble(&_dt, 1);

    // Read spacing of bounding box
    from.ReadAsDouble(&_origin._x, 1);
    from.ReadAsDouble(&_origin._y, 1);
    from.ReadAsDouble(&_origin._z, 1);

    // Read time domain
    from.ReadAsDouble(&_tMin, 1);
    from.ReadAsDouble(&_tMax, 1);

    // Initialize control points
    _xdata = Allocate(_xdata, _x, _y, _z, _t);
    _ydata = Allocate(_ydata, _x, _y, _z, _t);
    _zdata = Allocate(_zdata, _x, _y, _z, _t);

    // Initialize control points
    for (i = -2; i < _x+2; i++) {
        for (j = -2; j < _y+2; j++) {
            for (k = -2; k < _z+2; k++) {
                for (l = -2; l < _t+2; l++) {
                    _xdata[l][k][j][i] = 0;
                    _ydata[l][k][j][i] = 0;
                    _zdata[l][k][j][i] = 0;
                }
            }
        }
    }

    // Allocate temporary memory
    data = new double[3*_x*_y*_z*_t];

    // Read control point data
    from.ReadAsDouble(data, 3*_x*_y*_z*_t);

    // Convert data
    index = 0;
    for (i = 0; i < _x; i++) {
        for (j = 0; j < _y; j++) {
            for (k = 0; k < _z; k++) {
                for (l = 0; l < _t; l++) {
                    _xdata[l][k][j][i] = data[index];
                    _ydata[l][k][j][i] = data[index+1];
                    _zdata[l][k][j][i] = data[index+2];
                    index += 3;
                }
            }
        }
    }

    // Free temporary memory
    delete []data;

    // Initialize memory for control point status
    _status = new _Status[3*_x*_y*_z*_t];

    // Read control point status
    from.ReadAsInt((int *)_status, 3*_x*_y*_z*_t);

    // Update transformation matrix
    this->UpdateMatrix();

    return from;
}

irtkCofstream& irtkBSplineFreeFormTransformationPeriodic::Write(irtkCofstream& to)
{
    double *data;
    int i, j, k, l, index;
    unsigned int magic_no, trans_type;

    // Write magic no. for transformations
    magic_no = IRTKTRANSFORMATION_MAGIC;
    to.WriteAsUInt(&magic_no, 1);

    // Write transformation type
    trans_type = IRTKTRANSFORMATION_PERIODIC;
    to.WriteAsUInt(&trans_type, 1);

    // Write no of control points
    to.WriteAsInt(&_x, 1);
    to.WriteAsInt(&_y, 1);
    to.WriteAsInt(&_z, 1);
    to.WriteAsInt(&_t, 1);

    // Write orientation of bounding box
    to.WriteAsDouble(_xaxis, 3);
    to.WriteAsDouble(_yaxis, 3);
    to.WriteAsDouble(_zaxis, 3);

    // Write spacing of bounding box
    to.WriteAsDouble(&_dx, 1);
    to.WriteAsDouble(&_dy, 1);
    to.WriteAsDouble(&_dz, 1);
    to.WriteAsDouble(&_dt, 1);

    // Write spacing of bounding box
    to.WriteAsDouble(&_origin._x, 1);
    to.WriteAsDouble(&_origin._y, 1);
    to.WriteAsDouble(&_origin._z, 1);

    // Write time domain
    to.WriteAsDouble(&_tMin, 1);
    to.WriteAsDouble(&_tMax, 1);

    // Allocate temporary memory
    data = new double[3*_x*_y*_z*_t];

    // Convert data
    index = 0;
    for (i = 0; i < _x; i++) {
        for (j = 0; j < _y; j++) {
            for (k = 0; k < _z; k++) {
                for (l = 0; l < _t; l++) {
                    data[index]   = _xdata[l][k][j][i];
                    data[index+1] = _ydata[l][k][j][i];
                    data[index+2] = _zdata[l][k][j][i];
                    index += 3;
                }
            }
        }
    }

    // Write control point data
    to.WriteAsDouble(data, 3*_x*_y*_z*_t);

    // Free temporary memory
    delete []data;

    // Write control point status
    to.WriteAsInt((int *)_status, 3*_x*_y*_z*_t);

    return to;
}



irtkBSplineFreeFormTransformation3D *irtkBSplineFreeFormTransformationPeriodic::Compute3DFFD1(double time)
{
    int i ,j ,k, t, d, V, l;
    double v, B_L, CPx, CPy, CPz;

    // initialise 3D FFD
    double x1, y1, z1, x2, y2, z2;
    x1 = y1 = z1 = 0;
    this->LatticeToWorld(x1, y1, z1);
    x2 = _x-1; y2 = _y-1; z2 = _z-1;
    this->LatticeToWorld(x2, y2, z2);
    irtkBSplineFreeFormTransformation3D * transform = new irtkBSplineFreeFormTransformation3D(x1, y1, z1, x2, y2, z2, _dx, _dy, _dz, _xaxis, _yaxis, _zaxis);

    // periodic model:
    if (_periodic) {
        while (time < 0)
            time += double(_t-1);
        while (time >= _t-1)
            time -= double(_t-1);
    } else {
        if ((time < -2) || (time > _t+1)) {
            cerr << "time out of bounds" << endl;
            exit(1);
        }
    }

    d = (int)floor(time)-1;
    v = time-(d+1);
    V = round(LUTSIZE*v);

    // calculate 3D FFD
    //	irtkMatrix BS = new irtkMatrix(_x*_y*_z, _x*_y*_z);
    for (i = 0; i < _x; i++) {
        for (j = 0; j < _y; j++) {
            for (k = 0; k < _z; k++) {
                CPx = 0;
                CPy = 0;
                CPz = 0;

                for (l = 0; l < 4; l++) {
                    B_L = this->LookupTable[V][l];
                    CPx += B_L * _xdata[d+l][k][j][i];
                    CPy += B_L * _ydata[d+l][k][j][i];
                    CPz += B_L * _zdata[d+l][k][j][i];
                }

                transform->Put(i, j, k, CPx, CPy, CPz);
            }
        }
    }

    return transform;
}

irtkBSplineFreeFormTransformation3D *irtkBSplineFreeFormTransformationPeriodic::Compute3DFFD2(double time)
{
    cerr << "doesn't work (use Approximate function)" << endl;
    exit(1);

    int i, j, k, t, iC, jC, kC, m, n, o, a;
    double B_I, B_J, B_K;
    double **CP = new double *[3];
    double **CP_old = new double *[3];
    for (i = 0; i < 3; i++) {
        CP[i] = new double [_x*_y*_z];
        CP_old[i] = new double [_x*_y*_z];
    }

    // initialise 3D FFD
    double x1, y1, z1, x2, y2, z2;
    x1 = y1 = z1 = 0;
    this->LatticeToWorld(x1, y1, z1);
    x2 = _x-1; y2 = _y-1; z2 = _z-1;
    this->LatticeToWorld(x2, y2, z2);
    irtkBSplineFreeFormTransformation3D * transform = new irtkBSplineFreeFormTransformation3D(x1, y1, z1, x2, y2, z2, _dx, _dy, _dz, _xaxis, _yaxis, _zaxis);

    // periodic model:
    if (_periodic) {
        while (time < 0)
            time += double(_t-1);
        while (time >= _t-1)
            time -= double(_t-1);
    } else {
        if ((time < -2) || (time > _t+1)) {
            cerr << "time out of bounds" << endl;
            exit(1);
        }
    }

    // calculate 3D FFD conversion matrix (same for all 3 displacement components (x, y, z))
    cout<<"fill matrix"<<endl;
    irtkMatrix BS = irtkMatrix(_x*_y*_z, _x*_y*_z);
    // loops over all rows
    for (i = 0; i < _x; i++) {
        for (j = 0; j < _y; j++) {
            for (k = 0; k < _z; k++) {

                // loops over all columns in current row
                for (iC = 0; iC < _x; iC++) {
                    for (jC = 0; jC < _y; jC++) {
                        for (kC = 0; kC < _z; kC++) {
                            BS((k*_x*_y + j*_x + i), (kC*_x*_y + jC*_x + iC)) = 0;

                            // loops over all BSplines
                            for (m = -1; m < 3; m++) {
                                for (n = -1; n < 3; n++) {
                                    for (o = -1; o < 3; o++) {
                                        if (i+o == iC && j+n == jC && k+m == kC) {
                                            B_K = this->LookupTable[0][m+1];
                                            B_J = this->LookupTable[0][n+1] * B_K;
                                            B_I = this->LookupTable[0][o+1] * B_J;

                                            BS((k*_x*_y + j*_x + i), (kC*_x*_y + jC*_x + iC)) = B_I;
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
    }
    // invert matrix
    cout<<"invert matrix"<<endl;
    BS.Invert();

    // fill in diplacement vector
    cout<<"fill displacement vector"<<endl;
    double xDisp, yDisp, zDisp;
    for (i = 0; i < _x; i++) {
        for (j = 0; j < _y; j++) {
            for (k = 0; k < _z; k++) {
                xDisp = i;
                yDisp = j;
                zDisp = k;
                this->LocalDisplacement(xDisp, yDisp, zDisp, time);
                CP_old[0][k*_x*_y + j*_x + i] = xDisp;
                CP_old[1][k*_x*_y + j*_x + i] = yDisp;
                CP_old[2][k*_x*_y + j*_x + i] = zDisp;
            }
        }
    }

    // compute new control point values
    cout<<"compute new control point values"<<endl;
    for (i = 0; i < _x; i++) {
        for (j = 0; j < _y; j++) {
            for (k = 0; k < _z; k++) {
                CP[a][k*_x*_y + j*_x + i] = 0.;

                for (iC = 0; iC < _x; iC++) {
                    for (jC = 0; jC < _y; jC++) {
                        for (kC = 0; kC < _z; kC++) {
                            for (a = 0; a < 3; a++) {
                                CP[a][k*_x*_y + j*_x + i] += BS((k*_x*_y + j*_x + i), (kC*_x*_y + jC*_x + iC))*CP_old[a][kC*_x*_y + jC*_x + iC];
                            }
                        }
                    }
                }
            }
        }
    }

    // write to transformation
    cout<<"write new transformation"<<endl;
    for (i = 0; i < _x; i++) {
        for (j = 0; j < _y; j++) {
            for (k = 0; k < _z; k++) {
                transform->Put(i, j, k, CP[0][k*_x*_y + j*_x + i], CP[1][k*_x*_y + j*_x + i], CP[2][k*_x*_y + j*_x + i]);
            }
        }
    }

    return transform;
}


irtkBSplineFreeFormTransformation3D *irtkBSplineFreeFormTransformationPeriodic::Compute3DFFD3(double time)
{
    cerr << "doesn't work (use Approximate function)" << endl;
    exit(1);

    int i, j, k, t, iC, jC, kC, m, n, o, a;
    double B_I, B_J, B_K, delta, stepsize, stepsize_max, stepsize_min, metric_new, metric_old;
    double **CP = new double *[3];
    double **gradient = new double *[3];
    double **Ab = new double *[3];
    double **b = new double *[3];
    for (i = 0; i < 3; i++) {
        CP[i] = new double [_x*_y*_z];
        gradient[i] = new double [_x*_y*_z];
        Ab[i] = new double [_x*_y*_z];
        b[i] = new double [_x*_y*_z];
    }

    stepsize_max = stepsize = 2;
    stepsize_min = 0.01;

    // initialise 3D FFD
    double x1, y1, z1, x2, y2, z2;
    x1 = y1 = z1 = 0;
    this->LatticeToWorld(x1, y1, z1);
    x2 = _x-1; y2 = _y-1; z2 = _z-1;
    this->LatticeToWorld(x2, y2, z2);
    irtkBSplineFreeFormTransformation3D * transform = new irtkBSplineFreeFormTransformation3D(x1, y1, z1, x2, y2, z2, _dx, _dy, _dz, _xaxis, _yaxis, _zaxis);

    // periodic model:
    if (_periodic) {
        while (time < 0)
            time += double(_t-1);
        while (time >= _t-1)
            time -= double(_t-1);
    } else {
        if ((time < -2) || (time > _t+1)) {
            cerr << "time out of bounds" << endl;
            exit(1);
        }
    }

    // calculate 3D FFD conversion matrix (same for all 3 displacement components (x, y, z))
    cout<<"fill matrix (A)"<<endl;
    irtkMatrix BS = irtkMatrix(_x*_y*_z, _x*_y*_z);
    // loops over all rows
    for (i = 0; i < _x; i++) {
        for (j = 0; j < _y; j++) {
            for (k = 0; k < _z; k++) {

                // loops over all columns in current row
                for (iC = 0; iC < _x; iC++) {
                    for (jC = 0; jC < _y; jC++) {
                        for (kC = 0; kC < _z; kC++) {
                            BS((k*_x*_y + j*_x + i), (kC*_x*_y + jC*_x + iC)) = 0;

                            // loops over all BSplines
                            for (m = -1; m < 3; m++) {
                                for (n = -1; n < 3; n++) {
                                    for (o = -1; o < 3; o++) {
                                        if (i+o == iC && j+n == jC && k+m == kC) {
                                            B_K = this->LookupTable[0][m];
                                            B_J = this->LookupTable[0][n] * B_K;
                                            B_I = this->LookupTable[0][o] * B_J;

                                            BS((k*_x*_y + j*_x + i), (kC*_x*_y + jC*_x + iC)) = B_I;
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
    }

    // fill in diplacement vector
    cout<<"initialize gradient components (A*A), (A*b), (A) and (b)"<<endl;
    double xDisp, yDisp, zDisp;
    for (i = 0; i < _x; i++) {
        for (j = 0; j < _y; j++) {
            for (k = 0; k < _z; k++) {
                xDisp = i;
                yDisp = j;
                zDisp = k;
                this->LocalDisplacement(xDisp, yDisp, zDisp, time);
                gradient[0][k*_x*_y + j*_x + i] = b[0][k*_x*_y + j*_x + i] = xDisp;
                gradient[1][k*_x*_y + j*_x + i] = b[1][k*_x*_y + j*_x + i] = yDisp;
                gradient[2][k*_x*_y + j*_x + i] = b[2][k*_x*_y + j*_x + i] = zDisp;
            }
        }
    }
    // and multiply it with B-spline matrix (A)
    irtkMatrix BSTranspose = ~BS;
    for (i = 0; i < _x; i++) {
        for (j = 0; j < _y; j++) {
            for (k = 0; k < _z; k++) {
                Ab[0][k*_x*_y + j*_x + i] = 0;
                Ab[1][k*_x*_y + j*_x + i] = 0;
                Ab[2][k*_x*_y + j*_x + i] = 0;
                for (iC = 0; iC < _x; iC++) {
                    for (jC = 0; jC < _y; jC++) {
                        for (kC = 0; kC < _z; kC++) {
                            for (a = 0; a < 3; a++) {
                                Ab[a][k*_x*_y + j*_x + i] += BSTranspose((k*_x*_y + j*_x + i), (kC*_x*_y + jC*_x + iC))*gradient[a][kC*_x*_y + jC*_x + iC];
                            }
                        }
                    }
                }
            }
        }
    }
    // build matrix (A*A)
    BS = BSTranspose*BS;
    BSTranspose.Transpose();

    // get initialisation for new control point values
    cout<<"get initialisation for new control point values"<<endl;
    delete transform;
    transform = this->Compute3DFFD1(time);
    for (i = 0; i < _x; i++) {
        for (j = 0; j < _y; j++) {
            for (k = 0; k < _z; k++) {
                transform->Get(i, j, k, CP[0][k*_x*_y + j*_x + i], CP[1][k*_x*_y + j*_x + i], CP[2][k*_x*_y + j*_x + i]);
            }
        }
    }
    // build metric of initialisation
    metric_new = 0;
    for (i = 0; i < _x; i++) {
        for (j = 0; j < _y; j++) {
            for (k = 0; k < _z; k++) {
                for (a = 0; a < 3; a++) {
                    double tmp = 0;
                    for (iC = 0; iC < _x; iC++) {
                        for (jC = 0; jC < _y; jC++) {
                            for (kC = 0; kC < _z; kC++) {
                                tmp += BSTranspose((k*_x*_y + j*_x + i), (kC*_x*_y + jC*_x + iC))*CP[a][kC*_x*_y + jC*_x + iC];
                            }
                        }
                    }
                    metric_new += (tmp - b[a][k*_x*_y + j*_x + i]) * (tmp - b[a][k*_x*_y + j*_x + i]);
                }
            }
        }
    }
    metric_old = metric_new;

    // compute new control point values
    // x(n+1) = x(n) - 2*A*(A*x(n) - b)
    cout<<"compute new control point values"<<endl;
    delta = 1;
    do {

        // build gradient 2*(A*A-A*b)
        for (i = 0; i < _x; i++) {
            for (j = 0; j < _y; j++) {
                for (k = 0; k < _z; k++) {
                    gradient[0][k*_x*_y + j*_x + i] = 0;
                    gradient[1][k*_x*_y + j*_x + i] = 0;
                    gradient[2][k*_x*_y + j*_x + i] = 0;
                    for (iC = 0; iC < _x; iC++) {
                        for (jC = 0; jC < _y; jC++) {
                            for (kC = 0; kC < _z; kC++) {
                                for (a = 0; a < 3; a++) {
                                    gradient[a][k*_x*_y + j*_x + i] -= BS((k*_x*_y + j*_x + i), (kC*_x*_y + jC*_x + iC))*CP[a][kC*_x*_y + jC*_x + iC];
                                }
                            }
                        }
                    }
                    gradient[0][k*_x*_y + j*_x + i] -= Ab[0][k*_x*_y + j*_x + i];
                    gradient[1][k*_x*_y + j*_x + i] -= Ab[1][k*_x*_y + j*_x + i];
                    gradient[2][k*_x*_y + j*_x + i] -= Ab[2][k*_x*_y + j*_x + i];
                    gradient[0][k*_x*_y + j*_x + i] *= 2;
                    gradient[1][k*_x*_y + j*_x + i] *= 2;
                    gradient[2][k*_x*_y + j*_x + i] *= 2;
                }
            }
        }

        // update x(n+1)
        for (i = 0; i < _x; i++) {
            for (j = 0; j < _y; j++) {
                for (k = 0; k < _z; k++) {
                    for (a = 0; a < 3; a++) {
                        CP[a][k*_x*_y + j*_x + i] -= stepsize*gradient[a][k*_x*_y + j*_x + i];
                    }
                }
            }
        }

        // build metric
        metric_new = 0;
        for (i = 0; i < _x; i++) {
            for (j = 0; j < _y; j++) {
                for (k = 0; k < _z; k++) {
                    for (a = 0; a < 3; a++) {
                        double tmp = 0;
                        for (iC = 0; iC < _x; iC++) {
                            for (jC = 0; jC < _y; jC++) {
                                for (kC = 0; kC < _z; kC++) {
                                    tmp += BSTranspose((k*_x*_y + j*_x + i), (kC*_x*_y + jC*_x + iC))*CP[a][kC*_x*_y + jC*_x + iC];
                                }
                            }
                        }
                        metric_new += (tmp - b[a][k*_x*_y + j*_x + i]) * (tmp - b[a][k*_x*_y + j*_x + i]);
                    }
                }
            }
        }

        if (metric_new < metric_old) {
            delta = metric_old - metric_new;
            metric_old = metric_new;
            stepsize = stepsize_max;
        } else {
            for (i = 0; i < _x; i++) {
                for (j = 0; j < _y; j++) {
                    for (k = 0; k < _z; k++) {
                        for (a = 0; a < 3; a++) {
                            CP[a][k*_x*_y + j*_x + i] += stepsize*gradient[a][k*_x*_y + j*_x + i];
                        }
                    }
                }
            }
            stepsize *= 0.5;
            if (stepsize < stepsize_min)
                break;
        }

    } while (delta < 0.0001);

    // write to transformation
    cout<<"write new transformation"<<endl;
    for (i = 0; i < _x; i++) {
        for (j = 0; j < _y; j++) {
            for (k = 0; k < _z; k++) {
                transform->Put(i, j, k, CP[0][k*_x*_y + j*_x + i], CP[1][k*_x*_y + j*_x + i], CP[2][k*_x*_y + j*_x + i]);
            }
        }
    }

    return transform;
}

