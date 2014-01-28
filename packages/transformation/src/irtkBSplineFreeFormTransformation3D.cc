/*=========================================================================

Library   : Image Registration Toolkit (IRTK)
Module    : $Id$
Copyright : Imperial College, Department of Computing
Visual Information Processing (VIP), 2008 onwards
Date      : $Date$
Version   : $Revision$
Changes   : $Author$

Copyright (c) IXICO LIMITED
All rights reserved.
See COPYRIGHT for details

=========================================================================*/

#include <irtkTransformation.h>

#define LUTSIZE (double)(FFDLOOKUPTABLESIZE-1)

double irtkBSplineFreeFormTransformation3D::LookupTable   [FFDLOOKUPTABLESIZE][4];

double irtkBSplineFreeFormTransformation3D::LookupTable_I [FFDLOOKUPTABLESIZE][4];

double irtkBSplineFreeFormTransformation3D::LookupTable_II[FFDLOOKUPTABLESIZE][4];

irtkBSplineFreeFormTransformation3D::irtkBSplineFreeFormTransformation3D()
{
	int i;

	// Initialize control point domain
	_origin._x = 0.5;
	_origin._y = 0.5;
	_origin._z = 0.5;

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

	// Initialize control point spacing
	_dx = 1;
	_dy = 1;
	_dz = 1;

	// Initialize transformation matrix
	_matL2W = irtkMatrix(4, 4);
	_matW2L = irtkMatrix(4, 4);

	// Update transformation matrix
	this->UpdateMatrix();

	// Intialize memory for control point values
	_data = this->Allocate(_data, _x, _y, _z);

	// Initialize memory for control point status
	_status = new _Status[3*_x*_y*_z];
	if (_z == 1) {
		for (i = 0; i < 2*_x*_y*_z; i++) {
			_status[i] = _Active;
		}
		for (i = 2*_x*_y*_z; i < 3*_x*_y*_z; i++) {
			_status[i] = _Passive;
		}
	} else {
		for (i = 0; i < 3*_x*_y*_z; i++) {
			_status[i] = _Active;
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

irtkBSplineFreeFormTransformation3D::irtkBSplineFreeFormTransformation3D(irtkGenericImage<double> &image)
{
	int i, j, k;
	double x, y, z;

	// Get out the image attributes
	irtkImageAttributes attr = image.GetImageAttributes();

	// and copy
	_x         = attr._x;
	_y         = attr._y;
	_z         = attr._z;
	_dx        = attr._dx;
	_dy        = attr._dy;
	_dz        = attr._dz;
	_origin._x = attr._xorigin;
	_origin._y = attr._yorigin;
	_origin._z = attr._zorigin;

	// Copy axis
	image.GetOrientation(_xaxis, _yaxis, _zaxis);

	// Initialize transformation matrix
	_matL2W = irtkMatrix(4, 4);
	_matW2L = irtkMatrix(4, 4);

	// Update transformation matrix
	this->UpdateMatrix();

	// Intialize memory for control point values
	_data = this->Allocate(_data, _x, _y, _z);
	for (i = 0; i < _x; i++) {
		for (j = 0; j < _y; j++) {
			for (k = 0; k < _z; k++) {
				x = i;
				y = j;
				z = k;
				image.ImageToWorld(x, y, z);
				_data[k][j][i]._x = image(i, j, k, 0) - x;
				_data[k][j][i]._y = image(i, j, k, 1) - y;
				_data[k][j][i]._z = image(i, j, k, 2) - z;
			}
		}
	}

	// Initialize memory for control point status
	_status = new _Status[3*_x*_y*_z];
	if (_z == 1) {
		for (i = 0; i < 2*_x*_y*_z; i++) {
			_status[i] = _Active;
		}
		for (i = 2*_x*_y*_z; i < 3*_x*_y*_z; i++) {
			_status[i] = _Passive;
		}
	} else {
		for (i = 0; i < 3*_x*_y*_z; i++) {
			_status[i] = _Active;
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

irtkBSplineFreeFormTransformation3D::irtkBSplineFreeFormTransformation3D(irtkImageAttributes &attr, double dx, double dy, double dz)
{
    irtkMatrix matI2W = irtkBaseImage::GetImageToWorldMatrix(attr);
    irtkMatrix matW2I = irtkBaseImage::GetWorldToImageMatrix(attr);

    // Figure out FOV
    double x1 = matI2W(0, 3);
    double y1 = matI2W(1, 3);
    double z1 = matI2W(2, 3);
    double x2 = matI2W(0, 0)*(attr._x-1)+matI2W(0, 1)*(attr._y-1)+matI2W(0, 2)*(attr._z-1)+matI2W(0, 3);
    double y2 = matI2W(1, 0)*(attr._x-1)+matI2W(1, 1)*(attr._y-1)+matI2W(1, 2)*(attr._z-1)+matI2W(1, 3);
    double z2 = matI2W(2, 0)*(attr._x-1)+matI2W(2, 1)*(attr._y-1)+matI2W(2, 2)*(attr._z-1)+matI2W(2, 3);

    // Initialize control point domain
    _origin._x = (x2 + x1) / 2.0;
    _origin._y = (y2 + y1) / 2.0;
    _origin._z = (z2 + z1) / 2.0;

    // Initialize x-axis and y-axis
    _xaxis[0] = attr._xaxis[0];
    _xaxis[1] = attr._xaxis[1];
    _xaxis[2] = attr._xaxis[2];
    _yaxis[0] = attr._yaxis[0];
    _yaxis[1] = attr._yaxis[1];
    _yaxis[2] = attr._yaxis[2];
    _zaxis[0] = attr._zaxis[0];
    _zaxis[1] = attr._zaxis[1];
    _zaxis[2] = attr._zaxis[2];

    double a, b, c;

    a = x1 * _xaxis[0] + y1 * _xaxis[1] + z1 * _xaxis[2];
    b = x1 * _yaxis[0] + y1 * _yaxis[1] + z1 * _yaxis[2];
    c = x1 * _zaxis[0] + y1 * _zaxis[1] + z1 * _zaxis[2];
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

    // Initialize control point spacing
    _dx = (x2 - x1) / (_x - 1);
    _dy = (y2 - y1) / (_y - 1);
    _dz = _z > 1 ? (z2 - z1) / (_z - 1) : dz;

    // Initialize transformation matrix
    _matL2W = irtkMatrix(4, 4);
    _matW2L = irtkMatrix(4, 4);

    // Update transformation matrix
    this->UpdateMatrix();

    // Initialize memory for control point values
    _data = this->Allocate(_data, _x, _y, _z);

    // Initialize memory for control point status
    _status = new _Status[3*_x*_y*_z];
    for (int i = 0; i < 3*_x*_y*_z; i++) {
      _status[i] = _Active;
    }

    // Initialize lookup table
    for (int i = 0; i < FFDLOOKUPTABLESIZE; i++) {
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

irtkBSplineFreeFormTransformation3D::irtkBSplineFreeFormTransformation3D(irtkBaseImage &image, double dx, double dy, double dz)
{
    int i;
    double x1, y1, z1, x2, y2, z2;

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
    if (image.GetZ() > 1) {
        // Check if spacing is too large in any direction
        if (_x < 2) {
            cerr << "irtkBSplineFreeFormTransformation3D::irtkBSplineFreeFormTransformation3D: Requested control point spacing in the x-direction too large" << endl;
            exit(1);
        }
        if (_y < 2) {
            cerr << "irtkBSplineFreeFormTransformation3D::irtkBSplineFreeFormTransformation3D: Requested control point spacing in the y-direction too large" << endl;
            exit(1);
        }
        if (_z < 2) {
            cerr << "irtkBSplineFreeFormTransformation3D::irtkBSplineFreeFormTransformation3D: Requested control point spacing in the z-direction too large" << endl;
            exit(1);
        }
    } else {
        // Check if spacing is too large in any direction
        if (_x < 2) {
            cerr << "irtkBSplineFreeFormTransformation3D::irtkBSplineFreeFormTransformation3D: Requested control point spacing in the x-direction too large" << endl;
            exit(1);
        }
        if (_y < 2) {
            cerr << "irtkBSplineFreeFormTransformation3D::irtkBSplineFreeFormTransformation3D: Requested control point spacing in the y-direction too large" << endl;
            exit(1);
        }
        _z = 1;
    }

    // Initialize control point spacing
    _dx = (x2 - x1) / (_x - 1);
    _dy = (y2 - y1) / (_y - 1);
    if (_z > 1) {
        _dz = (z2 - z1) / (_z - 1);
    } else {
        _dz = dz;
    }

    // Initialize transformation matrix
    _matL2W = irtkMatrix(4, 4);
    _matW2L = irtkMatrix(4, 4);

    // Update transformation matrix
    this->UpdateMatrix();

    // Intialize memory for control point values
    _data = this->Allocate(_data, _x, _y, _z);

    // Initialize memory for control point status
    _status = new _Status[3*_x*_y*_z];
    if (_z == 1) {
        for (i = 0; i < 2*_x*_y*_z; i++) {
            _status[i] = _Active;
        }
        for (i = 2*_x*_y*_z; i < 3*_x*_y*_z; i++) {
            _status[i] = _Passive;
        }
    } else {
        for (i = 0; i < 3*_x*_y*_z; i++) {
            _status[i] = _Active;
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

irtkBSplineFreeFormTransformation3D::irtkBSplineFreeFormTransformation3D(double x1, double y1, double z1, double x2, double y2, double z2, double dx, double dy, double dz, double* xaxis, double* yaxis, double* zaxis)
{
	int i;

	// Initialize control point domain
	_origin._x = (x2 + x1) / 2.0;
	_origin._y = (y2 + y1) / 2.0;
	_origin._z = (z2 + z1) / 2.0;

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
	if (z2 > z1) {
		_z = round((z2 - z1) / dz) + 1;
		// Check if spacing is too large in any direction
		if (_x < 2) {
			cerr << "irtkBSplineFreeFormTransformation3D::irtkBSplineFreeFormTransformation3D: Requested control point spacing in the x-direction too large" << endl;
			exit(1);
		}
		if (_y < 2) {
			cerr << "irtkBSplineFreeFormTransformation3D::irtkBSplineFreeFormTransformation3D: Requested control point spacing in the y-direction too large" << endl;
			exit(1);
		}
		if (_z < 2) {
			cerr << "irtkBSplineFreeFormTransformation3D::irtkBSplineFreeFormTransformation3D: Requested control point spacing in the z-direction too large" << endl;
			exit(1);
		}
	} else {
		// Check if spacing is too large in any direction
		if (_x < 2) {
			cerr << "irtkBSplineFreeFormTransformation3D::irtkBSplineFreeFormTransformation3D: Requested control point spacing in the x-direction too large" << endl;
			exit(1);
		}
		if (_y < 2) {
			cerr << "irtkBSplineFreeFormTransformation3D::irtkBSplineFreeFormTransformation3D: Requested control point spacing in the y-direction too large" << endl;
			exit(1);
		}
		_z = 1;
	}

	// Initialize control point spacing
	_dx = (x2 - x1) / (_x - 1);
	_dy = (y2 - y1) / (_y - 1);
	if (z2 > z1) {
		_dz = (z2 - z1) / (_z - 1);
	} else {
		_dz = 1;
	}

	// Initialize transformation matrix
	_matL2W = irtkMatrix(4, 4);
	_matW2L = irtkMatrix(4, 4);

	// Update transformation matrix
	this->UpdateMatrix();

	// Intialize memory for control point values
	_data = this->Allocate(_data, _x, _y, _z);

	// Initialize memory for control point status
	_status = new _Status[3*_x*_y*_z];
	if (_z == 1) {
		for (i = 0; i < 2*_x*_y*_z; i++) {
			_status[i] = _Active;
		}
		for (i = 2*_x*_y*_z; i < 3*_x*_y*_z; i++) {
			_status[i] = _Passive;
		}
	} else {
		for (i = 0; i < 3*_x*_y*_z; i++) {
			_status[i] = _Active;
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

irtkBSplineFreeFormTransformation3D::irtkBSplineFreeFormTransformation3D(const irtkBSplineFreeFormTransformation3D &ffd) : irtkFreeFormTransformation3D(ffd)
{
	int i, j, k;

	// Initialize origin
	_origin = ffd._origin;

	// Initialize control point dimensions
	_x = ffd._x;
	_y = ffd._y;
	_z = ffd._z;

	// Intialize control point spacing
	_dx = ffd._dx;
	_dy = ffd._dy;
	_dz = ffd._dz;

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
	_data = this->Allocate(_data, _x, _y, _z);
	for (i = -4; i < _x+4; i++) {
		for (j = -4; j < _y+4; j++) {
			for (k = -4; k < _z+4; k++) {
				_data[k][j][i]._x = ffd._data[k][j][i]._x;
				_data[k][j][i]._y = ffd._data[k][j][i]._y;
				_data[k][j][i]._z = ffd._data[k][j][i]._z;
			}
		}
	}

	// Initialize memory for control point status
	_status = new _Status[3*_x*_y*_z];
	for (i = 0; i < 3*_x*_y*_z; i++) {
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

irtkBSplineFreeFormTransformation3D::~irtkBSplineFreeFormTransformation3D()
{
	// Free memory for control points if necessary
	if (_data != NULL) _data = this->Deallocate(_data, _x, _y, _z);

	_x = 0;
	_y = 0;
	_z = 0;
}

void irtkBSplineFreeFormTransformation3D::FFD2D(double &x, double &y) const
{
	// Check if there is some work to do
	if ((x < -2) || (y < -2) || (x > _x+1) || (y > _y+1)) {
		x = 0;
		y = 0;
		return;
	}

	irtkVector3D<double> *data;
	double s, t, B_J, xii, yii;
	int i, j, l, m, S, T;

	// Now calculate the real stuff
	l = (int)floor(x);
	m = (int)floor(y);
	s = x-l;
	t = y-m;

	// Calculate offset
	i = (_x + 8) * (_y + 4);
	x = 0;
	y = 0;
	S = round(LUTSIZE*s);
	T = round(LUTSIZE*t);
	data = &(_data[0][m-1][l-1]);
	for (j = 0; j < 4; j++) {
		B_J = this->LookupTable[T][j];

		// Inner most loop unrolled starts here
		xii  = data->_x * this->LookupTable[S][0];
		yii  = data->_y * this->LookupTable[S][0];
		data++;
		xii += data->_x * this->LookupTable[S][1];
		yii += data->_y * this->LookupTable[S][1];
		data++;
		xii += data->_x * this->LookupTable[S][2];
		yii += data->_y * this->LookupTable[S][2];
		data++;
		xii += data->_x * this->LookupTable[S][3];
		yii += data->_y * this->LookupTable[S][3];
		data++;
		// Inner most loop unrolled stops here

		x += xii * B_J;
		y += yii * B_J;
		data += _x + 4;
	}
}

void irtkBSplineFreeFormTransformation3D::FFD3D(double &x, double &y, double &z) const
{
	// Check if there is some work to do
	if ((x < -2) || (y < -2) || (z < -2) || (x > _x+1) || (y > _y+1) || (z > _z+1)) {
		x = 0;
		y = 0;
		z = 0;
		return;
	}

	irtkVector3D<double> *data;
	double s, t, u, B_J, B_K, xi, yi, zi, xii, yii, zii;
	int i, j, k, l, m, n, S, T, U;

	// Now calculate the real stuff
	l = (int)floor(x);
	m = (int)floor(y);
	n = (int)floor(z);
	s = x-l;
	t = y-m;
	u = z-n;

	// Calculate offset
	i = (_x + 8) * (_y + 4);
	x = 0;
	y = 0;
	z = 0;
	S = round(LUTSIZE*s);
	T = round(LUTSIZE*t);
	U = round(LUTSIZE*u);
	data = &(_data[n-1][m-1][l-1]);
	for (k = 0; k < 4; k++) {
		B_K = this->LookupTable[U][k];
		xi = 0;
		yi = 0;
		zi = 0;
		for (j = 0; j < 4; j++) {
			B_J = this->LookupTable[T][j];
			xii = 0; yii = 0; zii = 0;

			// Inner most loop unrolled starts here
			xii  += data->_x * this->LookupTable[S][0];
			yii  += data->_y * this->LookupTable[S][0];
			zii  += data->_z * this->LookupTable[S][0];
			data++;
			xii  += data->_x * this->LookupTable[S][1];
			yii  += data->_y * this->LookupTable[S][1];
			zii  += data->_z * this->LookupTable[S][1];
			data++;
			xii  += data->_x * this->LookupTable[S][2];
			yii  += data->_y * this->LookupTable[S][2];
			zii  += data->_z * this->LookupTable[S][2];
			data++;
			xii  += data->_x * this->LookupTable[S][3];
			yii  += data->_y * this->LookupTable[S][3];
			zii  += data->_z * this->LookupTable[S][3];
			data++;
			// Inner most loop unrolled stops here

			xi += xii * B_J;
			yi += yii * B_J;
			zi += zii * B_J;
			data += _x + 4;
		}
		x += xi * B_K;
		y += yi * B_K;
		z += zi * B_K;
		data += i;
	}
}

void irtkBSplineFreeFormTransformation3D::Jacobian(irtkMatrix &jac, double x, double y, double z, double t)
{
	this->LocalJacobian(jac, x, y, z, t);
}

void irtkBSplineFreeFormTransformation3D::GlobalJacobian(irtkMatrix &jac, double, double, double, double)
{
	// Jacobian matrix is 3 x 3
	jac.Initialize(3, 3);

	// Set matrix to identity
	jac(0, 0) = 1;
	jac(1, 1) = 1;
	jac(2, 2) = 1;
}

void irtkBSplineFreeFormTransformation3D::LocalJacobian(irtkMatrix &jac, double x, double y, double z, double)
{
	int i, j, k, l, m, n, I, J, K, S, T, U;
	double s, t, u, v, B_K, B_J, B_I, B_K_I, B_J_I, B_I_I;
	double z_k=0, y_k=0, x_k=0, z_j=0, y_j=0, x_j=0, z_i=0, y_i=0, x_i=0;

	// Convert to lattice coordinates
	this->WorldToLattice(x, y, z);

	// Compute derivatives
	l = (int)floor(x);
	m = (int)floor(y);
	n = (int)floor(z);
	s = x-l;
	t = y-m;
	u = z-n;
	S = round(LUTSIZE*s);
	T = round(LUTSIZE*t);
	U = round(LUTSIZE*u);
	for (k = 0; k < 4; k++) {
		K = k + n - 1;
		if ((K >= 0) && (K < _z)) {
			B_K   = this->LookupTable[U][k];
			B_K_I = this->LookupTable_I[U][k];
			for (j = 0; j < 4; j++) {
				J = j + m - 1;
				if ((J >= 0) && (J < _y)) {
					B_J   = this->LookupTable[T][j];
					B_J_I = this->LookupTable_I[T][j];
					for (i = 0; i < 4; i++) {
						I = i + l - 1;
						if ((I >= 0) && (I < _x)) {
							B_I   = this->LookupTable[S][i];
							B_I_I = this->LookupTable_I[S][i];
							v = B_I * B_J * B_K_I;
							x_k += _data[K][J][I]._x * v;
							y_k += _data[K][J][I]._y * v;
							z_k += _data[K][J][I]._z * v;
							v = B_I * B_J_I * B_K;
							x_j += _data[K][J][I]._x * v;
							y_j += _data[K][J][I]._y * v;
							z_j += _data[K][J][I]._z * v;
							v = B_I_I * B_J * B_K;
							x_i += _data[K][J][I]._x * v;
							y_i += _data[K][J][I]._y * v;
							z_i += _data[K][J][I]._z * v;
						}
					}
				}
			}
		}
	}

	// Jacobian matrix is 3 x 3
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

void irtkBSplineFreeFormTransformation3D::JacobianDetDerivative(irtkMatrix *detdev, int x, int y, int z)
{
	int i;
	double x_b,y_b,z_b,x_f,y_f,z_f,b_i,b_j,b_k;

	switch(x) {
	case -1:
		x_b=(1.0/6.0);
		x_f=(0.5);
		break;
	case 0:
		x_b=(2.0/3.0);
		x_f=(0.0);
		break;
	case 1:
		x_b=(1.0/6.0);
		x_f=(-0.5);
		break;
	default:
		x_b=0.0;
		x_f=0.0;
		break;
	}

	switch(y) {
	case -1:
		y_b=(1.0/6.0);
		y_f=(0.5);
		break;
	case 0:
		y_b=(2.0/3.0);
		y_f=(0.0);
		break;
	case 1:
		y_b=(1.0/6.0);
		y_f=(-0.5);
		break;
	default:
		y_b=0.0;
		y_f=0.0;
		break;
	}

	switch(z) {
	case -1:
		z_b=(1.0/6.0);
		z_f=(0.5);
		break;
	case 0:
		z_b=(2.0/3.0);
		z_f=(0.0);
		break;
	case 1:
		z_b=(1.0/6.0);
		z_f=(-0.5);
		break;
	default:
		z_b=0.0;
		z_f=0.0;
		break;
	}

	for(i = 0; i < 3; i++) {
		detdev[i].Initialize(3,3);
	}

	b_i = x_f*y_b*z_b;
	b_j = x_b*y_f*z_b;
	b_k = x_b*y_b*z_f;

	for(i = 0; i < 3; i++) {
		// with respect to ui
		detdev[i](i,0) = b_i; detdev[i](i,1) = b_j; detdev[i](i,2) = b_k;
		detdev[i] = detdev[i] * _matW2L(0, 0, 3, 3);
	}

}

void irtkBSplineFreeFormTransformation3D::JacobianDOFs(double jac[3], int dof, double x, double y, double z, double)
{
	int i, j, k;

	// Transforms index of control point location
	this->IndexToLattice(dof, i, j, k);

	// Compute lattice coordinates
	this->WorldToLattice(x, y, z);

	jac[0] = this->B(x - i) * this->B(y - j) * this->B(z - k);
	jac[1] = jac[0];
	jac[2] = jac[0];
}

double irtkBSplineFreeFormTransformation3D::Bending2D(int i, int j)
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
			y_ii += _data[0][J][I]._y * v;
			x_ii += _data[0][J][I]._x * v;

			v = B_I * B_J_II;
			y_jj += _data[0][J][I]._y * v;
			x_jj += _data[0][J][I]._x * v;

			v = B_I_I * B_J_I;
			y_ij += _data[0][J][I]._y * v;
			x_ij += _data[0][J][I]._x * v;
		}
	}

	// Compute bending
	return (x_ii*x_ii + x_jj*x_jj + y_ii*y_ii + y_jj*y_jj + 2*(x_ij*x_ij + y_ij*y_ij));
}

double irtkBSplineFreeFormTransformation3D::Bending3D(int i, int j, int k)
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
				z_kk += _data[K][J][I]._z * v;
				y_kk += _data[K][J][I]._y * v;
				x_kk += _data[K][J][I]._x * v;

				v = B_I * B_J_II * B_K;
				z_jj += _data[K][J][I]._z * v;
				y_jj += _data[K][J][I]._y * v;
				x_jj += _data[K][J][I]._x * v;

				v = B_I_II * B_J * B_K;
				z_ii += _data[K][J][I]._z * v;
				y_ii += _data[K][J][I]._y * v;
				x_ii += _data[K][J][I]._x * v;

				v = B_I_I * B_J_I * B_K;
				z_ij += _data[K][J][I]._z * v;
				y_ij += _data[K][J][I]._y * v;
				x_ij += _data[K][J][I]._x * v;

				v = B_I_I * B_J * B_K_I;
				z_ik += _data[K][J][I]._z * v;
				y_ik += _data[K][J][I]._y * v;
				x_ik += _data[K][J][I]._x * v;

				v = B_I * B_J_I * B_K_I;
				z_jk += _data[K][J][I]._z * v;
				y_jk += _data[K][J][I]._y * v;
				x_jk += _data[K][J][I]._x * v;
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

double irtkBSplineFreeFormTransformation3D::Bending()
{
	int i, j, k;
	double bending;

	bending = 0;
	if (_z == 1) {
		for (j = 0; j < _y; j++) {
			for (i = 0; i < _x; i++) {
				bending += this->Bending2D(i, j);
			}
		}
	} else {
		for (k = 0; k < _z; k++) {
			for (j = 0; j < _y; j++) {
				for (i = 0; i < _x; i++) {
					bending += this->Bending3D(i, j, k);
				}
			}
		}
	}
	return bending;
}

void irtkBSplineFreeFormTransformation3D::BendingGradient2D(double *gradient)
{
	int I, J, index, i, j, k, n = _x*_y*_z;
	double B_J, B_I, B_J_I, B_I_I, B_J_II, B_I_II, v, tmp[3];

	// Derivatives
	irtkVector3D<double> ***d_jj = NULL;
	d_jj = this->Allocate(d_jj, _x, _y, _z);
	irtkVector3D<double> ***d_ii = NULL;
	d_ii = this->Allocate(d_ii, _x, _y, _z);
	irtkVector3D<double> ***d_ij = NULL;
	d_ij = this->Allocate(d_ij, _x, _y, _z);

	// Values of the B-spline basis functions and its derivative (assuming that we compute the bending energy only at the control point location i, j, k)
	double b[3]    = {1.0/6.0, 2.0/3.0, 1.0/6.0};
	double b_i[3]  = {-0.5, 0, 0.5};
	double b_ii[3] = {1.0, -2.0, 1.0};

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
						d_jj[k][j][i]._y += 2 * _data[k][J][I]._y * v;
						d_jj[k][j][i]._x += 2 * _data[k][J][I]._x * v;

						v = B_I_II * B_J;
						d_ii[k][j][i]._y += 2 * _data[k][J][I]._y * v;
						d_ii[k][j][i]._x += 2 * _data[k][J][I]._x * v;

						v = B_I_I * B_J_I;
						d_ij[k][j][i]._y += 4 * _data[k][J][I]._y * v;
						d_ij[k][j][i]._x += 4 * _data[k][J][I]._x * v;

					}
				}
			}
		}
	}

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
						tmp[1] += d_jj[k][J][I]._y * v;
						tmp[0] += d_jj[k][J][I]._x * v;

						v = B_I_II * B_J;
						tmp[1] += d_ii[k][J][I]._y * v;
						tmp[0] += d_ii[k][J][I]._x * v;

						v = B_I_I * B_J_I;
						tmp[1] += d_ij[k][J][I]._y * v;
						tmp[0] += d_ij[k][J][I]._x * v;

					}
				}
				index = this->LatticeToIndex(i, j, k);
				gradient[index]     += -tmp[0];
				gradient[index+n]   += -tmp[1];
				gradient[index+2*n] += 0;
			}
		}
	}

	this->Deallocate(d_jj, _x, _y, _z);
	this->Deallocate(d_ii, _x, _y, _z);
	this->Deallocate(d_ij, _x, _y, _z);
}

void irtkBSplineFreeFormTransformation3D::BendingGradient3D(double *gradient)
{
	int I, J, K, index, i, j, k, n = _x*_y*_z;
	double B_K, B_J, B_I, B_K_I, B_J_I, B_I_I, B_K_II, B_J_II, B_I_II, v, tmp[3];

	// Derivatives
	irtkVector3D<double> ***d_kk = NULL;
	d_kk = this->Allocate(d_kk, _x, _y, _z);
	irtkVector3D<double> ***d_jj = NULL;
	d_jj = this->Allocate(d_jj, _x, _y, _z);
	irtkVector3D<double> ***d_ii = NULL;
	d_ii = this->Allocate(d_ii, _x, _y, _z);
	irtkVector3D<double> ***d_ij = NULL;
	d_ij = this->Allocate(d_ij, _x, _y, _z);
	irtkVector3D<double> ***d_ik = NULL;
	d_ik = this->Allocate(d_ik, _x, _y, _z);
	irtkVector3D<double> ***d_jk = NULL;
	d_jk = this->Allocate(d_jk, _x, _y, _z);

	// Values of the B-spline basis functions and its derivative (assuming that we compute the bending energy only at the control point location i, j, k)
	double b[3]    = {1.0/6.0, 2.0/3.0, 1.0/6.0};
	double b_i[3]  = {-0.5, 0, 0.5};
	double b_ii[3] = {1.0, -2.0, 1.0};

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
							d_kk[k][j][i]._z += 2 * _data[K][J][I]._z * v;
							d_kk[k][j][i]._y += 2 * _data[K][J][I]._y * v;
							d_kk[k][j][i]._x += 2 * _data[K][J][I]._x * v;

							v = B_I * B_J_II * B_K;
							d_jj[k][j][i]._z += 2 * _data[K][J][I]._z * v;
							d_jj[k][j][i]._y += 2 * _data[K][J][I]._y * v;
							d_jj[k][j][i]._x += 2 * _data[K][J][I]._x * v;

							v = B_I_II * B_J * B_K;
							d_ii[k][j][i]._z += 2 * _data[K][J][I]._z * v;
							d_ii[k][j][i]._y += 2 * _data[K][J][I]._y * v;
							d_ii[k][j][i]._x += 2 * _data[K][J][I]._x * v;

							v = B_I_I * B_J_I * B_K;
							d_ij[k][j][i]._z += 4 * _data[K][J][I]._z * v;
							d_ij[k][j][i]._y += 4 * _data[K][J][I]._y * v;
							d_ij[k][j][i]._x += 4 * _data[K][J][I]._x * v;

							v = B_I_I * B_J * B_K_I;
							d_ik[k][j][i]._z += 4 * _data[K][J][I]._z * v;
							d_ik[k][j][i]._y += 4 * _data[K][J][I]._y * v;
							d_ik[k][j][i]._x += 4 * _data[K][J][I]._x * v;

							v = B_I * B_J_I * B_K_I;
							d_jk[k][j][i]._z += 4 * _data[K][J][I]._z * v;
							d_jk[k][j][i]._y += 4 * _data[K][J][I]._y * v;
							d_jk[k][j][i]._x += 4 * _data[K][J][I]._x * v;
						}
					}
				}
			}
		}
	}

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
							tmp[2] += d_kk[K][J][I]._z * v;
							tmp[1] += d_kk[K][J][I]._y * v;
							tmp[0] += d_kk[K][J][I]._x * v;

							v = B_I * B_J_II * B_K;
							tmp[2] += d_jj[K][J][I]._z * v;
							tmp[1] += d_jj[K][J][I]._y * v;
							tmp[0] += d_jj[K][J][I]._x * v;

							v = B_I_II * B_J * B_K;
							tmp[2] += d_ii[K][J][I]._z * v;
							tmp[1] += d_ii[K][J][I]._y * v;
							tmp[0] += d_ii[K][J][I]._x * v;

							v = B_I_I * B_J_I * B_K;
							tmp[2] += d_ij[K][J][I]._z * v;
							tmp[1] += d_ij[K][J][I]._y * v;
							tmp[0] += d_ij[K][J][I]._x * v;

							v = B_I_I * B_J * B_K_I;
							tmp[2] += d_ik[K][J][I]._z * v;
							tmp[1] += d_ik[K][J][I]._y * v;
							tmp[0] += d_ik[K][J][I]._x * v;

							v = B_I * B_J_I * B_K_I;
							tmp[2] += d_jk[K][J][I]._z * v;
							tmp[1] += d_jk[K][J][I]._y * v;
							tmp[0] += d_jk[K][J][I]._x * v;
						}
					}
				}
				index = this->LatticeToIndex(i, j, k);
				gradient[index]     += -tmp[0];
				gradient[index+n]   += -tmp[1];
				gradient[index+2*n] += -tmp[2];
			}
		}
	}

	this->Deallocate(d_kk, _x, _y, _z);
	this->Deallocate(d_jj, _x, _y, _z);
	this->Deallocate(d_ii, _x, _y, _z);
	this->Deallocate(d_ij, _x, _y, _z);
	this->Deallocate(d_ik, _x, _y, _z);
	this->Deallocate(d_jk, _x, _y, _z);
}

void irtkBSplineFreeFormTransformation3D::BendingGradient(double *gradient)
{
	if (_z == 1) {
		this->BendingGradient2D(gradient);
	} else {
		this->BendingGradient3D(gradient);
	}
}

double irtkBSplineFreeFormTransformation3D::Bending2D(double x, double y)
{
	int i, j, l, m, I, J, S, T;
	double s, t, v, z;
	double B_J, B_I, B_J_I, B_I_I, B_J_II, B_I_II;
	double y_jj=0, x_jj=0, y_ii=0, x_ii=0, y_ij=0, x_ij=0;

	z = 0;
	this->WorldToLattice(x, y, z);

	// Compute BSpline derivatives
	l = (int)floor(x);
	m = (int)floor(y);
	s = x-l;
	t = y-m;
	S = round(LUTSIZE*s);
	T = round(LUTSIZE*t);
	for (j = 0; j < 4; j++) {
		J = j + m - 1;
		if ((J >= 0) && (J < _y)) {
			B_J    = this->LookupTable[T][j];
			B_J_I  = this->LookupTable_I[T][j];
			B_J_II = this->LookupTable_II[T][j];
			for (i = 0; i < 4; i++) {
				I = i + l - 1;
				if ((I >= 0) && (I < _x)) {
					B_I    = this->LookupTable[S][i];
					B_I_I  = this->LookupTable_I[S][i];
					B_I_II = this->LookupTable_II[S][i];

					v = B_I * B_J_II;
					y_jj += _data[0][J][I]._y * v;
					x_jj += _data[0][J][I]._x * v;

					v = B_I_II * B_J;
					y_ii += _data[0][J][I]._y * v;
					x_ii += _data[0][J][I]._x * v;

					v = B_I_I * B_J_I;
					y_ij += _data[0][J][I]._y * v;
					x_ij += _data[0][J][I]._x * v;
				}
			}
		}
	}

	// Compute bending
	return (x_ii*x_ii + x_jj*x_jj + y_ii*y_ii + y_jj*y_jj + 2*(x_ij*x_ij + y_ij*y_ij));
}

double irtkBSplineFreeFormTransformation3D::Bending3D(double x, double y, double z)
{
	int i, j, k, l, m, n, I, J, K, S, T, U;
	double s, t, u, v;
	double B_K, B_J, B_I, B_K_I, B_J_I, B_I_I, B_K_II, B_J_II, B_I_II;
	double z_kk=0, y_kk=0, x_kk=0, z_jj=0, y_jj=0, x_jj=0, z_ii=0, y_ii=0, x_ii=0;
	double z_ij=0, y_ij=0, x_ij=0, z_ik=0, y_ik=0, x_ik=0, z_jk=0, y_jk=0, x_jk=0;

	this->WorldToLattice(x, y, z);

	// Compute BSpline derivatives
	l = (int)floor(x);
	m = (int)floor(y);
	n = (int)floor(z);
	s = x-l;
	t = y-m;
	u = z-n;
	S = round(LUTSIZE*s);
	T = round(LUTSIZE*t);
	U = round(LUTSIZE*u);
	for (k = 0; k < 4; k++) {
		K = k + n - 1;
		if ((K >= 0) && (K < _z)) {
			B_K    = this->LookupTable[U][k];
			B_K_I  = this->LookupTable_I[U][k];
			B_K_II = this->LookupTable_II[U][k];
			for (j = 0; j < 4; j++) {
				J = j + m - 1;
				if ((J >= 0) && (J < _y)) {
					B_J    = this->LookupTable[T][j];
					B_J_I  = this->LookupTable_I[T][j];
					B_J_II = this->LookupTable_II[T][j];
					for (i = 0; i < 4; i++) {
						I = i + l - 1;
						if ((I >= 0) && (I < _x)) {
							B_I    = this->LookupTable[S][i];
							B_I_I  = this->LookupTable_I[S][i];
							B_I_II = this->LookupTable_II[S][i];

							v = B_I * B_J * B_K_II;
							z_kk += _data[K][J][I]._z * v;
							y_kk += _data[K][J][I]._y * v;
							x_kk += _data[K][J][I]._x * v;

							v = B_I * B_J_II * B_K;
							z_jj += _data[K][J][I]._z * v;
							y_jj += _data[K][J][I]._y * v;
							x_jj += _data[K][J][I]._x * v;

							v = B_I_II * B_J * B_K;
							z_ii += _data[K][J][I]._z * v;
							y_ii += _data[K][J][I]._y * v;
							x_ii += _data[K][J][I]._x * v;

							v = B_I_I * B_J_I * B_K;
							z_ij += _data[K][J][I]._z * v;
							y_ij += _data[K][J][I]._y * v;
							x_ij += _data[K][J][I]._x * v;

							v = B_I_I * B_J * B_K_I;
							z_ik += _data[K][J][I]._z * v;
							y_ik += _data[K][J][I]._y * v;
							x_ik += _data[K][J][I]._x * v;

							v = B_I * B_J_I * B_K_I;
							z_jk += _data[K][J][I]._z * v;
							y_jk += _data[K][J][I]._y * v;
							x_jk += _data[K][J][I]._x * v;
						}
					}
				}
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

int irtkBSplineFreeFormTransformation3D::CheckHeader(char *name)
{
	char buffer[255];

	ifstream from(name);

	if (!from) {
		cerr << "irtkBSplineFreeFormTransformation3D::CheckHeader: Can't open file "
			<< name << "\n";
		exit(1);
	}

	// Read keyword
	from >> buffer;
	if ((strcmp(buffer, "AFFD:") != 0) && (strcmp(buffer, "AFFD_BINARY:") != 0)) {
		return false;
	}

	return true;
}

void irtkBSplineFreeFormTransformation3D::ApproximateAsNew2D(const double *x1, const double *y1, const double *z1, double *x2, double *y2, double *z2, int no)
{
	int i, j, l, m, I, J, S, T, index;
	double s, t, x, y, z, B_I, B_J, basis, basis2, phi, norm;

	// Allocate memory for control points
	irtkVector3D<double> ***data = NULL;
	data = Allocate(data, _x, _y, _z);

	// Allocate memory for temporary storage
	double ***ds = NULL;
	ds = Allocate(ds, _x, _y, _z);

	// Initialize data structures
	for (j = -2; j < _y+2; j++) {
		for (i = -2; i < _x+2; i++) {
			data[0][j][i]._x = 0;
			data[0][j][i]._y = 0;
			data[0][j][i]._z = 0;
			ds[0][j][i]      = 0;
		}
	}

	// Initial loop: Calculate change of control points
	for (index = 0; index < no; index++) {
		x = x1[index];
		y = y1[index];
		z = z1[index];
		this->WorldToLattice(x, y, z);
		l = (int)floor(x);
		m = (int)floor(y);
		s = x-l;
		t = y-m;
		S = round(LUTSIZE*s);
		T = round(LUTSIZE*t);
		norm = 0;
		for (j = 0; j < 4; j++) {
			B_J = this->LookupTable[T][j];
			for (i = 0; i < 4; i++) {
				B_I = B_J * this->LookupTable[S][i];
				norm += B_I * B_I;
			}
		}
		for (j = 0; j < 4; j++) {
			J = j + m - 1;
			if ((J >= 0) && (J < _y)) {
				B_J = this->LookupTable[T][j];
				for (i = 0; i < 4; i++) {
					I = i + l - 1;
					if ((I >= 0) && (I < _x)) {
						B_I = B_J * this->LookupTable[S][i];
						basis = B_I / norm;
						basis2 = B_I * B_I;
						phi = x2[index] * basis;
						data[0][J][I]._x += basis2 * phi;
						phi = y2[index] * basis;
						data[0][J][I]._y += basis2 * phi;
						ds[0][J][I] += basis2;
					}
				}
			}
		}
	}

	// Final loop: Calculate new control points
	for (j = -2; j < _y+2; j++) {
		for (i = -2; i < _x+2; i++) {
			if(ds[0][j][i]>0){
				_data[0][j][i]._x = data[0][j][i]._x / ds[0][j][i];
				_data[0][j][i]._y = data[0][j][i]._y / ds[0][j][i];
			}
		}
	}

	// Calculate residual error
	for (index = 0; index < no; index++) {
		x = x1[index];
		y = y1[index];
		z = z1[index];
		this->LocalDisplacement(x, y, z);
		x2[index] -= x;
		y2[index] -= y;
		z2[index] -= z;
	}

	// Deallocate memory
	Deallocate(data, _x, _y, _z);
	Deallocate(ds, _x, _y, _z);
}

void irtkBSplineFreeFormTransformation3D::ApproximateAsNew3D(const double *x1, const double *y1, const double *z1, double *x2, double *y2, double *z2, int no)
{
	int i, j, k, l, m, n, I, J, K, S, T, U, index;
	double s, t, u, x, y, z, B_I, B_J, B_K, basis, basis2, phi, norm;;

	// Allocate memory for control points
	irtkVector3D<double> ***data = NULL;
	data = Allocate(data, _x, _y, _z);

	// Allocate memory for temporary storage
	double ***ds = NULL;
	ds = Allocate(ds, _x, _y, _z);

	// Initialize data structures
	for (k = -2; k < _z+2; k++) {
		for (j = -2; j < _y+2; j++) {
			for (i = -2; i < _x+2; i++) {
				data[k][j][i]._x = 0;
				data[k][j][i]._y = 0;
				data[k][j][i]._z = 0;
				ds[k][j][i] = 0;
			}
		}
	}
	// Initial loop: Calculate change of control points
	for (index = 0; index < no; index++) {
		x = x1[index];
		y = y1[index];
		z = z1[index];
		this->WorldToLattice(x, y, z);
		l = (int)floor(x);
		m = (int)floor(y);
		n = (int)floor(z);
		s = x-l;
		t = y-m;
		u = z-n;
		S = round(LUTSIZE*s);
		T = round(LUTSIZE*t);
		U = round(LUTSIZE*u);
		norm = 0;
		for (k = 0; k < 4; k++) {
			B_K = this->LookupTable[U][k];
			for (j = 0; j < 4; j++) {
				B_J = B_K * this->LookupTable[T][j];
				for (i = 0; i < 4; i++) {
					B_I = B_J * this->LookupTable[S][i];
					norm += B_I * B_I;
				}
			}
		}
		for (k = 0; k < 4; k++) {
			K = k + n - 1;
			if ((K >= 0) && (K < _z)) {
				B_K = this->LookupTable[U][k];
				for (j = 0; j < 4; j++) {        
					J = j + m - 1;
					if ((J >= 0) && (J < _y)) {
						B_J = B_K * this->LookupTable[T][j];
						for (i = 0; i < 4; i++) {
							I = i + l - 1;
							if ((I >= 0) && (I < _x)) {
								B_I = B_J * this->LookupTable[S][i];
								basis = B_I / norm;
								basis2 = B_I * B_I;
								phi = x2[index] * basis;
								data[K][J][I]._x += basis2 * phi;
								phi = y2[index] * basis;
								data[K][J][I]._y += basis2 * phi;
								phi = z2[index] * basis;
								data[K][J][I]._z += basis2 * phi;
								ds[K][J][I] += basis2;
							}
						}
					}
				}
			}
		}
	}

	// Final loop: Calculate new control points
	for (k = -2; k < _z+2; k++) {
		for (j = -2; j < _y+2; j++) {
			for (i = -2; i < _x+2; i++) {
				if (ds[k][j][i] > 0) {
					_data[k][j][i]._x = data[k][j][i]._x / ds[k][j][i];
					_data[k][j][i]._y = data[k][j][i]._y / ds[k][j][i];
					_data[k][j][i]._z = data[k][j][i]._z / ds[k][j][i];
				}
			}
		}
	}

	// Calculate residual error
	for (index = 0; index < no; index++) {
		x = x1[index];
		y = y1[index];
		z = z1[index];
		this->LocalDisplacement(x, y, z);
		x2[index] -= x;
		y2[index] -= y;
		z2[index] -= z;
	}

	// Deallocate memory
	Deallocate(data, _x, _y, _z);
	Deallocate(ds, _x, _y, _z);
}

void irtkBSplineFreeFormTransformation3D::ApproximateAsNew(const double *x1, const double *y1, const double *z1, double *x2, double *y2, double *z2, int no)
{
	if (_z == 1) {
		ApproximateAsNew2D(x1, y1, z1, x2, y2, z2, no);
	} else {
		ApproximateAsNew3D(x1, y1, z1, x2, y2, z2, no);
	}
}

void irtkBSplineFreeFormTransformation3D::ApproximateGradient(const double *x1, const double *y1, const double *z1, double *x2, double *y2, double *z2, int no, double *gradient)
{
	if (_z == 1) {
		ApproximateGradient2D(x1, y1, z1, x2, y2, z2, no, gradient);
	} else {
		ApproximateGradient3D(x1, y1, z1, x2, y2, z2, no, gradient);
	}
}

void irtkBSplineFreeFormTransformation3D::ApproximateGradient2D(const double *x1, const double *y1, const double *z1, double *x2, double *y2, double *z2, int no, double *gradient)
{

	int i, j, l, m, I, J, S, T, index;
	double s, t, x, y, z, B_I, B_J, basis;

	// Allocate memory for control points
	irtkVector3D<double> ***data = NULL;
	data = Allocate(data, _x, _y, _z);

	// Initialize data structures
	for (j = 0; j < _y; j++) {
		for (i = 0; i < _x; i++) {
			data[0][j][i]._x = 0;
			data[0][j][i]._y = 0;
			data[0][j][i]._z = 0;
		}
	}
	// Initial loop: Calculate change of control points
	for (index = 0; index < no; index++) {
		// Check if it is active
		if(x2[index] != 0 || y2[index] != 0 || z2[index] != 0){
			x = x1[index];
			y = y1[index];
			z = z1[index];
			this->WorldToLattice(x, y, z);
			l = (int)floor(x);
			m = (int)floor(y);
			s = x-l;
			t = y-m;
			S = round(LUTSIZE*s);
			T = round(LUTSIZE*t);
			for (j = 0; j < 4; j++) {
				J = j + m - 1;
				if ((J >= 0) && (J < _y)) {
					B_J = this->LookupTable[T][j];
					for (i = 0; i < 4; i++) {
						I = i + l - 1;
						if ((I >= 0) && (I < _x)) {
							B_I = B_J * this->LookupTable[S][i];
							basis = B_I;
							data[0][J][I]._x += x2[index] * basis;
							data[0][J][I]._y += y2[index] * basis;
						}
					}
				}
			}
		}
	}

	// Final loop: Calculate new control points
	for (j = 0; j < _y; j++) {
		for (i = 0; i < _x; i++) {
			index = this->LatticeToIndex(i,j,0);
			if(this->irtkTransformation::GetStatus(index) == _Active 
				||this->irtkTransformation::GetStatus(index+_x*_y*_z) == _Active){
					gradient[index] = data[0][j][i]._x;
					gradient[index + _x*_y*_z] = data[0][j][i]._y;
					gradient[index + 2*_x*_y*_z] = 0;
			}else{
				gradient[index] = 0;
				gradient[index + _x*_y*_z] = 0;
				gradient[index + 2*_x*_y*_z] = 0;
			}
		}
	}

	// Deallocate memory
	Deallocate(data, _x, _y, _z);
}

void irtkBSplineFreeFormTransformation3D::ApproximateGradient3D(const double *x1, const double *y1, const double *z1, double *x2, double *y2, double *z2, int no, double *gradient)
{

	int i, j, k, l, m, n, I, J, K, S, T, U, index;
	double s, t, u, x, y, z, B_I, B_J, B_K, basis;

	// Allocate memory for control points
	irtkVector3D<double> ***data = NULL;
	data = Allocate(data, _x, _y, _z);

	// Initialize data structures
	for (k = 0; k < _z; k++) {
		for (j = 0; j < _y; j++) {
			for (i = 0; i < _x; i++) {
				data[k][j][i]._x = 0;
				data[k][j][i]._y = 0;
				data[k][j][i]._z = 0;
			}
		}
	}
	// Initial loop: Calculate change of control points
	for (index = 0; index < no; index++) {
		// Check if it is active
		if(x2[index] != 0 || y2[index] != 0 || z2[index] != 0){
			x = x1[index];
			y = y1[index];
			z = z1[index];
			this->WorldToLattice(x, y, z);
			l = (int)floor(x);
			m = (int)floor(y);
			n = (int)floor(z);
			s = x-l;
			t = y-m;
			u = z-n;
			S = round(LUTSIZE*s);
			T = round(LUTSIZE*t);
			U = round(LUTSIZE*u);
			for (k = 0; k < 4; k++) {
				K = k + n - 1;
				if ((K >= 0) && (K < _z)) {
					B_K = this->LookupTable[U][k];
					for (j = 0; j < 4; j++) {
						J = j + m - 1;
						if ((J >= 0) && (J < _y)) {
							B_J = B_K * this->LookupTable[T][j];
							for (i = 0; i < 4; i++) {             
								I = i + l - 1;
								if ((I >= 0) && (I < _x)) {
									B_I = B_J * this->LookupTable[S][i];
									basis = B_I;
									data[K][J][I]._x += x2[index] * basis;
									data[K][J][I]._y += y2[index] * basis;
									data[K][J][I]._z += z2[index] * basis;
								}
							}
						}
					}
				}
			}
		}
	}

	// Final loop: Calculate new control points
	for (k = 0; k < _z; k++) {
		for (j = 0; j < _y; j++) {
			for (i = 0; i < _x; i++) {
				index = this->LatticeToIndex(i,j,k);
				if(this->irtkTransformation::GetStatus(index) == _Active 
					||this->irtkTransformation::GetStatus(index+_x*_y*_z) == _Active
					||this->irtkTransformation::GetStatus(index+2*_x*_y*_z) == _Active){
						gradient[index] = data[k][j][i]._x;
						gradient[index + _x*_y*_z] = data[k][j][i]._y;
						gradient[index + 2*_x*_y*_z] = data[k][j][i]._z;
				}else{
					gradient[index] = 0;
					gradient[index + _x*_y*_z] = 0;
					gradient[index + 2*_x*_y*_z] = 0;
				}
			}
		}
	}

	// Deallocate memory
	Deallocate(data, _x, _y, _z);
}

double irtkBSplineFreeFormTransformation3D::Approximate2D(const double *x1, const double *y1, const double *z1, double *x2, double *y2, double *z2, int no)
{
	int i, j, l, m, I, J, S, T, index;
	double s, t, x, y, z, B_I, B_J, basis, basis2, error, phi, norm;

	// Allocate memory for control points
	irtkVector3D<double> ***data = NULL;
	data = Allocate(data, _x, _y, _z);

	// Allocate memory for temporary storage
	double ***ds = NULL;
	ds = Allocate(ds, _x, _y, _z);

	// Subtract displacements which are approximated by current control points
	for (index = 0; index < no; index++) {
		x = x1[index];
		y = y1[index];
		z = z1[index];
		this->LocalDisplacement(x, y, z);
		x2[index] -= x;
		y2[index] -= y;
		z2[index] -= z;
	}

	// Initialize data structures
	for (j = -2; j < _y+2; j++) {
		for (i = -2; i < _x+2; i++) {
			data[0][j][i]._x = 0;
			data[0][j][i]._y = 0;
			data[0][j][i]._z = 0;
			ds[0][j][i]      = 0;
		}
	}

	// Initial loop: Calculate change of control points
	for (index = 0; index < no; index++) {
		x = x1[index];
		y = y1[index];
		z = z1[index];
		this->WorldToLattice(x, y, z);
		l = (int)floor(x);
		m = (int)floor(y);
		s = x-l;
		t = y-m;
		S = round(LUTSIZE*s);
		T = round(LUTSIZE*t);
		norm = 0;
		for (j = 0; j < 4; j++) {
			B_J = this->LookupTable[T][j];
			for (i = 0; i < 4; i++) {
				B_I = B_J * this->LookupTable[S][i];
				norm += B_I * B_I;
			}
		}
		for (j = 0; j < 4; j++) {
			B_J = this->LookupTable[T][j];
			J = j + m - 1;
			if ((J >= -2) && (J < _y+2)) {
				for (i = 0; i < 4; i++) {
					B_I = B_J * this->LookupTable[S][i];
					I = i + l - 1;
					if ((I >= -2) && (I < _x+2)) {
						basis = B_I / norm;
						basis2 = B_I * B_I;
						phi = x2[index] * basis;
						data[0][J][I]._x += basis2 * phi;
						phi = y2[index] * basis;
						data[0][J][I]._y += basis2 * phi;
						ds[0][J][I] += basis2;
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
		this->LocalDisplacement(x, y, z);
		x2[index] += x;
		y2[index] += y;
		z2[index] += z;
	}

	// Final loop: Calculate new control points
	for (j = -2; j < _y+2; j++) {
		for (i = -2; i < _x+2; i++) {
			if (ds[0][j][i] > 0) {
				_data[0][j][i]._x += data[0][j][i]._x / ds[0][j][i];
				_data[0][j][i]._y += data[0][j][i]._y / ds[0][j][i];
			}
		}
	}

	// Calculate residual error
	error = 0;
	for (index = 0; index < no; index++) {
		x = x1[index];
		y = y1[index];
		z = z1[index];
		this->LocalDisplacement(x, y, z);
		x2[index] -= x;
		y2[index] -= y;
		z2[index] -= z;
		// Calculate error
		error += sqrt(x2[index]*x2[index]+y2[index]*y2[index]+z2[index]*z2[index]);
	}
	error = error / (double)no;

	// Deallocate memory
	Deallocate(data, _x, _y, _z);
	Deallocate(ds, _x, _y, _z);

	// Return error
	return error;
}

double irtkBSplineFreeFormTransformation3D::Approximate3D(const double *x1, const double *y1, const double *z1, double *x2, double *y2, double *z2, int no)
{
	int i, j, k, l, m, n, I, J, K, S, T, U, index;
	double s, t, u, x, y, z, B_I, B_J, B_K, basis, basis2, error, phi, norm;

	// Allocate memory for control points
	irtkVector3D<double> ***data = NULL;
	data = Allocate(data, _x, _y, _z);

	// Allocate memory for temporary storage
	double ***ds = NULL;
	ds = Allocate(ds, _x, _y, _z);

	// Subtract displacements which are approximated by current control points
	for (index = 0; index < no; index++) {
		x = x1[index];
		y = y1[index];
		z = z1[index];
		this->LocalDisplacement(x, y, z);
		x2[index] -= x;
		y2[index] -= y;
		z2[index] -= z;
	}

	// Initialize data structures
	for (k = -2; k < _z+2; k++) {
		for (j = -2; j < _y+2; j++) {
			for (i = -2; i < _x+2; i++) {
				data[k][j][i]._x = 0;
				data[k][j][i]._y = 0;
				data[k][j][i]._z = 0;
				ds[k][j][i] = 0;
			}
		}
	}
	// Initial loop: Calculate change of control points
	for (index = 0; index < no; index++) {
		x = x1[index];
		y = y1[index];
		z = z1[index];
		this->WorldToLattice(x, y, z);
		l = (int)floor(x);
		m = (int)floor(y);
		n = (int)floor(z);
		s = x-l;
		t = y-m;
		u = z-n;
		S = round(LUTSIZE*s);
		T = round(LUTSIZE*t);
		U = round(LUTSIZE*u);
		norm = 0;
		for (k = 0; k < 4; k++) {
			B_K = this->LookupTable[U][k];
			for (j = 0; j < 4; j++) {
				B_J = B_K * this->LookupTable[T][j];
				for (i = 0; i < 4; i++) {
					B_I = B_J * this->LookupTable[S][i];
					norm += B_I * B_I;
				}
			}
		}
		for (k = 0; k < 4; k++) {
			B_K = this->LookupTable[U][k];
			K = k + n - 1;
			if ((K >= -2) && (K < _z+2)) {
				for (j = 0; j < 4; j++) {
					B_J = B_K * this->LookupTable[T][j];
					J = j + m - 1;
					if ((J >= -2) && (J < _y+2)) {
						for (i = 0; i < 4; i++) {
							B_I = B_J * this->LookupTable[S][i];
							I = i + l - 1;
							if ((I >= -2) && (I < _x+2)) {
								basis = B_I / norm;
								basis2 = B_I * B_I;
								phi = x2[index] * basis;
								data[K][J][I]._x += basis2 * phi;
								phi = y2[index] * basis;
								data[K][J][I]._y += basis2 * phi;
								phi = z2[index] * basis;
								data[K][J][I]._z += basis2 * phi;
								ds[K][J][I] += basis2;
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
		this->LocalDisplacement(x, y, z);
		x2[index] += x;
		y2[index] += y;
		z2[index] += z;
	}

	// Final loop: Calculate new control points
	for (k = -2; k < _z+2; k++) {
		for (j = -2; j < _y+2; j++) {
			for (i = -2; i < _x+2; i++) {
				if (ds[k][j][i] > 0) {
					_data[k][j][i]._x += data[k][j][i]._x / ds[k][j][i];
					_data[k][j][i]._y += data[k][j][i]._y / ds[k][j][i];
					_data[k][j][i]._z += data[k][j][i]._z / ds[k][j][i];
				}
			}
		}
	}

	// Calculate residual error
	error = 0;
	for (index = 0; index < no; index++) {
		x = x1[index];
		y = y1[index];
		z = z1[index];
		this->LocalDisplacement(x, y, z);
		x2[index] -= x;
		y2[index] -= y;
		z2[index] -= z;
		// Calculate error
		error += sqrt(x2[index]*x2[index]+y2[index]*y2[index]+z2[index]*z2[index]);
	}
	error = error / (double)no;

	// Deallocate memory
	Deallocate(data, _x, _y, _z);
	Deallocate(ds, _x, _y, _z);

	// Return error
	return error;
}

double irtkBSplineFreeFormTransformation3D::Approximate(const double *x1, const double *y1, const double *z1, double *x2, double *y2, double *z2, int no)
{
	if (_z == 1) {
		return Approximate2D(x1, y1, z1, x2, y2, z2, no);
	} else {
		return Approximate3D(x1, y1, z1, x2, y2, z2, no);
	}
}

void irtkBSplineFreeFormTransformation3D::Interpolate(const double* dxs, const double* dys, const double* dzs)
{
	irtkGenericImage<double> xCoeffs, yCoeffs, zCoeffs;

	if (_z == 1) {
		ComputeCoefficients2D(dxs, dys, dzs, xCoeffs, yCoeffs, zCoeffs);
	} else {
		ComputeCoefficients3D(dxs, dys, dzs, xCoeffs, yCoeffs, zCoeffs);

	}
	for (int z = 0; z < _z; z++) {
		for (int y = 0; y < _y; y++) {
			for (int x = 0; x < _x; x++) {
				_data[z][y][x]._x = xCoeffs(x, y, z);
				_data[z][y][x]._y = yCoeffs(x, y, z);
				_data[z][y][x]._z = zCoeffs(x, y, z);
			}
		}
	}
}

void irtkBSplineFreeFormTransformation3D::Subdivide()
{
	if (_z == 1) {
		this->Subdivide2D();
	} else {
		this->Subdivide3D();
	}
}

void irtkBSplineFreeFormTransformation3D::Subdivide2D()
{
	int i, j, i1, j1, i2, j2;

	// Weights for subdivision
	double w[2][3];
	w[1][0] = 0;
	w[1][1] = 1.0/2.0;
	w[1][2] = 1.0/2.0;
	w[0][0] = 1.0/8.0;
	w[0][1] = 6.0/8.0;
	w[0][2] = 1.0/8.0;

	// Allocate memory for new control points
	irtkVector3D<double> ***data = NULL;
	data = this->Allocate(data, 2*_x-1, 2*_y-1, _z);

	for (i = 0; i < _x; i++) {
		for (j = 0; j < _y; j++) {
			for (i1 = 0; i1 < 2; i1++) {
				for (j1 = 0; j1 < 2; j1++) {
					data[0][2*j+j1][2*i+i1]._x = 0;
					data[0][2*j+j1][2*i+i1]._y = 0;
					data[0][2*j+j1][2*i+i1]._z = 0;
					for (i2 = 0; i2 < 3; i2++) {
						for (j2 = 0; j2 < 3; j2++) {
							data[0][2*j+j1][2*i+i1]._x += w[i1][i2] * w[j1][j2] * _data[0][j+j2-1][i+i2-1]._x;
							data[0][2*j+j1][2*i+i1]._y += w[i1][i2] * w[j1][j2] * _data[0][j+j2-1][i+i2-1]._y;
							data[0][2*j+j1][2*i+i1]._z += w[i1][i2] * w[j1][j2] * _data[0][j+j2-1][i+i2-1]._z;
						}
					}
				}
			}
		}
	}

	// Deallocate points
	this->Deallocate(_data, _x, _y, _z);
	delete []_status;

	// Update pointers to control points
	_data = data;

	// Increase number of control points
	_x = 2*_x - 1;
	_y = 2*_y - 1;
	_z = 1;

	// Recalculate control point spacing
	_dx /= 2.0;
	_dy /= 2.0;

	// Update transformation matrix
	this->UpdateMatrix();

	// Initialize memory for control point status
	_status = new _Status[3*_x*_y*_z];
	for (i = 0; i < 2*_x*_y*_z; i++) {
		_status[i] = _Active;
	}
	for (i = 2*_x*_y*_z; i < 3*_x*_y*_z; i++) {
		_status[i] = _Passive;
	}
}

void irtkBSplineFreeFormTransformation3D::Subdivide3D()
{
	int i, j, k, i1, j1, k1, i2, j2, k2;

	// Weights for subdivision
	double w[2][3];
	w[1][0] = 0;
	w[1][1] = 1.0/2.0;
	w[1][2] = 1.0/2.0;
	w[0][0] = 1.0/8.0;
	w[0][1] = 6.0/8.0;
	w[0][2] = 1.0/8.0;

	// Allocate memory for new control points
	irtkVector3D<double> ***data = NULL;
	data = this->Allocate(data, 2*_x-1, 2*_y-1, 2*_z-1);

	for (i = 0; i < _x; i++) {
		for (j = 0; j < _y; j++) {
			for (k = 0; k < _z; k++) {
				for (i1 = 0; i1 < 2; i1++) {
					for (j1 = 0; j1 < 2; j1++) {
						for (k1 = 0; k1 < 2; k1++) {
							data[2*k+k1][2*j+j1][2*i+i1]._x = 0;
							data[2*k+k1][2*j+j1][2*i+i1]._y = 0;
							data[2*k+k1][2*j+j1][2*i+i1]._z = 0;
							for (i2 = 0; i2 < 3; i2++) {
								for (j2 = 0; j2 < 3; j2++) {
									for (k2 = 0; k2 < 3; k2++) {
										data[2*k+k1][2*j+j1][2*i+i1]._x += w[i1][i2] * w[j1][j2] * w[k1][k2] * _data[k+k2-1][j+j2-1][i+i2-1]._x;
										data[2*k+k1][2*j+j1][2*i+i1]._y += w[i1][i2] * w[j1][j2] * w[k1][k2] * _data[k+k2-1][j+j2-1][i+i2-1]._y;
										data[2*k+k1][2*j+j1][2*i+i1]._z += w[i1][i2] * w[j1][j2] * w[k1][k2] * _data[k+k2-1][j+j2-1][i+i2-1]._z;
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
	this->Deallocate(_data, _x, _y, _z);
	delete []_status;

	// Update pointers to control points
	_data = data;

	// Increase number of control points
	_x = 2*_x - 1;
	_y = 2*_y - 1;
	_z = 2*_z - 1;

	// Recalculate control point spacing
	_dx /= 2.0;
	_dy /= 2.0;
	_dz /= 2.0;

	// Update transformation matrix
	this->UpdateMatrix();

	// Initialize memory for control point status
	_status = new _Status[3*_x*_y*_z];
	for (i = 0; i < 3*_x*_y*_z; i++) {
		_status[i] = _Active;
	}
}

double irtkBSplineFreeFormTransformation3D::InitialAntiCausalCoefficient(double c[], int DataLength, double z)
{
	/* this initialization corresponds to mirror boundaries */
	return((z / (z * z - 1.0)) * (z * c[DataLength - 2] + c[DataLength - 1]));
}

double irtkBSplineFreeFormTransformation3D::InitialCausalCoefficient(double c[], int DataLength, double z, double Tolerance)
{
	double Sum, zn, z2n, iz;
	int n, Horizon;

	/* this initialization corresponds to mirror boundaries */
	Horizon = DataLength;
	if (Tolerance > 0.0) {
		Horizon = (int)ceil(log(Tolerance) / log(fabs(z)));
	}
	if (Horizon < DataLength) {
		/* accelerated loop */
		zn = z;
		Sum = c[0];
		for (n = 1; n < Horizon; n++) {
			Sum += zn * c[n];
			zn *= z;
		}
		return(Sum);
	} else {
		/* full loop */
		zn = z;
		iz = 1.0 / z;
		z2n = pow(z, (double)(DataLength - 1));
		Sum = c[0] + z2n * c[DataLength - 1];
		z2n *= z2n * iz;
		for (n = 1; n <= DataLength - 2; n++) {
			Sum += (zn + z2n) * c[n];
			zn *= z;
			z2n *= iz;
		}
		return(Sum / (1.0 - zn * zn));
	}
}

void irtkBSplineFreeFormTransformation3D::ConvertToInterpolationCoefficients(double* c, int DataLength, double* z, int NbPoles, double Tolerance)
{
	double Lambda = 1.0;
	int n, k;

	/* special case required by mirror boundaries */
	if (DataLength == 1) {
		return;
	}

	/* compute the overall gain */
	for (k = 0; k < NbPoles; k++) {
		Lambda = Lambda * (1.0 - z[k]) * (1.0 - 1.0 / z[k]);
	}

	/* apply the gain */
	for (n = 0; n < DataLength; n++) {
		c[n] *= Lambda;
	}

	/* loop over all poles */
	for (k = 0; k < NbPoles; k++) {
		/* causal initialization */
		c[0] = InitialCausalCoefficient(c, DataLength, z[k], Tolerance);
		/* causal recursion */
		for (n = 1; n < DataLength; n++) {
			c[n] += z[k] * c[n - 1];
		}
		/* anticausal initialization */
		c[DataLength - 1] = InitialAntiCausalCoefficient(c, DataLength, z[k]);
		/* anticausal recursion */
		for (n = DataLength - 2; 0 <= n; n--) {
			c[n] = z[k] * (c[n + 1] - c[n]);
		}
	}
}

void irtkBSplineFreeFormTransformation3D::ComputeCoefficients2D(const double* dxs, const double* dys, const double*, irtkGenericImage<double>& xCoeffs, irtkGenericImage<double>& yCoeffs, irtkGenericImage<double>& zCoeffs)
{
	int x, y, NbPoles = 1;
	double Pole[2];

	Pole[0] = sqrt(3.0) - 2.0;

	// Initialize coefficient images.
	xCoeffs = irtkGenericImage<double>(_x, _y, _z);
	yCoeffs = irtkGenericImage<double>(_x, _y, _z);
	zCoeffs = irtkGenericImage<double>(_x, _y, _z);

	// Convert the displacements into interpolation coefficients for each
	// direction.

	// In-place separable process, along x.
	double* xdata = new double[_x];
	double* ydata = new double[_x];
	for (y = 0; y < _y; y++) {
		for (x = 0; x < _x; x++) {
			int index = x + y*_x;

			xdata[x] = dxs[index];
			ydata[x] = dys[index];
		}

		ConvertToInterpolationCoefficients(xdata, _x, Pole, NbPoles,
			DBL_EPSILON);
		ConvertToInterpolationCoefficients(ydata, _x, Pole, NbPoles,
			DBL_EPSILON);

		for (x = 0; x < _x; x++) {
			xCoeffs(x, y, 0) = xdata[x];
			yCoeffs(x, y, 0) = ydata[x];
		}
	}

	delete[] xdata;
	delete[] ydata;

	// In-place separable process, along y.
	xdata = new double[_y];
	ydata = new double[_y];
	for (x = 0; x < _x; x++) {
		for (y = 0; y < _y; y++) {
			xdata[y] = xCoeffs(x, y, 0);
			ydata[y] = yCoeffs(x, y, 0);
		}

		ConvertToInterpolationCoefficients(xdata, _y, Pole, NbPoles,
			DBL_EPSILON);
		ConvertToInterpolationCoefficients(ydata, _y, Pole, NbPoles,
			DBL_EPSILON);

		for (y = 0; y < _y; y++) {
			xCoeffs(x, y, 0) = xdata[y];
			yCoeffs(x, y, 0) = ydata[y];
		}
	}
	delete[] xdata;
	delete[] ydata;
}

void irtkBSplineFreeFormTransformation3D::ComputeCoefficients3D(const double* dxs, const double* dys, const double* dzs, irtkGenericImage<double>& xCoeffs, irtkGenericImage<double>& yCoeffs, irtkGenericImage<double>& zCoeffs)
{
	int x, y, z, NbPoles = 1;
	double Pole[2];

	Pole[0] = sqrt(3.0) - 2.0;

	// Initialize coefficient images.
	xCoeffs = irtkGenericImage<double>(_x, _y, _z);
	yCoeffs = irtkGenericImage<double>(_x, _y, _z);
	zCoeffs = irtkGenericImage<double>(_x, _y, _z);

	// Convert the displacements into interpolation coefficients for each
	// direction.

	// In-place separable process, along x.
	double* xdata = new double[_x];
	double* ydata = new double[_x];
	double* zdata = new double[_x];
	for (z = 0; z < _z; z++) {
		for (y = 0; y < _y; y++) {
			for (x = 0; x < _x; x++) {
				int index = x + y*_x + z*_x*_y;

				xdata[x] = dxs[index];
				ydata[x] = dys[index];
				zdata[x] = dzs[index];
			}

			ConvertToInterpolationCoefficients(xdata, _x, Pole, NbPoles,
				DBL_EPSILON);
			ConvertToInterpolationCoefficients(ydata, _x, Pole, NbPoles,
				DBL_EPSILON);
			ConvertToInterpolationCoefficients(zdata, _x, Pole, NbPoles,
				DBL_EPSILON);

			for (x = 0; x < _x; x++) {
				xCoeffs(x, y, z) = xdata[x];
				yCoeffs(x, y, z) = ydata[x];
				zCoeffs(x, y, z) = zdata[x];
			}
		}
	}
	delete[] xdata;
	delete[] ydata;
	delete[] zdata;

	// In-place separable process, along y.
	xdata = new double[_y];
	ydata = new double[_y];
	zdata = new double[_y];
	for (z = 0; z < _z; z++) {
		for (x = 0; x < _x; x++) {
			for (y = 0; y < _y; y++) {
				xdata[y] = xCoeffs(x, y, z);
				ydata[y] = yCoeffs(x, y, z);
				zdata[y] = zCoeffs(x, y, z);
			}

			ConvertToInterpolationCoefficients(xdata, _y, Pole, NbPoles,
				DBL_EPSILON);
			ConvertToInterpolationCoefficients(ydata, _y, Pole, NbPoles,
				DBL_EPSILON);
			ConvertToInterpolationCoefficients(zdata, _y, Pole, NbPoles,
				DBL_EPSILON);

			for (y = 0; y < _y; y++) {
				xCoeffs(x, y, z) = xdata[y];
				yCoeffs(x, y, z) = ydata[y];
				zCoeffs(x, y, z) = zdata[y];
			}
		}
	}
	delete[] xdata;
	delete[] ydata;
	delete[] zdata;

	// In-place separable process, along z.
	xdata = new double[_z];
	ydata = new double[_z];
	zdata = new double[_z];
	for (y = 0; y < _y; y++) {
		for (x = 0; x < _x; x++) {
			for (z = 0; z < _z; z++) {
				xdata[z] = xCoeffs(x, y, z);
				ydata[z] = yCoeffs(x, y, z);
				zdata[z] = zCoeffs(x, y, z);
			}

			ConvertToInterpolationCoefficients(xdata, _z, Pole, NbPoles,
				DBL_EPSILON);
			ConvertToInterpolationCoefficients(ydata, _z, Pole, NbPoles,
				DBL_EPSILON);
			ConvertToInterpolationCoefficients(zdata, _z, Pole, NbPoles,
				DBL_EPSILON);

			for (z = 0; z < _z; z++) {
				xCoeffs(x, y, z) = xdata[z];
				yCoeffs(x, y, z) = ydata[z];
				zCoeffs(x, y, z) = zdata[z];
			}
		}
	}
	delete[] xdata;
	delete[] ydata;
	delete[] zdata;
}

void irtkBSplineFreeFormTransformation3D::BoundingBoxCP(int index, irtkPoint &p1, irtkPoint &p2, double fraction) const
{
	int i, j, k;

	if (index >= _x*_y*_z) {
		index -= _x*_y*_z;
		if (index >= _x*_y*_z) {
			index -= _x*_y*_z;
		}
	}
	i = index/(_y*_z);
	j = index%(_y*_z)/_z;
	k = index%(_y*_z)%_z;
	p1 = irtkPoint(i-2*fraction, j-2*fraction, k-2*fraction);
	p2 = irtkPoint(i+2*fraction, j+2*fraction, k+2*fraction);
	this->LatticeToWorld(p1);
	this->LatticeToWorld(p2);
}

void irtkBSplineFreeFormTransformation3D::BoundingBoxImage(irtkGreyImage *image, int index, int &i1, int &j1, int &k1, int &i2, int &j2, int &k2, double fraction) const
{
	irtkPoint p1, p2;
	double x1, y1, z1, x2, y2, z2;

	// Calculate bounding box in world coordinates
	this->BoundingBoxCP(index, p1, p2, fraction);
	x1 = p1._x;
	y1 = p1._y;
	z1 = p1._z;
	x2 = p2._x;
	y2 = p2._y;
	z2 = p2._z;

	// Transform world coordinates to image coordinates
	image->WorldToImage(x1, y1, z1);
	image->WorldToImage(x2, y2, z2);

	// Calculate bounding box in image coordinates
	i1 = (x1 < 0) ? 0 : int(x1)+1;
	j1 = (y1 < 0) ? 0 : int(y1)+1;
	k1 = (z1 < 0) ? 0 : int(z1)+1;
	i2 = (int(x2) >= image->GetX()) ? image->GetX()-1 : int(x2);
	j2 = (int(y2) >= image->GetY()) ? image->GetY()-1 : int(y2);
	k2 = (int(z2) >= image->GetZ()) ? image->GetZ()-1 : int(z2);
}

void irtkBSplineFreeFormTransformation3D::BoundingBoxImage(irtkRealImage *image, int index, int &i1, int &j1, int &k1, int &i2, int &j2, int &k2, double fraction) const
{
	irtkPoint p1, p2;
	double x1, y1, z1, x2, y2, z2;

	// Calculate bounding box in world coordinates
	this->BoundingBoxCP(index, p1, p2, fraction);
	x1 = p1._x;
	y1 = p1._y;
	z1 = p1._z;
	x2 = p2._x;
	y2 = p2._y;
	z2 = p2._z;

	// Transform world coordinates to image coordinates
	image->WorldToImage(x1, y1, z1);
	image->WorldToImage(x2, y2, z2);

	// Calculate bounding box in image coordinates
	i1 = (x1 < 0) ? 0 : int(x1)+1;
	j1 = (y1 < 0) ? 0 : int(y1)+1;
	k1 = (z1 < 0) ? 0 : int(z1)+1;
	i2 = (int(x2) >= image->GetX()) ? image->GetX()-1 : int(x2);
	j2 = (int(y2) >= image->GetY()) ? image->GetY()-1 : int(y2);
	k2 = (int(z2) >= image->GetZ()) ? image->GetZ()-1 : int(z2);
}

void irtkBSplineFreeFormTransformation3D::MultiBoundingBoxImage(irtkGreyImage *image, int index, int &i1, int &j1, int &k1, int &i2, int &j2, int &k2, double fraction) const
{
	irtkPoint p1, p2;
	double x1, y1, z1, x2, y2, z2;

	// Calculate bounding box in world coordinates
	this->BoundingBoxCP(index, p1, p2, fraction*2);
	x1 = p1._x;
	y1 = p1._y;
	z1 = p1._z;
	x2 = p2._x;
	y2 = p2._y;
	z2 = p2._z;

	// Transform world coordinates to image coordinates
	image->WorldToImage(x1, y1, z1);
	image->WorldToImage(x2, y2, z2);

	if(x2<x1) {
		i1 = x2; x2 = x1; x1 = i1;
	}
	if(y2<y1) {
		j1 = y2; y2 = y1; y1 = j1;
	}
	if(z2<z1) {
		k1 = z2; z2 = z1; z1 = k1;
	}

	// Calculate bounding box in image coordinates
	i1 = (x1 < 0) ? 0 : int(x1)+1;
	j1 = (y1 < 0) ? 0 : int(y1)+1;
	k1 = (z1 < 0) ? 0 : int(z1)+1;
	i2 = (int(x2) >= image->GetX()) ? image->GetX()-1 : int(x2);
	j2 = (int(y2) >= image->GetY()) ? image->GetY()-1 : int(y2);
	k2 = (int(z2) >= image->GetZ()) ? image->GetZ()-1 : int(z2);
}

void irtkBSplineFreeFormTransformation3D::MultiBoundingBoxImage(irtkRealImage *image, int index, int &i1, int &j1, int &k1, int &i2, int &j2, int &k2, double fraction) const
{
	irtkPoint p1, p2;
	double x1, y1, z1, x2, y2, z2;

	// Calculate bounding box in world coordinates
	this->BoundingBoxCP(index, p1, p2, fraction*2);
	x1 = p1._x;
	y1 = p1._y;
	z1 = p1._z;
	x2 = p2._x;
	y2 = p2._y;
	z2 = p2._z;

	// Transform world coordinates to image coordinates
	image->WorldToImage(x1, y1, z1);
	image->WorldToImage(x2, y2, z2);

	if(x2<x1) {
		i1 = x2; x2 = x1; x1 = i1;
	}
	if(y2<y1) {
		j1 = y2; y2 = y1; y1 = j1;
	}
	if(z2<z1) {
		k1 = z2; z2 = z1; z1 = k1;
	}

	// Calculate bounding box in image coordinates
	i1 = (x1 < 0) ? 0 : int(x1)+1;
	j1 = (y1 < 0) ? 0 : int(y1)+1;
	k1 = (z1 < 0) ? 0 : int(z1)+1;
	i2 = (int(x2) >= image->GetX()) ? image->GetX()-1 : int(x2);
	j2 = (int(y2) >= image->GetY()) ? image->GetY()-1 : int(y2);
	k2 = (int(z2) >= image->GetZ()) ? image->GetZ()-1 : int(z2);
}

istream& irtkBSplineFreeFormTransformation3D::Import(istream& is)
{
	char buffer[255];
	int i, j, k, n, index, no_attr;
	double x1, y1, z1, x2, y2, z2;

	// Free memory if necessary
	_data = Deallocate(_data, _x, _y, _z);
	delete []_status;
	_x = 0;
	_y = 0;
	_z = 0;

	// Read keyword
	is >> buffer;

	if (strcmp(buffer, "AFFD_BINARY:") == 0) {

		// Read no. of control points
		is >> buffer;
		is >> buffer;
		is >> buffer >> _x  >> buffer >> _y >> buffer >> _z;

		// Read bounding box for control points
		is >> buffer >> x1 >> y1 >> z1 >> buffer >> x2 >> y2 >> z2;

		// Calculate origin
		_origin._x = (x2 + x1) / 2.0;
		_origin._y = (y2 + y1) / 2.0;
		_origin._z = (z2 + z1) / 2.0;

		// Calculate control point spacing
		_dx = (x2 - x1) / (_x-1) ;
		_dy = (y2 - y1) / (_y-1) ;
		if (z2 > z1) {
			_dz = (z2 - z1) / (_z-1);
		} else {
			_dz = 1;
		}

		// Skip rest of line
		is.get(buffer, 255);
		is.clear();
		is.seekg(1, ios::cur);

	} else {
		if (strcmp(buffer, "FFD_BINARY:") == 0) {

			// Read no. of control points
			is >> buffer;
			is >> buffer;
			is >> buffer >> _x  >> buffer >> _y >> buffer >> _z;

			// Skip rest of line
			is.get(buffer, 255);
			is.clear();
			is.seekg(1, ios::cur);

			// Read orientation, origin and spacing
			double tmp;

			is.read((char *)&tmp, sizeof(double));
#ifndef WORDS_BIGENDIAN
			swap64((char *)&tmp, (char *)&tmp, 1);
#endif
			_xaxis[0] = tmp;
			is.read((char *)&tmp, sizeof(double));
#ifndef WORDS_BIGENDIAN
			swap64((char *)&tmp, (char *)&tmp, 1);
#endif
			_xaxis[1] = tmp;
			is.read((char *)&tmp, sizeof(double));
#ifndef WORDS_BIGENDIAN
			swap64((char *)&tmp, (char *)&tmp, 1);
#endif
			_xaxis[2] = tmp;
			is.read((char *)&tmp, sizeof(double));
#ifndef WORDS_BIGENDIAN
			swap64((char *)&tmp, (char *)&tmp, 1);
#endif
			_yaxis[0] = tmp;
			is.read((char *)&tmp, sizeof(double));
#ifndef WORDS_BIGENDIAN
			swap64((char *)&tmp, (char *)&tmp, 1);
#endif
			_yaxis[1] = tmp;
			is.read((char *)&tmp, sizeof(double));
#ifndef WORDS_BIGENDIAN
			swap64((char *)&tmp, (char *)&tmp, 1);
#endif
			_yaxis[2] = tmp;
			is.read((char *)&tmp, sizeof(double));
#ifndef WORDS_BIGENDIAN
			swap64((char *)&tmp, (char *)&tmp, 1);
#endif
			_dx = tmp;
			is.read((char *)&tmp, sizeof(double));
#ifndef WORDS_BIGENDIAN
			swap64((char *)&tmp, (char *)&tmp, 1);
#endif
			_dy = tmp;
			is.read((char *)&tmp, sizeof(double));
#ifndef WORDS_BIGENDIAN
			swap64((char *)&tmp, (char *)&tmp, 1);
#endif
			_dz = tmp;
			is.read((char *)&tmp, sizeof(double));
#ifndef WORDS_BIGENDIAN
			swap64((char *)&tmp, (char *)&tmp, 1);
#endif
			_origin._x = tmp;
			is.read((char *)&tmp, sizeof(double));
#ifndef WORDS_BIGENDIAN
			swap64((char *)&tmp, (char *)&tmp, 1);
#endif
			_origin._y = tmp;
			is.read((char *)&tmp, sizeof(double));
#ifndef WORDS_BIGENDIAN
			swap64((char *)&tmp, (char *)&tmp, 1);
#endif
			_origin._z = tmp;

		} else {
			cerr << "istream& operator>> (istream& is, irtkBSplineFreeFormTransformation3D &transformation): Invalid file format: " << buffer << endl;
			exit(1);
		}
	}

	// Compute z-axis.
	_zaxis[0] = _xaxis[1]*_yaxis[2] - _xaxis[2]*_yaxis[1];
	_zaxis[1] = _xaxis[2]*_yaxis[0] - _xaxis[0]*_yaxis[2];
	_zaxis[2] = _xaxis[0]*_yaxis[1] - _xaxis[1]*_yaxis[0];

	// Calculate number of control points
	n = _x*_y*_z;

	// Initialize control points
	_data = Allocate(_data, _x, _y, _z);

	// Initialize control points
	for (i = -2; i < _x+2; i++) {
		for (j = -2; j < _y+2; j++) {
			for (k = -2; k < _z+2; k++) {
				_data[k][j][i]._x = 0;
				_data[k][j][i]._y = 0;
				_data[k][j][i]._z = 0;
			}
		}
	}

	// Initialize memory for control point status
	_status = new _Status[3*n];

	// Initialize control point status to active
	for (i = 0; i < 3*n; i++) {
		_status[i] = _Active;
	}

	// Allocate temporary memory
	float *data = new float[3*n*sizeof(float)];

	// Read binary data for control point values
	is.read((char *)data, 3*n*sizeof(float));

#ifndef WORDS_BIGENDIAN
	swap32((char *)data, (char *)data, 3*n);
#endif

	// Convert data
	index = 0;
	for (i = 0; i < _x; i++) {
		for (j = 0; j < _y; j++) {
			for (k = 0; k < _z; k++) {
				_data[k][j][i]._x = data[index];
				_data[k][j][i]._y = data[index+1];
				_data[k][j][i]._z = data[index+2];
				index += 3;
			}
		}
	}

	// Free temporary memory
	delete []data;

	// Skip rest of line
	is.get(buffer, 255);
	is.clear();
	is.seekg(1, ios::cur);

	// Read number of attributes
	is >> buffer;

	if (strcmp(buffer, "ATTRIBUTES:") == 0) {
		is >> no_attr;
	} else {
		cerr << "istream& operator>> irtkBSplineFreeFormTransformation3D: "
			<< "Expecting attributes " << buffer << endl;
		exit(1);
	}

	for (i = 0; i < no_attr; i++) {

		// Read attribute
		is >> buffer;

		// Allocate temporary memory
		unsigned char *status = new unsigned char[n*sizeof(unsigned char)];

		if (strcmp(buffer, "STATUS_X:") == 0) {

			// Skip rest of line
			is.get(buffer, 255);
			is.clear();
			is.seekg(1, ios::cur);
			// Read binary data for control point status
			is.read((char *)&(status[0]), n*sizeof(unsigned char));

			// Convert status
			for (j = 0; j < n; j++) {
				_status[j] = (_Status) status[j];
			}
		} else {
			if (strcmp(buffer, "STATUS_Y:") == 0) {

				// Skip rest of line
				is.get(buffer, 255);
				is.clear();
				is.seekg(1, ios::cur);
				// Read binary data for control point status
				is.read((char *)&(status[0]), n*sizeof(unsigned char));

				// Convert status
				for (j = 0; j < n; j++) {
					_status[j+n] = (_Status) status[j];
				}
			} else {
				if (strcmp(buffer, "STATUS_Z:") == 0) {

					// Skip rest of line
					is.get(buffer, 255);
					is.clear();
					is.seekg(1, ios::cur);
					// Read binary data for control point status
					is.read((char *)&(status[0]), n*sizeof(unsigned char));

					// Convert status
					for (j = 0; j < n; j++) {
						_status[j+2*n] = (_Status) status[j];
					}
				} else {
					cerr << "istream& operator>> irtkBSplineFreeFormTransformation3D ";
					cerr << "Unknown attribute:" << buffer << " terminating" << endl;
					exit(1);
				}
			}
		}

		// Free temporary memory
		delete [] status;

		// Skip rest of line
		is.get(buffer, 255);
		is.clear();
		is.seekg(1, ios::cur);

	}
	// Update transformation matrix
	this->UpdateMatrix();

	return is;
}

ostream& irtkBSplineFreeFormTransformation3D::Export(ostream& os)
{
	int i, j, k, n, index;
	double x1, y1, z1, x2, y2, z2;
	double zaxis[3];

	cerr << "irtkBSplineFreeFormTransformation3D::Export: DO NOT USE. THIS METHOD IS OBSOLETE" << endl;

	// Compute z-axis.
	zaxis[0] = _xaxis[1]*_yaxis[2] - _xaxis[2]*_yaxis[1];
	zaxis[1] = _xaxis[2]*_yaxis[0] - _xaxis[0]*_yaxis[2];
	zaxis[2] = _xaxis[0]*_yaxis[1] - _xaxis[1]*_yaxis[0];

	if ((zaxis [0] != _zaxis[0]) || (zaxis [1] != _zaxis[1]) || (zaxis [2] != _zaxis[2])) {
		cerr << "irtkBSplineFreeFormTransformation3D::Export: z-axis has not default value. Writing transformation will omit z-axis." << endl;
	}

	// Check if we need to save orientation
	if ((_xaxis[0] != 1) || (_xaxis[1] != 0) ||
		(_xaxis[2] != 0) || (_yaxis[0] != 0) ||
		(_yaxis[1] != 1) || (_yaxis[2] != 0)) {

			// Write keyword and no. of DOFs
			os << "FFD_BINARY: " << NumberOfDOFs() << endl;

			// Write no. of control points
			os << "Control points: " << _x << " x " <<
				_y << " x " << _z << endl;

			// Write orientation, origin and spacing
			double tmp;

			tmp = _xaxis[0];
#ifndef WORDS_BIGENDIAN
			swap64((char *)&tmp, (char *)&tmp, 1);
#endif
			os.write((char *)&tmp, sizeof(double));
			tmp = _xaxis[1];
#ifndef WORDS_BIGENDIAN
			swap64((char *)&tmp, (char *)&tmp, 1);
#endif
			os.write((char *)&tmp, sizeof(double));
			tmp = _xaxis[2];
#ifndef WORDS_BIGENDIAN
			swap64((char *)&tmp, (char *)&tmp, 1);
#endif
			os.write((char *)&tmp, sizeof(double));
			tmp = _yaxis[0];
#ifndef WORDS_BIGENDIAN
			swap64((char *)&tmp, (char *)&tmp, 1);
#endif
			os.write((char *)&tmp, sizeof(double));
			tmp = _yaxis[1];
#ifndef WORDS_BIGENDIAN
			swap64((char *)&tmp, (char *)&tmp, 1);
#endif
			os.write((char *)&tmp, sizeof(double));
			tmp = _yaxis[2];
#ifndef WORDS_BIGENDIAN
			swap64((char *)&tmp, (char *)&tmp, 1);
#endif
			os.write((char *)&tmp, sizeof(double));
			tmp = _dx;
#ifndef WORDS_BIGENDIAN
			swap64((char *)&tmp, (char *)&tmp, 1);
#endif
			os.write((char *)&tmp, sizeof(double));
			tmp = _dy;
#ifndef WORDS_BIGENDIAN
			swap64((char *)&tmp, (char *)&tmp, 1);
#endif
			os.write((char *)&tmp, sizeof(double));
			tmp = _dz;
#ifndef WORDS_BIGENDIAN
			swap64((char *)&tmp, (char *)&tmp, 1);
#endif
			os.write((char *)&tmp, sizeof(double));
			tmp = _origin._x;
#ifndef WORDS_BIGENDIAN
			swap64((char *)&tmp, (char *)&tmp, 1);
#endif
			os.write((char *)&tmp, sizeof(double));
			tmp = _origin._y;
#ifndef WORDS_BIGENDIAN
			swap64((char *)&tmp, (char *)&tmp, 1);
#endif
			os.write((char *)&tmp, sizeof(double));
			tmp = _origin._z;
#ifndef WORDS_BIGENDIAN
			swap64((char *)&tmp, (char *)&tmp, 1);
#endif
			os.write((char *)&tmp, sizeof(double));

	} else {

		// Write keyword and no. of DOFs
		os << "AFFD_BINARY: " << this->NumberOfDOFs() << endl;

		// Write no. of control points
		os << "Control points: " << _x << " x " <<
			_y << " x " << _z << endl;

		// Calculate bounding box
		x1 = 0;
		y1 = 0;
		z1 = 0;
		this->LatticeToWorld(x1, y1, z1);
		x2 = _x-1;
		y2 = _y-1;
		z2 = _z-1;
		this->LatticeToWorld(x2, y2, z2);

		// Write bounding box
		os << "Point_1: " << x1 << " " << y1 << " " << z1 << " ";
		os << "Point_2: " << x2 << " " << y2 << " " << z2 << endl;
	}

	// Compute no. of control points
	n = _x*_y*_z;

	// Allocate temporary memory
	float *data = new float[3*n*sizeof(float)];

	// Convert data
	index = 0;
	for (i = 0; i < _x; i++) {
		for (j = 0; j < _y; j++) {
			for (k = 0; k < _z; k++) {
				data[index]   = _data[k][j][i]._x;
				data[index+1] = _data[k][j][i]._y;
				data[index+2] = _data[k][j][i]._z;
				index += 3;
			}
		}
	}

#ifndef WORDS_BIGENDIAN
	swap32((char *)data, (char *)data, 3*n);
#endif

	// Write binary data for control point values
	os.write((char *)data, 3*n*sizeof(float));

	// Free temporary memory
	delete []data;

	// Write end of line
	os << endl;

	// Cast status (later label) to unsigned char to avoid byte order problems
	unsigned char *status = new unsigned char [3*n*sizeof(unsigned char)];
	for (index = 0; index < 3*n; index++) {
		status[index] = _status[index];
	}

	// Write attributes
	os << "ATTRIBUTES: 3" << endl;

	// Write data for control point status
	os << "STATUS_X:" << endl;
	os.write((char *)&(status[0]), n*sizeof(unsigned char));
	os << endl;
	os << "STATUS_Y:" << endl;
	os.write((char *)&(status[n]), n*sizeof(unsigned char));
	os << endl;
	os << "STATUS_Z:" << endl;
	os.write((char *)&(status[2*n]), n*sizeof(unsigned char));

	// Free temporary memory
	delete []status;

	// Write end of line
	os << endl;

	return os;
}

void irtkBSplineFreeFormTransformation3D::Print()
{
	// Print keyword and no. of DOFs
	cout << "BSplineFFD: " << this->NumberOfDOFs() << endl;

	// Print no. of control points
	cout << "Control points: " << _x << " x " << _y << " x " << _z << endl;
	cout << "Spacing: " << _dx << " x " << _dy << " x " << _dz << endl;
	cout << "Origin: " << _origin._x << " " << _origin._y << " " << _origin._z << endl;
	// Print x-axis
	cout << "X-axis is " << _xaxis[0] << " " << _xaxis[1] << " " << _xaxis[2] << endl;
	// Print x-axis
	cout << "Y-axis is " << _yaxis[0] << " " << _yaxis[1] << " " << _yaxis[2] << endl;
	// Print x-axis
	cout << "Z-axis is " << _zaxis[0] << " " << _zaxis[1] << " " << _zaxis[2] << endl;
}

irtkCifstream& irtkBSplineFreeFormTransformation3D::Read(irtkCifstream& from)
{
    unsigned int magic_no, trans_type;

    // Read magic no. for transformations
    from.ReadAsUInt(&magic_no, 1);
    if (magic_no != IRTKTRANSFORMATION_MAGIC) {
        cerr << "irtkBSplineFreeFormTransformation3D::Read: Not a vaild transformation file" << endl;
        exit(1);
    }

    // Read transformation type
    from.ReadAsUInt(&trans_type, 1);
    if ((trans_type != IRTKTRANSFORMATION_BSPLINE_FFD) && (trans_type != IRTKTRANSFORMATION_BSPLINE_FFD_EXT1)) {
        cerr << "irtkBSplineFreeFormTransformation3D::Read: Not a vaild B-Spline FFD transformation" << endl;
        exit(1);
    }

    // Read control point data
    this->ReadCP(from, trans_type == IRTKTRANSFORMATION_BSPLINE_FFD_EXT1);

    // Update transformation matrix
    this->UpdateMatrix();

    return from;
}

irtkCifstream& irtkBSplineFreeFormTransformation3D::ReadCP(irtkCifstream& from, bool ext1)
{
    double *data;
    int i, j, k, index;

    // Free memory if necessary
    _data = Deallocate(_data, _x, _y, _z);
    delete []_status;

    // Read no of control points
    from.ReadAsInt(&_x, 1);
    from.ReadAsInt(&_y, 1);
    from.ReadAsInt(&_z, 1);

    // Read orientation of bounding box
    from.ReadAsDouble(_xaxis, 3);
    from.ReadAsDouble(_yaxis, 3);
    if (ext1) {
        from.ReadAsDouble(_zaxis, 3);
    } else {
        _zaxis[0] = _xaxis[1]*_yaxis[2] - _xaxis[2]*_yaxis[1];
        _zaxis[1] = _xaxis[2]*_yaxis[0] - _xaxis[0]*_yaxis[2];
        _zaxis[2] = _xaxis[0]*_yaxis[1] - _xaxis[1]*_yaxis[0];
    }

    // Read spacing of bounding box
    from.ReadAsDouble(&_dx, 1);
    from.ReadAsDouble(&_dy, 1);
    from.ReadAsDouble(&_dz, 1);

    // Read spacing of bounding box
    from.ReadAsDouble(&_origin._x, 1);
    from.ReadAsDouble(&_origin._y, 1);
    from.ReadAsDouble(&_origin._z, 1);

    // Initialize control points
    _data = Allocate(_data, _x, _y, _z);

    // Initialize control points
    for (i = -2; i < _x+2; i++) {
        for (j = -2; j < _y+2; j++) {
            for (k = -2; k < _z+2; k++) {
                _data[k][j][i]._x = 0;
                _data[k][j][i]._y = 0;
                _data[k][j][i]._z = 0;
            }
        }
    }

    // Allocate temporary memory
    data = new double[3*_x*_y*_z];

    // Read control point data
    from.ReadAsDouble(data, 3*_x*_y*_z);

    // Convert data
    index = 0;
    for (i = 0; i < _x; i++) {
        for (j = 0; j < _y; j++) {
            for (k = 0; k < _z; k++) {
                _data[k][j][i]._x = data[index];
                _data[k][j][i]._y = data[index+1];
                _data[k][j][i]._z = data[index+2];
                index += 3;
            }
        }
    }

    // Free temporary memory
    delete []data;

    // Initialize memory for control point status
    _status = new _Status[3*_x*_y*_z];

    // Read control point status
    from.ReadAsInt((int *)_status, 3*_x*_y*_z);

    return from;
}

irtkCofstream& irtkBSplineFreeFormTransformation3D::Write(irtkCofstream& to)
{
    unsigned int magic_no, trans_type;

    // Write magic no. for transformations
    magic_no = IRTKTRANSFORMATION_MAGIC;
    to.WriteAsUInt(&magic_no, 1);

    // Write transformation type
    trans_type = IRTKTRANSFORMATION_BSPLINE_FFD_EXT1;
    to.WriteAsUInt(&trans_type, 1);

    // Write control point data
    return this->WriteCP(to);
}

irtkCofstream& irtkBSplineFreeFormTransformation3D::WriteCP(irtkCofstream& to)
{
    double *data;
    int i, j, k, index;

    // Write no of control points
    to.WriteAsInt(&_x, 1);
    to.WriteAsInt(&_y, 1);
    to.WriteAsInt(&_z, 1);

    // Write orientation of bounding box
    to.WriteAsDouble(_xaxis, 3);
    to.WriteAsDouble(_yaxis, 3);
    to.WriteAsDouble(_zaxis, 3);

    // Write spacing of bounding box
    to.WriteAsDouble(&_dx, 1);
    to.WriteAsDouble(&_dy, 1);
    to.WriteAsDouble(&_dz, 1);

    // Write spacing of bounding box
    to.WriteAsDouble(&_origin._x, 1);
    to.WriteAsDouble(&_origin._y, 1);
    to.WriteAsDouble(&_origin._z, 1);

    // Allocate temporary memory
    data = new double[3*_x*_y*_z];

    // Convert data
    index = 0;
    for (i = 0; i < _x; i++) {
        for (j = 0; j < _y; j++) {
            for (k = 0; k < _z; k++) {
                data[index]   = _data[k][j][i]._x;
                data[index+1] = _data[k][j][i]._y;
                data[index+2] = _data[k][j][i]._z;
                index += 3;
            }
        }
    }

    // Write control point data
    to.WriteAsDouble(data, 3*_x*_y*_z);

    // Free temporary memory
    delete []data;

    // Write control point status
    to.WriteAsInt((int *)_status, 3*_x*_y*_z);

    return to;
}
