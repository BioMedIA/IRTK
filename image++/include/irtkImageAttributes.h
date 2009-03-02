#ifndef IRTKattrATTRIBUTES_H

#define IRTKattrATTRIBUTES_H

/**
 * Class which defines the attributes of the imaging geometry
 */

class irtkImageAttributes : public irtkObject
{

public:

  /// Image x-dimension (in voxels)
  int _x;
  
  /// Image y-dimension (in voxels)
  int _y;
  
  /// Image z-dimension (in voxels)
  int _z;
  
  /// Image t-dimension (in voxels)
  int _t;

  /// Voxel x-dimensions (in mm)
  double _dx;
  
  /// Voxel y-dimensions (in mm)
  double _dy;
  
  /// Voxel z-dimensions (in mm)
  double _dz;
  
  /// Voxel t-dimensions (in ms)
  double _dt;

  /// Image x-origin
  double _xorigin;

  /// Image y-origin
  double _yorigin;
  
  /// Image z-origin
  double _zorigin;

  /// Image origin (temporal)
  double _torigin;

  /// Direction of x-axis
  double _xaxis[3];

  /// Direction of y-axis
  double _yaxis[3];

  /// Direction of z-axis
  double _zaxis[3];

  /// Constructor
  irtkImageAttributes();

  /// Copy constructor
  irtkImageAttributes(const irtkImageAttributes &);

  /// Copy operator
  irtkImageAttributes& operator= (const irtkImageAttributes &);

  /// Comparison operator
  Bool operator==(const irtkImageAttributes &attr) const;

};

inline irtkImageAttributes::irtkImageAttributes()
{
  _x  = 0;
  _y  = 0;
  _z  = 1;
  _t  = 1;

  // Default voxel size
  _dx = 1;
  _dy = 1;
  _dz = 1;
  _dt = 1;

  // Default origin
  _xorigin = 0;
  _yorigin = 0;
  _zorigin = 0;
  _torigin = 0;

  // Default x-axis
  _xaxis[0] = 1;
  _xaxis[1] = 0;
  _xaxis[2] = 0;

  // Default y-axis
  _yaxis[0] = 0;
  _yaxis[1] = 1;
  _yaxis[2] = 0;

  // Default z-axis
  _zaxis[0] = 0;
  _zaxis[1] = 0;
  _zaxis[2] = 1;
}

inline irtkImageAttributes::irtkImageAttributes(const irtkImageAttributes &attr)
{
  _x  = attr._x;
  _y  = attr._y;
  _z  = attr._z;
  _t  = attr._t;

  // Default voxel size
  _dx = attr._dx;
  _dy = attr._dy;
  _dz = attr._dz;
  _dt = attr._dt;

  // Default origin
  _xorigin = attr._xorigin;
  _yorigin = attr._yorigin;
  _zorigin = attr._zorigin;
  _torigin = attr._torigin;

  // Default x-axis
  _xaxis[0] = attr._xaxis[0];
  _xaxis[1] = attr._xaxis[1];
  _xaxis[2] = attr._xaxis[2];

  // Default y-axis
  _yaxis[0] = attr._yaxis[0];
  _yaxis[1] = attr._yaxis[1];
  _yaxis[2] = attr._yaxis[2];

  // Default z-axis
  _zaxis[0] = attr._zaxis[0];
  _zaxis[1] = attr._zaxis[1];
  _zaxis[2] = attr._zaxis[2];
}

inline irtkImageAttributes& irtkImageAttributes::operator=(const irtkImageAttributes &attr)
{
	
  _x  = attr._x;
  _y  = attr._y;
  _z  = attr._z;
  _t  = attr._t;

  // Default voxel size
  _dx = attr._dx;
  _dy = attr._dy;
  _dz = attr._dz;
  _dt = attr._dt;

  // Default origin
  _xorigin = attr._xorigin;
  _yorigin = attr._yorigin;
  _zorigin = attr._zorigin;
  _torigin = attr._torigin;

  // Default x-axis
  _xaxis[0] = attr._xaxis[0];
  _xaxis[1] = attr._xaxis[1];
  _xaxis[2] = attr._xaxis[2];

  // Default y-axis
  _yaxis[0] = attr._yaxis[0];
  _yaxis[1] = attr._yaxis[1];
  _yaxis[2] = attr._yaxis[2];

  // Default z-axis
  _zaxis[0] = attr._zaxis[0];
  _zaxis[1] = attr._zaxis[1];
  _zaxis[2] = attr._zaxis[2];
  
  return *this;
}

inline Bool irtkImageAttributes::operator==(const irtkImageAttributes &attr) const
{
  return ((_x  == attr._x)  && (_y  == attr._y)  && (_z  == attr._z) && (_t  == attr._t) &&
          (_dx == attr._dx) && (_dy == attr._dy) && (_dz == attr._dz) && (_dt == attr._dt) &&
          (_xaxis[0] == attr._xaxis[0]) && (_xaxis[1] == attr._xaxis[1]) && (_xaxis[2] == attr._xaxis[2]) &&
          (_yaxis[0] == attr._yaxis[0]) && (_yaxis[1] == attr._yaxis[1]) && (_yaxis[2] == attr._yaxis[2]) &&
          (_xaxis[0] == attr._zaxis[0]) && (_zaxis[1] == attr._zaxis[1]) && (_zaxis[2] == attr._zaxis[2]) &&
          (_xorigin == attr._xorigin) && (_yorigin == attr._yorigin) && (_zorigin == attr._zorigin) && 
          (_torigin == attr._torigin));
}

#endif
