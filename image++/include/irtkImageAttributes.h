#ifndef IRTKIMAGEATTRIBUTES_H

#define IRTKIMAGEATTRIBUTES_H

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

  /// Image x-origin (in mm)
  double _xorigin;

  /// Image y-origin (in mm)
  double _yorigin;
  
  /// Image z-origin (in mm)
  double _zorigin;

  /// Image t-origin (in ms)
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
  bool operator==(const irtkImageAttributes &attr) const;

  /// Get Index from Lattice
  int LatticeToIndex(int i, int j, int k, int l = 0) const;

  /// Get Index from Lattice
  void IndexToLattice(int index, int *i, int *j, int *k, int *l = NULL);

  /// Print attributes
  void Print();

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

inline irtkImageAttributes::irtkImageAttributes(const irtkImageAttributes &attr) : irtkObject(attr)
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

inline bool irtkImageAttributes::operator==(const irtkImageAttributes &attr) const
{
  return ((_x  == attr._x)  && (_y  == attr._y)  && (_z  == attr._z) && (_t  == attr._t) &&
          (_dx == attr._dx) && (_dy == attr._dy) && (_dz == attr._dz) && (_dt == attr._dt) &&
          (_xaxis[0] == attr._xaxis[0]) && (_xaxis[1] == attr._xaxis[1]) && (_xaxis[2] == attr._xaxis[2]) &&
          (_yaxis[0] == attr._yaxis[0]) && (_yaxis[1] == attr._yaxis[1]) && (_yaxis[2] == attr._yaxis[2]) &&
          (_zaxis[0] == attr._zaxis[0]) && (_zaxis[1] == attr._zaxis[1]) && (_zaxis[2] == attr._zaxis[2]) &&
          (_xorigin == attr._xorigin) && (_yorigin == attr._yorigin) && (_zorigin == attr._zorigin) && 
          (_torigin == attr._torigin));
}

inline void irtkImageAttributes::Print()
{
	
  cerr<<_x<<" "<<_y<<" "<<_z<<" "<<_t<<endl;
  cerr<<_dx<<" "<<_dy<<" "<<_dz<<" "<<_dt<<endl;
  cerr<<_xorigin<<" "<<_yorigin<<" "<<_zorigin<<" "<<_torigin<<endl;
  cerr<<_xaxis[0]<<" "<<_xaxis[1]<<" "<<_xaxis[2]<<endl;
  cerr<<_yaxis[0]<<" "<<_yaxis[1]<<" "<<_yaxis[2]<<endl;
  cerr<<_zaxis[0]<<" "<<_zaxis[1]<<" "<<_zaxis[2]<<endl;
}

inline int irtkImageAttributes::LatticeToIndex(int i, int j, int k, int l) const
{
  return l*_z*_y*_x + k*_y*_x + j*_x + i;
}

inline void irtkImageAttributes::IndexToLattice(int index, int *i, int *j, int *k, int *l)
{
    if(l != NULL){
        *l = index/(_x*_y*_z);
    }
	*k = index%(_x*_y*_z)/(_y*_x);
	*j = index%(_x*_y*_z)%(_y*_x)/_x;
	*i = index%(_x*_y*_z)%(_y*_x)%_x;
}


#endif
