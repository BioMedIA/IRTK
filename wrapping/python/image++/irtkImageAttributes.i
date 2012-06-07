extern class irtkImageAttributes
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

  /// Get Index from Lattice
  int LatticeToIndex(int i, int j, int k, int l = 0) const;

%extend
{
  %apply int* OUTPUT {int* iOUTPUT};
  %apply int* OUTPUT {int* jOUTPUT};
  %apply int* OUTPUT {int* kOUTPUT};
  %apply int* OUTPUT {int* lOUTPUT};

  %feature("docstring", "IndexToLattice() -> (min, max)

  Returns a 2-tuple representing the min and max values of the image.");

  /// Get Index from Lattice
  virtual void IndexToLattice(int index, int* iOUTPUT, int* jOUTPUT, int* kOUTPUT, int* lOUTPUT)
  {
    int i, j, k, l=0;
    self->IndexToLattice(index, &i, &j, &k, &l);
    *iOUTPUT = i;
    *jOUTPUT = j;
    *kOUTPUT = k;
    *lOUTPUT = l;
  }
}

  /// Print attributes
  void Print();

};
