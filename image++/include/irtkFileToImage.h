/*=========================================================================

  Library   : Image Registration Toolkit (IRTK)
  Module    : $Id$
  Copyright : Imperial College, Department of Computing
              Visual Information Processing (VIP), 2008 onwards
  Date      : $Date$
  Version   : $Revision$
  Changes   : $Author$

=========================================================================*/

#ifndef _IRTKFILETOIMAGE_H

#define _IRTKFILETOIMAGE_H

#define VOXEL_UNKNOWN  0
#define VOXEL_CHAR     1
#define VOXEL_U_CHAR   2
#define VOXEL_SHORT    3
#define VOXEL_U_SHORT  4
#define VOXEL_INT      5
#define VOXEL_U_INT    6
#define VOXEL_FLOAT    7
#define VOXEL_DOUBLE   8
#define VOXEL_RGB      9
#define VOXEL_VECTOR_3D_CHAR   10
#define VOXEL_VECTOR_3D_SHORT  11
#define VOXEL_VECTOR_3D_FLOAT  12
#define VOXEL_VECTOR_3D_DOUBLE 13


/**
 * Abstract base class for any general file to image filter.
 *
 * This is the abstract base class which defines a common interface for all
 * filters which take a filename (referrencing an image file) as input and
 * produce an image as output. For each file format a derived class should
 * be created. Each derived class has to implement abstract member functions
 * ReadHeader() and CheckHeader().
 */

template <class VoxelType> class irtkFileToImage : protected irtkCifstream
{

  /// File name of image file
  char *_imagename;

  /// Pointer to image for output
  irtkGenericImage<VoxelType> *_output;

protected:

  /// Image dimensions (in x). Should be initialized by calling ReadHeader().
  int _x;

  /// Image dimensions (in y). Should be initialized by calling ReadHeader().
  int _y;

  /// Image dimensions (in z). Should be initialized by calling ReadHeader().
  int _z;

  /// Image dimensions (in t). Should be initialized by calling ReadHeader().
  int _t;

  /// Voxel dimensions (in x). Should be initialized by calling ReadHeader().
  double _xsize;

  /// Voxel dimensions (in y). Should be initialized by calling ReadHeader().
  double _ysize;

  /// Voxel dimensions (in z). Should be initialized by calling ReadHeader().
  double _zsize;

  /// Voxel dimensions (in t). Should be initialized by calling ReadHeader().
  double _tsize;

  /// Image origin (in x). Should be initialized by calling ReadHeader().
  double _xorigin;

  /// Image origin (in y). Should be initialized by calling ReadHeader().
  double _yorigin;

  /// Image origin (in z). Should be initialized by calling ReadHeader().
  double _zorigin;

  /// Image origin (in t). Should be initialized by calling ReadHeader().
  double _torigin;

  /// Image orientation (in x). Should be initialized by calling ReadHeader().
  double _xaxis[3];

  /// Image orientation (in y). Should be initialized by calling ReadHeader().
  double _yaxis[3];

  /// Image orientation (in y). Should be initialized by calling ReadHeader().
  double _zaxis[3];

  /// No. of bytes per voxel. Should be initialized by calling ReadHeader().
  int _bytes;

  /// Type of voxels. Should be initialized by calling ReadHeader().
  int _type;

  /// Intensity scaling parameter - slope (default: 1)
  float _slope;

  /// Intensity scaling parameter -  intercept (default: 0)
  float _intercept;

  /** Lookup table which data address (in bytes) for each slice. Should be
   *  initialized by calling ReadHeader().
   */
  int *_addr;

  /// Flag whether to reflect X axis
  int _reflectX;

  /// Flag whether to reflect Y axis
  int _reflectY;

  /// Flag whether to reflect Z axis
  int _reflectZ;

  /// Debug flag
  int _debug;

  /** Read header. This is an abstract function. Each derived class has to
   *  implement this function in order to initialize image dimensions, voxel
   *  dimensions, voxel type and a lookup table which the address for each
   *  slice.
   */
  virtual void ReadHeader() = 0;

public:

  /// Contructor
  irtkFileToImage();

  /// Destructor
  virtual ~irtkFileToImage();

  /** Static constructor. This constructor allocates a derived class which
   *  can be used to read the image file. This involves checking the file
   *  format of the file and dynamically allocating the corresonding derived
   *  class.
   */
  static irtkFileToImage *New(const char *);

  /// Set input
  virtual void SetInput (const char *);

  /// Set output
  virtual void SetOutput(irtkGenericImage<VoxelType> *);

  /// Get debug flag
  virtual int  GetDebugFlag();

  /// Put debug flag
  virtual void PutDebugFlag(int);

  /// Read entire image
  virtual void Run();

  /// Print image file information
  virtual void Print();

  /// Print  debugging info
  virtual void Debug(char *);

  // Returns the name of the class
  virtual const char *NameOfClass() = 0;
};

#include <irtkFilePGMToImage.h>
#include <irtkFileVTKToImage.h>
#include <irtkFileGIPLToImage.h>
#ifdef HAS_NIFTI
#include <irtkFileNIFTIToImage.h>
#endif
#include <irtkFileANALYZEToImage.h>

#endif
