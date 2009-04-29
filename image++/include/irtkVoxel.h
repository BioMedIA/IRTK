#ifndef IRTKVOXEL_H_

#define IRTKVOXEL_H_

typedef unsigned char  irtkBytePixel;
typedef short          irtkGreyPixel;
typedef float          irtkRealPixel;

#define MIN_GREY std::numeric_limits<short>::min()
#define MAX_GREY std::numeric_limits<short>::max()

#define IRTK_VOXEL_UNKNOWN  				0
#define IRTK_VOXEL_CHAR     				1
#define IRTK_VOXEL_UNSIGNED_CHAR    2
#define IRTK_VOXEL_SHORT    				3
#define IRTK_VOXEL_UNSIGNED_SHORT   4
#define IRTK_VOXEL_INT      				5
#define IRTK_VOXEL_UNSIGNED_INT     6
#define IRTK_VOXEL_FLOAT    				7
#define IRTK_VOXEL_DOUBLE   				8
#define IRTK_VOXEL_RGB      				9

template <typename T> struct voxel_limits {

public:
	
  static T min() throw();
  static T max() throw();
};

template <> struct voxel_limits<char> {

	static double min() throw() { return static_cast<char>(0x80); }
	static double max() throw() { return static_cast<char>(0x7f); }
	
};

template <> struct voxel_limits<unsigned char> {

	static double min() throw() { return static_cast<unsigned char>(0x80); }
	static double max() throw() { return static_cast<unsigned char>(0xffu); }
	
};

template <> struct voxel_limits<short> {

	static double min() throw() { return static_cast<short>(0x8000); }
	static double max() throw() { return static_cast<short>(0x7fff); }
	
};

template <> struct voxel_limits<unsigned short> {

	static double min() throw() { return static_cast<unsigned short>(0u); }
	static double max() throw() { return static_cast<unsigned short>(0xffffu); }
	
};

template <> struct voxel_limits<int> {

	static double min() throw() { return static_cast<int>(~(~0u >> 1)); }
	static double max() throw() { return static_cast<int>( ~0u >> 1); }
	
};

template <> struct voxel_limits<float> {

	static double min() throw() { return static_cast<float>(-1.0e+38f); }
	static double max() throw() { return static_cast<float>( 1.0e+38f); }
	
};

template <> struct voxel_limits<double> {

	static double min() throw() { return static_cast<double>(-1.0e+299); }
	static double max() throw() { return static_cast<double>( 1.0e+299); }
	
};

#include <irtkVector3D.h>

#endif /*IRTKVOXEL_H_*/
