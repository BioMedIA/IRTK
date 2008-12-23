/*=========================================================================

  Library   : Image Registration Toolkit (IRTK)
  Module    : $Id$
  Copyright : Imperial College, Department of Computing
              Visual Information Processing (VIP), 2008 onwards
  Date      : $Date$
  Version   : $Revision$
  Changes   : $Author$

=========================================================================*/

#ifndef _IRTKIMAGE_H

#define _IRTKIMAGE_H

/// Image contrast modes
typedef  enum { LightDark, DarkLight, AnyContrast } irtkContrastMode;

/// Image differentiation modes
typedef enum { F_x, F_y, F_z,
               F_x_f, F_y_f, F_z_f,
               F_x_b, F_y_b, F_z_b,
               F_xx, F_yy, F_zz,
               F_xy, F_yx, F_xz, F_zx, F_yz, F_zy,
               F_xxyy, F_yyxx, F_xxzz, F_zzxx, F_yyzz, F_zzyy,
               F_v, F_vv, F_w, F_ww, F_vw, F_wv,
               F_v_2D,  F_vv_2D, F_w_2D,  F_ww_2D, F_vw_2D, F_wv_2D,
               F_u_3D,  F_uu_3D, F_v_3D,  F_vv_3D, F_w_3D,  F_ww_3D,
               F_uv_3D, F_vu_3D, F_uw_3D, F_wu_3D, F_vw_3D, F_wv_3D,
               F_Gradient_2D, F_Gradient_3D,
               F_Laplace_2D,  F_Laplace_3D
             } irtkDiffMode;

// Basic pixel class
#include <irtkBytePixel.h>
#include <irtkGreyPixel.h>
#include <irtkRealPixel.h>
#include <irtkRGBPixel.h>
#include <irtkVector3D.h>

// Basic image class
template <class Type> class irtkGenericImage;

// Includes
#include <irtkGeometry.h>

#include <irtkBaseImage.h>
#include <irtkGenericImage.h>

/// Unsigned char image
typedef class irtkGenericImage<irtkBytePixel> irtkByteImage;
/// Short image
typedef class irtkGenericImage<irtkGreyPixel> irtkGreyImage;
/// Float image
typedef class irtkGenericImage<irtkRealPixel> irtkRealImage;
/// RGB image
typedef class irtkGenericImage<irtkRGBPixel> irtkRGBImage;
/// 3D vector images.
typedef class irtkGenericImage<irtkVector3D<char> > irtkVector3DCharImage;
typedef class irtkGenericImage<irtkVector3D<short> > irtkVector3DShortImage;
typedef class irtkGenericImage<irtkVector3D<float> > irtkVector3DFloatImage;
typedef class irtkGenericImage<irtkVector3D<double> > irtkVector3DDoubleImage;

#define IRTK_VOID            0
#define IRTK_BIT             1
#define IRTK_CHAR            2
#define IRTK_UNSIGNED_CHAR   3
#define IRTK_SHORT           4
#define IRTK_UNSIGNED_SHORT  5
#define IRTK_INT             6
#define IRTK_UNSIGNED_INT    7
#define IRTK_LONG            8
#define IRTK_UNSIGNED_LONG   9
#define IRTK_FLOAT          10
#define IRTK_DOUBLE         11

#ifndef _IMPLEMENTS_GENERICIMAGE_

extern template class irtkGenericImage<irtkBytePixel>;
extern template class irtkGenericImage<irtkGreyPixel>;
extern template class irtkGenericImage<irtkRealPixel>;
extern template class irtkGenericImage<irtkRGBPixel>;

#endif

#endif
