/*=========================================================================

  Library   : Image Registration Toolkit (IRTK)
  Module    : $Id$
  Copyright : Imperial College, Department of Computing
              Visual Information Processing (VIP), 2008 onwards
  Date      : $Date$
  Version   : $Revision$
  Changes   : $Author$

=========================================================================*/

#ifndef _IRTKRVIEWCONFIG_H

#define _IRTKRVIEWCONFIG_H

typedef struct {

  double xmin;
  double ymin;
  double xmax;
  double ymax;
  irtkViewerMode mode;

} irtkRViewConfig;

typedef enum { _View_XY, _View_XZ, _View_YZ,
               _View_XY_XZ_v, _View_XY_YZ_v, _View_XZ_YZ_v,
               _View_XY_XZ_h, _View_XY_YZ_h, _View_XZ_YZ_h,
               _View_XY_XZ_YZ
             } irtkConfigViewerMode;


/// RView configuration with single reslice plane in the X-Y plane
extern irtkRViewConfig View_XY[];

/// RView configuration with single reslice plane in the X-Z plane
extern irtkRViewConfig View_XZ[];

/// RView configuration with single reslice plane in the Y-Z plane
extern irtkRViewConfig View_YZ[];

/// RView configuration with two vertical reslice planes (X-Y and X-Z plane).
extern irtkRViewConfig View_XY_XZ_v[];

/// RView configuration with two vertical reslice planes (X-Y and Y-Z plane).
extern irtkRViewConfig View_XY_YZ_v[];

/// RView configuration with two vertical reslice planes (X-Z and Y-Z plane).
extern irtkRViewConfig View_XZ_YZ_v[];

/// RView configuration with two horizontal reslice planes (X-Y and Z-X plane).
extern irtkRViewConfig View_XY_XZ_h[];

/// RView configuration with two horizontal reslice planes (X-Y and Y-Z plane).
extern irtkRViewConfig View_XY_YZ_h[];

/// RView configuration with two horizontal reslice planes (X-Z and Y-Z plane).
extern irtkRViewConfig View_XZ_YZ_h[];

/// RView configuration with three reslice planes
extern irtkRViewConfig View_XY_XZ_YZ[];

#endif
