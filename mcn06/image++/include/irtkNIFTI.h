/*=========================================================================

  Library   : Image Registration Toolkit (IRTK)
  Module    : $Id$
  Copyright : Imperial College, Department of Computing
              Visual Information Processing (VIP), 2008 onwards
  Date      : $Date$
  Version   : $Revision$
  Changes   : $Author$

=========================================================================*/

#ifndef _IRTKNIFTI_H

#define _IRTKNIFTI_H

#ifdef HAS_NIFTI

#include <nifti1_io.h>

#define NIFTI_RADIOLOGICAL        -1
#define NIFTI_NEUROLOGICAL         1
#define NIFTI_INCONSISTENT         0

/** The NIFTI header class.

   This is a wrapper around the nifti_image struct.

*/
class irtkNIFTIHeader
{

public:

  /// The "NIFTI-1" image storage struct.
  nifti_image *nim;

  /// Constructor
  irtkNIFTIHeader();

  /// Destructor
  ~irtkNIFTIHeader();

  /// Read header
  void Read(const char *);

  // Initialize header with minimal set of fields and orientation info.
  void Initialize(int, int, int, int, double, double, double, double, int, irtkMatrix &i2w);

  /// Print header (for debugging purposes)
  void Print();

};

inline irtkNIFTIHeader::irtkNIFTIHeader()
{
  nim = nifti_simple_init_nim();
}

inline irtkNIFTIHeader::~irtkNIFTIHeader()
{
  if (nim != NULL) {
    nim->data = NULL; // (de)allocated elsewhere
    nifti_image_free(nim);
  }
  nim = NULL;
}

inline void irtkNIFTIHeader::Read(const char *filename)
{
  // Be good
  if (nim != NULL) nifti_image_free(nim);

  // Use nifti1_io to read header (not data)
  int read_data = 0;
  nim = nifti_image_read(filename, read_data);

#ifdef HAS_DEBUG
  // Just checking
  this->Print();
#endif
}

inline void irtkNIFTIHeader::Initialize(int x, int y, int z, int t,
                                        double xsize, double ysize, double zsize, double tsize,
                                        int datatype,
                                        irtkMatrix &i2w)
{
  int nbytepix = 0, ss = 0;
  int i, j;

  // Copy passed fields
  nim->datatype = datatype; // Will be NIFTI_TYPE_UINT8 | NIFTI_TYPE_INT16 | NIFTI_TYPE_FLOAT32
  nim->nx = x;
  nim->ny = y;
  nim->nz = z;
  nim->dx = fabs(xsize);   // Store only absolute pixel size values
  nim->dy = fabs(ysize);   // dito
  nim->dz = fabs(zsize);   // ...
  // Redundant fields ndim->dim/pixdim will be filled by nim2nhdr/nhdr2nim conversions in Write()

  // Default values for 3D
  if (t == 0) {
    nim->nifti_type = 1;     // 1==NIFTI-1 (1 file) - should set the magic in nhdr in Write()
    nim->ndim = 3;
    nim->nt   = 1;           // So that nvox can be computed correctly, see below
    nim->nu   = 1;           // dito
    nim->nv   = 1;           // ...
    nim->nw   = 1;           // ...
    nim->dt   = 1;           // or 0???
  } else {
    nim->nifti_type = 1;     // 1==NIFTI-1 (1 file) - should set the magic in nhdr in Write()
    nim->ndim = 4;
    nim->nt   = t;           // So that nvox can be computed correctly, see below
    nim->nu   = 1;           // dito
    nim->nv   = 1;           // ...
    nim->nw   = 1;           // ...
    nim->dt   = fabs(tsize);
  }

  // Derived values
  nifti_datatype_sizes(datatype, &nbytepix, &ss);
  nim->nbyper = nbytepix;
  nim->nvox   = nim->nx * nim->ny * nim->nz * nim->nt * nim->nu * nim->nv * nim->nw;

  // Compute qform
  mat44 mat_44;

  for (i = 0; i < 4; ++i) {
    for (j = 0; j < 4; ++j) {
      mat_44.m[i][j] = i2w(i, j);
    }
  }

  nim->qform_code = 1;
  nim->qto_xyz = mat_44;
  nifti_mat44_to_quatern(
    mat_44,
    &(nim->quatern_b), &(nim->quatern_c), &(nim->quatern_d),
    &(nim->qoffset_x), &(nim->qoffset_y), &(nim->qoffset_z),
    &(nim->dx)       , &(nim->dy)       , &(nim->dz)       ,
    &(nim->qfac));

  // Don't use sform
  nim->sform_code = 0;

  // Set the units to mm and seconds
  nim->xyz_units  = NIFTI_UNITS_MM;
  nim->time_units = NIFTI_UNITS_SEC;
  nim->toffset = 0;

  // Other stuff I could think of:
  nim->scl_slope = 1;
  nim->scl_inter = 0;
  nim->intent_code = 0;
  nim->intent_p1 = nim->intent_p2 = nim->intent_p3 = 0;

  // Data will be written out by ImageToFile class
  nim->data = NULL;
}

inline void irtkNIFTIHeader::Print()
{
  cerr << "irtkNIFTIHeader::Print() : \n";
  if (nim != NULL) nifti_image_infodump(nim);
}


#endif

#endif
