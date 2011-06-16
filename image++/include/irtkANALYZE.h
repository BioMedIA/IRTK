/*=========================================================================

  Library   : Image Registration Toolkit (IRTK)
  Module    : $Id$
  Copyright : Imperial College, Department of Computing
              Visual Information Processing (VIP), 2008 onwards
  Date      : $Date$
  Version   : $Revision$
  Changes   : $Author$

=========================================================================*/

#ifndef _IRTKANALYZE_H

#define _IRTKANALYZE_H

#define ANALYZE_NONE                   0
#define ANALYZE_UNKNOWN                0
#define ANALYZE_BINARY                 1
#define ANALYZE_UNSIGNED_CHAR          2
#define ANALYZE_SIGNED_SHORT           4
#define ANALYZE_SIGNED_INT             8
#define ANALYZE_FLOAT                  16
#define ANALYZE_COMPLEX                32
#define ANALYZE_DOUBLE                 64
#define ANALYZE_RGB                    128
#define ANALYZE_ALL                    255

class irtkANALYZEHeader
{

public:

  int   sizeof_hdr;
  char  padding1[28];
  int   extents;
  short padding2;
  char  regular;
  char  padding3;
  short dims[5];
  short padding4[11];
  short data_type;
  short bits;
  short padding5;
  float pixdims[4];
  float padding6[12];
  int   glmax;
  int   glmin;
  char  padding7[168];
  int   padding8[8];

  /// Constructor
  irtkANALYZEHeader();

  /// Write header
  void Write(char *);

  /// Read header
  void Read(char *);
};

inline irtkANALYZEHeader::irtkANALYZEHeader()
{
  int i;

  // Initialize header
  sizeof_hdr=348;
  extents=16384;
  regular=114;
  data_type=4;
  bits=16;
  glmin = 0;
  glmax = 0;

  for (i = 0;i < 4;i++)
    dims[i] = 0;
  for (i = 0;i < 4;i++)
    pixdims[i] = 0;

  for (i = 0; i < 28; i++) padding1[i] = 0;
  padding2 = 0;
  padding3 = 0;
  for (i = 0; i < 11; i++) padding4[i] = 0;
  padding5 = 0;
  for (i = 0; i < 12; i++) padding6[i] = 0;
  for (i = 0; i < 168; i++) padding7[i] = 0;
  for (i = 0; i < 8; i++) padding8[i] = 0;
}

inline void irtkANALYZEHeader::Write(char *filename)
{
  irtkCofstream to;

  // Write header
  to.Open(filename);
  to.WriteAsInt(sizeof_hdr, 0);
  to.WriteAsChar(padding1, 28, 4);
  to.WriteAsInt(extents, 32);
  to.WriteAsShort(padding2, 36);
  to.WriteAsChar(regular, 38);
  to.WriteAsChar(padding3, 39);
  to.WriteAsShort(dims, 4, 40);
  to.WriteAsShort(padding4, 11, 48);
  to.WriteAsShort(data_type, 70);
  to.WriteAsShort(bits, 72);
  to.WriteAsShort(padding5, 74);
  to.WriteAsFloat(pixdims, 4, 76);
  to.WriteAsFloat(padding6, 12, 92);
  to.WriteAsInt(glmax, 140);
  to.WriteAsInt(glmin, 144);
  to.WriteAsChar(padding7, 168, 148);
  to.WriteAsInt(padding8, 8, 316);
  to.Close();
}

inline void irtkANALYZEHeader::Read(char *filename)
{
  irtkCifstream from;

  // Write header
  from.Open(filename);
  from.ReadAsInt(&sizeof_hdr, 1, 0);
  if (sizeof_hdr != 348) {
    from.IsSwapped(!from.IsSwapped());
  }
  from.ReadAsChar(padding1, 28, 4);
  from.ReadAsInt(&extents, 1, 32);
  from.ReadAsShort(&padding2, 1, 36);
  from.ReadAsChar(&regular, 1, 38);
  from.ReadAsChar(&padding3, 1, 39);
  from.ReadAsShort(dims, 4, 40);
  from.ReadAsShort(padding4, 11, 48);
  from.ReadAsShort(&data_type, 1, 70);
  from.ReadAsShort(&bits, 1, 72);
  from.ReadAsShort(&padding5, 1, 74);
  from.ReadAsFloat(pixdims, 4, 76);
  from.ReadAsFloat(padding6, 12, 92);
  from.ReadAsInt(&glmax, 1, 140);
  from.ReadAsInt(&glmin, 1, 144);
  from.ReadAsChar(padding7, 168, 148);
  from.ReadAsInt(padding8, 8, 316);
  from.Close();
}
#endif
