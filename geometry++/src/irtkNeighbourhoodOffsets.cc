/*=========================================================================

  Library   : Image Registration Toolkit (IRTK)
  Module    : $Id$
  Copyright : Imperial College, Department of Computing
  Visual Information Processing (VIP), 2008 onwards
  Date      : $Date$
  Version   : $Revision$
  Changes   : $Author$

  =========================================================================*/


#include <irtkGeometry.h>

irtkNeighbourhoodOffsets::irtkNeighbourhoodOffsets()
{
  this->_Connectivity = CONNECTIVITY_26;
  for (int i = 0; i < 26; ++i)
    this->_Offsets[i] = 0;

  this->_Size = 0;
}

irtkNeighbourhoodOffsets::irtkNeighbourhoodOffsets(irtkBaseImage* image, irtkConnectivityType connectivity)
{
  this->Initialize(image, connectivity);
}

void irtkNeighbourhoodOffsets::Initialize(irtkBaseImage* image, irtkConnectivityType connectivity)
{
  int xdim, ydim;

  xdim = image->GetX();
  ydim = image->GetY();

  this->Initialize(xdim, ydim, connectivity);
}

void irtkNeighbourhoodOffsets::Initialize(int xdim, int ydim, irtkConnectivityType connectivity)
{
  switch (connectivity){
  case CONNECTIVITY_04:
    this->_Size = 4;
    break;
  case CONNECTIVITY_06:
    this->_Size = 6;
    break;
  case CONNECTIVITY_18:
    this->_Size = 18;
    break;
  case CONNECTIVITY_26:
    this->_Size = 26;
    break;
  default:
    cerr << "irtkNeighbourhoodOffsets: Invalid connectivty type." << endl;
    exit(1);
  }

  this->_Connectivity = connectivity;

  // Bare minimum, 4 connectivity.
  // Get the face neighbours along X and Y.
  this->_Offsets[0] =  1;
  this->_Offsets[1] = -1;
  this->_Offsets[2] = xdim * -1;
  this->_Offsets[3] = xdim;

  if (this->_Connectivity == CONNECTIVITY_06 ||
      this->_Connectivity == CONNECTIVITY_18 ||
      this->_Connectivity == CONNECTIVITY_26){
    // Add the face neighbours along Z.
    this->_Offsets[4] = xdim * ydim * -1;
    this->_Offsets[5] = xdim * ydim;
  }

  if (this->_Connectivity == CONNECTIVITY_18 ||
      this->_Connectivity == CONNECTIVITY_26){
    // Add the edge neighbours.
    this->_Offsets[6]  =  1 + this->_Offsets[2];
    this->_Offsets[7]  =  1 + this->_Offsets[3];
    this->_Offsets[8]  =  1 + this->_Offsets[4];
    this->_Offsets[9]  =  1 + this->_Offsets[5];

    this->_Offsets[10] = -1 + this->_Offsets[2];
    this->_Offsets[11] = -1 + this->_Offsets[3];
    this->_Offsets[12] = -1 + this->_Offsets[4];
    this->_Offsets[13] = -1 + this->_Offsets[5];

    this->_Offsets[14] = this->_Offsets[2] + this->_Offsets[4];
    this->_Offsets[15] = this->_Offsets[2] + this->_Offsets[5];
    this->_Offsets[16] = this->_Offsets[3] + this->_Offsets[4];
    this->_Offsets[17] = this->_Offsets[3] + this->_Offsets[5];
  }

  if (this->_Connectivity == CONNECTIVITY_26){
    // Add the vertex neighbours for the 26 neighbourhood.
    this->_Offsets[18] =  1 + this->_Offsets[14];
    this->_Offsets[19] =  1 + this->_Offsets[15];
    this->_Offsets[20] =  1 + this->_Offsets[16];
    this->_Offsets[21] =  1 + this->_Offsets[17];

    this->_Offsets[22] = -1 + this->_Offsets[14];
    this->_Offsets[23] = -1 + this->_Offsets[15];
    this->_Offsets[24] = -1 + this->_Offsets[16];
    this->_Offsets[25] = -1 + this->_Offsets[17];
  }
}


irtkNeighbourhoodOffsets::~irtkNeighbourhoodOffsets()
{
}
