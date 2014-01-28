/*=========================================================================

  Library   : Image Registration Toolkit (IRTK)
  Module    : $Id$
  Copyright : Imperial College, Department of Computing
              Visual Information Processing (VIP), 2008 onwards
  Date      : $Date$
  Version   : $Revision$
  Changes   : $Author$

Copyright (c) 1999-2014 and onwards, Imperial College London
All rights reserved.
See LICENSE for details

=========================================================================*/

#ifndef _IRTKPOINTTOIMAGE_H

#define _IRTKPOINTTOIMAGE_H

/**
 * Class for point to image filter.
 *
 * This class uses a point set to produce an image as output. The filter
 * loops through point of the input point set and sets the nearest voxel
 * in the output image to 1. All voxels of the output image are initially
 * set to 0.
 */

template <class VoxelType> class irtkPointToImage : public irtkObject
{

protected:

  /// Input for the image to point filter
  irtkPointSet *_input;

  /// Output for the image to point filter
  irtkGenericImage<VoxelType> *_output;

  /// Debugging flag
  bool _DebugFlag;

  /// Flag to use world or image coordinates for point locations
  bool _UseWorldCoordinates;

public:

  /// Constructor (using world coordinates by default)
  irtkPointToImage(bool = true);

  // Deconstuctor
  virtual ~irtkPointToImage();

  /// Set input
  virtual void SetInput (irtkPointSet  *);

  /// Set output
  virtual void SetOutput(irtkGenericImage<VoxelType> *);

  /// Run image to point filter
  virtual void Run();

  /// Set debugging flag
  SetMacro(DebugFlag, bool);

  /// Get debugging flag
  GetMacro(DebugFlag, bool);

  /// Print debugging messages if debugging is enabled
  virtual void Debug(char *);
};

#include <irtkFill.h>

#endif
