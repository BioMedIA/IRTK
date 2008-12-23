/*=========================================================================

  Library   : Image Registration Toolkit (IRTK)
  Module    : $Id$
  Copyright : Imperial College, Department of Computing
              Visual Information Processing (VIP), 2008 onwards
  Date      : $Date$
  Version   : $Revision$
  Changes   : $Author$

=========================================================================*/

#ifndef _IRTKUNIFORMNOISE_H

#define _IRTKUNIFORMNOISE_H

/**
 * Class for adding uniform noise to images
 *
 * This class implements a filter for adding uniformly distributed noise to
 * images.
 *
 */

template <class VoxelType> class irtkUniformNoise : public irtkNoise<VoxelType>
{

protected:

  /// Returns the name of the class
  virtual const char *NameOfClass();

public:

  /// Constructor
  irtkUniformNoise(double amplitude = 1);

  /// Destructor (empty).
  ~irtkUniformNoise() {};

  /// Run uniform noise filter
  virtual double Run(int, int, int, int);

};

#endif
