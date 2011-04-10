/*=========================================================================

  Library   : Image Registration Toolkit (IRTK)
  Module    : $Id: irtkHistogram_1D.h 2 2008-12-23 12:40:14Z dr $
  Copyright : Imperial College, Department of Computing
              Visual Information Processing (VIP), 2008 onwards
  Date      : $Date: 2008-12-23 12:40:14 +0000 (‰∫? 23 ÂçÅ‰∫åÊú?2008) $
  Version   : $Revision: 2 $
  Changes   : $Author: dr $

=========================================================================*/

#ifndef _IRTKIMAGEHISTOGRAM_1D_H

#define _IRTKIMAGEHISTOGRAM_1D_H

/** Class for 2D histograms.
 *  This class defines and implements 2D histograms.
 */

template <class VoxelType> class irtkImageHistogram_1D : public irtkHistogram_1D<double>
{

public:
	/// Evaluate the histogram from a given image with padding value
	virtual void Evaluate(irtkGenericImage<VoxelType> *, double padding = -1);
	/// Histogram Equalization
	virtual void Equalize(VoxelType min,VoxelType max);
	/// Back project the equalized histogram to image
	virtual void BackProject(irtkGenericImage<VoxelType> *); 
};

#endif
