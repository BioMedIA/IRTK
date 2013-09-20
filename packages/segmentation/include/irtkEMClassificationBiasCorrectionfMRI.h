/*=========================================================================

  Library   : Image Registration Toolkit (IRTK)
  Module    : $Id: irtkEMClassificationBiasCorrection.h 933 2013-07-03 15:33:19Z mm3 $
  Copyright : Imperial College, Department of Computing
              Visual Information Processing (VIP), 2008 onwards
  Date      : $Date: 2013-07-03 16:33:19 +0100 (Wed, 03 Jul 2013) $
  Version   : $Revision: 933 $
  Changes   : $Author: mm3 $

=========================================================================*/

#ifndef _IRTKEMCLASSIFICATIONBIASCORRECTIONFMRI_H

#define _IRTKEMCLASSIFICATIONBIASCORRECTIONFMRI_H

#include <irtkImage.h>

#include <irtkEMClassification.h>
#include <irtkEMClassificationBiasCorrection.h>


/*

EM maximisation algorithm with bias field estimation using Gaussian blurring and equal variance

*/

class irtkEMClassificationBiasCorrectionfMRI : public irtkEMClassificationBiasCorrection
{

protected:
public:

  irtkEMClassificationBiasCorrectionfMRI(double sigma);
  /// Estimates bias field
  virtual void BStep();
  void NormalizeInTime(irtkRealImage &image);

};
#endif
