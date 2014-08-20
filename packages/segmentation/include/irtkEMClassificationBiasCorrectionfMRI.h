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
