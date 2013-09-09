/*=========================================================================

  Library   : Image Registration Toolkit (IRTK)
  Module    : $Id$
  Copyright : Imperial College, Department of Computing
              Visual Information Processing (VIP), 2011 onwards
  Date      : $Date$
  Version   : $Revision$
  Changes   : $Author$

=========================================================================*/



#ifndef _IRTKBIASCORRECTIONMASK_H

#define _IRTKBIASCORRECTIONMASK_H

#include <irtkImage.h>

#include <irtkResampling.h>

#include <irtkTransformation.h>

#include <irtkBiasField.h>
#include <irtkBiasCorrection.h>

class irtkBiasCorrectionMask : public irtkBiasCorrection
{

protected:
  irtkRealImage *_mask;

public:

  /// Constructor
  irtkBiasCorrectionMask();

  /// Destructor
  virtual ~irtkBiasCorrectionMask();

  virtual void SetMask( irtkRealImage *);

  /// Runs the bias correction filter
  virtual void Run();

  /// Returns the name of the class
  virtual const char *NameOfClass();
};

inline const char *irtkBiasCorrectionMask::NameOfClass()
{
  return "irtkBiasCorrectionMask";
}

#endif
