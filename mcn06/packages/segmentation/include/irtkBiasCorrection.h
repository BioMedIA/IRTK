/*=========================================================================

  Library   : Image Registration Toolkit (IRTK)
  Module    : $Id$
  Copyright : Imperial College, Department of Computing
              Visual Information Processing (VIP), 2008 onwards
  Date      : $Date$
  Version   : $Revision$
  Changes   : $Author$

=========================================================================*/

#ifndef _IRTKBIASCORRECTION_H

#define _IRTKBIASCORRECTION_H

#include <irtkImage.h>

#include <irtkResampling.h>

#include <irtkTransformation.h>

#include <irtkBiasField.h>

class irtkBiasCorrection : public irtkObject
{

protected:

  irtkRealImage *_target;

  irtkRealImage *_reference;

  irtkRealImage *_weights;

  /// Output
  irtkBiasField *_biasfield;

  /// Padding value
  irtkGreyPixel _Padding;

  /// Initial set up for the registration
  virtual void Initialize();

  /// Final set up for the registration
  virtual void Finalize();

public:

  /// Constructor
  irtkBiasCorrection();

  /// Destructor
  virtual ~irtkBiasCorrection();

  /// Sets input for the bias correction filter
  virtual void SetInput (irtkRealImage *, irtkRealImage *);

  /// Sets weights for the bias correction filter
  virtual void SetWeights (irtkRealImage *);

  /// Sets output for the bias correction filter
  virtual void SetOutput(irtkBiasField *);

  /// Runs the bias correction filter
  virtual void Run();

  /// Apply bias correction to _input
  virtual void Apply(irtkRealImage &);

  /// Apply bias correction to any image
  virtual void ApplyToImage(irtkRealImage &);

  /// Apply bias correction to any image including logarithmic transform
  virtual void ApplyToImage(irtkGreyImage &);

  // Access parameters
  virtual SetMacro(Padding,      short);
  virtual GetMacro(Padding,      short);

  /// Returns the name of the class
  virtual const char *NameOfClass();
};

inline void irtkBiasCorrection::SetInput(irtkRealImage *target, irtkRealImage *reference)
{
  _target    = target;
  _reference = reference;
}

inline void irtkBiasCorrection::SetWeights(irtkRealImage *weights)
{
  _weights    = weights;
}

inline void irtkBiasCorrection::SetOutput(irtkBiasField *biasfield)
{
  _biasfield = biasfield;
}

inline const char *irtkBiasCorrection::NameOfClass()
{
  return "irtkBiasCorrection";
}

#endif
