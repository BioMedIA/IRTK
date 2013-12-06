/*=========================================================================

  Library   : Image Registration Toolkit (IRTK)
  Module    : $Id$
  Copyright : Imperial College, Department of Computing
              Visual Information Processing (VIP), 2011 onwards
  Date      : $Date$
  Version   : $Revision$
  Changes   : $Author$

=========================================================================*/


#include <stdio.h>
#include <iostream>
#include <cmath>

#include "irtkImage.h"

#include "GCoptimization.h"

typedef irtkRealPixel pixel_t;

typedef short LabelID;
typedef int SiteID;
typedef int EnergyTermType;
typedef long long EnergyType;

class irtkCRF : public irtkObject
{
private:
    irtkGenericImage<pixel_t> &_img;
    irtkGenericImage<LabelID> &_labels;
    irtkGenericImage<irtkRealPixel> &_proba;
    double _lambda;
    double _sigma;
    double _sigmaZ;

    size_t index( size_t i, size_t j, size_t k );
    
public:

    /// Constructor
    irtkCRF( irtkGenericImage<pixel_t> &img,
             irtkGenericImage<LabelID> &labels,
             irtkGenericImage<irtkRealPixel> &proba );

    /// Destructor
    ~irtkCRF();
    
    void SetLambda( double lambda );
    void SetSigma( double sigma );
    void SetSigmaZ( double sigmaZ );

    void Run();
};

inline size_t irtkCRF::index( size_t z, size_t y, size_t x ) {
    return x + _img.GetX()*( y + _img.GetY()*z );
}

inline void irtkCRF::SetLambda( double lambda ) {
    _lambda = lambda;
}

inline void irtkCRF::SetSigma( double sigma ) {
    _sigma = sigma;
}

inline void irtkCRF::SetSigmaZ( double sigmaZ ) {
    _sigmaZ = sigmaZ;
}
