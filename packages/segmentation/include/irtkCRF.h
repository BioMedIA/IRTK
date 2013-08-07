#include <stdio.h>
#include <iostream>
#include <cmath>

#include "irtkImage.h"

#include "GCoptimization.h"

typedef double pixel_t;

typedef short LabelID;
typedef int SiteID;
typedef int EnergyTermType;
typedef long long EnergyType;

class irtkCRF : public irtkObject
{
private:
    irtkGenericImage<pixel_t> &_img;
    irtkGenericImage<LabelID> &_labels;
    irtkGenericImage<double> &_proba;
    double _lambda;

    size_t index( size_t i, size_t j, size_t k );
    
public:

    /// Constructor
    irtkCRF( irtkGenericImage<pixel_t> &img,
             irtkGenericImage<LabelID> &labels,
             irtkGenericImage<double> &proba );

    /// Destructor
    ~irtkCRF();
    
    void SetLambda( double lambda );

    void Run();
};

inline size_t irtkCRF::index( size_t z, size_t y, size_t x ) {
    return x + _img.GetX()*( y + _img.GetY()*z );
}

inline void irtkCRF::SetLambda( double lambda ) {
    _lambda = lambda;
}
