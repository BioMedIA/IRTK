/*=========================================================================

  Library   : Image Registration Toolkit (IRTK)
  Module    : $Id$
  Copyright : Imperial College, Department of Computing
              Visual Information Processing (VIP), 2011 onwards
  Date      : $Date$
  Version   : $Revision$
  Changes   : $Author$

=========================================================================*/



#include "irtkCRF.h"

irtkCRF::irtkCRF( irtkGenericImage<pixel_t> &img,
                  irtkGenericImage<LabelID> &labels,
                  irtkGenericImage<irtkRealPixel> &proba )
 : _img(img), _labels(labels), _proba(proba)
{
    _lambda = 1.0;
}

irtkCRF::~irtkCRF()
{ }

void irtkCRF::Run()
{
    try {
        SiteID nb_pixels = _img.GetX() * _img.GetY() * _img.GetZ();
        LabelID nb_labels = 2;
        GCoptimizationGeneralGraph graph(nb_pixels, nb_labels);
    
        EnergyTermType* datacost = new EnergyTermType[nb_pixels*nb_labels];
        double epsilon = 0.00000001;
        double p1, p2;
        for ( int z = 0; z < _img.GetZ(); z++ )
            for ( int y = 0; y < _img.GetY(); y++ )
                for ( int x = 0; x < _img.GetX(); x++ ) {
                    SiteID id = this->index(z,y,x);
                    graph.setLabel( id, _labels(x,y,z) );
                    p1 = epsilon + 1 - _proba(x,y,z);
                    p2 = epsilon + _proba(x,y,z);
                    if (p1 > 1.0) p1 = 1.0;
                    if (p1 < epsilon) p1 = epsilon;
                    if (p2 > 1.0) p2 = 1.0;
                    if (p2 < epsilon) p2 = epsilon;
                    datacost[id*2+0] = -log( p1 )*1000;
                    datacost[id*2+1] = -log( p2 )*1000;
                }

        std::cout << "will set datacost\n";
               
        graph.setDataCost(datacost);

        std::cout << "done\n";
    
        EnergyTermType* labelcost = new EnergyTermType[nb_labels*nb_labels];
        for ( LabelID l1 = 0; l1 < nb_labels; l1++ )
            for ( LabelID l2 = 0; l2 < nb_labels; l2++ )
                labelcost[l1+nb_labels*l2] = double(l1 != l2) * _lambda;

        std::cout << "will set smooth cost\n";
        graph.setSmoothCost(labelcost);
        std::cout << "done\n";

        int x0, y0, z0;
        SiteID id, id0;
        double std2 = 2.0 * pow(_img.GetSD(),2);
        for ( int z = 0; z < _img.GetZ(); z++ )
            for ( int y = 0; y < _img.GetY(); y++ )
                for ( int x = 0; x < _img.GetX(); x++ ) 
                    for ( int a = -1; a <= 1; a++ )
                        for ( int b = -1; b <= 1; b++ )
                            for ( int c = -1; c <= 1; c++ ) {
                                z0 = z+a;
                                y0 = y+b;
                                x0 = x+c;
                                id = index(z,y,x);
                                id0 = index(z0,y0,x0);
                                if ( 0 <= z0 && z0 < _img.GetZ() &&
                                     0 <= y0 && y0 < _img.GetY() &&
                                     0 <= x0 && x0 < _img.GetX() &&
                                     id < id0 ) { 
                                    double dist = sqrt( pow(double(a),2.0) +
                                                        pow(double(b),2.0) +
                                                        pow(double(c),2.0) );
                                    double w = exp( -pow( double(_img(x,y,z) -
                                                          _img(x0,y0,z0)), 2.0)
                                                    / std2 )
                                        / dist;
                                    graph.setNeighbors( id, id0, w*1000 );
                                }
                            }

        std::cout << "energy before expansion: " << graph.compute_energy() << "\n";
        graph.expansion();
        std::cout << "energy after expansion: " << graph.compute_energy() << "\n";

        int l;
        for ( int z = 0; z < _img.GetZ(); z++ )
            for ( int y = 0; y < _img.GetY(); y++ )
                for ( int x = 0; x < _img.GetX(); x++ ) {
                    SiteID id = index(z,y,x);
                    l = graph.whatLabel(id);
                    if (l == 1)
                        _labels(x,y,z) = 1;
                    else
                        _labels(x,y,z) = 0;
                }

        // clean
        delete datacost;
        delete labelcost;

    }
    catch (GCException e){
        e.Report();
        exit(1);
    }    
    
    return;
}
