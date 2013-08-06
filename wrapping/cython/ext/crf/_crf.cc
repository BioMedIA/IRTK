#include "_crf.h"

inline double _pow( double x, int a ) { return x*x; }

void _crf( pixel_t* img,
           int shape0, int shape1,
           double std,
           LabelID* labels,
           double* proba,
           double lambda,
           int degree ) {

    SiteID nb_pixels = shape0 * shape1;
    LabelID nb_labels = 2;
    GCoptimizationGeneralGraph graph(nb_pixels, nb_labels);
    
    EnergyTermType* datacost = new EnergyTermType[nb_pixels*nb_labels];
    double epsilon = 0.00000001;
    for ( int y = 0; y < shape0; y++ )
        for ( int x = 0; x < shape1; x++ ) {
            SiteID id = index(y,x,shape0,shape1);
            graph.setLabel( id, labels[id] );
            datacost[id*2+0] = -log( epsilon + 1 - proba[id]);
            datacost[id*2+1] = -log( epsilon + proba[id]);
        } 
               
    graph.setDataCost(datacost);
    
    EnergyTermType* labelcost = new EnergyTermType[nb_labels*nb_labels];
    for ( LabelID l1 = 0; l1 < nb_labels; l1++ )
        for ( LabelID l2 = 0; l2 < nb_labels; l2++ )
            labelcost[l1+nb_labels*l2] = double(l1 != l2) * lambda;
        //     {
        //     if ( l1 != l2 )
        //         labelcost[l1+nb_labels*l2] = 1.0;
        //     else
        //         labelcost[l1+nb_labels*l2] = 0.0;
        // }
    graph.setSmoothCost(labelcost);

    int x0, y0;
    SiteID id, id0;
    double std2 = 2.0 * _pow(std,2);
    for ( int y = 0; y < shape0; y++ )
        for ( int x = 0; x < shape1; x++ ) 
            for ( int a = -degree; a <= degree; a++ )
                for ( int b = -degree; b <= degree; b++ ) {
                    y0 = y+a;
                    x0 = x+b;
                    id = index(y,x,shape0,shape1);
                    id0 = index(y0,x0,shape0,shape1);
                    if ( 0 <= y0  && y0 < shape0
                         && 0 <= x0 && x0 < shape1
                         && id < id0 ) { 
                        double dist = sqrt( _pow(double(a),2) + _pow(double(b),2) );
                        double w = exp( -_pow(img[id] - img[id0],2)/std2 )/dist;
                        graph.setNeighbors( id, id0, w );
                    }
                }

    std::cout << "energy before expansion: " << graph.compute_energy() << "\n";
    graph.expansion();
    std::cout << "energy after expansion: " << graph.compute_energy() << "\n";
    

    for ( int y = 0; y < shape0; y++ )
        for ( int x = 0; x < shape1; x++ ) {
            SiteID id = index(y,x,shape0,shape1);
            labels[id] = graph.whatLabel(id);
        }

    // clean
    delete datacost;
    delete labelcost;
    
    return;
}

void _crf3D( pixel_t* img,
          int shape0, int shape1, int shape2,
          double std,
          LabelID* labels,
          double* proba,
          double lambda ) {

    SiteID nb_pixels = shape0 * shape1 * shape2;
    LabelID nb_labels = 2;
    GCoptimizationGeneralGraph graph(nb_pixels, nb_labels);
    
    EnergyTermType* datacost = new EnergyTermType[nb_pixels*nb_labels];
    double epsilon = 0.00000001;
    for ( int z = 0; z < shape0; z++ )
        for ( int y = 0; y < shape1; y++ )
            for ( int x = 0; x < shape2; x++ ) {
                SiteID id = index(z,y,x,shape0,shape1,shape2);
                graph.setLabel( id, labels[id] );
                datacost[id*2+0] = -log( epsilon + 1 - proba[id]);
                datacost[id*2+1] = -log( epsilon + proba[id]);
            }

    std::cout << "will set datacost\n";
               
    graph.setDataCost(datacost);

    std::cout << "done\n";
    
    EnergyTermType* labelcost = new EnergyTermType[nb_labels*nb_labels];
    for ( LabelID l1 = 0; l1 < nb_labels; l1++ )
        for ( LabelID l2 = 0; l2 < nb_labels; l2++ )
            labelcost[l1+nb_labels*l2] = double(l1 != l2) * lambda;
        //     {
        //     if ( l1 != l2 )
        //         labelcost[l1+nb_labels*l2] = 1.0;
        //     else
        //         labelcost[l1+nb_labels*l2] = 0.0;
        // }
    std::cout << "will set smooth cost\n";
    graph.setSmoothCost(labelcost);
    std::cout << "done\n";

    int x0, y0, z0;
    SiteID id, id0;
    double std2 = 2.0 * _pow(std,2);
    for ( int z = 0; z < shape0; z++ )
        for ( int y = 0; y < shape1; y++ )
            for ( int x = 0; x < shape2; x++ ) 
                for ( int a = -1; a <= 1; a++ )
                    for ( int b = -1; b <= 1; b++ )
                        for ( int c = -1; c <= 1; c++ ) {
                            // connectivity: 6 or 26 ?
                            if (abs(a) + abs(b) + abs(c) == 1)
                                continue;
                            z0 = z+a;
                            y0 = y+b;
                            x0 = x+c;
                            id = index(z,y,x,shape0,shape1,shape2);
                            id0 = index(z0,y0,x0,shape0,shape1,shape2);
                            if ( 0 <= z0  && z0 < shape0
                                 && 0 <= y0  && y0 < shape1
                                 && 0 <= x0 && x0 < shape2
                                 && id < id0 ) { 
                                double dist = sqrt( _pow(double(a),2)
                                                    + _pow(double(b),2)
                                                    + _pow(double(c),2) );
                                double w = exp( -_pow(img[id] - img[id0],2)/std2 )/dist;
                                graph.setNeighbors( id, id0, w );
                            }
                        }

    std::cout << "energy before expansion: " << graph.compute_energy() << "\n";
    graph.expansion();
    std::cout << "energy after expansion: " << graph.compute_energy() << "\n";

    int l;
    for ( int z = 0; z < shape0; z++ )
        for ( int y = 0; y < shape1; y++ )
            for ( int x = 0; x < shape2; x++ ) {
                SiteID id = index(z,y,x,shape0,shape1,shape2);
                l = graph.whatLabel(id);
                if (l == 1)
                    labels[id] = 1;
                else
                    labels[id] = 0;
            }

    // clean
    delete datacost;
    delete labelcost;
    
    return;
}
