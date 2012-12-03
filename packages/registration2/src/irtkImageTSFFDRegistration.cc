/*=========================================================================

Library   : Image Registration Toolkit (IRTK)
Module    : $Id$
Copyright : Imperial College, Department of Computing
Visual Information Processing (VIP), 2008 onwards
Date      : $Date$
Version   : $Revision$
Changes   : $Author$

=========================================================================*/

#include <irtkRegistration2.h>

#ifdef WIN32
#include <time.h>
#else
#include <sys/resource.h>
#endif

#define MAX_NO_LINE_ITERATIONS 12
#define MAX_SSD 0
#define MAX_NMI 2

extern irtkGreyImage **tmp_ttarget, *tmp_source;

irtkImageTSFFDRegistration::irtkImageTSFFDRegistration()
{
    // Print debugging information
    this->Debug("irtkImageTSFFDRegistration::irtkImageTSFFDRegistration");

    // Default optimization
    _OptimizationMethod = GradientDescent;

    // Default parameters for non-rigid registration
    _Lambda1     = 0;
    _Lambda2     = 0;
    _Lambda3     = 0;
    _Lambda3off  = true;
    _LargestSpacing = 512;
    _FinestSpacing = 1;
    _Mode        = RegisterXYZ;
    _mffd = NULL;
    _affd = NULL;
    _NumberOfModels = 0;
    _currentgradient = NULL;
    _tmp = NULL;
    _mask = NULL;
    _g = NULL;
    _h = NULL;
}

void irtkImageTSFFDRegistration::GuessParameter()
{
    int i, slices = false;
    double xsize, ysize, zsize, spacing, minsize;

    if ((_target == NULL) || (_source == NULL)) {
        cerr << "irtkImageTSFFDRegistration::GuessParameter: Target and source image not found" << endl;
        exit(1);
    }

    // Default parameters for registration
    _NumberOfLevels     = 4;
    _NumberOfBins       = 64;

    // Default parameters for optimization
    _SimilarityMeasure  = SSD;
    _Epsilon            = 0.000001;

    _periodic = true;

    // Read target pixel size
    _target[0]->GetPixelSize(&xsize, &ysize, &zsize);
    // Default target parameters
    _TargetBlurring[0]      = 0.5;
    if(xsize < _TargetBlurring[0]*2){
        _TargetBlurring[0] = xsize / 2;
    }
    if(ysize < _TargetBlurring[0]*2){
        _TargetBlurring[0] = ysize / 2;
    }
    if(zsize < _TargetBlurring[0]*2){
        _TargetBlurring[0] = zsize / 2;
    }

    _SourceBlurring[0]      = 0.5;
    if(xsize < _SourceBlurring[0]*2){
        _SourceBlurring[0] = xsize / 2;
    }
    if(ysize < _SourceBlurring[0]*2){
        _SourceBlurring[0] = ysize / 2;
    }
    if(zsize < _SourceBlurring[0]*2){
        _SourceBlurring[0] = zsize / 2;
    }
    minsize = (xsize <= ysize)   ? xsize : ysize;
    minsize = (zsize <= minsize) ? zsize : minsize;

    // Use xsize as spacing
    spacing = xsize;

    for (i = 0; i < _N_target; i++) {
        if (_target[i]->GetZ()==1) {
            slices = true;
            break;
        }
    }

    // Default target parameters
    _TargetResolution[0][0] = xsize;
    _TargetResolution[0][1] = ysize;
    _TargetResolution[0][2] = zsize;

    for (i = 1; i < _NumberOfLevels; i++) {
        _TargetBlurring[i]      = _TargetBlurring[i-1] * 2;
        _TargetResolution[i][0] = _TargetResolution[i-1][0] * 2;
        _TargetResolution[i][1] = _TargetResolution[i-1][1] * 2;
        if (slices)
            _TargetResolution[i][2] = zsize;
        else 
            _TargetResolution[i][2] = _TargetResolution[i-1][2] * 2;
    }

    // Read source pixel size
    _source->GetPixelSize(&xsize, &ysize, &zsize);
    minsize = (xsize <= ysize)   ? xsize : ysize;
    minsize = (zsize <= minsize) ? zsize : minsize;

    // Default source parameters
    _SourceResolution[0][0] = xsize;
    _SourceResolution[0][1] = ysize;
    _SourceResolution[0][2] = zsize;

    for (i = 1; i < _NumberOfLevels; i++) {
        _SourceBlurring[i]      = _SourceBlurring[i-1] * 2;
        _SourceResolution[i][0] = _SourceResolution[i-1][0] * 2;
        _SourceResolution[i][1] = _SourceResolution[i-1][1] * 2;
        if (slices)
            _SourceResolution[i][2] = zsize;
        else 
            _SourceResolution[i][2] = _SourceResolution[i-1][2] * 2;
    }

    // Default parameters for non-rigid registration
    _Lambda1            = 0.00001;
    _Lambda3            = 0.04;
    _LargestSpacing     = 512;
    _FinestSpacing      = 1;

    // Remaining parameters
    for (i = 0; i < _NumberOfLevels; i++) {
        _NumberOfIterations[i] = 100;
        _MinStep[i]            = 0.01;
        _MaxStep[i]            = pow(2.0,i);
    }

    _TargetPadding = 0;
    _SourcePadding = 0;
}

void irtkImageTSFFDRegistration::Initialize()
{
    // Print debugging information
    this->Debug("irtkImageTSFFDRegistration::Initialize");

    // Initialize base class
    this->irtkTemporalImageRegistration::Initialize();

    // Pointer to multi-level FFD
    _mffd = (irtkMultiLevelFreeFormTransformation *)_transformation;

    if(_FinestSpacing >= _LargestSpacing){
        cout << "Finest spacing larger than largest spacing!" << endl;
        exit(1);
    }

    for(int i = 0; i < _N_target; i++){
        if(_target[0]->GetX() != _target[i]->GetX() || 
            _target[0]->GetY() != _target[i]->GetY() || 
            _target[0]->GetZ() != _target[i]->GetZ()){
                cout << "target image need to be the same image at different phase" << endl;
                exit(1);
        }
    }

    this->InitializeTransformation();

    if(_SimilarityMeasure == SSD){
        _MaxSimilarity = MAX_SSD;
    }else if(_SimilarityMeasure == NMI){
        _MaxSimilarity = MAX_NMI;
    }

    _Lambda2 = _Lambda3;
}

void irtkImageTSFFDRegistration::InitializeTransformation(){
    // Print debugging information
    this->Debug("irtkImageTSFFDRegistration::InitializeTransformation()");

    // Previous transformation
    bool hasprevioustransformation = false;
    if(_mffd->NumberOfLevels() > 0){
        hasprevioustransformation = true;
    }
    // TODO Store and pop all current transformations
    // if(hasprevioustransformation)
    //    this->Transformation2Image();
    // Intialize number of models
    int i,j,k,t;
    double dx,dy,dz,dt,tdx,tdy,tdz,tdt,odx,ody,odz,odt,interval;

    interval = 0.5/(this->_N_target);

    dx = _source->GetXSize()*_LargestSpacing;
    dy = _source->GetYSize()*_LargestSpacing;
    if(_source->GetZ() > 1)
        dz = _source->GetZSize()*_LargestSpacing;
    else
        dz = 1;
    dt = interval*_LargestSpacing;

    _NumberOfModels = 0;
    _NumberOfDofs = 0;
    odx = 0;
    ody = 0;
    odz = 0;
    odt = 0;
    while(dx > _source->GetXSize()*_FinestSpacing && dy > _source->GetYSize()*_FinestSpacing
        &&(dz > _source->GetZSize()*_FinestSpacing || _source->GetZ() == 1)
        &&(dt > interval*_FinestSpacing)){

            if(dx > _source->GetXSize()*_source->GetX()/3.0){
                tdx = _source->GetXSize()*_source->GetX()/3.0;
            }else{
                tdx = dx;
            }

            if(dy > _source->GetYSize()*_source->GetY()/3.0){
                tdy = _source->GetYSize()*_source->GetY()/3.0;
            }else{
                tdy = dy;
            }

            if(_source->GetZ() > 1){
                if(dz > _source->GetZSize()*_source->GetZ()/3.0){
                    tdz = _source->GetZSize()*_source->GetZ()/3.0;
                }else{
                    tdz = dz;
                }
            }else{
                tdz = dz;
            }

            if(dt > 0.25){
                tdt = 0.25;
            }else{
                tdt = dt;
            }

            // check new transformation is different from previous
            if(tdx != odx || tdy != ody || tdz != odz || tdt != odt){

                odx = tdx;
                ody = tdy;
                odz = tdz;
                odt = tdt;

                _affd = new irtkBSplineFreeFormTransformationPeriodic(*_source, tdx, tdy, tdz, tdt);
                _affd->PeriodicOn();

                _mffd->PushLocalTransformation(_affd);
                // TODO Padding of FFD
                irtkPadding(_target, this->_TargetPadding, _affd, _N_target, _t_real);
                // Register in the x-direction only
                if (_Mode == RegisterX) {
                    for (t = 0; t < _affd->GetT(); t++){
                        for (i = 0; i < _affd->GetX(); i++) {
                            for (j = 0; j < _affd->GetY(); j++) {
                                for (k = 0; k < _affd->GetZ(); k++) {
                                    _Status sx, sy, sz;
                                    _affd->GetStatusCP(i, j, k, t, sx, sy, sz);
                                    _affd->PutStatusCP(i, j, k, t, sx, _Passive, _Passive);
                                }
                            }
                        }
                    }
                }

                // Register in the y-direction only
                if (_Mode == RegisterY) {
                    for (t = 0; t < _affd->GetT(); t++){
                        for (i = 0; i < _affd->GetX(); i++) {
                            for (j = 0; j < _affd->GetY(); j++) {
                                for (k = 0; k < _affd->GetZ(); k++) {
                                    _Status sx, sy, sz;
                                    _affd->GetStatusCP(i, j, k, t, sx, sy, sz);
                                    _affd->PutStatusCP(i, j, k, t, _Passive, sy, _Passive);
                                }
                            }
                        }
                    }
                }

                // Register in the x- and y-direction only
                if (_Mode == RegisterXY) {
                    for (t = 0; t < _affd->GetT(); t++){
                        for (i = 0; i < _affd->GetX(); i++) {
                            for (j = 0; j < _affd->GetY(); j++) {
                                for (k = 0; k < _affd->GetZ(); k++) {
                                    _Status sx, sy, sz;
                                    _affd->GetStatusCP(i, j, k, t, sx, sy, sz);
                                    _affd->PutStatusCP(i, j, k, t, sx, sy, _Passive);
                                }
                            }
                        }
                    }
                }

                _NumberOfDofs += _affd->NumberOfDOFs();
                _NumberOfModels++;
            }

            dx /= 2; dy /= 2; dt /= 2;
            if(_source->GetZ() > 1) dz /= 2;
    }
}

void irtkImageTSFFDRegistration::InitializeCoordLut()
{
    int i, j, k, n;
    double x, y, z, *ptr2latt;
    // Allocate memory for lattice coordinates
    _latticeCoordLUT = new double*[_NumberOfModels];

    for (n = 0; n < _NumberOfModels; n++) {
        _affd = (irtkBSplineFreeFormTransformationPeriodic *)_mffd->GetLocalTransformation(n);
        _latticeCoordLUT[n] = new double[_target[0]->GetNumberOfVoxels() * 3];

        ptr2latt = _latticeCoordLUT[n];

        for (k = 0; k < _target[0]->GetZ(); k++) {
            for (j = 0; j < _target[0]->GetY(); j++) {
                for (i = 0; i < _target[0]->GetX(); i++) {
                    x = i;
                    y = j;
                    z = k;
                    _target[0]->ImageToWorld(x, y, z);
                    _affd->WorldToLattice(x, y, z);
                    ptr2latt[0] = x;
                    ptr2latt[1] = y;
                    ptr2latt[2] = z;
                    ptr2latt += 3;
                }
            }
        }
    }
}

void irtkImageTSFFDRegistration::InitilizeSparsityParmeter()
{
    // Initialize lambda3
    if(_Lambda2 > 0){
        _Lambda3 = 0;
        this->Update(true);
        double norm;
        int count, i, j;;

        // Compute _currentgradient with respect to displacements
        this->irtkTemporalImageRegistration::EvaluateGradient(&norm);

        if (_affd->GetZ() == 1) {
            this->EvaluateGradient2D();
        } else {
            this->EvaluateGradient3D();
        }

        if (this->_Lambda1 > 0) {
            this->SmoothnessPenaltyGradient();
        }

        norm = 0;
        count = 0;

        i = 0;
        for (j = 0; j < _NumberOfModels; j++){
            _affd = (irtkBSplineFreeFormTransformationPeriodic*)_mffd->GetLocalTransformation(j);
            // Move along _gradient direction
            for (int k = 0; k < _affd->NumberOfDOFs(); k++) {
                if(_affd->GetStatus(k) == Active && fabs(_currentgradient[i]) > 0){
                    count ++;
                }
                i++;
            }
        }

        for(i = 0; i < _NumberOfDofs; i++){
            norm += fabs(_currentgradient[i]);
        }

        if(norm > 0){
            cout << norm << " " << count  << endl;
            norm = norm/count;
            _Lambda3 = norm*_Lambda2;
            cout << "normalized sparsity penalty with respect to finite convergence property is:" << _Lambda3 << endl;
        }else{
            _Lambda3 = 0;
            cout << "current gradient is 0!" << endl;
        }
    }
}

void irtkImageTSFFDRegistration::Initialize(int level)
{
    // Print debugging information
    this->Debug("irtkImageTSFFDRegistration::Initialize(int)");

    // Initialize base class
    this->irtkTemporalImageRegistration::Initialize(level);

    this->InitializeTransformation(level);

    cout << "Number of models used: " << _NumberOfModels << endl;

    // Allocate memory for _currentgradient vector
    _currentgradient = new double[_NumberOfDofs];
    _tmp = new double[_NumberOfDofs];
    _mask = new bool[_NumberOfDofs];

    this->InitializeCoordLut();

    this->InitilizeSparsityParmeter();
}

void irtkImageTSFFDRegistration::Finalize()
{
    // Print debugging information
    this->Debug("irtkImageTSFFDRegistration::Finalize");

    // Finalize base class
    this->irtkTemporalImageRegistration::Finalize();

    if (_g != NULL) delete []_g;
    _g = NULL;
    if (_h != NULL) delete []_h;
    _h = NULL;
}

void irtkImageTSFFDRegistration::Finalize(int level)
{
    // Print debugging information
    this->Debug("irtkImageTSFFDRegistration::Finalize(int)");

    // Finalize base class
    this->irtkTemporalImageRegistration::Finalize(level);

    // Dellocate memory for _currentgradient vector
    if(_currentgradient != NULL){
        delete []_currentgradient;
    }
    _currentgradient = NULL;

    if(_tmp != NULL){
        delete []_tmp;
    }
    _tmp = NULL;

    if(_mask != NULL){
        delete []_mask;
    }
    _mask = NULL;

    if(_latticeCoordLUT != NULL){
        for (int n = 0; n < _NumberOfModels; n++) {
            delete []_latticeCoordLUT[n];
        }
        delete []_latticeCoordLUT;
    }
    _latticeCoordLUT = NULL;
}

double irtkImageTSFFDRegistration::SparsePenalty()
{
    int t,index;
    double sparsity,norm;

    norm = 0;

    sparsity = 0;
    for (t = 0; t < _NumberOfModels; t++){
        _affd = (irtkBSplineFreeFormTransformationPeriodic *)_mffd->GetLocalTransformation(t);
        // approximate using a b spline model.
        for(index = 0; index < _affd->NumberOfDOFs(); index++){
            if(_affd->GetStatus(index) == _Active){
                sparsity += fabs(_affd->Get(index));
                norm ++;
            }
        }
    }
    if(sparsity > 0){
        return sparsity/norm;
    }else{
        return 0;
        _Lambda3off = true;
    }
}

void irtkImageTSFFDRegistration::SparsePenaltyGradient()
{
    int i,t,index;

    double sparsity;

    // Sparsity gradient
    index = 0;
    for (t = 0; t < _NumberOfModels; t++){
        _affd = (irtkBSplineFreeFormTransformationPeriodic *)_mffd->GetLocalTransformation(t);
        // approximate using a b spline model.
        for (i = 0; i < _affd->NumberOfDOFs(); i++){
            if(_affd->Get(i) != 0){
                sparsity = _affd->Get(i)/fabs(_affd->Get(i));
            }else{
                sparsity = 0;
            }

            _currentgradient[index] -= _Lambda3*sparsity;

            if(_currentgradient[index] * sparsity >= 0)
                _mask[index] = false;
            else
                _mask[index] = true;

            index++;
        }
    }
}

void irtkImageTSFFDRegistration::UpdateSource()
{
    short *ptr1;
    double x, y, z, t, t1, t2, u1, u2, v1, v2;
    int a, b, c, i, j, k, n, offset1, offset2, offset3, offset4, offset5, offset6, offset7, offset8;
    double xl, yl, zl, xt, yt, zt;

#ifdef USE_TIMING
    // Start timing
    clock_t start, end;
    double cpu_time_used;
    start = clock();
#endif

    // Generate transformed tmp image
    for (n = 0; n < _N_target; n++) {
        _transformedSource[n] = *_target[n];
    }

    // Calculate offsets for fast pixel access
    offset1 = 0;
    offset2 = 1;
    offset3 = this->_source->GetX();
    offset4 = this->_source->GetX()+1;
    offset5 = this->_source->GetX()*this->_source->GetY();
    offset6 = this->_source->GetX()*this->_source->GetY()+1;
    offset7 = this->_source->GetX()*this->_source->GetY()+this->_source->GetX();
    offset8 = this->_source->GetX()*this->_source->GetY()+this->_source->GetX()+1;
    //cout<<"irtkImageTSFFDRegistration::UpdateSource start"<<endl;
    for (n = 0; n < _N_target; n++) {
        if ((_target[n]->GetZ() == 1) && (_source->GetZ() == 1)) {
            for (j = 0; j < _target[n]->GetY(); j++) {
                for (i = 0; i < _target[n]->GetX(); i++) {
                    if (_target[n]->Get(i, j, 0) >= 0) {
                        x = i;
                        y = j;
                        z = 0;
                        _target[n]->ImageToWorld(x, y, z);
                        _mffd->Transform(x, y, z, _t_real[n]);
                        _source->WorldToImage(x, y, z);

                        // Check whether transformed point is inside volume
                        if ((x > 0) && (x < _source->GetX()-1) &&
                            (y > 0) && (y < _source->GetY()-1)) {

                                if (_InterpolationMode == Interpolation_Linear) {
                                    // Calculated integer coordinates
                                    a  = int(x);
                                    b  = int(y);

                                    // Calculated fractional coordinates
                                    t1 = x - a;
                                    u1 = y - b;
                                    t2 = 1 - t1;
                                    u2 = 1 - u1;

                                    // Linear interpolation in source image
                                    ptr1 = (short *)_source->GetScalarPointer(a, b, 0);
                                    _transformedSource[n](i, j, 0) = t1 * (u2 * ptr1[offset2] + u1 * ptr1[offset4]) + t2 * (u2 * ptr1[offset1] + u1 * ptr1[offset3]);
                                } else {
                                    // Interpolation in source image
                                    _transformedSource[n](i, j, 0) = _interpolator->Evaluate(x, y, 0);
                                }
                        } else {
                            _transformedSource[n](i, j, 0) = -1;
                        }
                    } else {
                        _transformedSource[n](i, j, 0) = -1;
                    }
                }
            }
        } else {
            for (k = 0; k < _target[n]->GetZ(); k++) {
                for (j = 0; j < _target[n]->GetY(); j++) {
                    for (i = 0; i < _target[n]->GetX(); i++) {
                        if (_target[n]->Get(i, j, k) >= 0) {
                            x = i;
                            y = j;
                            z = k;
                            _target[n]->ImageToWorld(x, y, z);
                            _mffd->Transform(x, y, z, _t_real[n]);
                            _source->WorldToImage(x, y, z);

                            // Check whether transformed point is inside volume
                            if ((x > 0) && (x < _source->GetX()-1) &&
                                (y > 0) && (y < _source->GetY()-1) &&
                                (z > 0) && (z < _source->GetZ()-1)) {

                                    if (_InterpolationMode == Interpolation_Linear) {

                                        // Calculated integer coordinates
                                        a  = int(x);
                                        b  = int(y);
                                        c  = int(z);

                                        // Calculated fractional coordinates
                                        t1 = x - a;
                                        u1 = y - b;
                                        v1 = z - c;
                                        t2 = 1 - t1;
                                        u2 = 1 - u1;
                                        v2 = 1 - v1;

                                        // Linear interpolation in source image
                                        ptr1 = (short *)_source->GetScalarPointer(a, b, c);
                                        _transformedSource[n](i, j, k) = (t1 * (u2 * (v2 * ptr1[offset2] + v1 * ptr1[offset6]) +
                                            u1 * (v2 * ptr1[offset4] + v1 * ptr1[offset8])) +
                                            t2 * (u2 * (v2 * ptr1[offset1] + v1 * ptr1[offset5]) +
                                            u1 * (v2 * ptr1[offset3] + v1 * ptr1[offset7])));
                                    } else {
                                        // Interpolation in source image
                                        _transformedSource[n](i, j, k) = _interpolator->Evaluate(x, y, z);
                                    }
                            } else {
                                _transformedSource[n](i, j, k) = -1;
                            }
                        } else {
                            _transformedSource[n](i, j, k) = -1;
                        }
                    }
                }
            }
        }
    }
    //cout<<"irtkImageTSFFDRegistration::UpdateSource end"<<endl;
#ifdef USE_TIMING
    // Stop timing
    end = clock();
    cpu_time_used = ((double) (end - start)) / CLOCKS_PER_SEC;
    cout << "CPU time for irtkImageTSFFDRegistration::UpdateSource() = " << cpu_time_used << endl;
#endif

}

void irtkImageTSFFDRegistration::UpdateSourceAndGradient()
{
    short *ptr1;
    double x, y, z, t, t1, t2, u1, u2, v1, v2, *ptr2;
    int a, b, c, i, j, k, n, offset1, offset2, offset3, offset4, offset5, offset6, offset7, offset8;
    double xl, yl, zl, xt, yt, zt;

#ifdef USE_TIMING
    // Start timing
    clock_t start, end;
    double cpu_time_used;
    start = clock();
#endif

    // Generate transformed tmp image
    for (n = 0; n < _N_target; n++) {
        _transformedSource[n] = *_target[n];
    }

    // Calculate offsets for fast pixel access
    offset1 = 0;
    offset2 = 1;
    offset3 = this->_source->GetX();
    offset4 = this->_source->GetX()+1;
    offset5 = this->_source->GetX()*this->_source->GetY();
    offset6 = this->_source->GetX()*this->_source->GetY()+1;
    offset7 = this->_source->GetX()*this->_source->GetY()+this->_source->GetX();
    offset8 = this->_source->GetX()*this->_source->GetY()+this->_source->GetX()+1;
    //cout<<"irtkImageTSFFDRegistration::UpdateSourceAndGradient start"<<endl;
    for (n = 0; n < _N_target; n++) {
        if ((_target[n]->GetZ() == 1) && (_source->GetZ() == 1)) {
            for (j = 0; j < _target[n]->GetY(); j++) {
                for (i = 0; i < _target[n]->GetX(); i++) {
                    if (_target[n]->Get(i, j, 0) >= 0) {
                        x = i;
                        y = j;
                        z = 0;
                        _target[n]->ImageToWorld(x, y, z);
                        _mffd->Transform(x, y, z, _t_real[n]);
                        _source->WorldToImage(x, y, z);

                        // Check whether transformed point is inside volume
                        if ((x > 0) && (x < _source->GetX()-1) &&
                            (y > 0) && (y < _source->GetY()-1)) {

                                if (_InterpolationMode == Interpolation_Linear) {
                                    // Calculated integer coordinates
                                    a  = int(x);
                                    b  = int(y);

                                    // Calculated fractional coordinates
                                    t1 = x - a;
                                    u1 = y - b;
                                    t2 = 1 - t1;
                                    u2 = 1 - u1;

                                    // Linear interpolation in source image
                                    ptr1 = (short *)_source->GetScalarPointer(a, b, 0);
                                    _transformedSource[n](i, j, 0) = t1 * (u2 * ptr1[offset2] + u1 * ptr1[offset4]) + t2 * (u2 * ptr1[offset1] + u1 * ptr1[offset3]);

                                    // Linear interpolation in _currentgradient image
                                    ptr2 = _sourceGradient.GetPointerToVoxels(a, b, 0, 0);
                                    _transformedSourceGradient[n](i, j, 0, 0) = t1 * (u2 * ptr2[offset2] + u1 * ptr2[offset4]) + t2 * (u2 * ptr2[offset1] + u1 * ptr2[offset3]);
                                    ptr2 = _sourceGradient.GetPointerToVoxels(a, b, 0, 1);
                                    _transformedSourceGradient[n](i, j, 0, 1) = t1 * (u2 * ptr2[offset2] + u1 * ptr2[offset4]) + t2 * (u2 * ptr2[offset1] + u1 * ptr2[offset3]);

                                } else {
                                    // Interpolation in source image
                                    _transformedSource[n](i, j, 0) = _interpolator->Evaluate(x, y, 0);

                                    // Interpolation in _currentgradient image
                                    _transformedSourceGradient[n](i, j, 0, 0) = _interpolatorGradient->Evaluate(x, y, 0, 0);
                                    _transformedSourceGradient[n](i, j, 0, 1) = _interpolatorGradient->Evaluate(x, y, 0, 1);
                                }
                        } else {
                            _transformedSource[n](i, j, 0) = -1;
                            _transformedSourceGradient[n](i, j, 0, 0) = 0;
                            _transformedSourceGradient[n](i, j, 0, 1) = 0;
                        }
                    } else {
                        _transformedSource[n](i, j, 0) = -1;
                        _transformedSourceGradient[n](i, j, 0, 0) = 0;
                        _transformedSourceGradient[n](i, j, 0, 1) = 0;
                    }
                }
            }
        } else {
            for (k = 0; k < _target[n]->GetZ(); k++) {
                for (j = 0; j < _target[n]->GetY(); j++) {
                    for (i = 0; i < _target[n]->GetX(); i++) {
                        if (_target[n]->Get(i, j, k) >= 0) {
                            x = i;
                            y = j;
                            z = k;
                            _target[n]->ImageToWorld(x, y, z);
                            _mffd->Transform(x, y, z, _t_real[n]);
                            _source->WorldToImage(x, y, z);

                            // Check whether transformed point is inside volume
                            if ((x > 0) && (x < _source->GetX()-1) &&
                                (y > 0) && (y < _source->GetY()-1) &&
                                (z > 0) && (z < _source->GetZ()-1)) {

                                    if (_InterpolationMode == Interpolation_Linear) {
                                        // Calculated integer coordinates
                                        a  = int(x);
                                        b  = int(y);
                                        c  = int(z);

                                        // Calculated fractional coordinates
                                        t1 = x - a;
                                        u1 = y - b;
                                        v1 = z - c;
                                        t2 = 1 - t1;
                                        u2 = 1 - u1;
                                        v2 = 1 - v1;

                                        // Linear interpolation in source image
                                        ptr1 = (short *)_source->GetScalarPointer(a, b, c);
                                        _transformedSource[n](i, j, k) = (t1 * (u2 * (v2 * ptr1[offset2] + v1 * ptr1[offset6]) +
                                            u1 * (v2 * ptr1[offset4] + v1 * ptr1[offset8])) +
                                            t2 * (u2 * (v2 * ptr1[offset1] + v1 * ptr1[offset5]) +
                                            u1 * (v2 * ptr1[offset3] + v1 * ptr1[offset7])));

                                        // Linear interpolation in _currentgradient image
                                        ptr2 = _sourceGradient.GetPointerToVoxels(a, b, c, 0);
                                        _transformedSourceGradient[n](i, j, k, 0) = (t1 * (u2 * (v2 * ptr2[offset2] + v1 * ptr2[offset6]) +
                                            u1 * (v2 * ptr2[offset4] + v1 * ptr2[offset8])) +
                                            t2 * (u2 * (v2 * ptr2[offset1] + v1 * ptr2[offset5]) +
                                            u1 * (v2 * ptr2[offset3] + v1 * ptr2[offset7])));
                                        ptr2 = _sourceGradient.GetPointerToVoxels(a, b, c, 1);
                                        _transformedSourceGradient[n](i, j, k, 1) = (t1 * (u2 * (v2 * ptr2[offset2] + v1 * ptr2[offset6]) +
                                            u1 * (v2 * ptr2[offset4] + v1 * ptr2[offset8])) +
                                            t2 * (u2 * (v2 * ptr2[offset1] + v1 * ptr2[offset5]) +
                                            u1 * (v2 * ptr2[offset3] + v1 * ptr2[offset7])));
                                        ptr2 = _sourceGradient.GetPointerToVoxels(a, b, c, 2);
                                        _transformedSourceGradient[n](i, j, k, 2) = (t1 * (u2 * (v2 * ptr2[offset2] + v1 * ptr2[offset6]) +
                                            u1 * (v2 * ptr2[offset4] + v1 * ptr2[offset8])) +
                                            t2 * (u2 * (v2 * ptr2[offset1] + v1 * ptr2[offset5]) +
                                            u1 * (v2 * ptr2[offset3] + v1 * ptr2[offset7])));
                                    } else {
                                        // Interpolation in source image
                                        _transformedSource[n](i, j, k) = _interpolator->Evaluate(x, y, z);

                                        // Interpolation in _currentgradient image
                                        _transformedSourceGradient[n](i, j, k, 0) = _interpolatorGradient->Evaluate(x, y, z, 0);
                                        _transformedSourceGradient[n](i, j, k, 1) = _interpolatorGradient->Evaluate(x, y, z, 1);
                                        _transformedSourceGradient[n](i, j, k, 2) = _interpolatorGradient->Evaluate(x, y, z, 2);
                                    }
                            } else {
                                _transformedSource[n](i, j, k) = -1;
                                _transformedSourceGradient[n](i, j, k, 0) = 0;
                                _transformedSourceGradient[n](i, j, k, 1) = 0;
                                _transformedSourceGradient[n](i, j, k, 2) = 0;
                            }
                        } else {
                            _transformedSource[n](i, j, k) = -1;
                            _transformedSourceGradient[n](i, j, k, 0) = 0;
                            _transformedSourceGradient[n](i, j, k, 1) = 0;
                            _transformedSourceGradient[n](i, j, k, 2) = 0;
                        }
                    }
                }
            }
        }
    }
    //cout<<"irtkImageTSFFDRegistration::UpdateSourceAndGradient end"<<endl;
#ifdef USE_TIMING
    // Stop timing
    end = clock();
    cpu_time_used = ((double) (end - start)) / CLOCKS_PER_SEC;
    cout << "CPU time for irtkImageTSFFDRegistration::UpdateSourceAndGradient() = " << cpu_time_used << endl;
#endif

}

void irtkImageTSFFDRegistration::Update(bool updateGradient)
{
    // Print debugging information
    this->Debug("irtkImageTSFFDRegistration::Update()");

    // Finalize base class
    this->irtkTemporalImageRegistration::Update(updateGradient);
}

double irtkImageTSFFDRegistration::SmoothnessPenalty()
{
    int i;
    double penalty;
    penalty = 0;
    for(i = 0; i < _NumberOfModels; i++){
        _affd = (irtkBSplineFreeFormTransformationPeriodic*)_mffd->GetLocalTransformation(i);
        if (_affd->GetZ() == 1) {
            penalty += _affd->Bending() / double(_affd->GetX()*_affd->GetY()*_affd->GetT());
        } else {
            penalty += _affd->Bending() / double(_affd->GetX()*_affd->GetY()*_affd->GetZ()*_affd->GetT());
        }
    }
    return penalty;
}

void irtkImageTSFFDRegistration::SmoothnessPenaltyGradient()
{

    // TODO
    int i,t,index;
    double norm;

    index = 0;
    for(t = 0; t < _NumberOfModels; t++){
        _affd = (irtkBSplineFreeFormTransformationPeriodic*)_mffd->GetLocalTransformation(t);

        norm = _target[0]->GetNumberOfVoxels()*_N_target / double(_affd->GetX()*_affd->GetY()*_affd->GetZ()*_affd->GetT());

        // Allocate memory
        double *tmp_currentgradient = new double[_affd->NumberOfDOFs()];

        // and initialize memory
        for (i = 0; i < _affd->NumberOfDOFs(); i++) {
            tmp_currentgradient[i] = 0.0;
        }

        // Compute _currentgradient of smoothness term
        _affd->BendingGradient(tmp_currentgradient);

        // Add _currentgradient to existing _currentgradient
        for (i = 0; i < _affd->NumberOfDOFs(); i++) {
            _currentgradient[index] += this->_Lambda1 * tmp_currentgradient[i] * norm;
            index++;
        }

        // Free memory
        delete []tmp_currentgradient;
    }
}

void irtkImageTSFFDRegistration::InitializeTransformation(int level){
    // Print debugging information
    this->Debug("irtkImageTSFFDRegistration::InitializeTransformation(level)");

    int i,j,k,t;
    double dx,dy,dz,dt,tdx,tdy,tdz,tdt,odx,ody,odz,odt,interval;

    interval = 0.5/(this->_N_target);

    dx = tmp_source->GetXSize()*_LargestSpacing;
    dy = tmp_source->GetYSize()*_LargestSpacing;
    if(_source->GetZ() > 1)
        dz = _source->GetZSize()*_LargestSpacing;
    else
        dz = 1;
    dt = interval*_LargestSpacing;

    _NumberOfModels = 0;
    _NumberOfDofs = 0;
    odx = 0;
    ody = 0;
    odz = 0;
    odt = 0;
    while(dx > _source->GetXSize()*_FinestSpacing && dy > _source->GetYSize()*_FinestSpacing
        &&(dz > _source->GetZSize()*_FinestSpacing || _source->GetZ() == 1)
        &&(dt > interval*_FinestSpacing)){

            if(dx > _source->GetXSize()*_source->GetX()/3.0){
                tdx = _source->GetXSize()*_source->GetX()/3.0;
            }else{
                tdx = dx;
            }

            if(dy > _source->GetYSize()*_source->GetY()/3.0){
                tdy = _source->GetYSize()*_source->GetY()/3.0;
            }else{
                tdy = dy;
            }

            if(_source->GetZ() > 1){
                if(dz > _source->GetZSize()*_source->GetZ()/3.0){
                    tdz = _source->GetZSize()*_source->GetZ()/3.0;
                }else{
                    tdz = dz;
                }
            }else{
                tdz = dz;
            }

            if(dt > 0.25){
                tdt = 0.25;
            }else{
                tdt = dt;
            }

            // check new transformation is different from previous
            if(tdx != odx || tdy != ody || tdz != odz || tdt != odt){

                odx = tdx;
                ody = tdy;
                odz = tdz;
                odt = tdt;

                _affd = (irtkBSplineFreeFormTransformationPeriodic *)_mffd->GetLocalTransformation(_NumberOfModels);

                _NumberOfDofs += _affd->NumberOfDOFs();
                _NumberOfModels++;
            }

            dx /= 2; dy /= 2; dt /= 2;
            if(_source->GetZ() > 1) dz /= 2;
    }
}

double irtkImageTSFFDRegistration::Evaluate()
{
    double tmp, similarity;

    // Evaluate similarity
    similarity = this->irtkTemporalImageRegistration::Evaluate();

    // Current derivative base
    _CurrentSimilarity = _MaxSimilarity - similarity;

    // minimization
    similarity = _CurrentSimilarity;

    cout << "image: " << similarity;

    // Add penalty for volume preservation
    if (this->_Lambda1 > 0) {
        tmp = this->_Lambda1*this->SmoothnessPenalty();
        cout << " + Bending: " << tmp;
        similarity += tmp;
    }

    if (this->_Lambda3 > 0 && _Lambda3off == false) {
        tmp = this->_Lambda3*this->SparsePenalty();
        cout << " + Sparsity: " << tmp;
        similarity += tmp;
    }

    cout << endl;

    //Return similarity measure + penalty terms
    return similarity;
}

void irtkImageTSFFDRegistration::EvaluateGradient2D()
{
    double basis, pos[3], t1, t2, tt;
    int i, j, n, m, i1, i2, j1, j2, k1, k2, x, y, t, index, index1, index2, index3, globaloffset, offset;
    _Status stat[3];
    double dist, dist1, dist2;
    double x1, y1, z1, x2, y2, z2, xi, yi, zi;

    globaloffset = 0;

    for(m = 0; m < _NumberOfModels; m++){
        _affd = (irtkBSplineFreeFormTransformationPeriodic *)_mffd->GetLocalTransformation(m);
        // Initialize _currentgradient to zero
        for (i = 0; i < _affd->NumberOfDOFs(); i++) {
            _currentgradient[i+globaloffset] = 0;
        }

        offset = _affd->GetX()*_affd->GetY()*_affd->GetZ()*_affd->GetT();

        // Loop over control points
        for (t = 0; t < _affd->GetT(); t++) {
            for (y = 0; y < _affd->GetY(); y++) {
                for (x = 0; x < _affd->GetX(); x++) {

                    // Get status of DoFs corresponding to the control point
                    _affd->GetStatusCP(x, y, 0, t, stat[0], stat[1], stat[2]);
                    // Check if any DoF corresponding to the control point is active
                    if ((stat[0] == _Active) || (stat[1] == _Active) || (stat[2] == _Active)) {

                        // If so, calculate bounding box of control point in image coordinates
                        // Note: temporal bounding box [t1,t2] is in world time and doesn't correspond to the indices of the target images
                        index = _affd->LatticeToIndex(x, y, 0, t);

                        // loop over all target images
                        for (n = 0; n < _N_target; n++) {

                            // t1, t2 not in lattice coordinates at the moment!!!!!!!
                            //		    _affd->BoundingBoxImage(_target[n], index, i1, j1, k1, i2, j2, k2, t1, t2, 1.0);
                            // spatial coordinates in world system
                            _affd->BoundingBoxCP(index, x1, y1, z1, t1, x2, y2, z2, t2, 1.0);

                            // transform time point of current target image to lattice coordinates and check for periodicity
                            tt = _t_real[n];
                            // map time to relative time intervall [0,1]
                            while (tt < 0)
                                tt += 1.;
                            while (tt >= 1)
                                tt -= 1.;
                            tt = _affd->TimeToLattice(tt);

                            // check whether time point of current target image is in temporal bounding box (check in lattice coord.)
                            if (    ( (t1 >= 0) && (t2 < _affd->GetT()-1) &&  (tt >= t1) && (tt<=t2) )
                                || ( (t1 <  0) && 							  ( (tt <= t2) || (tt >= t1+_affd->GetT()-1) ) )
                                || ( (t2 >= _affd->GetT()-1) && 			  ( (tt >= t1) || (tt <= t2-_affd->GetT()+1) ) ) ) {

                                    // Loop over all voxels in the current target (reference) volume
                                    for (j = 0; j < _target[n]->GetY(); j++) {
                                        for (i = 0; i < _target[n]->GetX(); i++) {

                                            // check whether point is in bounding box
                                            xi = i;
                                            yi = j;
                                            zi = 0;
                                            _target[n]->ImageToWorld(xi, yi, zi);
                                            if (   (xi>=x1) && (xi<=x2)
                                                && (yi>=y1) && (yi<=y2)
                                                && (zi>=z1) && (zi<=z2)) {

                                                    // Check whether reference point is valid
                                                    if ((_target[n]->Get(i, j, 0) >= 0) && (_transformedSource[n](i, j, 0) >= 0)) {

                                                        // Convert position from voxel coordinates to world coordinates
                                                        pos[0] = i;
                                                        pos[1] = j;
                                                        pos[2] = 0;
                                                        _target[n]->ImageToWorld(pos[0], pos[1], pos[2]);

                                                        // Convert world coordinates into lattice coordinates
                                                        _affd->WorldToLattice(pos[0], pos[1], pos[2]);

                                                        // Compute B-spline tensor product at pos
                                                        basis = _affd->B(pos[0] - x) * _affd->B(pos[1] - y);
                                                        // with time:
                                                        dist1 = tt - t;
                                                        dist2 = (dist1<0) ? ((tt - _affd->GetT()+1) - t) : ((tt + _affd->GetT()-1) - t);
                                                        dist = (abs(dist1)<abs(dist2)) ? dist1 : dist2;
                                                        basis *= _affd->B(dist);

                                                        index1 = _affd->LatticeToIndex(x,y,0,t) + globaloffset;
                                                        index2 = index1 + offset;
                                                        index3 = index2 + offset;
                                                        // Convert voxel-based _currentgradient into _currentgradient with respect to parameters (chain rule)
                                                        // NOTE: This currently assumes that the control points displacements are aligned with the world coordinate displacements
                                                        _currentgradient[index1] += basis * _similarityGradient[n](i, j, 0, 0);
                                                        _currentgradient[index2] += basis * _similarityGradient[n](i, j, 0, 1);
                                                        _currentgradient[index3] += 0;
                                                    }
                                            }
                                        }
                                    }
                            }
                        }
                    }
                }
            }
        }
        globaloffset += _affd->NumberOfDOFs();
    }
}

void irtkImageTSFFDRegistration::EvaluateGradient3D()
{
    double basis, *ptr, t1, t2, tt;
    int i, j, k, n, m, i1, i2, j1, j2, k1, k2, x, y, z, t, index,index1, index2, index3, globaloffset, offset;
    _Status stat[3];
    double dist, dist1, dist2;
    double x1, y1, z1, x2, y2, z2;

    globaloffset = 0;

    cout<<"irtkImageTSFFDRegistration::EvaluateGradient3D start"<<endl;

    for(m = 0; m < _NumberOfModels; m++){
        _affd = (irtkBSplineFreeFormTransformationPeriodic *)_mffd->GetLocalTransformation(m);
        // Initialize _currentgradient to zero
        for (i = 0; i < _affd->NumberOfDOFs(); i++) {
            _currentgradient[i+globaloffset] = 0;
        }

        offset = _affd->GetX()*_affd->GetY()*_affd->GetZ()*_affd->GetT();
        int blablub = offset;
        cout<<"number of control points: "<<blablub<<";";
        cout.flush();
        // Loop over control points
        for (t = 0; t < _affd->GetT(); t++) {
            for (z = 0; z < _affd->GetZ(); z++) {
                for (y = 0; y < _affd->GetY(); y++) {
                    for (x = 0; x < _affd->GetX(); x++) {
                        // Get status of DoFs corresponding to the control point
                        _affd->GetStatusCP(x, y, z, t, stat[0], stat[1], stat[2]);
                        // Check if any DoF corresponding to the control point is active
                        if ((stat[0] == _Active) || (stat[1] == _Active) || (stat[2] == _Active)) {

                            // If so, calculate bounding box of control point in image coordinates
                            // Note: temporal bounding box [t1,t2] is in world time and doesn't correspond to the indices of the target images
                            index = _affd->LatticeToIndex(x, y, z, t);
                            index1 = index + globaloffset;
                            index2 = index1 + offset;
                            index3 = index2 + offset;
                            // t1, t2 not in lattice coordinates at the moment!!!!!!!
                            // spatial coordinates in world system
                            //		    _affd->BoundingBoxCP(index, x1, y1, z1, t1, x2, y2, z2, t2, 1.0);
                            // t1, t2 not in lattice coordinates at the moment!!!!!!!
                            _affd->BoundingBoxImage(_target[0], index, i1, j1, k1, i2, j2, k2, t1, t2, 1.0);

                            // loop over all target images
                            for (n = 0; n < _N_target; n++) {
                                // transform time point of current target image to lattice coordinates and check for periodicity
                                tt = _t_real[n];
                                // map time to relative time intervall [0,1]
                                while (tt < 0)
                                    tt += 1.;
                                while (tt >= 1)
                                    tt -= 1.;
                                tt = _affd->TimeToLattice(tt);

                                // check whether time point of current target image is in temporal bounding box (check in lattice coord.)
                                if (( (t1 >= 0) 
                                    && (t2 < _affd->GetT()-1) 
                                    &&  (tt >= t1) 
                                    && (tt<=t2) )
                                    || ( (t1 <  0) 
                                    && ( (tt <= t2) 
                                    || (tt >= t1+_affd->GetT()-1) ) )
                                    || ( (t2 >= _affd->GetT()-1) 
                                    && ( (tt >= t1) 
                                    || (tt <= t2-_affd->GetT()+1) ) ) ) {

                                        // Loop over all voxels in the target (reference) volume
                                        for (k = k1; k <= k2; k++) {
                                            for (j = j1; j <= j2; j++) {
                                                ptr = &(_latticeCoordLUT[m]
                                                [3 * (k * (_target[0]->GetX()*_target[0]->GetY())
                                                    + j * _target[0]->GetX() + i1)]);
                                                for (i = i1; i <= i2; i++) {

                                                    // Check whether reference point is valid
                                                    if ((_target[n]->Get(i, j, k) >= 0) && (_transformedSource[n](i, j, k) >= 0)) {

                                                        // Compute B-spline tensor product at current position
                                                        basis = _affd->B(ptr[0] - x) * _affd->B(ptr[1] - y) * _affd->B(ptr[2] - z);
                                                        // with time:
                                                        dist1 = tt - t;
                                                        dist2 = (dist1<0) ? ((tt - _affd->GetT()+1) - t) : ((tt + _affd->GetT()-1) - t);
                                                        dist = (abs(dist1)<abs(dist2)) ? dist1 : dist2;
                                                        basis *= _affd->B(dist);

                                                        // Convert voxel-based _currentgradient into _currentgradient with respect to parameters (chain rule)
                                                        // NOTE: This currently assumes that the control points displacements are aligned with the world coordinate displacements
                                                        _currentgradient[index1] += basis * _similarityGradient[n](i, j, k, 0);
                                                        _currentgradient[index2] += basis * _similarityGradient[n](i, j, k, 1);
                                                        _currentgradient[index3] += basis * _similarityGradient[n](i, j, k, 2);
                                                    }
                                                    ptr += 3;
                                                }
                                            }
                                        }
                                }
                            }
                        }
                    }
                }
            }
        }
        globaloffset += _affd->NumberOfDOFs();
        cout<<endl;
    }

    cout<<"irtkImageTSFFDRegistration::EvaluateGradient3D end"<<endl;
}

void irtkImageTSFFDRegistration::NormalizeGradient()
{
    int x,y,z,t,m,index1,index2,index3,offset,globaloffset;
    double norm,spacingnorm;

    globaloffset = 0;

    for (m = 0; m < _NumberOfModels; m++){
        _affd = (irtkBSplineFreeFormTransformationPeriodic *)_mffd->GetLocalTransformation(m);
        offset = _affd->GetX()*_affd->GetY()*_affd->GetZ()*_affd->GetT();

        spacingnorm = (double(_source->GetNumberOfVoxels()*_N_target) / double(offset));

        // approximate using a b spline model.
        for (t = 0; t < _affd->GetT(); t++){
            for (z = 0; z < _affd->GetZ(); z++) {
                for (y = 0; y < _affd->GetY(); y++) {
                    for (x = 0; x < _affd->GetX(); x++) {
                        index1 = _affd->LatticeToIndex(x,y,z,t) + globaloffset;
                        index2 = index1 + offset;
                        index3 = index2 + offset;
                        norm = pow(_currentgradient[index1],2.0)
                            + pow(_currentgradient[index2],2.0)
                            + pow(_currentgradient[index3],2.0);

                        //normalize
                        if(norm > 0){
                            norm = sqrt(norm);
                            _currentgradient[index1] = _currentgradient[index1]/(norm + spacingnorm*_Epsilon);
                            _currentgradient[index2] = _currentgradient[index2]/(norm + spacingnorm*_Epsilon);
                            _currentgradient[index3] = _currentgradient[index3]/(norm + spacingnorm*_Epsilon);
                        }
                    }
                }
            }
        }

        globaloffset += _affd->NumberOfDOFs();
    }
}

double irtkImageTSFFDRegistration::EvaluateGradient()
{
    double norm, max_length;
    int i, x, y, z, t, index, index2, index3;
    double gg, dgg, gamma;

    // Compute _currentgradient with respect to displacements
    this->irtkTemporalImageRegistration::EvaluateGradient(&norm);

    // Start timing
    clock_t start, end;
    //double cpu_time_used;
    //start = clock();

    if (_affd->GetZ() == 1) {
        this->EvaluateGradient2D();
    } else {
        this->EvaluateGradient3D();
    }

    if (this->_Lambda1 > 0) {
        this->SmoothnessPenaltyGradient();
    }

    if(_Lambda3 > 0  && _Lambda3off == false){
        this->SparsePenaltyGradient();
    }

    this->NormalizeGradient();

    // Update _currentgradient
    if (_CurrentIteration == 0) {
        // First iteration, so let's initialize
        if (_g != NULL) delete []_g;
        _g = new double [_NumberOfDofs];
        if (_h != NULL) delete []_h;
        _h = new double [_NumberOfDofs];
        for (i = 0; i < _NumberOfDofs; i++) {
            _g[i] = -_currentgradient[i];
            _h[i] = _g[i];
        }
    } else {
        // Update _currentgradient direction to be conjugate
        gg = 0;
        dgg = 0;
        for (i = 0; i < _NumberOfDofs; i++) {
            gg  += _g[i]*_h[i];
            dgg += (_currentgradient[i]+_g[i])*_currentgradient[i];
        }
        gamma = dgg/gg;
        for (i = 0; i < _NumberOfDofs; i++) {
            _g[i] = -_currentgradient[i];
            _h[i] = _g[i] + gamma*_h[i];
            _currentgradient[i] = -_h[i];
        }
    }

    // Calculate maximum of _currentgradient vector
    max_length = 0;
    int m = 0;
    for(int n = 0; n < _NumberOfModels; n++){
        _affd = (irtkBSplineFreeFormTransformationPeriodic*)_mffd->GetLocalTransformation(n);
        for(t = 0; t < _affd->GetT(); t++){
            for (z = 0; z < _affd->GetZ(); z++) {
                for (y = 0; y < _affd->GetY(); y++) {
                    for (x = 0; x < _affd->GetX(); x++) {
                        index  = m+_affd->LatticeToIndex(x, y, z, t);
                        index2 = index+_affd->GetX()*_affd->GetY()*_affd->GetZ()*_affd->GetT();
                        index3 = index2+_affd->GetX()*_affd->GetY()*_affd->GetZ()*_affd->GetT();
                        norm = sqrt(_currentgradient[index] * _currentgradient[index] 
                        + _currentgradient[index2] * _currentgradient[index2] 
                        + _currentgradient[index3] * _currentgradient[index3]);
                        if (norm > max_length) 
                            max_length = norm;
                    }
                }
            }
        }

        // Deal with active and passive control points
        for (i = 0; i < _affd->NumberOfDOFs(); i++) {
            if (_affd->irtkTransformation::GetStatus(i) == _Passive) {
                _currentgradient[m+i] = 0;
                _affd->Put(i,0);
            }
        }

        m += _affd->NumberOfDOFs();
    }

    return max_length;
}

void irtkImageTSFFDRegistration::Run()
{
    int i, j, k, index;
    char buffer[256];
    double gradient,delta, step, min_step, max_step, max_length, best_similarity, new_similarity, old_similarity, norm;

    // Print debugging information
    this->Debug("irtkImageTSFFDRegistration::Run");

    if (_source == NULL) {
        cerr << "irtkImageTSFFDRegistration::Run: Filter has no source input" << endl;
        exit(1);
    }

    if (_target == NULL) {
        cerr << "irtkImageTSFFDRegistration::Run: Filter has no target input" << endl;
        exit(1);
    }

    if (_transformation == NULL) {
        cerr << "irtkImageTSFFDRegistration::Run: Filter has no transformation output" << endl;
        exit(1);
    }

    // Do the initial set up for all levels
    this->Initialize();

    // Loop over levels
    for (_CurrentLevel = _NumberOfLevels-1; _CurrentLevel >= 0; _CurrentLevel--) {
        gradient = 1.0;
        // Initial step size
        min_step = _MinStep[_CurrentLevel];
        max_step = _MaxStep[_CurrentLevel];

        // Print resolution level
        cout << "Resolution level no. " << _CurrentLevel+1 << " (step sizes " << min_step << " to " << max_step  << ")\n";

        // Initialize for this level
        this->Initialize(_CurrentLevel);
        cout<<"Initialize("<<_CurrentLevel<<") done"<<endl;

        // Save pre-processed images if we are debugging
        //    _DebugFlag = true;
        sprintf(buffer, "source_%d.nii.gz", _CurrentLevel);
        if (_DebugFlag == true) _source->Write(buffer);
        for (int n = 0; n < _N_target; n++) {
            sprintf(buffer, "target_%d_%d.nii.gz", _CurrentLevel, n);
            if (_DebugFlag == true) _target[n]->Write(buffer);
        }
        _DebugFlag = false;

        // Run the registration filter at this resolution
        _CurrentIteration = 0;
        while (_CurrentIteration < _NumberOfIterations[_CurrentLevel]) {
            cout << "Iteration = " << _CurrentIteration + 1 << " (out of " << _NumberOfIterations[_CurrentLevel] << ")"<< endl;

            // Update source image
            this->Update(true);

            // Compute current metric value
            best_similarity = old_similarity = this->Evaluate();
            cout << "Current objective function value is " << best_similarity << endl;

            _Epsilon = fabs(best_similarity) * 0.0001;

            // Compute _currentgradient of similarity metric. The function EvaluateGradient() returns the maximum control point length in the _currentgradient
            max_length = this->EvaluateGradient();

            // Step along _currentgradient direction until no further improvement is necessary
            i = 0;
            delta = 0;

            if(max_length > 0){

                if(_Lambda3 > 0){
                    // Define steps
                    norm = 0;
                    int count = 0;
                    for (j = 0; j < _NumberOfModels; j++){
                        _affd = (irtkBSplineFreeFormTransformationPeriodic*)_mffd->GetLocalTransformation(j);
                        // Move along _gradient direction
                        for (k = 0; k < _affd->NumberOfDOFs(); k++) {
                            if(_affd->GetStatus(k) == Active){
                                norm += fabs(_affd->Get(k));
                                count ++;
                            }
                        }
                    }
                    norm = norm / count;

                    if(min_step < norm){
                        min_step = norm / 128;
                        max_step = norm;
                    }
                }
            }
            step = max_step;

            do {
                double current = step / max_length;

                index = 0;
                // Move along _currentgradient direction
                for (j = 0; j < _NumberOfModels; j++){
                    _affd = (irtkBSplineFreeFormTransformationPeriodic*)_mffd->GetLocalTransformation(j);
                    // Move along _gradient direction
                    for (k = 0; k < _affd->NumberOfDOFs(); k++) {
                        _tmp[index] = _affd->Get(k);
                        if(_currentgradient[index] != 0){
                            _affd->Put(k, _affd->Get(k) + current *  _currentgradient[index]);
                            //sign changed
                            if(_mask[index] == true && _tmp[index] * _affd->Get(k) <= 0)
                                _affd->Put(k,0);
                        }
                        index++;
                    }
                }

                // We have just changed the transformation parameters, so we need to update
                this->Update(false);

                cout << "Current metric value is ";

                new_similarity = this->Evaluate();

                if (new_similarity < best_similarity - _Epsilon) {
                    cout << " = " << new_similarity << " accepted; step = " << step << endl;
                    best_similarity = new_similarity;
                    delta += step;

                    _Lambda3off = false;

                    step = step * 1.1;
                } else {
                    // Last step was no improvement, so back track
                    cout << " = " << new_similarity << " rejected; step = " << step << endl;
                    index = 0;
                    for (j = 0; j < _NumberOfModels; j++){
                        _affd = (irtkBSplineFreeFormTransformationPeriodic*)_mffd->GetLocalTransformation(j);
                        // Move along _currentgradient direction
                        for (k = 0; k < _affd->NumberOfDOFs(); k++) {
                            if(_currentgradient[index] != 0)
                                _affd->Put(k, _tmp[index]);
                            index++;
                        }
                    }
                    step = step * 0.5;

                    if(delta > 0)
                        break;
                }
                i++;
            } while ((i < MAX_NO_LINE_ITERATIONS) && (step > min_step));

            _CurrentIteration++;

            // Check for convergence
            if (delta == 0) break;
        }

        // Do the final cleaning up for this level
        this->Finalize(_CurrentLevel);
    }

    // Do the final cleaning up for all levels
    this->Finalize();
}

bool irtkImageTSFFDRegistration::Read(char *buffer1, char *buffer2, int &level)
{
    int ok = false;

    if ((strstr(buffer1, "Bending penalty") != NULL)) {
        this->_Lambda1 = atof(buffer2);
        cout << "Bending penalty is ... " << this->_Lambda1 << endl;
        ok = true;
    }

    if ((strstr(buffer1, "Sparsity constrain") != NULL)) {
        this->_Lambda3 = atof(buffer2);
        cout << "Sparsity constrain lambda is ... " << this->_Lambda3 << endl;
        ok = true;
    }

    if ((strstr(buffer1, "Coarsest spacing") != NULL)) {
        this->_LargestSpacing = atof(buffer2);
        cout << "Coarsest grid spacing is ... " << this->_LargestSpacing << endl;
        ok = true;
    }

    if ((strstr(buffer1, "Finest spacing") != NULL)) {
        this->_FinestSpacing = atof(buffer2);
        cout << "Finest grid spacing is ... " << this->_FinestSpacing << endl;
        ok = true;
    }

    if (ok == false) {
        return this->irtkTemporalImageRegistration::Read(buffer1, buffer2, level);
    } else {
        return ok;
    }
}

void irtkImageTSFFDRegistration::Write(ostream &to)
{
    to << "\n#\n# Sparse non-rigid registration parameters\n#\n\n";
    to << "Bending penalty                     = " << this->_Lambda1 << endl;
    to << "Sparsity constrain                  = " << this->_Lambda3 << endl;
    to << "Coarsest spacing                    = " << this->_LargestSpacing << endl;
    to << "Finest spacing                      = " << this->_FinestSpacing << endl;

    this->irtkTemporalImageRegistration::Write(to);
}

