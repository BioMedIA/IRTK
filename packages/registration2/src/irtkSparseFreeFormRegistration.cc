/*=========================================================================

Library   : Image Registration Toolkit (IRTK)
Module    : $Id$
Copyright : Imperial College, Department of Computing
Visual Information Processing (VIP), 2009 onwards
Date      : $Date$
Version   : $Revision$
Changes   : $Author$

Copyright (c) 1999-2014 and onwards, Imperial College London
All rights reserved.
See LICENSE for details

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

extern irtkRealImage *tmp_target, *tmp_source;

irtkSparseFreeFormRegistration::irtkSparseFreeFormRegistration()
{
    // Print debugging information
    this->Debug("irtkSparseFreeFormRegistration::irtkSparseFreeFormRegistration");

    // Default optimization
    _OptimizationMethod = GradientDescent;

    // Default parameters for non-rigid registration
    _Lambda1     = 0;
    _Lambda2     = 0;
    _Lambda3     = 0;
    _LargestSpacing = 512;
    _FinestSpacing = 1;
    _Mode        = RegisterXYZ;
    _mffd = NULL;
    _affd = NULL;
    _mdffd = new irtkMultiLevelFreeFormTransformation;
    _NumberOfModels = 0;
    _currentgradient = NULL;
    _tmp = NULL;
    _mask = NULL;
    _g = NULL;
    _h = NULL;
}

irtkSparseFreeFormRegistration::~irtkSparseFreeFormRegistration()
{
    // Print debugging information
    this->Debug("irtkSparseFreeFormRegistration::~irtkSparseFreeFormRegistration");
}

void irtkSparseFreeFormRegistration::GuessParameter()
{
    int i,slices;
    double xsize, ysize, zsize;
    irtkRealPixel min,max;

    if ((_target == NULL) || (_source == NULL)) {
        cerr << "irtkSparseFreeFormRegistration::GuessParameter: Target and source image not found" << endl;
        exit(1);
    }

    slices = false;

    if (_target->GetZ()==1) {
        slices = true;
    }

    // Default parameters for registration
    _NumberOfLevels     = 4;
    // Initial guess of bins with range
    _target->GetMinMax(&min,&max);
    _NumberOfBins       = round((max - min)/5.0);
    // Initial guess of bins with number of voxels
    if(_NumberOfBins > round(_target->GetNumberOfVoxels() / 1000.0)){
        _NumberOfBins = round(_target->GetNumberOfVoxels() / 1000.0);
    }
    // Don't be too small
    if(_NumberOfBins < 16){
        _NumberOfBins = 16;
    }else if(_NumberOfBins > 64){
        _NumberOfBins = 64;
    }

    // Default parameters for optimization
    _SimilarityMeasure  = SSD;
    _Epsilon            = 0.00001;
    _InterpolationMode  = Interpolation_BSpline;

    // Read target pixel size
    _target->GetPixelSize(&xsize, &ysize, &zsize);

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

    _TargetResolution[0][0] = xsize;
    _TargetResolution[0][1] = ysize;
    if (_target->GetZ() > 1) {
        _TargetResolution[0][2] = zsize;
    } else {
        _TargetResolution[0][2] = 0;
    }

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

    // Default source parameters
    _SourceResolution[0][0] = xsize;
    _SourceResolution[0][1] = ysize;
    if (_source->GetZ() > 1) {
        _SourceResolution[0][2] = zsize;
    } else {
        _SourceResolution[0][2] = 0;
    }

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
    // recommended values
    _Lambda1            = 0.0001;
    _Lambda3            = 0.04;
    _LargestSpacing     = 1024;
    _FinestSpacing      = 1;

    // Remaining parameters
    for (i = 0; i < _NumberOfLevels; i++) {
        _NumberOfIterations[i] = 100;
        _MinStep[i]            = 0.01;
        _MaxStep[i]            = pow(2.0,i);
    }

    // Try to guess padding by looking at voxel values in all eight corners of the volume:
    // If all values are the same we assume that they correspond to the padding value
    _TargetPadding = -1;

    _SourcePadding = -1;
}

void irtkSparseFreeFormRegistration::Initialize()
{
    // Print debugging information
    this->Debug("irtkSparseFreeFormRegistration::Initialize");

    // Initialize base class
    this->irtkImageRegistration2::Initialize();

    // Pointer to multi-level FFD
    _mffd = (irtkMultiLevelFreeFormTransformation *)_transformation;

    if(_FinestSpacing >= _LargestSpacing){
        cout << "Finest spacing larger than largest spacing!" << endl;
        exit(1);
    }

    this->InitializeTransformation();

    if(_SimilarityMeasure == SSD){
        _MaxSimilarity = MAX_SSD;
    }else if(_SimilarityMeasure == NMI){
        _MaxSimilarity = MAX_NMI;
    }

    _Lambda3tmp = _Lambda3;

}

void irtkSparseFreeFormRegistration::Image2Transformation(){

    cout << "Image 2 Transformation" << endl;

    // approximate start point
    int i,j,k,t,index,numberofpoints;
    double *x_pos,*y_pos,*z_pos,*x_dis,*y_dis,*z_dis;

    numberofpoints =_similarityGradient.GetNumberOfVoxels()/3;

    // generate displacement to be approximated
    x_pos = new double[numberofpoints];
    y_pos = new double[numberofpoints];
    z_pos = new double[numberofpoints];

    x_dis = new double[numberofpoints];
    y_dis = new double[numberofpoints];
    z_dis = new double[numberofpoints];

    index = 0;
    for(k = 0; k <_similarityGradient.GetZ(); k++){
        for(j = 0; j <_similarityGradient.GetY(); j++){
            for(i = 0; i <_similarityGradient.GetX(); i++){
                x_pos[index] = i;
                y_pos[index] = j;
                z_pos[index] = k;
                _similarityGradient.ImageToWorld(x_pos[index],y_pos[index],z_pos[index]);
                x_dis[index] =_similarityGradient.GetAsDouble(i, j, k, 0);
                y_dis[index] =_similarityGradient.GetAsDouble(i, j, k, 1);
                if(_similarityGradient.GetZ()>1){
                    z_dis[index] =_similarityGradient.GetAsDouble(i, j, k, 2);
                }
                else{
                    z_dis[index] = 0;
                }
                index ++;
            }
        }
    }

    for (t = 0; t < _NumberOfModels; t++){
        _affd = (irtkBSplineFreeFormTransformation *)_mffd->GetLocalTransformation(t);
        // approximate using a b spline model.
        _affd->Approximate(x_pos,y_pos,z_pos,x_dis,y_dis,z_dis,numberofpoints);
    }

    cout << "Done" << endl;

    delete []x_pos;
    delete []y_pos;
    delete []z_pos;
    delete []x_dis;
    delete []y_dis;
    delete []z_dis;

}

void irtkSparseFreeFormRegistration::Transformation2Image(){

    cout << "Transformation 2 Image" << endl;

    // Initialize
    int i,j,k;
    double x1,y1,z1,x2,y2,z2;
    // Determine attributes of source image
    irtkImageAttributes attr = _target->GetImageAttributes();

    // Set up _currentgradient of similarity metric
    attr._t = 3;
    _similarityGradient.Initialize(attr);

    // save levels to image
    for (k = 0; k < _target->GetZ(); k++) {
        for (j = 0; j < _target->GetY(); j++) {
            for (i = 0; i < _target->GetX(); i++) {
                x1 = i;
                y1 = j;
                z1 = k;
                _target->ImageToWorld(x1,y1,z1);
                x2 = x1;
                y2 = y1;
                z2 = z1;
                _mffd->LocalDisplacement(x2,y2,z2);
                _similarityGradient(i,j,k,0) = x2;
                _similarityGradient(i,j,k,1) = y2;
                _similarityGradient(i,j,k,2) = z2;
            }
        }
    }

    // now saved delete all levels
    while(_mffd->NumberOfLevels() > 0){
        irtkFreeFormTransformation *ffd = _mffd->PopLocalTransformation();
        delete ffd;
    }
}

void irtkSparseFreeFormRegistration::InitializeTransformation(){
    // Print debugging information
    this->Debug("irtkSparseFreeFormRegistration::InitializeTransformation()");

    // Previous transformation
    bool hasprevioustransformation = false;
    if(_mffd->NumberOfLevels() > 0){
        hasprevioustransformation = true;
    }
    // Store and pop all current transformations
    if(hasprevioustransformation)
        this->Transformation2Image();
    // Intialize number of models
    int i,j,k;
    double dx,dy,dz,tdx,tdy,tdz,odx,ody,odz;

    dx = _target->GetXSize()*_LargestSpacing;
    dy = _target->GetYSize()*_LargestSpacing;
    if(_target->GetZ() > 1)
        dz = _target->GetZSize()*_LargestSpacing;
    else
        dz = 1;

    _NumberOfModels = 0;
    _NumberOfDofs = 0;
    odx = 0;
    ody = 0;
    odz = 0;
    while(dx > _target->GetXSize()*_FinestSpacing && dy > _target->GetYSize()*_FinestSpacing
        &&(dz > _target->GetZSize()*_FinestSpacing || _target->GetZ() == 1)){

            if(dx > _target->GetXSize()*_target->GetX()/3.0){
                tdx = _target->GetXSize()*_target->GetX()/3.0;
            }else{
                tdx = dx;
            }

            if(dy > _target->GetYSize()*_target->GetY()/3.0){
                tdy = _target->GetYSize()*_target->GetY()/3.0;
            }else{
                tdy = dy;
            }

            if(_target->GetZ() > 1){
                if(dz > _target->GetZSize()*_target->GetZ()/3.0){
                    tdz = _target->GetZSize()*_target->GetZ()/3.0;
                }else{
                    tdz = dz;
                }
            }else{
                tdz = dz;
            }

            // check new transformation is different from previous
            if(tdx != odx || tdy != ody || tdz != odz){

                odx = tdx;
                ody = tdy;
                odz = tdz;

                _affd = new irtkBSplineFreeFormTransformation3D(*_target, tdx, tdy, tdz);
                _mdffd->PushLocalTransformation(_affd);
                // Padding of FFD
                irtkPadding(*_target, this->_TargetPadding, _affd);
                // Register in the x-direction only
                if (_Mode == RegisterX) {
                    for (i = 0; i < _affd->GetX(); i++) {
                        for (j = 0; j < _affd->GetY(); j++) {
                            for (k = 0; k < _affd->GetZ(); k++) {
                                _Status sx, sy, sz;
                                _affd->GetStatusCP(i, j, k, sx, sy, sz);
                                _affd->PutStatusCP(i, j, k, sx, _Passive, _Passive);
                            }
                        }
                    }
                }

                // Register in the y-direction only
                if (_Mode == RegisterY) {
                    for (i = 0; i < _affd->GetX(); i++) {
                        for (j = 0; j < _affd->GetY(); j++) {
                            for (k = 0; k < _affd->GetZ(); k++) {
                                _Status sx, sy, sz;
                                _affd->GetStatusCP(i, j, k, sx, sy, sz);
                                _affd->PutStatusCP(i, j, k, _Passive, sy, _Passive);
                            }
                        }
                    }
                }

                // Register in the x- and y-direction only
                if (_Mode == RegisterXY) {
                    for (i = 0; i < _affd->GetX(); i++) {
                        for (j = 0; j < _affd->GetY(); j++) {
                            for (k = 0; k < _affd->GetZ(); k++) {
                                _Status sx, sy, sz;
                                _affd->GetStatusCP(i, j, k, sx, sy, sz);
                                _affd->PutStatusCP(i, j, k, sx, sy, _Passive);
                            }
                        }
                    }
                }

                _affd = new irtkBSplineFreeFormTransformation3D(*_target, tdx, tdy, tdz);
                _mffd->PushLocalTransformation(_affd);
                // Padding of FFD
                irtkPadding(*_target, this->_TargetPadding, _affd);
                // Register in the x-direction only
                if (_Mode == RegisterX) {
                    for (i = 0; i < _affd->GetX(); i++) {
                        for (j = 0; j < _affd->GetY(); j++) {
                            for (k = 0; k < _affd->GetZ(); k++) {
                                _Status sx, sy, sz;
                                _affd->GetStatusCP(i, j, k, sx, sy, sz);
                                _affd->PutStatusCP(i, j, k, sx, _Passive, _Passive);
                            }
                        }
                    }
                }

                // Register in the y-direction only
                if (_Mode == RegisterY) {
                    for (i = 0; i < _affd->GetX(); i++) {
                        for (j = 0; j < _affd->GetY(); j++) {
                            for (k = 0; k < _affd->GetZ(); k++) {
                                _Status sx, sy, sz;
                                _affd->GetStatusCP(i, j, k, sx, sy, sz);
                                _affd->PutStatusCP(i, j, k, _Passive, sy, _Passive);
                            }
                        }
                    }
                }

                // Register in the x- and y-direction only
                if (_Mode == RegisterXY) {
                    for (i = 0; i < _affd->GetX(); i++) {
                        for (j = 0; j < _affd->GetY(); j++) {
                            for (k = 0; k < _affd->GetZ(); k++) {
                                _Status sx, sy, sz;
                                _affd->GetStatusCP(i, j, k, sx, sy, sz);
                                _affd->PutStatusCP(i, j, k, sx, sy, _Passive);
                            }
                        }
                    }
                }

                _NumberOfDofs += _affd->NumberOfDOFs();
                _NumberOfModels++;
            }

            dx /= 2; dy /= 2; 
            if(_target->GetZ() > 1) dz /= 2;
    }

    // save transformation back to mffd
    if(hasprevioustransformation)
        this->Image2Transformation();
}

void irtkSparseFreeFormRegistration::InitializeTransformation(int level){
    // Print debugging information
    this->Debug("irtkSparseFreeFormRegistration::InitializeTransformation(int)");

    // Intialize number of models
    double dx,dy,dz,tdx,tdy,tdz,odx,ody,odz;

    dx = tmp_target->GetXSize()*_LargestSpacing;
    dy = tmp_target->GetYSize()*_LargestSpacing;
    if(_target->GetZ() > 1)
        dz = tmp_target->GetZSize()*_LargestSpacing;
    else
        dz = 1;

    _NumberOfModels = 0;
    _NumberOfDofs = 0;
    odx = 0;
    ody = 0;
    odz = 0;
    while(dx > _target->GetXSize() && dy > _target->GetYSize()
        &&(dz > _target->GetZSize() || _target->GetZ() == 1)
        && _NumberOfModels < _mffd->NumberOfLevels()){

            if(dx > _target->GetXSize()*_target->GetX()/3.0){
                tdx = _target->GetXSize()*_target->GetX()/3.0;
            }else{
                tdx = dx;
            }

            if(dy > _target->GetYSize()*_target->GetY()/3.0){
                tdy = _target->GetYSize()*_target->GetY()/3.0;
            }else{
                tdy = dy;
            }

            if(_target->GetZ() > 1){
                if(dz > _target->GetZSize()*_target->GetZ()/3.0){
                    tdz = _target->GetZSize()*_target->GetZ()/3.0;
                }else{
                    tdz = dz;
                }
            }else{
                tdz = dz;
            }

            // check new transformation is different from previous
            if(tdx != odx || tdy != ody || tdz != odz){

                odx = tdx;
                ody = tdy;
                odz = tdz;

                _affd = (irtkBSplineFreeFormTransformation3D *)_mffd->GetLocalTransformation(_NumberOfModels);

                _NumberOfDofs += _affd->NumberOfDOFs();
                _NumberOfModels++;
            }

            dx /= 2; dy /= 2; 
            if(_target->GetZ() > 1) dz /= 2;
    }
}

void irtkSparseFreeFormRegistration::Initialize(int level)
{

    // Print debugging information
    this->Debug("irtkSparseFreeFormRegistration::Initialize(int)");

    // Initialize base class
    this->irtkImageRegistration2::Initialize(level);

    this->InitializeTransformation(level);

    cout << "Number of models used: " << _NumberOfModels << endl;

    // Allocate memory for _currentgradient vector
    _currentgradient = new double[_NumberOfDofs];
    _tmp = new double[_NumberOfDofs];
    _mask = new bool[_NumberOfDofs];
    // Initialize lambda3
    if(_Lambda3tmp > 0){
        this->Update(true);
        double norm;
        int count;

        this->irtkImageRegistration2::Evaluate();

        // Compute _currentgradient with respect to displacements
        this->irtkImageRegistration2::EvaluateGradient(&norm);

        this->ParametricGradient();

        if (this->_Lambda1 > 0){
            this->SmoothnessPenaltyGradient();
        }

        norm = 0;
        count = 0;

        for (int j = 0; j < _NumberOfModels; j++){
            _affd = (irtkBSplineFreeFormTransformation*)_mffd->GetLocalTransformation(j);
            // Move along __currentgradient direction
            for (int k = 0; k < _affd->NumberOfDOFs(); k++) {
                if(_affd->GetStatus(k) == Active){
                    count ++;
                }
            }
        }

        for(int i = 0; i < _NumberOfDofs; i++){
            norm += fabs(_currentgradient[i]);
        }

        if(norm > 0){
            cout << norm << " " << count << endl;
            norm = norm/count;
            _Lambda3 = norm*_Lambda3tmp;
            cout << "normalized sparsity penalty with respect to finite convergence property is:" << _Lambda3 << endl;
        }else{
            _Lambda3 = 0;
            cout << "current _currentgradient is 0!" << endl;
        }
    }

    _Lambda2 = _Lambda3;
}

void irtkSparseFreeFormRegistration::Finalize()
{
    // Print debugging information
    this->Debug("irtkSparseFreeFormRegistration::Finalize");

    // Finalize base class
    this->irtkImageRegistration2::Finalize();

    delete _mdffd;

    if (_g != NULL) delete []_g;
    _g = NULL;
    if (_h != NULL) delete []_h;
    _h = NULL;
}

void irtkSparseFreeFormRegistration::Finalize(int level)
{
    // Print debugging information
    this->Debug("irtkSparseFreeFormRegistration::Finalize(int)");

    // Finalize base class
    this->irtkImageRegistration2::Finalize(level);

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

}

double irtkSparseFreeFormRegistration::SmoothnessPenalty()
{
    int i;
    double penalty;
    penalty = 0;
    for(i = 0; i < _NumberOfModels; i++){
        _affd = (irtkBSplineFreeFormTransformation*)_mffd->GetLocalTransformation(i);
        if (_affd->GetZ() == 1) {
            penalty += _affd->Bending() / double(_affd->GetX()*_affd->GetY());
        } else {
            penalty += _affd->Bending() / double(_affd->GetX()*_affd->GetY()*_affd->GetZ());
        }
    }
    return penalty;
}

void irtkSparseFreeFormRegistration::SmoothnessPenaltyGradient()
{
    int i,t,index;
    double norm;

    index = 0;
    for(t = 0; t < _NumberOfModels; t++){
        _affd = (irtkBSplineFreeFormTransformation*)_mffd->GetLocalTransformation(t);

        // Compute normalization factor
        if (_affd->GetZ() == 1) {
            norm = _target->GetNumberOfVoxels() / double(_affd->GetX()*_affd->GetY());
        } else {
            norm = _target->GetNumberOfVoxels() / double(_affd->GetX()*_affd->GetY()*_affd->GetZ());
        }

        // Allocate memory
        double *tmp__currentgradient = new double[_affd->NumberOfDOFs()];

        // and initialize memory
        for (i = 0; i < _affd->NumberOfDOFs(); i++) {
            tmp__currentgradient[i] = 0.0;
        }

        // Compute _currentgradient of smoothness term
        _affd->BendingGradient(tmp__currentgradient);

        // Add _currentgradient to existing _currentgradient
        for (i = 0; i < _affd->NumberOfDOFs(); i++) {
            _currentgradient[index] += this->_Lambda1 * tmp__currentgradient[i] * norm;
            index++;
        }

        // Free memory
        delete []tmp__currentgradient;
    }
}

double irtkSparseFreeFormRegistration::SparsePenalty()
{
    int t,index;
    double sparsity,norm;

    norm = 0;

    sparsity = 0;
    for (t = 0; t < _NumberOfModels; t++){
        _affd = (irtkBSplineFreeFormTransformation *)_mffd->GetLocalTransformation(t);
        // approximate using a b spline model.
        for(index = 0; index < _affd->NumberOfDOFs(); index++){
            if(_affd->GetStatus(index) == _Active){
                sparsity += fabs(_affd->Get(index));
                norm ++;
            }
        }
    }
    if(sparsity > 0)
        return sparsity/norm;
    else
        return 0;
}

void irtkSparseFreeFormRegistration::SparsePenaltyGradient()
{
    int i,t,index;

    double sparsity;

    // Sparsity _currentgradient
    index = 0;
    for (t = 0; t < _NumberOfModels; t++){
        _affd = (irtkBSplineFreeFormTransformation *)_mffd->GetLocalTransformation(t);
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

double irtkSparseFreeFormRegistration::Evaluate()
{
    // Print debugging information
    this->Debug("irtkSparseFreeFormRegistration::Evaluate()");

    double tmp, similarity;

    // Evaluate similarity
    similarity = this->irtkImageRegistration2::Evaluate();

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

    if (this->_Lambda3 > 0) {
        tmp = this->_Lambda3*this->SparsePenalty();
        cout << " + Sparsity: " << tmp;
        similarity += tmp;
    }

    //Return similarity measure + penalty terms
    return similarity;
}

void irtkSparseFreeFormRegistration::InitializeIteration()
{
    // Print debugging information
    this->Debug("irtkSparseFreeFormRegistration::InitializeIteration()");

    // Update source image
    this->Update(true);

    if(this->SparsePenalty() > 0){
        _Lambda3 = _Lambda2;
    }else{
        _Lambda3 = 0;
    }
}

void irtkSparseFreeFormRegistration::FinalizeIteration()
{
    // Print debugging information
    this->Debug("irtkSparseFreeFormRegistration::FinalizeIteration()");

}

void irtkSparseFreeFormRegistration::ParametricGradient()
{
    // Print debugging information
    this->Debug("irtkSparseFreeFormRegistration::ParametricGradient");

    cout << "Evaluation of similarity gradient" << endl;

    /*if (_affd->GetZ() == 1) {
    this->EvaluateGradient2D();
    } else {
    this->EvaluateGradient3D();
    }*/

    // approximate start point
    int i,j,k,t,index,numberofpoints;
    double *x_pos,*y_pos,*z_pos,*x_dis,*y_dis,*z_dis;

    numberofpoints =_similarityGradient.GetNumberOfVoxels()/3;

    // generate displacement to be approximated
    x_pos = new double[numberofpoints];
    y_pos = new double[numberofpoints];
    z_pos = new double[numberofpoints];

    x_dis = new double[numberofpoints];
    y_dis = new double[numberofpoints];
    z_dis = new double[numberofpoints];

    index = 0;
    for(k = 0; k <_similarityGradient.GetZ(); k++){
        for(j = 0; j <_similarityGradient.GetY(); j++){
            for(i = 0; i <_similarityGradient.GetX(); i++){
                x_pos[index] = i;
                y_pos[index] = j;
                z_pos[index] = k;
                _similarityGradient.ImageToWorld(x_pos[index],y_pos[index],z_pos[index]);
                x_dis[index] =_similarityGradient.GetAsDouble(i, j, k, 0);
                y_dis[index] =_similarityGradient.GetAsDouble(i, j, k, 1);
                if(_similarityGradient.GetZ()>1){
                    z_dis[index] =_similarityGradient.GetAsDouble(i, j, k, 2);
                }
                else{
                    z_dis[index] = 0;
                }
                index ++;
            }
        }
    }

    int offset = 0;

    for (t = 0; t < _NumberOfModels; t++){
        _affd = (irtkBSplineFreeFormTransformation *)_mdffd->GetLocalTransformation(t);
        // approximate using a b spline model.
        _affd->ApproximateGradient(x_pos,y_pos,z_pos,x_dis,y_dis,z_dis,numberofpoints,&_currentgradient[offset]);
        offset += _affd->NumberOfDOFs();
    }

    cout << "Done" << endl;

    delete []x_pos;
    delete []y_pos;
    delete []z_pos;
    delete []x_dis;
    delete []y_dis;
    delete []z_dis;
}

void irtkSparseFreeFormRegistration::EvaluateGradient2D()
{
    double basis, pos[3];
    int i, j, m, i1, i2, j1, j2, k1, k2, x, y, index, index2, index3, globaloffset;

    globaloffset = 0;
    for (m = 0; m < _NumberOfModels; m++){
        _affd = (irtkBSplineFreeFormTransformation *)_mffd->GetLocalTransformation(m);
        // Initialize _currentgradient to zero
        for (i = 0; i < _affd->NumberOfDOFs(); i++) {
            _currentgradient[i+globaloffset] = 0;
        }

        // Loop over control points
        for (y = 0; y < _affd->GetY(); y++) {
            for (x = 0; x < _affd->GetX(); x++) {

                // Compute DoFs corresponding to the control point
                index  = _affd->LatticeToIndex(x, y, 0);
                index2 = index+_affd->GetX()*_affd->GetY()*_affd->GetZ();
                index3 = index+2*_affd->GetX()*_affd->GetY()*_affd->GetZ();

                // Check if any DoF corresponding to the control point is active
                if ((_affd->irtkTransformation::GetStatus(index) == _Active) || (_affd->irtkTransformation::GetStatus(index2) == _Active) || (_affd->irtkTransformation::GetStatus(index3) == _Active)) {

                    // If so, calculate bounding box of control point in image coordinates
                    _affd->BoundingBoxImage(_target, index, i1, j1, k1, i2, j2, k2, 1.0);

                    // Loop over all voxels in the target (reference) volume
                    for (j = j1; j <= j2; j++) {
                        for (i = i1; i <= i2; i++) {

                            // Check whether reference point is valid
                        	if ((_distanceMask(i, j, 0) == 0) && (_transformedSource(i, j, 0) > _SourcePadding)) {

                                // Convert position from voxel coordinates to world coordinates
                                pos[0] = i;
                                pos[1] = j;
                                pos[2] = 0;
                                _target->ImageToWorld(pos[0], pos[1], pos[2]);

                                // Convert world coordinates into lattice coordinates
                                _affd->WorldToLattice(pos[0], pos[1], pos[2]);

                                // Compute B-spline tensor product at pos
                                basis = _affd->B(pos[0] - x) * _affd->B(pos[1] - y);

                                // Convert voxel-based _currentgradient into _currentgradient with respect to parameters (chain rule)
                                //
                                // NOTE: This currently assumes that the control points displacements are aligned with the world coordinate displacements
                                //
                                _currentgradient[index + globaloffset]  += basis * _similarityGradient(i, j, 0, 0);
                                _currentgradient[index2 + globaloffset] += basis * _similarityGradient(i, j, 0, 1);
                                _currentgradient[index3 + globaloffset] += 0;
                            }

                        }
                    }
                }
            }
        }
        globaloffset += _affd->NumberOfDOFs();
    }
}

void irtkSparseFreeFormRegistration::EvaluateGradient3D()
{
    double basis, pos[3];
    int i, j, k, m, i1, i2, j1, j2, k1, k2, x, y, z, index, index2, index3, globaloffset;
    globaloffset = 0;
    for (m = 0; m < _NumberOfModels; m++){
        _affd = (irtkBSplineFreeFormTransformation *)_mffd->GetLocalTransformation(m);
        // Initialize _currentgradient to zero
        for (i = 0; i < _affd->NumberOfDOFs(); i++) {
            _currentgradient[i+globaloffset] = 0;
        }

        // Loop over control points
        for (z = 0; z < _affd->GetZ(); z++) {
            for (y = 0; y < _affd->GetY(); y++) {
                for (x = 0; x < _affd->GetX(); x++) {

                    // Compute DoFs corresponding to the control point
                    index  = _affd->LatticeToIndex(x, y, z);
                    index2 = index+_affd->GetX()*_affd->GetY()*_affd->GetZ();
                    index3 = index+2*_affd->GetX()*_affd->GetY()*_affd->GetZ();

                    // Check if any DoF corresponding to the control point is active
                    if ((_affd->irtkTransformation::GetStatus(index) == _Active) || (_affd->irtkTransformation::GetStatus(index2) == _Active) || (_affd->irtkTransformation::GetStatus(index3) == _Active)) {

                        // If so, calculate bounding box of control point in image coordinates
                        _affd->BoundingBoxImage(_target, index, i1, j1, k1, i2, j2, k2, 1);

                        // Loop over all voxels in the target (reference) volume
                        //
                        // NOTE: This currently assumes that the control point lattice is aligned with the target image
                        //
                        for (k = k1; k <= k2; k++) {
                            for (j = j1; j <= j2; j++) {
                                for (i = i1; i <= i2; i++) {
                                    // Check whether reference point is valid
                                	if ((_distanceMask(i, j, k) == 0) && (_transformedSource(i, j, k) > _SourcePadding)) {
                                        pos[0] = i;
                                        pos[1] = j;
                                        pos[2] = k;
                                        _target->ImageToWorld(pos[0], pos[1], pos[2]);

                                        // Convert world coordinates into lattice coordinates
                                        _affd->WorldToLattice(pos[0], pos[1], pos[2]);
                                        // Compute B-spline tensor product at current position
                                        basis = _affd->B(pos[0] - x) * _affd->B(pos[1] - y) * _affd->B(pos[2] - z);

                                        // Convert voxel-based _currentgradient into _currentgradient with respect to parameters (chain rule)
                                        //
                                        // NOTE: This currently assumes that the control points displacements are aligned with the world coordinate displacements
                                        //
                                        _currentgradient[index+globaloffset]  += basis * _similarityGradient(i, j, k, 0);
                                        _currentgradient[index2+globaloffset] += basis * _similarityGradient(i, j, k, 1);
                                        _currentgradient[index3+globaloffset] += basis * _similarityGradient(i, j, k, 2);
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

void irtkSparseFreeFormRegistration::Update(bool updateGradient)
{
    // Print debugging information
    this->Debug("irtkSparseFreeFormRegistration::Update()");

    // Update
    if (updateGradient == true) {
        this->UpdateSourceAndGradient();
    } else {
        this->UpdateSource();
    }
}

void irtkSparseFreeFormRegistration::UpdateSource()
{
    short *ptr1;
    double x, y, z, t1, t2, u1, u2, v1, v2;
    int a, b, c, i, j, k, offset1, offset2, offset3, offset4, offset5, offset6, offset7, offset8;

    IRTK_START_TIMING();

    // Generate transformed tmp image
    _transformedSource = *_target;

    // Calculate offsets for fast pixel access
    offset1 = 0;
    offset2 = 1;
    offset3 = this->_source->GetX();
    offset4 = this->_source->GetX()+1;
    offset5 = this->_source->GetX()*this->_source->GetY();
    offset6 = this->_source->GetX()*this->_source->GetY()+1;
    offset7 = this->_source->GetX()*this->_source->GetY()+this->_source->GetX();
    offset8 = this->_source->GetX()*this->_source->GetY()+this->_source->GetX()+1;

    if ((_target->GetZ() == 1) && (_source->GetZ() == 1)) {
        for (j = 0; j < _target->GetY(); j++) {
            for (i = 0; i < _target->GetX(); i++) {
            	if (_distanceMask(i, j, 0) == 0) {
                    x = i;
                    y = j;
                    z = 0;
                    _target->ImageToWorld(x, y, z);
                    _transformation->Transform(x, y, z);
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
                                _transformedSource(i, j, 0) = t1 * (u2 * ptr1[offset2] + u1 * ptr1[offset4]) + t2 * (u2 * ptr1[offset1] + u1 * ptr1[offset3]);
                            } else {
                                // Interpolation in source image
                                _transformedSource(i, j, 0) = _interpolator->Evaluate(x, y, 0);
                            }
                    } else {
                        _transformedSource(i, j, 0) = -1;
                    }
                } else {
                    _transformedSource(i, j, 0) = -1;
                }
            }
        }
    } else {
        for (k = 0; k < _target->GetZ(); k++) {
            for (j = 0; j < _target->GetY(); j++) {
                for (i = 0; i < _target->GetX(); i++) {
                	if (_distanceMask(i, j, k) == 0) {
                        x = i;
                        y = j;
                        z = k;
                        _target->ImageToWorld(x, y, z);
                        _transformation->Transform(x, y, z);
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
                                    _transformedSource(i, j, k) = (t1 * (u2 * (v2 * ptr1[offset2] + v1 * ptr1[offset6]) +
                                        u1 * (v2 * ptr1[offset4] + v1 * ptr1[offset8])) +
                                        t2 * (u2 * (v2 * ptr1[offset1] + v1 * ptr1[offset5]) +
                                        u1 * (v2 * ptr1[offset3] + v1 * ptr1[offset7])));
                                } else {
                                    // Interpolation in source image
                                    _transformedSource(i, j, k) = _interpolator->Evaluate(x, y, z);
                                }
                        } else {
                            _transformedSource(i, j, k) = -1;
                        }
                    } else {
                        _transformedSource(i, j, k) = -1;
                    }
                }
            }
        }
    }

    IRTK_END_TIMING("irtkImageRegistration2::UpdateSource");
}

void irtkSparseFreeFormRegistration::UpdateSourceAndGradient()
{
    short *ptr1;
    irtkRealPixel *ptr2;
    double x, y, z, t1, t2, u1, u2, v1, v2;
    int a, b, c, i, j, k, offset1, offset2, offset3, offset4, offset5, offset6, offset7, offset8;

    IRTK_START_TIMING();

    // Generate transformed tmp image
    _transformedSource = *_target;

    // Calculate offsets for fast pixel access
    offset1 = 0;
    offset2 = 1;
    offset3 = this->_source->GetX();
    offset4 = this->_source->GetX()+1;
    offset5 = this->_source->GetX()*this->_source->GetY();
    offset6 = this->_source->GetX()*this->_source->GetY()+1;
    offset7 = this->_source->GetX()*this->_source->GetY()+this->_source->GetX();
    offset8 = this->_source->GetX()*this->_source->GetY()+this->_source->GetX()+1;

    if ((_target->GetZ() == 1) && (_source->GetZ() == 1)) {
        for (j = 0; j < _target->GetY(); j++) {
            for (i = 0; i < _target->GetX(); i++) {
            	if (_distanceMask(i, j, 0) == 0) {
                    x = i;
                    y = j;
                    z = 0;
                    _target->ImageToWorld(x, y, z);
                    _transformation->Transform(x, y, z);
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
                                _transformedSource(i, j, 0) = t1 * (u2 * ptr1[offset2] + u1 * ptr1[offset4]) + t2 * (u2 * ptr1[offset1] + u1 * ptr1[offset3]);

                                // Linear interpolation in _currentgradient image
                                ptr2 = _sourceGradient.GetPointerToVoxels(a, b, 0, 0);
                                _transformedSourceGradient(i, j, 0, 0) = t1 * (u2 * ptr2[offset2] + u1 * ptr2[offset4]) + t2 * (u2 * ptr2[offset1] + u1 * ptr2[offset3]);
                                ptr2 = _sourceGradient.GetPointerToVoxels(a, b, 0, 1);
                                _transformedSourceGradient(i, j, 0, 1) = t1 * (u2 * ptr2[offset2] + u1 * ptr2[offset4]) + t2 * (u2 * ptr2[offset1] + u1 * ptr2[offset3]);

                            } else {
                                // Interpolation in source image
                                _transformedSource(i, j, 0) = _interpolator->Evaluate(x, y, 0);

                                // Interpolation in _currentgradient image
                                _transformedSourceGradient(i, j, 0, 0) = _interpolatorGradient->Evaluate(x, y, 0, 0);
                                _transformedSourceGradient(i, j, 0, 1) = _interpolatorGradient->Evaluate(x, y, 0, 1);
                            }
                    } else {
                        _transformedSource(i, j, 0) = -1;
                        _transformedSourceGradient(i, j, 0, 0) = 0;
                        _transformedSourceGradient(i, j, 0, 1) = 0;
                    }
                } else {
                    _transformedSource(i, j, 0) = -1;
                    _transformedSourceGradient(i, j, 0, 0) = 0;
                    _transformedSourceGradient(i, j, 0, 1) = 0;

                }
            }
        }
    } else {
        for (k = 0; k < _target->GetZ(); k++) {
            for (j = 0; j < _target->GetY(); j++) {
                for (i = 0; i < _target->GetX(); i++) {
                	if (_distanceMask(i, j, k) == 0) {
                        x = i;
                        y = j;
                        z = k;
                        _target->ImageToWorld(x, y, z);
                        _transformation->Transform(x, y, z);
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
                                    _transformedSource(i, j, k) = (t1 * (u2 * (v2 * ptr1[offset2] + v1 * ptr1[offset6]) +
                                        u1 * (v2 * ptr1[offset4] + v1 * ptr1[offset8])) +
                                        t2 * (u2 * (v2 * ptr1[offset1] + v1 * ptr1[offset5]) +
                                        u1 * (v2 * ptr1[offset3] + v1 * ptr1[offset7])));

                                    // Linear interpolation in _currentgradient image
                                    ptr2 = _sourceGradient.GetPointerToVoxels(a, b, c, 0);
                                    _transformedSourceGradient(i, j, k, 0) = (t1 * (u2 * (v2 * ptr2[offset2] + v1 * ptr2[offset6]) +
                                        u1 * (v2 * ptr2[offset4] + v1 * ptr2[offset8])) +
                                        t2 * (u2 * (v2 * ptr2[offset1] + v1 * ptr2[offset5]) +
                                        u1 * (v2 * ptr2[offset3] + v1 * ptr2[offset7])));
                                    ptr2 = _sourceGradient.GetPointerToVoxels(a, b, c, 1);
                                    _transformedSourceGradient(i, j, k, 1) = (t1 * (u2 * (v2 * ptr2[offset2] + v1 * ptr2[offset6]) +
                                        u1 * (v2 * ptr2[offset4] + v1 * ptr2[offset8])) +
                                        t2 * (u2 * (v2 * ptr2[offset1] + v1 * ptr2[offset5]) +
                                        u1 * (v2 * ptr2[offset3] + v1 * ptr2[offset7])));
                                    ptr2 = _sourceGradient.GetPointerToVoxels(a, b, c, 2);
                                    _transformedSourceGradient(i, j, k, 2) = (t1 * (u2 * (v2 * ptr2[offset2] + v1 * ptr2[offset6]) +
                                        u1 * (v2 * ptr2[offset4] + v1 * ptr2[offset8])) +
                                        t2 * (u2 * (v2 * ptr2[offset1] + v1 * ptr2[offset5]) +
                                        u1 * (v2 * ptr2[offset3] + v1 * ptr2[offset7])));

                                } else {
                                    // Interpolation in source image
                                    _transformedSource(i, j, k) = _interpolator->Evaluate(x, y, z);

                                    // Interpolation in _currentgradient image
                                    _transformedSourceGradient(i, j, k, 0) = _interpolatorGradient->Evaluate(x, y, z, 0);
                                    _transformedSourceGradient(i, j, k, 1) = _interpolatorGradient->Evaluate(x, y, z, 1);
                                    _transformedSourceGradient(i, j, k, 2) = _interpolatorGradient->Evaluate(x, y, z, 2);
                                }
                        } else {
                            _transformedSource(i, j, k) = -1;
                            _transformedSourceGradient(i, j, k, 0) = 0;
                            _transformedSourceGradient(i, j, k, 1) = 0;
                            _transformedSourceGradient(i, j, k, 2) = 0;
                        }
                    } else {
                        _transformedSource(i, j, k) = -1;
                        _transformedSourceGradient(i, j, k, 0) = 0;
                        _transformedSourceGradient(i, j, k, 1) = 0;
                        _transformedSourceGradient(i, j, k, 2) = 0;
                    }
                }
            }
        }
    }

    IRTK_END_TIMING("irtkImageRegistration2::UpdateSourceAndGradient");
}

void irtkSparseFreeFormRegistration::NormalizeGradient()
{
    int x,y,z,t,index1,index2,index3,offset,globaloffset;
    double norm,spacingnorm;

    globaloffset = 0;

    for (t = 0; t < _NumberOfModels; t++){
        _affd = (irtkBSplineFreeFormTransformation *)_mffd->GetLocalTransformation(t);
        offset = _affd->GetX()*_affd->GetY()*_affd->GetZ();

        // Compute normalization factor
        if (_affd->GetZ() == 1) {
            spacingnorm = (double(_target->GetNumberOfVoxels()) / double(_affd->GetX()*_affd->GetY()));
        } else {
            spacingnorm = (double(_target->GetNumberOfVoxels()) / double(_affd->GetX()*_affd->GetY()*_affd->GetZ()));
        }

        // approximate using a b spline model.
        for (z = 0; z < _affd->GetZ(); z++) {
            for (y = 0; y < _affd->GetY(); y++) {
                for (x = 0; x < _affd->GetX(); x++) {
                    index1 = _affd->LatticeToIndex(x,y,z) + globaloffset;
                    index2 = index1 + offset;
                    index3 = index2 + offset;
                    norm = pow(_currentgradient[index1],2.0)
                        + pow(_currentgradient[index2],2.0)
                        + pow(_currentgradient[index3],2.0);

                    //normalize
                    if(norm > 0){
                        norm = sqrt(norm);
                        //_currentgradient[index1] = _currentgradient[index1]/(spacingnorm);
                        //_currentgradient[index2] = _currentgradient[index2]/(spacingnorm);
                        //_currentgradient[index3] = _currentgradient[index3]/(spacingnorm);
                        _currentgradient[index1] = _currentgradient[index1]/(norm + spacingnorm*_Epsilon);
                        _currentgradient[index2] = _currentgradient[index2]/(norm + spacingnorm*_Epsilon);
                        _currentgradient[index3] = _currentgradient[index3]/(norm + spacingnorm*_Epsilon);
                    }
                }
            }
        }

        globaloffset += _affd->NumberOfDOFs();
    }


}

double irtkSparseFreeFormRegistration::EvaluateGradient(double *gradient)
{
    double norm, max_length, gg, dgg, gamma;
    int i, x, y, z, n, m, index, index2, index3;

    // Compute _currentgradient with respect to displacements
    this->irtkImageRegistration2::EvaluateGradient(gradient);

    if(_DebugFlag){
        _target->Write("target.nii.gz");
        _transformedSource.Write("source.nii.gz");
        _transformedSourceGradient.Write("source_currentgradient.nii.gz");
    }

    this->ParametricGradient();

    if (this->_Lambda1 > 0){
        this->SmoothnessPenaltyGradient();
    }

    if (this->_Lambda3 > 0){
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
    m = 0;
    for(n = 0; n < _NumberOfModels; n++){
        _affd = (irtkBSplineFreeFormTransformation*)_mffd->GetLocalTransformation(n);
        for (z = 0; z < _affd->GetZ(); z++) {
            for (y = 0; y < _affd->GetY(); y++) {
                for (x = 0; x < _affd->GetX(); x++) {
                    index  = m+_affd->LatticeToIndex(x, y, z);
                    index2 = index+_affd->GetX()*_affd->GetY()*_affd->GetZ();
                    index3 = index+2*_affd->GetX()*_affd->GetY()*_affd->GetZ();
                    norm = sqrt(_currentgradient[index] * _currentgradient[index] 
                    + _currentgradient[index2] * _currentgradient[index2] 
                    + _currentgradient[index3] * _currentgradient[index3]);
                    if (norm > max_length) 
                        max_length = norm;
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

void irtkSparseFreeFormRegistration::Run()
{
    int i, j, k, index;
    char buffer[256];
    double gradient, delta, step, min_step, max_step, max_length, best_similarity, new_similarity, old_similarity,norm;

    // Print debugging information
    this->Debug("irtkSparseFreeFormRegistration::Run");

    if (_source == NULL) {
        cerr << "irtkSparseFreeFormRegistration::Run: Filter has no source input" << endl;
        exit(1);
    }

    if (_target == NULL) {
        cerr << "irtkSparseFreeFormRegistration::Run: Filter has no target input" << endl;
        exit(1);
    }

    if (_transformation == NULL) {
        cerr << "irtkSparseFreeFormRegistration::Run: Filter has no transformation output" << endl;
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

        // Save pre-processed images if we are debugging
        sprintf(buffer, "source_%d.nii.gz", _CurrentLevel);
        if (_DebugFlag == true) _source->Write(buffer);
        sprintf(buffer, "target_%d.nii.gz", _CurrentLevel);
        if (_DebugFlag == true) _target->Write(buffer);

        _CurrentIteration = 0;
        while (_CurrentIteration < _NumberOfIterations[_CurrentLevel]) {
            cout << "Iteration = " << _CurrentIteration + 1 << " (out of " << _NumberOfIterations[_CurrentLevel] << ")"<< endl;

            this->InitializeIteration();

            // Compute current metric value
            best_similarity = old_similarity = this->Evaluate();
            cout << " Current best metric value is " << best_similarity << endl;

            _Epsilon = fabs(best_similarity) * 0.0001;

            // Compute _currentgradient of similarity metric. The function EvaluateGradient() returns the maximum control point length in the _currentgradient
            max_length = this->EvaluateGradient(&gradient);

            cout.flush();

            // Step along _currentgradient direction until no further improvement is necessary
            i = 0;
            delta = 0;

            if(max_length > 0){

                if(_Lambda3 > 0){
                    // Define steps
                    norm = 0;
                    int count = 0;
                    for (j = 0; j < _NumberOfModels; j++){
                        _affd = (irtkBSplineFreeFormTransformation*)_mffd->GetLocalTransformation(j);
                        // average of dof
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
                step = max_step;
                do {
                    double current = step/max_length;

                    index = 0;
                    for (j = 0; j < _NumberOfModels; j++){
                        _affd = (irtkBSplineFreeFormTransformation*)_mffd->GetLocalTransformation(j);
                        // Move along gradient direction
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

                    // Compute new similarity

                    cout << "Current metric value is ";

                    new_similarity = this->Evaluate();

                    if (new_similarity < best_similarity - _Epsilon) {
                        cout << " = " << new_similarity << " accepted; step = " << step << endl;
                        best_similarity = new_similarity;
                        delta += step;

                        step = step * 1.1;
                    } else {
                        // Last step was no improvement, so back track
                        cout << " = " << new_similarity << " rejected; step = " << step << endl;
                        index = 0;
                        for (j = 0; j < _NumberOfModels; j++){
                            _affd = (irtkBSplineFreeFormTransformation*)_mffd->GetLocalTransformation(j);
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
            }

            // Check for convergence
            if (delta == 0) {
                break;
            }

            this->FinalizeIteration();
        }

        // Do the final cleaning up for this level
        this->Finalize(_CurrentLevel);
    }

    // Do the final cleaning up for all levels
    this->Finalize();
}

bool irtkSparseFreeFormRegistration::Read(char *buffer1, char *buffer2, int &level)
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
        return this->irtkImageRegistration2::Read(buffer1, buffer2, level);
    } else {
        return ok;
    }
}

void irtkSparseFreeFormRegistration::Write(ostream &to)
{
    to << "\n#\n# Sparse non-rigid registration parameters\n#\n\n";
    to << "Bending penalty                     = " << this->_Lambda1 << endl;
    to << "Sparsity constrain                  = " << this->_Lambda3 << endl;
    to << "Coarsest spacing                    = " << this->_LargestSpacing << endl;
    to << "Finest spacing                      = " << this->_FinestSpacing << endl;

    this->irtkImageRegistration2::Write(to);
}
