/*=========================================================================

Library   : Image Registration Toolkit (IRTK)
Module    : $Id$
Copyright : Imperial College, Department of Computing
Visual Information Processing (VIP), 2009 onwards
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

extern irtkGreyImage **tmp_tsource, *tmp_target;

irtkImageTFFDRegistration::irtkImageTFFDRegistration()
{
    // Print debugging information
    this->Debug("irtkImageTFFDRegistration::irtkImageTFFDRegistration");

    // Default optimization
    _OptimizationMethod = GradientDescent;

    // Default parameters for non-rigid registration
    _Lambda1     = 0;
    _Lambda2     = 0;
    _Lambda3     = 0;
    _DX          = 20;
    _DY          = 20;
    _DZ          = 20;
    _DT		   = 0.05;
    _Subdivision = true;
    _Mode        = RegisterXYZ;
    _MFFDMode    = false;
    _adjugate    = NULL;
    _determinant = NULL;
    _periodic    = true;
}

void irtkImageTFFDRegistration::GuessParameter()
{
    int i, slices = false;
    double xsize, ysize, zsize, spacing, minsize;

    if ((_target == NULL) || (_source == NULL)) {
        cerr << "irtkImageTFFDRegistration::GuessParameter: Target and source image not found" << endl;
        exit(1);
    }

    // Default parameters for registration
    _NumberOfLevels     = 3;
    _NumberOfBins       = 64;

    // Default parameters for optimization
    _SimilarityMeasure  = SSD;
    _Epsilon            = 0.0001;

    _periodic = true;

    // Read target pixel size
    _target->GetPixelSize(&xsize, &ysize, &zsize);
    minsize = (xsize <= ysize)   ? xsize : ysize;
    minsize = (zsize <= minsize) ? zsize : minsize;

    // Use xsize as spacing
    spacing = xsize;

	if (_target->GetZ()==1) {
		slices = true;
	}
    
    // Default target parameters
    _TargetResolution[0][0] = xsize;
    _TargetResolution[0][1] = ysize;
    _TargetResolution[0][2] = zsize;

    for (i = 1; i < _NumberOfLevels; i++) {
        _TargetResolution[i][0] = _TargetResolution[i-1][0] * 2;
        _TargetResolution[i][1] = _TargetResolution[i-1][1] * 2;
        if (slices)
            _TargetResolution[i][2] = zsize;
        else 
            _TargetResolution[i][2] = _TargetResolution[i-1][2] * 2;
    }

    // Read source pixel size
	slices = false;
    _source[0]->GetPixelSize(&xsize, &ysize, &zsize);
    minsize = (xsize <= ysize)   ? xsize : ysize;
    minsize = (zsize <= minsize) ? zsize : minsize;

	for(int i = 0; i < _N_source; i++){
		if (_source[i]->GetZ()==1) {
			slices = true;
			break;
		}
	}

    // Default source parameters
    _SourceResolution[0][0] = xsize;
    _SourceResolution[0][1] = ysize;
    _SourceResolution[0][2] = zsize;

    for (i = 1; i < _NumberOfLevels; i++) {
        _SourceResolution[i][0] = _SourceResolution[i-1][0] * 2;
        _SourceResolution[i][1] = _SourceResolution[i-1][1] * 2;
        if (slices)
            _SourceResolution[i][2] = zsize;
        else 
            _SourceResolution[i][2] = _SourceResolution[i-1][2] * 2;
    }

    // Default parameters for non-rigid registration
    _Lambda1            = 0; //recommended value 0.0001
    _Lambda2            = 0; //recommended value 1
    _Lambda3            = 0;
    _DX = xsize * pow(2.0,_NumberOfLevels-1);
    _DY = ysize * pow(2.0,_NumberOfLevels-1);
    if (slices) {
        _DZ               = 1;
    } else {
        _DZ               = zsize * pow(2.0,_NumberOfLevels-1);
    }
    _DT				  = 0.05;
    _Subdivision        = true;

    // Remaining parameters
    for (i = 0; i < _NumberOfLevels; i++) {
        _NumberOfIterations[i] = 40;
        _MinStep[i]            = 0.01;
        _MaxStep[i]            = 1.0;
    }

    _TargetPadding = MIN_GREY;
}

void irtkImageTFFDRegistration::Initialize()
{
    // Print debugging information
    this->Debug("irtkImageTFFDRegistration::Initialize");

    // Initialize base class
    this->irtkTemporalImageRegistration::Initialize();

    // Pointer to multi-level FFD
    _mffd = (irtkMultiLevelFreeFormTransformation *)_transformation;

    // Create FFD
    if (_MFFDMode == false) {
        if (_mffd->NumberOfLevels() == 0) {
            _affd = new irtkBSplineFreeFormTransformationPeriodic(*_target, this->_DX, this->_DY, this->_DZ, this->_DT);
            _affd->PeriodicOn();
        } else {
            _affd = (irtkBSplineFreeFormTransformationPeriodic *)_mffd->PopLocalTransformation();
            _affd->PeriodicOn();
        }
    } else {
        _affd = new irtkBSplineFreeFormTransformationPeriodic(*_target, this->_DX, this->_DY, this->_DZ, this->_DT);
        _affd->PeriodicOn();
    }

    _adjugate = new irtkMatrix[_affd->NumberOfDOFs()/3];
    _determinant = new double[_affd->NumberOfDOFs()/3];
}

void irtkImageTFFDRegistration::Initialize(int level)
{
    int i, j, k, l, n;
    double x, y, z, t, *ptr2latt, *ptr2disp;

    // Print debugging information
    this->Debug("irtkImageTFFDRegistration::Initialize(int)");

    // Initialize base class
    this->irtkTemporalImageRegistration::Initialize(level);

    // Padding of FFD
    //  for (int n=0; n<_N_source; n++)
    //    irtkPadding(*_target[n], this->_TargetPadding, _affd);

    // Register in the x-direction only
    if (_Mode == RegisterX) {
        for (i = 0; i < _affd->GetX(); i++) {
            for (j = 0; j < _affd->GetY(); j++) {
                for (k = 0; k < _affd->GetZ(); k++) {
                    for (l = 0; l < _affd->GetT(); l++) {
                        _Status sx, sy, sz;
                        _affd->GetStatusCP(i, j, k, l, sx, sy, sz);
                        _affd->PutStatusCP(i, j, k, l, sx, _Passive, _Passive);
                    }
                }
            }
        }
    }

    // Register in the y-direction only
    if (_Mode == RegisterY) {
        for (i = 0; i < _affd->GetX(); i++) {
            for (j = 0; j < _affd->GetY(); j++) {
                for (k = 0; k < _affd->GetZ(); k++) {
                    for (l = 0; l < _affd->GetT(); l++) {
                        _Status sx, sy, sz;
                        _affd->GetStatusCP(i, j, k, l, sx, sy, sz);
                        _affd->PutStatusCP(i, j, k, l, _Passive, sy, _Passive);
                    }
                }
            }
        }
    }

    // Register in the x- and y-direction only
    if (_Mode == RegisterXY) {
        for (i = 0; i < _affd->GetX(); i++) {
            for (j = 0; j < _affd->GetY(); j++) {
                for (k = 0; k < _affd->GetZ(); k++) {
                    for (l = 0; l < _affd->GetT(); l++) {
                        _Status sx, sy, sz;
                        _affd->GetStatusCP(i, j, k, l, sx, sy, sz);
                        _affd->PutStatusCP(i, j, k, l, sx, sy, _Passive);
                    }
                }
            }
        }
    }

    // Allocate memory for global displacements
    _displacementLUT = new double*[_N_source];

    // Allocate memory for lattice coordinates
    _latticeCoordLUT = new double*[_N_source];

    for (n = 0; n < _N_source; n++) {
        _displacementLUT[n] = new double[_target->GetNumberOfVoxels() * 3];
        _latticeCoordLUT[n] = new double[_target->GetNumberOfVoxels() * 3];

        ptr2disp = _displacementLUT[n];
        ptr2latt = _latticeCoordLUT[n];

        for (k = 0; k < _target->GetZ(); k++) {
            for (j = 0; j < _target->GetY(); j++) {
                for (i = 0; i < _target->GetX(); i++) {
                    if (_target->Get(i, j, k) >= 0) {
                        x = i;
                        y = j;
                        z = k;
                        _target->ImageToWorld(x, y, z);
                        t = _t_real[n];
                        _mffd->Transform(x, y, z, t);
                        ptr2disp[0] = x;
                        ptr2disp[1] = y;
                        ptr2disp[2] = z;
                        x = i;
                        y = j;
                        z = k;
                        _target->ImageToWorld(x, y, z);
                        _affd->WorldToLattice(x, y, z);
                        ptr2latt[0] = x;
                        ptr2latt[1] = y;
                        ptr2latt[2] = z;
                    } else {
                        ptr2disp[0] = 0;
                        ptr2disp[1] = 0;
                        ptr2disp[2] = 0;
                        ptr2latt[0] = 0;
                        ptr2latt[1] = 0;
                        ptr2latt[2] = 0;
                    }
                    ptr2disp += 3;
                    ptr2latt += 3;
                }
            }
        }
    }
}

void irtkImageTFFDRegistration::Finalize()
{
    // Push local transformation back on transformation stack
    _mffd->PushLocalTransformation(_affd);

    // Print debugging information
    this->Debug("irtkImageTFFDRegistration::Finalize");

    // Finalize base class
    this->irtkTemporalImageRegistration::Finalize();
    delete []_adjugate;
    _adjugate = NULL;
    delete []_determinant;
    _determinant = NULL;
}

void irtkImageTFFDRegistration::Finalize(int level)
{
    // Print debugging information
    this->Debug("irtkImageTFFDRegistration::Finalize(int)");

    // Finalize base class
    this->irtkTemporalImageRegistration::Finalize(level);

    // Check if we are not at the lowest level of resolution
    if (level != 0) {
        if (this->_Subdivision == true) {
            //cout << "irtkImageTFFDRegistration::Finalize: Subdivide start" << endl;
            switch (_SubdivisionDim[level])
            {
            case 4:
                _affd->Subdivide4D();
                break;
            case 3:
                _affd->Subdivide3D();
                break;
            case 2:
                _affd->Subdivide2D();
                break;
            case 0:
                cout << "no subdivision!!!" << endl;
                break;
            default:
                cerr<<"wrong subdivision mode/dimension for level "<<level<<endl;
                exit(1);
            }
            //cout << "irtkImageTFFDRegistration::Finalize: Subdivide start" << endl;
        } else {

            // Push local transformation back on transformation stack
            _mffd->PushLocalTransformation(_affd);

            // Create new FFD
            _affd = new irtkBSplineFreeFormTransformationPeriodic(*_target,
                this->_DX / pow(2.0, this->_NumberOfLevels-level),
                this->_DY / pow(2.0, this->_NumberOfLevels-level),
                this->_DZ / pow(2.0, this->_NumberOfLevels-level),
                this->_DT);
            _affd->PeriodicOn();

            // Push local transformation back on transformation stack
            _mffd->PushLocalTransformation(_affd);
        }
    }
    delete []_adjugate;
    _adjugate = NULL;
    _adjugate = new irtkMatrix[_affd->NumberOfDOFs()/3];
    delete []_determinant;
    _determinant = NULL;
    _determinant = new double[_affd->NumberOfDOFs()/3];
    delete []_displacementLUT;
    _displacementLUT = NULL;
    delete []_latticeCoordLUT;
    _latticeCoordLUT = NULL;
}

void irtkImageTFFDRegistration::UpdateSource()
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
    for (n = 0; n < _N_source; n++) {
        _transformedSource[n] = *_target;
    }

    //cout<<"irtkImageTFFDRegistration::UpdateSource start"<<endl;
    for (n = 0; n < _N_source; n++) {
        double *ptr2disp = _displacementLUT[n];
        double *ptr2latt = _latticeCoordLUT[n];
		// Calculate offsets for fast pixel access
		offset1 = 0;
		offset2 = 1;
		offset3 = this->_source[n]->GetX();
		offset4 = this->_source[n]->GetX()+1;
		offset5 = this->_source[n]->GetX()*this->_source[n]->GetY();
		offset6 = this->_source[n]->GetX()*this->_source[n]->GetY()+1;
		offset7 = this->_source[n]->GetX()*this->_source[n]->GetY()+this->_source[n]->GetX();
		offset8 = this->_source[n]->GetX()*this->_source[n]->GetY()+this->_source[n]->GetX()+1;
        if ((_target->GetZ() == 1) && (_source[n]->GetZ() == 1)) {
            for (j = 0; j < _target->GetY(); j++) {
                for (i = 0; i < _target->GetX(); i++) {
                    if (_target->Get(i, j, 0) >= 0) {
                        x = ptr2latt[0];
                        y = ptr2latt[1];
                        z = ptr2latt[2];
                        t = _affd->TimeToLattice(_t_real[n]);
                        _affd->FFD1(x, y, z, t);
                        z  = 0;
                        x += ptr2disp[0];
                        y += ptr2disp[1];
                        z += ptr2disp[2];
                        _source[n]->WorldToImage(x, y, z);

                        // Check whether transformed point is inside volume
                        if ((x > 0) && (x < _source[n]->GetX()-1) &&
                            (y > 0) && (y < _source[n]->GetY()-1)) {

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
                                    ptr1 = (short *)_source[n]->GetScalarPointer(a, b, 0);
                                    _transformedSource[n](i, j, 0) = t1 * (u2 * ptr1[offset2] + u1 * ptr1[offset4]) + t2 * (u2 * ptr1[offset1] + u1 * ptr1[offset3]);
                                } else {
                                    // Interpolation in source image
                                    _transformedSource[n](i, j, 0) = _interpolator[n]->Evaluate(x, y, 0);
                                }
                        } else {
                            _transformedSource[n](i, j, 0) = -1;
                        }
                    } else {
                        _transformedSource[n](i, j, 0) = -1;
                    }
                    ptr2disp += 3;
                    ptr2latt += 3;
                }
            }
        } else {
            for (k = 0; k < _target->GetZ(); k++) {
                for (j = 0; j < _target->GetY(); j++) {
                    for (i = 0; i < _target->GetX(); i++) {
                        if (_target->Get(i, j, k) >= 0) {
                            x = ptr2latt[0];
                            y = ptr2latt[1];
                            z = ptr2latt[2];
                            t = _affd->TimeToLattice(_t_real[n]);
                            _affd->FFD1(x, y, z, t);
                            x += ptr2disp[0];
                            y += ptr2disp[1];
                            z += ptr2disp[2];
                            _source[n]->WorldToImage(x, y, z);

                            // Check whether transformed point is inside volume
                            if ((x > 0) && (x < _source[n]->GetX()-1) &&
                                (y > 0) && (y < _source[n]->GetY()-1) &&
                                (z > 0) && (z < _source[n]->GetZ()-1)) {

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
                                        ptr1 = (short *)_source[n]->GetScalarPointer(a, b, c);
                                        _transformedSource[n](i, j, k) = (t1 * (u2 * (v2 * ptr1[offset2] + v1 * ptr1[offset6]) +
                                            u1 * (v2 * ptr1[offset4] + v1 * ptr1[offset8])) +
                                            t2 * (u2 * (v2 * ptr1[offset1] + v1 * ptr1[offset5]) +
                                            u1 * (v2 * ptr1[offset3] + v1 * ptr1[offset7])));
                                    } else {
                                        // Interpolation in source image
                                        _transformedSource[n](i, j, k) = _interpolator[n]->Evaluate(x, y, z);
                                    }
                            } else {
                                _transformedSource[n](i, j, k) = -1;
                            }
                        } else {
                            _transformedSource[n](i, j, k) = -1;
                        }
                        ptr2disp += 3;
                        ptr2latt += 3;
                    }
                }
            }
        }
    }
    //cout<<"irtkImageTFFDRegistration::UpdateSource end"<<endl;
#ifdef USE_TIMING
    // Stop timing
    end = clock();
    cpu_time_used = ((double) (end - start)) / CLOCKS_PER_SEC;
    cout << "CPU time for irtkImageTFFDRegistration::UpdateSource() = " << cpu_time_used << endl;
#endif

}

void irtkImageTFFDRegistration::UpdateSourceAndGradient()
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
    for (n = 0; n < _N_source; n++) {
        _transformedSource[n] = *_source[n];
    }

    //cout<<"irtkImageTFFDRegistration::UpdateSourceAndGradient start"<<endl;
    for (n = 0; n < _N_source; n++) {
        double *ptr2disp = _displacementLUT[n];
        double *ptr2latt = _latticeCoordLUT[n];
		// Calculate offsets for fast pixel access
		offset1 = 0;
		offset2 = 1;
		offset3 = this->_source[n]->GetX();
		offset4 = this->_source[n]->GetX()+1;
		offset5 = this->_source[n]->GetX()*this->_source[n]->GetY();
		offset6 = this->_source[n]->GetX()*this->_source[n]->GetY()+1;
		offset7 = this->_source[n]->GetX()*this->_source[n]->GetY()+this->_source[n]->GetX();
		offset8 = this->_source[n]->GetX()*this->_source[n]->GetY()+this->_source[n]->GetX()+1;
        if ((_target->GetZ() == 1) && (_source[n]->GetZ() == 1)) {
            for (j = 0; j < _target->GetY(); j++) {
                for (i = 0; i < _target->GetX(); i++) {
                    if (_target->Get(i, j, 0) >= 0) {
                        x = ptr2latt[0];
                        y = ptr2latt[1];
                        z = ptr2latt[2];
                        t = _affd->TimeToLattice(_t_real[n]);
                        _affd->FFD1(x, y, z, t);
                        z  = 0;
                        x += ptr2disp[0];
                        y += ptr2disp[1];
                        z += ptr2disp[2];
                        _source[n]->WorldToImage(x, y, z);

                        // Check whether transformed point is inside volume
                        if ((x > 0) && (x < _source[n]->GetX()-1) &&
                            (y > 0) && (y < _source[n]->GetY()-1)) {

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
                                    ptr1 = (short *)_source[n]->GetScalarPointer(a, b, 0);
                                    _transformedSource[n](i, j, 0) = t1 * (u2 * ptr1[offset2] + u1 * ptr1[offset4]) + t2 * (u2 * ptr1[offset1] + u1 * ptr1[offset3]);

                                    // Linear interpolation in gradient image
                                    ptr2 = _sourceGradient[n].GetPointerToVoxels(a, b, 0, 0);
                                    _transformedSourceGradient[n](i, j, 0, 0) = t1 * (u2 * ptr2[offset2] + u1 * ptr2[offset4]) + t2 * (u2 * ptr2[offset1] + u1 * ptr2[offset3]);
                                    ptr2 = _sourceGradient[n].GetPointerToVoxels(a, b, 0, 1);
                                    _transformedSourceGradient[n](i, j, 0, 1) = t1 * (u2 * ptr2[offset2] + u1 * ptr2[offset4]) + t2 * (u2 * ptr2[offset1] + u1 * ptr2[offset3]);

                                } else {
                                    // Interpolation in source image
                                    _transformedSource[n](i, j, 0) = _interpolator[n]->Evaluate(x, y, 0);

                                    // Interpolation in gradient image
                                    _transformedSourceGradient[n](i, j, 0, 0) = _interpolatorGradient[n]->Evaluate(x, y, 0, 0);
                                    _transformedSourceGradient[n](i, j, 0, 1) = _interpolatorGradient[n]->Evaluate(x, y, 0, 1);
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
                    ptr2disp += 3;
                    ptr2latt += 3;
                }
            }
        } else {
            for (k = 0; k < _target->GetZ(); k++) {
                for (j = 0; j < _target->GetY(); j++) {
                    for (i = 0; i < _target->GetX(); i++) {
                        if (_target->Get(i, j, k) >= 0) {
                            x = ptr2latt[0];
                            y = ptr2latt[1];
                            z = ptr2latt[2];
                            t = _affd->TimeToLattice(_t_real[n]);
                            _affd->FFD1(x, y, z, t);
                            x += ptr2disp[0];
                            y += ptr2disp[1];
                            z += ptr2disp[2];
                            _source[n]->WorldToImage(x, y, z);

                            // Check whether transformed point is inside volume
                            if ((x > 0) && (x < _source[n]->GetX()-1) &&
                                (y > 0) && (y < _source[n]->GetY()-1) &&
                                (z > 0) && (z < _source[n]->GetZ()-1)) {

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
                                        ptr1 = (short *)_source[n]->GetScalarPointer(a, b, c);
                                        _transformedSource[n](i, j, k) = (t1 * (u2 * (v2 * ptr1[offset2] + v1 * ptr1[offset6]) +
                                            u1 * (v2 * ptr1[offset4] + v1 * ptr1[offset8])) +
                                            t2 * (u2 * (v2 * ptr1[offset1] + v1 * ptr1[offset5]) +
                                            u1 * (v2 * ptr1[offset3] + v1 * ptr1[offset7])));

                                        // Linear interpolation in gradient image
                                        ptr2 = _sourceGradient[n].GetPointerToVoxels(a, b, c, 0);
                                        _transformedSourceGradient[n](i, j, k, 0) = (t1 * (u2 * (v2 * ptr2[offset2] + v1 * ptr2[offset6]) +
                                            u1 * (v2 * ptr2[offset4] + v1 * ptr2[offset8])) +
                                            t2 * (u2 * (v2 * ptr2[offset1] + v1 * ptr2[offset5]) +
                                            u1 * (v2 * ptr2[offset3] + v1 * ptr2[offset7])));
                                        ptr2 = _sourceGradient[n].GetPointerToVoxels(a, b, c, 1);
                                        _transformedSourceGradient[n](i, j, k, 1) = (t1 * (u2 * (v2 * ptr2[offset2] + v1 * ptr2[offset6]) +
                                            u1 * (v2 * ptr2[offset4] + v1 * ptr2[offset8])) +
                                            t2 * (u2 * (v2 * ptr2[offset1] + v1 * ptr2[offset5]) +
                                            u1 * (v2 * ptr2[offset3] + v1 * ptr2[offset7])));
                                        ptr2 = _sourceGradient[n].GetPointerToVoxels(a, b, c, 2);
                                        _transformedSourceGradient[n](i, j, k, 2) = (t1 * (u2 * (v2 * ptr2[offset2] + v1 * ptr2[offset6]) +
                                            u1 * (v2 * ptr2[offset4] + v1 * ptr2[offset8])) +
                                            t2 * (u2 * (v2 * ptr2[offset1] + v1 * ptr2[offset5]) +
                                            u1 * (v2 * ptr2[offset3] + v1 * ptr2[offset7])));
                                    } else {
                                        // Interpolation in source image
                                        _transformedSource[n](i, j, k) = _interpolator[n]->Evaluate(x, y, z);

                                        // Interpolation in gradient image
                                        _transformedSourceGradient[n](i, j, k, 0) = _interpolatorGradient[n]->Evaluate(x, y, z, 0);
                                        _transformedSourceGradient[n](i, j, k, 1) = _interpolatorGradient[n]->Evaluate(x, y, z, 1);
                                        _transformedSourceGradient[n](i, j, k, 2) = _interpolatorGradient[n]->Evaluate(x, y, z, 2);
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
                        ptr2disp += 3;
                        ptr2latt += 3;
                    }
                }
            }
        }
    }
    //cout<<"irtkImageTFFDRegistration::UpdateSourceAndGradient end"<<endl;
#ifdef USE_TIMING
    // Stop timing
    end = clock();
    cpu_time_used = ((double) (end - start)) / CLOCKS_PER_SEC;
    cout << "CPU time for irtkImageTFFDRegistration::UpdateSourceAndGradient() = " << cpu_time_used << endl;
#endif

}

void irtkImageTFFDRegistration::Update(bool updateGradient)
{
    // Print debugging information
    this->Debug("irtkImageTFFDRegistration::Update()");

    // Finalize base class
    this->irtkTemporalImageRegistration::Update(updateGradient);

    if(_Lambda2 > 0) {
        int i, j, k, l;
        double x,y,z,t,jacobian;
        // Update Jacobian and jacobian determinants;
        for (i = 0; i < _affd->GetX(); i++) {
            for (j = 0; j < _affd->GetY(); j++) {
                for (k = 0; k < _affd->GetZ(); k++) {
                    for (l = 0; l < _affd->GetT(); l++) {
                        x = i;
                        y = j;
                        z = k;
                        t = l;
                        _affd->LatticeToWorld(x, y, z);
                        _affd->LatticeToTime(t);
                        irtkMatrix jac;
                        jac.Initialize(3, 3);
                        _mffd->LocalJacobian(jac, x, y, z, t);
                        jac.Adjugate(jacobian);
                        if (jacobian < 0.0000001) jacobian = 0.0000001;
                        _determinant[((l*_affd->GetZ() + k)*_affd->GetY() + j)*_affd->GetX() + i] = jacobian;
                        _adjugate[((l*_affd->GetZ() + k)*_affd->GetY() + j)*_affd->GetX() + i] = jac;
                    }
                }
            }
        }
    }
}

double irtkImageTFFDRegistration::SmoothnessPenalty()
{
    cerr << "irtkImageTFFDRegistration::SmoothnessPenalty: Not implemented yet" << endl;
    exit(1);

    //  if (_affd->GetZ() == 1) {
    //    return -_affd->Bending() / double(2.0*_affd->GetX()*_affd->GetY());
    //  } else {
    //    return -_affd->Bending() / double(3.0*_affd->GetX()*_affd->GetY()*_affd->GetZ());
    //  }
}

void irtkImageTFFDRegistration::SmoothnessPenaltyGradient(double *gradient)
{
    cerr << "irtkImageTFFDRegistration::SmoothnessPenaltyGradient: Not implemented yet" << endl;
    exit(1);

    //  int i;
    //  double norm;
    //
    //  // Compute normalization factor
    //  if (_affd->GetZ() == 1) {
    //    norm = (double(_source->GetNumberOfVoxels()) / double(2.0*_affd->GetX()*_affd->GetY()));
    //  } else {
    //    norm = (double(_source->GetNumberOfVoxels()) / double(3.0*_affd->GetX()*_affd->GetY()*_affd->GetZ()));
    //  }
    //
    //  // Allocate memory
    //  double *tmp_gradient = new double[_affd->NumberOfDOFs()];
    //
    //  // and initialize memory (thanks to Stefan for his bug fix)
    //  for (i = 0; i < _affd->NumberOfDOFs(); i++) {
    //    tmp_gradient[i] = 0.0;
    //  }
    //
    //  // Compute gradient of smoothness term
    //  _affd->BendingGradient(tmp_gradient);
    //
    //  // Add gradient to existing gradient
    //  for (i = 0; i < _affd->NumberOfDOFs(); i++) {
    //    gradient[i] += this->_Lambda1 * tmp_gradient[i] * norm;
    //  }
    //
    //  // Free memory
    //  delete []tmp_gradient;
}

double irtkImageTFFDRegistration::VolumePreservationPenalty()
{
    cerr << "irtkImageTFFDRegistration::VolumePreservationPenalty: Not implemented yet" << endl;
    exit(1);

    //  int k;
    //  double penalty, jacobian;
    //
    //  penalty = 0;
    //  for (k = 0; k < _affd->NumberOfDOFs()/3; k++) {
    //    // Determinant of Jacobian of deformation derivatives
    //    jacobian = _determinant[k];
    //    penalty += pow(log(jacobian), 2.0);
    //  }
    //
    //  // Normalize sum by number of DOFs
    //  return -penalty / (double) _affd->NumberOfDOFs();
}

void irtkImageTFFDRegistration::VolumePreservationPenaltyGradient(double *gradient)
{
    cerr << "irtkImageTFFDRegistration::VolumePreservationPenaltyGradient: Not implemented yet" << endl;
    exit(1);

    //  int i, j, k, l, m, n, o, i1, j1, k1, i2, j2, k2, x, y, z, count, index, index1, index2, index3;
    //  double jacobian, drv[3];
    //  irtkMatrix jac, det_drv[3];
    //
    //  double norm;
    //
    //  // Compute normalization factor
    //  if (_affd->GetZ() == 1) {
    //      norm = (double(_target->GetNumberOfVoxels()) / double(2.0*_affd->GetX()*_affd->GetY()));
    //  } else {
    //      norm = (double(_target->GetNumberOfVoxels()) / double(3.0*_affd->GetX()*_affd->GetY()*_affd->GetZ()));
    //  }
    //
    //  // Loop over control points
    //  for (z = 0; z < _affd->GetZ(); z++) {
    //    for (y = 0; y < _affd->GetY(); y++) {
    //      for (x = 0; x < _affd->GetX(); x++) {
    //
    //        // Compute DoFs corresponding to the control point
    //        index  = _affd->LatticeToIndex(x, y, z);
    //        index2 = index+_affd->GetX()*_affd->GetY()*_affd->GetZ();
    //        index3 = index+2*_affd->GetX()*_affd->GetY()*_affd->GetZ();
    //
    //        // Check if any DoF corresponding to the control point is active
    //        if ((_affd->irtkTransformation::GetStatusCP(index)  == _Active) ||
    //            (_affd->irtkTransformation::GetStatusCP(index2) == _Active) ||
    //            (_affd->irtkTransformation::GetStatusCP(index3) == _Active)) {
    //
    //          _affd->IndexToLattice(index, i, j, k);
    //          l = i;
    //          m = j;
    //          n = k;
    //          count = 0;
    //          k1 = (k-1) > 0 ? (k-1) : 0;
    //          j1 = (j-1) > 0 ? (j-1) : 0;
    //          i1 = (i-1) > 0 ? (i-1) : 0;
    //          k2 = (k+2) < _affd->GetZ() ? (k+2) : _affd->GetZ();
    //          j2 = (j+2) < _affd->GetY() ? (j+2) : _affd->GetY();
    //          i2 = (i+2) < _affd->GetX() ? (i+2) : _affd->GetX();
    //
    //          for (i = 0; i < 3; i++) drv[i] = 0;
    //
    //          for (k = k1; k < k2; k++) {
    //            for (j = j1; j < j2; j++) {
    //              for (i = i1; i < i2; i++) {
    //                if(k != n || j != m || i != l) {
    //                  index1 = _affd->LatticeToIndex(i,j,k);
    //                  // Torsten Rohlfing et al. MICCAI'01 (w/o scaling correction):
    //                  // Calculate jacobian
    //                  irtkMatrix jac = _adjugate[index1];
    //
    //                  // find jacobian derivatives
    //                  jacobian = _determinant[index1];
    //
    //                  // if jacobian < 0
    //                  jacobian = (2.0*log(jacobian))/jacobian;
    //
    //                  _affd->JacobianDetDerivative(det_drv,i-l,j-m,k-n);
    //
    //                  double tmpdrv[3];
    //
    //                  for(o = 0; o < 3; o++) {
    //                    // trace * adjugate * derivative
    //                    tmpdrv[o] = jac(0,o)*det_drv[o](o,0) + jac(1,o)*det_drv[o](o,1) + jac(2,o)*det_drv[o](o,2);
    //                    // * rest of the volume preservation derivative
    //                    drv[o] += (jacobian*tmpdrv[o]);
    //                  }
    //                  count ++;
    //                }
    //              }
    //            }
    //          }
    //
    //          for (l = 0; l < 3; l++) drv[l] = -drv[l]/count;
    //
    //          gradient[index]  += this->_Lambda2  * drv[0] * norm;
    //          gradient[index2] += this->_Lambda2  * drv[1] * norm;
    //          gradient[index3] += this->_Lambda2  * drv[2] * norm;
    //        }
    //      }
    //    }
    //  }

    return;
}

double irtkImageTFFDRegistration::Evaluate()
{
    double tmp, similarity;

    // Evaluate similarity
    similarity = this->irtkTemporalImageRegistration::Evaluate();
    cout << "Similarity = " << similarity << "\t";

    // Add penalty for smoothness
    if (this->_Lambda1 > 0) {
        tmp = this->_Lambda1*this->SmoothnessPenalty();
        cout << "Bending = " << tmp << "\t";
        similarity += tmp;
    }
    // Add penalty for volume preservation
    if (this->_Lambda2 > 0) {
        tmp = this->_Lambda2*this->VolumePreservationPenalty();
        cout << "Volume = " << tmp;
        similarity += tmp;
    }
    if ((this->_Lambda1 > 0) || (this->_Lambda2 > 0)) cout << endl;

    //Return similarity measure + penalty terms
    return similarity;
}

void irtkImageTFFDRegistration::EvaluateGradient2D(double *gradient)
{
    double basis, pos[3], t1, t2, tt;
    int i, j, n, i1, i2, j1, j2, k1, k2, x, y, t, index;
    _Status stat[3];
    double dist, dist1, dist2;
    double x1, y1, z1, x2, y2, z2, xi, yi, zi;

    // Initialize gradient to zero
    for (i = 0; i < _affd->NumberOfDOFs(); i++) {
        gradient[i] = 0;
    }

    // Loop over control points
    for (y = 0; y < _affd->GetY(); y++) {
        for (x = 0; x < _affd->GetX(); x++) {
            for (t = 0; t < _affd->GetT(); t++) {

                // Get status of DoFs corresponding to the control point
                _affd->GetStatusCP(x, y, 0, t, stat[0], stat[1], stat[2]);
                // Check if any DoF corresponding to the control point is active
                if ((stat[0] == _Active) || (stat[1] == _Active) || (stat[2] == _Active)) {

                    // If so, calculate bounding box of control point in image coordinates
                    // Note: temporal bounding box [t1,t2] is in world time and doesn't correspond to the indices of the target images
                    index = _affd->LatticeToIndex(x, y, 0, t);

                    // loop over all target images
                    for (n = 0; n < _N_source; n++) {

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
                                for (j = 0; j < _target->GetY(); j++) {
                                    for (i = 0; i < _target->GetX(); i++) {

                                        // check whether point is in bounding box
                                        xi = i;
                                        yi = j;
                                        zi = 0;
                                        _target->ImageToWorld(xi, yi, zi);
                                        if (   (xi>=x1) && (xi<=x2)
                                            && (yi>=y1) && (yi<=y2)
                                            && (zi>=z1) && (zi<=z2)) {

                                                // Check whether reference point is valid
                                                if ((_target->Get(i, j, 0) >= 0) && (_transformedSource[n](i, j, 0) >= 0)) {

                                                    // Convert position from voxel coordinates to world coordinates
                                                    pos[0] = i;
                                                    pos[1] = j;
                                                    pos[2] = 0;
                                                    _target->ImageToWorld(pos[0], pos[1], pos[2]);

                                                    // Convert world coordinates into lattice coordinates
                                                    _affd->WorldToLattice(pos[0], pos[1], pos[2]);

                                                    // Compute B-spline tensor product at pos
                                                    basis = _affd->B(pos[0] - x) * _affd->B(pos[1] - y);
                                                    // with time:
                                                    dist1 = tt - t;
                                                    dist2 = (dist1<0) ? ((tt - _affd->GetT()+1) - t) : ((tt + _affd->GetT()-1) - t);
                                                    dist = (abs(dist1)<abs(dist2)) ? dist1 : dist2;
                                                    basis *= _affd->B(dist);

                                                    // Convert voxel-based gradient into gradient with respect to parameters (chain rule)
                                                    // NOTE: This currently assumes that the control points displacements are aligned with the world coordinate displacements
                                                    gradient[(((0*_affd->GetT() + t)*_affd->GetZ() + 0)*_affd->GetY() + y)*_affd->GetX() + x] += basis * _similarityGradient[n](i, j, 0, 0);
                                                    gradient[(((1*_affd->GetT() + t)*_affd->GetZ() + 0)*_affd->GetY() + y)*_affd->GetX() + x] += basis * _similarityGradient[n](i, j, 0, 1);
                                                    gradient[(((2*_affd->GetT() + t)*_affd->GetZ() + 0)*_affd->GetY() + y)*_affd->GetX() + x] += 0;
                                                }
                                        }
                                    }
                                }
                                //		      // Loop over all voxels in the current target (reference) volume
                                //		      for (j = j1; j <= j2; j++) {
                                //			    for (i = i1; i <= i2; i++) {
                                //
                                //				  // Check whether reference point is valid
                                //				  if ((_target[n]->Get(i, j, 0) >= 0) && (_transformedSource[n](i, j, 0) >= 0)) {
                                //
                                //					// Convert position from voxel coordinates to world coordinates
                                //					pos[0] = i;
                                //					pos[1] = j;
                                //					pos[2] = 0;
                                //					_target[n]->ImageToWorld(pos[0], pos[1], pos[2]);
                                //
                                //					// Convert world coordinates into lattice coordinates
                                //					_affd->WorldToLattice(pos[0], pos[1], pos[2]);
                                //
                                //					// Compute B-spline tensor product at pos
                                //					basis = _affd->B(pos[0] - x) * _affd->B(pos[1] - y);
                                //					// with time:
                                //					dist1 = tt - t;
                                //					dist2 = (dist1<0) ? ((tt - _affd->GetT()+1) - t) : ((tt + _affd->GetT()-1) - t);
                                //					dist = (abs(dist1)<abs(dist2)) ? dist1 : dist2;
                                //					basis *= _affd->B(dist);
                                //
                                //					// Convert voxel-based gradient into gradient with respect to parameters (chain rule)
                                //					// NOTE: This currently assumes that the control points displacements are aligned with the world coordinate displacements
                                //					gradient[(((0*_affd->GetT() + t)*_affd->GetZ() + 0)*_affd->GetY() + y)*_affd->GetX() + x] += basis * _similarityGradient[n](i, j, 0, 0);
                                //					gradient[(((1*_affd->GetT() + t)*_affd->GetZ() + 0)*_affd->GetY() + y)*_affd->GetX() + x] += basis * _similarityGradient[n](i, j, 0, 1);
                                //					gradient[(((2*_affd->GetT() + t)*_affd->GetZ() + 0)*_affd->GetY() + y)*_affd->GetX() + x] += 0;
                                //				  }
                                //			    }
                                //			  }
                        }
                    }
                }
            }
        }
    }
}

void irtkImageTFFDRegistration::EvaluateGradient3D(double *gradient)
{
    double basis, *ptr, t1, t2, tt;
    int i, j, k, n, i1, i2, j1, j2, k1, k2, x, y, z, t, index;
    _Status stat[3];
    double dist, dist1, dist2;
    double x1, y1, z1, x2, y2, z2, xi, yi, zi;
    double pos[3];

    // Initialize gradient to zero
    for (i = 0; i < _affd->NumberOfDOFs(); i++) {
        gradient[i] = 0.;
    }
    cout<<"irtkImageTFFDRegistration::EvaluateGradient3D start"<<endl;
    int bla = 0, blablub = _affd->NumberOfDOFs()/3;
    cout<<"number of control points: "<<blablub<<";";
    cout.flush();
    //cout<<bla<<"/"<<blablub<<";";
    // Loop over control points
    for (z = 0; z < _affd->GetZ(); z++) {
        for (y = 0; y < _affd->GetY(); y++) {
            for (x = 0; x < _affd->GetX(); x++) {
                for (t = 0; t < _affd->GetT(); t++) {
                    //cout<<"\r";
                    //cout<<bla<<"/"<<blablub<<";"; bla++;

                    // Get status of DoFs corresponding to the control point
                    _affd->GetStatusCP(x, y, z, t, stat[0], stat[1], stat[2]);
                    // Check if any DoF corresponding to the control point is active
                    if ((stat[0] == _Active) || (stat[1] == _Active) || (stat[2] == _Active)) {

                        // If so, calculate bounding box of control point in image coordinates
                        // Note: temporal bounding box [t1,t2] is in world time and doesn't correspond to the indices of the target images
                        index = _affd->LatticeToIndex(x, y, z, t);
                        // t1, t2 not in lattice coordinates at the moment!!!!!!!
                        // spatial coordinates in world system
                        //		    _affd->BoundingBoxCP(index, x1, y1, z1, t1, x2, y2, z2, t2, 1.0);

                        // loop over all target images
                        for (n = 0; n < _N_source; n++) {

                            // t1, t2 not in lattice coordinates at the moment!!!!!!!
                            _affd->BoundingBoxImage(_target, index, i1, j1, k1, i2, j2, k2, t1, t2, 1.0);

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
                                || ( (t1 <  0) && 							( (tt <= t2) || (tt >= t1+_affd->GetT()-1) ) )
                                || ( (t2 >= _affd->GetT()-1) && 			    ( (tt >= t1) || (tt <= t2-_affd->GetT()+1) ) ) ) {

                                    //		    	// Loop over all voxels in the target (reference) volume
                                    //		    	ptr = _latticeCoordLUT[n];
                                    //		    	for (k = 0; k < _target[n]->GetZ(); k++) {
                                    //				  for (j = 0; j < _target[n]->GetY(); j++) {
                                    //		    	    for (i = 0; i < _target[n]->GetX(); i++) {
                                    //
                                    //		    	      // check whether point is in bounding box
                                    //		    	      xi = i;
                                    //		    	      yi = j;
                                    //		    	      zi = k;
                                    //		    	      _target[n]->ImageToWorld(xi, yi, zi);
                                    //		    	      if (   ((x1<=xi)&&(xi<=x2)) || ((x2<=xi)&&(xi<=x1))
                                    //		    	    	  && ((y1<=yi)&&(yi<=y2)) || ((y2<=yi)&&(yi<=y1))
                                    //		    	    	  && ((z1<=zi)&&(zi<=z2)) || ((z2<=zi)&&(zi<=z1))) {
                                    //					    // Check whether reference point is valid
                                    //					    if ((_target[n]->Get(i, j, k) >= 0) && (_transformedSource[n](i, j, k) >= 0)) {
                                    //
                                    //						  // Compute B-spline tensor product at pos
                                    //						  basis = _affd->B(ptr[0] - x) * _affd->B(ptr[1] - y) * _affd->B(ptr[2] - z);
                                    //						  // with time:
                                    //						  dist1 = tt - t;
                                    //						  dist2 = (dist1<0) ? ((tt - _affd->GetT()+1) - t) : ((tt + _affd->GetT()-1) - t);
                                    //						  dist = (abs(dist1)<abs(dist2)) ? dist1 : dist2;
                                    //						  basis *= _affd->B(dist);
                                    //
                                    //						  // Convert voxel-based gradient into gradient with respect to parameters (chain rule)
                                    //						  // NOTE: This currently assumes that the control points displacements are aligned with the world coordinate displacements
                                    //						  gradient[(((0*_affd->GetT() + t)*_affd->GetZ() + z)*_affd->GetY() + y)*_affd->GetX() + x] += basis * _similarityGradient[n](i, j, k, 0);
                                    //						  gradient[(((1*_affd->GetT() + t)*_affd->GetZ() + z)*_affd->GetY() + y)*_affd->GetX() + x] += basis * _similarityGradient[n](i, j, k, 1);
                                    //						  gradient[(((2*_affd->GetT() + t)*_affd->GetZ() + z)*_affd->GetY() + y)*_affd->GetX() + x] += basis * _similarityGradient[n](i, j, k, 2);
                                    //					    }
                                    //		    	      }
                                    //					  ptr += 3;
                                    //					}
                                    //				  }
                                    //				}
                                    // Loop over all voxels in the target (reference) volume
                                    for (k = k1; k <= k2; k++) {
                                        for (j = j1; j <= j2; j++) {
                                            ptr = &(_latticeCoordLUT[n][3 * (k * (_target->GetX()*_target->GetY()) + j * _target->GetX() + i1)]);
                                            for (i = i1; i <= i2; i++) {

                                                // Check whether reference point is valid
                                                if ((_target->Get(i, j, k) >= 0) && (_transformedSource[n](i, j, k) >= 0)) {

                                                    // Convert world coordinates into lattice coordinates
                                                    pos[0] = i;
                                                    pos[1] = j;
                                                    pos[2] = k;
                                                    _target->ImageToWorld(pos[0], pos[1], pos[2]);
                                                    _affd->WorldToLattice(pos[0], pos[1], pos[2]);

                                                    // Compute B-spline tensor product at pos
                                                    //					    basis = _affd->B(ptr[0] - x) * _affd->B(ptr[1] - y) * _affd->B(ptr[2] - z);
                                                    basis = _affd->B(pos[0] - x) * _affd->B(pos[1] - y) * _affd->B(pos[2] - z);
                                                    // with time:
                                                    dist1 = tt - t;
                                                    dist2 = (dist1<0) ? ((tt - _affd->GetT()+1) - t) : ((tt + _affd->GetT()-1) - t);
                                                    dist = (abs(dist1)<abs(dist2)) ? dist1 : dist2;
                                                    basis *= _affd->B(dist);

                                                    // Convert voxel-based gradient into gradient with respect to parameters (chain rule)
                                                    // NOTE: This currently assumes that the control points displacements are aligned with the world coordinate displacements
                                                    gradient[(((0*_affd->GetT() + t)*_affd->GetZ() + z)*_affd->GetY() + y)*_affd->GetX() + x] += basis * _similarityGradient[n](i, j, k, 0);
                                                    gradient[(((1*_affd->GetT() + t)*_affd->GetZ() + z)*_affd->GetY() + y)*_affd->GetX() + x] += basis * _similarityGradient[n](i, j, k, 1);
                                                    gradient[(((2*_affd->GetT() + t)*_affd->GetZ() + z)*_affd->GetY() + y)*_affd->GetX() + x] += basis * _similarityGradient[n](i, j, k, 2);
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
    cout<<endl;
    //  // gradient[_t-1] = gradient[0]
    //  for (z = 0; z < _affd->GetZ(); z++) {
    //  	for (y = 0; y < _affd->GetY(); y++) {
    //  	  for (x = 0; x < _affd->GetX(); x++) {
    //  		gradient[(((0*_affd->GetT() + _affd->GetT()-1)*_affd->GetZ() + z)*_affd->GetY() + y)*_affd->GetX() + x] = gradient[(((0*_affd->GetT() + 0)*_affd->GetZ() + z)*_affd->GetY() + y)*_affd->GetX() + x];
    //		gradient[(((1*_affd->GetT() + _affd->GetT()-1)*_affd->GetZ() + z)*_affd->GetY() + y)*_affd->GetX() + x] = gradient[(((1*_affd->GetT() + 0)*_affd->GetZ() + z)*_affd->GetY() + y)*_affd->GetX() + x];
    //		gradient[(((2*_affd->GetT() + _affd->GetT()-1)*_affd->GetZ() + z)*_affd->GetY() + y)*_affd->GetX() + x] = gradient[(((2*_affd->GetT() + 0)*_affd->GetZ() + z)*_affd->GetY() + y)*_affd->GetX() + x];
    //  	  }
    //  	}
    //  }
    cout<<"irtkImageTFFDRegistration::EvaluateGradient3D end"<<endl;
}

void irtkImageTFFDRegistration::NormalizeGradient(double *gradient)
{
    int x,y,z,t, index,index1,index2;
    double norm,spacingnorm;

    // Compute normalization factor
    if (_affd->GetZ() == 1) {
        spacingnorm = (double(_target->GetNumberOfVoxels()) / double(2.0*_affd->GetX()*_affd->GetY()));
    } else {
        spacingnorm = (double(_target->GetNumberOfVoxels()) / double(3.0*_affd->GetX()*_affd->GetY()*_affd->GetZ()));
    }

    for (z = 0; z < _affd->GetZ(); z++) {
        for (y = 0; y < _affd->GetY(); y++) {
            for (x = 0; x < _affd->GetX(); x++) {
                for (t = 0; t < _affd->GetT(); t++) {
                    //weired to me
                    index  = (((0*_affd->GetT() + t)*_affd->GetZ() + z)*_affd->GetY() + y)*_affd->GetX() + x;
                    index1 = (((1*_affd->GetT() + t)*_affd->GetZ() + z)*_affd->GetY() + y)*_affd->GetX() + x;
                    index2 = (((2*_affd->GetT() + t)*_affd->GetZ() + z)*_affd->GetY() + y)*_affd->GetX() + x;
                    norm = pow(gradient[index ],2.0)
                        + pow(gradient[index1],2.0)
                        + pow(gradient[index2],2.0);

                    //normalize
                    if(norm > 0){
                        norm = sqrt(norm);
                        gradient[index ] = gradient[index ]/pow((norm + _Epsilon*spacingnorm),_Lambda3);
                        gradient[index1] = gradient[index1]/pow((norm + _Epsilon*spacingnorm),_Lambda3);
                        gradient[index2] = gradient[index2]/pow((norm + _Epsilon*spacingnorm),_Lambda3);
                    }
                }
            }
        }
    }
}

double irtkImageTFFDRegistration::EvaluateGradient(double *gradient)
{
    double norm, max_length;
    int i, x, y, z, t, index, index2, index3;
    static double *g = NULL, *h = NULL, gg, dgg, gamma;

    // Compute gradient with respect to displacements
    this->irtkTemporalImageRegistration::EvaluateGradient(gradient);

    // Start timing
    clock_t start, end;
    double cpu_time_used;
    start = clock();

    if (_affd->GetZ() == 1) {
        this->EvaluateGradient2D(gradient);
    } else {
        this->EvaluateGradient3D(gradient);
    }

    // Update gradient
    if (_CurrentIteration == 0) {
        // First iteration, so let's initialize
        if (g != NULL) delete []g;
        g = new double [_affd->NumberOfDOFs()];
        if (h != NULL) delete []h;
        h = new double [_affd->NumberOfDOFs()];
        for (i = 0; i < _affd->NumberOfDOFs(); i++) {
            g[i] = -gradient[i];
            h[i] = g[i];
        }
    } else {
        // Update gradient direction to be conjugate
        gg = 0;
        dgg = 0;
        for (i = 0; i < _affd->NumberOfDOFs(); i++) {
            gg  += g[i]*h[i];
            dgg += (gradient[i]+g[i])*gradient[i];
        }
        gamma = dgg/gg;
        for (i = 0; i < _affd->NumberOfDOFs(); i++) {
            g[i] = -gradient[i];
            h[i] = g[i] + gamma*h[i];
            gradient[i] = -h[i];
        }
    }

    if (this->_Lambda1 > 0) {
        this->SmoothnessPenaltyGradient(gradient);
    }

    if (this->_Lambda2 > 0) {
        this->VolumePreservationPenaltyGradient(gradient);
    }

    if(_Lambda3 > 0){
        this->NormalizeGradient(gradient);
    }

    // Calculate maximum of gradient vector
    max_length = 0;
    for (z = 0; z < _affd->GetZ(); z++) {
        for (y = 0; y < _affd->GetY(); y++) {
            for (x = 0; x < _affd->GetX(); x++) {
                for (t = 0; t < _affd->GetT(); t++) {
                    index  = (((0*_affd->GetT() + t)*_affd->GetZ() + z)*_affd->GetY() + y)*_affd->GetX() + x;
                    index2 = (((1*_affd->GetT() + t)*_affd->GetZ() + z)*_affd->GetY() + y)*_affd->GetX() + x;
                    index3 = (((2*_affd->GetT() + t)*_affd->GetZ() + z)*_affd->GetY() + y)*_affd->GetX() + x;
                    norm = sqrt(gradient[index] * gradient[index] + gradient[index2] * gradient[index2] + gradient[index3] * gradient[index3]);
                    if (norm > max_length) max_length = norm;
                }
            }
        }
    }

    // Deal with active and passive control points
    for (i = 0; i < _affd->NumberOfDOFs(); i++) {
        if (_affd->irtkTransformation::GetStatus(i) == _Passive) {
            gradient[i] = 0;
        }
    }

    // Stop timing
    end = clock();
    cpu_time_used = ((double) (end - start)) / CLOCKS_PER_SEC;
    //cout << "CPU time for irtkImageTFFDRegistration::EvaluateGradient() = " << cpu_time_used << endl;

    return max_length;
}

void irtkImageTFFDRegistration::Run()
{
    int i, k;
    char buffer[256];
    double *gradient, delta, step, min_step, max_step, max_length, best_similarity, new_similarity, old_similarity;

    // Print debugging information
    this->Debug("irtkImageTFFDRegistration::Run");

    if (_source == NULL) {
        cerr << "irtkImageTFFDRegistration::Run: Filter has no source input" << endl;
        exit(1);
    }

    if (_target == NULL) {
        cerr << "irtkImageTFFDRegistration::Run: Filter has no target input" << endl;
        exit(1);
    }

    if (_transformation == NULL) {
        cerr << "irtkImageTFFDRegistration::Run: Filter has no transformation output" << endl;
        exit(1);
    }

    // Do the initial set up for all levels
    this->Initialize();

    // Loop over levels
    for (_CurrentLevel = _NumberOfLevels-1; _CurrentLevel >= 0; _CurrentLevel--) {

        // Initial step size
        min_step = _MinStep[_CurrentLevel];
        max_step = _MaxStep[_CurrentLevel];

        // Print resolution level
        cout << "Resolution level no. " << _CurrentLevel+1 << " (step sizes " << min_step << " to " << max_step  << ")\n";
        // Initialize for this level
        this->Initialize(_CurrentLevel);
        cout<<"Initialize("<<_CurrentLevel<<") done";

        // Save pre-processed images if we are debugging
        //    _DebugFlag = true;
        sprintf(buffer, "target_%d.nii.gz", _CurrentLevel);
        if (_DebugFlag == true) _target->Write(buffer);
        for (int n = 0; n < _N_source; n++) {
            sprintf(buffer, "source_%d_%d.nii.gz", _CurrentLevel, n);
            if (_DebugFlag == true) _source[n]->Write(buffer);
        }

        // Allocate memory for gradient vector
        gradient = new double[_affd->NumberOfDOFs()];

        // Run the registration filter at this resolution
        _CurrentIteration = 0;
        while (_CurrentIteration < _NumberOfIterations[_CurrentLevel]) {
            cout << "Iteration = " << _CurrentIteration + 1 << " (out of " << _NumberOfIterations[_CurrentLevel] << ")"<< endl;

            // Update source image
            this->Update(true);

            //      for (int n = 0; n < _N_source; n++) {
            //    	  char * filename = NULL;
            //    	  char buffert1[5];
            //		  sprintf(buffert1, "%d", _CurrentLevel);
            //		  char buffert2[5];
            //		  sprintf(buffert2, "%d", n);
            //		  filename = new char[100];
            //		  strcpy(filename,"tmp/transformedSource_");
            //		  strcat(filename,buffert1);
            //		  strcat(filename,"_");
            //		  strcat(filename,buffert2);
            //		  strcat(filename,".nii.gz");
            //		  _transformedSource[n].Write(filename);
            //		  delete filename;
            //      }
            //      for (int n = 0; n < _N_source; n++) {
            //    	  char * filename = NULL;
            //		  char buffert1[5];
            //		  sprintf(buffert1, "%d", _CurrentLevel);
            //		  char buffert2[5];
            //		  sprintf(buffert2, "%d", n);
            //		  filename = new char[100];
            //		  strcpy(filename,"tmp/transformedSourceGradient_");
            //		  strcat(filename,buffert1);
            //		  strcat(filename,"_");
            //		  strcat(filename,buffert2);
            //		  strcat(filename,".nii.gz");
            //		  _transformedSourceGradient[n].Write(filename);
            //		  delete filename;
            //      }

            // Compute current metric value
            best_similarity = old_similarity = this->Evaluate();
            cout << "Current objective function value is " << best_similarity << endl;

            // Compute gradient of similarity metric. The function EvaluateGradient() returns the maximum control point length in the gradient
            max_length = this->EvaluateGradient(gradient);

            //      for (int n = 0; n < _N_source; n++) {
            //		  char * filename = NULL;
            //		  char buffert1[5];
            //		  sprintf(buffert1, "%d", _CurrentLevel);
            //		  char buffert2[5];
            //		  sprintf(buffert2, "%d", n);
            //		  filename = new char[100];
            //		  strcpy(filename,"tmp/similarityGradient_");
            //		  strcat(filename,buffert1);
            //		  strcat(filename,"_");
            //		  strcat(filename,buffert2);
            //		  strcat(filename,".nii.gz");
            //		  _similarityGradient[n].Write(filename);
            //		  delete filename;
            //      }

            // Step along gradient direction until no further improvement is necessary
            i = 0;
            delta = 0;
            step = max_step;
            do {
                double current = step / max_length;

                // Move along gradient direction
                for (k = 0; k < _affd->NumberOfDOFs(); k++) {
                    _affd->Put(k, _affd->Get(k) + current * gradient[k]);
                }

                // We have just changed the transformation parameters, so we need to update
                this->Update(false);

                // Compute new similarity
                new_similarity = this->Evaluate();

                if (new_similarity > best_similarity + _Epsilon) {
                    cout << "New objective value function is " << new_similarity << "; step = " << step << endl;
                    best_similarity = new_similarity;
                    delta += step;
                    step = step * 1.1;
                    if (step > max_step) step = max_step;

                } else {
                    // Last step was no improvement, so back track
                    cout << "Rejected objective function value is " << new_similarity << "; step = " << step << endl;
                    for (k = 0; k < _affd->NumberOfDOFs(); k++) {
                        _affd->Put(k, _affd->Get(k) - current * gradient[k]);
                    }
                    step = step * 0.5;
                }
                i++;
            } while ((i < MAX_NO_LINE_ITERATIONS) && (step > min_step));

            _CurrentIteration++;

            // Check for convergence
            if (delta == 0) break;
        }

        // Delete gradient
        delete gradient;

        // Do the final cleaning up for this level
        this->Finalize(_CurrentLevel);
    }

    // Do the final cleaning up for all levels
    this->Finalize();
}

bool irtkImageTFFDRegistration::Read(char *buffer1, char *buffer2, int &level)
{
    int ok = false;

    if ((strstr(buffer1, "Lambda ") != NULL) ||
        (strstr(buffer1, "Lambda1") != NULL)) {
            this->_Lambda1 = atof(buffer2);
            cout << "Lambda 1 is ... " << this->_Lambda1 << endl;
            ok = true;
    }
    if (strstr(buffer1, "Lambda2") != NULL) {
        this->_Lambda2 = atof(buffer2);
        cout << "Lambda 2 is ... " << this->_Lambda2 << endl;
        ok = true;
    }
    if (strstr(buffer1, "Lambda3") != NULL) {
        this->_Lambda3 = atof(buffer2);
        cout << "Lambda 3 is ... " << this->_Lambda3 << endl;
        ok = true;
    }
    if (strstr(buffer1, "MFFDMode") != NULL) {
        if ((strcmp(buffer2, "False") == 0) || (strcmp(buffer2, "No") == 0)) {
            this->_MFFDMode = false;
            cout << "MFFDMode is ... false" << endl;
        } else {
            if ((strcmp(buffer2, "True") == 0) || (strcmp(buffer2, "Yes") == 0)) {
                this->_MFFDMode = true;
                cout << "MFFDMode is ... true" << endl;
            } else {
                cerr << "Can't read boolean value = " << buffer2 << endl;
                exit(1);
            }
        }
        ok = true;
    }
    if (strstr(buffer1, "Control point spacing in X") != NULL) {
        this->_DX = atof(buffer2);
        cout << "Control point spacing in X is ... " << this->_DX << endl;
        ok = true;
    }
    if (strstr(buffer1, "Control point spacing in Y") != NULL) {
        this->_DY = atof(buffer2);
        cout << "Control point spacing in Y is ... " << this->_DY << endl;
        ok = true;
    }
    if (strstr(buffer1, "Control point spacing in Z") != NULL) {
        this->_DZ = atof(buffer2);
        cout << "Control point spacing in Z is ... " << this->_DZ << endl;
        ok = true;
    }
    if (strstr(buffer1, "Control point spacing in T") != NULL) {
        this->_DT = atof(buffer2);
        cout << "Control point spacing in T is ... " << this->_DT << endl;
        ok = true;
    }
    if (strstr(buffer1, "Subdivision") != NULL) {
        if ((strcmp(buffer2, "False") == 0) || (strcmp(buffer2, "No") == 0)) {
            this->_Subdivision = false;
            cout << "Subdivision is ... false" << endl;
        } else {
            if ((strcmp(buffer2, "True") == 0) || (strcmp(buffer2, "Yes") == 0)) {
                this->_Subdivision = true;
                cout << "Subdivision is ... true" << endl;
            } else {
                cerr << "Can't read boolean value = " << buffer2 << endl;
                exit(1);
            }
        }
        ok = true;
    }

    if (ok == false) {
        return this->irtkTemporalImageRegistration::Read(buffer1, buffer2, level);
    } else {
        return ok;
    }
}

void irtkImageTFFDRegistration::Write(ostream &to)
{
    to << "\n#\n# Non-rigid registration parameters\n#\n\n";
    to << "Lambda1                           = " << this->_Lambda1 << endl;
    to << "Lambda2                           = " << this->_Lambda2 << endl;
    to << "Lambda3                           = " << this->_Lambda3 << endl;
    to << "Control point spacing in X        = " << this->_DX << endl;
    to << "Control point spacing in Y        = " << this->_DY << endl;
    to << "Control point spacing in Z        = " << this->_DZ << endl;
    to << "Control point spacing in T        = " << this->_DT << endl;
    if (_Subdivision == true) {
        to << "Subdivision                       = True" << endl;
    } else {
        to << "Subdivision                       = False" << endl;
    }
    if (_MFFDMode == true) {
        to << "MFFDMode                          = True" << endl;
    } else {
        to << "MFFDMode                          = False" << endl;
    }

    this->irtkTemporalImageRegistration::Write(to);
}

