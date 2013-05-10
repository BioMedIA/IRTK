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

irtkImageFreeFormRegistration2::irtkImageFreeFormRegistration2()
{
    // Print debugging information
    this->Debug("irtkImageFreeFormRegistration2::irtkImageFreeFormRegistration2");

    // Default optimization
    _OptimizationMethod = GradientDescent;

    // Default parameters for non-rigid registration
    _Lambda1     = 0;
    _Lambda2     = 0;
    _Lambda3     = 0;
    _DX          = 20;
    _DY          = 20;
    _DZ          = 20;
    _Subdivision = true;
    _Mode        = RegisterXYZ;
    _MFFDMode    = true;
    _adjugate    = NULL;
    _determinant = NULL;
}

void irtkImageFreeFormRegistration2::GuessParameter()
{
    int i;
    double xsize, ysize, zsize, spacing;

    if ((_target == NULL) || (_source == NULL)) {
        cerr << "irtkImageFreeFormRegistration2::GuessParameter: Target and source image not found" << endl;
        exit(1);
    }

    // Default parameters for registration
    _NumberOfLevels     = 4;
    _NumberOfBins       = 64;

    // Default parameters for optimization
    _SimilarityMeasure  = NMI;
    _Epsilon            = 0.0001;

    // Read target pixel size
    _target->GetPixelSize(&xsize, &ysize, &zsize);

    // Default parameters for non-rigid registration
    _Lambda1            = 0; //recommended value 0.0001
    _Lambda2            = 0; //recommended value 1
    _Lambda3            = 0;
    _DX                 = _target->GetX() * xsize / 10.0;
    _DY                 = _target->GetY() * ysize / 10.0;
    if (_target->GetZ() > 1) {
        _DZ               = _target->GetZ() * zsize / 10.0;
    } else {
        _DZ               = 1;
    }
    _Subdivision        = true;

    // Default target parameters
    _TargetBlurring[0]      = GuessResolution(xsize, ysize, zsize) / 2.0;
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
        _TargetResolution[i][2] = _TargetResolution[i-1][2] * 2;
    }

    // Read source pixel size
    _source->GetPixelSize(&xsize, &ysize, &zsize);

    // Default source parameters
    _SourceBlurring[0]      = GuessResolution(xsize, ysize, zsize) / 2.0;
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
        _SourceResolution[i][2] = _SourceResolution[i-1][2] * 2;
    }

    // Remaining parameters
    for (i = 0; i < _NumberOfLevels; i++) {
        _NumberOfIterations[i] = 40;
        _MinStep[i]            = 0.01;
        _MaxStep[i]            = 2.0;
    }

    // Try to guess padding by looking at voxel values in all eight corners of the volume:
    // If all values are the same we assume that they correspond to the padding value
    _TargetPadding = MIN_GREY;
    if ((_target->Get(_target->GetX()-1, 0, 0)                                 == _target->Get(0, 0, 0)) &&
        (_target->Get(0, _target->GetY()-1, 0)                                 == _target->Get(0, 0, 0)) &&
        (_target->Get(0, 0, _target->GetZ()-1)                                 == _target->Get(0, 0, 0)) &&
        (_target->Get(_target->GetX()-1, _target->GetY()-1, 0)                 == _target->Get(0, 0, 0)) &&
        (_target->Get(0, _target->GetY()-1, _target->GetZ()-1)                 == _target->Get(0, 0, 0)) &&
        (_target->Get(_target->GetX()-1, 0, _target->GetZ()-1)                 == _target->Get(0, 0, 0)) &&
        (_target->Get(_target->GetX()-1, _target->GetY()-1, _target->GetZ()-1) == _target->Get(0, 0, 0))) {
            _TargetPadding = _target->Get(0, 0, 0);
    }

    _SourcePadding = MIN_GREY;
    if ((_source->Get(_source->GetX()-1, 0, 0)                                 == _source->Get(0, 0, 0)) &&
        (_source->Get(0, _source->GetY()-1, 0)                                 == _source->Get(0, 0, 0)) &&
        (_source->Get(0, 0, _source->GetZ()-1)                                 == _source->Get(0, 0, 0)) &&
        (_source->Get(_source->GetX()-1, _source->GetY()-1, 0)                 == _source->Get(0, 0, 0)) &&
        (_source->Get(0, _source->GetY()-1, _source->GetZ()-1)                 == _source->Get(0, 0, 0)) &&
        (_source->Get(_source->GetX()-1, 0, _source->GetZ()-1)                 == _source->Get(0, 0, 0)) &&
        (_source->Get(_source->GetX()-1, _source->GetY()-1, _source->GetZ()-1) == _source->Get(0, 0, 0))) {
            _SourcePadding = _source->Get(0, 0, 0);
    }
}

void irtkImageFreeFormRegistration2::Initialize()
{
    // Print debugging information
    this->Debug("irtkImageFreeFormRegistration2::Initialize");

    // Initialize base class
    this->irtkImageRegistration2::Initialize();

    // Pointer to multi-level FFD
    _mffd = (irtkMultiLevelFreeFormTransformation *)_transformation;

    // Create FFD
    if (_MFFDMode == false) {
        if (_mffd->NumberOfLevels() == 0) {
            _affd = new irtkBSplineFreeFormTransformation3D(*_target, this->_DX, this->_DY, this->_DZ);
        } else {
            _affd = (irtkBSplineFreeFormTransformation *)_mffd->PopLocalTransformation();
        }
    } else {
        _affd = new irtkBSplineFreeFormTransformation3D(*_target, this->_DX, this->_DY, this->_DZ);
    }

    _adjugate = new irtkMatrix[_affd->NumberOfDOFs()/3];
    _determinant = new double[_affd->NumberOfDOFs()/3];
}

void irtkImageFreeFormRegistration2::Initialize(int level)
{
    int i, j, k;
    double x, y, z, *ptr2latt, *ptr2disp;

    // Print debugging information
    this->Debug("irtkImageFreeFormRegistration2::Initialize(int)");

    // Initialize base class
    this->irtkImageRegistration2::Initialize(level);

    // Padding of FFD
    irtkPadding(*_target, -1, _affd);

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

    // Allocate memory for global displacements
    _displacementLUT = new double[_target->GetNumberOfVoxels() * 3];

    // Allocate memory for lattice coordinates
    _latticeCoordLUT = new double[_target->GetNumberOfVoxels() * 3];

    ptr2disp = _displacementLUT;
    ptr2latt = _latticeCoordLUT;
    for (k = 0; k < _target->GetZ(); k++) {
        for (j = 0; j < _target->GetY(); j++) {
            for (i = 0; i < _target->GetX(); i++) {
                if (_target->Get(i, j, k) >= 0) {
                    x = i;
                    y = j;
                    z = k;
                    _target->ImageToWorld(x, y, z);
                    _mffd->Transform(x, y, z);
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

void irtkImageFreeFormRegistration2::Finalize()
{
    // Push local transformation back on transformation stack
    _mffd->PushLocalTransformation(_affd);

    // Print debugging information
    this->Debug("irtkImageFreeFormRegistration2::Finalize");

    // Finalize base class
    this->irtkImageRegistration2::Finalize();
    delete []_adjugate;
    delete []_determinant;
}

void irtkImageFreeFormRegistration2::Finalize(int level)
{
    // Print debugging information
    this->Debug("irtkImageFreeFormRegistration2::Finalize(int)");

    // Finalize base class
    this->irtkImageRegistration2::Finalize(level);

    // Check if we are not at the lowest level of resolution
    if (level != 0) {
        if (this->_Subdivision == true) {
            _affd->Subdivide();
        } else {

            // Push local transformation back on transformation stack
            _mffd->PushLocalTransformation(_affd);

            // Create new FFD
            _affd = new irtkBSplineFreeFormTransformation(*_target,
                this->_DX / pow(2.0, this->_NumberOfLevels-level),
                this->_DY / pow(2.0, this->_NumberOfLevels-level),
                this->_DZ / pow(2.0, this->_NumberOfLevels-level));

            // Push local transformation back on transformation stack
            _mffd->PushLocalTransformation(_affd);
        }
    }
    delete []_adjugate;
    delete []_determinant;
    _adjugate = new irtkMatrix[_affd->NumberOfDOFs()/3];
    _determinant = new double[_affd->NumberOfDOFs()/3];
    delete []_displacementLUT;
    delete []_latticeCoordLUT;
}

void irtkImageFreeFormRegistration2::UpdateSource()
{
    short *ptr1;
    double x, y, z, t1, t2, u1, u2, v1, v2;
    int a, b, c, i, j, k, offset1, offset2, offset3, offset4, offset5, offset6, offset7, offset8;

#ifdef USE_TIMING
    // Start timing
    clock_t start, end;
    double cpu_time_used;
    start = clock();
#endif

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

    double *ptr2disp = _displacementLUT;
    double *ptr2latt = _latticeCoordLUT;
    if ((_target->GetZ() == 1) && (_source->GetZ() == 1)) {
        for (j = 0; j < _target->GetY(); j++) {
            for (i = 0; i < _target->GetX(); i++) {
                if (_target->Get(i, j, 0) >= 0) {
                    x = ptr2latt[0];
                    y = ptr2latt[1];
                    z = ptr2latt[2];
                    _affd->FFD2D(x, y);
                    z  = 0;
                    x += ptr2disp[0];
                    y += ptr2disp[1];
                    z += ptr2disp[2];
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
                        _affd->FFD3D(x, y, z);
                        x += ptr2disp[0];
                        y += ptr2disp[1];
                        z += ptr2disp[2];
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
                    ptr2disp += 3;
                    ptr2latt += 3;
                }
            }
        }
    }

#ifdef USE_TIMING
    // Stop timing
    end = clock();
    cpu_time_used = ((double) (end - start)) / CLOCKS_PER_SEC;
    cout << "CPU time for irtkImageFreeFormRegistration2::UpdateSource() = " << cpu_time_used << endl;
#endif

}

void irtkImageFreeFormRegistration2::UpdateSourceAndGradient()
{
    short *ptr1;
    double x, y, z, t1, t2, u1, u2, v1, v2, *ptr2;
    int a, b, c, i, j, k, offset1, offset2, offset3, offset4, offset5, offset6, offset7, offset8;

#ifdef USE_TIMING
    // Start timing
    clock_t start, end;
    double cpu_time_used;
    start = clock();
#endif

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

    double *ptr2disp = _displacementLUT;
    double *ptr2latt = _latticeCoordLUT;
    if ((_target->GetZ() == 1) && (_source->GetZ() == 1)) {
        for (j = 0; j < _target->GetY(); j++) {
            for (i = 0; i < _target->GetX(); i++) {
                if (_target->Get(i, j, 0) >= 0) {
                    x = ptr2latt[0];
                    y = ptr2latt[1];
                    z = ptr2latt[2];
                    _affd->FFD2D(x, y);
                    z  = 0;
                    x += ptr2disp[0];
                    y += ptr2disp[1];
                    z += ptr2disp[2];
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

                                // Linear interpolation in gradient image
                                ptr2 = _sourceGradient.GetPointerToVoxels(a, b, 0, 0);
                                _transformedSourceGradient(i, j, 0, 0) = t1 * (u2 * ptr2[offset2] + u1 * ptr2[offset4]) + t2 * (u2 * ptr2[offset1] + u1 * ptr2[offset3]);
                                ptr2 = _sourceGradient.GetPointerToVoxels(a, b, 0, 1);
                                _transformedSourceGradient(i, j, 0, 1) = t1 * (u2 * ptr2[offset2] + u1 * ptr2[offset4]) + t2 * (u2 * ptr2[offset1] + u1 * ptr2[offset3]);

                            } else {
                                // Interpolation in source image
                                _transformedSource(i, j, 0) = _interpolator->Evaluate(x, y, 0);

                                // Interpolation in gradient image
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
                        _affd->FFD3D(x, y, z);
                        x += ptr2disp[0];
                        y += ptr2disp[1];
                        z += ptr2disp[2];
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

                                    // Linear interpolation in gradient image
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

                                    // Interpolation in gradient image
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
                    ptr2disp += 3;
                    ptr2latt += 3;
                }
            }
        }
    }

#ifdef USE_TIMING
    // Stop timing
    end = clock();
    cpu_time_used = ((double) (end - start)) / CLOCKS_PER_SEC;
    cout << "CPU time for irtkImageFreeFormRegistration2::UpdateSourceAndGradient() = " << cpu_time_used << endl;
#endif

}

void irtkImageFreeFormRegistration2::Update(bool updateGradient)
{
    // Print debugging information
    this->Debug("irtkImageFreeFormRegistration2::Update()");

    // Finalize base class
    this->irtkImageRegistration2::Update(updateGradient);

    if(_Lambda2 > 0) {
        int index, index2, index3, i, j, k;
        double x,y,z,jacobian;
        // Update Jacobian and jacobian determinants;
        for (index = 0; index < _affd->NumberOfDOFs()/3; index++) {
            index2 = _affd->NumberOfDOFs()/3 + index;
            index3 = _affd->NumberOfDOFs()/3*2 + index;
            _affd->IndexToLattice(index,i,j,k);
            x = i;
            y = j;
            z = k;
            _affd->LatticeToWorld(x,y,z);
            irtkMatrix jac;
            jac.Initialize(3,3);
            _mffd->LocalJacobian(jac,x,y,z);
            jac.Adjugate(jacobian);
            if(jacobian < 0.0000001) jacobian = 0.0000001;
            _determinant[index] = jacobian;
            _adjugate[index] = jac;
        }
    }
}

double irtkImageFreeFormRegistration2::SmoothnessPenalty()
{
    if (_affd->GetZ() == 1) {
        return -_affd->Bending() / double(2.0*_affd->GetX()*_affd->GetY());
    } else {
        return -_affd->Bending() / double(3.0*_affd->GetX()*_affd->GetY()*_affd->GetZ());
    }
}

void irtkImageFreeFormRegistration2::SmoothnessPenaltyGradient(double *gradient)
{
    int i;
    double norm;

    // Compute normalization factor
    if (_affd->GetZ() == 1) {
        norm = (double(_target->GetNumberOfVoxels()) / double(2.0*_affd->GetX()*_affd->GetY()));
    } else {
        norm = (double(_target->GetNumberOfVoxels()) / double(3.0*_affd->GetX()*_affd->GetY()*_affd->GetZ()));
    }

    // Allocate memory
    double *tmp_gradient = new double[_affd->NumberOfDOFs()];

    // and initialize memory (thanks to Stefan for his bug fix)
    for (i = 0; i < _affd->NumberOfDOFs(); i++) {
        tmp_gradient[i] = 0.0;
    }

    // Compute gradient of smoothness term
    _affd->BendingGradient(tmp_gradient);

    // Add gradient to existing gradient
    for (i = 0; i < _affd->NumberOfDOFs(); i++) {
        gradient[i] += this->_Lambda1 * tmp_gradient[i] * norm;
    }

    // Free memory
    delete []tmp_gradient;
}

double irtkImageFreeFormRegistration2::VolumePreservationPenalty()
{
    int k;
    double penalty, jacobian;

    penalty = 0;
    for (k = 0; k < _affd->NumberOfDOFs()/3; k++) {
        // Determinant of Jacobian of deformation derivatives
        jacobian = _determinant[k];
        penalty += pow(log(jacobian), 2.0);
    }

    // Normalize sum by number of DOFs
    return -penalty / (double) _affd->NumberOfDOFs();
}

void irtkImageFreeFormRegistration2::VolumePreservationPenaltyGradient(double *gradient)
{
    int i, j, k, l, m, n, o, i1, j1, k1, i2, j2, k2, x, y, z, count, index, index1, index2, index3;
    double jacobian, drv[3];
    irtkMatrix jac, det_drv[3];

    double norm;

    // Compute normalization factor
    if (_affd->GetZ() == 1) {
        norm = (double(_target->GetNumberOfVoxels()) / double(2.0*_affd->GetX()*_affd->GetY()));
    } else {
        norm = (double(_target->GetNumberOfVoxels()) / double(3.0*_affd->GetX()*_affd->GetY()*_affd->GetZ()));
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
                if ((_affd->irtkTransformation::GetStatus(index)  == _Active) ||
                    (_affd->irtkTransformation::GetStatus(index2) == _Active) ||
                    (_affd->irtkTransformation::GetStatus(index3) == _Active)) {

                        _affd->IndexToLattice(index, i, j, k);
                        l = i;
                        m = j;
                        n = k;
                        count = 0;
                        k1 = (k-1) > 0 ? (k-1) : 0;
                        j1 = (j-1) > 0 ? (j-1) : 0;
                        i1 = (i-1) > 0 ? (i-1) : 0;
                        k2 = (k+2) < _affd->GetZ() ? (k+2) : _affd->GetZ();
                        j2 = (j+2) < _affd->GetY() ? (j+2) : _affd->GetY();
                        i2 = (i+2) < _affd->GetX() ? (i+2) : _affd->GetX();

                        for (i = 0; i < 3; i++) drv[i] = 0;

                        for (k = k1; k < k2; k++) {
                            for (j = j1; j < j2; j++) {
                                for (i = i1; i < i2; i++) {
                                    if(k != n || j != m || i != l) {
                                        index1 = _affd->LatticeToIndex(i,j,k);
                                        // Torsten Rohlfing et al. MICCAI'01 (w/o scaling correction):
                                        // Calculate jacobian
                                        irtkMatrix jac = _adjugate[index1];

                                        // find jacobian derivatives
                                        jacobian = _determinant[index1];

                                        // if jacobian < 0
                                        jacobian = (2.0*log(jacobian))/jacobian;

                                        _affd->JacobianDetDerivative(det_drv,i-l,j-m,k-n);

                                        double tmpdrv[3];

                                        for(o = 0; o < 3; o++) {
                                            // trace * adjugate * derivative
                                            tmpdrv[o] = jac(0,o)*det_drv[o](o,0) + jac(1,o)*det_drv[o](o,1) + jac(2,o)*det_drv[o](o,2);
                                            // * rest of the volume preservation derivative
                                            drv[o] += (jacobian*tmpdrv[o]);
                                        }
                                        count ++;
                                    }
                                }
                            }
                        }

                        for (l = 0; l < 3; l++) drv[l] = -drv[l]/count;

                        gradient[index]  += this->_Lambda2  * drv[0] * norm;
                        gradient[index2] += this->_Lambda2  * drv[1] * norm;
                        gradient[index3] += this->_Lambda2  * drv[2] * norm;
                }
            }
        }
    }

    return;
}

double irtkImageFreeFormRegistration2::Evaluate()
{
    double tmp, similarity;

    // Evaluate similarity
    similarity = this->irtkImageRegistration2::Evaluate();
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

void irtkImageFreeFormRegistration2::EvaluateGradient2D(double *gradient)
{
    double basis, pos[3], offset;
    int i, j, i1, i2, j1, j2, k1, k2, x, y, index, index2;

    // Initialize gradient to zero
    for (i = 0; i < _affd->NumberOfDOFs(); i++) {
        gradient[i] = 0;
    }

    // Loop over control points
    for (y = 0; y < _affd->GetY(); y++) {
        for (x = 0; x < _affd->GetX(); x++) {

            // Compute DoFs corresponding to the control point
            index  = _affd->LatticeToIndex(x, y, 0);
            index2 = index+_affd->GetX()*_affd->GetY()*_affd->GetZ();

            // Check if any DoF corresponding to the control point is active
            if ((_affd->irtkTransformation::GetStatus(index) == _Active) || (_affd->irtkTransformation::GetStatus(index2) == _Active)) {

                // If so, calculate bounding box of control point in image coordinates
                _affd->BoundingBoxImage(_target, index, i1, j1, k1, i2, j2, k2, 1.0);

                // Loop over all voxels in the target (reference) volume
                for (j = j1; j <= j2; j++) {
                    for (i = i1; i <= i2; i++) {

                        // Check whether reference point is valid
                        if ((_target->Get(i, j, 0) >= 0) && (_transformedSource(i, j, 0) >= 0)) {

                            // Convert position from voxel coordinates to world coordinates
                            pos[0] = i;
                            pos[1] = j;
                            pos[2] = 0;
                            _target->ImageToWorld(pos[0], pos[1], pos[2]);

                            // Convert world coordinates into lattice coordinates
                            _affd->WorldToLattice(pos[0], pos[1], pos[2]);

                            // Compute B-spline tensor product at pos
                            basis = _affd->B(pos[0] - x) * _affd->B(pos[1] - y);

                            // Convert voxel-based gradient into gradient with respect to parameters (chain rule)
                            //
                            // NOTE: This currently assumes that the control points displacements are aligned with the world coordinate displacements
                            //
                            gradient[index]  += basis * _similarityGradient(i, j, 0, 0);
                            gradient[index2] += basis * _similarityGradient(i, j, 0, 1);
                        }
                    }
                }
            }
        }
    }
}

void irtkImageFreeFormRegistration2::EvaluateGradient3D(double *gradient)
{
    double basis, *ptr;
    int i, j, k, i1, i2, j1, j2, k1, k2, x, y, z, index, index2, index3;

    // Initialize gradient to zero
    for (i = 0; i < _affd->NumberOfDOFs(); i++) {
        gradient[i] = 0;
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
                            ptr = &(_latticeCoordLUT[3 * (k * (_target->GetX()*_target->GetY()) + j * _target->GetX() + i1)]);
                            for (i = i1; i <= i2; i++) {
                                // Check whether reference point is valid
                                if ((_target->Get(i, j, k) >= 0) && (_transformedSource(i, j, k) >= 0)) {
                                    // Compute B-spline tensor product at current position
                                    basis = _affd->B(ptr[0] - x) * _affd->B(ptr[1] - y) * _affd->B(ptr[2] - z);

                                    // Convert voxel-based gradient into gradient with respect to parameters (chain rule)
                                    //
                                    // NOTE: This currently assumes that the control points displacements are aligned with the world coordinate displacements
                                    //
                                    gradient[index]  += basis * _similarityGradient(i, j, k, 0);
                                    gradient[index2] += basis * _similarityGradient(i, j, k, 1);
                                    gradient[index3] += basis * _similarityGradient(i, j, k, 2);
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

void irtkImageFreeFormRegistration2::NormalizeGradient(double *gradient)
{
    int x,y,z,index1,index2,index3,offset;
    double norm,spacingnorm;

    // Compute normalization factor
    if (_affd->GetZ() == 1) {
        spacingnorm = (double(_target->GetNumberOfVoxels()) / double(2.0*_affd->GetX()*_affd->GetY()));
    } else {
        spacingnorm = (double(_target->GetNumberOfVoxels()) / double(3.0*_affd->GetX()*_affd->GetY()*_affd->GetZ()));
    }

    offset = _affd->GetX()*_affd->GetY()*_affd->GetZ();

    for (z = 0; z < _affd->GetZ(); z++) {
        for (y = 0; y < _affd->GetY(); y++) {
            for (x = 0; x < _affd->GetX(); x++) {
                index1 = _affd->LatticeToIndex(x,y,z);
                index2 = index1 + offset;
                index3 = index2 + offset;
                norm = pow(gradient[index1],2.0)
                    + pow(gradient[index2],2.0)
                    + pow(gradient[index3],2.0);

                //normalize
                if(norm > 0){
                    norm = sqrt(norm);
                    gradient[index1] = gradient[index1]/pow((norm + _Epsilon*spacingnorm),_Lambda3);
                    gradient[index2] = gradient[index2]/pow((norm + _Epsilon*spacingnorm),_Lambda3);
                    gradient[index3] = gradient[index3]/pow((norm + _Epsilon*spacingnorm),_Lambda3);
                }
            }
        }
    }
}

double irtkImageFreeFormRegistration2::EvaluateGradient(double *gradient)
{
    double norm, max_length;
    int i, x, y, z, index, index2, index3;
    static double *g = NULL, *h = NULL, gg, dgg, gamma;

    // Compute gradient with respect to displacements
    this->irtkImageRegistration2::EvaluateGradient(gradient);

    // Start timing
    //clock_t start = clock();

    if (_affd->GetZ() == 1) {
        this->EvaluateGradient2D(gradient);
    } else {
        this->EvaluateGradient3D(gradient);
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

    // Calculate maximum of gradient vector
    max_length = 0;
    for (z = 0; z < _affd->GetZ(); z++) {
        for (y = 0; y < _affd->GetY(); y++) {
            for (x = 0; x < _affd->GetX(); x++) {
                index  = _affd->LatticeToIndex(x, y, z);
                index2 = index+_affd->GetX()*_affd->GetY()*_affd->GetZ();
                index3 = index+2*_affd->GetX()*_affd->GetY()*_affd->GetZ();
                norm = sqrt(gradient[index] * gradient[index] + gradient[index2] * gradient[index2] + gradient[index3] * gradient[index3]);
                if (norm > max_length) max_length = norm;
            }
        }
    }

    // Deal with active and passive control points
    for (i = 0; i < _affd->NumberOfDOFs(); i++) {
        if (_affd->irtkTransformation::GetStatus(i) == _Passive) {
            gradient[i] = 0;
            _affd->Put(i,0);
        }
    }

    // Stop timing
    //clock_t end = clock();
    //double cpu_time_used = ((double) (end - start)) / CLOCKS_PER_SEC;
    //cout << "CPU time for irtkImageFreeFormRegistration2::EvaluateGradient() = " << cpu_time_used << endl;

    return max_length;
}

void irtkImageFreeFormRegistration2::Run()
{
    int i, k;
    char buffer[256];
    double *gradient, delta, step, min_step, max_step, max_length, best_similarity, new_similarity, old_similarity;

    // Print debugging information
    this->Debug("irtkImageFreeFormRegistration2::Run");

    if (_source == NULL) {
        cerr << "irtkImageFreeFormRegistration2::Run: Filter has no source input" << endl;
        exit(1);
    }

    if (_target == NULL) {
        cerr << "irtkImageFreeFormRegistration2::Run: Filter has no target input" << endl;
        exit(1);
    }

    if (_transformation == NULL) {
        cerr << "irtkImageFreeFormRegistration2::Run: Filter has no transformation output" << endl;
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

        // Save pre-processed images if we are debugging
        sprintf(buffer, "source_%d.nii.gz", _CurrentLevel);
        if (_DebugFlag == true) _source->Write(buffer);
        sprintf(buffer, "target_%d.nii.gz", _CurrentLevel);
        if (_DebugFlag == true) _target->Write(buffer);

        // Allocate memory for gradient vector
        gradient = new double[_affd->NumberOfDOFs()];

        // Run the registration filter at this resolution
        _CurrentIteration = 0;
        while (_CurrentIteration < _NumberOfIterations[_CurrentLevel]) {
            cout << "Iteration = " << _CurrentIteration + 1 << " (out of " << _NumberOfIterations[_CurrentLevel] << ")"<< endl;

            // Update source image
            this->Update(true);


            if (_DebugFlag == true){
                sprintf(buffer, "transformedsource_%d_%03d.nii.gz", _CurrentLevel, _CurrentIteration);
                _transformedSource.Write(buffer);
                sprintf(buffer, "transformedgradient_%d_%03d.nii.gz", _CurrentLevel, _CurrentIteration);
                _transformedSourceGradient.Write(buffer);
            }

            // Compute current metric value
            best_similarity = old_similarity = this->Evaluate();
            cout << "Current objective function value is " << best_similarity << endl;

            // Compute gradient of similarity metric. The function EvaluateGradient() returns the maximum control point length in the gradient
            max_length = this->EvaluateGradient(gradient);

            if (_DebugFlag == true){
                sprintf(buffer, "similaritygradient_%d_%03d.nii.gz", _CurrentLevel, _CurrentIteration);
                _similarityGradient.Write(buffer);
                irtkImageAttributes attr;
                attr._x  = _affd->GetX();
                attr._y  = _affd->GetY();
                attr._z  = _affd->GetZ();
                attr._t  = 3;
                attr._dx = _affd->GetXSpacing();
                attr._dy = _affd->GetYSpacing();
                attr._dz = _affd->GetZSpacing();
                irtkGenericImage<double> g(attr);
                int idx = 0;
                for (int k = 0; k < g.GetZ(); k++) {
                  for (int j = 0; j < g.GetY(); j++) {
                    for (int i = 0; i < g.GetX(); i++) {
                      int index1 = _affd->LatticeToIndex(i, j, k);
                      int index2 = index1+_affd->GetX()*_affd->GetY()*_affd->GetZ();
                      int index3 = index2+_affd->GetX()*_affd->GetY()*_affd->GetZ();
                      g(i, j, k, 0) = gradient[index1];
                      g(i, j, k, 1) = gradient[index2];
                      g(i, j, k, 2) = gradient[index3];
                    }
                  }
                }
                sprintf(buffer, "gradient_%d_%03d.nii.gz", _CurrentLevel, _CurrentIteration);
                g.Write(buffer);
            }

            // Step along gradient direction until no further improvement is necessary
            i = 0;
            delta = 0;
            step = max_step;

            if(max_length > 0){
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
            }

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

bool irtkImageFreeFormRegistration2::Read(char *buffer1, char *buffer2, int &level)
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
        return this->irtkImageRegistration2::Read(buffer1, buffer2, level);
    } else {
        return ok;
    }
}

void irtkImageFreeFormRegistration2::Write(ostream &to)
{
    to << "\n#\n# Non-rigid registration parameters\n#\n\n";
    to << "Lambda1                           = " << this->_Lambda1 << endl;
    to << "Lambda2                           = " << this->_Lambda2 << endl;
    to << "Lambda3                           = " << this->_Lambda3 << endl;
    to << "Control point spacing in X        = " << this->_DX << endl;
    to << "Control point spacing in Y        = " << this->_DY << endl;
    to << "Control point spacing in Z        = " << this->_DZ << endl;
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

    this->irtkImageRegistration2::Write(to);
}
