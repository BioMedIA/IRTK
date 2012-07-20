/*=========================================================================

Library   : Image Registration Toolkit (IRTK)
Module    : $Id$
Copyright : Imperial College, Department of Computing
Visual Information Processing (VIP), 2009 onwards
Date      : $Date$
Version   : $Revision$
Changes   : $Author$

=========================================================================*/

#ifndef _irtkSparseFreeFormRegistration_H

#define _irtkSparseFreeFormRegistration_H

/**
* Filter for non-rigid registration based on voxel similarity measures.
*
* This class implements a registration filter for the non-rigid registration
* of two images. The basic algorithm is described in Rueckert et al., IEEE
* Transactions on Medical Imaging, In press.
*
*/

class irtkSparseFreeFormRegistration : public irtkImageRegistration2
{

protected:

    /// Pointer to the b spline gradient
    irtkMultiLevelFreeFormTransformation *_mdffd;

    /// Pointer to the global transformation which is constant
    irtkMultiLevelFreeFormTransformation *_mffd;

    /// Pointer to the local bspline transformation which is currently optimized
    irtkBSplineFreeFormTransformation *_affd;

    /// current gradient
    double *_currentgradient;

    double *_g;
    double *_h;
    double *_tmp;
    bool *_mask;

    int _NumberOfDofs;

    /// Bending energy
    double _Lambda1;

    /// Sparsity constraint lambda
    double _Lambda3;

    /// temp for sparsity constrain
    double _Lambda3tmp;

    /// sparsity log barrier
    double _Lambda2;

    /// Largest ffd grid spacing
    double _LargestSpacing;

    /// Marxium similarity
    double _MaxSimilarity;

    /// Current similarity
    double _CurrentSimilarity;

    /// Registration mode
    irtkImageFreeFormRegistrationMode _Mode;

    /// Number of sparse model levels
    int _NumberOfModels;

    /// Initial set up for the registration
    virtual void Initialize();

    /// Initial set up for the registration
    virtual void Initialize(int);

    /// Initialize transformation for registration level
    virtual void InitializeTransformation();

    /// Transformation to image (store current transformation in image)
    virtual void Transformation2Image();

    /// Image to transformation (approximate image with transformation)
    virtual void Image2Transformation();

    /// Initialize transformation for registration level
    virtual void InitializeTransformation(int);

    /// Final set up for the registration
    virtual void Finalize();

    /// Final set up for the registration
    virtual void Finalize(int);

    /// Update state of the registration based on current transformation estimate
    virtual void Update(bool);

    /// Update state of the registration based on current transformation estimate (source image)
    virtual void UpdateSource();

    /// Update state of the registration based on current transformation estimate (source image and source image gradient)
    virtual void UpdateSourceAndGradient();

    /** Evaluates the smoothness preservation term. */
    virtual double SmoothnessPenalty();

    /** Evaluates the gradient of the smoothness term. */
    virtual void SmoothnessPenaltyGradient();

    /** Evaluates the sparse preservation term. */
    virtual double SparsePenalty();

    /** Evaluates the gradient of the sparse term. */
    virtual void SparsePenaltyGradient();

    /** Evaluates the registration. This function evaluates the registration by
    *  looping over the target image and interpolating the transformed source
    *  image while filling the joint histogram. This function returns the value
    *  of the similarity measure using Similarity().
    */
    virtual double Evaluate();

    /// Evaluate parametric gradient of the MFFD model
    virtual void ParametricGradient();

    /// Evaluate the gradient of the similarity measure for the current transformation.
    virtual double EvaluateGradient(double *);

    /// Normalize similarity gradient
    virtual void NormalizeGradient();

    /// Initialize Iteration
    virtual void InitializeIteration();

    /// Finalize Iteration
    virtual void FinalizeIteration();

public:

    /// Constructor
    irtkSparseFreeFormRegistration();

    ~irtkSparseFreeFormRegistration();

    /// Set output for the registration filter
    virtual void SetOutput(irtkTransformation *);

    /// Runs the registration filter
    virtual void Run();

    /// Returns the name of the class
    virtual const char *NameOfClass();

    /// Print some information
    virtual void Print();

    /// Guess parameters
    virtual void GuessParameter();

    /// Read single line of registration parameters
    virtual bool Read(char *, char *, int &);

    /// Write registration parameters to file
    virtual void Write(ostream &);

    // Access parameters for registration mode
    virtual SetMacro(Mode, irtkImageFreeFormRegistrationMode);
    virtual GetMacro(Mode, irtkImageFreeFormRegistrationMode);
    virtual SetMacro(Lambda1, double);
    virtual GetMacro(Lambda1, double);
    virtual SetMacro(Lambda3, double);
    virtual GetMacro(Lambda3, double);

};

inline void irtkSparseFreeFormRegistration::SetOutput(irtkTransformation *transformation)
{
    // Print debugging information
    this->Debug("irtkSparseFreeFormRegistration::SetOutput");

    if (strcmp(transformation->NameOfClass(),
        "irtkMultiLevelFreeFormTransformation") != 0) {
            cerr << "irtkImageFreeFormRegistration::SetOutput: Transformation must be "
                << "irtkMultiLevelFreeFormTransformation" << endl;
            exit(0);
    }
    _transformation = transformation;
}

inline const char *irtkSparseFreeFormRegistration::NameOfClass()
{
    return "irtkFreeFormRegistration";
}

inline void irtkSparseFreeFormRegistration::Print()
{}

#endif
