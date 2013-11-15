/*=========================================================================

  Library   : Image Registration Toolkit (IRTK)
  Module    : $Id: irtkReconstructionCuda.cc 1 2013-11-15 14:36:30 bkainz $
  Copyright : Imperial College, Department of Computing
  Visual Information Processing (VIP), 2011 onwards
  Date      : $Date: 2013-11-15 14:36:30 +0100 (Fri, 15 Nov 2013) $
  Version   : $Revision: 1 $
  Changes   : $Author: bkainz $

  =========================================================================*/

#ifndef _irtkReconstructionCuda_H
#define _irtkReconstructionCuda_H

#include <irtkImage.h>
#include <irtkTransformation.h>
#include <irtkGaussianBlurring.h>
#include "reconstruction_cuda.cuh"

#include <vector>
using namespace std;

/*

  Reconstruction of volume from 2D slices

*/

/*struct POINT3D
{
    short x;
    short y;
    short z;
    double value;
};

typedef std::vector<POINT3D> VOXELCOEFFS; 
typedef std::vector<std::vector<VOXELCOEFFS> > SLICECOEFFS;*/

class irtkReconstructionCuda : public irtkObject
{

 protected:

    //Structures to store the matrix of transformation between volume and slices
    std::vector<SLICECOEFFS> _volcoeffs;

    //SLICES
    /// Slices
    vector<irtkRealImage> _slices;
    vector<irtkRealImage> _simulated_slices;
    vector<irtkRealImage> _simulated_weights;
    vector<irtkRealImage> _simulated_inside;
  
    /// Transformations
    vector<irtkRigidTransformation> _transformations;
    /// Indicator whether slice has an overlap with volumetric mask
    vector<bool> _slice_inside;
  
    //VOLUME
    /// Reconstructed volume
    irtkRealImage _reconstructed;
    /// Flag to say whether the template volume has been created
    bool _template_created;
    /// Volume mask
    irtkRealImage _mask;
    irtkRealImage _brain_probability;
    vector<irtkRealImage> _probability_maps;
  
    /// Flag to say whether we have a mask
    bool _have_mask;
    /// Weights for Gaussian reconstruction
    irtkRealImage _volume_weights;
    /// Weights for regularization
    irtkRealImage _confidence_map;
  
    //EM algorithm
    /// Variance for inlier voxel errors
    double _sigma;
    /// Proportion of inlier voxels
    double _mix;
    /// Uniform distribution for outlier voxels
    double _m;
    /// Mean for inlier slice errors
    double _mean_s;
    /// Variance for inlier slice errors
    double _sigma_s;
    /// Mean for outlier slice errors
    double _mean_s2;
    /// Variance for outlier slice errors
    double _sigma_s2;
    /// Proportion of inlier slices
    double _mix_s;
    /// Step size for likelihood calculation
    double _step;
    /// Voxel posteriors
    vector<irtkRealImage> _weights;
    ///Slice posteriors
    vector<double> _slice_weight;
   
    //Bias field
    ///Variance for bias field
    double _sigma_bias;
    /* /// Blurring object for bias field */
    /* irtkGaussianBlurring<irtkRealPixel>* _gb; */
    /// Slice-dependent bias fields
    vector<irtkRealImage> _bias;

    ///Slice-dependent scales
    vector<double> _scale;
  
    ///Quality factor - higher means slower and better
    double _quality_factor;
    ///Intensity min and max
    double _max_intensity;
    double _min_intensity;
  
    //Gradient descent and regulatization parameters
    ///Step for gradient descent
    double _alpha;
    ///Determine what is en edge in edge-preserving smoothing
    double _delta;
    ///Amount of smoothing
    double _lambda;
    ///Average voxel wights to modulate parameter alpha
    double _average_volume_weight;

    
    //global bias field correction
    ///global bias correction flag
    bool _global_bias_correction;
    ///low intensity cutoff for bias field estimation
    double _low_intensity_cutoff;
  
    //to restore original signal intensity of the MRI slices
    vector<double> _stack_factor;
    double _average_value;
    vector<int> _stack_index;
  
    //forced excluded slices
    vector<int> _force_excluded;

    //slices identify as too small to be used
    vector<int> _small_slices;

    /// use adaptive or non-adaptive regularisation (default:false)
    bool _adaptive;
    
    //utility
    ///Debug mode
    bool _debug;

  
    //Probability density functions
    ///Zero-mean Gaussian PDF
    inline double G(double x,double s);
    ///Uniform PDF
    inline double M(double m);

    int _directions[13][3];
  
	Reconstruction reconstructionGPU;

 public:

    ///Constructor
    irtkReconstructionCuda();
    ///Destructor
    ~irtkReconstructionCuda();

	Matrix4 toMatrix4(irtkMatrix mat);

	void UpdateGPUTranformationMatrices();

    ///Create zero image as a template for reconstructed volume
    double CreateTemplate( irtkRealImage stack,
                           double resolution=0 );
    double CreateLargeTemplate( vector<irtkRealImage>& stacks,
                                vector<irtkRigidTransformation>& stack_transformations,
                                double resolution,
                                double smooth_mask,
                                double threshold_mask,
                                double expand=0 );
  
    ///If template image has been masked instead of creating the mask in separate
    ///file, this function can be used to create mask from the template image
    irtkRealImage CreateMask(irtkRealImage image);
  
    ///Remember volumetric mask and smooth it if necessary
    void SetMask(irtkRealImage * mask, double sigma, double threshold=0.5 );
  
    ///Create mask from black background if the flag is set
    void CreateMaskFromBlackBackground( vector<irtkRealImage>& stacks,
                                        vector<irtkRigidTransformation>& stack_transformations,
                                        double smooth_mask );
    void CreateMaskFromAllMasks( vector<irtkRealImage> &stacks,
                                 vector<irtkRigidTransformation> &stack_transformations,
                                 double smooth_mask,
                                 double threshold_mask );
    void UpdateMaskFromAllMasks( double smooth_mask,
                                 double threshold_mask );

    void UpdateProbabilityMap();
    void SaveProbabilityMap( int i );

    void crf3DMask( double smooth_mask,
                    double threshold_mask );
  
    void CenterStacks( vector<irtkRealImage>& stacks,
                       vector<irtkRigidTransformation>& stack_transformations,
                       int templateNumber );

    //Create average image from the stacks and volumetric transformations
    irtkRealImage CreateAverage( vector<irtkRealImage>& stacks,
                                 vector<irtkRigidTransformation>& stack_transformations );

    ///Crop image according to the mask
    void CropImage( irtkRealImage& image,
                    irtkRealImage& mask );
  
    /// Transform and resample mask to the space of the image
    void TransformMask( irtkRealImage& image,
                        irtkRealImage& mask,
                        irtkRigidTransformation& transformation );

    /// Rescale image ignoring negative values
    void Rescale( irtkRealImage &img, double max);
    
    ///Calculate initial registrations
    void StackRegistrations( vector<irtkRealImage>& stacks,
                             vector<irtkRigidTransformation>& stack_transformations,
                             int templateNumber);
  
    ///Create slices from the stacks and slice-dependent transformations from
    ///stack transformations
    void CreateSlicesAndTransformations( vector<irtkRealImage>& stacks,
                                         vector<irtkRigidTransformation>& stack_transformations,
                                         vector<double>& thickness,
                                         const vector<irtkRealImage> &probability_maps=vector<irtkRealImage>() );
    void SetSlicesAndTransformations( vector<irtkRealImage>& slices,
                                      vector<irtkRigidTransformation>& slice_transformations,
                                      vector<int>& stack_ids,
                                      vector<double>& thickness );
    void ResetSlices( vector<irtkRealImage>& stacks,
                      vector<double>& thickness );

    ///Update slices if stacks have changed
    void UpdateSlices(vector<irtkRealImage>& stacks, vector<double>& thickness);
  
    void GetSlices( vector<irtkRealImage>& second_stacks );  
  
    ///Invert all stack transformation
    void InvertStackTransformations( vector<irtkRigidTransformation>& stack_transformations );

    ///Match stack intensities
    void MatchStackIntensities (vector<irtkRealImage>& stacks,
                                vector<irtkRigidTransformation>& stack_transformations,
                                double averageValue,
                                bool together=false);
 
    ///Match stack intensities with masking
    void MatchStackIntensitiesWithMasking (vector<irtkRealImage>& stacks,
                                vector<irtkRigidTransformation>& stack_transformations,
                                double averageValue,
                                bool together=false);
 
    ///Mask all slices
    void MaskSlices();
  
    ///Calculate transformation matrix between slices and voxels
    void CoeffInit();
  
    ///Reconstruction using weighted Gaussian PSF
    void GaussianReconstruction();
  
    ///Initialise variables and parameters for EM
    void InitializeEM();
  
    ///Initialise values of variables and parameters for EM
    void InitializeEMValues();
  
    ///Initalize robust statistics
    void InitializeRobustStatistics();
  
    ///Perform E-step 
    void EStep();
  
    ///Calculate slice-dependent scale
    void Scale();
  
	//sync functions GPU CPU
	void SyncGPU();
	void SyncCPU();

    ///Calculate slice-dependent bias fields
    void Bias();
    void NormaliseBias(int iter);
  
    ///Superresolution
    void Superresolution(int iter);
  
    ///Calculation of voxel-vise robust statistics
    void MStep(int iter);
  
    ///Edge-preserving regularization
    void Regularization(int iter);
  
    ///Edge-preserving regularization with confidence map
    //void AdaptiveRegularization(int iter, irtkRealImage& original);
  
    ///Slice to volume registrations
    void SliceToVolumeRegistration();
  
    ///Correct bias in the reconstructed volume
    void BiasCorrectVolume(irtkRealImage& original);
  
    ///Mask the volume
    void MaskVolume();
    void MaskImage( irtkRealImage& image, double padding=-1);
  
    ///Save slices
    void SaveSlices();
    void SlicesInfo( const char* filename, vector<string> &stack_files );
  
    ///Save weights
    void SaveWeights();
  
    ///Save transformations
    void SaveTransformations();
    void GetTransformations( vector<irtkRigidTransformation> &transformations );
    void SetTransformations( vector<irtkRigidTransformation> &transformations );
  
    ///Save confidence map
    void SaveConfidenceMap();
  
    ///Save bias field
    void SaveBiasFields();
  
    ///Remember stdev for bias field
    inline void SetSigma( double sigma );
  
    ///Return reconstructed volume
    inline irtkRealImage GetReconstructed();
    void SetReconstructed(irtkRealImage &reconstructed);
  
    ///Return resampled mask
    inline irtkRealImage GetMask();
  
    ///Set smoothing parameters
    inline void SetSmoothingParameters( double delta, double lambda );
  
    ///Use faster lower quality reconstruction
    inline void SpeedupOn();
  
    ///Use slower better quality reconstruction
    inline void SpeedupOff();
  
    ///Switch on global bias correction
    inline void GlobalBiasCorrectionOn();
  
    ///Switch off global bias correction
    inline void GlobalBiasCorrectionOff();
  
    ///Set lower threshold for low intensity cutoff during bias estimation
    inline void SetLowIntensityCutoff( double cutoff );
  
    ///Set slices which need to be excluded by default
    inline void SetForceExcludedSlices( vector<int>& force_excluded );

	irtkRealImage GetReconstructedGPU();

    //utility
    ///Save intermediate results
    inline void DebugOn();
    ///Do not save intermediate results
    inline void DebugOff();

    inline void UseAdaptiveRegularisation();
    
    ///Write included/excluded/outside slices
    void Evaluate( int iter );
  
    /// Read Transformations
    void ReadTransformation( char* folder );
  
    //To recover original scaling
    ///Restore slice intensities to their original values
    void RestoreSliceIntensities();
    ///Scale volume to match the slice intensities
    void ScaleVolume();
  
    ///To compare how simulation from the reconstructed volume matches the original stacks
    void SimulateStacks(vector<irtkRealImage>& stacks);

    void SimulateSlices();
  
    ///Puts origin of the image into origin of world coordinates
    void ResetOrigin( irtkGreyImage &image,
                      irtkRigidTransformation& transformation);
  
    ///Packages to volume registrations
    void PackageToVolume( vector<irtkRealImage>& stacks,
                          vector<int> &pack_num,
                          bool evenodd=false,
                          bool half=false,
                          int half_iter=1);
  
    ///Splits stacks into packages
    void SplitImage( irtkRealImage image,
                     int packages,
                     vector<irtkRealImage>& stacks );
    ///Splits stacks into packages and each package into even and odd slices
    void SplitImageEvenOdd( irtkRealImage image,
                            int packages,
                            vector<irtkRealImage>& stacks );
    ///Splits image into top and bottom half roi according to z coordinate
    void HalfImage( irtkRealImage image,
                    vector<irtkRealImage>& stacks );
    ///Splits stacks into packages and each package into even and odd slices and top and bottom roi
    void SplitImageEvenOddHalf( irtkRealImage image,
                                int packages,
                                vector<irtkRealImage>& stacks,
                                int iter=1);

    friend class ParallelStackRegistrations;
    friend class ParallelSliceToVolumeRegistration;
    friend class ParallelCoeffInit;

    friend class ParallelAverage;
    friend class ParallelSliceAverage;
};

inline double irtkReconstructionCuda::G(double x,double s)
{
    return _step*exp(-x*x/(2*s))/(sqrt(6.28*s));
}

inline double irtkReconstructionCuda::M(double m)
{
    return m*_step;
}

inline irtkRealImage irtkReconstructionCuda::GetReconstructed()
{
    return _reconstructed;
}

inline irtkRealImage irtkReconstructionCuda::GetMask()
{
    return _mask;
}

inline void irtkReconstructionCuda::DebugOn()
{
    _debug=true;
    cout<<"Debug mode."<<endl;
}

inline void irtkReconstructionCuda::UseAdaptiveRegularisation()
{
    _adaptive = true;
}

inline void irtkReconstructionCuda::DebugOff()
{
    _debug=false;
}

inline void irtkReconstructionCuda::SetSigma(double sigma)
{
    _sigma_bias=sigma;
}


inline void irtkReconstructionCuda::SpeedupOn()
{
    _quality_factor=1;
}

inline void irtkReconstructionCuda::SpeedupOff()
{
    _quality_factor=2;
}

inline void irtkReconstructionCuda::GlobalBiasCorrectionOn()
{
    _global_bias_correction=true;
}

inline void irtkReconstructionCuda::GlobalBiasCorrectionOff()
{
    _global_bias_correction=false;
}

inline void irtkReconstructionCuda::SetLowIntensityCutoff(double cutoff)
{
    if (cutoff>1) cutoff=1;
    if (cutoff<0) cutoff=0;
    _low_intensity_cutoff = cutoff;
    //cout<<"Setting low intensity cutoff for bias correction to "<<_low_intensity_cutoff<<" of the maximum intensity."<<endl;
}


inline void irtkReconstructionCuda::SetSmoothingParameters(double delta, double lambda)
{
    _delta=delta;
    _lambda=lambda*delta*delta;
    _alpha = 0.05/lambda;
    if (_alpha>1) _alpha= 1;
    cout<<"delta = "<<_delta<<" lambda = "<<lambda<<" alpha = "<<_alpha<<endl;
}

inline void irtkReconstructionCuda::SetForceExcludedSlices(vector<int>& force_excluded)
{
    _force_excluded = force_excluded;  
}

#endif
