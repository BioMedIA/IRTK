#include "reconstruction.h"

void _reconstruct(
                 // input stacks
                 float* img,
                 double* pixelSize,
                 double* xAxis,
                 double* yAxis,
                 double* zAxis,
                 double* origin,
                 int* dim,

                 // number of stacks
                 int n,

                 // stack ids: which stack each slice
                 // comes from
                 int* _stack_ids,

                 // number of reconstruction iterations to run
                 int iterations,

                 // initial transformations
                 double* tx,
                 double* ty,
                 double* tz,
                 double* rx,
                 double* ry,
                 double* rz,

                 // slice thickness
                 double* _thickness,

                 // mask (header same as template)
                 float* mask_img,
                 
                 // output: reconstructed image                 
                 float* reconstructed_img,
                 double* reconstructed_pixelSize,
                 double* reconstructed_xAxis,
                 double* reconstructed_yAxis,
                 double* reconstructed_zAxis,
                 double* reconstructed_origin,
                 int* reconstructed_dim ) {

    /// Slices
    std::vector< irtkGenericImage<float> > slices;
    pyList2irtkVector( slices,
                       img,
                       pixelSize,
                       xAxis,
                       yAxis,
                       zAxis,
                       origin,
                       dim,
                       n );

    /// Stack transformation
    std::vector<irtkRigidTransformation> slice_transformations;
    pyList2rigidVector( slice_transformations,
                        tx, ty, tz,
                        rx, ry, rz,
                        n );

    // convert thickness and stack_ids to std::vector
    std::vector<double> thickness(n);
    std::vector<int> stack_ids(n);
    for ( int i = 0; i < n; i++ ) {
        thickness[i] = _thickness[i];
        stack_ids[i] = _stack_ids[i];
    }

    // template
    irtkGenericImage<float> reconstructed;
    py2irtk<float>( reconstructed,
                    reconstructed_img,
                    reconstructed_pixelSize,
                    reconstructed_xAxis,
                    reconstructed_yAxis,
                    reconstructed_zAxis,
                    reconstructed_origin,
                    reconstructed_dim );

    // mask
    irtkGenericImage<float> mask;
    py2irtk<float>( mask,
                    mask_img,
                    reconstructed_pixelSize,
                    reconstructed_xAxis,
                    reconstructed_yAxis,
                    reconstructed_zAxis,
                    reconstructed_origin,
                    reconstructed_dim );

    std::cout << "now creating reconstruction object\n";
    
    //Create reconstruction object
    irtkReconstruction reconstruction;
    int levels = 3; // Number of resolution levels
    double sigma = 20; // Stdev for bias field
    double lambda = 0.02;
    double delta = 150;
    double lastIterLambda = 0.01;
    int rec_iterations;

    //Set debug
    reconstruction.DebugOn();
    
    //Set low intensity cutoff for bias estimation
    reconstruction.SetLowIntensityCutoff(0.01);

    //Set global bias correction flag
    reconstruction.GlobalBiasCorrectionOff();    

    //Rescale intensities of the stacks to have the same average
    // reconstruction.MatchStackIntensities(stacks,stack_transformations,700);  
  
    reconstruction.SetSlicesAndTransformations( slices,
                                                slice_transformations,
                                                stack_ids,
                                                thickness );

    std::cout << "slices have been read\n";

    reconstructed.Print();

    reconstruction.SetReconstructed( reconstructed );

    reconstruction.SetMask( &mask, 0 );

    std::cout << " ready to go...\n";

    //return;

    // //Initialise data structures for EM
    // reconstruction.InitializeEM();    

    // //Calculate matrix of transformation between voxels of slices and volume
    // reconstruction.CoeffInit();
    
    // //Initialize reconstructed image with Gaussian weighted reconstruction
    // reconstruction.GaussianReconstruction();
        
    //Initialise data structures for EM
    reconstruction.InitializeEM();

    //
    // Interleaved registration-reconstruction iterations
    //
    for ( int iter = 0; iter < iterations; iter++ ) {
        std::cout<<"Iteration "<<iter<<".\n";


        //perform slice-to-volume registrations - skip the first iteration 
        if ( iter > 0 )      
            reconstruction.SliceToVolumeRegistration(); // UPSAMPLING
   
        //Set smoothing parameters 
        //amount of smoothing (given by lambda) is decreased with improving alignment
        //delta (to determine edges) stays constant throughout
        if ( iter == iterations-1 )
            reconstruction.SetSmoothingParameters(delta,lastIterLambda);
        else {
            double l = lambda;
            for ( int i = 0; i < levels; i++ ) {
                if ( iter == iterations*(levels-i-1)/levels )
                    reconstruction.SetSmoothingParameters(delta, l);
                l*=2;
            }
        }
    
        //Use faster reconstruction during iterations and slower for final reconstruction
        if ( iter < iterations-1 )
            reconstruction.SpeedupOn();
        else 
            reconstruction.SpeedupOff();
    
        //Initialise values of weights, scales and bias fields
        reconstruction.InitializeEMValues();
    
        //Calculate matrix of transformation between voxels of slices and volume
        reconstruction.CoeffInit();
    
        //Initialize reconstructed image with Gaussian weighted reconstruction
        reconstruction.GaussianReconstruction();
    
        //Initialize robust statistics parameters
        reconstruction.InitializeRobustStatistics();
    
        //EStep
        reconstruction.EStep();

        //number of reconstruction iterations
        if ( iter == iterations-1 ) 
            rec_iterations = 5;      
        else 
        rec_iterations = 3;//5;
        
        //reconstruction iterations
        for ( int i = 0; i < rec_iterations; i++ ) {
            cout<<endl<<"  Reconstruction iteration "<<i<<". "<<endl;
      
            // // if (intensity_matching) {
            //     //calculate bias fields
            //     if (sigma>0)
            //         reconstruction.Bias();
            //     //calculate scales
            //     reconstruction.Scale();
            // }
      
            //MStep and update reconstructed volume
            reconstruction.Superresolution( i+1 );
            // if ( sigma > 0 )
            //     reconstruction.NormaliseBias(i);

            cout << "starting MStep\n";
            reconstruction.MStep(i+1);
            cout << "finished MStep\n";
            
            //E-step
            cout << "starting EStep\n";
            reconstruction.EStep();
            cout << "finished EStep\n";
            
        } //end of reconstruction iterations
    
        // Mask reconstructed image to ROI given by the mask
        // reconstruction.MaskVolume();
   
    }// end of interleaved registration-reconstruction iterations

    // final result
    //reconstruction.RestoreSliceIntensities();
    reconstruction.ScaleVolume();
    //reconstructed.Print();

    //irtkGenericImage<float> reconstructed2;
   
    reconstructed = reconstruction.GetReconstructed();
    //reconstructed.Write( "toto.nii");

    // reconstructed2.Print();
    
    irtk2py<float>( reconstructed,
                    reconstructed_img,
                    reconstructed_pixelSize,
                    reconstructed_xAxis,
                    reconstructed_yAxis,
                    reconstructed_zAxis,
                    reconstructed_origin,
                    reconstructed_dim );

    reconstruction.GetTransformations( slice_transformations );

    rigidVector2pyList( slice_transformations,
                        tx, ty, tz,
                        rx, ry, rz );

}
