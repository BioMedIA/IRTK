/*=========================================================================

  Library   : Image Registration Toolkit (IRTK)
  Module    : $Id: irtkImageRigidRegistration.h 2 2008-12-23 12:40:14Z dr $
  Copyright : Imperial College, Department of Computing
              Visual Information Processing (VIP), 2008 onwards
  Date      : $Date: 2008-12-23 12:40:14 +0000 (Tue, 23 Dec 2008) $
  Version   : $Revision: 2 $
  Changes   : $Author: dr $

=========================================================================*/

#ifndef _IRTKIMAGERIGIDREGISTRATIONGPU_H

#define _IRTKIMAGERIGIDREGISTRATIONGPU_H


/**
 * Filter for rigid registration based on voxel similarity measures.
 *
 * This class implements a registration filter for the rigid registration of
 * two images. The basic algorithm is described in Studholme, Medical Image
 * Analysis, Vol. 1, No. 2, 1996.
 *
 */

class irtkImageRigidRegistrationGPU : public irtkImageRegistration
{

protected:

	int* outBI;
	int* outNI;
	int* outIMG;
	int lastxyz;

	int *d_output; //output array//
	int *d_results; //host companion//


	int *d_samplesArray; //output array//
	int *h_samplesArray; //host companion//

	float *d_matrix; //output array//

	//For use in checking voxel results
	float *d_sourceFloats;
	float *h_sourceFloats;

	  // Pointer to reference data
  irtkGreyPixel *ptr2target;
  irtkGreyPixel *ptr2targetTest;
  irtkGreyPixel *ptr2source;

	irtkMatrix matrix;

	double _clockAcc;

  /// Evaluate the similarity measure for a given transformation.
  virtual double Evaluate();

  /// Perform interpolation in class flatly
  virtual double EvaluateInside(double x, double y, double z, double t,short * ptr);

  virtual void LinearMetricAddition(int i, int j, int k,int* _bins,int *_nsamp,
													  int SX, int SY, int SZ,
													  int TX, int TY, int TZ,
													  short * ptr2target,
													  short * ptr2source,
													  float* matrix) ;

  /// Initial set up for the GPU parallesiation
  virtual void SetupGPU();

  /// Initial set up for the registration
  virtual void Initialize();

  /// Final set up for the registration
  virtual void Finalize();

  virtual void checkFloats();

  virtual float* ConstructGPUMatrix();

public:

  /** Sets the output for the registration filter. The output must be a rigid
   *  transformation. The current parameters of the rigid transformation are
   *  used as initial guess for the rigid registration. After execution of the
   *  filter the parameters of the rigid transformation are updated with the
   *  optimal transformation parameters.
   */
  virtual void SetOutput(irtkTransformation *);

  /// Returns the name of the class
  virtual const char *NameOfClass();

  /// Print information about the progress of the registration
  virtual void Print();

  /// Guess parameters
  virtual void GuessParameter();
};

inline void irtkImageRigidRegistrationGPU::SetOutput(irtkTransformation *transformation)
{
  if (strcmp(transformation->NameOfClass(), "irtkRigidTransformation") != 0) {
    cerr << "irtkImageRigidRegistration::SetOutput: Transformation must be rigid"
         << endl;
    exit(0);
  }
  _transformation = transformation;
}

inline const char *irtkImageRigidRegistrationGPU::NameOfClass()
{
  return "irtkImageRigidRegistrationGPU";
}

inline void irtkImageRigidRegistrationGPU::Print()
{
  _transformation->Print();
}

#endif
