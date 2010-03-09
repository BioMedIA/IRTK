/*=========================================================================

  Library   : Image Registration Toolkit (IRTK)
  Module    : $Id: irtkImageRigidRegistration.cc 2 2008-12-23 12:40:14Z dr $
  Copyright : Imperial College, Department of Computing
              Visual Information Processing (VIP), 2008 onwards
  Date      : $Date: 2008-12-23 12:40:14 +0000 (Tue, 23 Dec 2008) $
  Version   : $Revision: 2 $
  Changes   : $Author: dr $

=========================================================================*/

#include <irtkRegistration.h>

#include <irtkHomogeneousTransformationIterator.h>

#include <irtkImageRigidRegistrationGPU.h>

#include <ctime>

#include <cutil_inline.h>
#include <driver_functions.h>


typedef unsigned char uchar;

////////////////////////////////////////////////////////////////////////////////
// declaration, forward
extern "C" int runTest(const int argc,
                        char* data, int2* data_int2, unsigned int len);
extern "C" void initCuda(short *h_target, short *h_source, cudaExtent volumeSize);
extern "C" int runSimpleTexture(dim3 gridSize, dim3 blockSize, int * d_output, float * d_matrix,int imageX,int imageY,int imageZ,int modX, int modY);
extern "C" int runSlowTexture(dim3 gridSize, dim3 blockSize, int * d_output,float * d_matrix,int imageX,int imageY,int imageZ,int modX, int modY);
extern "C" int runRegistrationWithSamples(dim3 gridSize, dim3 blockSize, int * d_output,int * d_samplesArray, float * d_matrix,int imageX,int imageY,int imageZ,int modX, int modY);
extern "C" int runOutputNormSourceVals(dim3 gridSize, dim3 blockSize, float * d_output,float * d_matrix,int imageX,int imageY,int imageZ,int modX, int modY);

void irtkImageRigidRegistrationGPU::GuessParameter()
{
  int i;
  double xsize, ysize, zsize;

  if ((_target == NULL) || (_source == NULL)) {
    cerr << "irtkImageRigidRegistration::GuessParameter: Target and source image not found" << endl;
    exit(1);
  }

  // Default parameters for registration
  _NumberOfLevels     = 3;
  _NumberOfBins       = 64;

  // Default parameters for optimization
  _SimilarityMeasure  = NMI;
  _OptimizationMethod = GradientDescent;
  _Epsilon            = 0.0001;

  // Read target pixel size
  _target->GetPixelSize(&xsize, &ysize, &zsize);

  // Default target parameters
  _TargetBlurring[0]      = GuessResolution(xsize, ysize, zsize) / 2.0;
  _TargetResolution[0][0] = GuessResolution(xsize, ysize, zsize);
  _TargetResolution[0][1] = GuessResolution(xsize, ysize, zsize);
  _TargetResolution[0][2] = GuessResolution(xsize, ysize, zsize);

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
  _SourceResolution[0][0] = GuessResolution(xsize, ysize, zsize);
  _SourceResolution[0][1] = GuessResolution(xsize, ysize, zsize);
  _SourceResolution[0][2] = GuessResolution(xsize, ysize, zsize);

  for (i = 1; i < _NumberOfLevels; i++) {
    _SourceBlurring[i]      = _SourceBlurring[i-1] * 2;
    _SourceResolution[i][0] = _SourceResolution[i-1][0] * 2;
    _SourceResolution[i][1] = _SourceResolution[i-1][1] * 2;
    _SourceResolution[i][2] = _SourceResolution[i-1][2] * 2;
  }

  // Remaining parameters
  for (i = 0; i < _NumberOfLevels; i++) {
    _NumberOfIterations[i] = 20;
    _NumberOfSteps[i]      = 4;
    _LengthOfSteps[i]      = 2 * pow(2.0, i);
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
}



void irtkImageRigidRegistrationGPU::SetupGPU()
{
		cudaSetDevice( cutGetMaxGflopsDeviceId() );
			cutilCheckMsg("Device set failed");

}

void irtkImageRigidRegistrationGPU::Initialize()
{
  // Call base class
  this->irtkImageRegistration::Initialize();

  lastxyz = 0;
  // Invert rigid transformation (to be backwards compatible)
  ((irtkRigidTransformation *)_transformation)->Invert();
  ((irtkRigidTransformation *)_transformation)->UpdateParameter();
  SetupGPU();
}

void irtkImageRigidRegistrationGPU::Finalize()
{
  // Call base class
  this->irtkImageRegistration::Finalize();

  // Invert rigid transformation (to be backwards compatible)
  ((irtkRigidTransformation *)_transformation)->Invert();
  ((irtkRigidTransformation *)_transformation)->UpdateParameter();
}

double irtkImageRigidRegistrationGPU::EvaluateInside(double x, double y, double z, double time,short * ptr)
{
	return 0;
}



void irtkImageRigidRegistrationGPU::LinearMetricAddition(int i, int j, int k,int* _bins,int* _nsamp,
													  int SX, int SY, int SZ,
													  int TX, int TY, int TZ,
													  short * ptr2target,
													  short * ptr2source,
													  float* matrix) {

}

float* irtkImageRigidRegistrationGPU::ConstructGPUMatrix(){
	// Transformation matrix
	matrix = _source->GetWorldToImageMatrix() * ((irtkHomogeneousTransformation *)_transformation)->GetMatrix() * _target->GetImageToWorldMatrix();

	//Construct matrix as a non special data type

	float * outputMatrix = (float *)malloc(sizeof(float) *4 * 4 *4);

	for (int mi = 0; mi < 4; mi++){
		for (int mj = 0; mj < 4; mj++){
			outputMatrix[mi*4 + mj] = (float)matrix(mi,mj);
		}
	}

	for (int mi = 15; mi < 64; mi++){
		outputMatrix[mi] = 0;
	}

	return outputMatrix;
}

void irtkImageRigidRegistrationGPU::checkFloats(){

	  // Invert transformation
	((irtkRigidTransformation *)_transformation)->Invert();

	// Create iterator
	irtkHomogeneousTransformationIterator
		iterator((irtkHomogeneousTransformation *)_transformation);

	// Initialize metric
	_metric->Reset();

	// Pointer to voxels in target image
	ptr2target = _target->GetPointerToVoxels();
	ptr2source = _source->GetPointerToVoxels();
	ptr2targetTest = _target->GetPointerToVoxels();

	size_t imageX = (size_t)_target->GetX();
	size_t imageY = (size_t)_target->GetY();
	size_t imageZ = (size_t)_target->GetZ();

	size_t modX = imageX;
	size_t modY = imageY;

	size_t blockX = 32;
	size_t blockY = 8;


	//Work out necessary size for an x multiple of 32!
	if(imageX % blockX != 0){
		modX = ((imageX / blockX) + 1) * blockX;
	}
	size_t gridX = modX / 32;

	if(imageY % blockY != 0){
		modY = ((imageY / blockY) + 1) * blockY;
	}
	size_t gridY = modY / blockY;


	if(lastxyz !=  _source->GetX()*_source->GetY()*_source->GetZ())
	{
		lastxyz =  _source->GetX()*_source->GetY()*_source->GetZ();

		const cudaExtent volumeSize = make_cudaExtent(imageX, imageY, imageZ);
		initCuda(ptr2target,ptr2source,volumeSize);

		cutilSafeCall( cudaMalloc((void **)&d_sourceFloats, 2*modX*modY*imageZ*sizeof(float)) ); // create output array
		cutilSafeCall( cudaMemset(d_sourceFloats,-1,2*modX*modY*imageZ*sizeof(float)));
	}

	h_sourceFloats = (float*) malloc(modX*modY*imageZ*sizeof(float));


	float* transformationMatrix = ConstructGPUMatrix();


	cutilSafeCall( cudaMalloc((void **)&d_matrix, 64*sizeof(float)) ); // the matrix
	cudaMemcpy(d_matrix,transformationMatrix,64*sizeof(float),cudaMemcpyHostToDevice);


	const dim3 gridSize(gridX * imageZ, gridY, 1);
	const dim3 blockSize(blockX,blockY,1);

	runOutputNormSourceVals(gridSize,blockSize,d_sourceFloats,d_matrix,imageX,imageY,imageZ,modX,modY);

	cutilSafeCall( cudaThreadSynchronize() );
	cutilCheckMsg("render_kernel failed");

	cutilSafeCall( cudaMemcpy( h_sourceFloats, d_sourceFloats,modX*modY*imageZ*sizeof(int), cudaMemcpyDeviceToHost) );

	for(int k = 0; k < imageZ; k++){
		for(int j = 0; j < imageY; j++){
			for (int i = 0; i < imageX; i++){
				if(h_sourceFloats[i + j * imageX + k * imageX*imageY] > -1){
					printf("x: %d y: %d z: %d val: %f\n",i,j,k,h_sourceFloats[i + j * imageX + k * imageX*imageY]);
				}
			}
		}
	}

	cutilSafeCall(cudaFree(d_matrix));
	cutilSafeCall(cudaFree(d_sourceFloats));
	free(h_sourceFloats);
	free(transformationMatrix);
}

double irtkImageRigidRegistrationGPU::Evaluate()
{

	//checkFloats();

		if (_clockAcc < 0) _clockAcc = 0;
		std::clock_t start = std::clock();

	int **_bins;
	int _nsamp = 0;
	float _x,_y,_z;

  // Print debugging information
  this->Debug("irtkImageRigidRegistration::Evaluate");

  // Invert transformation
  ((irtkRigidTransformation *)_transformation)->Invert();

  // Create iterator
  irtkHomogeneousTransformationIterator
  iterator((irtkHomogeneousTransformation *)_transformation);

  irtkGreyPixel target_min, target_max, target_nbins;
  irtkGreyPixel source_min, source_max, source_nbins;

  // Initialize metric
  _metric->Reset();
  
  _bins = _metric->GetBins();


  // Pointer to voxels in target image
  ptr2target = _target->GetPointerToVoxels();
  ptr2source = _source->GetPointerToVoxels();
  ptr2targetTest = _target->GetPointerToVoxels();

  if (_transformation == NULL) {
    cout << "irtkHomogeneousTransformationIterator::Initialize(): Transformation has not been set." << endl;
    exit(1);
  }

  	size_t imageX = (size_t)_target->GetX();
	size_t imageY = (size_t)_target->GetY();
	size_t imageZ = (size_t)_target->GetZ();

	size_t modX = imageX;
	size_t modY = imageY;

	size_t blockX = 32;
	size_t blockY = 8;


	//Work out necessary size for an x multiple of 32!
	if(imageX % blockX != 0){
		modX = ((imageX / blockX) + 1) * blockX;
	}
	size_t gridX = modX / 32;

	if(imageY % blockY != 0){
		modY = ((imageY / blockY) + 1) * blockY;
	}
	size_t gridY = modY / blockY;

		
	if(lastxyz !=  _source->GetX()*_source->GetY()*_source->GetZ())
	{
		lastxyz =  _source->GetX()*_source->GetY()*_source->GetZ();

		const cudaExtent volumeSize = make_cudaExtent(imageX, imageY, imageZ);
			initCuda(ptr2target,ptr2source,volumeSize);

		cutilSafeCall( cudaMalloc((void **)&d_output, 2*modX*modY*imageZ*sizeof(int)) ); // create output array
		
	}
	cutilSafeCall( cudaMemset(d_output,0,2*modX*modY*imageZ*sizeof(float)));

	cutilSafeCall( cudaMalloc((void **)&d_samplesArray, 1000*sizeof(int)) );

	d_results = (int*) malloc(modX*modY*imageZ*sizeof(int));
	h_samplesArray = (int*) malloc(1000*sizeof(int));


	float* transformationMatrix = ConstructGPUMatrix();

	cutilSafeCall( cudaMalloc((void **)&d_matrix, 64*sizeof(float)) ); // the matrix
	cudaMemcpy(d_matrix,transformationMatrix,64*sizeof(float),cudaMemcpyHostToDevice);

	
		const dim3 gridSize(gridX * imageZ, gridY, 1);
		const dim3 blockSize(blockX,blockY,1);

		int result = runSlowTexture(gridSize,blockSize,d_output,d_matrix,imageX,imageY,imageZ,modX,modY);

	cutilSafeCall( cudaThreadSynchronize() );
	cutilCheckMsg("render_kernel failed");

	cutilSafeCall( cudaMemcpy( d_results, d_output,modX*modY*imageZ*sizeof(int), cudaMemcpyDeviceToHost) );
	//cutilSafeCall( cudaMemcpy( h_samplesArray, d_samplesArray,1000*sizeof(int), cudaMemcpyDeviceToHost) );
	//cout << "max(short): " << numeric_limits<short>::max() << endl;
/*	int testingSamp = 0;

		for (int nsampi = 0; nsampi < modX*modY*imageZ; nsampi++){
		  if(d_results[nsampi] > -1){
			  //printf("result: %d\n",d_results[nsampi]);
			(*_bins)[d_results[nsampi]] += 1;
			_nsamp += 1;
			//testingSamp += d_results[nsampi];
		  }
	  }
	 for (int nsampi = 0; nsampi < 1000; nsampi++){
		  if(h_samplesArray[nsampi] > 0){
			testingSamp += h_samplesArray[nsampi];
		  }
	  }
*/
	for (int nsampi = 0; nsampi < 63*64; nsampi++){
		if(d_results[nsampi] > -1){
			//printf("result: %d\n",d_results[nsampi]);
			(*_bins)[nsampi] = d_results[nsampi];
			_nsamp += d_results[nsampi];
			//testingSamp += d_results[nsampi];
		}
	}

	  _metric->PutBins(_bins);
	  _metric->PutNSamp(_nsamp);

	//cutilSafeCall(cudaFree(d_output));
	cutilSafeCall(cudaFree(d_matrix));
	cutilSafeCall(cudaFree(d_samplesArray));
	free(h_samplesArray);
	free(d_results);
	free(transformationMatrix);
	//cudaThreadExit();

std::cout<< ( ( std::clock() - start ) / (double)CLOCKS_PER_SEC ) <<'\n';
  _clockAcc += (( std::clock() - start ) / (double)CLOCKS_PER_SEC);
	std::cout<< ( _clockAcc ) <<'\n';

	return _metric->Evaluate();


}
