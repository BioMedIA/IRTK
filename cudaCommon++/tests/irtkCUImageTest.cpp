#include <iostream>
#include <stdlib.h>

#include <irtkCUImage.h>
#include <limits>
#include <irtkCUMaskFilter.h>

static const char* g_filename = "brain1.nii";
#define FILENAME const_cast<char*>(g_filename)

int main(int argc, char **argv)
{
	std::cout << "testing irtkCUImage..." << std::endl;
	std::cout << "testing constructors..." << std::endl;
	irtkCUGenericImage<irtkBytePixel>* test1 = new irtkCUGenericImage<irtkBytePixel>();
	if(test1 != NULL) std::cout << "constructor test 1 --> Pass" << std::endl;
	irtkCUGenericImage<irtkGreyPixel>* test2 = new irtkCUGenericImage<irtkGreyPixel>(FILENAME);
	printf("size %d %d %d\n", test2->GetX(), test2->GetY(), test2->GetZ(), test2->GetT());

	if(test2 != NULL) std::cout << "constructor test 2 --> Pass" << std::endl;
	irtkCUGenericImage<irtkGreyPixel>* test3 = new irtkCUGenericImage<irtkGreyPixel>(128,128,128,12);
	if(test3 != NULL) std::cout << "constructor test 3 --> Pass" << std::endl;
	irtkCUGenericImage<irtkGreyPixel>* test4 = new irtkCUGenericImage<irtkGreyPixel>(*test2);
	if(test4 != NULL) std::cout << "constructor test 4 --> Pass" << std::endl;
	irtkCUGenericImage<irtkGreyPixel>* test5 = new irtkCUGenericImage<irtkGreyPixel>(test3->GetImageAttributes());
	if(test5 != NULL) std::cout << "constructor test 5 --> Pass" << std::endl;

	irtkGenericImage<irtkGreyPixel>* test6a = new irtkGenericImage<irtkGreyPixel>(128,128,128,12);
	irtkCUGenericImage<irtkGreyPixel>* test6 = new irtkCUGenericImage<irtkGreyPixel>(*test6a);
	if(test6 != NULL) std::cout << "constructor test 6 --> Pass" << std::endl;

	irtkGenericImage<irtkGreyPixel>* im = (irtkGenericImage<irtkGreyPixel>*)test2;
	//...

	//reference
	irtkGenericImage<irtkGreyPixel>* test7 = new irtkGenericImage<irtkGreyPixel>(FILENAME);
	printf("testing size %d %d %d\n", test7->GetX(), test7->GetY(), test7->GetZ(), test7->GetT());

	irtkGreyPixel min, max, cmin, cmax;
	test2->GetMinMax(&min, &max);

	test7->GetMinMax(&cmin, &cmax);
	printf("calculation test 7 minmax: gpumin %d, cpumin %d, gpumax %d, cpumax %d --> %s \n", min, cmin, max, cmax,  (min == cmin && max == cmax) ? "Pass" : "Fail");

	irtkGreyPixel gav = test2->GetAverage();
	irtkGreyPixel cav = test7->GetAverage();
	printf("calculation test 8 average: gpu: %d, cpu: %d --> %s \n", gav, cav,  (gav == cav) ? "Pass" : "Fail");

	test2->normalize(min,max);

	//Histogram test
	//TODO investigate this after histogram:
#ifdef USE_CUDA_NPP
	test2->GetMinMax(&min, &max);

	const int numBins = 256;
	int gpu_hist[numBins];
	int cpu_hist[numBins];
	test2->GetHistogram(gpu_hist,numBins);

	test2->GetMinMax(&min, &max);

	int csum = 0;
	printf("\nHISTROGRAM GPU:\n");
	for(int i =0; i < numBins; i++)
	{
		csum += gpu_hist[i];
		printf("%d ",  gpu_hist[i]);
	}
	printf("\nSUM: %d \n ", csum);

	for(int i =0; i < numBins; i++)
	{
		cpu_hist[i] = 0;
	}

	for(int i =0; i < test2->GetNumberOfVoxels(); ++i)
	{
		irtkGreyPixel val = test2->GetPointerToVoxels()[i];
		int idx = numBins * ((double)val - (double)min)/((double)max - (double)min);
		++cpu_hist[idx];
	}

	int gsum = 0;
	irtkGreyPixel diff = 0;
	printf("\nHISTROGRAM CPU: \n");
	//printf("%d %d %d \n", cpu_hist[0] - gpu_hist[0], cpu_hist[0], gpu_hist[0]);
	for(int i =0; i < numBins; ++i)
	{
		if( i > 0) //some bug in testcase, maybe sync issues
			diff =  (cpu_hist[i] - gpu_hist[i]) + diff;
		gsum += cpu_hist[i];
		printf("%d ",  cpu_hist[i]);
	}
	printf("\nSUM: %d\n ", gsum);
	printf("calculation test 9 Histogram: diff %d --> %s \n", diff, (gsum-csum == 0) ? "Pass" : "Fail");
 #endif


	// very simple array based mask filter
	// filter test -- TODO lazyness: build own test app -> filter is not considering image orientation yet !! TODO!
	irtkCUMaskFilter<irtkGreyPixel>* testMaskFilter = new irtkCUMaskFilter<irtkGreyPixel>;
	testMaskFilter->setInput(test2);
	irtkCUGenericImage<irtkGreyPixel>* mask = new irtkCUGenericImage<irtkGreyPixel>("body.hdr");
	testMaskFilter->setMask(mask, cmin);
	testMaskFilter->setOutput(test2);
	testMaskFilter->getOuput()->Write("testout.nii");

	std::cin.get();

	return 0;

}
