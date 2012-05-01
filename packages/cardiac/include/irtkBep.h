#ifndef _IRTKBEP_H

#define _IRTKBEP_H

#include <irtkRegistration.h>
#include <irtkGaussianBlurring.h>
#include <irtkSegmentationFunction.h>

#ifdef HAS_VTK

class irtkBep : public irtkObject
{
private:
	vtkPolyData *bepsurface;
	vtkPolyData *datasurface;
	char* outputfilename;
    int numberofsegments;

public:

	irtkBep();

	virtual void SetInput(vtkPolyData *bepsurface, vtkPolyData *datasurface, int);

	virtual void SetOutput(char *outputfilename);

	virtual void Initialize();

	virtual void Finalize();

    /// Evaluate
    virtual void Bullseyeplot();

    /// Segment infarction
	virtual void EvaluateInfarction(irtkRealImage *, irtkRealImage **, int number = 1,double timeweight = 1,int cutmode = 3,int connectmode = 0, double regionweight = 4);

};

#endif

#endif
