#ifndef _IRTKBEP_H

#define _IRTKBEP_H

#ifdef HAS_VTK

#include <irtkRegistration.h>
#include <irtkGaussianBlurring.h>
#include <irtkSegmentationFunction.h>


class irtkBep : public irtkObject
{
private:
	vtkPolyData *bepsurface;
	vtkPolyData *datasurface;
	char* outputfilename;

public:

	irtkBep();

	virtual void SetInput(vtkPolyData *bepsurface, vtkPolyData *datasurface);

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