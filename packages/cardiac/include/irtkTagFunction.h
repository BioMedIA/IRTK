#ifndef _IRTKTAGFUNCTION_H

#define _IRTKTAGFUNCTION_H

#include <irtkRegistration.h>

#ifdef HAS_OPENCV

#include <highgui.h>

#endif

class irtkTagFunction : public irtkObject
{
protected:
	// interest region image
	irtkGreyImage *_interest;
	// interest theshold
	irtkRealImage *_threshold;
	// pointset
	irtkPointSet _p1,_p2;
	// interest region affinetransformation
	irtkAffineTransformation* _interestaffine;
	// target and source images
	irtkGreyImage *_target,*_source;
	// transformation and registration;
	irtkAffineTransformation *_transformation;
	// registration
	irtkImageAffineRegistration *_registration;
    //switch sampled or not
	int sampled;
	irtkGreyImage* _maxmap;
	// exemption image for sampling tag
	irtkGreyImage* _exemption;

	virtual void Sample(void);
    virtual void Initialize(void);
	virtual void Finalize(void);
	virtual void InitializeTransformation(irtkGreyImage& target, irtkGreyImage &source);
	virtual void FinalizeTransformation(irtkGreyImage &target, irtkGreyImage &source);

public:
	irtkTagFunction(void);
	~irtkTagFunction(void);
    // SetInput of the function
	virtual void SetInput(irtkGreyImage* interest, irtkRealImage* threshold, irtkGreyImage* target, irtkGreyImage* source, irtkAffineTransformation* interestaffine);
	// SetOutput pointset
	virtual irtkPointSet GetOutput(void);
	// toggle 0 no center toggle 2 cente
	virtual void Track(int toggle = 2, int ds = 2);
	virtual void SetPointSet(irtkPointSet& p1);
	virtual irtkPointSet GetPointSet(void);

};

#ifdef HAS_OPENCV

void cvShiftDFT(CvArr * src_arr, CvArr * dst_arr );

#endif
#endif