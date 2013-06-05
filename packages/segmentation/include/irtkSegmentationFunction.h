#ifndef _IRTKSEGMENTATIONFUNCTION_H

#define _IRTKSEGMENTATIONFUNCTION_H

#include <irtkRegistration.h>

#include <irtkEMClassification.h>

#include <irtkImageGraphCut.h>

#include <irtkMultiImageGraphCut.h>

#include <irtkGradientImageFilter.h>

#include <irtkPatchMatch.h>

#include <irtkPatchMatchSegmentation.h>

#ifdef HAS_OPENCV

#include <highgui.h>

#endif

class irtkSegmentationFunction : public irtkObject
{


public:

	// Evaluate Heart Center Location
#ifdef HAS_OPENCV
	/// detect object with threshold
	virtual irtkAffineTransformation DetectObject(irtkRealImage *, irtkGreyImage *, CvHaarClassifierCascade *, int, double fscale = 1.1, int size = 50);
	/// detect object without threshold
	virtual irtkAffineTransformation DetectObject(irtkGreyImage *, CvHaarClassifierCascade *, int, double fscale = 1.1, int size = 50);
#else
	virtual irtkAffineTransformation DetectObject(irtkRealImage *, irtkGreyImage *, int, double fscale = 1.1);
#endif
  virtual void GenerateBox(irtkAffineTransformation& interestregion, irtkGreyImage* interest, irtkPoint& ip1, irtkPoint& ip2, irtkGreyImage* utarget);

  virtual void EvaluateThreshold(irtkRealImage *, irtkGreyImage *, int padding = 0, int n = 3);

  virtual void EvaluateGraphCut(irtkGreyImage *, irtkGreyImage **,irtkRealImage **, int numberofinput, int numberofatlas, double timeweight = 1,int cutmode = 2, int connectmode = 0, double regionweight = 0.9, int padding = 0, int *c = NULL);

  virtual void EvaluateGraphCut(irtkRealImage **, irtkRealImage **,irtkRealImage **, int numberofinput, int numberofatlas, double timeweight = 1,int cutmode = 2, int connectmode = 0, double regionweight = 0.9, int padding = 0, int *c = NULL);

  virtual void EvaluateGraphCut(irtkGreyImage *, irtkGreyImage **, int number = 1,double timeweight = 1,int cutmode = 2,int connectmode = 0, double regionweight = 0.9,int padding = 0, int n = 2);

  virtual void RegionGrow(irtkPoint& center, irtkGreyImage* target, irtkRealImage* threshold, int minT, int value = 2);

  virtual void FixThreshold(irtkRealImage* threshold, irtkRealImage* region);

  virtual void Normalize(irtkRealImage* target, irtkRealImage* source, int n = 3);

  virtual void DetectBackGround(irtkRealImage* target, irtkRealImage* source, int n = 3);

};

#ifdef HAS_OPENCV

void on_mouse( int event, int x, int y, int flags, void* param );

#endif

#endif