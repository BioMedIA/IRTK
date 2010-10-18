/*=========================================================================

  Library   : Image Registration Toolkit (IRTK)
  Module    : $Id$
  Copyright : Imperial College, Department of Computing
              Visual Information Processing (VIP), 2009 onwards
  Date      : $Date$
  Version   : $Revision$
  Changes   : $Author$

=========================================================================*/

#ifndef _IRTKMODELREGISTRATION_H

#define _IRTKMODELREGISTRATION_H

#ifdef HAS_VTK

#include <irtkModelSimilarityMetric.h>

/**
 * Generic for registering a model (surface) to an image.
 *
 * This class implements a registration filter which takes a model and an 
 * images and calculates the transformation which maps the model into the 
 * coordinate system of the image. This is the abstract base class which 
 * defines a common interface for arbitrary registration filters. Each 
 * derived class has to implement all abstract member functions.
 *
 */

class irtkModelRegistration : public irtkRegistration
{

  /// Interface to input file stream
  friend istream& operator>> (istream&, irtkModelRegistration*);

  /// Interface to output file stream
  friend ostream& operator<< (ostream&, const irtkModelRegistration*);

protected:

  /// Image
  irtkGreyImage *_image;

  /// Model
  vtkPolyData *_model;
  
  /// Output
  irtkTransformation *_transformation;

  /// Metric
  irtkModelSimilarityMetric *_metric;

  /// Optimizer
  irtkOptimizer *_optimizer;

  /// Blurring of image (in mm)
  double _ImageBlurring[MAX_NO_RESOLUTIONS];

  /// Resolution of image (in mm)
  double _ImageResolution[MAX_NO_RESOLUTIONS][3];

  /// Number of step sizes
  int    _NumberOfSteps[MAX_NO_RESOLUTIONS];

  /// Length of steps
  double _LengthOfSteps[MAX_NO_RESOLUTIONS];

  /// Max. number of iterations per step size
  int    _NumberOfIterations[MAX_NO_RESOLUTIONS];

  /// Number of levels of multiresolution pyramid
  int    _NumberOfLevels;

  /// Optimization method for registration
  irtkOptimizationMethod _OptimizationMethod;

  /// Convergence parameter for optimization based on change in similarity.
  double _Epsilon;

  /// Convergence parameter for optimization based on change in the transformation.
  double _Delta[MAX_NO_RESOLUTIONS];

  /// Debugging flag
  int    _DebugFlag;

  /// Initial set up for the registration
  virtual void Initialize();

  /// Final set up for the registration
  virtual void Finalize();

  /// Initial set up for the registration at a multiresolution level
  virtual void Initialize(int);

  /// Final set up for the registration at a multiresolution level
  virtual void Finalize(int);

public:

  /// Constructor
  irtkModelRegistration();

  /// Destructor
  virtual ~irtkModelRegistration();

  /// Sets input for the registration filter
  virtual void SetInput (vtkPolyData *, irtkGreyImage *);

  /// Sets output for the registration filter
  virtual void SetOutput(irtkTransformation *) = 0;

  /// Runs the registration filter
  virtual void Run();

  /** Evaluates the similarity metric. This function evaluates the similarity
   *  metric of the registration by looping over the model points, transforming
   *  the model points into the image and computing the fit of the model in the
   *  image. This function returns the value of the similarity measure using 
   *  Similarity().
   */
  virtual double Evaluate() = 0;

  /// Evaluates the gradient of the similarity metric. This function
  virtual double EvaluateGradient(float, float *);

  /// Returns the name of the class
  virtual const char *NameOfClass() = 0;

  /// Prints debugging messages if debugging is enabled
  virtual void Debug(string);

  /// Prints information about the progress of the registration
  virtual void Print() = 0;

  /// Guess parameters
  virtual void GuessParameter() = 0;

  /// Read registration parameters from file
  virtual void Read (char *);

  /// Parse parameter line
  virtual bool Read(char *, char *, int &);

  /// Write registration parameters to file
  virtual void Write(char *);

  /// Write parameters to stream
  virtual void Write(ostream &);

  // Access parameters
  virtual SetMacro(DebugFlag, int);
  virtual GetMacro(DebugFlag, int);
  virtual SetMacro(OptimizationMethod, irtkOptimizationMethod);
  virtual GetMacro(OptimizationMethod, irtkOptimizationMethod);

};

inline void irtkModelRegistration::SetInput(vtkPolyData *model, irtkGreyImage *image)
{
  _model = model;
  _image = image;
}

inline void irtkModelRegistration::Debug(string message)
{
  if (_DebugFlag == true) cout << message << endl;
}

#include <irtkModelRigidRegistration.h>
#include <irtkModelAffineRegistration.h>
#include <irtkModelFreeFormRegistration.h>

#endif

#endif
