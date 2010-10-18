/*=========================================================================

  Library   : Image Registration Toolkit (IRTK)
  Module    : $Id$
  Copyright : Imperial College, Department of Computing
              Visual Information Processing (VIP), 2008 onwards
  Date      : $Date$
  Version   : $Revision$
  Changes   : $Author$

=========================================================================*/

#ifndef _IRTKDEMONSREGISTRATION_H

#define _IRTKDEMONSREGISTRATION_H

#define MAX_NO_RESOLUTIONS 10

#include <irtkImage.h>

#include <irtkResampling.h>

#include <irtkImageFunction.h>

#include <irtkTransformation.h>

#include <irtkUtil.h>

/**
 * Generic for image registration based on voxel similarity measures.
 *
 * This class implements a registration filter which takes two input images
 * and calculates the transformation which maps the second image (denotes as
 * source image) into the coordinate system of the first image (denotes as
 * target image).  This is the abstract base class which defines a common
 * interface for arbitrary registration filters. Each derived class has to
 * implement all abstract member functions.
 *
 */

class irtkDemonsRegistration
{

protected:

  /** First input image. This image is denoted as target image and its
   *  coordinate system defines the frame of reference for the registration.
   */
  irtkRealImage *_target;

  /** Second input image. This image is denoted as source image. The goal of
   *  the registration is to find the transformation which maps the source
   *  image into the coordinate system of the target image.
   */
  irtkRealImage *_source;

  /** Temporary target image
   */
  irtkRealImage _targetTmp;

  /** Temporary source image
   */
  irtkRealImage _sourceTmp;

  /** Gradient of the source image.
   */
  irtkRealImage _sourceGradient;

  /** Gradient of the target image.
   */
  irtkRealImage _targetGradient;

  /** Local displacement field.
   */
  irtkGenericImage<double> _local1, _local2;

  /// Output
  irtkMultiLevelFreeFormTransformation *_transformation1;
  irtkMultiLevelFreeFormTransformation *_transformation2;

  irtkMultiLevelFreeFormTransformation *_transhist1;
  irtkMultiLevelFreeFormTransformation *_transhist2;

  /// Linear free-form transformation
  irtkLinearFreeFormTransformation *_ffd1;
  irtkLinearFreeFormTransformation *_ffd2;

  /// Transformation filters
  irtkImageTransformation _imagetransformation1;
  irtkImageTransformation _imagetransformation2;

  /// Interpolator
  irtkInterpolateImageFunction *_interpolator1;
  irtkInterpolateImageFunction *_interpolator2;

  /// Blurring of target image (in mm)
  double _TargetBlurring[MAX_NO_RESOLUTIONS];

  /// Resolution of target image (in mm)
  double _TargetResolution[MAX_NO_RESOLUTIONS][3];

  /// Blurring of source image (in mm)
  double _SourceBlurring[MAX_NO_RESOLUTIONS];

  /// Resolution of source image (in mm)
  double _SourceResolution[MAX_NO_RESOLUTIONS][3];

  /// Padding value of target image
  short  _TargetPadding;

  /// Padding value of source image
  short  _SourcePadding;

  /// Smoothing of deformation field at every iteration
  double _Smoothing[MAX_NO_RESOLUTIONS];

  /// Number of levels of multiresolution pyramid
  int    _NumberOfLevels;

  /// Max. number of iterations per step size
  int    _NumberOfIterations;

  /// Step size
  double _StepSize;

  /// Epsilon
  double _Epsilon;

  /// Regridding
  int _Regridding;

  /// Symmetric demons?
  bool _Symmetric;

  /// Debugging flag
  int    _DebugFlag;

  /// Source image domain which can be interpolated fast
  double _source_x1, _source_y1, _source_z1;
  double _source_x2, _source_y2, _source_z2;

  /// Target image domain which can be interpolated fast
  double _target_x1, _target_y1, _target_z1;
  double _target_x2, _target_y2, _target_z2;

  /// Initial set up for the registration
  virtual void Initialize();

  /// Initial set up for specific level
  virtual void Initialize(int);

  /// Final set up for the registration
  virtual void Finalize();

  /// Final set up for specific level
  virtual void Finalize(int);

  /// Compute force
  virtual double Force();
  virtual double Force2();

  /// Compute smoothing
  virtual void Smooth(double);

  /// Compute update
  virtual void Update();

public:

  /// Constructor
  irtkDemonsRegistration();

  /// Destructor
  virtual ~irtkDemonsRegistration();

  /// Sets input for the registration filter
  virtual void SetInput (irtkRealImage *, irtkRealImage *);

  /// Sets output for the registration filter
  virtual void SetOutput(irtkMultiLevelFreeFormTransformation *, irtkMultiLevelFreeFormTransformation *);

  /// Runs the registration filter
  virtual void Run();

  /// Read registration parameters from file
  virtual void Read (char *);

  /// Parse parameter line
  virtual bool Read(char *, char *, int &);

  /// Write registration parameters to file
  virtual void Write(char *);

  /// Write parameters to stream
  virtual void Write(ostream &);

  /// Guess parameters
  virtual void GuessParameter();

  // Access parameters
  virtual SetMacro(TargetPadding,      short);
  virtual GetMacro(TargetPadding,      short);
  virtual SetMacro(SourcePadding,      short);
  virtual GetMacro(SourcePadding,      short);
  virtual SetMacro(NumberOfLevels,     int);
  virtual GetMacro(NumberOfLevels,     int);
  virtual SetMacro(NumberOfIterations, int);
  virtual GetMacro(NumberOfIterations, int);
  virtual SetMacro(DebugFlag, int);
  virtual GetMacro(DebugFlag, int);


};

inline void irtkDemonsRegistration::SetInput(irtkRealImage *target, irtkRealImage *source)
{
  _target = target;
  _source = source;
}

inline void irtkDemonsRegistration::SetOutput(irtkMultiLevelFreeFormTransformation *transformation1,
    irtkMultiLevelFreeFormTransformation *transformation2)
{
  _transformation1 = transformation1;
  _transformation2 = transformation2;
}

// Function to read in a line from an istream, to be used by derived classes
extern int read_line(istream &, char *, char *&);

#endif
