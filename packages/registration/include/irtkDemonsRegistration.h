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

#include <irtkImage.h>

#include <irtkResampling.h>

#include <irtkImageFunction.h>

#include <irtkTransformation.h>

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

  /// Interface to input file stream
  friend istream& operator>> (istream&, irtkDemonsRegistration*);

  /// Interface to output file stream
  friend ostream& operator<< (ostream&, const irtkDemonsRegistration*);

protected:

  /** First input image. This image is denoted as target image and its
   *  coordinate system defines the frame of reference for the registration.
   */
  irtkGreyImage *_target;

  /** Second input image. This image is denoted as source image. The goal of
   *  the registration is to find the transformation which maps the source
   *  image into the coordinate system of the target image.
   */
  irtkGreyImage *_source;

  /** Temporary source image
   */
  irtkGreyImage _tmp;

  /** Gradient of the target image.
   */
  irtkGreyImage _targetGradientX, _targetGradientY, _targetGradientZ;

  /** Local displacement field.
   */
  irtkRealImage _localDX, _localDY, _localDZ;

  /** Total displacement field.
   */
  irtkRealImage _globalDX, _globalDY, _globalDZ;

  /// Output
  irtkMultiLevelFreeFormTransformation *_transformation;

  /// Linear free-form transformation
  irtkLinearFreeFormTransformation *_ffd;

  /// Interpolator
  irtkInterpolateImageFunction *_interpolator;

  /// Blurring of target image (in mm)
  double _TargetBlurring;

  /// Resolution of target image (in mm)
  double _TargetResolution;

  /// Padding value of target image
  short  _TargetPadding;

  /// Blurring of source image (in mm)
  double _SourceBlurring;

  /// Resolution of source image (in mm)
  double _SourceResolution;

  /// Padding value of source image
  short  _SourcePadding;

  /// Number of levels of multiresolution pyramid
  int    _NumberOfLevels;

  /// Max. number of iterations per step size
  int    _NumberOfIterations;

  /// Reduction factor in multi-resolution pyramid
  double _ReductionFactor;

  /// Smoothing of deformation field at every iteration
  double _Smoothing;

  /// Interpolation mode to use during resampling and registration
  irtkInterpolationMode _InterpolationMode;

  /// Source image domain which can be interpolated fast
  double _source_x1, _source_y1, _source_z1;
  double _source_x2, _source_y2, _source_z2;

  /// Initial set up for the registration
  virtual void Initialize();

  /// Initial set up for the registration
  virtual void Initialize(int);

  /// Final set up for the registration
  virtual void Finalize();

  /// Final set up for the registration
  virtual void Finalize(int);

public:

  /// Constructor
  irtkDemonsRegistration();

  /// Destructor
  virtual ~irtkDemonsRegistration();

  /// Sets input for the registration filter
  virtual void SetInput (irtkGreyImage *, irtkGreyImage *);

  /// Sets output for the registration filter
  virtual void SetOutput(irtkMultiLevelFreeFormTransformation *);

  /// Runs the registration filter
  virtual void Run();

  /// Runs the registration filter
  virtual void Run(irtkGreyImage, irtkGreyImage, int);

  /// Copy parameters
  virtual void SetParameter(const irtkDemonsRegistration *r);

  /// Read registration parameters from file
  virtual void Read (char *);

  /// Write registration parameters to file
  virtual void Write(char *);

  // Access parameters
  virtual SetMacro(TargetBlurring,     double);
  virtual GetMacro(TargetBlurring,     double);
  virtual SetMacro(TargetResolution,   double);
  virtual GetMacro(TargetResolution,   double);
  virtual SetMacro(TargetPadding,      short);
  virtual GetMacro(TargetPadding,      short);
  virtual SetMacro(SourceBlurring,     double);
  virtual GetMacro(SourceBlurring,     double);
  virtual SetMacro(SourceResolution,   double);
  virtual GetMacro(SourceResolution,   double);
  virtual SetMacro(SourcePadding,      short);
  virtual GetMacro(SourcePadding,      short);
  virtual SetMacro(NumberOfLevels,     int);
  virtual GetMacro(NumberOfLevels,     int);
  virtual SetMacro(NumberOfIterations, int);
  virtual GetMacro(NumberOfIterations, int);
  virtual SetMacro(ReductionFactor,    double);
  virtual GetMacro(ReductionFactor,    double);
  virtual SetMacro(Smoothing,          double);
  virtual GetMacro(Smoothing,          double);
  virtual SetMacro(InterpolationMode,  irtkInterpolationMode);
  virtual GetMacro(InterpolationMode,  irtkInterpolationMode);

};

inline void irtkDemonsRegistration::SetInput(irtkGreyImage *target, irtkGreyImage *source)
{
  _target = target;
  _source = source;
}

inline void irtkDemonsRegistration::SetOutput(irtkMultiLevelFreeFormTransformation *transformation)
{
  _transformation = transformation;
}

// Function to read in a line from an istream, to be used by derived classes
extern int read_line(istream &, char *, char *&);

#endif
