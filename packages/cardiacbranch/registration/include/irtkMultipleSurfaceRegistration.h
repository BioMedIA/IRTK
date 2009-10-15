/*=========================================================================

  Library   : Image Registration Toolkit (IRTK)
  Module    : $Id$
  Copyright : Imperial College, Department of Computing
              Visual Information Processing (VIP), 2008 onwards
  Date      : $Date$
  Version   : $Revision$
  Changes   : $Author$

=========================================================================*/

#ifdef HAS_VTK

#ifndef _IRTKMULTIPLESURFACEREGISTRATION_H

#define _IRTKMULTIPLESURFACEREGISTRATION_H

#include <vtkPointData.h>
#include <vtkPolyData.h>
#include <vtkPolyDataWriter.h>
#include <vtkPolyDataReader.h>

#include <irtkRegistration.h>

#include <irtkUtil.h>
#include <irtkLocator.h>
#include <irtkOptimizer.h>

/**
 * Generic for surface registration based on the ICP algorithm.
 *
 * This class implements a registration filter which takes two input surfaces
 * and calculates the transformation which maps the second surface (denoted as
 * source surface) into the coordinate system of the first surface (denoted as
 * target surface).  This is the abstract base class which defines a common
 * interface for arbitrary registration filters. Each derived class has to
 * implement all abstract member functions.
 *
 */

class irtkMultipleSurfaceRegistration : public irtkObject
{

protected:

  /// Target point set
  vtkPolyData **_target;

  /// Source point set
  vtkPolyData **_source;

  /// Number of surfaces
  int _NumberOfSurfaces;

  /// Output
  irtkTransformation *_transformation;

  /// Locator for target points
  irtkLocator **_target_locator;

  /// Locator for source points
  irtkLocator **_source_locator;

  /// Pointer to point registration
  irtkPointRegistration *_preg;

  /// Initial set up for the registration
  virtual void Initialize();

  /// Final set up for the registration
  virtual void Finalize();

  /// Optimize registration
  virtual void Optimize();

  /// Convergence factor of registration
  double _Epsilon;

  /// Number of iterations
  int _NumberOfIterations;

  /// Flag whether to use symmetric distance function
  bool _UseSymmetricDistance;

public:

  /// Constructor
  irtkMultipleSurfaceRegistration();

  /// Destructor
  virtual ~irtkMultipleSurfaceRegistration();

  /// Sets input for the registration filter
  virtual void SetInput (vtkPolyData **, vtkPolyData **);

  /// Sets output for the registration filter
  virtual void SetOutput(irtkTransformation *);

  /// Run the filter
  virtual void Run();

  /// Returns the name of the class
  virtual const char *NameOfClass();

  /// Set locator
  virtual void SetTargetLocator(irtkLocator **);

  /// Set locator
  virtual void SetSourceLocator(irtkLocator **);

  /// Set epsilon
  virtual void SetEpsilon(float);

  /// Get epsilon
  virtual double GetEpsilon();

  /// Use symmetric distance function
  void UseSymmetricDistance();

  // Access parameters
  virtual SetMacro(NumberOfSurfaces, int);
  virtual GetMacro(NumberOfSurfaces, int);
  virtual SetMacro(NumberOfIterations, int);
  virtual GetMacro(NumberOfIterations, int);
};

inline void irtkMultipleSurfaceRegistration::SetTargetLocator(irtkLocator **locator)
{
  _target_locator = locator;
}

inline void irtkMultipleSurfaceRegistration::SetSourceLocator(irtkLocator **locator)
{
  _source_locator = locator;
}

inline void irtkMultipleSurfaceRegistration::SetEpsilon(float _e)
{
  _Epsilon = _e;
}

inline double irtkMultipleSurfaceRegistration::GetEpsilon()
{
  return _Epsilon;
}

inline void irtkMultipleSurfaceRegistration::UseSymmetricDistance()
{
  _UseSymmetricDistance = 1;
}

#include <irtkMultipleSurfaceRigidRegistration.h>
#include <irtkMultipleSurfaceAffineRegistration.h>
#include <irtkMultipleSurfaceFreeFormRegistration.h>

#endif

#endif
