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

#ifndef __vtkKDTreePointLocator_h
#define __vtkKDTreePointLocator_h

#include "vtkLocator.h"
#include "vtkPoints.h"
#include "vtkIdList.h"
#include <irtkImage.h>

class vtkKDTreeNode;
class PointData;

class vtkKDTreePointLocator : public vtkLocator
{
public:
  ~vtkKDTreePointLocator();
  static vtkKDTreePointLocator *New();
  const char *GetClassName() {
    return "vtkKDTreePointLocator";
  };
  //void PrintSelf(ostream& os, vtkIndent indent);

  vtkTypeRevisionMacro(vtkKDTreePointLocator,vtkLocator);

  // Description:
  // Construct with 10 records per terminal bucket
  vtkKDTreePointLocator();

  // Description:
  // Specify the average number of points in each bucket.
  vtkSetClampMacro(NumberOfPointsPerBucket,int,1,VTK_LARGE_INTEGER);
  vtkGetMacro(NumberOfPointsPerBucket,int);

  // Description:
  // Given a position x, return the id of the point closest to it. The
  // dimensionality of x should be the same as GetDataDimension
  // Calls BuildLocator, so this method may throw an exception.
  virtual int FindClosestPoint(double *x);

  // Description:
  // Given a position x when the input data dimensionality is 3,
  // returns the id of the point closest to it.  If GetDataDimension != 3,
  // returns -1.
  // Calls BuildLocator, so this method may throw an exception.
  int FindClosestPoint(double x, double y, double z);

  // Description:
  // Find the closest N points to a position. This returns the closest
  // N points to a position. A faster method could be created that returned
  // N close points to a position, but neccesarily the exact N closest.
  // The returned points are sorted from closest to farthest.
  /*  virtual void FindClosestNPoints(int N, float x[3], vtkIdList *result)
    virtual void FindClosestNPoints(int N, float x, float y, float z,
  				  vtkIdList *result);
  */
  // Description:
  // Find all points within a specified radius R of position x.
  // The result is not sorted in any specific manner.
  /*  virtual void FindPointsWithinRadius(float R, float x[3], vtkIdList *result);
    virtual void FindPointsWithinRadius(float R, float x, float y, float z,
  				      vtkIdList *result);
  */
  // Description:
  // See vtkLocator interface documentation.
  void FreeSearchStructure();

  // Description:
  // Builds the internal search structure.  If no input data set has been
  // assigned or if its dimensionality is indeterminate, this method
  // throws an exception.
  // Note: this method is NOT thread-safe
  void BuildLocator();

  // This method is unimplemented.  The interface description is vague
  // enough that I don't know what it's supposed to do.
  void GenerateRepresentation(int level, vtkPolyData *pd);

#ifdef _DEBUG
  // Description:
  // Test method that verifies internal data structure consistency.
  // xmlOutputFName is the filename of an XML file to be written that
  // encodes the KD-Tree.  If xmlOutputFName == NULL, no file is written.
  // TreeIsConsistent returns false if:
  //    1) Any points are found in the wrong KDTree bins
  //    2) Any points were left out of the KDTree
  //    3) Any points were included in two tree bins
  //    4) Any illegal point IDs are stored in the tree
  // Calls BuildLocator, so this method may throw an exception.
  bool TreeIsConsistent(const char *xmlOutputFName = NULL);
#endif

  // Description:
  // Determines the dimensionality of the input data.  Note: this only works
  // if the data set is a vtkPointSet, vtkImageData, or vtkRectilinearGrid.
  // For all other data types, -1 is returned and vtkWarningMacro is called.
  virtual int GetDataDimension();

protected:
  // Description:
  // Recursively builds a KD-tree from the points given.
  // Precondition: GetDataDimension should return a valid result
  // Note: this method is NOT thread-safe.
  vtkKDTreeNode *BuildTree(int numPoints, PointData *pts);

  int NumberOfPointsPerBucket; //Used with previous boolean to control subdivide
  vtkKDTreeNode *Root;

  // Description:
  // Recursive method that finds the single closest data set point to
  // a given point in 3D.  Based on [Friedman, pg. 224-5].
  // Assumes that the tree structure has been built before calling.
  bool SearchForClosestPoint(
    double *x, vtkKDTreeNode *node,
    int &closestPointId, double &closestPointDistanceSquared,
    double *closestPoint, double extent[6]);

  // Description:
  // Based on [Friedman, pg. 225].  Slightly adapted.  We used the Euclidean
  // distance for the DISSIM and COORDINATE DISTANCE functions.  Instead of
  // taking the square root (e.g. doing DISSIM) on the sum, we squared the
  // PQD[1] number (e.g. doing inverse(DISSIM)).
  // Assumes that the tree structure has been built before calling.
  bool BoundsOverlapBall(double *x,
                         double closestPointDistanceSquared, double extent[6]);
};

#endif

#endif
