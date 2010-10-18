/*=========================================================================

  Library   : packages
  Module    : $RCSfile: vtkKDTreePointLocator.cxx,v $
  Authors   : Daniel Rueckert
  Copyright : Imperial College, Department of Computing
              Visual Information Processing (VIP), 2000-2001
  Purpose   :
  Date      : $Date: 2006-03-21 14:52:40 $
  Version   : $Revision: 1.6 $
  Changes   : $Locker:  $
	      $Log: vtkKDTreePointLocator.cxx,v $
	      Revision 1.6  2006-03-21 14:52:40  dr
	      Added compatibility to vtk 5.0

	      Revision 1.5  2005/12/22 14:43:17  dr
	      Updated surface registration class hierarchy
	
	      Revision 1.4  2005/04/04 14:13:29  dr
	      Bug fix
	
	      Revision 1.3  2004/08/12 10:06:34  dr
	      Added compatibility for VTK 4.4 and higher
	
	      Revision 1.2  2004/07/16 11:30:45  dr
	      Added #ifdef HAS_VTK
	
	      Revision 1.1  2004/06/15 18:39:14  tp500
	      Class written by 3d party needed for ND registrations
	

=========================================================================*/
/*
  Program:   Visualization Toolkit
    Language:  C++

  Copyright (c) 2001-2003 Gerald Dalley
  All rights reserved.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notice for more information.

 =========================================================================*/

#ifdef HAS_VTK

#include "vtkKDTreePointLocator.h"
#include "vtkMath.h"
#include "vtkIntArray.h"
#include "vtkPointSet.h"
#include "vtkImageData.h"
#include "vtkRectilinearGrid.h"
#if ((VTK_MAJOR_VERSION == 5) || ((VTK_MAJOR_VERSION == 4) && (VTK_MINOR_VERSION > 3)))
#include "vtkDoubleArray.h"
#else
#include "vtkFloatArray.h"
#endif
#include "vtkObjectFactory.h"
#include <assert.h>

vtkCxxRevisionMacro(vtkKDTreePointLocator, "$Revision: 1.6 $");
vtkStandardNewMacro(vtkKDTreePointLocator);

class PointData {
public:
    int PointId;
    double *coords;
    ~PointData() { delete [] coords; };
};

struct TerminalNodeData {
    vtkIdList *PointIDs;  // partition
};

struct InternalNodeData {
    double SplitValue;      // Xq
    int SplitAxis;         // discriminator

    vtkKDTreeNode *LeftChild;
    vtkKDTreeNode *RightChild;
};

class vtkKDTreeNode {
public:
    ~vtkKDTreeNode();

    bool Terminal;

    union {
        TerminalNodeData t;
        InternalNodeData i;
    } NodeData;
};

vtkKDTreeNode::~vtkKDTreeNode()
{
    if (this->Terminal) {
        this->NodeData.t.PointIDs->Delete();
    } else {
        delete this->NodeData.i.LeftChild;
        delete this->NodeData.i.RightChild;
    }
}

// =======================================================
//    M A I N   M E T H O D S
// =======================================================

vtkKDTreePointLocator::vtkKDTreePointLocator()
{
    this->NumberOfPointsPerBucket = 10;
    this->Root = NULL;
}

vtkKDTreePointLocator::~vtkKDTreePointLocator()
{
    this->FreeSearchStructure();
}

void vtkKDTreePointLocator::FreeSearchStructure()
{
    if (this->Root != NULL) 
    {
        delete Root;
        Root = NULL;
    }
    this->Modified();
}

void vtkKDTreePointLocator::BuildLocator()
{
    // Skip the build if the tree is already up-to-date
    if ( (this->BuildTime > this->MTime) && 
        (this->BuildTime > this->DataSet->GetMTime()) )
    {
        return;
    }

    // Get rid of old structure if it exists
    this->FreeSearchStructure();

    // Make sure we have a usable data set
    int numDims = this->GetDataDimension();
    if (numDims <= 0)
    {
        vtkErrorMacro(<< "input data set dimension cannot be determined "
            "(make sure you have some input data before calling "
            "BuildLocator)");
        throw;
    }

    // Transfer points to a PointData array
    double *point;
    int numPoints = this->DataSet->GetNumberOfPoints();
    PointData *pts = new PointData[numPoints]; 
    for (int i=0; i<numPoints; i++) {
        point = this->DataSet->GetPoint(i);
        pts[i].PointId = i;
        pts[i].coords = new double[numDims];
        for (int d=0; d<numDims; d++) 
        {            
            pts[i].coords[d] = point[d];
        }
    }

    // Build the tree
    this->Root = this->BuildTree(numPoints, pts);

    // Cleanup, etc.
    delete [] pts;
    this->BuildTime.Modified();
}

// Note: because I wanted to use the standard C qsort function and it can't
// be parameterized, I'm using a non-threadsafe method of doing the 
// comparisons.  You could make this threadsafe with mutexes or something
// like that if you want.
static int compareDimension = 0;

static int CompareDim(const void *e1, const void *e2) 
{
    if (((PointData*)e1)->coords[compareDimension] > 
        ((PointData*)e2)->coords[compareDimension]) 
    {
        return 1;
    } else {
        return -1;
    }
}

/* Unfortunately, the median position isn't always the best place to have
   a split point.  Consider the 1D point sequence:
        1.2, 1.2, 1.2, 1.2, 1.3, 1.4, 1.5
    If we selected the fourth entry to indicate the split point, we would
    use 1.2, but this would actually make the split happen before the first
    point.  In other words, we'd have no splitting happening at all.  This
    function is designed to search for the closest split point to the median
    position so that the split actually occurs.  Precondition: not all of the
    values are along the axis dimension have identical values.  Returns -1
    if an error is detected.  */
static inline int findSplitPoint(PointData *pts, int numPoints, int axis) 
{
    int medianPosition = int(ceil(numPoints/2.0));

    // First see if the median position is okay
    if (pts[medianPosition].coords[axis] == 
        pts[medianPosition-1].coords[axis])
    { 
        // Nope it's not, so search left and right from the median
        int leftSplitIdx = medianPosition - 1; 
        while ((leftSplitIdx >= 0) && 
            (pts[leftSplitIdx].coords[axis] == 
             pts[medianPosition].coords[axis])) 
        { 
            leftSplitIdx--; 
        } 
        int rightSplitIdx = medianPosition + 1; 
        while ((rightSplitIdx < numPoints) && 
            (pts[rightSplitIdx].coords[axis] == 
                pts[medianPosition].coords[axis])) 
        { 
            rightSplitIdx++; 
        } 
        if (leftSplitIdx >= 0) 
        { 
            if (rightSplitIdx < numPoints) 
            { 
                if ((medianPosition - leftSplitIdx) < 
                    (medianPosition - rightSplitIdx)) 
                { 
                    return leftSplitIdx+1; 
                } 
                else 
                { 
                    return rightSplitIdx; 
                } 
            } 
            else 
            { 
                return leftSplitIdx+1; 
            } 
        } 
        else 
        { 
            if (rightSplitIdx < numPoints) 
            { 
                return rightSplitIdx; 
            } 
            else 
            { 
                vtkGenericWarningMacro(
                    << "The precondition that all points "
                    << "do not have the same coordinates along the axis "
                    << "dimension was violated.");
                return -1;
            } 
        } 
    } 
    else 
    { 
        return medianPosition; 
    }
}

vtkKDTreeNode *vtkKDTreePointLocator::BuildTree(int numPoints, PointData *pts) 
{
    int d; // Dimension index variable
    int numDims = this->GetDataDimension();
    assert(numDims > 0); // Check precondition in debug mode

    // Create KD-Tree node
    vtkKDTreeNode *node = new vtkKDTreeNode();

    // Jump out if we just need to stuff the results 
    // in a terminal bucket
    if (numPoints <= this->NumberOfPointsPerBucket) {
        node->Terminal = true;

        vtkIdList *bucket = vtkIdList::New();
        bucket->SetNumberOfIds(numPoints);

        for (int i=0; i<numPoints; i++) 
        { 
            bucket->SetId(i, pts[i].PointId);
        }

        node->NodeData.t.PointIDs = bucket;

        return node;
    } 

    node->Terminal = false;

    // Find which axis has the greatest spread.  Then sort the points
    // based on the values for that axis.
    int bestDim = -1;
    double bestSpread = 0.0f;
    for (d = 0; d < numDims; d++)
    {
        double minVal =  VTK_FLOAT_MAX;
        double maxVal = -VTK_FLOAT_MAX;
        for (int i = 0; i < numPoints; i++)
        {
            if (pts[i].coords[d] < minVal) minVal = pts[i].coords[d];
            if (pts[i].coords[d] > maxVal) maxVal = pts[i].coords[d];
        }

        if ((maxVal - minVal) > bestSpread)
        {
            bestSpread = maxVal - minVal;
            bestDim = d;
        }
    }

    // Check for degenerate data sets where all spreads 
    // are zero and just dump the results into a single bin
    if (bestSpread == 0.0f)
    {
        node->Terminal = true;
        vtkIdList *bucket = vtkIdList::New();
        bucket->SetNumberOfIds(numPoints);
        for (int i=0; i<numPoints; i++) {
            bucket->SetId(i, pts[i].PointId);
        }
        node->NodeData.t.PointIDs = bucket;
        vtkWarningMacro(<< numPoints << " points, all with exactly the same "
            "coordinates were found.  They are being placed in one kd-tree "
            "bin.");
        return node;
    }

    // To make this threadsafe, lock a mutex here
    compareDimension = bestDim; 
    qsort((void*)pts, numPoints, sizeof(PointData), CompareDim);
    // To make this threadsafe, unlock the mutex here
    
    int splitPosition = findSplitPoint(pts, numPoints, bestDim);
    assert(splitPosition > 0); // should be true since max spread > 0

    // make split value be between points to make the searches be able to 
    // be slightly faster in some rare cases (otherwise just assign
    // pts[splitPosition].coords[bestDim] to the split value).
    node->NodeData.i.SplitValue = 
        (pts[splitPosition].coords[bestDim] + 
        pts[splitPosition-1].coords[bestDim]) / 2.0f; 
    node->NodeData.i.SplitAxis = bestDim;

    node->NodeData.i.LeftChild = this->BuildTree(splitPosition, pts); 
    node->NodeData.i.RightChild = this->BuildTree(                    
        numPoints - splitPosition, pts + splitPosition);             

    return node;
}

void vtkKDTreePointLocator::GenerateRepresentation(
    int vtkNotUsed(level), vtkPolyData *)
{
    vtkErrorMacro(<<"I don't understand what this method is "
        "supposed to do, so it's not implemented.");
}

#ifdef _DEBUG

static bool TreeIsConsistentRecurse(ofstream &treefile, 
    vtkKDTreeNode *node, bool *pointIsPresent, 
    vtkDataSet *ds, double *nodeExtent, int numDims)
{    
    int d;

    if (node->Terminal)
    {
        // Print the leaf tree node info to the xml file
        treefile << "<leaf ";
        for (d = 0; d < numDims; d++)
        {
            treefile << "min" << d << "=\"" << nodeExtent[d*2]   << "\" ";
            treefile << "max" << d << "=\"" << nodeExtent[d*2+1] << "\" ";
        }
        treefile << ">" << endl;

        // Mark of all points that are contained in this node.  In the 
        // process, make sure no other nodes contain this point and make 
        // sure that the point actually is within the nodeExtent
        int i;
        for (i=0; i < node->NodeData.t.PointIDs->GetNumberOfIds(); i++)
        {
            int dsID = node->NodeData.t.PointIDs->GetId(i);

            // Make sure that the point has not already been used and that the
            // point ID is a valid one
            if (dsID < 0)
            {
                vtkGenericWarningMacro(<< "Somehow a negative point ID "
                    << "was inserted in the tree: " << dsID);
                return false;
            }
            else if (dsID >= ds->GetNumberOfPoints())
            {
                vtkGenericWarningMacro(<< "Somehow point ID " << dsID 
                    << " was inserted into the tree, but the tree only has "
                    << ds->GetNumberOfPoints() << " points.");
                return false;
            }
            else if (pointIsPresent[dsID]) 
            {
                vtkGenericWarningMacro(<< "Point " << dsID 
                    << " is present in more than one node!");
                return false;
            }
            else 
            {
                // Mark the point as existing in the tree
                pointIsPresent[dsID] = true;
            }

            // Make sure the point is in the right node
            double *pt = ds->GetPoint(dsID);
            for (d = 0; d < numDims; d++)
            {
                if ((pt[d] < nodeExtent[d*2]) || (pt[d] >= nodeExtent[d*2+1]))
                {
                    vtkGenericWarningMacro(<< "Point " << dsID 
                        << " is in the wrong tree node! (dimension " 
                        << d << ")"); 
                    return false;
                }
            }

            // Write out the point info to the xml file
            treefile << "<point id=\"" << dsID << "\" ";
            for (d = 0; d < numDims; d++)
            {
                treefile << "x" << d << "=\"" << pt[d]   << "\" ";
            }
            treefile << "/>" << endl;
        }
        treefile << "</leaf>" << endl;

        return true;
    }
    else
    {
        treefile << "<node axis=\"" << node->NodeData.i.SplitAxis << "\" ";
        treefile << "splitValue=\"" << node->NodeData.i.SplitValue << "\" ";
        for (d = 0; d < numDims; d++)
        {
            treefile << "min" << d << "=\"" << nodeExtent[d*2]   << "\" ";
            treefile << "max" << d << "=\"" << nodeExtent[d*2+1] << "\" ";
        }
        treefile << ">" << endl;

        double tmpExtent;

        if (node->NodeData.i.SplitValue >= 
            nodeExtent[node->NodeData.i.SplitAxis*2 + 1])
        {
            vtkGenericWarningMacro(
                << "Attempting to split an internal node on "
                << "dimension " << node->NodeData.i.SplitAxis 
                << " at " << node->NodeData.i.SplitValue << ", but "
                << "this node's max extent is only " 
                << nodeExtent[node->NodeData.i.SplitAxis*2 + 1]);
            return false;
        }
        tmpExtent = nodeExtent[node->NodeData.i.SplitAxis*2 + 1];  
        nodeExtent[node->NodeData.i.SplitAxis*2 + 1] = 
            node->NodeData.i.SplitValue;
        if (!TreeIsConsistentRecurse(treefile, node->NodeData.i.LeftChild, 
            pointIsPresent, ds, nodeExtent, numDims))
        {
            // No vtkErrorMacro needed here because TreeIsConsistentRecurse
            // should have already printed the error
            return false;
        }
        nodeExtent[node->NodeData.i.SplitAxis*2 + 1] = tmpExtent;

        if (node->NodeData.i.SplitValue <= 
            nodeExtent[node->NodeData.i.SplitAxis*2])
        {
            vtkGenericWarningMacro(
                << "Attempting to split an internal node on "
                << "dimension " << node->NodeData.i.SplitAxis 
                << " at " << node->NodeData.i.SplitValue << ", but "
                << "this node's min extent is already " 
                << nodeExtent[node->NodeData.i.SplitAxis*2]);
            return false;
        }
        tmpExtent = nodeExtent[node->NodeData.i.SplitAxis*2];  
        nodeExtent[node->NodeData.i.SplitAxis*2] = node->NodeData.i.SplitValue;
        if (!TreeIsConsistentRecurse(treefile, node->NodeData.i.RightChild, 
            pointIsPresent, ds, nodeExtent, numDims))
        {
            // No vtkErrorMacro needed here because TreeIsConsistentRecurse
            // should have already printed the error
            return false;
        }
        nodeExtent[node->NodeData.i.SplitAxis*2] = tmpExtent;

        treefile << "</node>" << endl;

        return true;
    }
}

bool vtkKDTreePointLocator::TreeIsConsistent(const char *xmlOutputFName)
{
    this->BuildLocator();

    // Opens up the XML output file.  If xmlOutputFName is NULL, 
    // zero length, etc. ofstream should send the output to /dev/null
    ofstream treefile(xmlOutputFName);

    bool *pointIsPresent = new bool[this->DataSet->GetNumberOfPoints()];
    memset(pointIsPresent, 0, 
        sizeof(bool) * this->DataSet->GetNumberOfPoints());

    double *extent = new double[this->GetDataDimension()*2];
    for (int d = 0; d < this->GetDataDimension(); d++)
    {
        extent[d*2]     = -VTK_FLOAT_MAX;
        extent[d*2+1]   =  VTK_FLOAT_MAX;
    }

    if (!TreeIsConsistentRecurse(treefile,
        this->Root, pointIsPresent, this->DataSet, extent, 
        this->GetDataDimension()))
    {
        delete [] extent;
        return false;
    }
    delete [] extent;

    for (int i=0; i<this->DataSet->GetNumberOfPoints(); i++)
    {
        if (!pointIsPresent[i])
        {
            vtkErrorMacro(<< "Point " << i << " is missing from the tree!");
            return false;
        }
    }

    delete pointIsPresent;

    return true;
}
#endif

// Given a position x-y-z, return the id of the point closest to it.
int vtkKDTreePointLocator::FindClosestPoint(double x, double y, double z)
{
    double xyz[3];

    if (this->GetDataDimension() != 3)
    {
        vtkErrorMacro(<< "Attempting to find the closest point in 3D "
            << "coordinates, but the actual data dimensionality is "
            << this->GetDataDimension());
        return -1;
    }

    xyz[0] = x; xyz[1] = y; xyz[2] = z;
    return this->FindClosestPoint(xyz);
}

int vtkKDTreePointLocator::FindClosestPoint(double *x)
{
    this->BuildLocator();

    int closestPointId = -1;
    double closestPointDistanceSquared = VTK_FLOAT_MAX;
    int numDims = this->GetDataDimension();
    double *extent = new double[numDims*2];
    double *closestPoint = new double[numDims];
    for (int d = 0; d < numDims; d++)
    {
        extent[d*2]     = -VTK_FLOAT_MAX;
        extent[d*2+1]   =  VTK_FLOAT_MAX;
        closestPoint[d] =  VTK_FLOAT_MAX;
    }

    this->SearchForClosestPoint(x, this->Root, 
        closestPointId, closestPointDistanceSquared, closestPoint, extent);

    delete [] extent;
    delete [] closestPoint;

    return closestPointId;
}

// Some globals for performance measurement hacks and debugging.
#ifdef DEBUG
long numExecutes = 0;
long numTerminals = 0;
long numDrips = 0;
long numDoubleChecks = 0;
#endif

// Retrieves a linear point list from the data set if possible (and we
// know how to do so).  Otherwise returns NULL.
inline static double *getLinearPointListIfPossible(vtkDataSet *ds)
{
    if (ds == NULL)
    {
        vtkGenericWarningMacro(<<"No input data set has been assigned to this "
            "locator.");
        return NULL;
    }
    else if (ds->IsA("vtkPointSet"))
    {
        vtkDataArray *arr = ((vtkPointSet*)ds)->GetPoints()->GetData();
#if ((VTK_MAJOR_VERSION == 5) || ((VTK_MAJOR_VERSION == 4) && (VTK_MINOR_VERSION > 3)))
       if (arr->IsA("vtkDoubleArray"))
        {
            return ((vtkDoubleArray*)arr)->GetPointer(0);
        }
        else
        {
            return NULL;
        }
#else
        if (arr->IsA("vtkFloatArray"))
        {
            return ((vtkFloatArray*)arr)->GetPointer(0);
        }
        else
        {
            return NULL;
        }
#endif
    } 
    else if (ds->IsA("vtkImageData"))
    {
        // We'd have to build it the hard way, so don't bother
        return NULL;
    }
    else if (ds->IsA("vtkRectilinearGrid"))
    {
        // We'd have to build it the hard way, so don't bother
        return NULL;
    }
    else
    {        
        // We'd have to build it the hard way, so don't bother
        return NULL;
    }
}

bool vtkKDTreePointLocator::SearchForClosestPoint(
    double *x, vtkKDTreeNode *node, 
    int &closestPointId, double &closestPointDistanceSquared, 
    double *closestPoint, double *extent)
{
    int numDims = this->GetDataDimension();

#ifdef DEBUG
    numExecutes++;
#endif
    if (node->Terminal)
    {
#ifdef DEBUG
        numTerminals++;
#endif
        // If we have hit a terminal node, look for the closest closestPoint
        // in that node to x.

        int numIds = node->NodeData.t.PointIDs->GetNumberOfIds();
        int *ids = (int*) node->NodeData.t.PointIDs->GetPointer(0);
#ifdef DEBUG
        numDrips += numIds;
#endif
        
        // If we can use direct array indexing, we can save the overhead
        // of calling GetPoint all the time.
        double *pointList = getLinearPointListIfPossible(this->DataSet);

        if (pointList == NULL)
        {
            // If you modify anything here, make sure you make the
            // corresponding changes to the "else" block below
            for (int idNum=0; idNum<numIds; idNum++)
            {
                double *tmpClosestPoint = this->DataSet->GetPoint(ids[idNum]);
                // k-dimensional inline version of
                // vtkMath::Distance2BetweenPoints(x, tmpClosestPoint);
                double distSquared = 0.0f;
                for (int d=0; d<numDims; d++)
                {
                    distSquared += (x[d] - tmpClosestPoint[d]) * 
                        (x[d] - tmpClosestPoint[d]);
                }

                if (distSquared < closestPointDistanceSquared)
                {
                    closestPointDistanceSquared = distSquared;
                    closestPointId = ids[idNum];
                    closestPoint = tmpClosestPoint;
                }
            }
        }
        else 
        {           
            // If you modify anything here, make sure you make the
            // corresponding changes to the "then" block above
            for (int idNum=0; idNum<numIds; idNum++)
            {
                double *tmpClosestPoint = &pointList[ids[idNum]*numDims];
                // k-dimensional inline version of
                //vtkMath::Distance2BetweenPoints(x, tmpClosestPoint);
                double distSquared = 0.0f;
                for (int d=0; d<numDims; d++)
                {
                    distSquared += (x[d] - tmpClosestPoint[d]) * 
                        (x[d] - tmpClosestPoint[d]);
                }

                if (distSquared < closestPointDistanceSquared)
                {
                    closestPointDistanceSquared = distSquared;
                    closestPointId = ids[idNum];
                    closestPoint = tmpClosestPoint;
                }
            }
        }
    } 
    else 
    {
        double tmpExtent;
        bool done;

        // We're not at a terminal node, so we need to figure out
        // which subtree to search
        if (x[node->NodeData.i.SplitAxis] < node->NodeData.i.SplitValue)
        {
            // First case: x is contained in the left subtree
            // Search for the closest point in that subtree.
            tmpExtent = extent[node->NodeData.i.SplitAxis*2+1];
            extent[node->NodeData.i.SplitAxis*2+1] = 
                node->NodeData.i.SplitValue;
            done = this->SearchForClosestPoint(x, node->NodeData.i.LeftChild,
                closestPointId, closestPointDistanceSquared, 
                closestPoint, extent);
            extent[node->NodeData.i.SplitAxis*2+1] = tmpExtent;

            // Closest point definitively found in the subtree.  We
            // can jump out and return here.
            if (done) return true;

            // Second case: we did not definitively find x in the
            // closer subtree, so we'll check the farther one.

            // Check other subtree if there could be a closer point
            // in it.
            tmpExtent = extent[node->NodeData.i.SplitAxis*2];
            extent[node->NodeData.i.SplitAxis*2] = node->NodeData.i.SplitValue;
            if (this->BoundsOverlapBall(x, 
                closestPointDistanceSquared, extent))
            {
#ifdef DEBUG
                numDoubleChecks++;
#endif
                this->SearchForClosestPoint(x, node->NodeData.i.RightChild,
                    closestPointId, closestPointDistanceSquared, 
                    closestPoint, extent);
            }
            extent[node->NodeData.i.SplitAxis*2] = tmpExtent;
        }
        else 
        {
            // First case: x is contained in the right subtree
            // Search for the closest point in that subtree.
            tmpExtent = extent[node->NodeData.i.SplitAxis*2];
            extent[node->NodeData.i.SplitAxis*2] = 
                node->NodeData.i.SplitValue;
            done = this->SearchForClosestPoint(x, node->NodeData.i.RightChild,
                closestPointId, closestPointDistanceSquared, 
                closestPoint, extent);
            extent[node->NodeData.i.SplitAxis*2] = tmpExtent;

            // Closest point definitively found in the subtree.  We
            // can jump out and return here.
            if (done) return true;


            // Second case: we did not definitively find x in the
            // closer subtree, so we'll check the farther one.

            // Check other subtree
            tmpExtent = extent[node->NodeData.i.SplitAxis*2+1];
            extent[node->NodeData.i.SplitAxis*2+1] = 
                node->NodeData.i.SplitValue;
            if (this->BoundsOverlapBall(x, 
                closestPointDistanceSquared, extent))
            {
#ifdef DEBUG
                numDoubleChecks++;
#endif
                this->SearchForClosestPoint(x, node->NodeData.i.LeftChild,
                    closestPointId, closestPointDistanceSquared, 
                    closestPoint, extent);
            }
            extent[node->NodeData.i.SplitAxis*2+1] = tmpExtent;
        }
    }

    // BALL_WITHIN_BOUNDS ...

    // Finally, see if we can guarantee that the closest closestPoint is
    // closer to x than any of the edges of the current kD-tree
    // node.  If it is not, then some other node may contain a
    // closestPoint closer than the one just found.  We return false if
    // more work needs to be done to guarantee that we found the closest
    // point.
    double xToPoint2;
    double xToBoundary2;

    xToPoint2 = 0.0f;
    for (int d=0; d<numDims; d++)
    {
        xToPoint2 += (closestPoint[d] - x[d]) * (closestPoint[d] - x[d]);
    }

    double *ePtr = extent;
    double *xPtr = x;
    for (int d=0; d<numDims; d++)
    {
        // Check min boundary
        xToBoundary2 = (*xPtr) - (*ePtr); ePtr++;
        //xToBoundary2 *= xToBoundary2;
        if (xToBoundary2 <= xToPoint2) return false;

        // Check max boundary
        xToBoundary2 = (*xPtr) - (*ePtr); ePtr++;
        //xToBoundary2 *= xToBoundary2;
        if (xToBoundary2 <= xToPoint2) return false;

        xPtr++;
    }
    return true;
}

bool vtkKDTreePointLocator::BoundsOverlapBall(
    double *x, double closestPointDistanceSquared, 
    double *extent)
{
    double xToBoundary;
    double *ePtr = extent;
    double *xPtr = x;

    double xToExtentCorner = 0.0f;

    int numDims = this->GetDataDimension();
    for (int d=0; d<numDims; d++) 
    {
        // Check min boundary
        xToBoundary = (*xPtr) - (*ePtr); ePtr++;
        if (xToBoundary < 0.0f) {
            xToExtentCorner +=  xToBoundary * xToBoundary;
            if (xToExtentCorner > closestPointDistanceSquared) return false;
        }

        // Check max boundary
        xToBoundary = (*xPtr) - (*ePtr); ePtr++;
        if (xToBoundary > 0.0f) {
            xToExtentCorner +=  xToBoundary * xToBoundary;
            if (xToExtentCorner > closestPointDistanceSquared) return false;
        }

        xPtr++;        
    }

    return true; 
}

int vtkKDTreePointLocator::GetDataDimension()
{
    if (this->DataSet == NULL)
    {
        vtkWarningMacro(<<"No input data set has been assigned to this "
            "locator.  Uncomputable data dimensionality.");
        return -1;
    }
    else if (this->DataSet->IsA("vtkPointSet"))
    {
        return ((vtkPointSet*)this->DataSet)->
            GetPoints()->GetData()->GetNumberOfComponents();
    } 
    else if (this->DataSet->IsA("vtkImageData"))
    {
        return ((vtkImageData*)this->DataSet)->
            GetDataDimension();
    }
    else if (this->DataSet->IsA("vtkRectilinearGrid"))
    {
        return ((vtkRectilinearGrid*)this->DataSet)->
            GetDataDimension();
    }
    else
    {        
        vtkWarningMacro(<< "Uncomputable data dimension. " 
            << "Unsupported DataSet type: " << this->DataSet->GetClassName());
        return -1;
    }
}

#endif
