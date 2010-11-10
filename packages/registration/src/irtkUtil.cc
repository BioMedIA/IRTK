/*=========================================================================

  Library   : Image Registration Toolkit (IRTK)
  Module    : $Id$
  Copyright : Imperial College, Department of Computing
              Visual Information Processing (VIP), 2008 onwards
  Date      : $Date$
  Version   : $Revision$
  Changes   : $Author$

=========================================================================*/

#include <irtkRegistration.h>

double combine_similarity(irtkSimilarityMetric *s1, irtkSimilarityMetric *s2, double p1, double p2)
{
	double combined;

	double factor;

	double sv1,sv2;

	//combined similarity = (sv1*p1 + sv2*p2)/(p1 + p2)

	sv1 = s1->Evaluate();
	sv2 = s2->Evaluate();
	factor = p1 + p2;

	if(factor > 0)
		combined = (sv1*p1 + sv2*p2)/factor;
	else
		combined = 0;

	return combined;
}

void irtkPadding(irtkGreyImage &image, irtkGreyPixel padding)
{
  int i, j, k, l, m, n, p, t;

  // Calculate padding
  m = 0;
  n = 0;

  for (t = 0; t < image.GetT(); t++) {
    for (i = 0; i < image.GetX(); i++) {
      for (j = 0; j < image.GetY(); j++) {
        for (k = 0; k < image.GetZ(); k++) {
          if (image(i, j, k, t) < 0) {
            // Count no. of padded voxels
            n++;
          } else {
            // Count no. of unpadded voxels
            m++;
          }
        }
      }
    }
  }

  if (n > 0) {

    // Print some padding information
    cout << "Padding value = " << padding << endl;
    cout << "Padding ratio = " << double(100*n)/double(m+n) << " %" << endl;

    // Calculate distances
    for (t = 0; t < image.GetT(); t++) {
      for (k = 0; k < image.GetZ(); k++) {
        for (j = 0; j < image.GetY(); j++) {
          for (i = 0; i < image.GetX(); i++) {
            if (image(i, j, k, t) == -1) {
              for (l = i; l < image.GetX(); l++) {
                if (image(l, j, k, t) != -1) {
                  break;
                }
              }
              for (p = i; p < l; p++) {
                image(p, j, k, t) = - (l - p);
              }
              i = l - 1;
            }
          }
        }
      }
    }
  }
}

void irtkPadding(irtkGreyImage &image, irtkGreyPixel padding, irtkFreeFormTransformation3D *ffd)
{
  int i, j, k, x, y, z, x1, y1, z1, x2, y2, z2, ok, index;
  int t;

  // Calculate number of active and passive control points
  for (i = 0; i < ffd->GetX(); i++) {
    for (j = 0; j < ffd->GetY(); j++) {
      for (k = 0; k < ffd->GetZ(); k++) {
        // Convert control points to index
        index = ffd->LatticeToIndex(i, j, k);

        // Calculate bounding box of control point in voxels
        ffd->BoundingBox(&image, index, x1, y1, z1, x2, y2, z2);

        ok = false;
        for (t = 0; t < image.GetT(); t++) {
          for (z = z1; z <= z2; z++) {
            for (y = y1; y <= y2; y++) {
              for (x = x1; x <= x2; x++) {
                if (image(x, y, z, t) > padding) {
                  ok = true;
                }
              }
            }
          }
        }
        if (ok == false) {
          ffd->PutStatus(i, j, k, _Passive);
        }
      }
    }
  }
}

void irtkPadding(irtkGreyImage **image, irtkGreyPixel padding, irtkFreeFormTransformation3D *ffd, int numberOfImages)
{
  int i, j, k, x, y, z, x1, y1, z1, x2, y2, z2, ok, index;
  int t, n;

  // Calculate number of active and passive control points
  for (i = 0; i < ffd->GetX(); i++) {
    for (j = 0; j < ffd->GetY(); j++) {
      for (k = 0; k < ffd->GetZ(); k++) {
        // Convert control points to index
        index = ffd->LatticeToIndex(i, j, k);

		ok = false;

		for (n = 0; n < numberOfImages; n++){
			// Calculate bounding box of control point in voxels
			ffd->BoundingBox(image[n], index, x1, y1, z1, x2, y2, z2);
			for (t = 0; t < image[n]->GetT(); t++) {
				for (z = z1; z <= z2; z++) {
					for (y = y1; y <= y2; y++) {
						for (x = x1; x <= x2; x++) {
							if (image[n]->GetAsDouble(x, y, z, t) > padding) {
								ok = true;
							}
						}
					}
				}
			}
		}
        if (ok == false) {
          ffd->PutStatus(i, j, k, _Passive);
        }
      }
    }
  }
}

double GuessResolution(double xsize, double ysize, double zsize)
{
  if ((xsize <= ysize) && (xsize <= zsize)) return xsize;
  if ((ysize <= xsize) && (ysize <= zsize)) return ysize;
  return zsize;
}

double GuessResolution(double xsize, double ysize)
{
  if (xsize > ysize) return xsize;
  return ysize;
}

int GuessPadding(irtkGreyImage &image)
{
  bool padding = true;

  if (image(0, 0, 0) != image(image.GetX() - 1, 0, 0)) padding = false;
  if (image(0, 0, 0) != image(0, image.GetY() - 1, 0)) padding = false;
  if (image(0, 0, 0) != image(0, 0, image.GetZ() - 1)) padding = false;
  if (image(0, 0, 0) != image(image.GetX() - 1, image.GetY() - 1, 0)) padding = false;
  if (image(0, 0, 0) != image(0, image.GetY() - 1, image.GetZ() - 1)) padding = false;
  if (image(0, 0, 0) != image(image.GetX() - 1, 0, image.GetZ() - 1)) padding = false;
  if (image(0, 0, 0) != image(image.GetX() - 1, image.GetY() - 1, image.GetZ() - 1)) padding = false;

  if (padding == false) {
    return MIN_GREY;
  } else {
    return image(0, 0, 0);
  }
}

int irtkCalculateNumberOfBins(irtkGreyImage *image, int maxbin, int min, int max)
{
  int i, nbins, range, width;
  irtkGreyPixel *ptr;

  // Find the min and max intensities
  range = max - min + 1;
  nbins = max - min + 1;
  width = 1;

  // Calculate number of bins to use
  if (maxbin > 0) {
    while (int(ceil(range/(double)width)) > maxbin) {
      width++;
    }
    nbins = int(ceil(range/(double)width));

    // Print out number of bins
    cout << "Using " << nbins << " out of " << maxbin << " bin(s) with width "
         << width << endl;
  } else {
    // Print out number of bins
    cout << "Using " << nbins << " bin(s) with width " 	 << width << endl;
  }

  // Rescale intensities to the number of bins
  ptr = image->GetPointerToVoxels();
  for (i = 0; i < image->GetNumberOfVoxels(); i++) {
    if (*ptr > 0) {
      *ptr = int(*ptr/(double)width);
    }
    ptr++;
  }

  // Return number of bins
  return nbins;
}

int irtkCalculateNumberOfBins(irtkGreyImage **image, int maxbin, int min, int max, int n)
{
  int i, j, nbins, range, width;
  irtkGreyPixel *ptr;

  // Find the min and max intensities
  range = max - min + 1;
  nbins = max - min + 1;
  width = 1;

  // Calculate number of bins to use
  if (maxbin > 0) {
    while (int(ceil(range/(double)width)) > maxbin) {
      width++;
    }
    nbins = int(ceil(range/(double)width));

    // Print out number of bins
    cout << "Using " << nbins << " out of " << maxbin << " bin(s) with width "
         << width << endl;
  } else {
    // Print out number of bins
    cout << "Using " << nbins << " bin(s) with width " 	 << width << endl;
  }

  // Rescale intensities to the number of bins
  for (j = 0; j < n; j++) {
    ptr = image[j]->GetPointerToVoxels();
    for (i = 0; i < image[j]->GetNumberOfVoxels(); i++) {
      if (*ptr > 0) {
        *ptr = int(*ptr/(double)width);
      }
      ptr++;
    }
  }

  // Return number of bins
  return nbins;
}

int read_line(istream &in, char *buffer1, char *&buffer2)
{
  char c;

  do {
    if (in.eof() == true) return 0;
    in.getline(buffer1, 255);
    c = buffer1[0];
  } while ((strlen(buffer1) == 0) || (c == '#') || (c == 13));

  if ((buffer2 = strchr(buffer1, '=')) == NULL) {
    cerr << "No valid line format\n";
    exit(1);
  }
  do {
    buffer2++;
  } while ((*buffer2 == ' ') || (*buffer2 == '\t'));

  return strlen(buffer1);
}

#ifdef HAS_VTK

#include <vtkPolyData.h>
#include <vtkCleanPolyData.h>
#include <vtkFloatArray.h>
#include <vtkPointData.h>
#include <vtkFeatureEdges.h>
#include <vtkPointLocator.h>
#include <vtkGenericCell.h>

void MarkBoundary(vtkPolyData *polydata)
{
  double x[3];
  int i, id, no_polydata, no_edges;

  // Detect edges
  vtkFeatureEdges *edges = vtkFeatureEdges::New();
  edges->BoundaryEdgesOn();
  edges->FeatureEdgesOff();
  edges->ManifoldEdgesOff();
  edges->NonManifoldEdgesOff();
  edges->SetColoring(1);
  edges->SetInput(polydata);
  edges->Update();

  // Calculate number of points
  no_polydata = polydata->GetNumberOfPoints();
  no_edges    = edges->GetOutput()->GetNumberOfPoints();

  // Setup locator
  vtkPointLocator *locator = vtkPointLocator::New();
  locator->SetDataSet(polydata);
  locator->BuildLocator();

  // Create scalar array to indicate if point is on edge or not
  vtkFloatArray *scalars = vtkFloatArray::New();
  scalars->SetNumberOfTuples(no_polydata);

  for (i = 0; i < no_edges; i++) {
    edges->GetOutput()->GetPoint(i, x);
    id = locator->FindClosestPoint(x);
    // If point is an edge, set scalar for point to 0
    scalars->InsertTuple1(id, 0);
  }

  for (i = 0; i < no_polydata; i++) {
    // If point has no been marked as edge, set scalar for point to 1
    if (scalars->GetValue(i) != 0) {
      scalars->InsertTuple1(i, 1);
    }
  }

  // Set up name
  scalars->SetName("EDGEPOINTS");
  polydata->GetPointData()->SetScalars(scalars);
}

void GetConnectedVertices(vtkSmartPointer<vtkPolyData> mesh, int seed, vtkSmartPointer<vtkIdList> connectedVertices)
{

	//get all cells that vertex 'seed' is a part of
	vtkSmartPointer<vtkIdList> cellIdList =
		vtkSmartPointer<vtkIdList>::New();
	mesh->GetPointCells(seed, cellIdList);

	//cout << "There are " << cellIdList->GetNumberOfIds() << " cells that use point " << seed << endl;
		//loop through all the cells that use the seed point
		for(vtkIdType i = 0; i < cellIdList->GetNumberOfIds(); i++)
		{

			vtkCell* cell = mesh->GetCell(cellIdList->GetId(i));
			//cout << "The cell has " << cell->GetNumberOfEdges() << " edges." << endl;

			//if the cell doesn't have any edges, it is a line
			if(cell->GetNumberOfEdges() <= 0)
			{
				continue;
			}

			for(vtkIdType e = 0; e < cell->GetNumberOfEdges(); e++)
			{
				vtkCell* edge = cell->GetEdge(e);

				vtkIdList* pointIdList = edge->GetPointIds();
				//cout << "This cell uses " << pointIdList->GetNumberOfIds() <<" points" << endl;
				/*
				for(vtkIdType p = 0; p < pointIdList->GetNumberOfIds(); p++)
				{
				cout << "Edge " << i << " uses point " << pointIdList->GetId(p) << endl;
				}
				*/
				if(pointIdList->GetId(0) == seed || pointIdList->GetId(1) == seed)
				{
					if(pointIdList->GetId(0) == seed)
					{
						connectedVertices->InsertUniqueId(pointIdList->GetId(1));
					}
					else
					{
						connectedVertices->InsertUniqueId(pointIdList->GetId(0));
					}
				}
			}
		}
		//cout << "There are " << connectedVertices->GetNumberOfIds() << " points connected to point " << seed << endl;
}

#endif
