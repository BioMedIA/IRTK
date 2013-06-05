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

double combine_mysimilarity(double sv1, double sv2, double p1, double p2)
{
    double combined;

    double factor;

    factor = p1 + p2;

    if(factor > 0)
        combined = (sv1*p1 + sv2*p2)/factor;
    else
        combined = 0;

    return combined;
}

double combine_mysimilarity(irtkSimilarityMetric **s, double *weight, double number)
{
    double combined, factor;
    int i;

    factor = 0; combined = 0;

    for (i = 0; i < number; i++) {
        if ( weight[i] > 0 ) {
            combined += s[i]->Evaluate() * weight[i];
            factor += weight[i];
        }
    }

    if(factor > 0)
        combined = combined/factor;
    else
        combined = 0;

    return combined;
}

void irtkPadding(irtkRealImage &image, irtkRealPixel padding, irtkGreyImage *result)
{
    int i, j, k, l, m, n, p, t;

    // Calculate padding
    m = 0;
    n = 0;

    for (t = 0; t < image.GetT(); t++) {
        for (i = 0; i < image.GetX(); i++) {
            for (j = 0; j < image.GetY(); j++) {
                for (k = 0; k < image.GetZ(); k++) {
                    if (image(i, j, k, t) <= padding) {
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
        cout << "Padding ratio = " << 100*double(n)/(double(m)+double(n)) << " %" << endl;

        // Calculate distances
        for (t = 0; t < image.GetT(); t++) {
            for (k = 0; k < image.GetZ(); k++) {
                for (j = 0; j < image.GetY(); j++) {
                    for (i = 0; i < image.GetX(); i++) {
                        if (image(i, j, k, t) <= padding) {
                            for (l = i; l < image.GetX(); l++) {
                                if (image(l, j, k, t) > padding) {
                                    break;
                                }
                            }
                            for (p = i; p < l; p++) {
                            	result->Put(p, j, k, t, l - p);
                            }
                            i = l - 1;
                        }
                    }
                }
            }
        }
    }
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
        cout << "Padding ratio = " << 100*double(n)/(double(m)+double(n)) << " %" << endl;

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

void irtkPadding(irtkRealImage &image, irtkRealPixel padding, irtkFreeFormTransformation3D *ffd)
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
                ffd->BoundingBoxImage(&image, index, x1, y1, z1, x2, y2, z2, 0.5);

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
                    ffd->PutStatusCP(i, j, k, _Passive, _Passive, _Passive);
                }
            }
        }
    }
}

int irtkGetBinIndex(irtkRealPixel pix, int min, int max, int nbins) {
  int val;

  val = (int)(((pix - min)/(max - min)) * (double)(nbins));

  if (val < 0) val = 0;
  if (val > nbins - 1) val = nbins - 1;

  return val;
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
                ffd->BoundingBoxImage(&image, index, x1, y1, z1, x2, y2, z2, 0.5);

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
                    ffd->PutStatusCP(i, j, k, _Passive, _Passive, _Passive);
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

                for (n = 0; n < numberOfImages; n++) {
                    // Calculate bounding box of control point in voxels
                    ffd->BoundingBoxImage(image[n], index, x1, y1, z1, x2, y2, z2);
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
                    ffd->PutStatusCP(i, j, k, _Passive, _Passive, _Passive);
                }
            }
        }
    }
}

void irtkPadding(irtkGreyImage **image, irtkGreyPixel padding, irtkBSplineFreeFormTransformationPeriodic *ffd, int numberOfImages, double *time)
{
    int i, j, k, x, y, z, i1, j1, k1, i2, j2, k2, ok, index;
    int t, n, number;
    double t1, t2, tt;

    number = 0;

    // Loop over control points
    for (t = 0; t < ffd->GetT(); t++) {
        for (z = 0; z < ffd->GetZ(); z++) {
            for (y = 0; y < ffd->GetY(); y++) {
                for (x = 0; x < ffd->GetX(); x++) {
                    ok = false;
                    index = ffd->LatticeToIndex(x, y, z, t);
                    // t1, t2 not in lattice coordinates at the moment!!!!!!!
                    ffd->BoundingBoxImage(image[0], index, i1, j1, k1, i2, j2, k2, t1, t2, 1.0);

                    // loop over all target images
                    for (n = 0; n < numberOfImages; n++) {
                        // transform time point of current target image to lattice coordinates and check for periodicity
                        tt = time[n];
                        // map time to relative time intervall [0,1]
                        while (tt < 0)
                            tt += 1.;
                        while (tt >= 1)
                            tt -= 1.;
                        tt = ffd->TimeToLattice(tt);

                        // check whether time point of current target image is in temporal bounding box (check in lattice coord.)
                        if (( (t1 >= 0) 
                            && (t2 < ffd->GetT()-1) 
                            &&  (tt >= t1) 
                            && (tt<=t2) )
                            || ( (t1 <  0) 
                            && ( (tt <= t2) 
                            || (tt >= t1+ffd->GetT()-1) ) )
                            || ( (t2 >= ffd->GetT()-1) 
                            && ( (tt >= t1) 
                            || (tt <= t2-ffd->GetT()+1) ) ) ) {

                                // Loop over all voxels in the target (reference) volume
                                for (k = k1; k <= k2; k++) {
                                    for (j = j1; j <= j2; j++) {
                                        for (i = i1; i <= i2; i++) {
                                            if (image[n]->GetAsDouble(i, j, k) > padding) {
                                                ok = true;
                                            }
                                        }
                                    }
                                }
                        }
                    }
                    if (ok == false) {
                        ffd->PutStatusCP(x, y, z, t, _Passive, _Passive, _Passive);
                        number++;
                    }
                }
            }
        }
    }
    cout << "Number of CP padded: " << number << " out of " << ffd->NumberOfDOFs() << endl;
}

void irtkPadding(irtkGreyImage *image, irtkGreyPixel padding, irtkBSplineFreeFormTransformationPeriodic *ffd, int numberOfImages, double *time)
{
    int i, j, k, x, y, z, i1, j1, k1, i2, j2, k2, ok, index;
    int t, n, number;
    double t1, t2, tt;

    number = 0;

    // Loop over control points
    for (t = 0; t < ffd->GetT(); t++) {
        for (z = 0; z < ffd->GetZ(); z++) {
            for (y = 0; y < ffd->GetY(); y++) {
                for (x = 0; x < ffd->GetX(); x++) {
                    ok = false;
                    index = ffd->LatticeToIndex(x, y, z, t);
                    // t1, t2 not in lattice coordinates at the moment!!!!!!!
                    ffd->BoundingBoxImage(image, index, i1, j1, k1, i2, j2, k2, t1, t2, 1.0);

                    // loop over all target images
                    for (n = 0; n < numberOfImages; n++) {
                        // transform time point of current target image to lattice coordinates and check for periodicity
                        tt = time[n];
                        // map time to relative time intervall [0,1]
                        while (tt < 0)
                            tt += 1.;
                        while (tt >= 1)
                            tt -= 1.;
                        tt = ffd->TimeToLattice(tt);

                        // check whether time point of current target image is in temporal bounding box (check in lattice coord.)
                        if (( (t1 >= 0) 
                            && (t2 < ffd->GetT()-1) 
                            &&  (tt >= t1) 
                            && (tt<=t2) )
                            || ( (t1 <  0) 
                            && ( (tt <= t2) 
                            || (tt >= t1+ffd->GetT()-1) ) )
                            || ( (t2 >= ffd->GetT()-1) 
                            && ( (tt >= t1) 
                            || (tt <= t2-ffd->GetT()+1) ) ) ) {

                                // Loop over all voxels in the target (reference) volume
                                for (k = k1; k <= k2; k++) {
                                    for (j = j1; j <= j2; j++) {
                                        for (i = i1; i <= i2; i++) {
                                            if (image->GetAsDouble(i, j, k) > padding) {
                                                ok = true;
                                            }
                                        }
                                    }
                                }
                        }
                    }
                    if (ok == false) {
                        ffd->PutStatusCP(x, y, z, t, _Passive, _Passive, _Passive);
                        number++;
                    }
                }
            }
        }
    }
    cout << "Number of CP padded: " << number << " out of " << ffd->NumberOfDOFs() << endl;
}

double GuessResolution(double xsize, double ysize, double zsize)
{
    if ((xsize >= ysize) && (xsize >= zsize)) return xsize;
    if ((ysize >= xsize) && (ysize >= zsize)) return ysize;
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

#endif
