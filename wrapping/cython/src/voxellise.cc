#include "voxellise.h"

void voxellise( vtkPolyData *poly, // input mesh (must be closed)
                irtkGreyImage &image,
                double value ) {

    double pt[3];

    int noOfPts = poly->GetNumberOfPoints();
    for (int i = 0; i < noOfPts; ++i){
        poly->GetPoints()->GetPoint (i, pt);
        image.WorldToImage( pt[0], pt[1], pt[2] );
        poly->GetPoints()->SetPoint(i, pt);
    }    

    vtkSmartPointer<vtkImageData> whiteImage = vtkSmartPointer<vtkImageData>::New();    
    double spacing[3]; 
    spacing[0] = 1.0;
    spacing[1] = 1.0;
    spacing[2] = 1.0;
    whiteImage->SetSpacing(spacing);
    whiteImage->SetExtent( 0, image.GetX() - 1,
                           0, image.GetY() - 1,
                           0, image.GetZ() - 1 );
    double origin[3];
    origin[0] = 0;
    origin[1] = 0;
    origin[2] = 0;
    whiteImage->SetOrigin(origin);
    whiteImage->SetScalarTypeToUnsignedChar();
    whiteImage->AllocateScalars();

    for ( int i = 0; i < whiteImage->GetNumberOfPoints(); i++ )
        whiteImage->GetPointData()->GetScalars()->SetTuple1(i, value);  

    // polygonal data --> image stencil:
    vtkSmartPointer<vtkPolyDataToImageStencil> pol2stenc = 
        vtkSmartPointer<vtkPolyDataToImageStencil>::New();
    pol2stenc->SetInput( poly );
    pol2stenc->SetOutputOrigin(origin);
    pol2stenc->SetOutputSpacing(spacing);
    pol2stenc->SetTolerance(0.0);
    pol2stenc->SetOutputWholeExtent(whiteImage->GetExtent());
    pol2stenc->Update();

    // cut the corresponding white image and set the background:
    vtkSmartPointer<vtkImageStencil> imgstenc = 
        vtkSmartPointer<vtkImageStencil>::New();
    imgstenc->SetInput(whiteImage);
    imgstenc->SetStencil(pol2stenc->GetOutput());
    imgstenc->ReverseStencilOff();
    imgstenc->SetBackgroundValue(0);
    imgstenc->Update();

    vtkSmartPointer<vtkImageData> vtkimageOut =
        vtkSmartPointer<vtkImageData>::New();
    vtkimageOut = imgstenc->GetOutput();
    vtkimageOut->Modified();
    vtkimageOut->Update();

    // Retrieve the output in IRTK
    int n    = image.GetNumberOfVoxels();
    short* ptr1 = image.GetPointerToVoxels();
    unsigned char* ptr2 = (unsigned char*)vtkimageOut->GetScalarPointer();
    for ( int i = 0; i < n; i++ ){
        *ptr1 = *ptr2;
        ptr1++;
        ptr2++;
    }    
}

void create_polydata( double* points,
                      int npoints,
                      int* triangles,
                      int ntriangles,
                      vtkPolyData *poly ) {

    int i;

    // Setup points
    vtkSmartPointer<vtkPoints> vtk_points = vtkSmartPointer<vtkPoints>::New();
    for ( i = 0; i < npoints; i++ )
        vtk_points->InsertNextPoint( points[index(i, 0, npoints, 3)],
                                 points[index(i, 1, npoints, 3)],
                                 points[index(i, 2, npoints, 3)] );
 
    // Setup triangles
    vtkSmartPointer<vtkCellArray> vtk_triangles = vtkSmartPointer<vtkCellArray>::New();
    for ( i = 0; i < ntriangles; i++ ) {
        vtkSmartPointer<vtkTriangle> vtk_triangle = vtkSmartPointer<vtkTriangle>::New();
        vtk_triangle->GetPointIds()->SetId(0, triangles[index(i, 0, ntriangles, 3)]);
        vtk_triangle->GetPointIds()->SetId(1, triangles[index(i, 1, ntriangles, 3)]);
        vtk_triangle->GetPointIds()->SetId(2, triangles[index(i, 2, ntriangles, 3)]);
        vtk_triangles->InsertNextCell(vtk_triangle);
    }

    poly->SetPoints(vtk_points);
    poly->SetPolys(vtk_triangles);
    poly->Update();
}

#include <vtkXMLPolyDataWriter.h>
#include <vtkPolyDataWriter.h>
void _voxellise( double* points,
                 int npoints,
                 int* triangles,
                 int ntriangles,
                 short* img, // irtkGreyImage
                 double* pixelSize,
                 double* xAxis,
                 double* yAxis,
                 double* zAxis,
                 double* origin,
                 int* dim ) {
    vtkSmartPointer<vtkPolyData> poly = vtkSmartPointer<vtkPolyData>::New();

    create_polydata( points,
                     npoints,
                     triangles,
                     ntriangles,
                     poly );

    // Write the file
    vtkSmartPointer<vtkPolyDataWriter> writer =  
        vtkSmartPointer<vtkPolyDataWriter>::New();
    writer->SetFileName("debug.vtk");
    writer->SetInput(poly);
    //writer->SetDataModeToAscii();
    writer->Write();

    irtkGenericImage<short> irtk_image;
    py2irtk<short>( irtk_image,
                    img,
                    pixelSize,
                    xAxis,
                    yAxis,
                    zAxis,
                    origin,
                    dim );
    
    voxellise( poly,
               irtk_image,
               1 );

    irtk2py<short>( irtk_image,
                    img,
                    pixelSize,
                    xAxis,
                    yAxis,
                    zAxis,
                    origin,
                    dim );    
}
