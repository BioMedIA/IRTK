"""
VTK shortcut functions

For more about VTK, visit:
    http://www.vtk.org/
    http://www.vtk.org/Wiki/VTK/Examples/Python
"""

import vtk
import sys
import numpy

from vtk.util.colors import peacock

def vtk_point_cloud(points, colors=[], point_size=2):
    """
    Represent a point cloud in VTK
    
    Parameters
    ----------
    points :  numpy array, each row is a point
    colors : list of colors, one per point
    point_size : rendering size for the points
    
    Returns
    -------
    actor : vtkActor representing the point cloud
    """     
    nb = len(points);
    vtk_points = vtk.vtkPoints();
    vtk_verts = vtk.vtkCellArray();
    if colors:
        vtk_colors = vtk.vtkUnsignedCharArray();
        vtk_colors.SetNumberOfComponents(3);
        vtk_colors.SetName( "Colors");
        
    for i in range(0,nb):
        
        p = points[i]
        if len(p) >= 3:
            print "ok",p
            coords = [p[0],p[1],p[2]]
        elif len(p) == 2:
            coords = [p[0],p[1],0]
        elif len(p) == 1:
            coords = [p[0],0,0]
        else:
            print "**ERROR** wrong dimension"
            sys.exit(1)
        
        id = vtk_points.InsertNextPoint( *coords )
        vtk_verts.InsertNextCell(1)
        vtk_verts.InsertCellPoint(id)
        if colors:
            vtk_colors.InsertNextTuple3( *colors[i] )
    
    poly = vtk.vtkPolyData()
    poly.SetPoints(vtk_points)
    poly.SetVerts(vtk_verts)
    if colors:
        poly.GetPointData().SetScalars(vtk_colors)
    poly.Update()

    mapper = vtk.vtkPolyDataMapper()
    mapper.SetInput(poly)

    actor = vtk.vtkActor()
    actor.SetMapper(mapper)
    actor.GetProperty().SetRepresentationToPoints
    actor.GetProperty().SetPointSize( point_size )

    return actor

def vtk_basic( actors, save="", magnification=3 ):
    """
    Create a window, renderer, interactor, add the actors and start the thing
    
    Parameters
    ----------
    actors :  list of vtkActors
    
    Returns
    -------
    nothing
    """     
    
    # create a rendering window and renderer
    ren = vtk.vtkRenderer()
    renWin = vtk.vtkRenderWindow()
    renWin.AddRenderer(ren)
    renWin.SetSize(600,600)
    # ren.SetBackground( 1, 1, 1)
 
    # create a renderwindowinteractor
    iren = vtk.vtkRenderWindowInteractor()
    iren.SetRenderWindow(renWin)

    for a in actors:
        # assign actor to the renderer
        ren.AddActor(a )

    # style = vtk.vtkInteractorStyleTerrain() 
    # iren.SetInteractorStyle( style )
    
    # render
    renWin.Render()
   
    if save:
        if not save.endswith('.png'):
            save += '.png'
        grabber = vtk.vtkWindowToImageFilter()
        grabber.SetInput( renWin )
        grabber.SetMagnification( magnification )
        grabber.Update()
        
        writer = vtk.vtkPNGWriter()
        writer.SetInput( grabber.GetOutput() )
        writer.SetFileName( save )
        writer.Write()

    else:
        # enable user interface interactor
        iren.Initialize()
        iren.Start()

def vtk_Nviews( actors, split='h' ):
    """
    Create a window, an interactor and one renderer per actor
    
    Parameters
    ----------
    actors :  list of vtkActors
    
    Returns
    -------
    nothing
    """    
    N = len(actors)
    
    # create a rendering window and renderers
    renderers = [vtk.vtkRenderer() for i in range(N)]
    renWin = vtk.vtkRenderWindow()
    renWin.SetSize( 600, 600 )
    for i in range(N):
        # split the viewport
        if split == 'h':
            renderers[i].SetViewport(0,float(N-i-1)/N,1,float(N-i)/N)
        else:
            renderers[i].SetViewport(float(N-i-1)/N,0,float(N-i)/N,1)
        renderers[i].SetBackground( 1, 1, 1)
        renderers[i].AddActor( actors[i] )
        renWin.AddRenderer(renderers[i])
 
    # create a renderwindowinteractor
    iren = vtk.vtkRenderWindowInteractor()
    iren.SetRenderWindow(renWin)
    
    #enable user interface interactor
    iren.Initialize()
    renWin.Render()
    iren.Start()

def vtk_show_points( points, colors=[] ):
    """
    Display a point cloud
    
    Parameters
    ----------
    points :  numpy array, each row is a point
    colors : list of colors, one per point
    
    Returns
    -------
    nothing
    """     
    point_cloud = vtk_point_cloud(points,colors)
    vtk_basic( [point_cloud] )


def vtk_colored_graph(points, edges, colors=[], line_width=2):
    """
    Represent a graph in VTK
    
    Parameters
    ----------
    points :  numpy array, each row is a point
    edges : numpy array of edges, each row is of the form
            [ point_1, point_2, distance ]      
    colors : list of colors, one per point
    line_width : rendering size for the lines
    
    Returns
    -------
    actor : vtkActor representing the graph
    """
    nb_points = len(points)    
    vtk_points = vtk.vtkPoints()
    vtk_lines = vtk.vtkCellArray()
    vtk_colors = vtk.vtkUnsignedCharArray()
    vtk_colors.SetNumberOfComponents(3)
    vtk_colors.SetName( "Colors")

    if (len(colors) ==0):
        for i in range(0,len(edges)):
            colors.append((0, 164, 180))
    
    for i in range(0,nb_points):
        
        p = points[i]
        if len(p) >= 3:
            coords = [p[0],p[1],p[2]]
        elif len(p) == 2:
            coords = [p[0],p[1],0]
        elif len(p) == 1:
            coords = [p[0],0,0]
        else:
            print "**ERROR** wrong dimension"
            sys.exit(1)
        
        id = vtk_points.InsertNextPoint( *coords )

    for i in range(0,len(edges)):
        line = vtk.vtkLine()
        line.GetPointIds().SetId(0,edges[i][0])
        line.GetPointIds().SetId(1,edges[i][1])
        vtk_lines.InsertNextCell(line)
        vtk_colors.InsertNextTuple3( *colors[i] )
        
    poly = vtk.vtkPolyData()
    poly.SetPoints(vtk_points)
    poly.SetLines(vtk_lines)
    poly.GetCellData().SetScalars(vtk_colors);
    poly.Update()
    
    cleaner = vtk.vtkCleanPolyData()
    cleaner.SetInput(poly)

    tubes = vtk.vtkTubeFilter()
    tubes.SetInputConnection(cleaner.GetOutputPort())
    tubes.SetRadius(0.1)
    tubes.SetNumberOfSides(6)
    mapEdges = vtk.vtkPolyDataMapper()
    mapEdges.SetInputConnection(tubes.GetOutputPort())
    edgeActor = vtk.vtkActor()
    edgeActor.SetMapper(mapEdges)
    edgeActor.GetProperty().SetSpecularColor(1, 1, 1)
    edgeActor.GetProperty().SetSpecular(0.3)
    edgeActor.GetProperty().SetSpecularPower(20)
    edgeActor.GetProperty().SetAmbient(0.2)
    edgeActor.GetProperty().SetDiffuse(0.8)
    return edgeActor

def vtk_triangles(points, triangles, colors=[]):
    """
    Display triangles in VTK
    
    Parameters
    ----------
    points :  numpy array, each row is a point
    triangle : numpy array of vertices, each row is of the form
            [ point_1, point_2, point_3 ]      
    colors : list of colors, one per triangle
    
    Returns
    -------
    actor : vtkActor representing the triangles
    """    
    nb_points = len(points)
    vtk_points = vtk.vtkPoints()
    vtk_triangles = vtk.vtkCellArray()
    vtk_colors = vtk.vtkUnsignedCharArray()
    vtk_colors.SetNumberOfComponents(3)
    vtk_colors.SetName( "Colors")

    if (len(colors) ==0):
        for i in range(0,nb_points):
            vtk_colors.InsertNextTuple3(0, 164, 180)
    else:
        for i in range(0,nb_points):
            vtk_colors.InsertNextTuple3( *colors[i] )
    
    for i in range(0,nb_points):
        
        p = points[i]
        if len(p) >= 3:
            coords = [p[0],p[1],p[2]]
        elif len(p) == 2:
            coords = [p[0],p[1],0]
        elif len(p) == 1:
            coords = [p[0],0,0]
        else:
            print "**ERROR** wrong dimension"
            sys.exit(1)
        
        id = vtk_points.InsertNextPoint( *coords )

    for i in range(0,len(triangles)):
        triangle = vtk.vtkTriangle()
        triangle.GetPointIds().SetId(0,triangles[i][0])
        triangle.GetPointIds().SetId(1,triangles[i][1])
        triangle.GetPointIds().SetId(2,triangles[i][2])
        vtk_triangles.InsertNextCell(triangle)
        
    poly = vtk.vtkPolyData()
    poly.SetPoints(vtk_points)
    poly.SetPolys(vtk_triangles)
    poly.GetPointData().SetScalars(vtk_colors)
    poly.Update()

    cleaner = vtk.vtkCleanPolyData()
    cleaner.SetInput(poly)

    mapper = vtk.vtkPolyDataMapper()
    mapper.SetInputConnection(cleaner.GetOutputPort())

    actor = vtk.vtkActor()
    actor.SetMapper(mapper)

    return actor    

from vtk.util.vtkConstants import *

def numpy2VTK(img,spacing=[1.0,1.0,1.0]):
    importer = vtk.vtkImageImport()
    
    img_data = img.astype('uint8')
    img_string = img_data.tostring() # type short
    dim = img.shape
    
    importer.CopyImportVoidPointer(img_string, len(img_string))
    importer.SetDataScalarType(VTK_UNSIGNED_CHAR)
    importer.SetNumberOfScalarComponents(1)
    
    extent = importer.GetDataExtent()
    importer.SetDataExtent(extent[0], extent[0] + dim[2] - 1,
                           extent[2], extent[2] + dim[1] - 1,
                           extent[4], extent[4] + dim[0] - 1)
    importer.SetWholeExtent(extent[0], extent[0] + dim[2] - 1,
                            extent[2], extent[2] + dim[1] - 1,
                            extent[4], extent[4] + dim[0] - 1)
    importer.SetDataSpacing( spacing[2], spacing[1], spacing[0])
    importer.SetDataOrigin( 0,0,0 )

    # flip = vtk.vtkImageFlip()
    # flip.SetFilteredAxis(2)
    # flip.SetInput(importer.GetOutput())

    return importer

def volumeRender(img, tf=[],spacing=[1.0,1.0,1.0], box=False):
    importer = numpy2VTK(img,spacing)

    # Transfer Functions
    opacity_tf = vtk.vtkPiecewiseFunction()
    color_tf = vtk.vtkColorTransferFunction()

    if len(tf) == 0:
        tf.append([img.min(),0,0,0,0])
        tf.append([img.max(),1,1,1,1])

    for p in tf:
        color_tf.AddRGBPoint(p[0], p[1], p[2], p[3])
        opacity_tf.AddPoint(p[0], p[4])

    volMapper = vtk.vtkGPUVolumeRayCastMapper()
    # compositeFunction = vtk.vtkVolumeRayCastCompositeFunction()
    # compositeFunction.SetCompositeMethodToInterpolateFirst()
    # volMapper.SetVolumeRayCastFunction(compositeFunction)

    # pix_diag = 5.0
    # volMapper.SetSampleDistance(pix_diag / 5.0)
    volMapper.SetInputConnection(importer.GetOutputPort())

    # The property describes how the data will look
    volProperty =  vtk.vtkVolumeProperty()
    volProperty.SetColor(color_tf)
    volProperty.SetScalarOpacity(opacity_tf)
    volProperty.ShadeOn()
    volProperty.SetInterpolationTypeToLinear()
    #volProperty.SetScalarOpacityUnitDistance(pix_diag)

    vol = vtk.vtkVolume()
    vol.SetMapper(volMapper)
    vol.SetProperty(volProperty)
    
    if not box:
        return [vol]
    else:
        bbox = outline(importer)
        return [vol,bbox]

def marchingCubes(img,spacing=[1.0,1.0,1.0],contours=[]):
    importer = numpy2VTK(img,spacing)

    if len(contours) == 0:
        contours = [[img.max(),1.0,1.0,1.0,1.0]]

    actors = []

    for c in contours:
        mc = vtk.vtkMarchingCubes()
        mc.ComputeScalarsOff()
        mc.ComputeGradientsOff()
        mc.ComputeNormalsOff()
        mc.SetValue( 0, c[0] )
        mc.SetInput( importer.GetOutput())

        # connectivityFilter = vtk.vtkPolyDataConnectivityFilter()
        # connectivityFilter.SetInput(mc.GetOutput())
        # connectivityFilter.ColorRegionsOff()       
        # connectivityFilter.SetExtractionModeToLargestRegion()

        # tris = vtk.vtkTriangleFilter()
        # tris.SetInput(mc.GetOutput())
        # tris.GetOutput().ReleaseDataFlagOn()
        # tris.Update()
        # strip = vtk.vtkStripper()
        # strip.SetInput(tris.GetOutput())
        # strip.GetOutput().ReleaseDataFlagOn()
        # strip.Update()
        
        mapper = vtk.vtkDataSetMapper()
        mapper.SetInput(mc.GetOutput() )
        mapper.ImmediateModeRenderingOn()

#        mapper.SetInput( connectivityFilter.GetOutput() )

        actor = vtk.vtkActor()
        actor.SetMapper( mapper)
        actor.GetProperty().SetColor(c[1],c[2],c[3])
        actor.GetProperty().SetOpacity(c[4])
        actor.GetProperty().SetRepresentationToSurface()

        actors.append(actor)

    return actors

def contours(img,spacing=[1.0,1.0,1.0],contours=[]):
    importer = numpy2VTK(img,spacing)

    if len(contours) == 0:
        contours = [[img.max(),1.0,1.0,1.0,1.0]]

    actors = []

    for c in contours:
        contourExtractor = vtk.vtkContourFilter()
        contourExtractor.SetInputConnection(importer.GetOutputPort())
        contourExtractor.SetValue(0, c[0])
        
        # contourNormals = vtk.vtkPolyDataNormals()
        # contourNormals.SetInputConnection(contourExtractor.GetOutputPort())
        # contourNormals.SetFeatureAngle(60.0)
        
        # contourStripper = vtk.vtkStripper()
        # contourStripper.SetInputConnection(contourNormals.GetOutputPort())

        deci = vtk.vtkDecimatePro()
        deci.SetInputConnection(contourExtractor.GetOutputPort())
        deci.SetTargetReduction(0.99)
        deci.PreserveTopologyOn ()
        # smoother = vtk.vtkSmoothPolyDataFilter()
        # smoother.SetInputConnection(deci.GetOutputPort())
        # smoother.SetNumberOfIterations(50)
        normals = vtk.vtkPolyDataNormals()
        normals.SetInputConnection(deci.GetOutputPort())
        normals.FlipNormalsOn()

        # smoothFilter = vtk.vtkSmoothPolyDataFilter()
        # smoothFilter.SetInputConnection(contourStripper.GetOutputPort())
        # smoothFilter.SetNumberOfIterations(50)
        # smoothFilter.Update()
        
        contourMapper = vtk.vtkPolyDataMapper()
#        contourMapper.SetInputConnection(contourStripper.GetOutputPort())
        contourMapper.SetInputConnection(normals.GetOutputPort())        
        contourMapper.ScalarVisibilityOff()
        
        actor = vtk.vtkActor()
        actor.SetMapper( contourMapper)
        actor.GetProperty().SetColor(c[1],c[2],c[3])
        actor.GetProperty().SetOpacity(c[4])
        actor.GetProperty().SetRepresentationToSurface()        

   # # An outline provides context around the data.
   # outlineData = vtk.vtkOutlineFilter()
   # outlineData.SetInputConnection(v16.GetOutputPort())
   # mapOutline = vtk.vtkPolyDataMapper()
   # mapOutline.SetInputConnection(outlineData.GetOutputPort())
   # outline = vtk.vtkActor()
   # outline.SetMapper(mapOutline)
   # outline.GetProperty().SetColor(0, 0, 0)
        
 

        actors.append(actor)

    return actors

def axes(width=10):
    pts = vtk.vtkPoints()
    pts.InsertNextPoint([0.0, 0.0, 0.0])
    pts.InsertNextPoint([width*10, 0.0, 0.0])
    pts.InsertNextPoint([0.0, width*10, 0.0])
    pts.InsertNextPoint([0.0, 0.0, width*10])

    # http://vtk.org/gitweb?p=VTK.git;a=blob;f=Examples/Annotation/Tcl/textOrigin.tcl
    # Xactor = vtk.vtkTextActor()
    # Xactor.GetTextProperty().SetFontSize( 24 )
    # Xactor.SetPosition( [width*10, 0.0, 0.0] )
    # Xactor.SetInput( "X" )
    # Xactor.GetTextProperty().SetColor( 1.0,0.0,0.0 )

    # Setup the colors array
    colors = vtk.vtkUnsignedCharArray()
    colors.SetNumberOfComponents(3)
    colors.SetName("Colors")
 
    # Add the colors we created to the colors array
    colors.InsertNextTuple3(255, 0, 0)
    colors.InsertNextTuple3(0, 255, 0)
    colors.InsertNextTuple3(0, 0, 255)

    # Create the first line (between Origin and P0)
    line0 = vtk.vtkLine()
    line0.GetPointIds().SetId(0,0) #the second 0 is the index of the Origin in the vtkPoints
    line0.GetPointIds().SetId(1,1) #the second 1 is the index of P0 in the vtkPoints
 
    # Create the second line (between Origin and P1)
    line1 = vtk.vtkLine()
    line1.GetPointIds().SetId(0,0) #the second 0 is the index of the Origin in the vtkPoints
    line1.GetPointIds().SetId(1,2) #2 is the index of P1 in the vtkPoints

    # Create the third line (between Origin and P2)
    line2 = vtk.vtkLine()
    line2.GetPointIds().SetId(0,0) #the second 0 is the index of the Origin in the vtkPoints
    line2.GetPointIds().SetId(1,3) #3 is the index of P2 in the vtkPoints
 
    # Create a cell array to store the lines in and add the lines to it
    lines = vtk.vtkCellArray()
    lines.InsertNextCell(line0)
    lines.InsertNextCell(line1)
    lines.InsertNextCell(line2)

    # Create a polydata to store everything in
    linesPolyData = vtk.vtkPolyData()
 
    # Add the points to the dataset
    linesPolyData.SetPoints(pts)
 
    # Add the lines to the dataset
    linesPolyData.SetLines(lines)
 
    # Color the lines - associate the first component (red) of the
    # colors array with the first component of the cell array (line 0)
    # and the second component (green) of the colors array with the
    # second component of the cell array (line 1)
    linesPolyData.GetCellData().SetScalars(colors)
 
    # Visualize
    mapper = vtk.vtkPolyDataMapper()
    mapper.SetInput(linesPolyData)
 
    actor = vtk.vtkActor()
    actor.SetMapper(mapper)
    actor.GetProperty().SetLineWidth(width)



    return actor#, Xactor

def move(actor, matrix):
    transfo_mat = vtk.vtkMatrix4x4()
    # for i in range(0,4):
    #     for j in range(0,4):
    #         transfo_mat.SetElement(i,j, matrix[i,j])
    for k in range(0,4):
        transfo_mat.SetElement(k,3, matrix[k,3])
    for k in range(0,4):
        transfo_mat.SetElement(k,0, matrix[k,2])
    for k in range(0,4):
        transfo_mat.SetElement(k,1, matrix[k,1])
    for k in range(0,4):
        transfo_mat.SetElement(k,2, matrix[k,0])
        
    print matrix
    print transfo_mat

    actor.SetUserMatrix(transfo_mat)

def outline(source):
    outline_filter = vtk.vtkOutlineFilter()
    outline_filter.SetInput(source.GetOutput())
    mapper = vtk.vtkPolyDataMapper()
    mapper.SetInput(outline_filter.GetOutput())
    actor = vtk.vtkActor()
    actor.SetMapper(mapper)
    return actor
