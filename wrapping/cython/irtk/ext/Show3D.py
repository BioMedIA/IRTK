__all__ = [ "render" ]

import vtk
import numpy as np
import sys

from scipy.stats.mstats import mquantiles
import scipy.interpolate

from vtk.util.vtkConstants import *

import matplotlib.pyplot as plt

def numpy2VTK(img,spacing=[1.0,1.0,1.0]):
    # evolved from code from Stou S.,
    # on http://www.siafoo.net/snippet/314
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

    importer.SetDataSpacing( spacing[0], spacing[1], spacing[2])
    importer.SetDataOrigin( 0,0,0 )

    return importer

def volumeRender(img, tf=[],spacing=[1.0,1.0,1.0]):
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
    volMapper.SetInputConnection(importer.GetOutputPort())

    # The property describes how the data will look
    volProperty =  vtk.vtkVolumeProperty()
    volProperty.SetColor(color_tf)
    volProperty.SetScalarOpacity(opacity_tf)
    volProperty.ShadeOn()
    volProperty.SetInterpolationTypeToLinear()

    vol = vtk.vtkVolume()
    vol.SetMapper(volMapper)
    vol.SetProperty(volProperty)
    
    return [vol]


def vtk_basic( actors ):
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
    
    # render
    renWin.Render()
   
    # enable user interface interactor
    iren.Initialize()
    iren.Start()

    
def render( img, tf=None, colors=None, opacity=None ):
    """
    tf: transfert function of the form [[voxel_value, r, g, b, opacity]]
    """
    data = img.rescale().get_data(dtype='uint8')
    if colors is not None:
        tf = plt.cm.get_cmap(colors, 256)(range(256))*255
        tf = np.hstack( (np.arange(256).reshape(256,1), tf[:,:3]) )
        if opacity is not None:
            opacity = np.array(opacity)
            x = opacity[:,0]
            y = opacity[:,1]
            f = scipy.interpolate.interp1d(x, y, kind='linear')
            tf = np.hstack( (tf, f(range(256)).reshape((256,1)) ) )
        else:
            tf = np.hstack( (tf,
                             np.linspace(0,1,len(tf)).reshape(len(tf),1),
                             np.linspace(0,1,len(tf)).reshape(len(tf),1) ) )
    if tf is None:
        q = mquantiles(data.flatten(),[0.5,0.98])
        q[0]=max(q[0],1)
        q[1] = max(q[1],1)
        tf=[[0,0,0,0,0],[q[0],0,0,0,0],[q[1],1,1,1,0.5],[data.max(),1,1,1,1]]

    actor_list = volumeRender(data, tf=tf, spacing=img.header['pixelSize'][:3])
    vtk_basic(actor_list)
