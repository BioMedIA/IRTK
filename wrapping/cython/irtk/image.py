from __future__ import division

"""
This module is a wrapper around IRTK (http://www.doc.ic.ac.uk/~dr/software/)
written using Cython and some automated scripts to simulate templated code.

iPython tricks:
# automatically reload modified modules:
%load_ext autoreload
%autoreload 2

Using Notebook:
ipython notebook --pylab inline

Using Qtconsole:
ipython qtconsole --pylab=inline

"""

__all__ = [ "Image",
            "imread",
            "imwrite",
            "imshow",
            "zeros",
            "ones",
            "rview",
            "display",
            "new_header",
            "PointsToImage",
            "drawSphere",
            "sphere",
            "fix_orientation" ]

import numpy as np
import pprint
import copy
import subprocess

from scipy.stats.mstats import mquantiles

import _irtk
import registration
import utils

import IPython.core.display

import tempfile
import os

garbage = []
import atexit
def _cleanup():
    for f in garbage:
        os.remove(f)
atexit.register(_cleanup)

def new_header( pixelSize=[1,1,1,1],
                orientation=np.eye( 3, dtype='float64'),
                origin=[0,0,0,0],
                dim=None ):
    if dim is None:
        raise ValueError( "at least dim must be specified" )
    return { 'pixelSize' : np.array( pixelSize, dtype='float64'),
             'orientation' : orientation,
             'origin' : np.array( origin, dtype='float64'),
             'dim' : np.array( dim, dtype='int32') }

class Image(np.ndarray):

    def __new__(cls,
                img=None,
                header=None,
                copy=True):
        if img is None:
            img = np.zeros(1)
        # else:
        #     img = np.array( img, copy=copy, order='C' )
            
        if header is None:
            dim = img.shape
            if len(dim) < 4:
                dim = list(dim)[::-1] + [1]*(4-len(dim))
            header = {
                'pixelSize' : np.ones(4,dtype='float64'),
                'orientation' : np.eye(3,dtype='float64'),                   
                'origin' : np.zeros(4,dtype='float64'),
                'dim' : np.array(dim,dtype='int32')
                }

        img = np.squeeze(img)
        obj = np.asarray(img).view(cls)

        for i in xrange(4):
            if header['dim'][i] == 0:
                header['dim'][i] = 1
        
        obj.header = header
        obj.I2W = np.eye(4,dtype='float64')
        obj.W2I = np.eye(4,dtype='float64')

        obj.__updateMatrix()
        
        return obj

    def __repr__(self):
        text = "Image data:\n" + repr(self.view(np.ndarray)) + "\n\n"
        text += "Image shape: " + repr(self.shape) + "\n\n"
        if self.header is not None:
            text += "Header:\n" + pprint.pformat(self.header,indent=4)
        else:
            text += "No header information\n"
        return text

    def __str__(self):
        if len(self.shape) == 0:
            return str(self.view(np.ndarray))
        else:
            return self.__repr__()

    def __array_finalize__(self, parent):
        # print 'In __array_finalize__:'
        # print self.shape, parent.shape
        # print '   self is %s' % repr(self)
        # print '   obj is %s' % repr(obj)
        if parent is None: return
        if self.shape != parent.shape:
            self.header = None
            self.I2W = None
            self.W2I = None
        else:
            self.header = getattr(parent, 'header', None)
            self.I2W = getattr(parent, 'I2W', None)
            self.W2I = getattr(parent, 'W2I', None)

    # def __array_wrap__(self, out_arr, context=None):
    #     print 'In __array_wrap__:'
    #     print self.shape, out_arr.shape
    #     # print '   self is %s' % repr(self)
    #     # print '   arr is %s' % repr(out_arr)
    #     #print context
    #     # then just call the parent
    #     return np.ndarray.__array_wrap__(self, out_arr, context)

    def __getitem__(self, index):
        # TODO: time shift
        #print 'In __getitem__:', index
        # Fix index, handling ellipsis and incomplete slices.
        if not isinstance(index, tuple): index = (index,)
        # for i in index:
        #     print type(i),i 
            
        shift = []
        for i,slice_ in enumerate(index):
            if slice_ is Ellipsis:
                raise "use of Ellipsis not implemented"
            elif isinstance(slice_, slice):
                if slice_.start >= 0:
                    shift.append((i,slice_.start,slice_.stop))
            elif isinstance(slice_, (int, long)):
                shift.append((i,slice_,slice_+1))
            else:
                # result is no longer and Image
                # raise "wrong index:", index
                return np.ndarray.__getitem__(self.view(np.ndarray), index)
                
        header = self.get_header()
        mapping = []
        axis_names = ['x','y','z','t']
        axis = {'x':0,'y':1,'z':2,'t':3}
        for i in xrange(4):
            if header['dim'][i] > 1:
                mapping.append(axis_names[i])

        mapping = mapping[::-1]
        o1 = np.zeros(3,dtype='float64')
        for i, start, stop in shift:

            if stop < start:
                raise "negative slices not allowed"
            header['dim'][axis[mapping[i]]] = stop - start
            if mapping[i] == 'x':
                o1[0] = float(start)
            if mapping[i] == 'y':
                o1[1] = float(start)
            if mapping[i] == 'z':
                o1[2] = float(start)
        
        data = np.ndarray.__getitem__(self.view(np.ndarray), index)#.copy('C')

        # Calculate position of first voxel in roi in original image
        o1 = self.ImageToWorld( o1 )

        # Calculate position of first voxel in roi in new image
        tmp_img = Image(header=header)
        o2 = tmp_img.ImageToWorld( np.zeros(3) )

        # Shift origin of new image accordingly
        header['origin'][:3] += o1 - o2
  
        return Image(data,header, copy=False)

    def __reduce__(self):
        """
        Required for pickling/unpickling, which is used for instance
        in joblib Parallel.
        An example implementation can be found in numpy/ma/core.py .
        """        
        return ( _ImageReconstruct,
                 (self.get_data(), self.get_header()) )

    def __updateMatrix(self):
        #  Update image to world coordinate system matrix
        self.I2W = np.eye(4,dtype='float64')
        
        xAxis,yAxis,zAxis = self.header['orientation']
        self.I2W[0:3,0] = xAxis
        self.I2W[0:3,1] = yAxis
        self.I2W[0:3,2] = zAxis

        x,y,z,t = map(float,self.header['dim'])
        tmp1 = np.eye(4,dtype='float64')
        tmp1[0:3,3] = [-(x-1)/2,-(y-1)/2,-(z-1)/2]

        dx,dy,dz,dt = self.header['pixelSize']
        tmp2 = np.eye(4,dtype='float64')*[dx,dy,dz,1]

        tmp3 = np.eye(4,dtype='float64')
        tmp3[0:3,3] = self.origin()

        self.I2W = np.dot(tmp3, np.dot(self.I2W, np.dot(tmp2, tmp1)))

        # Update world to image coordinate system matrix
        self.W2I = np.eye(4,dtype='float64')

        self.W2I[0,0:3] = xAxis
        self.W2I[1,0:3] = yAxis
        self.W2I[2,0:3] = zAxis

        tmp1 = np.eye(4,dtype='float64')
        tmp1[0:3,3] = [(x-1)/2,(y-1)/2,(z-1)/2]

        tmp2 = np.eye(4,dtype='float64')*[1/dx,1/dy,1/dz,1]

        tmp3 = np.eye(4,dtype='float64')
        tmp3[0:3,3] = - self.origin()

        self.W2I = np.dot(tmp1, np.dot(tmp2, np.dot(self.W2I, tmp3)))

    def origin(self):
        return self.header['origin'][:3]

    def nb_pixels(self):
        dim = self.header['dim']
        return dim[0]*dim[1]*dim[2]*dim[3]
    
    def ImageToWorld( self, pt ):
        """
        BEWARE: images are indexed as img[t,z,y,x],
        so for instance:
        tx,ty,tz = img.ImageToWorld( [(img.shape[2]-1)/2,
                              (img.shape[1]-1)/2,
                              (img.shape[0]-1)/2] )
        """
        tmp_pt = np.array( pt, dtype='float64' ).copy()
        if len(tmp_pt.shape) == 1:
            tmp_pt = np.hstack((tmp_pt,[1])).astype('float64')
            return np.dot( self.I2W, tmp_pt )[:3]
        else:
            tmp_pt = _irtk.transform_points( self.I2W, tmp_pt )
            return tmp_pt            
            #tmp_pt = np.hstack((pt,[[1]]*pt.shape[0])).astype('float64')
            #return np.transpose( np.dot( self.I2W, np.transpose(tmp_pt) ) )[:,:3]
            #return np.array( map( self.ImageToWorld, pt) )
            # return np.dot( self.I2W, tmp_pt )[:,:3]
    
    def WorldToImage( self, pt ):
        """
        Returns integers
        """
        tmp_pt = np.array( pt, dtype='float64' ).copy()
        if len(tmp_pt.shape) == 1:
            tmp_pt = np.hstack((tmp_pt,[1])).astype('float64')
            return np.rint( np.dot( self.W2I, tmp_pt )[:3] ).astype('int32')
        else:
            tmp_pt = np.rint( _irtk.transform_points( self.W2I, tmp_pt ) ).astype('int32')
            return tmp_pt
            #tmp_pt = np.hstack((pt,[[1]]*pt.shape[0])).astype('float64')
            #return np.rint( np.transpose( np.dot( self.W2I, np.transpose(tmp_pt) ) )[:,:3])
            #return np.array( map( self.WorldToImage, pt) )
                 
    def get_header(self):
        return copy.deepcopy(self.header)
    def set_header( self, new_header ):
        self.header = copy.deepcopy(new_header)
        self.__updateMatrix()
    def get_data(self, dtype=None, purpose="python" ):
        if dtype is None:
            dtype = self.dtype
        if purpose == "cython":
            shape = self.header['dim'][::-1]
            data = np.reshape( self.view(np.ndarray),
                               shape ).astype( dtype,
                                               copy=True,
                                               order='C')
        else:
            data = np.array( self, dtype=dtype, copy=True)
            
        return data

    def copy( self, dtype=None ):
        return Image( self.get_data(dtype), self.get_header() )

    def saturate( self, q0=0.01, q1=0.99 ):
        data = self.get_data('float32')
        if q0 is None:
            q0 = 0
        if q1 is None:
            q1 = 1
        q = mquantiles(data.flatten(),[q0,q1])
        data[data<q[0]] = q[0]
        data[data>q[1]] = q[1]
        return Image( data, self.header )
    
    def rescale( self, min=0, max=255 ):
        """ Stretch contrast """
        data = self.get_data('float32')
        data -= data.min()
        data /= data.max()
        data = data * (max - min) + min
        return Image( data, self.header )

    def resample( self,
                  pixelSize=[1,1,1],
                  interpolation='linear',
                  gaussian_parameter=1.0 ):
        if isinstance(pixelSize,tuple):
            pixelSize = list(pixelSize)
        if not isinstance(pixelSize,list):
            pixelSize = [pixelSize]*3
        pixelSize = np.asarray(pixelSize, dtype='float64')
        new_img, new_header = _irtk.resample( self.get_data('float32', 'cython'),
                                              self.get_header(),
                                              pixelSize,
                                              interpolation,
                                              gaussian_parameter )
        return Image( new_img, new_header )

    def resample2D( self,
                    pixelSize=[1,1],
                    interpolation='linear',
                    gaussian_parameter=1.0 ):
        # FIXME: do this in C++
        if isinstance(pixelSize,tuple):
            pixelSize = list(pixelSize)
        if not isinstance(pixelSize, list):
            pixelSize = [pixelSize]
        if len(pixelSize) == 1:
            pixelSize = pixelSize*2
        if len(pixelSize) > 2:
            raise ValueError( "this function is for XY resampling only" )
        pixelSize += [1]
        
        header = self.get_header()
        # Determine the new dimensions of the image
        dim = header['dim']
        new_x = int(float(dim[0]) * header['pixelSize'][0] / pixelSize[0])
        new_y = int(float(dim[1]) * header['pixelSize'][1] / pixelSize[1])
        header['pixelSize'][0] = pixelSize[0]
        header['pixelSize'][1] = pixelSize[1]
        header['dim'][0] = new_x
        header['dim'][1] = new_y
        res = zeros( header, dtype='float32' )
        if len(res.shape) == 4:
            for t in xrange(self.shape[0]):
                for z in xrange(self.shape[1]):
                    res[t,z] = self[t,z].resample( pixelSize,
                                                    interpolation=interpolation,
                                                    gaussian_parameter=gaussian_parameter)
        elif len(res.shape) == 3:
            for z in xrange(self.shape[0]):
                res[z] = self[z].resample( pixelSize,
                                            interpolation=interpolation,
                                            gaussian_parameter=gaussian_parameter)
        else:
            res = self.resample( pixelSize,
                                  interpolation=interpolation,
                                  gaussian_parameter=gaussian_parameter)
        return res    
    
    def register( self, target, transformation=None ):
        return registration.registration_rigid( self, target, transformation )

    def transform( self,
                   transformation=registration.RigidTransformation(),
                   target=None,
                   interpolation='linear',
                   gaussian_parameter=1.0):
        return transformation.apply( self,
                                     target,
                                     interpolation,
                                     gaussian_parameter )

    def mask( self, mask, threshold=0  ):
        # transform mask to current image coordinate system
        tmp_mask = mask.transform( target=self )
        new_img = self.copy()
        new_img[ tmp_mask <= threshold ] = 0
        return new_img

    def reset_header( self, header=None ):
        return Image( self.view(np.ndarray), header )

    def thumbnail( self,
                   seg=None,
                   overlay=None,
                   colors=None,
                   opacity=0.5,
                   step=2,
                   unroll=False ):
        img = self.copy('int')
        if overlay is not None and seg is not None:
            print "You cannot specify both seg and overlay"
            return
        if seg is not None:
            seg = seg.astype('int')
        if overlay is not None and colors is None:
            colors = 'jet' # default colormap
        if colors is not None and seg is None and overlay is None:
            overlay = img.copy() # we want to see img in colors
            img = zeros(img.get_header(), dtype='uint8')
            opacity = 1.0
        if isinstance( colors, str ):
            colors = utils.get_colormap( colors )
        if seg is not None and colors is None:
            colors = utils.random_colormap(seg.max())

        # now, seg == overlay
        if overlay is not None:
            seg = overlay

        if len(img.shape) == 2:
            if seg is not None:
                data = img.get_data('uint8')
                data = data.reshape(data.shape[0],
                                    data.shape[1],
                                    1)
                data = np.concatenate( [ data,
                                         data,
                                         data ], axis=2 )
                rgb_overlay = utils.remap( seg, colors )
                op = (1-opacity) * (seg != 0).astype('float') + (seg == 0).astype('float')
                op = op.reshape(op.shape[0],
                                op.shape[1],
                                1)
                op = np.concatenate( [ op,
                                       op,
                                       op ], axis=2 )
                img = op * data  + opacity*rgb_overlay

            return img.astype('uint8')
        
        elif len(img.shape) == 3:
            if not unroll:
                shape = np.array(img.shape).max()
                if seg is None:
                    output = np.ones((shape,shape*3+step*(3-1)))*255
                else:
                    output = np.ones((shape,shape*3+step*(3-1),3))*255

                offset1 = (shape - img.shape[1])/2
                offset2 = (shape - img.shape[2])/2
                if seg is None:
                    tmp_img = img[img.shape[0]/2,:,:]
                else:
                    tmp_img = Image(img[img.shape[0]/2,:,:]).thumbnail( seg[img.shape[0]/2,:,:], colors=colors, opacity=opacity )
                output[offset1:offset1+img.shape[1],
                       offset2:offset2+img.shape[2]] = tmp_img

                offset1 = (shape - img.shape[0])/2
                offset2 = shape + step + (shape - img.shape[2])/2
                if seg is None:
                    tmp_img = img[:,img.shape[1]/2,:]
                else:
                    tmp_img = Image(img[:,img.shape[1]/2,:]).thumbnail( seg[:,img.shape[1]/2,:], colors=colors, opacity=opacity)
                output[offset1:offset1+img.shape[0],
                       offset2:offset2+img.shape[2]] = tmp_img

                offset1 = (shape - img.shape[0])/2
                offset2 = 2*shape + 2*step + (shape - img.shape[1])/2
                if seg is None:
                    tmp_img = img[:,:,img.shape[2]/2]
                else:
                    tmp_img = Image(img[:,:,img.shape[2]/2]).thumbnail( seg[:,:,img.shape[2]/2], colors=colors, opacity=opacity )
                output[offset1:offset1+img.shape[0],
                       offset2:offset2+img.shape[1]] = tmp_img

                return output.astype('uint8')
            else: # unroll is True
                if seg is None:
                    output = np.ones( ( self.shape[1],
                                        self.shape[2]*self.shape[0]+2*(self.shape[0]-1) )
                                      ) * 255
                else:
                    output = np.ones( ( self.shape[1],
                                        self.shape[2]*self.shape[0]+2*(self.shape[0]-1),
                                        3 )
                                      ) * 255
                for k in xrange(self.shape[0]):
                    if seg is None:
                        tmp_img = img[k,:,:]
                    else:
                        tmp_img = Image(img[k,:,:]).thumbnail( seg[k,:,:], colors=colors, opacity=opacity )
                    output[:,
                        k*self.shape[2]+2*k:(k+1)*self.shape[2]+2*k] = tmp_img
                return output.astype('uint8')
        else:
            raise "Wrong number of dimensions for thumbnail: " + str(len(self.shape)) 

    # We need to override certain functions so that they do not return an image
    # https://github.com/numpy/numpy/blob/master/numpy/matrixlib/defmatrix.py
    # def max(self, **kwargs):
    #     return self.view(np.ndarray).max(**kwargs)
    # def min(self, **kwargs):
    #     return self.view(np.ndarray).min(**kwargs)
    # def sum(self, **kwargs):
    #     return self.view(np.ndarray).sum(**kwargs)

    # def unroll( self ):
    #     output = np.ones( ( self.shape[1],
    #                         self.shape[2]*self.shape[0]+2*(self.shape[0]-1) )
    #                       ) * 255
    #     for k in xrange(self.shape[0]):
    #         output[:,
    #             k*self.shape[2]+2*k:(k+1)*self.shape[2]+2*k] = self[k,:,:]

    #     return output


    def bbox( self, crop=False, world=False ):
        pts = np.transpose(np.nonzero(self>0))
        if world:
            pts = self.ImageToWorld( pts[:,::-1] )
            x_min, y_min, z_min = pts.min(axis=0)
            x_max, y_max, z_max = pts.max(axis=0)
            return (x_min, y_min, z_min, x_max, y_max, z_max )
        if len(self.shape) == 2:
            y_min, x_min = pts.min(axis=0)
            y_max, x_max = pts.max(axis=0)
            if not crop:
                return (x_min, y_min, x_max, y_max )
            else:
                return self[y_min:y_max+1,
                            x_min:x_max+1]        
        if len(self.shape) == 3:
            z_min, y_min, x_min = pts.min(axis=0)
            z_max, y_max, x_max = pts.max(axis=0)
            if not crop:
                return (x_min, y_min, z_min, x_max, y_max, z_max )
            else:
                return self[z_min:z_max+1,
                            y_min:y_max+1,
                            x_min:x_max+1]
        if len(self.shape) == 4:
            t_min, z_min, y_min, x_min = pts.min(axis=0)
            t_max, z_max, y_max, x_max = pts.max(axis=0)
            if not crop:
                return (x_min, y_min, z_min, t_min, x_max, y_max, z_max, t_max )
            else:
                return self[t_min:t_max+1,
                            z_min:z_max+1,
                            y_min:y_max+1,
                            x_min:x_max+1]            

    def order( self ):
          if np.linalg.det(self.I2W[:3,:3]) < 0:
              return "radiological"
          else:
              return "neurological"
    def neurological( self ):
        return self.order() == "neurological"
    def radiological( self ):
        return self.order() == "radiological"
    def orientation( self ):
        return _irtk.orientation( self.get_header() )



    
def imread( filename, dtype=None ):
    reader = {
        "int8" : _irtk.imread_char,
        "uint8" : _irtk.imread_uchar,
        "int16" : _irtk.imread_short,
        "uint16" : _irtk.imread_ushort,
        "int32" : _irtk.imread_int,
        "uint32" : _irtk.imread_uint,
        "float32" : _irtk.imread_float,
        "float64" : _irtk.imread_double,
        }

    header, original_dtype = _irtk.get_header( filename )
    if header is None:
        return False
    if dtype is None:
        dtype = original_dtype
    img = reader[dtype]( filename, header )
    if img is False:
        return False
    else:
        return Image(img,header)

def imwrite( filename, img ):
    if not isinstance(img,Image):
        img = Image(img)
    if img.dtype.name == "bool":
        img = img.astype("uint8")
    writer = {
        "int8" : _irtk.imwrite_char,
        "uint8" : _irtk.imwrite_uchar,
        "int16" : _irtk.imwrite_short,
        "uint16" : _irtk.imwrite_ushort,
        "int32" : _irtk.imwrite_int,
        "uint32" : _irtk.imwrite_uint,
        "float32" : _irtk.imwrite_float,
        "float64" : _irtk.imwrite_double,
        }

    shape = img.header['dim'][::-1]
    data = np.reshape( img.view(np.ndarray), shape ).copy('C')
    success = writer[img.dtype.name]( filename, data, img.header )
    return success

def zeros( header, dtype='float32' ):
    return Image( np.zeros( header['dim'][::-1], dtype=dtype ), header )

def ones( header, dtype='float32' ):
    return Image( np.ones( header['dim'][::-1], dtype=dtype ), header )

def rview( img, seg=None, overlay=None, colors=None, binary="rview" ):
    if isinstance( img, Image ):
        if overlay is not None and seg is not None:
            print "You cannot specify both seg and overlay"
            return
        if overlay is not None and colors is None:
            colors = 'jet' # default colormap
        if colors is not None and overlay is None:
            overlay = img.rescale() # we want to see img in colors
            img = zeros(img.get_header(), dtype='int16')
        if isinstance( colors, str ):
            colors = utils.get_colormap( colors )
        if seg is not None and colors is None:
            colors = utils.random_colormap(seg.max())

        # now, seg == overlay
        if overlay is not None:
            seg = overlay

        if seg is not None and not isinstance(seg, Image):
            seg = Image( seg, img.get_header() )
    
        filename = tmp_dir + "/rview.nii"
        imwrite( filename, img.astype('int16') ) # rview reads in "short" anyway
        args = [ binary, filename ]
        if seg is not None:
            seg_filename = tmp_dir + "/rview_seg.nii"
            colors_filename = tmp_dir + "/rview_colors.txt"
            imwrite( seg_filename, seg.astype('int16') )
            args.extend( ["-seg", seg_filename])
            utils.colormap( colors, colors_filename )
            args.extend( ["-lut", colors_filename])
        if overlay is not None:
            args.append("-labels")

    # launch rview as an independant process
    subprocess.Popen( args )

def display( img, seg=None, overlay=None, colors=None ):
    rview(img, seg, overlay, colors, binary="display")


    

def write_list( img_list ):
    n = len(img_list)
    nb_pixels = 0
    for img in img_list:
        nb_pixels += img.nb_pixels()

    pixelData = np.zeros( nb_pixels, dtype='float32' )
    pixelSize = np.zeros( 4*n, dtype='float64' )
    xAxis = np.zeros( 3*n, dtype='float64' )
    yAxis = np.zeros( 3*n, dtype='float64' )
    zAxis = np.zeros( 3*n, dtype='float64' )
    origin = np.zeros( 4*n, dtype='float64' )
    dim = np.zeros( 4*n, dtype='int32' )
    
    offset = 0
    for i, img in enumerate(img_list):
        pixelData[offset:offset+img.nb_pixels()] = img.flatten()
        offset += img.nb_pixels()
        pixelSize[4*i:4*(i+1)] = img.header['pixelSize']
        xAxis[3*i:3*(i+1)] = img.header['orientation'][0]
        yAxis[3*i:3*(i+1)] = img.header['orientation'][1]
        zAxis[3*i:3*(i+1)] = img.header['orientation'][2]
        origin[4*i:4*(i+1)] = img.header['origin']
        dim[4*i:4*(i+1)] = img.header['dim']

    print pixelData.shape

    _irtk.write_list( pixelData, pixelSize, xAxis, yAxis, zAxis, origin, dim, n )

def PointsToImage( pts, header ):
    if isinstance(header,Image):
        header = header.get_header()
    pts = np.array(pts, dtype="float64")
    img = _irtk.points_to_image( pts, header )
    return Image( img, header )

def drawSphere( img, (x,y,z), rx, ry=None, rz=None ):
    if ry is None:
        ry = rx
    if rz is None:
        rz = rx
    _irtk.drawSphere( img, x, y, z, rx, ry, rz )

def sphere( (x, y, z), r, header ):
    img = zeros( header, dtype='uint8' )
    x0, y0, z0 = img.WorldToImage( (x, y, z) )
    rx = int(round( r / header['pixelSize'][0] ))
    ry = int(round( r / header['pixelSize'][1] ))
    rz = int(round( r / header['pixelSize'][2] ))
    print "Sphere:",  x0, y0, z0, rx, ry, rz
    _irtk.drawSphere( img, x0, y0, z0, rx, ry, rz )
    return img #Image( img, header )

def _ImageReconstruct( data, header ):
    """Internal function that builds a new Image from the
    information stored in a pickle.
    """
    return Image( data, header )

def fix_orientation( img1, img2, verbose=False ):
    order1 = img1.order()
    order2 = img2.order()
    if verbose:
        print "Image1 is in " + order1 + " order"
        print "Image2 is in " + order2 + " order"
    if order1 == order2:
        if verbose:
            print "nothing to fix"
        return img2
    i1, j1, k1 = img1.orientation()
    i2, j2, k2 = img2.orientation()

    header = img2.get_header()
    data = img2.get_data()

    if i1 != i2:
        data = data[:,:,::-1]
        header['orientation'][0] *= -1
        if verbose:
            print "flipping X axis"

    if j1 != j2:
        data = data[:,::-1,:]
        header['orientation'][1] *= -1
        if verbose:
            print "flipping Y axis"

    if k1 != k2:
        data = data[::-1,:,:]
        header['orientation'][2] *= -1
        if verbose:
            print "flipping Z axis"

    return Image( data.copy(), header )

# function to check if ipython is running
def run_from_ipython():
    try:
        __IPYTHON__
        return True
    except NameError:
        return False

def rgb(img):
    new_img = np.zeros((3,img.shape[0],img.shape[1]),dtype="uint8")
    for i in xrange(3):
        new_img[i] = img[:,:,i]
    return new_img

def imshow( img,
            seg=None,
            overlay=None,
            colors=None,
            opacity=0.5,
            filename=None,
            unroll=False ):
    """
    Display a 2D image
    """
    if len(img.shape) not in [2,3]:
        print "Invalid image shape: ", img.shape
        return
    
    if filename is None:
        handle,filename = tempfile.mkstemp(suffix=".png")
        garbage.append(filename)
        write_only = False
    else:
        write_only = True

    img = img.rescale()
    img = img.thumbnail( seg=seg,
                         overlay=overlay,
                         colors=colors,
                         opacity=opacity,
                         unroll=unroll ).astype("uint8")

    if len(img.shape) == 3 and img.shape[2] == 3:
        img = rgb(img)

    if write_only:
        return imwrite(filename, img )

    elif run_from_ipython():
        if imwrite(filename, img ):
            return IPython.core.display.Image(filename=filename)

    return False
