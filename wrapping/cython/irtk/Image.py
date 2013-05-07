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

import numpy as np
import pprint
import copy
import subprocess

from scipy.stats.mstats import mquantiles
import IPython.core.display
import cv2

import _irtk
import Transformation
import Utils

# A temporary directory is created when the module is imported
# and is deleted when the module is exited.
# This is the directory used to save NIFTI files for rview, display and show
import tempfile
import shutil
import os
tmp_dir_location = "/tmp/irtk.Image/"
if not os.path.exists( tmp_dir_location ):
    os.makedirs( tmp_dir_location )

try:
    tmp_dir = tempfile.mkdtemp(dir=tmp_dir_location)
except Exception:
    tmp_dir = tempfile.mkdtemp(dir="./")

import atexit
def _cleanup():
    shutil.rmtree(tmp_dir)
atexit.register(_cleanup)

# function to check if ipython is running
def run_from_ipython():
    try:
        __IPYTHON__
        return True
    except NameError:
        return False

class Image(np.ndarray):

    def __new__(cls,
                img=None,
                header=None):
        if img is None:
            img = np.zeros(1)
        else:
            img = np.array( img, copy=True, order='C' )
            
        if header is None:
            dim = img.shape
            if len(dim) < 4:
                dim = list(dim) + [1]*(4-len(dim))
            header = {
                'pixelSize' : np.ones(4,dtype='float64'),
                'orientation' : (
                    np.array([1,0,0],dtype='float64'),
                    np.array([0,1,0],dtype='float64'),
                    np.array([0,0,1],dtype='float64')),
                'origin' : np.zeros(4,dtype='float64'),
                'dim' : np.array(dim,dtype='int32')
                }

        img = np.squeeze(img)
        obj = np.asarray(img).view(cls)
        
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
        
        data = np.ndarray.__getitem__(self.view(np.ndarray), index).copy('C')

        # Calculate position of first voxel in roi in original image
        o1 = self.ImageToWorld( o1 )

        # Calculate position of first voxel in roi in new image
        tmp_img = Image(header=header)
        o2 = tmp_img.ImageToWorld( np.zeros(3) )

        # Shift origin of new image accordingly
        header['origin'][:3] += o1 - o2
  
        return Image(data,header)

    # def __getslice__(self, i,j):
    #     print 'In __getslice__:'
    #     print i,j
    #     return self.__getitem__(slice(i,j,None))

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
    
    def ImageToWorld( self, pt ):
        """
        BEWARE: images are indexed as img[t,z,y,x],
        so for instance:
        tx,ty,tz = img.ImageToWorld( [(img.shape[2]-1)/2,
                              (img.shape[1]-1)/2,
                              (img.shape[0]-1)/2] )
        """
        #print pt.shape
        pt = np.array(pt)
        if len(pt.shape) == 1:
            tmp_pt = np.hstack((pt,[1])).astype('float64')
            return np.dot( self.I2W, tmp_pt )[:3]
        else:
            #tmp_pt = np.hstack((pt,[[1]]*pt.shape[0])).astype('float64')
            return np.array( map( self.ImageToWorld, pt) )
            # return np.dot( self.I2W, tmp_pt )[:,:3]
    
    def WorldToImage( self, pt ):
        tmp_pt = np.hstack((pt,[1])).astype('float64')
        return np.dot( self.W2I, tmp_pt )[:3]
    
    def get_header(self):
        return copy.deepcopy(self.header)
    
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

    def resample( self, pixelSize=[1,1,1], interpolation='linear', gaussian_parameter=1.0 ):
        pixelSize = np.asarray(pixelSize, dtype='float64')
        new_img, new_header = _irtk.resample( self.get_data('float32', 'cython'),
                                              self.get_header(),
                                              pixelSize,
                                              interpolation,
                                              gaussian_parameter )
        return Image( new_img, new_header )

    def register( self, target, transformation=None ):
        return Transformation.registration_rigid( self, target, transformation )

    def transform( self,
                   transformation=Transformation.RigidTransformation(),
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

    def thumbnail( self, step=2 ):
        if len(self.shape) == 3:
            data = self.rescale().get_data()
            shape = np.array(data.shape).max()
            output = np.ones((shape,shape*3+step*(3-1)))*255

            offset1 = (shape - data.shape[1])/2
            offset2 = (shape - data.shape[2])/2
            output[offset1:offset1+data.shape[1],
                   offset2:offset2+data.shape[2]] = data[data.shape[0]/2,:,:]
    
            offset1 = (shape - data.shape[0])/2
            offset2 = shape + step + (shape - data.shape[2])/2
            output[offset1:offset1+data.shape[0],
                   offset2:offset2+data.shape[2]] = data[:,data.shape[1]/2,:]
    
            offset1 = (shape - data.shape[0])/2
            offset2 = 2*shape + 2*step + (shape - data.shape[1])/2
            output[offset1:offset1+data.shape[0],
                   offset2:offset2+data.shape[1]] = data[:,:,data.shape[2]/2]

            return output
        else:
            raise "Wrong number of dimensions for thumbnail: " + str(len(self.shape)) 

    # We need to override certain functions so that they do not return an image
    #https://github.com/numpy/numpy/blob/master/numpy/matrixlib/defmatrix.py
    # def max(self, **kwargs):
    #     return self.view(np.ndarray).max(**kwargs)
    # def min(self, **kwargs):
    #     return self.view(np.ndarray).min(**kwargs)
    # def sum(self, **kwargs):
    #     return self.view(np.ndarray).sum(**kwargs)

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
    if dtype is None:
        dtype = original_dtype
    img = reader[dtype]( filename, header )
    return Image(img,header)

def imwrite( filename, img ):
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
    writer[img.dtype.name]( filename, data, img.header )
    return

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
            colors = Utils.get_colormap( colors )
        if seg is not None and colors is None:
            colors = Utils.random_colormap(seg.max())

        # now, seg == overlay
        if overlay is not None:
            seg = overlay
    
        filename = tmp_dir + "/rview.nii"
        imwrite( filename, img.astype('int16') ) # rview reads in "short" anyway
        args = [ binary, filename ]
        if seg is not None:
            seg_filename = tmp_dir + "/rview_seg.nii"
            colors_filename = tmp_dir + "/rview_colors.txt"
            imwrite( seg_filename, seg.astype('int16') )
            args.extend( ["-seg", seg_filename])
            Utils.colormap( colors, colors_filename )
            args.extend( ["-lut", colors_filename])
        if overlay is not None:
            args.append("-labels")

    # launch rview as an independant process
    subprocess.Popen( args )

def display( img, seg=None, overlay=None, colors=None ):
    rview(img, seg, overlay, colors, binary="display")

def imshow( img ):
    """
    Display a 2D image
    """
    filename = tmp_dir + "/show.png"
    if len(img.shape) == 2:
        cv2.imwrite(filename, img)
    elif len(img.shape) == 3:
        cv2.imwrite(filename, img.thumbnail())
    else:
        print "Invalid image shape: ", img.shape
        return
    return IPython.core.display.Image(filename=filename) 
    

