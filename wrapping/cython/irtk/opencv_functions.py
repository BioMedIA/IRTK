from image import Image
import cv2
import numpy as np
import IPython.core.display

__all__ = [ "imshow" ]

def resample2D( self, new_pixelSize, interpolation='linear' ):
    if not isinstance(new_pixelSize, tuple): new_pixelSize = (new_pixelSize,)
    if len(new_pixelSize) == 1:
        new_pixelSize = new_pixelSize*2
    if len(new_pixelSize) > 2:
        raise ValueError( "this function is for XY resampling only" )

    interpolation_cv2 = {
        'nearest' : cv2.INTER_NEAREST, # nearest-neighbor interpolation
        'linear' : cv2.INTER_LINEAR, # bilinear interpolation
        'area' : cv2.INTER_AREA, # resampling using pixel area relation. It
        # may be a preferred method for image decimation, as it gives
        # moire-free results. But when the image is zoomed, it is similar
        # to the INTER_NEAREST method.
        'cubic' : cv2.INTER_CUBIC, # bicubic interpolation over 4x4 pixel neighborhood
        'lanczos' : cv2.INTER_LANCZOS4 # Lanczos interpolation over 8x8
        # pixel neighborhood
        }
        
    data = self.get_data(purpose='cython')
    header = self.get_header()
    # Determine the new dimensions of the image
    dim = header['dim']
    pixelSize = header['pixelSize']
    new_x = int(float(dim[0]) * pixelSize[0] / new_pixelSize[0])
    new_y = int(float(dim[1]) * pixelSize[1] / new_pixelSize[1])
    header['pixelSize'][0] = new_pixelSize[0]
    header['pixelSize'][1] = new_pixelSize[1]
    header['dim'][0] = new_x
    header['dim'][1] = new_y
    new_data = np.zeros( header['dim'][::-1], dtype=data.dtype )
    for t in xrange(data.shape[0]):
        for z in xrange(data.shape[1]):
            new_data[t,z] = cv2.resize( data[t,z], (new_x, new_y),
                                        interpolation=interpolation_cv2[interpolation])
    return Image(new_data, header)

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
        filename = tmp_dir + "/show.png"
        write_only = False
    else:
        write_only = True

    img = img.rescale()
    img = img.thumbnail( seg=seg,
                         overlay=overlay,
                         colors=colors,
                         opacity=opacity,
                         unroll=unroll )
    
    if len(img.shape) == 3:
        # OpenCV uses BGR and not RGB
        img = img[:,:,::-1]

    if write_only:
        cv2.imwrite(filename, img )
        return

    elif run_from_ipython():
        cv2.imwrite(filename, img ) 
        return IPython.core.display.Image(filename=filename)

    else:
        cv2.imshow( "imshow", img )
        cv2.waitKey(0)
