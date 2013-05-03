import numpy as np
import cv2

from .. import _irtk

predefined_colors = np.array( [[0,0,0],
                               [0,0,255],
                               [0,255,0],
                               [255,0,0]], dtype='uint8' )

# def remap( a, colors=None ):
#     palette, index = np.unique(a, return_inverse = True)
#     if colors is None:
#         colors = np.random.random_integers( 0, 255, len(palette) )
#     return colors[index].reshape(a.shape)

def remap( a, colors=None, nb_colors=None ):
    if colors is None:
        colors = np.random.random_integers( 0, 255, nb_colors ).astype('uint8')
        colors[0] = 0
    return colors[a.flatten()].reshape(a.shape)

    # tmp_a = a.copy('C').astype('int32')
    # print tmp_a.shape, colors.shape
    # _irtk.remap(tmp_a, colors.copy('C') ) #colors[a.flatten()].reshape(a.shape)
    # return tmp_a



class WindowsXYZ:
    def __init__(self, data, seg=None, overlay=None, colormap=None, opacity=0.5 ):
        data = data.astype('uint8')

        if len(data.shape) == 3:
            data = data.reshape(data.shape[0],
                                data.shape[1],
                                data.shape[2],
                                1)
            data = np.concatenate( [ data,
                                     data,
                                     data ], axis=3 )

            if seg is not None:
                seg = seg.astype('int')
                nb_colors = seg.max()+1
                if colormap is None and nb_colors <= predefined_colors.shape[0]:
                    colormap = predefined_colors
                if colormap is None:
                    rgb_overlay = np.concatenate( [ remap(seg).reshape(seg.shape[0],
                                                                       seg.shape[1],
                                                                       seg.shape[2],
                                                                       1).astype('uint8'),
                                                    remap(seg).reshape(seg.shape[0],
                                                                       seg.shape[1],
                                                                       seg.shape[2],
                                                                       1).astype('uint8'),
                                                    remap(seg).reshape(seg.shape[0],
                                                                       seg.shape[1],
                                                                       seg.shape[2],
                                                                       1).astype('uint8')
                                                    ],
                                                  axis=3
                                                  )
                else:
                    rgb_overlay = colormap[seg.flatten()].reshape(seg.shape[0],
                                                                seg.shape[1],
                                                                seg.shape[2],
                                                                3).astype('uint8')
                    
                data = (1-opacity) * data  + opacity*rgb_overlay

            if overlay is not None:
                rgb_overlay = np.concatenate(
                    [ remap(overlay,
                            colormap[:,0]).reshape(overlay.shape[0],
                                                overlay.shape[1],
                                                overlay.shape[2],
                                                1).astype('uint8'),
                      remap(overlay,
                            colormap[:,1]).reshape(overlay.shape[0],
                                                overlay.shape[1],
                                                overlay.shape[2],
                                                1).astype('uint8'),
                      remap(overlay,
                            colormap[:,2]).reshape(overlay.shape[0],
                                                overlay.shape[1],
                                                overlay.shape[2],
                                                1).astype('uint8')
                      ],
                    axis=3
                    )
                data = (1-opacity) * data  + opacity*rgb_overlay
                

        # now, data is an RGB image
        self.data = data.astype('uint8')
        self.zyx = np.array( data.shape, dtype=int ) / 2

        cv2.startWindowThread()
        cv2.namedWindow("XY")
        cv2.namedWindow("XZ")
        cv2.namedWindow("YZ")

        cv2.createTrackbar('Z', 'XY', self.zyx[0],
                           self.data.shape[0]-1, self.XY_trackbar_callback)
        cv2.createTrackbar('Y', 'XZ', self.zyx[1],
                           self.data.shape[1]-1, self.XZ_trackbar_callback)
        cv2.createTrackbar('X', 'YZ', self.zyx[2],
                           self.data.shape[2]-1, self.YZ_trackbar_callback)

        cv2.setMouseCallback("XY", self.XY_callback)
        cv2.setMouseCallback("XZ", self.XZ_callback)
        cv2.setMouseCallback("YZ", self.YZ_callback)      

        self.show()
        
    def show(self):
        cv2.imshow( "XY", self.data[self.zyx[0],:,:] )
        cv2.imshow( "XZ", self.data[:,self.zyx[1],:] )
        cv2.imshow( "YZ", self.data[:,:,self.zyx[2]] )

    def XY_trackbar_callback(self,pos):
        self.zyx[0] = pos # update z
        self.show()
    def XZ_trackbar_callback(self,pos):
        self.zyx[1] = pos # update y
        self.show()
    def YZ_trackbar_callback(self,pos):
        self.zyx[2] = pos # update x
        self.show()

    def XY_callback(self, event, x, y, flags, param):
        if event == cv2.EVENT_LBUTTONDOWN:
            self.zyx[1] = y
            self.zyx[2] = x
            self.show()

    def XZ_callback(self, event, x, y, flags, param):
        if event == cv2.EVENT_LBUTTONDOWN:
            self.zyx[0] = y
            self.zyx[2] = x
            self.show()

    def YZ_callback(self, event, x, y, flags, param):
        if event == cv2.EVENT_LBUTTONDOWN:
            self.zyx[0] = y
            self.zyx[1] = x
            self.show()            

   
def imshow( img, seg=None, overlay=None, colormap=None ):
    if overlay is not None:
        overlay = overlay.view(np.ndarray)
        if colormap is None:
            colormap = 'jet'
        if isinstance( colormap, str ):
            colormap = get_colormap( colormap )
    if seg is not None:
        seg = seg.view(np.ndarray)
    windows = WindowsXYZ( img.view(np.ndarray),
                          seg,
                          overlay, colormap )

    while True:
        ch = 0xFF & cv2.waitKey()         
        if ( ch == 27 or ch == ord('q') # ESC
             or ch == ord(' ') # SPACE
             or ch == ord('\n') ): # ENTER
            break

    cv2.destroyAllWindows()
