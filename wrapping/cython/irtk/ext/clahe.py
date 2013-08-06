import numpy as np
import _clahe

__all__ = ["CLAHE"]

def CLAHE( img, radius, clip, nbins=256, rescale=True ):
    """
import irtk.ext.clahe as clahe
import cv2
img = cv2.imread("/home/kevin/Imperial/PhD/ITK/wrapping/shu_small.png",0)
cv2.imwrite( "test.png", clahe.CLAHE(img, 29,3.0))
    """
    if rescale:
        # rescale image
        img = img.astype('float')
        img -= img.min()
        img /= img.max()
        img *= 255
        
    img = img.astype('uint8').copy()

    # pad image if necessary
    padding0 = radius - (img.shape[0] % radius)
    padding1 = radius - (img.shape[1] % radius)
    if padding0 == radius:
        padding0 = 0
    if padding1 == radius:
        padding1 = 0
        
    left0 = padding0 / 2
    right0 = padding0 - left0
    left1 = padding1 / 2
    right1 = padding1 - left1
    if right0 == 0:
        right0 = -img.shape[0]
    if right1 == 0:
        right1 = -img.shape[1]

    new_img = np.zeros( ( img.shape[0] + padding0,
                          img.shape[1] + padding1 ),
                        dtype='uint8' )

    new_img[left0:-right0,
            left1:-right1] = img

    print (radius, img.shape, new_img.shape,
           new_img.shape[1]/radius,
           new_img.shape[0]/radius)
    new_img = _clahe.CLAHE( new_img,
                            new_img.shape[1]/radius,
                            new_img.shape[0]/radius,
                            nbins,
                            clip )
                            
    if new_img is None:
        return None
    else:
        return new_img[left0:-right0,
                       left1:-right1].copy()
    
