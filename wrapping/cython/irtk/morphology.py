import numpy as np

__all__ = [ "disk" ]

def disk( r ):
    """
    2D disk for morphological operations.
    """
    mask = np.zeros( (2*r+1,2*r+1), dtype='uint8' )
    for x in range(-r,r+1):
        for y in range(-r,r+1):
            if x**2 + y**2 <= r**2:
                mask[y+r,x+r] = 1
    return mask
