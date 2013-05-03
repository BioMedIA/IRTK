import _slic

def slic( img, size=300, compactness=10 ):
    """
    SLIC supervoxels
    """
    return _slic.slic3D( img.get_data('uint8'),
                         size,
                         compactness,
                         spacing=img.header['pixelSize'][:3] )
