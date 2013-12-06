import numpy as np
import random

def HSVtoRGB( h, s, v ):
    """
    from irtk/packages/rview/src/irtkColor.cc
    """
    if s == 0:
        if h < 0:
            r = int(255.0*v)
            g = int(255.0*v);
            b = int(255.0*v);
        else:
            raise ValueError( "irtkColor::HSV: Undefined HSV color" )

    else:
        if h == 1:
            h = 0
        h = h * 6
        i = int(h)
        f = h - i
        p = v * (1 - s)
        q = v * (1 - (s * f))
        t = v * (1 - (s * (1 - f)))

        if i == 0:
            r = int(255.0*v)
            g = int(255.0*t)
            b = int(255.0*p)
        elif i == 1:
            r = int(255.0*q)
            g = int(255.0*v)
            b = int(255.0*p)
        elif i == 2:
            r = int(255.0*p)
            g = int(255.0*v)
            b = int(255.0*t)
        elif i == 3:
            r = int(255.0*p)
            g = int(255.0*q)
            b = int(255.0*v)
        elif i == 4:
            r = int(255.0*t)
            g = int(255.0*p)
            b = int(255.0*v)
        elif i == 5:
            r = int(255*v)
            g = int(255*p)
            b = int(255*q)

    return r, g, b
        
def get_colormap( name, min=0, max=255 ):
    if name in [ "rainbow", "jet" ]:
        colormap = map( lambda x: HSVtoRGB(float(max - x)/float(max - min)*2.0/3.0,
                                           1.0,
                                           1.0 ),
                        range(min,max+1) )
        colormap = np.array( colormap, dtype='uint8' )
    else:
        raise ValueError( "Unknown colormap: " + name )
    color_dict = {}
    for i in xrange(min,max+1):
        color_dict[i] = colormap[i]
    return color_dict

def colormap( color_dict, filename ):
    f = open(filename,'w')

    f.write('irtkSegmentTable: ' + str(len(color_dict)) + "\n")    
    for i, (r,g,b) in color_dict.iteritems():
        f.write(str(i) + ' ')
        f.write(str(r) + ' ')
        f.write(str(g) + ' ')
        f.write(str(b) + ' ')
        f.write(str(1) + ' ')
        f.write(str(1) + ' ')
        f.write(str(i) + "\n")

    f.close()
        
def random_colormap(N):
    # http://cloford.com/resources/colours/500col.htm
    default_colors = np.array([[255,0,0],     # green
                               [0,255,0],     # red
                               [0,245,255],   # turquoise
                               [143,188,143], # darkseagreen
                               [250,128,114], # salmon                               
                               [218,112,214], # orchid
                               [255,182,193], # lightpink
                               [237,145,33],  # carrot
                               [61,89,171],   # cobalt
                               ], dtype='uint8' )
    color_dict = {}
    for i in xrange(1,N+2):
        if i <= default_colors.shape[0]:
            color_dict[i] = default_colors[i-1]
        else:
            color_dict[i] = np.random.random_integers( 0, 255, 3 )
    return color_dict

def remap( a, color_dict=None, colors=None ):
    if colors is None:
        colors = np.zeros( (np.max(color_dict.keys())+1,3) )
        for i, c in color_dict.iteritems():
            colors[i] = c
    return colors[a.flatten()].reshape(a.shape[0], a.shape[1],3).astype('uint8')

    

