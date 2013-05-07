import numpy as np
import random
import matplotlib.pyplot as plt

def list_colormaps():
    maps = [m for m in plt.cm.datad if not m.endswith("_r")]
    maps.sort()
    for m in maps:
        print m
        
def get_colormap( name, min=0, max=255 ):
    colormap = plt.cm.get_cmap(name, max+1-min)(range(min,max+1))*255
    colormap = colormap[:,:3].astype('uint8')
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

    

