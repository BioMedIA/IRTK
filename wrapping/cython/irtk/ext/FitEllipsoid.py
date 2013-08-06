import numpy as np
import scipy.ndimage as nd
import SimpleITK as sitk
from scipy.stats.mstats import mquantiles
import itertools

r = [-1,0, 1]
neighbourhood = []
for a,b,c in itertools.product(r,r,r):
    if abs(a)+abs(b)+abs(c) == 3:
        neighbourhood.append((a,b,c))

def fit_ellipsoid( points, spacing=[1.0,1.0,1.0] ):
    """ Adapted from cvFitEllipse2 in
        OpenCV-2.4.3/modules/imgproc/src/shapedescr.cpp
        
        Ax**2 + By**2 + Cz**2 + Dxy + Exz + Fyz + Gx + Hy + Iz = 1
    
    """

    if len(points) < 5:
        raise ValueError("Number of points should be >= 5")

    points = points.astype('float')
    points[:,0] *= spacing[2]
    points[:,1] *= spacing[1]
    points[:,2] *= spacing[0]
    #points = np.dot(points.astype('float'), np.array(spacing[::-1]))
    
    # first fit for parameters A - E
    center = points.mean(axis=0)
    #print "First center:", center
    
    A = np.zeros( (len(points),9), dtype='float' )
    A[:,0] = - (points[:,0]-center[0]) * (points[:,0]-center[0]) # x**2
    A[:,1] = - (points[:,1]-center[1]) * (points[:,1]-center[1]) # x**2
    A[:,2] = - (points[:,2]-center[2]) * (points[:,2]-center[2]) # x**2
    A[:,3] = - (points[:,0]-center[0]) * (points[:,1]-center[1]) # xy
    A[:,4] = - (points[:,0]-center[0]) * (points[:,2]-center[2]) # xz
    A[:,5] = - (points[:,1]-center[1]) * (points[:,2]-center[2]) # yz
    A[:,6] = (points[:,0]-center[0])  # x
    A[:,7] = (points[:,1]-center[1])  # y
    A[:,8] = (points[:,2]-center[2])  # z

    B = np.ones( (len(points),1), dtype='float' ) * 10000 

    x = np.linalg.lstsq( A, B )[0].flatten()

    # M = np.array([[x[0],x[3]/2,x[4]/2],
    #               [x[3]/2,x[1],x[5]/2],
    #               [x[4]/2,x[5]/2,x[2]]],dtype='float')

    # print M,M.shape

    # u, s, v = np.linalg.svd(M)

    # return center,s,v

    # now use general-form parameters A - E to find the ellipse center:
    # differentiate general form wrt x/y to get two equations for cx and cy

    A = np.zeros( (3,3), dtype='float' )
    B = np.zeros( (3,1), dtype='float' )
    A[0,0] = 2*x[0] # 2*A
    A[0,1] = x[3] # D
    A[0,2] = x[4] # E
    B[0] = x[6] # G

    A[1,0] = x[3] # D
    A[1,1] = 2*x[1] # 2*B
    A[1,2] = x[5] # F
    B[1] = x[7] # G

    A[2,0] = x[4] # E
    A[2,1] = x[5] # F
    A[2,2] = 2*x[2] # 2*C
    B[2] = x[8] # H

    x = np.linalg.lstsq( A, B )[0].flatten()
    center += x

    #print "Second center:", center

    # re-fit for parameters A - C with those center coordinates
    A = np.zeros( (len(points),6), dtype='float' )
    A[:,0] = (points[:,0]-center[0]) * (points[:,0]-center[0]) # x**2
    A[:,1] = (points[:,1]-center[1]) * (points[:,1]-center[1]) # y**2
    A[:,2] = (points[:,2]-center[2]) * (points[:,2]-center[2]) # z**2
    A[:,3] = (points[:,0]-center[0]) * (points[:,1]-center[1]) # xy
    A[:,4] = (points[:,0]-center[0]) * (points[:,2]-center[2]) # xz
    A[:,5] = (points[:,1]-center[1]) * (points[:,2]-center[2]) # yz
    # A[:,6] = points[:,0]-center[0]  # x
    # A[:,7] = points[:,1]-center[1]  # y
    # A[:,8] = points[:,2]-center[2]  # z

    B = np.ones( (len(points),1), dtype='float' )

    x = np.linalg.lstsq( A, B )[0].flatten()

    M = np.array([[x[0],x[3]/2,x[4]/2],
                  [x[3]/2,x[1],x[5]/2],
                  [x[4]/2,x[5]/2,x[2]]],dtype='float')

    #print M,M.shape

    u, s, v = np.linalg.svd(M)

    return center,1/np.sqrt(s[::-1]),v[::-1]

def fit_ellipsoid2( points, spacing=[1.0,1.0,1.0] ):
    points = points.astype('float')
    points[:,0] *= spacing[2]
    points[:,1] *= spacing[1]
    points[:,2] *= spacing[0]

    center = points.mean(axis=0)
    points -= center
    U, S, V = np.linalg.svd(points, full_matrices=False)
    
    # http://stackoverflow.com/questions/4801259/whats-wrong-with-my-pca
    #S /= len(points) - 1

    #print "ellipse",center, S, V
    #S = 1/S#**2
    S *= 1.96
    S /= np.sqrt(len(points)-1)
    return center,S,V
    #return center, 1/S[::-1]**2, V[::-1]

def draw_ellipsoid( img, (z0,y0,x0,z1,y1,x1), center, s, v, l, spacing=[1.0,1.0,1.0] ):
    spacing = np.array(spacing[::-1],dtype='float')
    bbox = np.array([ center - np.dot(v,(s*np.array(offset,dtype='float')).T)
             for offset in neighbourhood],dtype='int')
    z0 = max(0,np.min(bbox[:,0]))
    z1 = min(img.shape[0],np.max(bbox[:,0])+1)
    y0 = max(0,np.min(bbox[:,1]))
    y1 = min(img.shape[1],np.max(bbox[:,1])+1)
    x0 = max(0,np.min(bbox[:,2]))
    x1 = min(img.shape[2],np.max(bbox[:,2])+1)
    for z in range(z0, z1):
        for y in range(y0,y1):
            for x in range(x0,x1):
                p = np.array([z,y,x],dtype='float')*spacing-center
                p = np.dot( v, p.T)
                #if s[0]*p[0]**2 + s[1]*p[1]**2 + s[2]*p[2]**2 <= 1:
                if (p[0]/s[0])**2 + (p[1]/s[1])**2 + (p[2]/s[2])**2 <= 1:
                    img[z,y,x] = l
    return img

# def get_voxels(labels,raw_spacing=[1.0,1.0,1.0],selected=None):
#     v = [[] for i in range(np.max(labels)+1)]
#     obj = nd.find_objects(labels)
#     for i,o in enumerate(obj):
#         l = i+1
#         if selected is not None and l not in selected:
#             #print l, "is not in", selected
#             continue
#         if o is None:
#             #print "o is None"
#             continue
#         for z in range(o[0].start,o[0].stop):
#             for y in range(o[1].start,o[1].stop):
#                 for x in range(o[2].start,o[2].stop):
#                     if labels[z,y,x] == l:
#                         v[l].append((z,y,x))

#     return v

def get_voxels(labels,max_label):
    count = np.bincount(labels.flatten())

    tmp = np.unravel_index(np.argsort(labels,axis=None),labels.shape)
    v = [np.zeros((count[i],3),dtype='int') for i in range(max_label+1)]
    s = 0
    for i in range(max_label+1):
        v[i][:,0] = tmp[0][s:s+count[i]]
        v[i][:,1] = tmp[1][s:s+count[i]]
        v[i][:,2] = tmp[2][s:s+count[i]]
        s += count[i]
           
    return v

def adjacent_labels(labels):
    m1 = nd.maximum_filter(labels, size=3)
    m1[m1==labels] = 0
    m2 = nd.minimum_filter(labels, size=3)
    m2[m2==labels] = 0
    m1[m2>0] = m2[m2>0]
    return m1

def boundaries(labels,keep_labels=True):
    border_pixels = (nd.convolve(labels,np.ones((3,3,3),dtype='float')/27) - labels) != 0
    if not keep_labels:
        return border_pixels
    borders = labels.copy()
    borders[border_pixels==0] = 0
    return borders

# def intensity_histogram( img, (z0,y0,x0,z1,y1,x1), center, s, v ):
#     hist = np.zeros(3,dtype='float')
#     R_max = 4*max(np.sqrt(1/s[1]),np.sqrt(1/s[2]))+1
#     R0 = int(R_max)
#     count = np.zeros(10,dtype='float')
#     for z in range(max(0,z0-R0), min(img.shape[0],z1+R0)):
#         for y in range(max(0,y0-R0),min(img.shape[1],y1+R0)):
#             for x in range(max(0,x0),min(img.shape[2],x1+R0)):
#                 p = np.array([z,y,x],dtype='float') - center
#                 p = np.dot( v, p.T)
#                 if ( p[0]**2/(1/s[0])
#                      + p[1]**2/(4**2* 1/s[1]) + p[2]**2/(4**2* 1/s[2])
#                      <= 1 ):
#                     R = np.sqrt(p[1]**2 + p[2]**2)
#                 #if R < R_max:
#                     index = int(R/R_max*3)
#                     hist[index] += img[z,y,x]
#                     count[index] += 1
#     for i in range(len(hist)):
#         if count[i] > 0:
#             hist[i] /= count[i]
#     hist /= np.linalg.norm(hist)
#     return hist

# def intensity_diff( img, (z0,y0,x0,z1,y1,x1), center, s, v ):
#     inside = 0.0
#     outside = 0.0
#     nb_in = 0
#     nb_out = 0
#     R_max = 2*max(np.sqrt(1/s[1]),np.sqrt(1/s[2]))
#     R_in = np.mean([np.sqrt(1/s[1]),np.sqrt(1/s[2])])
#     R0 = int(R_max)
#     for z in range(max(0,z0-R0), min(img.shape[0],z1+R0)):
#         for y in range(max(0,y0-R0),min(img.shape[1],y1+R0)):
#             for x in range(max(0,x0-R0),min(img.shape[2],x1+R0)):
#                 p = np.array([z,y,x],dtype='float') - center
#                 p = np.dot( v, p.T)
#                 if ( p[0]**2/(1/s[0])
#                      + p[1]**2/(2**2* 1/s[1]) + p[2]**2/(2**2* 1/s[2])
#                      <= 1 ):                
#                     if ( p[0]**2/(1/s[0])
#                          + p[1]**2/(1/s[1]) + p[2]**2/ 1/s[2])
#                          <= 1 ):
#                 R = np.sqrt(p[1]**2 + p[2]**2)
#                 if R < R_max:
#                     if R <= R_in:
#                         inside += img[z,y,x]
#                         nb_in += 1
#                     else:
#                         outside +=  img[z,y,x]
#                         nb_out += 1

#     inside /= nb_in
#     outside /= nb_out

#     return (inside, outside)

def features( img, labels, l, (z0,y0,x0,z1,y1,x1), center, s, v, points, n_slic ):
    outside = img[z0:z1,y0:y1,x0:x1].mean()
    d = z1 - z0
    h = y1 - y0
    w = x1 - x0
    zc = (z0+z1)/2
    yc = (y0+y1)/2
    xc = (x0+x1)/2
    outside2 = img[max(0,zc - d/2):min(img.shape[0],zc+d/2),
                   max(0,yc - h/2):min(img.shape[1],yc+h/2),
                   max(0,xc - w/2):min(img.shape[2],xc+w/2)].mean()
    voxels = img[points[:,0],
                points[:,1],
                points[:,2]]
    gravity_center = points.mean(axis=0)
    
    # for z in range(z0,z1):
    #     for y in range(y0,y1):
    #         for x in range(x0,x1):
    #             #p = np.array([z,y,x],dtype='float') - center
    #             # p = np.dot( v, p.T)             
    #             # if s[0]*p[0]**2 + s[1]*p[1]**2 + s[2]*p[2]**2 <= 1:
    #             #     voxels.append(img[z,y,x])
    #             if labels[z,y,x] == l:
    #                 voxels.append(img[z,y,x])

    voxels = np.array(voxels)
    inside = voxels.mean()

    a,b,c = s#(1/np.sqrt(s))[::-1] # a,b,c are the lenghts of the semi-axes
    V = 4.0/3.0*np.pi*a*b*c

    # e1 = np.sqrt(1-(b/a)**2)
    # e2 = np.sqrt(1-(c/b)**2)

    # ratio = np.sqrt(1 - ((b+c)/(a+b+c))**2) # percentage of explained variance?

    nb_points = float(len(points))
    return ( a/(a+b+c),b/a,c/a,c/b,
             (inside-outside)/(inside+0.0001),(inside-outside2)/(inside+0.0001),
             voxels.std(),
             nb_points/n_slic,
             nb_points/V,
             nb_points/(d*h*w),
             V/n_slic,
             d*h*w/V,
             np.linalg.norm(gravity_center-center)/c)

if __name__=="__main__":
    test = "output/1671_slic.nii"
    sitk_img = sitk.ReadImage( "output/1671_resampled.nii" )
    img = sitk.GetArrayFromImage( sitk_img ).astype("float")

    ## Contrast-stretch with saturation
    q = mquantiles(img.flatten(),[0.01,0.99])
    img[img<q[0]] = q[0]
    img[img>q[1]] = q[1]
    img -= img.min()
    img /= img.max()
    
    sitk_img = sitk.ReadImage( test )
    data = sitk.GetArrayFromImage( sitk_img )

# print nd.find_objects(data)

# data = data[nd.find_objects(data)[4]]
# print data

    labels = np.unique( data )
    voxels = get_voxels(data)
    print "voxels done"
    objects = nd.find_objects(data)

    cyl = np.zeros(data.shape,dtype='uint8')
    
    for l in labels:
        if l == 0:
            continue
        print l
        points = np.array(voxels[l])
        center,s,v = fit_ellipsoid(points)
        obj = objects[l-1]
        cyl = draw_ellipsoid(cyl, (obj[0].start,obj[1].start,obj[2].start,
                                   obj[0].stop,obj[1].stop,obj[2].stop),
                             center,s,v,l)
        print "HIST:"
        print features( img,
                                   (obj[0].start,obj[1].start,obj[2].start,
                                    obj[0].stop,obj[1].stop,obj[2].stop),
                                   center, s, v, points,1000 )

    sitk_seg = sitk.GetImageFromArray( cyl )
    sitk_seg.SetSpacing( sitk_img.GetSpacing() )
    sitk_seg.SetOrigin( sitk_img.GetOrigin() )
    sitk_seg.SetDirection( sitk_img.GetDirection() )
    sitk.WriteImage( sitk_seg, "ellipsoids.nii" )
