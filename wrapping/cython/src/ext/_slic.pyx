import numpy as np
cimport numpy as np

from libcpp.vector cimport vector

np.import_array()

from libc.stdlib cimport malloc, free

ctypedef unsigned char pixeltype

cdef extern from "SLIC.h":
    cdef cppclass SLIC:
        SLIC()
        # img, width, heights, returned labels, num_labels, superpixel size, compactness
        #void DoSuperpixelSegmentation_ForGivenSuperpixelSize(unsigned int*, int, int, int *, int, int, double)
        void DoSupervoxelSegmentation(pixeltype**, int, int, int,
                                      double, double, double, int**, int, int, double)
#(unsigned int**, int, int, int, int **, int, int, double)
        # img, width, heights, returned labels, num_labels, superpixel number, compactness
        #void DoSuperpixelSegmentation_ForGivenNumberOfSuperpixels(unsigned int*, int, int, int *, int, int, double)
        # img, labeling, width, heights, color for contours
        #void DrawContoursAroundSegments(unsigned int** ,int*, int, int, unsigned int)


# def slic_s(np.ndarray[np.uint8_t, ndim=3] img, superpixel_size=300, compactness=10):
#     """SLIC Superpixels for fixed superpixel size.

#     Parameters
#     ----------
#     img : numpy array, dtype=uint8
#         Original image, ARGB (or AXXX) format, A channel is ignored.
#         Needs to be C-Contiguous.
#     superpixel_size: int, default=300
#         Desired size for superpixel
#     compactness: douple, default=10
#         Degree of compactness of superpixels.

#     Returns
#     -------
#     labels : numpy array

#     """

#     if (img.shape[2] != 3):
#         raise ValueError("Image needs to have 3 channels.")
#     if np.isfortran(img):
#         raise ValueError("The input image is not C-contiguous")
#     cdef np.ndarray[np.uint8_t, ndim=3] img_ = np.empty((img.shape[0], img.shape[1], 4), dtype=np.uint8)
#     img_[:, :, 1:] = img
#     cdef int h = img.shape[0]
#     cdef int w = img.shape[1]
#     cdef int * labels
#     cdef int n_labels
#     cdef SLIC* slic = new SLIC()
#     slic.DoSuperpixelSegmentation_ForGivenSuperpixelSize(<unsigned int *>img_.data, w, h,
#             labels, n_labels, superpixel_size, compactness)
#     cdef np.npy_intp shape[2]
#     shape[0] = h
#     shape[1] = w
#     label_array = np.PyArray_SimpleNewFromData(2, shape, np.NPY_INT32, <void*> labels)
#     return label_array


# def slic_n(np.ndarray[np.uint8_t, ndim=3] img, n_superpixels=500, compactness=10):
#     """SLIC Superpixels for fixed number of superpixels.

#     Parameters
#     ----------
#     img : numpy array, dtype=uint8
#         Original image RGB.
#         Needs to be C-Contiguous
#     n_superpixels: int, default=500
#         Desired number of superpixels.
#     compactness: douple, default=10
#         Degree of compactness of superpixels.

#     Returns
#     -------
#     labels : numpy array

#     """
#     if (img.shape[2] != 3):
#         raise ValueError("Image needs to have 3 channels.")
#     if np.isfortran(img):
#         raise ValueError("The input image is not C-contiguous")
#     cdef np.ndarray[np.uint8_t, ndim=3] img_ = np.empty((img.shape[0], img.shape[1], 4), dtype=np.uint8)
#     img_[:, :, :-1] = img
#     cdef int h = img.shape[0]
#     cdef int w = img.shape[1]
#     cdef int * labels
#     cdef int n_labels
#     cdef SLIC* slic = new SLIC()
#     slic.DoSuperpixelSegmentation_ForGivenNumberOfSuperpixels(<unsigned int *>img_.data, w, h,
#             labels, n_labels, n_superpixels, compactness)
#     cdef np.npy_intp shape[2]
#     shape[0] = h
#     shape[1] = w
#     label_array = np.PyArray_SimpleNewFromData(2, shape, np.NPY_INT32, <void*> labels)
#     return label_array


# def contours(np.ndarray[np.uint8_t, ndim=3] img, np.ndarray[np.int32_t, ndim=2] labels, color=10):
#     """Draw contours of superpixels into original image.
#     Destoys original!

#     Parameters
#     ----------
#     img : numpy array, dtype=uint8
#         Original image.
#         Needs to be uint, RGB.
#         Needs to be C-Contiguous
#     lables: numpy array, dtype=int
#         Same width and height as image, 
#     color: int
#         color for boundaries.
#     """
#     cdef SLIC* slic = new SLIC()
#     assert(img.shape[2] == 3)
#     cdef int h = img.shape[0]
#     cdef int w = img.shape[1]
#     cdef int n_labels
#     cdef np.ndarray[np.uint8_t, ndim=3] img_ = np.empty((img.shape[0], img.shape[1], 4), dtype=np.uint8)
#     img_[:, :, :-1] = img
#     slic.DrawContoursAroundSegments(<unsigned int **>&img_.data, <int*>labels.data, w, h,
#             color)
#     return img_


def slic3D(np.ndarray[np.uint8_t, ndim=3] img, supervoxel_size=300,
           compactness=10, spacing=[1.0,1.0,1.0]):
    """SLIC Superpixels for fixed superpixel size.

    Parameters
    ----------
    img : numpy array, dtype=uint8
        Original image, ARGB (or AXXX) format, A channel is ignored.
        Needs to be C-Contiguous.
    superpixel_size: int, default=300
        Desired size for superpixel
    compactness: douple, default=10
        Degree of compactness of superpixels.

    Returns
    -------
    labels : numpy array

    """

    # if (img.shape[2] != 3):
    #     raise ValueError("Image needs to have 3 channels.")
    if np.isfortran(img):
        raise ValueError("The input image is not C-contiguous")
    if not (<object>img).flags["C_CONTIGUOUS"]:
        img = img.copy('C')
    # cdef np.ndarray[np.uint8_t, ndim=4] img_ = np.empty((img.shape[0],
    #                                                      img.shape[1],
    #                                                      img.shape[2], 4), dtype=np.uint8)
    # img_[:, :, :, :-1] = img
    cdef int d = img.shape[0]
    cdef int h = img.shape[1]
    cdef int w = img.shape[2]
    cdef int ** labels
    cdef int n_labels
    cdef SLIC* slic = new SLIC()

    cdef pixeltype **data
    data = <pixeltype **>malloc(d*sizeof(pixeltype *))
    #cdef unsigned int* data = new unsigned int*[d]
    cdef int i
    cdef int j
    cdef int k
    #cdef np.ndarray[np.uint8_t, ndim=3] img_ = np.empty((h, w, 4), dtype=np.uint8)
    for i in range(d):
        data[i] = <pixeltype *>malloc(h*w*sizeof(pixeltype))
        for j in range(h):
            for k in range(w):
                data[i][j*w+k] = <pixeltype>img[i,j,k]
                    
        # img_[:, :, :-1] = img[d,:,:,:]
        # data[i] = &(<unsigned int *>img_.data)
        # for j in range(h):
        #     for k in range(w):
        #         for l in range(4):
        #             [j*w*4+k*4+l] = <unsigned int>img_[i,j,k,l]
        #data[i] = &(<unsigned int *>img_.data)[i * cols]
        #data[i] = &(<unsigned int *>img_.data)[i*h*w*4]
        #data[i] = <unsigned int *>&img_[i,0,0,0]
    print "OK"
    cdef double spacingx = spacing[0]
    cdef double spacingy = spacing[1]
    cdef double spacingz = spacing[2]
    slic.DoSupervoxelSegmentation(data, w, h, d,
                                  spacingx, spacingy, spacingz,
            labels, n_labels, supervoxel_size, compactness)

    #free(data)

    print "done"
    
    cdef int * plabels
    plabels = <int *>malloc(d*h*w*sizeof(int))
    
    for i in range(d):
        for j in range(h*w):
            plabels[i*h*w+j] = labels[i][j]
            
    cdef np.npy_intp shape[3]
    shape[0] = d
    shape[1] = h
    shape[2] = w
    label_array = np.PyArray_SimpleNewFromData(3, shape, np.NPY_INT32, <void*> plabels)
    #free(plabels)
    return label_array

# def get_voxels(np.ndarray[np.int32_t, ndim=3] labels, max_label):
#     cdef int d = labels.shape[0]
#     cdef int h = labels.shape[1]
#     cdef int w = labels.shape[2]
#     cdef vector[vector[vector[int]]] v
#     v.resize(max_label+1)
#     cdef int i
#     cdef int j
#     cdef int k    
#     for i in range(d):
#         for j in range(h):
#             for k in range(w):
#                 if labels[i,j,k] > 0:
#                     v[labels[i,j,k]].push_back((i,j,k))    
#     return v

# def get_voxels(np.ndarray[np.int32_t, ndim=3] labels, max_label):
#     cdef int i
#     cdef int j
#     cdef int k      
#     cdef int d = labels.shape[0]
#     cdef int h = labels.shape[1]
#     cdef int w = labels.shape[2]    
#     cdef vector[int] count
#     count.resize(max_label+1,0)
#     for i in range(labels.shape[0]):
#         for j in range(labels.shape[1]):
#             for k in range(labels.shape[2]):
#                 if labels[i,j,k] > 0:
#                     count[labels[i,j,k]] += 1

#     cdef vector[vector[vector[int]]] v
#     v.resize(max_label+1)
#     for i in range(max_label+1):
#         v[i].resize(count[i])
#         for j in range(count[i]):
#             v[i][j].resize(3,0)
#         count[i] = 0
        
#     for i in range(labels.shape[0]):
#         for j in range(labels.shape[1]):
#             for k in range(labels.shape[2]):
#                 if labels[i,j,k] > 0:
#                     v[labels[i,j,k]][count[labels[i,j,k]]] = (i,j,k)
#                     count[labels[i,j,k]] += 1
                                     
#     return v
