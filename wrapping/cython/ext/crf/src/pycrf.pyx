import numpy as np
cimport numpy as np

np.import_array()

ctypedef int LabelID
ctypedef int SiteID
ctypedef float EnergyTermType
ctypedef float EnergyType

# This software is for minimizing energy functions of the form:
# E(l) = sum_p D(p,l_p)  + sum_{p,q} Vpq(l_p,l_q)

cdef extern from "GCoptimization.h":
    cdef cppclass GCoptimizationGeneralGraph:
        GCoptimizationGeneralGraph(SiteID, LabelID) # num_sites, num_labels
        void setLabel(SiteID, LabelID)
        void setDataCost(EnergyTermType*)
        void setSmoothCost(EnergyTermType*)
        void setNeighbors(SiteID, SiteID, EnergyTermType)
        EnergyType compute_energy()
        void expansion()
        void swap()
        LabelID whatLabel(SiteID)

cdef class CRF:
    cdef GCoptimizationGeneralGraph *thisptr      # hold a C++ instance which we're wrapping
    cdef SiteID num_pixels
    cdef LabelID num_labels
    def __cinit__(self, SiteID num_pixels, LabelID num_labels):
        self.thisptr = new GCoptimizationGeneralGraph(num_pixels, num_labels)
        self.num_pixels = num_pixels
        self.num_labels = num_labels
    def __dealloc__(self):
        del self.thisptr
    def setLabel(self, SiteID i, LabelID cl):
        self.thisptr.setLabel(i,cl)
    def setLabels(self, np.ndarray[np.int32_t, ndim=1] labels):
        for i in xrange(labels.shape[0]):
            self.thisptr.setLabel(i,labels[i])  
    def setDataCost(self, np.ndarray[np.float32_t, ndim=1] datacost):
        # datacost is an array s.t. the data cost for pixel p and  label l is
        # stored at datacost[pixel*num_labels+l]
        if not (<object>datacost).flags["C_CONTIGUOUS"]:
            datacost = datacost.copy('C')
        self.thisptr.setDataCost(<EnergyTermType*>datacost.data)
    def setSmoothCost(self ,np.ndarray[np.float32_t, ndim=1] labelcost):
        # labelcost is an array of smoothness costs, such that
        # V_pq(label1,label2)  is stored at labelcost[label1+num_labels*label2]
        if not (<object>labelcost).flags["C_CONTIGUOUS"]:
            labelcost = labelcost.copy('C')        
        self.thisptr.setSmoothCost(<EnergyTermType*>labelcost.data)
    def setNeighbors(self, SiteID site1, SiteID site2, EnergyTermType weight=1):
        self.thisptr.setNeighbors(site1,site2,weight)
    def compute_energy(self):
        return self.thisptr.compute_energy()
    def expansion(self):
        self.thisptr.expansion()
    def swap(self):
        self.thisptr.swap()        
    def whatLabel(self, SiteID i):
        return self.thisptr.whatLabel(i)
    def getLabels(self):
        cdef np.npy_intp labels_shape[1]
        labels_shape[0] = self.num_pixels
        cdef np.ndarray[np.int32_t, ndim=1] labels = np.PyArray_SimpleNew(1,
                                                                          labels_shape,
                                                                          np.NPY_INT32)
        cdef int * labels_ptr = <int*>labels.data
        for i in xrange(self.num_pixels):
            labels_ptr[i] = self.thisptr.whatLabel(i)
        return labels
