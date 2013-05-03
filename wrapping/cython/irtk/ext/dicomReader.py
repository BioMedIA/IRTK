import dicom
import numpy as np

def load_dcm(filename):
    dcm = dicom.ReadFile(filename)

    # http://www.creatis.insa-lyon.fr/pipermail/dcmlib/2005-September/002141.html
    # http://www.dabsoft.ch/dicom/3/C.7.6.2.1.1/
    # http://xmedcon.sourceforge.net/Docs/OrientationDicomStandard
    img0 = dcm[(0x5200, 0x9230)][0] # Per-frame Functional Groups Sequence 
    spacingXY = img0[(0x0028, 0x9110)][0][(0x0028, 0x0030)].value # Pixel Spacing

    pos0 = img0[(0x0020, 0x9113)][0][(0x0020, 0x0032)].value # Image Position (Patient)
    orientation = img0[(0x0020, 0x9116)][0][(0x0020, 0x0037)].value # Image Orientation (Patient)    

    i = 1
    spacingZ = 0
    while spacingZ == 0:
        img1 = dcm[(0x5200, 0x9230)][i] # Per-frame Functional Groups Sequence
        pos1 = img1[(0x0020, 0x9113)][0][(0x0020, 0x0032)].value # Image Position (Patient)

        spacingZ = np.linalg.norm(np.array(pos1) - np.array(pos0) )
        i +=1

    dcm_spacing = [spacingXY[0],spacingXY[1],spacingZ]
    #dcm_spacing = [spacingZ,spacingXY[1],spacingXY[0]]
    
    print 'DICOM  - spacing', dcm_spacing
    
    X = np.array([orientation[0],orientation[1],orientation[2]])
    Y = np.array([orientation[3],orientation[4],orientation[5]])
    Z = np.cross(X,Y)

    X *= dcm_spacing[0]
    Y *= dcm_spacing[1]
    Z *= dcm_spacing[2]
    scaling = np.diag(dcm_spacing + [1] )
    dcm_affine = np.array([[X[0],Y[0],Z[0],pos0[0]],
                           [X[1],Y[1],Z[1],pos0[1]],
                           [X[2],Y[2],Z[2],pos0[2]],
                           [0,0,0,1]])

    #dcm_affine = np.dot( dcm_affine, scaling)
    

    data = dcm.pixel_array.astype('float')            

    return data, dcm_affine,dcm_spacing

def dcm2volume(filename):
    dcm = dicom.ReadFile(filename)

    data = dcm.pixel_array.astype('float')#[::-1,:,:]
    vol = []
    
    # http://www.creatis.insa-lyon.fr/pipermail/dcmlib/2005-September/002141.html
    # http://www.dabsoft.ch/dicom/3/C.7.6.2.1.1/
    # http://xmedcon.sourceforge.net/Docs/OrientationDicomStandard

    frames = dcm[(0x5200, 0x9230)].value # Per-frame Functional Groups Sequence
    cursor = 0
    length = 1
    #print len(frames.value)
    for i in range(len(frames)-1):
        current = frames[i]

        if length == 1:
            pos0 = current[(0x0020, 0x9113)][0][(0x0020, 0x0032)].value # Image Position (Patient)
            spacingXY = current[(0x0028, 0x9110)][0][(0x0028, 0x0030)].value # Pixel Spacing
            pos = current[(0x0020, 0x9113)][0][(0x0020, 0x0032)].value # Image Position (Patient)
            orientation = current[(0x0020, 0x9116)][0][(0x0020, 0x0037)].value # Image Orientation (Patient)    
        
            next = frames[i+1]
            next_pos = next[(0x0020, 0x9113)][0][(0x0020, 0x0032)].value # Image Position (Patient)
            next_orientation = next[(0x0020, 0x9116)][0][(0x0020, 0x0037)].value # Image Orientation (Patient)  

            spacingZ = np.linalg.norm(np.array(next_pos) - np.array(pos) )

            if spacingZ == 0:
                raise ValueError('null spacing: not implemented')

        orientation = current[(0x0020, 0x9116)][0][(0x0020, 0x0037)].value # Image Orientation (Patient)    
        next = frames[i+1]
        next_orientation = next[(0x0020, 0x9116)][0][(0x0020, 0x0037)].value # Image Orientation (Patient)
        
        if orientation == next_orientation and i < len(frames)-2:
            length += 1
        else:
            if i == len(frames)-2:
                length += 1
            dcm_spacing = [spacingXY[0],spacingXY[1],spacingZ]    
            X = np.array([orientation[0],orientation[1],orientation[2]])
            Y = np.array([orientation[3],orientation[4],orientation[5]])
            Z = np.cross(X,Y)

            X *= dcm_spacing[0]
            Y *= dcm_spacing[1]
            Z *= dcm_spacing[2]
            scaling = np.diag(dcm_spacing + [1] )
            dcm_affine = np.array([[X[0],Y[0],Z[0],pos0[0]],
                                   [X[1],Y[1],Z[1],pos0[1]],
                                   [X[2],Y[2],Z[2],pos0[2]],
                                   [0,0,0,1]])
            
            vol.append((data[:,:,cursor:cursor+length,],dcm_affine))
            cursor += length
            length = 1
          
    return vol

def load_movie(filename):
    dcm = dicom.ReadFile(filename)

    # http://www.creatis.insa-lyon.fr/pipermail/dcmlib/2005-September/002141.html
    # http://www.dabsoft.ch/dicom/3/C.7.6.2.1.1/
    # http://xmedcon.sourceforge.net/Docs/OrientationDicomStandard
    frames = dcm[(0x5200, 0x9230)].value # Per-frame Functional Groups Sequence
    
    img0 = frames[0] # Per-frame Functional Groups Sequence 
    spacingXY = img0[(0x0028, 0x9110)][0][(0x0028, 0x0030)].value # Pixel Spacing

    pos0 = img0[(0x0020, 0x9113)][0][(0x0020, 0x0032)].value # Image Position (Patient)
    orientation = img0[(0x0020, 0x9116)][0][(0x0020, 0x0037)].value # Image Orientation (Patient)    

    for i in range(1,len(frames)):
        img = frames[i]
        pos_i = img[(0x0020, 0x9113)][0][(0x0020, 0x0032)].value # Image Position (Patient)
        orientation_i = img[(0x0020, 0x9116)][0][(0x0020, 0x0037)].value # Image Orientation (Patient) 
        if pos_i != pos0 or orientation_i != orientation:
            raise ValueError('this is not a movie')

    dcm_spacing = [spacingXY[0],spacingXY[1],0]
    
    print 'DICOM  - spacing', dcm_spacing
    
    X = np.array([orientation[0],orientation[1],orientation[2]])
    Y = np.array([orientation[3],orientation[4],orientation[5]])

    X *= dcm_spacing[0]
    Y *= dcm_spacing[1]
    scaling = np.diag(dcm_spacing + [1] )
    dcm_affine = np.array([[X[0],Y[0],0,pos0[0]],
                           [X[1],Y[1],0,pos0[1]],
                           [X[2],Y[2],0,pos0[2]],
                           [0,0,0,1]])

    data = dcm.pixel_array.transpose()#.astype('int')            

    return data, dcm_affine,dcm_spacing

