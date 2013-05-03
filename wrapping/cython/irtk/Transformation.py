from __future__ import division

import numpy as np
from math import cos, sin, pi, asin, atan2

import _irtk
import irtk

class RigidTransformation:
    def __init__( self, filename=None, matrix=None,
                  tx=0, ty=0, tz=0,
                  rx=0, ry=0, rz=0 ):
        if filename is not None:
            ( tx, ty, tz,
              rx, ry, rz ) = _irtk.read_rigid( filename )
        elif matrix is not None:
            ( tx, ty, tz,
              rx, ry, rz ) = self.__from_matrix( matrix )
            
        ( self.tx, self.ty, self.tz,
          self.rx, self.ry, self.rz ) = map( float,
                                             ( tx, ty, tz,
                                               rx, ry, rz ) )

    def get_parameters( self ):
        return ( self.tx, self.ty, self.tz,
                 self.rx, self.ry, self.rz )

    def __repr__( self ):
        txt = "tx: " + str(self.tx) + " (mm)\n"
        txt += "ty: " + str(self.ty) + " (mm)\n"
        txt += "tz: " + str(self.tz) + " (mm)\n"
        txt += "rx: " + str(self.rx) + "(degrees)\n"
        txt += "ry: " + str(self.ry) + "(degrees)\n"
        txt += "rz: " + str(self.rz)+ " (degrees)"
        return txt

    def matrix( self ):
        cosrx = cos(self.rx*(pi/180.0))
        cosry = cos(self.ry*(pi/180.0))
        cosrz = cos(self.rz*(pi/180.0))
        sinrx = sin(self.rx*(pi/180.0))
        sinry = sin(self.ry*(pi/180.0))
        sinrz = sin(self.rz*(pi/180.0))

        # Create a transformation whose transformation matrix is an identity matrix
        m = np.eye( 4, dtype='float64' )

         # Add other transformation parameters to transformation matrix
        m[0,0] = cosry*cosrz
        m[0,1] = cosry*sinrz
        m[0,2] = -sinry
        m[0,3] = self.tx

        m[1,0] = (sinrx*sinry*cosrz-cosrx*sinrz)
        m[1,1] = (sinrx*sinry*sinrz+cosrx*cosrz)
        m[1,2] = sinrx*cosry
        m[1,3] = self.ty

        m[2,0] = (cosrx*sinry*cosrz+sinrx*sinrz)
        m[2,1] = (cosrx*sinry*sinrz-sinrx*cosrz)
        m[2,2] = cosrx*cosry
        m[2,3] = self.tz
        m[3,3] = 1.0

        return m

    def __from_matrix( self, m ):
        m = m.astype('float64')
        TOL = 0.000001

        tx = m[0,3]
        ty = m[1,3]
        tz = m[2,3]

        tmp = asin(-1 * m[0,2])

        # asin returns values for tmp in range -pi/2 to +pi/2, i.e. cos(tmp) >=
        # 0 so the division by cos(tmp) in the first part of the if clause was
        # not needed.
        if abs(cos(tmp)) > TOL:
            rx = atan2(m[1,2], m[2,2])
            ry = tmp
            rz = atan2(m[0,1], m[0,0])
        else:
            # m[0,2] is close to +1 or -1
           rx = atan2(-1.0*m[0,2]*m[1,0], -1.0*m[0,2]*m[2,0])
           ry = tmp
           rz = 0
           
        # Convert to degrees.
        rx *= 180.0/pi
        ry *= 180.0/pi
        rz *= 180.0/pi

        return ( tx, ty, tz,
                 rx, ry, rz )

    def write( self, filename ):
        _irtk.write_rigid( filename,
                           self.tx, self.ty, self.tz,
                           self.rx, self.ry, self.rz )

    def apply( self, img, target_header=None,
               interpolation='linear', gaussian_parameter=1.0 ):
        if target_header is None:
            target_header = img.get_header()
        if isinstance( target_header, irtk.Image ):
            target_header = target_header.get_header()
        data = img.get_data('cython', dtype='float32')
        new_data = _irtk.transform_rigid( self.tx, self.ty, self.tz,
                                          self.rx, self.ry, self.rz,
                                          data,
                                          img.get_header(),
                                          target_header,
                                          interpolation,
                                          gaussian_parameter )
        return irtk.Image( new_data, target_header )

    def invert( self ):
        return RigidTransformation( matrix=np.linalg.inv(self.matrix()) )

    def __mul__(self, other):
        """
        Overloading multiply operator to simulate the composition of transformations:
        returns self * other
        """
        return RigidTransformation( matrix=np.dot(self.matrix(),
                                             other.matrix()) )


def registration_rigid( source, target, transformation=None ):
    if transformation is None:
        transformation = RigidTransformation()
    tx, ty, tz, rx, ry, rz = transformation.get_parameters()
    tx, ty, tz, rx, ry, rz = _irtk.registration_rigid( source.get_data('cython',
                                                                       dtype='int16'),
                                                       source.get_header(),
                                                       target.get_data('cython',
                                                                       dtype='int16'),
                                                       target.get_header(),
                                                       tx, ty, tz,
                                                       rx, ry, rz )
    return RigidTransformation( tx=tx, ty=ty, tz=tz,
                                rx=rx, ry=ry, rz=rz )
                                                       
