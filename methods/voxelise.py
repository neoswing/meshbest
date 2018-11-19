#!/usr/bin/env python
'''
Created on Jun 13, 2018

@author: melnikov
'''
import numpy
from scipy import ndimage


def stretch(where3D):
    stretched_where = (numpy.concatenate((where3D[0]*10, where3D[0]*10+1, where3D[0]*10+2, where3D[0]*10+3, \
                                          where3D[0]*10+4, where3D[0]*10+5, where3D[0]*10+6, where3D[0]*10+7, \
                                          where3D[0]*10+8, where3D[0]*10+9)), \
                       numpy.concatenate((where3D[1]*10, where3D[1]*10+1, where3D[1]*10+2, where3D[1]*10+3, \
                                          where3D[1]*10+4, where3D[1]*10+5, where3D[1]*10+6, where3D[1]*10+7, \
                                          where3D[1]*10+8, where3D[1]*10+9)), \
                       numpy.concatenate((where3D[2]*10, where3D[2]*10+1, where3D[2]*10+2, where3D[2]*10+3, \
                                          where3D[2]*10+4, where3D[2]*10+5, where3D[2]*10+6, where3D[2]*10+7, \
                                          where3D[2]*10+8, where3D[2]*10+9)))
    
    
    return stretched_where


def voxelCoordinates(XY_0, XZ_0, XY_1, XZ_1, fi_list):
#    checking equality of X steps numbers
    if numpy.shape(XY_0)[1]==numpy.shape(XZ_0)[1] and numpy.shape(XY_1)[1]==numpy.shape(XZ_1)[1] and \
    numpy.shape(XY_0)[1]==numpy.shape(XY_1)[1]:
        offset1 = ndimage.measurements.center_of_mass(XY_0)
        offset2 = ndimage.measurements.center_of_mass(XZ_0)
        offset3 = ndimage.measurements.center_of_mass(XY_1)
        offset4 = ndimage.measurements.center_of_mass(XZ_1)

        fi3 = fi_list[2] - fi_list[0]
        fi4 = fi_list[3] - fi_list[0]

        L, W, H = numpy.shape(XY_0)[1], numpy.shape(XY_0)[0], numpy.shape(XZ_0)[0]
        
        voxel = numpy.ones(((L, W, H)))
        XY_0_indicator = XY_0>0.3
        XZ_0_indicator = XZ_0>0.3
        XY_1_indicator = XY_1>0.3
        XZ_1_indicator = XZ_1>0.3
        
        where = numpy.where(XY_0_indicator!=1)
        voxel[where[1], where[0], :] = 0
        where = numpy.where(XZ_0_indicator!=1)
        voxel[where[1], :, where[0]] = 0
        
        where = stretch(numpy.where(voxel==1))
        
        X = where[0]
        
        Y1 = numpy.cos(3.14*fi3/180.0)*(where[1]-offset1[0]*10) + \
             numpy.sin(3.14*fi3/180.0)*(where[2]-offset2[0]*10) + offset3[0]*10
        
        Z1 = -numpy.cos(3.14*fi4/180.0)*(where[1]-offset1[0]*10) + \
             numpy.sin(3.14*fi4/180.0)*(where[2]-offset2[0]*10) + offset4[0]*10
        
        Check = numpy.multiply((XY_1_indicator[(Y1/10.0-0.5).astype('int'), (X/10.0-0.5).astype('int')]==1), \
                               (XZ_1_indicator[(Z1/10.0-0.5).astype('int'), (X/10.0-0.5).astype('int')]==1))
        
        coordinates = (where[0][Check]/10.0, where[1][Check]/10.0, where[2][Check]/10.0)
        
        Darray = numpy.min([XY_0[coordinates[1].astype('int'), coordinates[0].astype('int')],\
        XZ_0[coordinates[2].astype('int'), coordinates[0].astype('int')],\
        XY_1[(numpy.cos(3.14*fi3/180.0)*(coordinates[1]-offset1[0]) + \
             numpy.sin(3.14*fi3/180.0)*(coordinates[2]-offset2[0]) + offset3[0]).astype('int'), coordinates[0].astype('int')],\
        XZ_1[(-numpy.cos(3.14*fi4/180.0)*(coordinates[1]-offset1[0]) + \
             numpy.sin(3.14*fi4/180.0)*(coordinates[2]-offset2[0]) + offset4[0]).astype('int'), coordinates[0].astype('int')]], axis=0)
        
        CenterOfMass = numpy.array([numpy.sum(numpy.multiply(Darray, coordinates[0])),\
                                    numpy.sum(numpy.multiply(Darray, coordinates[1])),\
                                    numpy.sum(numpy.multiply(Darray, coordinates[2]))])/numpy.sum(Darray)
        
        return coordinates
    
    else:
        print 'Error: mesh scans have different lengths'
        return numpy.nan

def CenterOfMass(XY_0, XZ_0, XY_1, XZ_1, fi_list):
    
    coordinates = voxelCoordinates(XY_0, XZ_0, XY_1, XZ_1, fi_list)
    if coordinates:
        offset1 = ndimage.measurements.center_of_mass(XY_0)
        offset2 = ndimage.measurements.center_of_mass(XZ_0)
        offset3 = ndimage.measurements.center_of_mass(XY_1)
        offset4 = ndimage.measurements.center_of_mass(XZ_1)
        
        fi3 = fi_list[2] - fi_list[0]
        fi4 = fi_list[3] - fi_list[0]
        
        Darray = numpy.min([XY_0[coordinates[1].astype('int'), coordinates[0].astype('int')],\
        XZ_0[coordinates[2].astype('int'), coordinates[0].astype('int')],\
        XY_1[(numpy.cos(3.14*fi3/180.0)*(coordinates[1]-offset1[0]) + \
             numpy.sin(3.14*fi3/180.0)*(coordinates[2]-offset2[0]) + offset3[0]).astype('int'), coordinates[0].astype('int')],\
        XZ_1[(-numpy.cos(3.14*fi4/180.0)*(coordinates[1]-offset1[0]) + \
             numpy.sin(3.14*fi4/180.0)*(coordinates[2]-offset2[0]) + offset4[0]).astype('int'), coordinates[0].astype('int')]], axis=0)
        
        CenterOfMass = numpy.array([numpy.sum(numpy.multiply(Darray, coordinates[0])),\
                                    numpy.sum(numpy.multiply(Darray, coordinates[1])),\
                                    numpy.sum(numpy.multiply(Darray, coordinates[2]))])/numpy.sum(Darray)
        return CenterOfMass
    else:
        return None
