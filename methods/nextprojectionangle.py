#!/usr/bin/env python
'''
Created on Sep 14, 2018

@author: melnikov
'''

import numpy
import json
import base64
from scipy import ndimage

def FWHM(array_1D):
    if numpy.all(array_1D==0):
        return 0
    else:
        C = ndimage.measurements.center_of_mass(array_1D)
        H = 2.355*numpy.sqrt(numpy.sum(array_1D*(numpy.arange(0, len(array_1D))-C)**2)/numpy.sum(array_1D))
        return H

def determine(json1, json2):
    
    f = open(json1, 'r')
    jdata1 = json.load(f)
    f.close()
    f = open(json2, 'r')
    jdata2 = json.load(f)
    f.close()
    
    Dtable1 = numpy.fromstring(base64.b64decode(jdata1['MeshBest']['Dtable']))
    Dtable1 = numpy.reshape(Dtable1, (jdata1['grid_info']['steps_y'], jdata1['grid_info']['steps_x']))

    Dtable2 = numpy.fromstring(base64.b64decode(jdata2['MeshBest']['Dtable']))
    Dtable2 = numpy.reshape(Dtable2, (jdata2['grid_info']['steps_y'], jdata2['grid_info']['steps_x']))

    Hmax1 = max([FWHM(column) for column in Dtable1.T])
    Hmax2 = max([FWHM(column) for column in Dtable2.T])
    
    angle = numpy.arctan(Hmax2/Hmax1)
    
    return (angle, -angle)
