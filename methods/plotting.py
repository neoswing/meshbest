#!/usr/bin/env python
'''
Created on May 15, 2018

@author: melnikov
'''

import numpy
import base64
from matplotlib import colors, pyplot as plt
from matplotlib.patches import Circle, Ellipse as El
from scipy import ndimage, signal
import random

try:
    from workflow_lib import workflow_logging
    logger = workflow_logging.getLogger()
except:
    import logging
    logger = logging.getLogger("MeshBest")
    

def ConstructColorlist(array):
    basecolors = ['#FF0101', '#F5A26F', '#668DE5', '#E224DE', '#04FEFD', '#00CA02', '#FEFE00', '#0004AF', '#B5FF06']

    N = int(numpy.max(array))
    AdjacentArray = numpy.identity(N)

    t = numpy.ones((3, 3), dtype='int32')
    for j in range(1, N + 1):
        cut = ndimage.measurements.label(array == j)
        c = signal.convolve2d(cut[0], t, mode='same')
        adjacentvalues = numpy.unique(array[numpy.where(c != 0)]).astype('int')

        for i in adjacentvalues:
            if i == -1 or i == -2:
                pass
            else:
                AdjacentArray[j - 1, i - 1] = 1
    t = [1, 2, 3, 4, 5, 6, 7, 8, 9]
    Adjacent_NULL = numpy.tril(AdjacentArray, k= -1)
    AdjacentArray = Adjacent_NULL

    ColorVector = numpy.ones(N)
    for i in xrange(N):
        BannedColors = numpy.unique(AdjacentArray[i, :])[1:]
        for item in BannedColors:
            t.remove(item)
        ColorVector[i] = t[0]
        AdjacentArray = ColorVector * Adjacent_NULL
        t = [1, 2, 3, 4, 5, 6, 7, 8, 9]
    ColorVector = ColorVector.astype(int)

#    random.shuffle(basecolors)
    Colors = ColorVector.astype('a7')
    for i in xrange(N):
        Colors[i] = basecolors[ColorVector[i] - 1]

    clrs = [(0.3, 0.3, 0.3), 'black', 'black']
    clrs.extend(Colors.tolist())

    return clrs





def MainPlot(jsondata, ax, addPositions=True):
    
    try:
        row, col = jsondata['grid_info']['steps_y'], jsondata['grid_info']['steps_x']
        difminpar = jsondata['MeshBest']['difminpar']
    except KeyError:
        logger.error('Experiment parameters are not communicated in the JSON')
        return None
    try:
        Dtable = numpy.fromstring(base64.b64decode(jsondata['MeshBest']['Dtable']))
        Dtable = numpy.reshape(Dtable, (row, col))
        Ztable = numpy.fromstring(base64.b64decode(jsondata['MeshBest']['Ztable']))
        Ztable = numpy.reshape(Ztable, (row, col))
    except KeyError:
        logger.error('Plotting: No data to work with in the JSON')
        return None
    
    if numpy.all(Ztable==-1):
        logger.error('Plotting: Crystal map is empty')
        return None
    elif numpy.all(Ztable<0):
        logger.info('Plotting: Only crystal overlaps were found')
        return None
    
    clrs = ConstructColorlist(Ztable)
#    logger.debug('Checkpoint for colormap function:', (time.time()-start_time))
    Ncolors = len(clrs)
    cmap = colors.ListedColormap(clrs)
    cmap2 = colors.LinearSegmentedColormap.from_list('my_cmap2', [(0, 0, 0), (0, 0, 0.001)], 256)
    cmap2._init()
    alphas = numpy.linspace(-0.01, -1, cmap2.N + 3)
    cmap2._lut[:, -1] = alphas

    bounds = numpy.linspace(-2, Ncolors - 2, Ncolors + 1)
    norm = colors.BoundaryNorm(bounds, ncolors=Ncolors)
    
    CrystImage = plt.imshow(Ztable, cmap=cmap, norm=norm, interpolation='nearest', origin='upper', \
                            extent=[0.5, (col + 0.5), (row + 0.5), 0.5])
    OpacityImage = plt.imshow(Dtable, cmap=cmap2, interpolation='nearest', origin='upper', vmin=difminpar, \
                              vmax=numpy.percentile(Dtable[Dtable>difminpar], 99), extent=[0.5, (col + 0.5), (row + 0.5), 0.5])
    
    
    opacity = 1+cmap2(Dtable)[:, :, 3]
    rgb = cmap(norm(Ztable))
    
    hexArray = numpy.empty(numpy.shape(opacity), dtype='U7')
    for i, v in numpy.ndenumerate(opacity):
        hexArray[i] = colors.rgb2hex(v*rgb[i][:3])
    
    jsondata['MeshBest']['hexArray'] = base64.b64encode(hexArray)

    
    for (j, i) in numpy.ndindex((row, col)):
        if Ztable[j, i] > 0:
            if (j, i + 1) in numpy.ndindex((row, col)):
                if Ztable[j, i + 1] != Ztable[j, i]:
                    line = plt.Line2D((i + 1.5, i + 1.5), (j + 0.5, j + 1.5), lw=3, color='white')
                    plt.gca().add_line(line)
            if (j, i - 1) in numpy.ndindex((row, col)):
                if Ztable[j, i - 1] != Ztable[j, i]:
                    line = plt.Line2D((i + 0.5, i + 0.5), (j + 0.5, j + 1.5), lw=3, color='white')
                    plt.gca().add_line(line)
            if (j + 1, i) in numpy.ndindex((row, col)):
                if Ztable[j + 1, i] != Ztable[j, i]:
                    line = plt.Line2D((i + 0.5, i + 1.5), (j + 1.5, j + 1.5), lw=3, color='white')
                    plt.gca().add_line(line)
            if (j - 1, i) in numpy.ndindex((row, col)):
                if Ztable[j - 1, i] != Ztable[j, i]:
                    line = plt.Line2D((i + 0.5, i + 1.5), (j + 0.5, j + 0.5), lw=3, color='white')
                    plt.gca().add_line(line)
    
    
    if addPositions:
        
        try:
            EllipseArray = numpy.fromstring(base64.b64decode(jsondata['MeshBest']['EllipseArray']))
            EllipseArray = numpy.reshape(EllipseArray, (numpy.size(EllipseArray)/6, 6))
            for line in xrange(numpy.size(EllipseArray)/6):
                edgecolor = 'orange'
                linewidth = 2
                ax.add_patch(El((EllipseArray[line, 0], EllipseArray[line, 1]), EllipseArray[line, 3],
                            EllipseArray[line, 2], EllipseArray[line, 4], linewidth=linewidth,
                            fill=False, edgecolor=edgecolor, zorder=5))
        except KeyError:
            try:
                BestPositions = numpy.fromstring(base64.b64decode(jsondata['MeshBest']['BestPositions']))
                BestPositions = numpy.reshape(BestPositions, (numpy.size(BestPositions)/4, 4))
                
                for position in BestPositions:
                    ax.add_artist(Circle((position[0], position[1]), position[2]/2.0, zorder=5, clip_on=False,
                                         linewidth=3, edgecolor='orange', fill=False))
            except KeyError:
                logger.error('addPositions error: Best Positions are not communicated in the JSON')
    




def LinePlot(jsondata):
    
    try:
        row, col = jsondata['grid_info']['steps_y'], jsondata['grid_info']['steps_x']
        difminpar = jsondata['MeshBest']['difminpar']
    except KeyError:
        logger.error('Experiment parameters are not communicated in the JSON')
        return None
    try:
        Dtable = numpy.fromstring(base64.b64decode(jsondata['MeshBest']['Dtable']))
        Dtable = numpy.reshape(Dtable, (row, col))
        Ztable = numpy.fromstring(base64.b64decode(jsondata['MeshBest']['Ztable']))
        Ztable = numpy.reshape(Ztable, (row, col))
    except KeyError:
        logger.error('Plotting: No data to work with in the JSON')
        return None
    
    basecolors = ['#FF0101', '#F5A26F', '#668DE5', '#E224DE', '#04FEFD', '#00CA02', '#FEFE00', '#0004AF', '#B5FF06']

    Zunique = numpy.unique(Ztable[Ztable>0])
    clrs = []
    for cycle in xrange(1+len(Zunique)/9):
        clrs = clrs + basecolors

    for index in xrange(len(Zunique)):
        plt.plot(Dtable*(Ztable==Zunique[index]), color=clrs[index])

