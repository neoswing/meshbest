'''
Created on May 15, 2018

@author: melnikov
'''
import numpy
from scipy import ndimage
import scipy.cluster.hierarchy as hca
import base64
from matplotlib import pyplot as plt
import multiprocessing as mp
import ctypes


try:
    from workflow_lib import workflow_logging
    logger = workflow_logging.getLogger()
except:
    import logging
    logger = logging.getLogger("MeshBest")



def From64ToSpotArray(string64):
    try:
        array = numpy.fromstring(base64.b64decode(string64))
        array = numpy.reshape(array, (numpy.size(array) / 5, 5))
    except TypeError:
        array = numpy.array([])
        print string64
    return array


def DistanceCalc_MP(queue):
    while True:
        string = queue.get()
        if string == None:
            break
#        logger.debug('MP_checkpoint_1')
        string1, string2, bufcoord = string
        
        if string1 != '0' and string2 != '0':
#            logger.debug('MP_checkpoint_2')
            array1 = From64ToSpotArray(string1)
            array2 = From64ToSpotArray(string2)
            
            a = numpy.size(array1) / 5
            b = numpy.size(array2) / 5
            
            RealCoords1 = numpy.zeros((a, 3))
            RealCoords2 = numpy.zeros((b, 3))

            for i in xrange(a):
                x = (array1[i, 1] - BeamCenter[0]) * DetectorPixel
                y = (array1[i, 2] - BeamCenter[1]) * DetectorPixel
                divider = numpy.sqrt(x ** 2 + y ** 2 + DetectorDistance ** 2)
                RealCoords1[i, 0] = x / divider
                RealCoords1[i, 1] = y / divider
                RealCoords1[i, 2] = DetectorDistance / divider

            for i in xrange(b):
                x = (array2[i, 1] - BeamCenter[0]) * DetectorPixel
                y = (array2[i, 2] - BeamCenter[1]) * DetectorPixel
                divider = numpy.sqrt(x ** 2 + y ** 2 + DetectorDistance ** 2)
                RealCoords2[i, 0] = x / divider
                RealCoords2[i, 1] = y / divider
                RealCoords2[i, 2] = DetectorDistance / divider


            output = []
            
            thrsh = 0.1*3.14159/180.0
            
            if a >= b:
                for i in xrange(b):
                    tarr = RealCoords1[:, :] - RealCoords2[i, :]
                    X = numpy.min(numpy.sum(numpy.multiply(tarr, tarr), axis=1))
                    X = min(X, thrsh**2)
                    output.append(X)
            else:
                for i in xrange(a):
                    tarr = RealCoords2[:, :] - RealCoords1[i, :]
                    X = numpy.min(numpy.sum(tarr**2, axis=1))
                    X = min(X, thrsh**2)
                    output.append(X)
            
            output = numpy.array(output)
            F = numpy.sqrt(numpy.mean(output))*180/3.14159
            
            Buffer[int(bufcoord)] = F

        else:
            logger.error('Problem of referencing - spot lists are not found')



def CalculateZone(jsondata, keys):
    global Buffer
    L = len(keys[0])

    Buffer = mp.RawArray(ctypes.c_double, L * L)
    nCPU = mp.cpu_count()
    queue = mp.Queue()
    positionReference = jsondata['MeshBest']['positionReference']
    
    for i in xrange(L):
        for j in range(i + 1, L):
            keyi, keyj = positionReference[keys[0][i], keys[1][i]], positionReference[keys[0][j], keys[1][j]]
            if numpy.abs(jsondata['meshPositions'][keyi]['omega'] - \
                         jsondata['meshPositions'][keyj]['omega'])>=0.5:
                Buffer[j + L * i] = numpy.nan
            else:
                try:
                    array1 = jsondata['meshPositions'][keyi]['dozorSpotList_saltremoved']
                except KeyError:
                    array1 = jsondata['meshPositions'][keyi]['dozorSpotList']
                except KeyError:
                    array1 = '0'
                try:
                    array2 = jsondata['meshPositions'][keyj]['dozorSpotList_saltremoved']
                except KeyError:
                    array2 = jsondata['meshPositions'][keyj]['dozorSpotList']
                except KeyError:
                    array2 = '0'
                query = [array1, array2, str(j + L * i)]
                queue.put(query)
            

    for i in xrange(nCPU):
        queue.put(None)
    workers = []
    for item in xrange(nCPU):
        worker = mp.Process(target=DistanceCalc_MP, args=(queue,))
        workers.append(worker)
        worker.start()
    for worker in workers:
        worker.join()


    T = numpy.frombuffer(Buffer)
    T = numpy.reshape(T, (L, L))
    
    Buffer = 0
    
    T = T+numpy.transpose(T)
    for key, value in numpy.ndenumerate(T):
        if value!=value:
            ar = numpy.maximum(T[key[0], :], T[:, key[1]])
            if numpy.all(numpy.isnan(ar)):
                T[key] = 0.1
            else:
                T[key] = numpy.nanmean(numpy.maximum(T[key[0], :], T[:, key[1]]))
    
    
    q = []
    for key in numpy.ndenumerate(T):
        if key[0][1] > key[0][0]:
            q.append(key[1])
    T = numpy.array(q)

    Linkages = hca.linkage(T, method='average')
    
    hca.dendrogram(Linkages)
    
    name = 'Dendrogram%d.png' % L
    plt.savefig(name, dpi=300, transparent=True, bbox_inches='tight')  # , pad_inches=0)
    plt.close()

    Clusters = hca.fcluster(Linkages, t=0.093, criterion='distance')
    Clusters = Clusters.astype('a5')


    ClusterIndex = numpy.unique(Clusters, return_inverse=True)

    Base = numpy.max(Ztable)+1
    for i in xrange(len(Clusters)):
        Ztable[keys[0][i], keys[1][i]] = ClusterIndex[1][i] + Base

def PerformCrystalRecognition(jsondata):
    global Wavelength, DetectorDistance, BeamCenter, DetectorPixel, Ztable
    
    Wavelength = jsondata['inputDozor']['wavelength']
    DetectorDistance = jsondata['inputDozor']['detectorDistance'] * 1000
    BeamCenter = (jsondata['inputDozor']['orgx'], jsondata['inputDozor']['orgy'])
    DetectorPixel = jsondata['beamlineInfo']['detectorPixelSize']*1000
    
    Ztable = jsondata['MeshBest']['Ztable']

    if numpy.size(Ztable[Ztable == 0]) == 1:
        Ztable[Ztable == 0] = 1
        
    elif numpy.size(Ztable[Ztable == 0]) > 1:

        ZoneMap = ndimage.measurements.label((Ztable > -1), structure=numpy.ones((3, 3)))
    
        for i in range(1, 1+ZoneMap[1]):
            keys = numpy.where(ZoneMap[0]==i)
            if len(keys[0])==1:
                M = numpy.max(Ztable)
                Ztable[keys[0][0], keys[1][0]] = M+1
            else:
                CalculateZone(jsondata, keys)
    else:
        logger.info('Only multi-crystal diffraction found')


    jsondata['MeshBest']['Ztable'] = Ztable
    numpy.savetxt('Ztable.txt', Ztable, fmt='%d')
