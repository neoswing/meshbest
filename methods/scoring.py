'''
Created on May 15, 2018

@author: melnikov
'''
import numpy
from scipy import ndimage, spatial
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
    array = numpy.fromstring(base64.b64decode(string64))
    array = array.reshape((int(array.size/5), 5))

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
            
            count1 = int(a)#1000 if a>1000 else int(a)
            count2 = int(b)#1000 if b>1000 else int(b)
            
            
#            RealCoords1 = None
#            RealCoords2 = None
#            
#            if len(numpy.atleast_1d(array1)) > 5:
#                RealCoords1 = numpy.zeros((numpy.shape(array1)[0], 5))
#        
#        
#                x = (array1[:, 1] - BeamCenter[0]) * DetectorPixel
#                y = (array1[:, 2] - BeamCenter[1]) * DetectorPixel
#                divider = numpy.sqrt(x ** 2 + y ** 2 + DetectorDistance ** 2)
#                RealCoords1[:, 0] = x / divider
#                RealCoords1[:, 1] = y / divider
#                RealCoords1[:, 2] = 1 - DetectorDistance/divider
#                RealCoords1[:, 3] = array1[:, 0]
#                RealCoords1[:, 4] = array1[:, 3]# / array1[:, 4]
#        
#            if len(numpy.atleast_1d(array2)) > 5:
#                RealCoords2 = numpy.zeros((numpy.shape(array2)[0], 5))
#        
#        
#                x = (array2[:, 1] - BeamCenter[0]) * DetectorPixel
#                y = (array2[:, 2] - BeamCenter[1]) * DetectorPixel
#                divider = numpy.sqrt(x ** 2 + y ** 2 + DetectorDistance ** 2)
#                RealCoords2[:, 0] = x / divider
#                RealCoords2[:, 1] = y / divider
#                RealCoords2[:, 2] = 1 - DetectorDistance/divider
#                RealCoords2[:, 3] = array2[:, 0]
#                RealCoords2[:, 4] = array2[:, 3]# / array2[:, 4]
#            
#            
#            
#            if isinstance(RealCoords1, numpy.ndarray) and isinstance(RealCoords2, numpy.ndarray):
#            
#                c = spatial.distance.cdist(RealCoords1[:, :3], RealCoords2[:, :3], metric='euclidean')
#            
#                I = RealCoords1[:, 4, numpy.newaxis] * RealCoords2[:, 4, numpy.newaxis].T
#                
#                thr = 0.1*3.14/180.0
#                
#                
#                F = numpy.sum(I[c<thr])
#            else:
#                F = 0.0
            
            
            
            
            RealCoords1 = numpy.zeros((count1, 3))
            RealCoords2 = numpy.zeros((count2, 3))
            

            x = (array1[:, 1] - BeamCenter[0]) * DetectorPixel
            y = (array1[:, 2] - BeamCenter[1]) * DetectorPixel
            divider = numpy.sqrt(x ** 2 + y ** 2 + DetectorDistance ** 2)
            RealCoords1[:, 0] = x / divider
            RealCoords1[:, 1] = y / divider
            RealCoords1[:, 2] = 1 - DetectorDistance/divider


            x = (array2[:, 1] - BeamCenter[0]) * DetectorPixel
            y = (array2[:, 2] - BeamCenter[1]) * DetectorPixel
            divider = numpy.sqrt(x ** 2 + y ** 2 + DetectorDistance ** 2)
            RealCoords2[:, 0] = x / divider
            RealCoords2[:, 1] = y / divider
            RealCoords2[:, 2] = 1 - DetectorDistance/divider
                

            output = []
            
            thrsh = 0.1*3.14159/180.0
            
#            logger.debug('start calculating distance matrix')
            
            DistanceMatrix = spatial.distance.cdist(RealCoords1, RealCoords2, metric='euclidean')

            if a >= b:
                output = numpy.min(DistanceMatrix, axis=0)
            else:
                output = numpy.min(DistanceMatrix, axis=1)
            
            
            output[output>thrsh] = thrsh

            
            output = numpy.array(output)

            F = numpy.sqrt(numpy.mean(output**2))*180/3.14159

            Buffer[int(bufcoord)] = F

#            logger.debug('finish calculating distance matrix')
        else:
            logger.error('Problem of referencing - spot lists are not found')

def C_Calc_MP(queue):
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
            
            count1 = 1000 if a>1000 else int(a)
            count2 = 1000 if b>1000 else int(b)
            
            
            RealCoords1 = None
            RealCoords2 = None
            
            if len(numpy.atleast_1d(array1)) > 5:
                RealCoords1 = numpy.zeros((numpy.shape(array1)[0], 5))
        
        
                x = (array1[:, 1] - BeamCenter[0]) * DetectorPixel
                y = (array1[:, 2] - BeamCenter[1]) * DetectorPixel
                divider = numpy.sqrt(x ** 2 + y ** 2 + DetectorDistance ** 2)
                RealCoords1[:, 0] = x / divider
                RealCoords1[:, 1] = y / divider
                RealCoords1[:, 2] = 1 - DetectorDistance/divider
                RealCoords1[:, 3] = array1[:, 0]
                RealCoords1[:, 4] = array1[:, 3]# / array1[:, 4]
        
            if len(numpy.atleast_1d(array2)) > 5:
                RealCoords2 = numpy.zeros((numpy.shape(array2)[0], 5))
        
        
                x = (array2[:, 1] - BeamCenter[0]) * DetectorPixel
                y = (array2[:, 2] - BeamCenter[1]) * DetectorPixel
                divider = numpy.sqrt(x ** 2 + y ** 2 + DetectorDistance ** 2)
                RealCoords2[:, 0] = x / divider
                RealCoords2[:, 1] = y / divider
                RealCoords2[:, 2] = 1 - DetectorDistance/divider
                RealCoords2[:, 3] = array2[:, 0]
                RealCoords2[:, 4] = array2[:, 3]# / array2[:, 4]
            
            
            
            if isinstance(RealCoords1, numpy.ndarray) and isinstance(RealCoords2, numpy.ndarray):
            
                c = spatial.distance.cdist(RealCoords1[:, :3], RealCoords2[:, :3], metric='euclidean')
            
                I = RealCoords1[:, 4, numpy.newaxis] * RealCoords2[:, 4, numpy.newaxis].T
                div = numpy.sqrt(numpy.sum(RealCoords1[:, 4]**2) * numpy.sum(RealCoords2[:, 4]**2))
                
                thr = 0.1*3.14/180.0
                
                
                F = numpy.sum(I[c<thr])/div
            else:
                F = 0.0
            
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
#    logger.debug('start')
    for i in range(L):
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
            
#    logger.debug('filled')
    for i in range(nCPU):
        queue.put(None)
    workers = []
    for item in range(nCPU):
        worker = mp.Process(target=DistanceCalc_MP, args=(queue,))
        workers.append(worker)
        worker.start()
    for worker in workers:
        worker.join()
    
#    logger.debug('calculated distances')

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

#    logger.debug('constructed distance matrix')



    plt.imshow(T, interpolation='nearest', cmap='hot')
    plt.colorbar()
    plt.savefig('Dmatrix.png', dpi=300, bbox_inches='tight')
    plt.close()
    
    numpy.savetxt('Dmatrix.txt', T, fmt='%.3f')

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
    for i in range(len(Clusters)):
        Ztable[keys[0][i], keys[1][i]] = ClusterIndex[1][i] + Base
        
#    logger.debug('finished clustering')

def Cmatrixcalc(jsondata, keys):
    global Buffer
    L = len(keys[0])

    Buffer = mp.RawArray(ctypes.c_double, L * L)
    nCPU = mp.cpu_count()
    queue = mp.Queue()
    positionReference = jsondata['MeshBest']['positionReference']
#    logger.debug('start')
    for i in range(L):
        for j in range(i, L):
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
            
#    logger.debug('filled')
    for i in range(nCPU):
        queue.put(None)
    workers = []
    for item in range(nCPU):
        worker = mp.Process(target=C_Calc_MP, args=(queue,))
        workers.append(worker)
        worker.start()
    for worker in workers:
        worker.join()
    
#    logger.debug('calculated distances')

    T = numpy.frombuffer(Buffer)
    T = numpy.reshape(T, (L, L))

    Buffer = 0
    
    T = T+numpy.transpose(T)-numpy.diag(numpy.diag(T))
    for key, value in numpy.ndenumerate(T):
        if value!=value:
            ar = numpy.maximum(T[key[0], :], T[:, key[1]])
            if numpy.all(numpy.isnan(ar)):
                T[key] = 0.1
            else:
                T[key] = numpy.nanmean(numpy.maximum(T[key[0], :], T[:, key[1]]))

#    logger.debug('constructed distance matrix')



    plt.imshow(T, interpolation='nearest', cmap='hot')
    plt.colorbar()
    plt.savefig('Cmatrix.png', dpi=300, bbox_inches='tight')
    plt.close()
    
    numpy.savetxt('Cmatrix.txt', T, fmt='%.3f')

        
#    logger.debug('finished clustering')


def PerformCrystalRecognition(jsondata):
    global Wavelength, DetectorDistance, BeamCenter, DetectorPixel, Ztable
    
    Wavelength = jsondata['wavelength']
    DetectorDistance = jsondata['detector_distance'] * 1000
    BeamCenter = (jsondata['orgx'], jsondata['orgy'])
    DetectorPixel = jsondata['detectorPixelSize']*1000
    
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

def getCmatrix(jsondata):
    global Wavelength, DetectorDistance, BeamCenter, DetectorPixel, Ztable
    
    Wavelength = jsondata['wavelength']
    DetectorDistance = jsondata['detector_distance'] * 1000
    BeamCenter = (jsondata['orgx'], jsondata['orgy'])
    DetectorPixel = jsondata['detectorPixelSize']*1000
    
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
                Cmatrixcalc(jsondata, keys)
    else:
        logger.info('Only multi-crystal diffraction found')



