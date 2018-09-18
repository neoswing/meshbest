'''
Created on May 15, 2018

@author: melnikov
'''


import numpy
import base64
from scipy import signal
import multiprocessing as mp



try:
    from workflow_lib import workflow_logging
    logger = workflow_logging.getLogger()
except:
    import logging
    logger = logging.getLogger("MeshBest")
    
    
    

def SpotMatrix(n):

    profile = numpy.ones((2*n+1, 2*n+1))

    for index in numpy.ndindex(numpy.shape(profile)):
        j, i = index
        if ((i-n)**2+(j-n)**2)>=n**2:
            profile[index] = 0
    return profile/numpy.sum(profile)**(0.75)

def SpotMatrix_t(n1, n2):

    profile = numpy.ones((2*n1+1, 2*n2+1))

    for index in numpy.ndindex(numpy.shape(profile)):
        j, i = index
        if (((i-n1)/n1)**2+((j-n2)/n2)**2)>=1:
            profile[index] = 0
    return profile/numpy.sum(profile)**(0.75)



def GenerateGaussian2D(x, y, mu, sigma):
    gauss = 1/(6.28*sigma**2)*numpy.exp(-( ((x-mu[0])**2+(y-mu[1])**2) / ( 2.0 * sigma**2 ) ) )
    return gauss



def BestCorr_MP(queue):
    

    dx = size_x*1000
    dy = size_y*1000
    
    mind = min(dx, dy)
    maxd = max(dx, dy)

    
    Set = [float(x) for x in AvailableApertures if x>=mind]

    while True:
        Value = queue.get()
        if Value == None:
            break

        
        where = numpy.where(Ztable == Value)
        
        if len(where[0]) * len(where[1]) == 1:
            # simple case of single image
            
            xy = (float(where[1]) + 1, float(where[0]) + 1)
            size = 1.0
            score = (3.14/4)*mind**2*Dtable[where[0], where[1]]
            
            ResultArray = numpy.array([[xy[0], xy[1], size, score]])
            
            Buffer[int(Value)] = base64.b64encode(ResultArray)
            
        else:
            Xlimits = (numpy.min(where[1]), numpy.max(where[1]))
            Ylimits = (numpy.min(where[0]), numpy.max(where[0]))
            
            submatrix = numpy.zeros((Ylimits[1] - Ylimits[0] + 3, Xlimits[1] - Xlimits[0] + 3))
            submatrix[1:-1, 1:-1] = Dtable[Ylimits[0]:Ylimits[1] + 1, Xlimits[0]:Xlimits[1] + 1] * \
                                (Ztable[Ylimits[0]:Ylimits[1] + 1, Xlimits[0]:Xlimits[1] + 1]==Value)

    
            stretched_submatrix = numpy.zeros((numpy.shape(submatrix)[0] * 10, numpy.shape(submatrix)[1] * 10))
            
            for index in numpy.ndindex(numpy.shape(stretched_submatrix)):
                j, i = index
                j = int(float(j) / int(10*dy/mind))
                i = int(float(i) / int(10*dx/mind))
                stretched_submatrix[index] = submatrix[j, i]
                
            silhouette = numpy.array(stretched_submatrix>difminpar, dtype=float)
            
            x = numpy.arange(numpy.shape(stretched_submatrix)[1])
            y = numpy.arange(numpy.shape(stretched_submatrix)[0])
            x, y = numpy.meshgrid(x, y)


            #------------------------------------------------------------------------------------   
            Correlations = {}
            gaussparameter = 0.25
            peaks = []
            for i in xrange(2):
                maximums = []
                
                for aperture in Set:

                    n = int(5*aperture/mind)
                    Corr =  signal.convolve2d(silhouette, SpotMatrix(n), mode='same')

                    Correlations[aperture] = Corr

                    maximums.append(numpy.max(Corr))
#                    print aperture
#                    print numpy.max(Corr)
#                    plt.imshow(Corr, interpolation='nearest', cmap='hot')
#                    plt.colorbar()
#                    plt.show()
    
                factor = 10.0
                #peaks (X, Y, size)

                while True:
                    num = maximums.index(max(maximums))
                    factor = maximums[num]

                    if factor<2.5:
                        break
                    Peak = numpy.where(Correlations[Set[num]]==maximums[num])
                    peaks.append((Peak[1][0], Peak[0][0], Set[num]))
    
                    silhouette = silhouette - ((x-Peak[1][0])**2+(y-Peak[0][0])**2<=((10/mind)*Set[num]/2.0)**2)
                    silhouette[silhouette!=1] = 0
#                    plt.imshow(silhouette, interpolation='nearest', cmap='hot')
#                    plt.colorbar()
#                    plt.show()
                    
                    for aperture in Set:
                        
                        subtr = maximums[num]*6.28*((Set[num]/(gaussparameter*mind))**2)*GenerateGaussian2D(x, y, (Peak[1][0], Peak[0][0]), sigma=Set[num]/(gaussparameter*mind))
                        Correlations[aperture] = Correlations[aperture]-subtr
    
                        maximums.append(numpy.max(Correlations[aperture]))
#                        plt.imshow(Correlations[aperture], interpolation='nearest', cmap='hot')
#                        plt.colorbar()
#                        plt.show()
                    
                    maximums = maximums[len(Set):]
    
                
                
                
            #-----------------------------------------------------------------------------

            
            ResultArray = numpy.zeros((len(peaks), 4))
            for n in xrange(len(peaks)):
                X = peaks[n][0]/(10.0*(dx/mind)) + Xlimits[0] - 0.5
                Y = peaks[n][1]/(10.0*(dy/mind)) + Ylimits[0] - 0.5
                size = peaks[n][2]/mind

                ResultArray[n, :3] = numpy.array([X, Y, size])
                ResultArray[n, 3] = (peaks[n][2]/10.0)**2*numpy.sum(numpy.multiply(stretched_submatrix, ((x-peaks[n][0])**2+(y-peaks[n][1])**2<=((10/mind)*peaks[n][2]/2.0)**2)))
                if Value==-2:
                    ResultArray[n, 3] = ResultArray[n, 3]/100.0
                
            Buffer[int(Value)] = base64.b64encode(ResultArray)



    
def GetAllPositions(jsondata):
    global Buffer, Dtable, Ztable, difminpar, AvailableApertures, size_x, size_y
    
    try:
        Dtable = jsondata['MeshBest']['Dtable']
        Ztable = jsondata['MeshBest']['Ztable']
    except KeyError:
        logger.error('jsondata misses precedent MeshBest steps')
        return None

    difminpar = jsondata['MeshBest']['difminpar']
    size_x = jsondata['grid_info']['beam_width']
    size_y = jsondata['grid_info']['beam_height']

    AvailableApertures = jsondata['beamlineInfo']['beamlineApertures']


    
    manager = mp.Manager()
    Buffer = manager.dict()
    nCPU = mp.cpu_count()
    queue = mp.Queue()

    for Value in numpy.unique(Ztable[Ztable!=-1]):
        queue.put(Value)
    for item in xrange(nCPU):
        queue.put(None)

    workers = []
    for item in xrange(nCPU):
        worker = mp.Process(target=BestCorr_MP, args=(queue,))
        workers.append(worker)
        worker.start()
    for worker in workers:
        worker.join()
    
    
    BestPositions = numpy.empty((0, 4), float)
    
    
    for Value in numpy.unique(Ztable[Ztable!=-1]):
        ApertureCorrArray = numpy.fromstring(base64.b64decode(Buffer[Value]))
        ApertureCorrArray = ApertureCorrArray.reshape(len(ApertureCorrArray)/4, 4)
        BestPositions = numpy.append(BestPositions, ApertureCorrArray, axis=0)
#        print '*******\n', BestPositions, '*******\n'

    BestPositions = BestPositions[BestPositions[:, 3].argsort()][::-1]
    
    jsondata['MeshBest']['BestPositions'] = base64.b64encode(numpy.ascontiguousarray(BestPositions))
    numpy.savetxt('Result_BestPositions.txt', BestPositions, fmt='%0.2f')
