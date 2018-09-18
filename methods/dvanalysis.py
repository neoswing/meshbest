'''
Created on May 15, 2018

@author: melnikov
'''

import numpy
from scipy import stats, signal
import base64
import random
import matplotlib
from matplotlib import pyplot as plt
import multiprocessing as mp
import ctypes

try:
    from workflow_lib import workflow_logging
    logger = workflow_logging.getLogger()
except:
    import logging
    logger = logging.getLogger("MeshBest")





def triangle(x0, y0, length):
    x = numpy.linspace(0, 99, 100)
    array = y0 - 2 * y0 * numpy.abs(x - x0) / length
    array = array * (array > 0)
    return array


def From64ToSpotArray(string64):
    array = numpy.fromstring(base64.b64decode(string64))
    array = numpy.reshape(array, (numpy.size(array) / 5, 5))
    return array


def AMPDiter(array):


    L = int(len(array) / 2)
    matrix = numpy.zeros((L, len(array)))

    for k in range(1, L + 1):
        for i in range(1, len(array) + 1):

            if i >= k + 2 and i < len(array) - k + 2:
#                W = 2 * k
                if array[i - 2] > array[i - k - 2] and array[i - 2] > array[i + k - 2]:
                    matrix[k - 1, i - 1] = 0
                else:
                    matrix[k - 1, i - 1] = 1 + random.random() / 2
            else:
                matrix[k - 1, i - 1] = 1 + random.random() / 2
    gammas = numpy.sum(matrix, axis=1)
#    logger.debug(gammas)



    Lambda = numpy.where(gammas == numpy.min(gammas))[0][0] + 1
#    logger.debug(Lambda)

    matrix = matrix[:Lambda, :]
    Sigma = numpy.std(matrix, axis=0)

    peaks = []
    for i in xrange(len(Sigma)):
        if Sigma[i] == 0:
            if (i - 1) / float(len(array)) > 0.00:
                peaks.append(i - 1)
    peaks = numpy.array(peaks, dtype=int)
#    logger.debug('AMPD-result_PEAKS: ', peaks)
    return peaks, Lambda


def AMPD(array_orig):

    fullpeaklist = numpy.array([], dtype=int)

    for cycle in xrange(10):
        M = numpy.mean(array_orig)
#        SD = numpy.std(array_orig)
        X = numpy.arange(0, len(array_orig))
        linfit = stats.linregress(X, array_orig)
        array = array_orig - (linfit[0] * X + linfit[1])
        array = array * (array > 0)
        MAX = numpy.max(array) - M

        allpeaks = numpy.array([], dtype=int)

        while True:
            substract = numpy.zeros(len(array_orig))
            peaks, Lambda = AMPDiter(array)

            peaks = peaks[(array_orig[peaks] - M > MAX / 5)]

            if len(peaks) > 0:
                pass
            else:
                break

            allpeaks = numpy.append(allpeaks, peaks)

            for peak in peaks:
                substract += triangle(peak, array[peak], Lambda)[:len(array_orig)]
            array = array - substract
            array = array * (array > 0)


        if len(numpy.atleast_1d(allpeaks)) == 0:
            break

        allpeaks = numpy.sort(allpeaks)
        allpeaks = allpeaks.astype(int)

        dels = []
        for i in xrange(len(allpeaks)):
            peak = allpeaks[i]

            if peak > 1 and peak < (len(array_orig) - 1):
                if array_orig[peak] < array_orig[peak + 1] or array_orig[peak] < array_orig[peak - 1]:
                    dels.append(i)
        allpeaks = numpy.delete(allpeaks, dels)
        fullpeaklist = numpy.append(fullpeaklist, allpeaks)


    fullpeaklist = numpy.unique(fullpeaklist)
#    fig = plt.plot(array_orig)
#    sc = plt.scatter(fullpeaklist, array_orig[fullpeaklist], color='red')
#    plt.show()
    
    if len(fullpeaklist)>1:
        return fullpeaklist
    else:
        return None

    
    
    
def CalcSeff(array):
#    N = numpy.size(array) / 5
    
    rarray = numpy.sqrt((array[:, 1] - BeamCenter[0])**2+(array[:, 2] - BeamCenter[1])**2)
    rmax = min(int(BeamCenter[0]), int(BeamCenter[1]))
    
    HIST = stats.histogram(rarray, numbins=50, defaultlimits=(0, rmax))
    
    Rspace = numpy.linspace(0, rmax, 50)
    density = numpy.zeros(50)
    for i in xrange(50):
        density[i] = HIST[0][i]/(2*3.14159*HIST[2]*HIST[2]*(i+0.5))
    
    Sdetec = numpy.diff(Rspace**2, 1)
    Sdetec = numpy.append(Sdetec, [0])
    
    Sreal = (DetectorPixel**2)*(Sdetec*DetectorDistance/(Wavelength**2)) / (numpy.sqrt(DetectorDistance**2+(DetectorPixel*Rspace)**2))**3
    Seff = numpy.sum(Sreal*(density>0))


#    plt.plot(density)
#    plt.show()
    return Seff


def SaltRingCheck_MP(queue, BeamCenter, Buffer):

    while True:
        spot = queue.get()
        if spot == None:
            break


        try:
            string = spot['dozorSpotList']
            array = From64ToSpotArray(string)

            if len(numpy.atleast_1d(array)) > 250:
                xsize = int(BeamCenter[0]) * 2 + 1
                ysize = int(BeamCenter[1]) * 2 + 1

                detector = numpy.zeros((ysize, xsize))

                for i in xrange(numpy.shape(array)[0]):
                    detector[int(array[i, 2]), int(array[i, 1])] = 1

                radius_array = numpy.sqrt((array[:, 1] - BeamCenter[0]) ** 2 + (array[:, 2] - BeamCenter[1]) ** 2)
                density = numpy.zeros(numpy.size(radius_array))
                for i in xrange(numpy.shape(array)[0]):
                    x0, y0 = int(array[i, 1]), int(array[i, 2])
                    density[i] = numpy.mean(detector[y0 - 5:y0 + 5, x0 - 5:x0 + 5])

                HIST = stats.histogram(radius_array, numbins=400, defaultlimits=(10, min(BeamCenter)), weights=(density))[0]

#---SALT_RING_CHECK_ANALYSIS---

                M = numpy.mean(HIST)

                # remove_outliers
                X = numpy.arange(400)

                p = numpy.polyfit(X[(HIST < 3 * M)], HIST[X[(HIST < 3 * M)]], 2)
                V = numpy.poly1d(p)

                based = HIST - 2 * V(X) * (V(X) > 0)
                based = based * (based > 0)

                M1 = numpy.mean(based[X[(HIST < 3 * M)]])
                SD = numpy.std(based[X[(HIST < 3 * M)]])
                peaks = X[(based > M1 + 10 * SD) * (HIST > 0.1)]


#---EXCLUDE_BAD_REGIONS---

                Excluded_regions = [((peak - 2) * (min(BeamCenter) - 10) / 399 + 10, (peak + 2) * (min(BeamCenter) - 10) / 399 + 10) for peak in peaks]

                r = numpy.zeros(numpy.shape(array)[0])

                for ring in Excluded_regions:
                    r += (radius_array > ring[0]) * (radius_array < ring[1])



                array = numpy.delete(array, numpy.where(r), 0)

                if len(array) > 5:
                    newstring = base64.b64encode(array)
                    Buffer[spot['index']] = newstring


        except KeyError:
            pass






def MakeHistogram_MP(queue):
    limits = (0.001, 0.04)
    while True:
        spot = queue.get()
        if spot == None:
            break


        if 'dozorSpotList_saltremoved' in spot.keys():
            string = spot['dozorSpotList_saltremoved']
            array = From64ToSpotArray(string)
        elif 'dozorSpotList' in spot.keys():
            string = spot['dozorSpotList']
            array = From64ToSpotArray(string)
        else:
            string = False

        if string != False:
            if len(numpy.atleast_1d(array)) > 5:
                RealCoords = numpy.zeros((numpy.shape(array)[0], 5))

                for i in xrange(numpy.shape(array)[0]):

                    x = (array[i, 1] - BeamCenter[0]) * DetectorPixel
                    y = (array[i, 2] - BeamCenter[1]) * DetectorPixel
                    divider = Wavelength * numpy.sqrt(x ** 2 + y ** 2 + DetectorDistance ** 2)
                    RealCoords[i, 0] = x / divider
                    RealCoords[i, 1] = y / divider
                    RealCoords[i, 2] = DetectorDistance / divider
                    RealCoords[i, 3] = array[i, 0]
                    RealCoords[i, 4] = float(array[i, 3]) / float(array[i, 4])

                array = numpy.array([])

                for i in xrange(len(RealCoords[:, 0])):
                    for j in range(i + 1, len(RealCoords[:, 0])):
                        if numpy.abs(RealCoords[i, 3] - RealCoords[j, 3]) < 15:
                            L = numpy.sqrt((RealCoords[i, 0] - RealCoords[j, 0]) ** 2 + (RealCoords[i, 1] - RealCoords[j, 1]) ** 2 + (RealCoords[i, 2] - RealCoords[j, 2]) ** 2)
                            if L < limits[1] and L >= limits[0]:
                                array = numpy.append(array, L)

                if len(array) > 1:
                    histogr = stats.histogram(array, numbins=100, defaultlimits=limits)
                    string = base64.b64encode(histogr[0])
                    Buffer[spot['index']] = string





def OverlapCheck_MP(queue, base_regions):

    while True:
        item = queue.get()
        if item == None:
            break
        
        Buffer[item['index']] = 0
        
        if item['DVHistogram'] != None:
            try:
                array = From64ToSpotArray(item['dozorSpotList_saltremoved'])
            except KeyError:
                array = From64ToSpotArray(item['dozorSpotList'])
            N = numpy.size(array) / 5
    
            Smax = CalcSeff(array)
            
            H0 = numpy.fromstring(base64.b64decode(item['DVHistogram']))
            H = H0[base_regions]
            X = numpy.linspace(0.001, 0.04, 100)[base_regions]
            
            if numpy.all(H==0):
                k = 0
            else:
                k = numpy.mean(X*H)/numpy.mean(X*X)
                std = numpy.sqrt(numpy.mean((k*X-H)**2))
                weights = numpy.exp((N/450.0)*(k*X-H)/std)
                
                k = numpy.mean(X*H*weights)/numpy.mean(X*X*weights)
            
            
            undersq = N**2-2*k*Smax/(0.039/99.0)
            
            if undersq>0:
                N2est = int(0.5*(N-numpy.sqrt(undersq)))
                N1est = int(0.5*(N+numpy.sqrt(undersq)))
                SatelliteRatio = float(N2est)/float(N1est)
            else:
                N2est = N/2
                N1est = N/2
                SatelliteRatio = 1

            if SatelliteRatio>0.3:
                if N2est>50:
                    Buffer[item['index']] = 1



def EliminateSaltRings(jsondata):
    
    BeamCenter = (jsondata['inputDozor']['orgx'], jsondata['inputDozor']['orgy'])
    
    manager = mp.Manager()
    Buffer = manager.dict()
    queue = mp.Queue()
    nCPU = mp.cpu_count()

    for item in jsondata['meshPositions']:
        queue.put(item)

    for item in xrange(nCPU):
        queue.put(None)

    workers = []
    for item in xrange(nCPU):
        worker = mp.Process(target=SaltRingCheck_MP, args=(queue, BeamCenter, Buffer, ))
        workers.append(worker)
        worker.start()
    for worker in workers:
        worker.join()

    for item in jsondata['meshPositions']:
        try:
            item['dozorSpotList_saltremoved'] = Buffer[item['index']]
        except KeyError:
            pass

    Buffer = 0


def DetermineMCdiffraction(jsondata):
    global Wavelength, DetectorDistance, BeamCenter, DetectorPixel, Buffer
    manager = mp.Manager()
    Buffer = manager.dict()
    nCPU = mp.cpu_count()
    queue = mp.Queue()

    Wavelength = jsondata['inputDozor']['wavelength']
    DetectorDistance = jsondata['inputDozor']['detectorDistance'] * 1000
    BeamCenter = (jsondata['inputDozor']['orgx'], jsondata['inputDozor']['orgy'])
    DetectorPixel = jsondata['beamlineInfo']['detectorPixelSize']*1000
    row, col = jsondata['grid_info']['steps_y'], jsondata['grid_info']['steps_x']


    for item in jsondata['meshPositions']:
        queue.put(item)
    for item in xrange(nCPU):
        queue.put(None)


    workers = []
    for item in xrange(nCPU):
        worker = mp.Process(target=MakeHistogram_MP, args=(queue, ))
        workers.append(worker)
        worker.start()
    for worker in workers:
        worker.join()


    HIST = numpy.zeros(100)
    for item in jsondata['meshPositions']:
        try:
            item['DVHistogram'] = Buffer[item['index']]
        except KeyError:
            item['DVHistogram'] = None

        if item['DVHistogram'] != None:
            HIST += numpy.fromstring(base64.b64decode(item['DVHistogram']))
    
    Buffer = 0

# filtering
    window = signal.general_gaussian(M=5, p=1, sig=1)
    filtered = signal.fftconvolve(HIST, window, mode='same')
    numpy.roll(filtered, -2)

    plt.plot(filtered)
    plt.savefig('Cumulative_DVhistogram_filter.png', dpi=300, transparent=True, bbox_inches='tight')
    plt.close()
    
    jsondata['MeshBest']['Cumulative_DVHistogram'] = base64.b64encode(HIST)

    init_base = AMPD(numpy.max(filtered)-filtered)
    
    if init_base is not None:
        X = numpy.linspace(0.001, 0.04, 100)[init_base]
        init_array = HIST[init_base]
    
        k = numpy.mean(X*init_array)/numpy.mean(X*X)
        
        std = numpy.sqrt(numpy.mean((k*X-init_array)**2))
        weights = numpy.exp((k*X-init_array)/std)
    #    print weights
        k = numpy.mean(X*init_array*weights)/numpy.mean(X*X*weights)
        
        X = numpy.linspace(0.001, 0.04, 100)
        fit = k*X
    
    #    base_regions = (HIST - fit)<0
        base_regions = (HIST - fit)/std<=0.5
        
#check for overlap
        
        Buffer = mp.RawArray(ctypes.c_double, row * col)
    
        for item in jsondata['meshPositions']:
            queue.put(item)
        for item in xrange(nCPU):
            queue.put(None)
    
        workers = []
        for item in xrange(nCPU):
            worker = mp.Process(target=OverlapCheck_MP, args=(queue, base_regions,))
            workers.append(worker)
            worker.start()
        for worker in workers:
            worker.join()
        
        
        MXDiff = numpy.frombuffer(Buffer)
        Buffer = 0
        jsondata['MeshBest']['MXDiff'] = base64.b64encode(MXDiff)
        
        for key, value in numpy.ndenumerate(jsondata['MeshBest']['positionReference']):
            if MXDiff[value]:
                jsondata['MeshBest']['Ztable'][key] = -2
        
    else:
        logger.error('AMPD: Too little data to estimate histogram baseline')
        base_regions = numpy.array([])
        jsondata['MeshBest']['MXDiff'] = base64.b64encode(numpy.zeros((row, col)))
        
    
    
    jsondata['MeshBest']['DVBase_Regions'] = base_regions.tolist()
        
    fig1 = plt.figure()
    
    ax1 = fig1.add_subplot(111)
    X = numpy.linspace(0.001, 0.04, 100)
    plt.plot(X, HIST, 'k')
#    ax1.plot(X, fit)
    if len(numpy.atleast_1d(base_regions)):
        xbase = X[base_regions]
        plt.scatter(xbase, HIST[base_regions], color='red')
    ax1.set_xlim([0.001, 0.04])
    ax1.set_ylim([0.0, 1.1*numpy.max(HIST)])
    plt.title('Histogram of inter-spot distances in the reciprocal space')
    plt.xlabel(r'$d^{-1}$, $\AA$')

    ax1.xaxis.set_major_formatter(matplotlib.ticker.FuncFormatter(lambda x, pos: '{0:0.1f}'.format(1/x) if x!=0 else 0))
    
    plt.savefig('Cumulative_DVHistogram.png', dpi=300, transparent=True, bbox_inches='tight')
    plt.close()
#    plt.show()


