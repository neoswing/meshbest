'''
Created on May 15, 2018

@author: melnikov
'''


import numpy
import ctypes
import multiprocessing as mp
import base64
from scipy.optimize import differential_evolution

try:
    from workflow_lib import workflow_logging
    logger = workflow_logging.getLogger()
except:
    import logging
    logger = logging.getLogger("MeshBest")



def Ellipsoid((x, y), X0, Y0, a, b, fi, h):
    fi = fi * 3.14159 / 180
    undersqrt = 1 - ((x - X0) * numpy.cos(fi) - (y - Y0) * numpy.sin(fi)) ** 2 / (a ** 2) - \
                ((x - X0) * numpy.sin(fi) + (y - Y0) * numpy.cos(fi)) ** 2 / (b ** 2)
    undersqrt[undersqrt < 0] = 0
    result = h * numpy.sqrt(undersqrt)
    return result.ravel()


def Residue(p, x, y, stretched_submatrix):
    X0, Y0, a, b, fi, h = p[0], p[1], p[2], p[3], p[4], p[5]
    
    if a<b:
        temp = b
        b = a
        a = temp
    
#    Penalty = (4*3.14*a*b*h/0.03) #Full ellipsoid volume**2

    E = Ellipsoid((x, y), X0, Y0, a, b, fi, h).reshape(numpy.shape(stretched_submatrix))
    MASK = (stretched_submatrix == 0)
    difarray = E - stretched_submatrix
    Focus1 = (numpy.sqrt(a**2-b**2)*numpy.cos(fi)+X0, numpy.sqrt(a**2-b**2)*numpy.sin(fi)+Y0)
    Focus2 = (-numpy.sqrt(a**2-b**2)*numpy.cos(fi)+X0, -numpy.sqrt(a**2-b**2)*numpy.sin(fi)+Y0)
    
    limits = numpy.shape(stretched_submatrix)
#    if Focus1[0]>10 and Focus1[0]<limits[1]-10 and Focus1[1]>10 and Focus1[1]<limits[0]-10 and Focus2[0]>10 and Focus2[0]<limits[1]-10 and Focus2[1]>10 and Focus2[1]<limits[0]-10:
#        Penalty = 0

    return numpy.sum(((difarray) * (5 + MASK)) ** 2)



def FitEllipse_MP(queue):

    while True:
        Value = queue.get()
        if Value == None:
            break

        where = numpy.where(improvedZ == Value)
#        logger.debug('---', where, '---')
        if len(where[0]) * len(where[1]) == 1:
    # simple case of single image
            xy = (float(where[1]), float(where[0]))
            width = 1
            height = 1
            angle = 0
            h = Dtable[int(xy[1]) - 1, int(xy[0]) - 1]
            Buffer[(Value - 1) * 7:Value * 7] = [Value, xy[0], xy[1], height, width, angle, h]

        else:
            dx = (numpy.min(where[1]), numpy.max(where[1]))
            dy = (numpy.min(where[0]), numpy.max(where[0]))
            submatrix = improvedD[dy[0] - 1:dy[1] + 2, dx[0] - 1:dx[1] + 2]
            submatrixZ = improvedZ[dy[0] - 1:dy[1] + 2, dx[0] - 1:dx[1] + 2]

            stretched_submatrix = numpy.zeros((numpy.shape(submatrix)[0] * 10, numpy.shape(submatrix)[1] * 10))

            for index in numpy.ndindex(numpy.shape(stretched_submatrix)):
                j, i = index
                j = int(float(j) / 10)
                i = int(float(i) / 10)

                if submatrixZ[j, i] == Value:
                    stretched_submatrix[index] = submatrix[j, i]
#                elif submatrixZ[j, i]==-1:
#                    stretched_submatrix[index] = 0
                else:
                    stretched_submatrix[index] = 0




            x = numpy.arange(0, numpy.shape(submatrix)[1], 0.1)
            y = numpy.arange(0, numpy.shape(submatrix)[0], 0.1)
            x, y = numpy.meshgrid(x, y)


            a = float(max(dx[1] - dx[0], dy[1] - dy[0]) + 1)
            b = float(min(dx[1] - dx[0], dy[1] - dy[0]) + 1)
            c = numpy.sqrt(a ** 2 + b ** 2)
            p0 = numpy.array([ (0.5, numpy.shape(submatrix)[1]-1.5), (0.5, numpy.shape(submatrix)[0]-1.5), (b / 2, c / 1), (0.1, b), (-90, 90), (numpy.mean(submatrix), numpy.max(submatrix))])

            xopt = differential_evolution(Residue, p0, args=(x, y, stretched_submatrix), tol=0.1)

    # creating ellipse

            xy = (dx[0] - 1.5 + xopt.x[0], dy[0] - 1.5 + xopt.x[1])

            if xopt.x[3] <= xopt.x[2]:
                width = xopt.x[3] * 2
                height = xopt.x[2] * 2
                angle = -xopt.x[4] + 90
#                logger.debug(height, width)
            else:
                width = xopt.x[2] * 2
                height = xopt.x[3] * 2
                angle = -xopt.x[4]
#                logger.debug(height, width)
            h = xopt.x[5]

            Buffer[(Value - 1) * 7:Value * 7] = [Value, xy[0], xy[1], height, width, angle, h]




def DoEllipseFit(jsondata):
    global improvedZ, improvedD, Ztable, Dtable, Buffer
    try:
        Dtable = jsondata['MeshBest']['Dtable']
        Ztable = jsondata['MeshBest']['Ztable']
    except KeyError:
        logger.error('jsondata misses precedent MeshBest steps')
        return None

    
    improvedZ = numpy.zeros((numpy.shape(Ztable)[0] + 2, numpy.shape(Ztable)[1] + 2))
    improvedD = numpy.zeros((numpy.shape(Ztable)[0] + 2, numpy.shape(Ztable)[1] + 2))
    improvedZ[:, :] = -1
    improvedZ[1:-1, 1:-1] = Ztable
    improvedD[1:-1, 1:-1] = Dtable

    Buffer = mp.RawArray(ctypes.c_double, int(numpy.max(Ztable)) * 7)
    nCPU = mp.cpu_count()
    queue = mp.Queue()
    
    for Value in range(1, int(numpy.max(Ztable)) + 1):
        queue.put(Value)
    for item in xrange(nCPU):
        queue.put(None)

    workers = []
    for item in xrange(nCPU):
        worker = mp.Process(target=FitEllipse_MP, args=(queue,))
        workers.append(worker)
        worker.start()
    for worker in workers:
        worker.join()

    temp = numpy.frombuffer(Buffer)

    temp = temp.reshape(int(numpy.max(Ztable)), 7)

    
    EllipseArray = numpy.array(temp)[:, 1:]
    EllipseArray[:, 4] = numpy.array(temp)[:, 5]
    EllipseArray[:, 5] = numpy.pi*temp[:, 3]*temp[:, 4]*temp[:, 6]/6.0
    EllipseArray = EllipseArray[EllipseArray[:, 5].argsort()][::-1]

    jsondata['MeshBest']['EllipseArray'] = base64.b64encode(numpy.ascontiguousarray(EllipseArray))
    numpy.savetxt('Result_Ellipses.txt', EllipseArray, fmt='%0.2f')
    
    BestPositions = numpy.array(EllipseArray)[:, :4]
    BestPositions[:, 2] = EllipseArray[:, 3]
    BestPositions[:, 3] = EllipseArray[:, 5]*EllipseArray[:, 3]/EllipseArray[:, 2]
    
    BestPositions = BestPositions[BestPositions[:, 3].argsort()][::-1]
    
    jsondata['MeshBest']['BestPositions'] = base64.b64encode(numpy.ascontiguousarray(BestPositions))
    numpy.savetxt('Result_BestPositions.txt', BestPositions, fmt='%0.2f')
