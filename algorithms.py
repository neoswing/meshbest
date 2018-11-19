#!/usr/bin/env python

import numpy
import os
import json
import time
import base64
from scipy import ndimage
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt




from meshbest.methods import jsoncheck, dvanalysis, scoring, plotting, ellipticfit, sizecorr, voxelise

try:
    from workflow_lib import workflow_logging
    logger = workflow_logging.getLogger()
except:
    import logging
    logger = logging.getLogger("MeshBest")
    logger.setLevel(logging.DEBUG)
    ch = logging.StreamHandler()
    ch.setLevel(logging.DEBUG)
    # create formatter
    formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')
    # add formatter to ch
    ch.setFormatter(formatter)
    # add ch to logger
    logger.addHandler(ch)

# numpy.set_printoptions(threshold='nan')

start_time = time.time()





def simple(jsonFilePath, resultsPath=None):
#    returns only crystal map

    if resultsPath!=None:
        os.chdir(resultsPath)
    
    logger.debug('Checkpoint: Start - {0}s'.format('%0.3f') % (time.time() - start_time))
    
    jsondata = jsoncheck.check(jsonFilePath, jobtype='simple')
    
    if jsondata==False:
        logger.error('Input json file not accepted, see check.json for details')
        return False
    
    logger.debug('Checkpoint: JsonCheck - {0}s'.format('%0.3f') % (time.time() - start_time))
    
    
    
    row, col = jsondata['grid_info']['steps_y'], jsondata['grid_info']['steps_x']
    Dtable = numpy.zeros((row, col))

    for item in jsondata['meshPositions']:
        i = item['indexY']
        j = item['indexZ']
        Dtable[j, i] = item['dozor_score']
    
    jsondata['MeshBest']['Dtable'] = Dtable
    numpy.savetxt('Dtable.txt', Dtable, fmt='%0.2f')
    plt.imshow(Dtable, cmap='hot', interpolation='nearest', origin='upper', extent=[0.5, (col + 0.5), (row + 0.5), 0.5])
    plt.colorbar()
    plt.savefig('Dtable.png', dpi=300, transparent=True, bbox_inches='tight')  # , pad_inches=0)
    plt.close()
    
    
    
    difminpar = jsondata['MeshBest']['difminpar']
    Ztable = numpy.zeros((row, col))
    Ztable[Dtable<difminpar] = -1
    jsondata['MeshBest']['Ztable'] = Ztable
    if numpy.all(Ztable==-1):
        logger.info('Diffraction signal is very weak for MeshBest')

        return jsondata
    
    logger.debug('Checkpoint: Initial data acquired - {0}s'.format('%0.3f') % (time.time() - start_time))
    
    dvanalysis.EliminateSaltRings(jsondata)
    logger.debug('Checkpoint: SaltRing Analysis - {0}s'.format('%0.3f') % (time.time() - start_time))

    dvanalysis.DetermineMCdiffraction(jsondata)
    logger.debug('Checkpoint: DV Analysis - {0}s'.format('%0.3f') % (time.time() - start_time))
    
    scoring.PerformCrystalRecognition(jsondata)
    logger.debug('Checkpoint: Crystal recognition - {0}s'.format('%0.3f') % (time.time() - start_time))

    jsondata['MeshBest']['Dtable'] = base64.b64encode(jsondata['MeshBest']['Dtable'])
    jsondata['MeshBest']['Ztable'] = base64.b64encode(jsondata['MeshBest']['Ztable'])
    jsondata['MeshBest']['positionReference'] = base64.b64encode(jsondata['MeshBest']['positionReference'])
    
    fig = plt.figure()
    ax = fig.add_subplot(111)
    
    plotting.MainPlot(jsondata, ax, addPositions=False)
    plt.savefig('CrystalMesh.png', dpi=150, transparent=True, bbox_inches='tight', pad_inches=0)
        
    with open('MeshResults.json', 'w') as outfile:
        json.dump(jsondata, outfile, sort_keys=True, indent=4, ensure_ascii=False)
    logger.debug('Checkpoint: Finish {0}s'.format('%0.3f') % (time.time() - start_time))
    
    return jsondata




def xraycentering(jsonFilePath, resultsPath=None):
    
    if resultsPath!=None:
        os.chdir(resultsPath)
    
    logger.debug('Checkpoint: Start - {0}s'.format('%0.3f') % (time.time() - start_time))
    
    jsondata = jsoncheck.check(jsonFilePath, jobtype='xraycentering')
    
    if jsondata==False:
        logger.error('Input json file not accepted, see check.json for details')
        return False
    
    logger.debug('Checkpoint: JsonCheck - {0}s'.format('%0.3f') % (time.time() - start_time))



    row, col = jsondata['grid_info']['steps_y'], jsondata['grid_info']['steps_x']
    Dtable = numpy.zeros((row, col))

    for item in jsondata['meshPositions']:
        i = item['indexY']
        j = item['indexZ']
        Dtable[j, i] = item['dozor_score']

    jsondata['MeshBest']['Dtable'] = Dtable
    numpy.savetxt('Dtable.txt', Dtable, fmt='%0.2f')
    plt.imshow(Dtable, cmap='hot', interpolation='nearest', origin='upper', extent=[0.5, (col + 0.5), (row + 0.5), 0.5])
    plt.colorbar()
    plt.savefig('Dtable.png', dpi=300, transparent=True, bbox_inches='tight')  # , pad_inches=0)
    plt.close()



    difminpar = jsondata['MeshBest']['difminpar']
    Ztable = numpy.zeros((row, col))
    Ztable[Dtable<difminpar] = -1
    jsondata['MeshBest']['Ztable'] = Ztable
    if numpy.all(Ztable==-1):
        logger.info('Diffraction signal is very weak for MeshBest')
        
        numpy.savetxt('Result_BestPositions.txt', [])

        jsondata['MeshBest']['Dtable'] = base64.b64encode(Dtable)
        jsondata['MeshBest']['Ztable'] = base64.b64encode(Ztable)
        
        jsondata['MeshBest']['BestPositions'] = base64.b64encode(numpy.ascontiguousarray(numpy.empty((0, 4), float)))
        with open('MeshResults.json', 'w') as outfile:
            json.dump(jsondata, outfile, sort_keys=True, indent=4, ensure_ascii=False)
            logger.debug('Checkpoint: Finish {0}s'.format('%0.3f') % (time.time() - start_time))
        return jsondata

    logger.debug('Checkpoint: Initial data acquired - {0}s'.format('%0.3f') % (time.time() - start_time))
    
    dvanalysis.EliminateSaltRings(jsondata)
    logger.debug('Checkpoint: SaltRing Analysis - {0}s'.format('%0.3f') % (time.time() - start_time))

    dvanalysis.DetermineMCdiffraction(jsondata)
    logger.debug('Checkpoint: DV Analysis - {0}s'.format('%0.3f') % (time.time() - start_time))

    scoring.PerformCrystalRecognition(jsondata)
    logger.debug('Checkpoint: Crystal recognition - {0}s'.format('%0.3f') % (time.time() - start_time))
    
    if numpy.all(jsondata['MeshBest']['Ztable'] < 0):
        logger.warning('Only multi-pattern diffraction found in the scanned area')
        
        sizecorr.GetAllPositions(jsondata)

    else:
        ellipticfit.DoEllipseFit(jsondata)
        
    logger.debug('Checkpoint: Shape Fit {0}s'.format('%0.3f') % (time.time() - start_time))

    jsondata['MeshBest']['Dtable'] = base64.b64encode(jsondata['MeshBest']['Dtable'])
    jsondata['MeshBest']['Ztable'] = base64.b64encode(jsondata['MeshBest']['Ztable'])
    jsondata['MeshBest']['positionReference'] = base64.b64encode(jsondata['MeshBest']['positionReference'])
    fig = plt.figure()
    ax = fig.add_subplot(111)
        
    plotting.MainPlot(jsondata, ax)
    plt.savefig('CrystalMesh.png', dpi=150, transparent=True, bbox_inches='tight')  # , pad_inches=0)
        
    with open('MeshResults.json', 'w') as outfile:
        json.dump(jsondata, outfile, sort_keys=True, indent=4, ensure_ascii=False)
    logger.debug('Checkpoint: Finish {0}s'.format('%0.3f') % (time.time() - start_time))
    
    return jsondata




def meshandcollect(jsonFilePath, resultsPath=None):
    
    if resultsPath!=None:
        os.chdir(resultsPath)
    
    logger.debug('Checkpoint: Start - {0}s'.format('%0.3f') % (time.time() - start_time))
    
    jsondata = jsoncheck.check(jsonFilePath, jobtype='meshandcollect')
    
    if jsondata==False:
        logger.error('Input json file not accepted, see check.json for details')
        return False
    
    logger.debug('Checkpoint: JsonCheck - {0}s'.format('%0.3f') % (time.time() - start_time))



    row, col = jsondata['grid_info']['steps_y'], jsondata['grid_info']['steps_x']
    Dtable = numpy.zeros((row, col))

    for item in jsondata['meshPositions']:
        i = item['indexY']
        j = item['indexZ']
        Dtable[j, i] = item['dozor_score']
    
    jsondata['MeshBest']['Dtable'] = Dtable
    numpy.savetxt('Dtable.txt', Dtable, fmt='%0.2f')
    plt.imshow(Dtable, cmap='hot', interpolation='nearest', origin='upper', extent=[0.5, (col + 0.5), (row + 0.5), 0.5])
    plt.colorbar()
    plt.savefig('Dtable.png', dpi=300, transparent=True, bbox_inches='tight')  # , pad_inches=0)
    plt.close()
    
    

    difminpar = jsondata['MeshBest']['difminpar']
    Ztable = numpy.zeros((row, col))
    Ztable[Dtable<difminpar] = -1
    jsondata['MeshBest']['Ztable'] = Ztable
    if numpy.all(Ztable==-1):
        logger.info('Diffraction signal is very weak for MeshBest')
        
        numpy.savetxt('Result_BestPositions.txt', [])
        
        jsondata['MeshBest']['Dtable'] = base64.b64encode(Dtable)
        jsondata['MeshBest']['Ztable'] = base64.b64encode(Ztable)
        
        jsondata['MeshBest']['BestPositions'] = base64.b64encode(numpy.ascontiguousarray(numpy.empty((0, 4), float)))
        
        return jsondata
    
    logger.debug('Checkpoint: Initial data acquired - {0}s'.format('%0.3f') % (time.time() - start_time))
    
    dvanalysis.EliminateSaltRings(jsondata)
    logger.debug('Checkpoint: SaltRing Analysis - {0}s'.format('%0.3f') % (time.time() - start_time))

    dvanalysis.DetermineMCdiffraction(jsondata)
    logger.debug('Checkpoint: DV Analysis - {0}s'.format('%0.3f') % (time.time() - start_time))

    scoring.PerformCrystalRecognition(jsondata)
    logger.debug('Checkpoint: Crystal recognition - {0}s'.format('%0.3f') % (time.time() - start_time))
    
    if numpy.all(jsondata['MeshBest']['Ztable'] < 0):
        logger.warning('Only multi-pattern diffraction found in the scanned area')

    sizecorr.GetAllPositions(jsondata)
    
    logger.debug('Checkpoint: Calculating positions {0}s'.format('%0.3f') % (time.time() - start_time))

    jsondata['MeshBest']['Dtable'] = base64.b64encode(jsondata['MeshBest']['Dtable'])
    jsondata['MeshBest']['Ztable'] = base64.b64encode(jsondata['MeshBest']['Ztable'])
    jsondata['MeshBest']['positionReference'] = base64.b64encode(jsondata['MeshBest']['positionReference'])
    fig = plt.figure()
    ax = fig.add_subplot(111)
    
    plotting.MainPlot(jsondata, ax)
    plt.savefig('CrystalMesh.png', dpi=150, transparent=True, bbox_inches='tight')  # , pad_inches=0)
        
    with open('MeshResults.json', 'w') as outfile:
        json.dump(jsondata, outfile, sort_keys=True, indent=4, ensure_ascii=False)
    logger.debug('Checkpoint: Finish {0}s'.format('%0.3f') % (time.time() - start_time))
    
    return jsondata




def linescan(jsonFilePath, resultsPath=None):
    
    if resultsPath!=None:
        os.chdir(resultsPath)
    
    logger.debug('Checkpoint: Start - {0}s'.format('%0.3f') % (time.time() - start_time))
    
    jsondata = jsoncheck.check(jsonFilePath, jobtype='linescan')
    
    if jsondata==False:
        logger.error('Input json file not accepted, see check.json for details')
        return False
    
    logger.debug('Checkpoint: JsonCheck - {0}s'.format('%0.3f') % (time.time() - start_time))
    
    
    
    row, col = jsondata['grid_info']['steps_y'], jsondata['grid_info']['steps_x']
    Dtable = numpy.zeros((row, col))
    for item in jsondata['meshPositions']:
        i = item['indexY']
        j = item['indexZ']
        Dtable[j, i] = item['dozor_score']

    jsondata['MeshBest']['Dtable'] = Dtable
    numpy.savetxt('Dtable.txt', Dtable, fmt='%0.2f')
    plt.imshow(Dtable, cmap='hot', interpolation='nearest', origin='upper',\
               extent=[0.5, (col + 0.5), (row + 0.5), 0.5])
    plt.colorbar()
    plt.savefig('Dtable.png', dpi=300, transparent=True, bbox_inches='tight')  # , pad_inches=0)
    plt.close()
    
    

    difminpar = jsondata['MeshBest']['difminpar']
    Ztable = numpy.zeros((row, col))
    Ztable[Dtable<difminpar] = -1
    jsondata['MeshBest']['Ztable'] = Ztable
    if numpy.all(Ztable==-1):
        logger.info('Diffraction signal is very weak for MeshBest line scan')
        
        numpy.savetxt('Result_BestPositions.txt', [])
        
        jsondata['MeshBest']['Dtable'] = base64.b64encode(Dtable)
        jsondata['MeshBest']['Ztable'] = base64.b64encode(Ztable)
        
        jsondata['MeshBest']['BestPositions'] = base64.b64encode(numpy.ascontiguousarray(numpy.empty((0, 4), float)))
        
        return jsondata

    logger.debug('Checkpoint: Initial data acquired - {0}s'.format('%0.3f') % (time.time() - start_time))
    
    dvanalysis.EliminateSaltRings(jsondata)
    logger.debug('Checkpoint: SaltRing Analysis - {0}s'.format('%0.3f') % (time.time() - start_time))

    dvanalysis.DetermineMCdiffraction(jsondata)
    logger.debug('Checkpoint: DV Analysis - {0}s'.format('%0.3f') % (time.time() - start_time))

    scoring.PerformCrystalRecognition(jsondata)
    logger.debug('Checkpoint: Crystal recognition - {0}s'.format('%0.3f') % (time.time() - start_time))
    
    if numpy.all(jsondata['MeshBest']['Ztable'] < 0):
        logger.info('Only multi-pattern diffraction found in the scanned area')
        
        C = 1 + numpy.sum(numpy.where(Ztable==-2)[0]*Dtable[Ztable==-2])/numpy.sum(Dtable[Ztable==-2])
        width = numpy.mean([numpy.size(Dtable[Dtable>level])\
                                    for level in numpy.linspace(0, 0.9*numpy.max(Dtable), 10)])
        
        BestPositions = numpy.array([[1.0, C, width, numpy.sum(Dtable)/100.0]])

    else:
        BestPositions = numpy.empty((0, 4), float)
        
        for v in numpy.unique(Ztable[Ztable>0]):
            C = 1 + numpy.sum(numpy.where(Ztable==v)[0]*Dtable[Ztable==v])/numpy.sum(Dtable[Ztable==v])
            
            eachArray = (Dtable*(Ztable==v))
            width = numpy.mean([numpy.size(eachArray[eachArray>level])\
                                    for level in numpy.linspace(0, 0.9*numpy.max(eachArray), 10)])
            BestPositions = numpy.append(BestPositions,\
                                         numpy.array([[1.0, C, width, numpy.sum(Dtable[Ztable==v])]]), axis=0)
        
        BestPositions = BestPositions[BestPositions[:, 3].argsort()][::-1]
        
    
    jsondata['MeshBest']['Dtable'] = base64.b64encode(jsondata['MeshBest']['Dtable'])
    jsondata['MeshBest']['Ztable'] = base64.b64encode(jsondata['MeshBest']['Ztable'])
    jsondata['MeshBest']['positionReference'] = base64.b64encode(jsondata['MeshBest']['positionReference'])
    jsondata['MeshBest']['BestPositions'] = base64.b64encode(numpy.ascontiguousarray(BestPositions))
    numpy.savetxt('Result_BestPositions.txt', BestPositions, fmt='%0.2f')
        
    
    fig = plt.figure()
    ax = fig.add_subplot(111)
    
    plotting.MainPlot(jsondata, ax, addPositions=False)
    plt.savefig('CrystalMesh.png', dpi=150, transparent=True, bbox_inches='tight', pad_inches=0)
    plt.clf()
    
    plotting.LinePlot(jsondata)
    plt.savefig('LineScan.png', dpi=150, transparent=True, bbox_inches='tight', pad_inches=0)
    plt.clf()
    

    with open('MeshResults.json', 'w') as outfile:
        json.dump(jsondata, outfile, sort_keys=True, indent=4, ensure_ascii=False)
    logger.debug('Checkpoint: Finish {0}s'.format('%0.3f') % (time.time() - start_time))
    
    return jsondata




#def hamburg(jsonFilePath, process_data=False, resultsPath=None):
#    
#    if resultsPath!=None:
#        os.chdir(resultsPath)
#
#    logger.debug('Checkpoint: Start - {0}s'.format('%0.3f') % (time.time() - start_time))
#    
#    jsondata = jsoncheck.check(jsonFilePath, jobtype='hamburg')
#    
#    if jsondata==False:
#        logger.error('Input json file not accepted, see check.json for details')
#        return False
#    
#    logger.debug('Checkpoint: JsonCheck - {0}s'.format('%0.3f') % (time.time() - start_time))
#    
#    
#    
#    row, col = jsondata['grid_info']['steps_y'], jsondata['grid_info']['steps_x']
#    Dtable = numpy.zeros((row, col))
#
#    for item in jsondata['meshPositions']:
#        i = item['indexY']
#        j = item['indexZ']
#        Dtable[j, i] = item['dozor_score']
#
#    jsondata['MeshBest']['Dtable'] = Dtable
#    numpy.savetxt('Dtable.txt', Dtable, fmt='%0.2f')
#    jsondata['MeshBest']['Dtable'] = base64.b64encode(Dtable)
#    plt.imshow(Dtable, cmap='hot', interpolation='nearest', origin='upper',\
#               extent=[0.5, (col + 0.5), (row + 0.5), 0.5])
#    plt.colorbar()
#    plt.savefig('Dtable.png', dpi=300, transparent=True, bbox_inches='tight')  # , pad_inches=0)
#    plt.close()
#    
#    difminpar = jsondata['MeshBest']['difminpar']
#    Ztable = numpy.zeros((row, col))
#    Ztable[Dtable<difminpar] = -1
#    jsondata['MeshBest']['Ztable'] = Ztable
#    if numpy.all(Ztable==-1):
#        logger.info('Diffraction signal is very weak for MeshBest')
#
#        jsondata['MeshBest']['Ztable'] = base64.b64encode(Ztable)
#        
#        return jsondata
#
#    logger.debug('Checkpoint: Initial data acquired - {0}s'.format('%0.3f') % (time.time() - start_time))
#    
#    dvanalysis.EliminateSaltRings(jsondata)
#    logger.debug('Checkpoint: SaltRing Analysis - {0}s'.format('%0.3f') % (time.time() - start_time))
#
#    dvanalysis.DetermineMCdiffraction(jsondata)
#    logger.debug('Checkpoint: DV Analysis - {0}s'.format('%0.3f') % (time.time() - start_time))
#    
#    scoring.PerformCrystalRecognition(jsondata)
#    logger.debug('Checkpoint: Crystal recognition - {0}s'.format('%0.3f') % (time.time() - start_time))
#    
#    Ztable = jsondata['MeshBest']['Ztable']
#    if numpy.all(jsondata['MeshBest']['Ztable'] < 0):
#        logger.info('Only multi-pattern diffraction found in the scanned area')
#        
#        print 'hello'
#    jsondata['MeshBest']['Ztable'] = base64.b64encode(jsondata['MeshBest']['Ztable'])
#
#    fig = plt.figure()
#    ax = fig.add_subplot(111)
#    plotting.MainPlot(jsondata, ax, addPositions=False)
#    plt.savefig('HambMesh.png', dpi=150, transparent=True, bbox_inches='tight')  # , pad_inches=0)
#    plt.close()
#
#    listOfEdges = []
#    Z = 0
#    first = -1
#    last = -1
#    for i in jsondata['meshPositions']:
#        if Ztable[i['indexZ'], i['indexY']]==Z:
#            last = i['index']
#        else:
#            if last!=-1:
#                if last!=first:
#                    listOfEdges.append((first, last, Z))
#                first = -1
#                last = -1
#            Z = 0
#            if Ztable[i['indexZ'], i['indexY']]>0:
#                first = i['index']
#                last = i['index']
#                Z = Ztable[i['indexZ'], i['indexY']]
#
#
#    jsondata['MeshBest']['listOfEdges'] = listOfEdges
#    jsondata['MeshBest']['positionReference'] = base64.b64encode(jsondata['MeshBest']['positionReference'])
#    
#    with open('MeshResults.json', 'w') as outfile:
#        json.dump(jsondata, outfile, sort_keys=True, indent=4, ensure_ascii=False)
#        
#    if process_data:
#        
#    
#    
#    
#    
#    
#    
#    
#
#    logger.debug('Checkpoint: Finish {0}s'.format('%0.3f') % (time.time() - start_time))
#    
#    return listOfEdges




def threedcentering(jsonFilePath1, jsonFilePath2, jsonFilePath3=None, jsonFilePath4=None, resultsPath=None):
    
    #all mesh scan have to be already processed with 'simple' algorithm
    
    if resultsPath!=None:
        os.chdir(resultsPath)
    else:
        resultsPath = os.getcwd()

    f1 = open(jsonFilePath1, 'r')
    jsondata1 = json.load(f1)
    f1.close()
    
    Dtable1 = numpy.fromstring(base64.b64decode(jsondata1['MeshBest']['Dtable']))
    Dtable1 = numpy.reshape(Dtable1, (jsondata1['grid_info']['steps_y'], jsondata1['grid_info']['steps_x']))
    if numpy.all(Dtable1<=0.3):
        return None
    
    fi1 = jsondata1['meshPositions'][0]['omega']
    
    f2 = open(jsonFilePath2, 'r')
    jsondata2 = json.load(f2)
    f2.close()
    
    Dtable2 = numpy.fromstring(base64.b64decode(jsondata2['MeshBest']['Dtable']))
    Dtable2 = numpy.reshape(Dtable2, (jsondata2['grid_info']['steps_y'], jsondata2['grid_info']['steps_x']))
    if numpy.all(Dtable2<=0.3):
        return None
    
    fi2 = jsondata2['meshPositions'][0]['omega']

    if jsonFilePath3!=None and jsonFilePath4!=None:
        f3 = open(jsonFilePath3, 'r')
        jsondata3 = json.load(f3)
        f3.close()
        
        Dtable3 = numpy.fromstring(base64.b64decode(jsondata3['MeshBest']['Dtable']))
        Dtable3 = numpy.reshape(Dtable3, (jsondata3['grid_info']['steps_y'], jsondata3['grid_info']['steps_x']))
        if numpy.all(Dtable3<=0.3):
            return None
        
        fi3 = jsondata3['meshPositions'][0]['omega']

        f4 = open(jsonFilePath4, 'r')
        jsondata4 = json.load(f4)
        f4.close()
        
        Dtable4 = numpy.fromstring(base64.b64decode(jsondata4['MeshBest']['Dtable']))
        Dtable4 = numpy.reshape(Dtable4, (jsondata4['grid_info']['steps_y'], jsondata4['grid_info']['steps_x']))
        if numpy.all(Dtable4<=0.3):
            return None
        
        fi4 = jsondata4['meshPositions'][0]['omega']
        
        fi_list = [fi1, fi2, fi3, fi4]
        
        Cmass = voxelise.CenterOfMass(Dtable1, Dtable2, Dtable3, Dtable4, fi_list)
        return Cmass
    
    else:
        C1 = ndimage.measurements.center_of_mass(Dtable1)
        C2 = ndimage.measurements.center_of_mass(Dtable2)
        return C1, C2

    

    
    
    
    
    
    
    
    
    
    
    
    
    
