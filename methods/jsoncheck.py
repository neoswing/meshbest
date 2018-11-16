'''
Created on Sep 11, 2018

@author: melnikov
'''

import numpy
import json
import os



try:
    from workflow_lib import workflow_logging
    logger = workflow_logging.getLogger()
except:
    import logging
    logger = logging.getLogger("MeshBest")





def check(jsonFilePath, jobtype='simple'):
    
    check = True
    runreport = {}
    difminpar = 0.2
    jsondata = {}
    
    if os.path.isfile(jsonFilePath):
        json_file = open(jsonFilePath, 'r')
        jsondata = json.load(json_file)
        json_file.close()
        runreport['File'] = jsonFilePath
    else:
        runreport['File'] = 'NOT OK'
        logger.error('Json File not found')
        check = False
    
    try:
        runreport['Ncolumns'] = jsondata['grid_info']['steps_x']
        runreport['Nrows'] = jsondata['grid_info']['steps_y']
    except KeyError:
        runreport['Ncolumns'] = 'NOT OK'
        runreport['Nrows'] = 'NOT OK'
        logger.error('Number of columns/rows not specified')
        check = False
    
    try:
        runreport['DetectorDistance'] = jsondata['inputDozor']['detectorDistance'] * 1000
        runreport['DetectorPixelSize'] = jsondata['beamlineInfo']['detectorPixelSize']*1000
    except KeyError:
        runreport['DetectorDistance'] = 'NOT OK'
        runreport['DetectorPixelSize'] = 'NOT OK'
        logger.error('Detector parameters not specified (detector distance, detector pixel size)')
        check = False
    
    try:
        runreport['Wavelength'] = jsondata['inputDozor']['wavelength']
        runreport['BeamCentreX'] = jsondata['inputDozor']['orgx']
        runreport['BeamCentreY'] = jsondata['inputDozor']['orgy']
        if jobtype!='linescan':
            runreport['BeamSizeX'] = int(jsondata['grid_info']['beam_width']*1000)
            runreport['BeamSizeY'] = int(jsondata['grid_info']['beam_height']*1000)
    except KeyError:
        runreport['Wavelength'] = 'NOT OK'
        runreport['BeamCentreX'] = 'NOT OK'
        runreport['BeamCentreY'] = 'NOT OK'
        runreport['BeamSizeX'] = 'NOT OK'
        runreport['BeamSizeY'] = 'NOT OK'
        logger.error('Beam parameters not specified (wavelength,\
        beam centre X and Y, beam size X and Y)')
        check = False
    
    try:
        runreport['Total N images'] = len(jsondata['meshPositions'])
    except KeyError:
        runreport['Total N images'] = 'NOT OK'
        logger.error('Json lacks meshPositions section')
        check = False
    
    try:
        for i in jsondata['meshPositions']:
            temp = i['index']
            temp = i['indexY']
            temp = i['indexZ']
            temp = i['omega']
            temp = i['dozor_score']
            if temp>difminpar:
                try:
                    temp2 = i['dozorSpotList']
#                    temp2 = i['dozorSpotListShape']
                except KeyError:
                    logger.error('Spot lists not found despite diffraction signal')
                    check = False
        runreport['meshPositions'] = 'OK'
    except KeyError:
        runreport['meshPositions'] = 'NOT OK'
        logger.error('meshPositions section not well structured (lacks index, Y, Z, omega)')
        check = False
    
    
    
    if jobtype=='meshandcollect' and check:
        try:
            runreport['AvailableApertures'] = jsondata['beamlineInfo']['beamlineApertures']
        except KeyError:
            logger.warning('List of available apertures not specified. Using the default.')
            jsondata['beamlineInfo']['beamlineApertures'] = [str(jsondata['grid_info']['beam_width']*1000)]
            runreport['AvailableApertures'] = 'default '+str(jsondata['grid_info']['beam_width']*1000)
    
    if jobtype=='linescan' and check:
        if runreport['Ncolumns'] != 1:
            logger.error('MeshBest linescan method is used not in line scan')
            check = False
    
    if jobtype=='hamburg' and check:
        if round(jsondata['meshPositions'][1]['omega']-jsondata['meshPositions'][0]['omega'], 2)==0:
            logger.error('Omega is not changed during the mesh scan')
            check = False
    
    
    
    
    with open('check.json', 'w') as outfile:
        json.dump(runreport, outfile, sort_keys=True, indent=4, ensure_ascii=False)
    
    if check:
        positionReference = numpy.empty((runreport['Nrows'], runreport['Ncolumns']), dtype='int')
        for i in jsondata['meshPositions']:
            positionReference[i['indexZ'], i['indexY']] = i['index']
            
        jsondata['MeshBest'] = {'type': jobtype, 'check': True,\
                                'difminpar': difminpar, 'positionReference': positionReference}
        return jsondata
    else:
        return False
    
    