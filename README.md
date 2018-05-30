# meshbest



    This is Mesh Best software module.

#        -------Tasks-------

    MeshBest performs recognition of the crystal samples based on a mesh scan (or 2D X-ray diffraction scan routine).
    It says how many crystals are in the scanned 2D area, how well they diffract and also reports their dispositions
    and size approximations.

#        -------Input-------

    *MeshResults*.json dictionary file containing:

    -Numbers of columns and rows of the mesh scan (col, row)

    -X-ray wavelength used, Detector distance, Centre coordinates of the beam (in detector px), Beam dimensions (mm),
    detector pixel size (microns) and a list of available aperture sizes (in microns) to predict aperture choice

    -Lists of detected spot coordinates and for each image of the mesh scan: in detector coordinates, in base64 string
    format

    -Overall diffraction score for each image

    MeshBest uses .json file containing all experiment parameters (of a mesh scan) and output given by Dozor pre-
    analysis. The output from Dozor should contain diffraction score evaluation for every image of the mesh scan and
    as well a list of detected diffraction spot coordinates in base64 string format.
    
#        -------Output-------
    
    MeshBest produces a 2D colour map indicating different crystal zones found in the sample area.
    
    If MeshBest detected less than 3 crystals in the sample area the output given is elliptic approximation of the
    crystals. If more crystals are detected then the output gives best positions and corresponding aperture choices
    for multicrystal data collection. The structure of the output array is independent of the case and is the
    following:
    
    column 0: X coordinate of the centre
    column 1: Y coordinate of the centre
    column 2: Suggested aperture size
    column 3: Integral diffraction quality
    
    Every row represents a particular position/crystal for data collection. "Result_BestPositions.txt" is the text
    file with an output array.
    
    The output also gives a text file with parameters of ellipses - "Result_Ellipses.txt", where the architecture
    is similar:
    
    column 0: X coordinate of the centre
    column 1: Y coordinate of the centre
    column 2: Ellipse long axis
    column 3: Ellipse short axis
    column 4: Angle between the long axis and mesh scan X-axis
    column 5: Integral diffraction quality

#        -------Usage-------
    
    To proceed with the classic algorithm of mesh scan analysis one should call the function classic() imported
    from meshbest.algorithms. One should pass the json file to the function with all experiment parameters and
    pre-analysis by Dozor. If working directory is not specified as a second argument to classic(), MeshBest
    will work in the CWD and produce related output files there.
    









#Changes_log

 v9
many changes
newDistanceFunction
new method for overlap treating
v9.2
polished structure
calculation now takes place only in each diffraction zone
v9.3
the distance between the images is no more counted in detector-dependent pixels but angular difference
v9.4
overlap detection was linked with number of satellite spots instead of histogram slope
v9.5
fixed bug with diagonal ellipse fit to the mesh scan squares
v10
Changed the way ellipse fit is organised: only for ~single crystals it is used; for multiple crystals in the sample
area we try to implement zone-size correlation to adapt the appropriate aperture size
Input parameters are assembled together
Modified linkage for clustering
