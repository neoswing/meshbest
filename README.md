# meshbest



    This is Mesh Best software module.

##        -------Tasks-------

MeshBest performs recognition of the crystal samples based on a mesh scan (or 2D X-ray diffraction scan routine).
It says how many crystals are in the scanned 2D area, how well they diffract and also reports their dispositions
and size approximations.

##        -------Input-------

*MeshResults*.json dictionary file containing:

* Numbers of columns and rows of the mesh scan (col, row)

* X-ray wavelength used, Detector distance, Centre coordinates of the beam (in detector px), Beam dimensions (mm),
detector pixel size (mm) and a list of available aperture sizes (in microns) to predict aperture choice

* Lists of detected spot coordinates and for each image of the mesh scan: in detector coordinates, in base64 string
format

* Overall diffraction score for each image

MeshBest uses .json file containing all experiment parameters (of a mesh scan) and output given by Dozor pre-
analysis. The output from Dozor should contain diffraction score evaluation for every image of the mesh scan and
as well a list of detected diffraction spot coordinates in base64 string format.

##        -------Output-------

MeshBest produces a 2D colour map indicating different crystal zones found in the sample area. This map appears
under the name **CrystalMesh.png** in the MeshBest working directory. Most of MeshBest output can be found in
**MeshResults.json file**, in **['MeshBest']** partition of the dictionary. The dictionary is also returned by
**simple()** method function.


### Best Positions

For multicrystal MeshAndCollect experiment pipeline MeshBest does pre-analysis where estimates positions and
appropriate beam sizes (depends on the available aperture sizes on the beamline) for collecting small wedges
of data.The output array is returned in **['MeshBest']['BestPositions']** dictionary partition in base64 string
format. The structure of the output array is the following:

*    column 0: X coordinate of the centre
*    column 1: Y coordinate of the centre
*    column 2: Suggested aperture size
*    column 3: Integral diffraction quality

Every row represents a particular position/crystal for data collection. "Result_BestPositions.txt" is the text
file with an output array.


### Elliptic Fit

For X-ray Centering pipeline MeshBest does crystal approximations on the map to elliptic shapes.

If elliptic approximations have been made to crystal areas, the output dictionary contains ellipse parameters in
**['MeshBest']['EllipseArray']** in base64 string format.

The output also gives a text file with parameters of ellipses - "Result_Ellipses.txt", where the architecture
is similar:

*    column 0: X coordinate of the centre
*    column 1: Y coordinate of the centre
*    column 2: Ellipse long axis
*    column 3: Ellipse short axis
*    column 4: Angle between the long axis and mesh scan X-axis
*    column 5: Integral diffraction quality

##        -------Usage-------

To proceed with the simple algorithm where returned is only a crystal map one should call the function **simple()** imported
from meshbest.algorithms. One should pass the json file to the function with all experiment parameters and
pre-analysis by Dozor. If working directory is not specified as a second argument to **simple()**, MeshBest
will work in the CWD and produce related output files there.
