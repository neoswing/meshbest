#!/usr/bin/env python

# v9
#many changes
#newDistanceFunction
#new method for overlap treating
#v9.2
#polished structure
#calculation now takes place only in each diffraction zone
#v9.3
#the distance between the images is no more counted in detector-dependent pixels but angular difference
#v9.4
#overlap detection was linked with number of satellite spots instead of histogram slope
#v9.5
#fixed bug with diagonal ellipse fit to the mesh scan squares
#v10
#Changed the way ellipse fit is organised: only for ~single crystals it is used; for multiple crystals
#    in the sample area we try to implement zone-size correlation to adapt the appropriate aperture size
#Input parameters are assembled together
#Modified linkage for clustering
#v11
#moved to python 3 syntax
#code optimisation to speed up calculations


'''
###################################################################################################################
###################################################################################################################
###################################################################################################################
###################################################################################################################
###################################################################################################################
##############.___  ___.  _______     _______. __    __     .______    _______     _______.___________.############
##############|   \/   | |   ____|   /       ||  |  |  |    |   _  \  |   ____|   /       |           |############
##############|  \  /  | |  |__     |   (----`|  |__|  |    |  |_)  | |  |__     |   (----`---|  |----`############
##############|  |\/|  | |   __|     \   \    |   __   |    |   _  <  |   __|     \   \       |  |     ############
##############|  |  |  | |  |____.----)   |   |  |  |  |    |  |_)  | |  |____.----)   |      |  |     ############
##############|__|  |__| |_______|_______/    |__|  |__|    |______/  |_______|_______/       |__|     ############
###################################################################################################################
#################################################            ######################################################
#############################################      METHODS      ###################################################
#################################################            ######################################################
###################################################################################################################
'''
from meshbest.methods import dvanalysis, scoring, plotting, ellipticfit, sizecorr

