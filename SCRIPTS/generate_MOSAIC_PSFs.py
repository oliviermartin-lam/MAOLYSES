#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue May  4 14:45:01 2021

@author: omartin
"""

#%% IMPORTING LIBRARIES
from aoSystem.fourierModel import fourierModel


#%% MANAGING PATHS
path_maolyses = '/home/omartin/Projects/MOSAIC/AOsimulations/MAOLYSES'
path_p3       = '/home/omartin/Projects/P3'
path_ini      = path_maolyses + '/INI/MosaicGLAOParams_onaxis_singlewvl.ini'

#%% INSTANTIATING THE FOURIER MODEL

fao = fourierModel(path_ini, path_root = path_p3,\
                   nyquistSampling=True,display=False,verbose=True,\
                   getErrorBreakDown=False,getFWHM=True,getEncircledEnergy=True)