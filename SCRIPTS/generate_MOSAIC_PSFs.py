#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue May  4 14:45:01 2021

@author: omartin
"""

#%% IMPORTING LIBRARIES
import numpy as np
from aoSystem.fourierModel import fourierModel


#%% MANAGING PATHS
path_maolyses = '/home/omartin/Projects/MOSAIC/AOsimulations/MAOLYSES'
path_p3       = '/home/omartin/Projects/P3'

#%% ON-AXIS CASE - SINGLE WAVELENGTH

path_ini = path_maolyses + '/INI/MosaicGLAOParams_onaxis_singlewvl.ini'
fao = fourierModel(path_ini, path_root = path_p3,\
                   nyquistSampling=True,display=False,verbose=True,\
                   getErrorBreakDown=False,getFWHM=True,getEncircledEnergy=True)

FWHM = np.mean(fao.FWHM[:,0,0],axis=0)
EE  = fao.EncE[:,0,0]
#%% ON-AXIS CASE - MULTI WAVELENGTH

path_ini = path_maolyses + '/INI/MosaicGLAOParams_onaxis_multiwvl.ini'
fao = fourierModel(path_ini, path_root = path_p3,\
                   nyquistSampling=True,display=False,verbose=True,\
                   getErrorBreakDown=False,getFWHM=True,getEncircledEnergy=True)

FWHM = np.squeeze(fao.FWHM).mean(axis=0)
EE  = np.squeeze(fao.EncE[:,:,0])

#%% OFF-AXIS CASE - SINGLE WAVELENGTH

path_ini = path_maolyses + '/INI/MosaicGLAOParams_grid_singlewvl.ini'
fao = fourierModel(path_ini, path_root = path_p3,\
                   nyquistSampling=True,display=False,verbose=True,\
                   getErrorBreakDown=False,getFWHM=True,getEncircledEnergy=True)

FWHM = np.squeeze(fao.FWHM).mean(axis=0)
EE  = np.squeeze(fao.EncE[:,:,0])

#%% OFF-AXIS CASE - MULTI WAVELENGTH

path_ini = path_maolyses + '/INI/MosaicGLAOParams_grid_multiwvl.ini'
fao = fourierModel(path_ini, path_root = path_p3,\
                   nyquistSampling=True,display=False,verbose=True,\
                   getErrorBreakDown=False,getFWHM=True,getEncircledEnergy=True)

FWHM = fao.FWHM.mean(axis=0)
EE  = np.squeeze(fao.EncE)