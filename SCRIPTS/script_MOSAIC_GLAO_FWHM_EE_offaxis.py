#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu May  6 10:11:46 2021

@author: omartin
"""

#%% Import librairies
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
from astropy.io import fits
from configparser import ConfigParser

from aoSystem.fourierModel import fourierModel
import maolyses_utils as utils

#%% PATHS
path_maolyses = '/home/omartin/Projects/MOSAIC/AOsimulations/MAOLYSES'
path_p3       = '/home/omartin/Projects/P3'
path_ini      = path_maolyses + '/INI/MosaicGLAOParams_grid_multiwvl_dev.ini'
pathCn2       = path_maolyses + '/DATA/'
fileCn2       = 'profil_turbulent_eso.fits'
path_save     = '/home/omartin/Projects/MOSAIC/AOsimulations/RESULTS/GLAO_offaxis/'
Cn2_eso       = fits.getdata(pathCn2+fileCn2)

parser = ConfigParser()
parser.optionxform = str
parser.read(path_ini)

#%% RUN THE SIMULATION

def GetContours(path_ini,EEradius,fovDiameterInArcsec,nSamples,wvl=[500e-9,640e-9,890e-9,1.2e-6,1.65e-6,2.12e-6],profile='Median'):
    
    # SELECT THE PROFILE
    nameProf   = np.array(['JQ1','JQ2','JQ3','JQ4','Median'])
    if profile not in nameProf:
        raise ValueError('The profile name is not recognized; should be JQ1, JQ2, JQ3, JQ4 or Median')
    else:
        idP  = np.argwhere(np.array(nameProf)==profile)[0][0]
        name = nameProf[idP]
        zenith_angle, seeing, Cn2Weights, Cn2Heights, WindSpeed, WindDir \
        = utils.FromProfileToAtmosphereClass(Cn2_eso,name)
    
    # DEFINING THE SAMPLING OF THE FOV
    x = np.linspace(-fovDiameterInArcsec//2,fovDiameterInArcsec//2,nSamples)
    X,Y = np.meshgrid(x,x)
    zS  = list(np.hypot(X,Y).reshape(-1))
    aS  = list(180/np.pi * np.arctan2(Y,X).reshape(-1))
    
    # UPDATE THE .INI FILE
    parser.set('sources_science','Wavelength',str(wvl))
    parser.set('sources_science','Zenith',str(zS))
    parser.set('sources_science','Azimuth',str(aS))
    parser.set('telescope','ZenithAngle',str(zenith_angle))
    parser.set('atmosphere','Cn2Weights',str(Cn2Weights))
    parser.set('atmosphere','Cn2Heights',str(Cn2Heights))
    parser.set('atmosphere','WindSpeed',str(WindSpeed))
    parser.set('atmosphere','WindDirection',str(WindDir))
    parser.set('atmosphere','Seeing',str(seeing))  
    with open(path_ini, 'w') as configfile:
        parser.write(configfile)
            
    # SIMULATE THE PSF AND MEASURE THE METRICS
    fao = fourierModel(path_ini, path_root = path_p3,\
                   nyquistSampling=True,display=False,verbose=True,\
                   getErrorBreakDown=False,getFWHM=True,getEncircledEnergy=True)

    # DISPLAY CONTOURS
    for n in range(fao.freq.nWvl):
        fao.displayPsfMetricsContours(EEradius,wvlIndex=n)
        
    return fao    


#%%
fao = GetContours(path_ini,200,420,7,wvl=[1.65e-6],profile='Median')