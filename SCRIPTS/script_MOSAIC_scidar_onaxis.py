#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jun  2 08:36:24 2021

@author: omartin
"""



# !!!!! IMPORTANT NOTES !!!!!!

# The deault L0 is 25 m
# PSFs/FWHM/EE400 mas are derived on-axis at 1.65µm
# The PSf is Nyquist-sampled (can be changed in specifying the pixel scale in the .ini file and setting nyquistSampling to False when instantiating fourierModel)
# The Cn2 profiles are the new Scidar profiles, assumed to be at zenith
# The default zenith angle is 30 deg (can be changed in the .ini file)
# No additionnal residual jitter (can be changed in the .ini file)
# The wind speed.wind directions are obtained from ESO 35-layers profiles. The wind speed is then mean-rescaled to the average wind speed from the Scidar 
# 4 cases: NOAO, GLAO, MOAO with 5'-diameter NGS constellation of 3 stars MOAO 2.5'- 3 NGS
#%% Import librairies
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
from astropy.io import fits
from configparser import ConfigParser
import time
from aoSystem.fourierModel import fourierModel
import aoSystem.FourierUtils as FourierUtils

#%% PATHS
path_ini    = '/home/omartin/Projects/MOSAIC/AOsimulations/MAOLYSES/INI/Mosaic_scidar.ini'
path_scicar = '/home/omartin/Projects/MOSAIC/AOsimulations/MAOLYSES/DATA/stereoSCIDAR2021/'
path_eso    = '/home/omartin/Projects/MOSAIC/AOsimulations/MAOLYSES/DATA/'
path_save   = '/home/omartin/Projects/MOSAIC/AOsimulations/RESULTS/GLAO_versus_MOAO/'

parser = ConfigParser()
parser.optionxform = str
parser.read(path_ini)

#%% SCIDAR PROFILES

#    • profilesTime.fits → Array[19356] Julian Date of profiles 
#    • profilesCn2.fits → Array[100, 19356] Cn2 profiles 
#    • profilesWind.fits → Array[100, 19356] wind profiles 
#    • profilesParams.fits → Array[19356, 3] seeing, theta0, and average wind speed (this is computed considering only layers where wind speed is >0) 
#    • referenceProfiles35layers.fits → Array[35, 6] altitude vector (35 layers) + 5 reference profiles
#    • referenceProfiles100layers.fits → Array[100, 6] altitude vector (100 layers) + 5 reference profiles
    
Ceso   = fits.getdata(path_eso + 'profil_turbulent_eso.fits')
Cn2    = fits.getdata(path_scicar + 'profilesCn2.fits')
wind   = fits.getdata(path_scicar + 'profilesWind.fits')
param  = fits.getdata(path_scicar + 'profilesParams.fits')
seeing = param[0]
wMean  = param[2]
wSpeed = Ceso[:,2]
wSpeed = wSpeed/wSpeed.mean()
wDir   = Ceso[:,8]*180/np.pi

nProf  = Cn2.shape[0]
nL     = Cn2.shape[1]

#%% SYSTEM CONFIG
savePSF = False

# atmosphere
L0 = 25
nEqLayers = len(wDir) # the profiles will be compressed over 35 layers
scidar_heights = np.linspace(0,24750,num=100)
# note: the wind direction is hardcoded

# Ensquared Energy
wvl     = 1.65e-6
EEbox   = 400 # diameter of the box in mas
psInMas = 3600 * 180 * 1000/np.pi * wvl/39/2
nntrue  = EEbox/psInMas/2
nn2     = int(nntrue)

# config
aoConfig_all = ['NOAO' , 'GLAO' , 'MOAO_poor' , 'MOAO_good']
nCases = len(aoConfig_all)
# instantiating outputs
nPup = 128
nPsf = 2 * nPup
parser.set('telescope','Resolution',str(nPup))
parser.set('sensor_science','FiedOfView',str(nPsf))
with open(path_ini, 'w') as configfile:
            parser.write(configfile)
            
if savePSF:
    psfSimu = np.zeros((nCases,nProf,nPsf,nPsf))  
FWHM_all    = np.zeros((nCases,nProf))
EE400_all   = np.zeros((nCases,nProf))

#%% LOOPS ON AO CONFIGURATION - CN2 CONDITIONS - 
t0 = time.time()
for k in range(nCases):
    # UPDATING THE CONFIGURATION
    aoConfig = aoConfig_all[k]   
    
    ## NOAO Case
    if aoConfig == 'NOAO':
        loopGain = 0
    else:
        loopGain = 0.5
        
        
    ## LGS ASTERISM 
    ###############
    LGSFoV              = 7.4*60 #2.4*60 #arcsec
    nLGS                = 4
    GuideStarZenith     = [] 
    GuideStarAz         = [] 
    for i in range(nLGS):
        GuideStarZenith.append(LGSFoV/2)
        GuideStarAz.append(i*360/nLGS)
        
        
    ## MOAO case ##
    ###############
    if aoConfig == 'MOAO_poor': # THE NGS ARE LOCATED AT 5' FROM ON-AXIS
        #add of  NGS 
        nbNGS= 4
        fovNGS = 5 # in arcmin
        if nbNGS !=0 :
            for istar in range(nbNGS) :
                GuideStarZenith.append(fovNGS*60/2)
                GuideStarAz.append(istar*360/nbNGS)
    
    if aoConfig == 'MOAO_good': # THE NGS ARE LOCATED AT 1' FROM ON-AXIS
        #add of  NGS 
        nbNGS= 4
        fovNGS = 2.5 # in arcmin
        if nbNGS !=0 :
            for istar in range(nbNGS) :
                GuideStarZenith.append(fovNGS*60/2)
                GuideStarAz.append(istar*360/nbNGS)
                
    ## OPTIMISATION DIRECTION 
    #########################
    optimRad     = []
    optimAz      = []
    optim        = []
    
    if aoConfig == 'MOAO_poor' or aoConfig == 'MOAO_good':  # WE OPTIMIZE ON-AXIS
        optimRad.append([0])
        optimAz.append([0])
        optim.append([1])
        NumberLenslets = [64,64,64,64,64,64,64,64]
        SizeLenslets   = [0.6,0.6,0.6,0.6,0.6,0.6,0.6,0.6]
        NumberPhotons  = [500,500,500,500,500,500,500,500] 
        nRec           = 9
    else : # WE OPTIMIZE in THE FOV
        NumberLenslets = [64,64,64,64]
        SizeLenslets   = [0.6,0.6,0.6,0.6]
        NumberPhotons  = [500,500,500,500] 
        nRec           = 1
        
        OptimFoV     = 5*60 #arcsec
        nopt         = 5 #nb point lineaire pour l'optim (au total nopt²) 
        optimX       = np.linspace(-OptimFoV/2,OptimFoV/2,num=nopt)#+1e-6
        optimY       = np.linspace(-OptimFoV/2,OptimFoV/2,num=nopt)#+1e-6
      
        for i in range(nopt):
            for j in range(nopt): 
                a  = np.sqrt( optimX[i] * optimX[i] + optimY[j] * optimY[j] )
                if optimY[j] == 0 : 
                    if optimX[i] < 0 : 
                        b = -90
                    else : 
                        b = 90
                else :  
                    b  = 180*np.arctan2(optimX[i],optimY[j])/np.pi
                optimRad.append(a)
                optimAz.append(b)
                optim.append(1)
                    
    # UPDATE .INI FILE
    parser.set('sources_HO','Zenith',str(GuideStarZenith))
    parser.set('sources_HO','Azimuth',str(GuideStarAz))
    parser.set('RTC','LoopGain_HO',str(loopGain))
    parser.set('DM','OptimizationZenith',str(optimRad))
    parser.set('DM','OptimizationAzimuth',str(optimAz))
    parser.set('DM','OptimizationWeight',str(optim))
    parser.set('DM','NumberReconstructedLayers',str(nRec))
    parser.set('sensor_HO','NumberPhotons',str(NumberPhotons))
    parser.set('sensor_HO','SizeLenslets',str(SizeLenslets))
    parser.set('sensor_HO','NumberLenslets',str(NumberLenslets))
    parser.set('atmosphere','L0',str(L0))
    
    
    for i  in range(nProf):
             
        # HANDLE THE PROFILES
        
        # select the ith profile
        Cn2_tmp = Cn2[i,:]
        wind_tmp = wind[i,:]
        # remove NAN values
        Cn2_tmp[Cn2_tmp!=Cn2_tmp] = 0
        idG = Cn2_tmp > 1e-20
        # compress the profile to speed up the simulation
        Cn2Weights, Cn2Heights  = FourierUtils.eqLayers(Cn2_tmp[idG]/Cn2_tmp[idG].sum(),scidar_heights[idG],nEqLayers)
        
        # UPDATE .INI FILE
        parser.set('telescope','ZenithAngle',str(0))
        parser.set('atmosphere','Cn2Weights',str(list(Cn2Weights)))
        parser.set('atmosphere','Cn2Heights',str(list(Cn2Heights)))
        parser.set('atmosphere','WindSpeed',str(list(wSpeed * wMean[i] )))
        parser.set('atmosphere','WindDirection',str(list(wDir)))
        parser.set('atmosphere','Seeing',str(seeing[i]))  
        with open(path_ini, 'w') as configfile:
            parser.write(configfile)
        
        # PSF SIMULATIONS
        fao = fourierModel(path_ini,calcPSF=True,verbose=False,display=False,\
                           getErrorBreakDown=False,getFWHM=True,getEncircledEnergy=True,\
                           getEnsquaredEnergy=False, displayContour=False,nyquistSampling=True) 
        
        # CONCATENATING PSFS
        if savePSF:
            psfSimu[k,i] = np.squeeze(np.transpose(fao.PSF,axes=(2,3,0,1)))
        
        # CONCATENATING FWHM
        FWHM_all[k,i] = np.squeeze(np.mean(fao.FWHM,axis=0))
        
        # CONCATENATING EE
        EE400_all[k,i] = (nntrue - nn2)*fao.EncE[nn2+1] + (nn2+1-nntrue)*fao.EncE[nn2]

tsimu = time.time() - t0
#%% SAVING RESULTS
#if path_save != '':
#    if savePSF:
#        # PSF simu 
#        hdu = fits.PrimaryHDU(psfSimu)
#        hdul= fits.HDUList(hdu)
#        hdul.writeto(path_save+'PSFsimu_on-axis_NOAO_GLAO_MOAOpoor_MOAOgood_scidar.fits',overwrite=True)
#    
#    # FWHM/EE 
#    hdu1 = fits.PrimaryHDU(wvl)
#    hdu2 = fits.ImageHDU(FWHM_all)
#    hdu3 = fits.ImageHDU(EE400_all)
#    hdul= fits.HDUList([hdu1, hdu2, hdu3])
#    hdul.writeto(path_save+'wvl_FWHM_EE400mas_on-axis_NOAO_GLAO_MOAOpoor_MOAOgood_scidar.fits',overwrite=True)