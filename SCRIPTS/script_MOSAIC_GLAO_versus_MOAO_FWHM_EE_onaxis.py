#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed May  5 20:43:52 2021

@author: omartin
"""

#%% Import librairies
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
from astropy.io import fits
from configparser import ConfigParser

from aoSystem.fourierModel import fourierModel

#%% PATHS
parfile = '/home/omartin/Projects/MOSAIC/AOsimulations/MAOLYSES/INI/Mosaic_scidar.ini'
pathCn2 = '/home/omartin/Projects/MOSAIC/AOsimulations/MAOLYSES/DATA/'
fileCn2 = 'profil_turbulent_eso.fits'
path_save = '/home/omartin/Projects/MOSAIC/AOsimulations/RESULTS/GLAO_versus_MOAO/'

parser = ConfigParser()
parser.optionxform = str
parser.read(parfile)

#%% OBSERVING CONDITIONS
C          = fits.open(pathCn2+fileCn2)
Cn2        = C[0].data
h          = Cn2[:,1]
v          = Cn2[:,2]
vDir       = Cn2[:,8]*180/np.pi
L0         = 25
nameProf   = ['JQ1','JQ2','JQ3','JQ4','Median']
ZenithProf = np.array([30, 30, 30, 45, 30])
wvl        = np.array([500e-09, 640e-09, 890e-9, 1.2e-06, 1.65e-06, 2.17e-06])
nWvl       = len(wvl)
nProf      = len(nameProf)

# Ensquared Energy
EEbox   = 400 # diameter of the box in mas
psInMas = 3600 * 180 * 1000/np.pi * wvl/39/2
nntrue  = EEbox/psInMas/2
nn2     = nntrue.astype('int')

# config
aoConfig_all = ['NOAO' , 'GLAO' , 'MOAO_poor' , 'MOAO_good']
nCases = len(aoConfig_all)
# instantiating outputs
nPup = 1024
nPsf = 2 * nPup
parser.set('telescope','Resolution',str(nPup))
parser.set('sensor_science','FiedOfView',str(nPsf))
with open(parfile, 'w') as configfile:
            parser.write(configfile)
            
psfSimu     = np.zeros((nCases,nProf,nWvl,nPsf,nPsf))  
FWHM_all    = np.zeros((nCases,nProf,nWvl))
EE400_all   = np.zeros((nCases,nProf,nWvl))


#The GLAO system shall deliver a PSF in H band (1.65 microns) with 
#50% encircled energy within 0.4 arcsec under median atmospheric conditions (TBC)
# anywhere in the patrol field.

#%% LOOPS ON AO CONFIGURATION - CN2 CONDITIONS - 
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
    
        # UPDATING THE ATMOSPHERE
        name = nameProf[i]
        if name == 'JQ1' : 
            Cn2Val = Cn2[:,4]/100
            r0z         = 0.234 
        if name == 'JQ2' : 
            Cn2Val  = Cn2[:,5]/100
            r0z         = 0.178 
        if name == 'JQ3' : 
            Cn2Val  = Cn2[:,6]/100
            r0z         = 0.139
        if name == 'JQ4' : 
            Cn2Val  = Cn2[:,7]/100
            r0z         = 0.097 
        if name == 'Median' : 
            Cn2Val  = Cn2[:,3]/100
            r0z         = 0.157

        zenithAngle = ZenithProf[i]
        seeing = 0.987*0.5/r0z/4.85
        Cn2Weights  = []
        Cn2Heights  = []
        wSpeed      = []
        wDir        = []
        for ii in range(len(Cn2[:,0])):
            Cn2Weights.append(Cn2Val[ii]) 
            Cn2Heights.append(h[ii]) 
            wSpeed.append(v[ii]) 
            wDir.append(vDir[ii])
         
         
        # UPDATE .INI FILE
        parser.set('telescope','ZenithAngle',str(zenithAngle))
        parser.set('atmosphere','Cn2Weights',str(Cn2Weights))
        parser.set('atmosphere','Cn2Heights',str(Cn2Heights))
        parser.set('atmosphere','WindSpeed',str(wSpeed))
        parser.set('atmosphere','WindDirection',str(wDir))
        parser.set('atmosphere','Seeing',str(seeing))  
        with open(parfile, 'w') as configfile:
            parser.write(configfile)
        
        # PSF SIMULATIONS
        fao = fourierModel(parfile,calcPSF=True,verbose=True,display=False,\
                           getErrorBreakDown=False,getFWHM=True,getEncircledEnergy=True,\
                           getEnsquaredEnergy=False, displayContour=False,nyquistSampling=True) 
        
        # CONCATENATING PSFS
        psfSimu[k,i] = np.squeeze(np.transpose(fao.PSF,axes=(2,3,0,1)))
        
        # CONCATENATING FWHM
        FWHM_all[k,i,:] = np.squeeze(np.mean(fao.FWHM,axis=0))
        
        # CONCATENATING EE
        EE       = np.squeeze(fao.EncE)
        EE_multi = (nntrue - nn2)*EE[nn2+1,:] + (nn2+1-nntrue)*EE[nn2,:]
        EE400_all[k,i,:] = np.diag(EE_multi)


#%% SAVING RESULTS
if path_save != '':
    # PSF simu 
    hdu = fits.PrimaryHDU(psfSimu)
    hdul= fits.HDUList(hdu)
    hdul.writeto(path_save+'PSFsimu_on-axis_NOAO_GLAO_MOAOpoor_MOAOgood_5profiles_6wvl.fits',overwrite=True)
    
    # FWHM/EE 
    hdu1 = fits.PrimaryHDU(wvl)
    hdu2 = fits.ImageHDU(FWHM_all)
    hdu3 = fits.ImageHDU(EE400_all)
    hdul= fits.HDUList([hdu1, hdu2, hdu3])
    hdul.writeto(path_save+'wvl_FWHM_EE400mas_on-axis_NOAO_GLAO_MOAOpoor_MOAOgood_5profiles_6wvl.fits',overwrite=True)
        
#%% DISPLAYING RESULTS



plt.close('all')
mpl.rcParams['font.size'] = 22
usetex = True
plt.rcParams.update({
    "text.usetex": usetex,
    "font.family": "serif",
    "font.serif": ["Palatino"],
})  
    
# FWHM vs wvl - Median conditions
fig = plt.figure()
ax = fig.add_subplot(1, 1, 1)
ax.set_xticks(wvl*1e6)
ax.set_xticks(np.sort(list(wvl) + list((wvl[1:] + wvl[:-1]) / 2))*1e6, minor=True)
ax.set_yticks(np.linspace(0,700,21),minor=True)
ax.set_yticks(np.linspace(0,700,11))

plt.plot(wvl*1e6,FWHM_all[0,-1,:],'bs--',label='NOAO')
plt.plot(wvl*1e6,FWHM_all[1,-1,:],'rs--',label='GLAO 4 LGSs')
plt.plot(wvl*1e6,FWHM_all[2,-1,:],'gs--',label='MOAO 4 LGS - 4 NGS at 5\'')
plt.plot(wvl*1e6,FWHM_all[3,-1,:],'ks--',label='MOAO 4 LGS - 4 NGS at 2.5\'')
plt.plot(wvl*1e6,1e3*0.987*wvl*1e6/(0.157*(wvl/500e-9)**1.2)/4.85,'k:',label='Seeing')
plt.plot(wvl*1e6,206264.8*1e3*wvl/39,'b:',label='Diffraction')
plt.xlabel('Wavelength [$\mu$m]')
plt.ylabel('FWHM  - Median profile [mas]')
plt.legend(loc=1,prop={'size': 14})
ax.grid(which='minor', alpha=0.2)
ax.grid(which='major', alpha=0.5)


#%
# EE vs wvl - Median conditions
fig = plt.figure()
ax = fig.add_subplot(1, 1, 1)
ax.set_xticks(wvl*1e6)
ax.set_xticks(np.sort(list(wvl) + list((wvl[1:] + wvl[:-1]) / 2))*1e6, minor=True)
ax.set_yticks(np.linspace(10,70,25),minor=True)
ax.set_yticks(np.linspace(10,70,13))

plt.plot(wvl*1e6,EE400_all[0,-1,:],'bs--',label='NOAO')
plt.plot(wvl*1e6,EE400_all[1,-1,:],'rs--',label='GLAO 4 LGSs')
plt.plot(wvl*1e6,EE400_all[2,-1,:],'gs--',label='MOAO 4 LGS - 4 NGS at 5\'')
plt.plot(wvl*1e6,EE400_all[3,-1,:],'ks--',label='MOAO 4 LGS - 4 NGS at 2.5\'')
plt.xlabel('Wavelength [$\mu$m]')
plt.ylabel('Encircled Energy at '+str(EEbox)+' mas  - Median profile [\%]')
plt.legend(loc=2,prop={'size': 14})
ax.grid(which='minor', alpha=0.2)
ax.grid(which='major', alpha=0.5)

#%
# FWHM vs profile - Median conditions
fig = plt.figure()
ax = fig.add_subplot(1, 1, 1)
ax.set_xticks(wvl*1e6)
ax.set_xticks(np.sort(list(wvl) + list((wvl[1:] + wvl[:-1]) / 2))*1e6, minor=True)
ax.set_yticks(np.linspace(0,1000,21),minor=True)
ax.set_yticks(np.linspace(0,1000,11))
plt.plot(wvl*1e6,FWHM_all[1,0,:],'bs--',label='JQ1')
plt.plot(wvl*1e6,FWHM_all[1,1,:],'rs--',label='JQ2')
plt.plot(wvl*1e6,FWHM_all[1,2,:],'gs--',label='JQ3')
plt.plot(wvl*1e6,FWHM_all[1,3,:],'ms--',label='JQ4')
plt.plot(wvl*1e6,FWHM_all[1,4,:],'ks--',label='Median')
plt.xlabel('Wavelength [$\mu$m]')
plt.ylabel('GLAO FWHM [mas]')
plt.legend()
plt.legend(loc=1,prop={'size': 14})
ax.grid(which='minor', alpha=0.2)
ax.grid(which='major', alpha=0.5)
#%
# FWHM vs profile - Median conditions
fig = plt.figure()
ax = fig.add_subplot(1, 1, 1)
ax.set_xticks(wvl*1e6)
ax.set_xticks(np.sort(list(wvl) + list((wvl[1:] + wvl[:-1]) / 2))*1e6, minor=True)
ax.set_yticks(np.linspace(5,70,27),minor=True)
ax.set_yticks(np.linspace(5,70,14))
plt.plot(wvl*1e6,EE400_all[1,0,:],'bs--',label='JQ1')
plt.plot(wvl*1e6,EE400_all[1,1,:],'rs--',label='JQ2')
plt.plot(wvl*1e6,EE400_all[1,2,:],'gs--',label='JQ3')
plt.plot(wvl*1e6,EE400_all[1,3,:],'ms--',label='JQ4')
plt.plot(wvl*1e6,EE400_all[1,4,:],'ks--',label='Median')
plt.xlabel('Wavelength [$\mu$m]')
plt.ylabel('GLAO Encircled Energy at '+str(EEbox)+' mas [\%]')
plt.legend()
plt.legend(loc=2,prop={'size': 14})
ax.grid(which='minor', alpha=0.2)
ax.grid(which='major', alpha=0.5)