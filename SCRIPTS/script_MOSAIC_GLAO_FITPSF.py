#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Feb 12 17:45:06 2021

@author: omartin
"""

#%% Import librairies
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
from astropy.io import fits
from configparser import ConfigParser

from aoSystem.fourierModel import fourierModel
import aoSystem.FourierUtils as FourierUtils
#from psfFitting.psfFitting import psfFitting
from maoppy.psfmodel import psffit, Psfao, Moffat
from maoppy.instrument import mosaic
#%% PATHS
parfile = '/home/omartin/Projects/MOSAIC/AOsimulations/MAOLYSES/INI/MosaicGLAOParams_dev.ini'
pathCn2 = '/home/omartin/Projects/MOSAIC/AOsimulations/MAOLYSES/DATA/'
fileCn2 = 'profil_turbulent_eso.fits'
#%% GRAB INI FILE

parser = ConfigParser()
parser.optionxform = str
parser.read(parfile)

#%% OBSERVING CONDITIONS


C          = fits.open(pathCn2+fileCn2)
Cn2        = C[0].data
h          = Cn2[:,1]
v          = Cn2[:,2]
vDir       = Cn2[:,8]*180/np.pi
L0         = 50
nameProf   = ['JQ1','JQ2','JQ3','JQ4','Median']
ZenithProf = np.array([30, 30, 30, 45, 30])
wvl        = np.array([500e-09, 640e-09, 890e-9, 1.2e-06, 1.65e-06, 2.17e-06])
nWvl       = len(wvl)
nProf      = len(nameProf)

# config
aoConfig_all = ['NOAO' , 'GLAO' , 'MOAO_poor' , 'MOAO_good']
nCases = len(aoConfig_all)
# instantiating outputs
nPup = 512
nPsf = 2 * nPup
parser.set('telescope','Resolution',str(nPup))
parser.set('sensor_science','FiedOfView',str(nPsf))
with open(parfile, 'w') as configfile:
            parser.write(configfile)
            
psfSimu     = np.zeros((nCases,nProf,nWvl,nPsf,nPsf))  
psfMoff     = np.zeros((nCases,nProf,nWvl,nPsf,nPsf)) 
psfFet      = np.zeros((nCases,nProf,nWvl,nPsf,nPsf)) 
FWHM_all    = np.zeros((nCases,nProf,nWvl,3))
EE300_all   = np.zeros((nCases,nProf,nWvl,3))
xMoff       = np.zeros((nCases,nProf,nWvl,4))
xFet        = np.zeros((nCases,nProf,nWvl,7))

for k in range(1):
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
    
    else : # WE OPTIMIZE in THE FOV
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
    parser.set('atmosphere','L0',str(L0))
    
    
    for i  in range(1):
    
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
        fao = fourierModel(parfile,calcPSF=True,verbose=False,display=False,\
                           getErrorBreakDown=False,getFWHM=False,getEncircledEnergy=False,\
                           getEnsquaredEnergy=False, displayContour=False,nyquistSampling=True) 
        
        # CONCATENATING PSF
        psfSimu[k,i] = np.transpose(fao.PSF[:,:,0,:],axes=(2,0,1))
        
        # FITTING + GRABBING FWHM and EE IN 300 MAS BOX
        for l in range(1):
            # param
            psInMas = 206264.8*1e3*wvl[l]/39/2
            nntrue = 300/psInMas/2
            nn2    = int(nntrue)
            # Moffat fitting
            sx = 1e3*seeing* (wvl[l]/0.5e-6)**(-1/6)/psInMas/2
            x0  = [sx,sx,0,3]
            #res = psfFitting(psfSimu[k,i,l],psfao,x0,fixed=(False,False,False,False,False,False,False,True,True,True,False,False,False,True))
            res = psffit(psfSimu[k,i,l],Moffat,x0,flux_bck=(True,False),fixed=(False,False,False,False))
            psfMoff[k,i,l] = res.psf/res.psf.sum()
            xMoff[k,i,l] = res.x
            # PSFAO19 fitting
            opd = 0.1*(39/(r0z * (wvl[l]/0.5e-6)**1.2))**(5/3)/np.sqrt(k+1)
            x0  = [r0z*(wvl[l]/0.5e-6)**1.2,1e-5,(2*np.pi*opd*1e-9/wvl[l])**2 ,1e-2,1,0,1.8]
            res = psffit(psfSimu[k,i,l],Psfao,x0,flux_bck=(True,False),system=mosaic,samp=2,fixed=(False,False,False,False,False,False,False),Lext=L0)
            psfFet[k,i,l] = res.psf/res.psf.sum()
            xFet[k,i,l] = res.x
            # FWHM
            FWHM_all[k,i,l,0] = FourierUtils.getFWHM(psfSimu[k,i,l],psInMas,nargout=1)
            FWHM_all[k,i,l,1] = FourierUtils.getFWHM(psfMoff[k,i,l],psInMas,nargout=1)
            FWHM_all[k,i,l,2] = FourierUtils.getFWHM(psfFet[k,i,l],psInMas,nargout=1)
            # ENCIRCLED ENERGY
            EE = FourierUtils.getEncircledEnergy(psfSimu[k,i,l])
            EE300_all[k,i,l,0] = (nntrue - nn2)*EE[nn2+1] + (nn2+1-nntrue)*EE[nn2]
            EE = FourierUtils.getEncircledEnergy(psfMoff[k,i,l])
            EE300_all[k,i,l,1] = (nntrue - nn2)*EE[nn2+1] + (nn2+1-nntrue)*EE[nn2]
            EE = FourierUtils.getEncircledEnergy(psfFet[k,i,l])
            EE300_all[k,i,l,2] = (nntrue - nn2)*EE[nn2+1] + (nn2+1-nntrue)*EE[nn2]

# to get an instance of the model : pp = Psfao((nPup,nPup),system=mosaic,samp=2,Lext=L0)
# check the pupil in PSFAO19 : plt.imshow(abs(pp.system.pupil((nPup,nPup))))
#%% SAVING RESULTS
path_save = '/home/omartin/Projects/MOSAIC/AOsimulations/RESULTS/FIT/'
# PSF simu 
hdu = fits.PrimaryHDU(psfSimu)
hdul= fits.HDUList(hdu)
hdul.writeto(path_save+'PSFsimu_on-axis_NOAO_GLAO_MOAO_5profiles_6wvl.fits')
# PSF Moffat 
hdu1 = fits.PrimaryHDU(psfMoff)
hdu2 = fits.ImageHDU(xMoff)
hdul= fits.HDUList([hdu1,hdu2])
hdul.writeto(path_save+'PSFmoffat_on-axis_NOAO_GLAO_MOAO_5profiles_6wvl.fits')
# PSF Fetick 
hdu1 = fits.PrimaryHDU(psfFet)
hdu2 = fits.ImageHDU(xFet)
hdul= fits.HDUList([hdu1,hdu2])
hdul.writeto(path_save+'PSFfetick_on-axis_NOAO_GLAO_MOAO_5profiles_6wvl.fits')
# FWHM/EE 
hdu1 = fits.PrimaryHDU(wvl)
hdu2 = fits.ImageHDU(FWHM_all)
hdu3 = fits.ImageHDU(EE300_all)
hdul= fits.HDUList([hdu1, hdu2, hdu3])
hdul.writeto(path_save+'wvl_FWHM_EE300mas_on-axis_NOAO_GLAO_MOAO_5profiles_6wvl.fits')

#%% LOADING
redo = False
if redo:
    path_save = '/home/omartin/Projects/MOSAIC/AOsimulations/RESULTS/FIT/'
    psfSimu = fits.open(path_save+'PSFsimu_on-axis_NOAO_GLAO_MOAO_5profiles_6wvl_save.fits')[0].data
    psfSimu_NOAO = fits.open(path_save+'PSFsimu_on-axis_NOAO_5profiles_6wvl.fits')[0].data

    tmp = fits.open(path_save+'PSFmoffat_on-axis_NOAO_GLAO_MOAO_5profiles_6wvl_save.fits')
    psfMoff  = tmp[0].data
    xMoff   = tmp[1].data
    tmp = fits.open(path_save+'PSFmoffat_on-axis_NOAO_5profiles_6wvl.fits')
    psfMoff_NOAO  = tmp[0].data
    xMoff_NOAO   = tmp[1].data
    
    
    tmp = fits.open(path_save+'PSFfetick_on-axis_NOAO_GLAO_MOAO_5profiles_6wvl_save.fits')
    psfFet  = tmp[0].data
    xFet   = tmp[1].data
    tmp = fits.open(path_save+'PSFfetick_on-axis_NOAO_5profiles_6wvl.fits')
    psfFet_NOAO  = tmp[0].data
    xFet_NOAO   = tmp[1].data
    
    
    tmp    = fits.open(path_save+'wvl_FWHM_EE300mas_on-axis_NOAO_GLAO_MOAO_5profiles_6wvl_save.fits')
    wvl    = tmp[0].data
    FWHM_all    = tmp[1].data
    EE300_all    = tmp[2].data
    tmp    = fits.open(path_save+'wvl_FWHM_EE300mas_on-axis_NOAO_5profiles_6wvl.fits')
    FWHM_NOAO    = tmp[1].data
    EE300_NOAO    = tmp[2].data

    psfSimu[0]  = psfSimu_NOAO
    psfMoff[0]  = psfMoff_NOAO
    psfFet[0]   = psfFet_NOAO
    xFet[0]     = xFet_NOAO[0]
    xMoff[0]     = xMoff_NOAO[0]
    FWHM_all[0] = FWHM_NOAO
    EE300_all[0] = EE300_NOAO
    

    hdu = fits.PrimaryHDU(psfSimu)
    hdul= fits.HDUList(hdu)
    hdul.writeto(path_save+'PSFsimu_on-axis_NOAO_GLAO_MOAO_5profiles_6wvl.fits')
    # PSF Moffat 
    hdu1 = fits.PrimaryHDU(psfMoff)
    hdu2 = fits.ImageHDU(xMoff)
    hdul= fits.HDUList([hdu1,hdu2])
    hdul.writeto(path_save+'PSFmoffat_on-axis_NOAO_GLAO_MOAO_5profiles_6wvl.fits')
    # PSF Fetick 
    hdu1 = fits.PrimaryHDU(psfFet)
    hdu2 = fits.ImageHDU(xFet)
    hdul= fits.HDUList([hdu1,hdu2])
    hdul.writeto(path_save+'PSFfetick_on-axis_NOAO_GLAO_MOAO_5profiles_6wvl.fits')
    # FWHM/EE 
    hdu1 = fits.PrimaryHDU(wvl)
    hdu2 = fits.ImageHDU(FWHM_all)
    hdu3 = fits.ImageHDU(EE300_all)
    hdul= fits.HDUList([hdu1, hdu2, hdu3])
    hdul.writeto(path_save+'wvl_FWHM_EE300mas_on-axis_NOAO_GLAO_MOAO_5profiles_6wvl.fits')
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

plt.plot(wvl*1e6,FWHM_all[0,-1,:,0],'bs--',label='NOAO')
plt.plot(wvl*1e6,FWHM_all[1,-1,:,0],'rs--',label='GLAO 4 LGSs')
plt.plot(wvl*1e6,FWHM_all[2,-1,:,0],'gs--',label='MOAO 4 LGS - 4 NGS at 5\'')
plt.plot(wvl*1e6,FWHM_all[3,-1,:,0],'ks--',label='MOAO 4 LGS - 4 NGS at 2.5\'')
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

plt.plot(wvl*1e6,1e2*EE300_all[0,-1,:,0],'bs--',label='NOAO')
plt.plot(wvl*1e6,1e2*EE300_all[1,-1,:,0],'rs--',label='GLAO 4 LGSs')
plt.plot(wvl*1e6,1e2*EE300_all[2,-1,:,0],'gs--',label='MOAO 4 LGS - 4 NGS at 5\'')
plt.plot(wvl*1e6,1e2*EE300_all[3,-1,:,0],'ks--',label='MOAO 4 LGS - 4 NGS at 2.5\'')
plt.xlabel('Wavelength [$\mu$m]')
plt.ylabel('Encircled Energy at 300 mas  - Median profile [\%]')
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
plt.plot(wvl*1e6,FWHM_all[1,0,:,0],'bs--',label='JQ1')
plt.plot(wvl*1e6,FWHM_all[1,1,:,0],'rs--',label='JQ2')
plt.plot(wvl*1e6,FWHM_all[1,2,:,0],'gs--',label='JQ3')
plt.plot(wvl*1e6,FWHM_all[1,3,:,0],'ms--',label='JQ4')
plt.plot(wvl*1e6,FWHM_all[1,4,:,0],'ks--',label='Median')
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
plt.plot(wvl*1e6,1e2*EE300_all[1,0,:,0],'bs--',label='JQ1')
plt.plot(wvl*1e6,1e2*EE300_all[1,1,:,0],'rs--',label='JQ2')
plt.plot(wvl*1e6,1e2*EE300_all[1,2,:,0],'gs--',label='JQ3')
plt.plot(wvl*1e6,1e2*EE300_all[1,3,:,0],'ms--',label='JQ4')
plt.plot(wvl*1e6,1e2*EE300_all[1,4,:,0],'ks--',label='Median')
plt.xlabel('Wavelength [$\mu$m]')
plt.ylabel('GLAO Encircled Energy at 300 mas [\%]')
plt.legend()
plt.legend(loc=2,prop={'size': 14})
ax.grid(which='minor', alpha=0.2)
ax.grid(which='major', alpha=0.5)