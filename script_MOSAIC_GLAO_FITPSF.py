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
from fourierModel import fourierModel
import FourierUtils
from maoppy.psfmodel import psffit, Psfao, Moffat
from maoppy.instrument import mosaic

#%% GRAB INI FILE
parfile = '/home/omartin/Projects/MOSAIC/AOsimulations/CODES/MosaicGLAOParams_onaxis_multiwvl.ini'
parser = ConfigParser()
parser.optionxform = str
parser.read(parfile)
#%% OBSERVING CONDITIONS
pathCn2 =  '/home/omartin/Projects/MOSAIC/AOsimulations/DATA/'
fileCn2 ='profil_turbulent_eso.fits'

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
nPup = 1000
nPsf = 2 * nPup
parser.set('telescope','resolution',str(nPup))
parser.set('PSF_DIRECTIONS','psf_FoV',str(nPsf))
with open(parfile, 'w') as configfile:
            parser.write(configfile)
            
psfSimu     = np.zeros((nCases,nProf,nWvl,nPsf,nPsf))  
psfMoff     = np.zeros((nCases,nProf,nWvl,nPsf,nPsf)) 
psfFet      = np.zeros((nCases,nProf,nWvl,nPsf,nPsf)) 
FWHM_all    = np.zeros((nCases,nProf,nWvl,3))
EE300_all   = np.zeros((nCases,nProf,nWvl,3))
xMoff       = np.zeros((nCases,nProf,nWvl,4))
xFet        = np.zeros((nCases,nProf,nWvl,7))

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
        fovNGS = 1 # in arcmin
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
                    
    technical_FoV = 1.5*LGSFoV            
    
        
    # UPDATE .INI FILE
    parser.set('GUIDESTARS_HO','GuideStarZenith_HO',str(GuideStarZenith))
    parser.set('GUIDESTARS_HO','GuideStarAzimuth_HO',str(GuideStarAz))
    parser.set('SENSOR_HO','loopGain',str(loopGain))
    parser.set('DM','OptimizationZenith',str(optimRad))
    parser.set('DM','OptimizationAzimuth',str(optimAz))
    parser.set('DM','OptimizationWeight',str(optim))
    parser.set('atmosphere','L0',str(L0))
    parser.set('PSF_DIRECTIONS','technical_FoV',str(technical_FoV))
    
    
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
        parser.set('telescope','zenithAngle',str(zenithAngle))
        parser.set('atmosphere','Cn2Weights',str(Cn2Weights))
        parser.set('atmosphere','Cn2Heights',str(Cn2Heights))
        parser.set('atmosphere','wSpeed',str(wSpeed))
        parser.set('atmosphere','wDir',str(wDir))
        parser.set('atmosphere','seeing',str(seeing))  
        with open(parfile, 'w') as configfile:
            parser.write(configfile)
        
        # PSF SIMULATIONS
        fao = fourierModel(parfile,calcPSF=True,verbose=False,display=False,getErrorBreakDown=False,getFWHM=False,getEncircledEnergy=True,getEnsquaredEnergy=False, displayContour=False) 
        
        # CONCATENATING PSF
        psfSimu[k,i] = np.transpose(fao.PSF[:,:,0,:],axes=(2,0,1))
        
        # FITTING + GRABBING FWHM and EE IN 300 MAS BOX
        for l in range(nWvl):
            # param
            psInMas = 206264.8*1e3*wvl[l]/39/2
            nntrue = 300/psInMas/2
            nn2    = int(nntrue)
            # Moffat fitting
            sx = 1e3*seeing* (wvl[l]/0.5e-6)**(-1/6)/psInMas/2
            x0  = [sx,sx,0,3]
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
plt.figure()
plt.plot(wvl*1e6,FWHM_all[0,-1,:,0],'b',label='NOAO')
plt.plot(wvl*1e6,FWHM_all[1,-1,:,0],'r',label='GLAO')
plt.plot(wvl*1e6,FWHM_all[2,-1,:,0],'g',label='MOAO')
#plt.plot(wvl*1e6,206264.8*1e3*wvl/39,'k--',label='Diffraction')
plt.xlabel('Wavelength [$\mu$m]')
plt.ylabel('FWHM [mas]')
plt.legend()

# FWHM vs profile - Median conditions
plt.figure()
plt.plot(wvl*1e6,FWHM_all[0,0,:,0],label='JQ1')
plt.plot(wvl*1e6,FWHM_all[0,1,:,0],label='JQ2')
plt.plot(wvl*1e6,FWHM_all[0,2,:,0],label='JQ3')
plt.plot(wvl*1e6,FWHM_all[0,3,:,0],label='JQ4')
plt.plot(wvl*1e6,FWHM_all[0,4,:,0],label='Median')
#plt.plot(wvl*1e6,206264.8*1e3*wvl/39,'k--',label='Diffraction')
plt.xlabel('Wavelength [$\mu$m]')
plt.ylabel('NOAO FWHM [mas]')
plt.legend()
