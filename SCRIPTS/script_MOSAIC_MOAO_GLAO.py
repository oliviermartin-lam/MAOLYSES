#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Nov  4 13:38:45 2020

@author: TFusco 
"""

#%% Import librairies
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
from astropy.io import fits
from configparser import ConfigParser
from fourierModel import fourierModel
import time
# seed the pseudorandom number generator
import random
#%matplotlib inline

#from distutils.spawn import find_executable
#find_executable('tex')

#text.usetex:False

# Define parFiles path
path = 'C:/Users/tfusco/Nextcloud/07-PYTHON/fourierPSF/'
file =  'MosaicParams_dev.ini'
# Load the parFile
parfile = '/home/omartin/Projects/fourierPSF/parFile/MosaicGLAOParams.ini'
parser = ConfigParser()
parser.optionxform = str
parser.read(parfile)

########################################
## GESTIION DES CONDITIONS DE TURBULENCE 
########################################
pathCn2 =  'C:/Users/tfusco/Nextcloud/04-DATA/HARMONI/'
fileCn2 ='profil_turbulent_eso.fits'

C          = fits.open(pathCn2+fileCn2)
Cn2        = C[0].data
h          = Cn2[:,1]
v          = Cn2[:,2]
vDir       = Cn2[:,8]*180/np.pi
L0         = 50
nameProf   = ['JQ1','JQ2','JQ3','JQ4','Median']
ZenithProf = np.array([30, 30, 30, 45, 30])
 
ScienceWavelength = 1650e-9

for i  in [4] :
    # for i in range(0,4) : 
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
       
    
    #################################################
    ## GESTION DES ASPECTS SYSTEMES & OBSERVATIONNELS  
    #################################################
    
    MOAO = 0 # ou O 
    
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
    if MOAO == 1 : 
        #add of  NGS 
        nbNGS= 4
        if nbNGS !=0 :
            for istar in range(nbNGS) :
                random.seed()
                GuideStarZenith.append(random.random()*5*60/2)
                random.seed()
                GuideStarAz.append(random.random()*360)
        print(GuideStarZenith)
        print(GuideStarAz)
    
    parser.set('GUIDESTARS_HO','GuideStarZenith_HO',str(GuideStarZenith))
    parser.set('GUIDESTARS_HO','GuideStarAzimuth_HO',str(GuideStarAz))
    
    
    ## OPTIMISATION DIRECTION 
    #########################
    optimRad     = []
    optimAz      = []
    optim        = []
    
    if MOAO == 1 : 
        random.seed()
        aopt   = random.random()*5/2*60 #arcsec
        random.seed()
        bopt   = random.random()*360 #deg
        print(aopt,bopt)
        optimRad.append(aopt)
        optimAz.append(bopt)
        optim.append(1)
    
    else : 
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
     
            
    ## EVAL DIRECTION 
    #################
    EvalFoV     = 5*60 #arcsec  
    if MOAO == 0 :
        nEval       = 5 #nb point lineaire pour la perf (au total nEval²)
        EvalX       = np.linspace(-EvalFoV/2,EvalFoV/2,num=nEval)#+1e-6 #1/4 de champ because symétrie 
        EvalY       = np.linspace(-EvalFoV/2,EvalFoV/2,num=nEval)#+1e-6
        
        EvalRad     = []
        EvalAz      = []
        
        for i in range(nEval):
            for j in range(nEval): 
                a  = np.sqrt( EvalX[i] * EvalX[i] + EvalY[j] * EvalY[j] )
                if EvalY[j] == 0 : 
                   if EvalX[i] < 0 : 
                       b = -90
                   else : 
                       b = 90
                else :  
                   b  = 180*np.arctan2(EvalX[i],EvalY[j])/np.pi
                EvalRad.append(a)
                EvalAz.append(b)
    else : 
        EvalRad     = optimRad
        EvalAz      = optimAz

    
    
    
    technical_FoV = 1.5*max([EvalFoV, LGSFoV])
    
    nLenslet_HO = 4 #64
    #DmPitchs = [39/nLenslet_HO]
    if MOAO == 1 :
        DmPitchs = [39/4]
    else : 
        DmPitchs = [30/4] # 0.6
    
    SensorFrameRate_HO = 0.1 #500.0
    SensorFrameRate_LO = SensorFrameRate_HO
    
    ###############################
    # Ecriture dans le fichier init
    ###############################
    
    parser.set('SENSOR_HO','nLenslet_HO',str(nLenslet_HO))
    parser.set('SENSOR_HO','SensorFrameRate_HO',str(SensorFrameRate_HO))
    parser.set('SENSOR_LO','SensorFrameRate_HO',str(SensorFrameRate_HO))
    
    parser.set('DM','DmPitchs',str(DmPitchs))
    
    parser.set('telescope','zenithAngle',str(zenithAngle))
    parser.set('atmosphere','Cn2Weights',str(Cn2Weights))
    parser.set('atmosphere','Cn2Heights',str(Cn2Heights))
    parser.set('atmosphere','wSpeed',str(wSpeed))
    parser.set('atmosphere','wDir',str(wDir))
    parser.set('atmosphere','seeing',str(seeing))
    parser.set('atmosphere','L0',str(L0))

    
    parser.set('PSF_DIRECTIONS','technical_FoV',str(technical_FoV))
    
    
    parser.set('DM','OptimizationZenith',str(optimRad))
    parser.set('DM','OptimizationAzimuth',str(optimAz))
    parser.set('DM','OptimizationWeight',str(optim))
    
    parser.set('PSF_DIRECTIONS','ScienceWavelength',str(ScienceWavelength))
    parser.set('PSF_DIRECTIONS','ScienceZenith',str(EvalRad))
    parser.set('PSF_DIRECTIONS','ScienceAzimuth',str(EvalAz))
    
    with open(parfile, 'w') as configfile:
        parser.write(configfile)
    
    #fao = fourierModel(parfile,calcPSF=True,verbose=True,display=True,getErrorBreakDown=False)
    
    t0 = time.time()
    fao = fourierModel(parfile,calcPSF=True,verbose=False,display=False,getErrorBreakDown=False,getFWHM=True,getEncircledEnergy=True,getEnsquaredEnergy=False, displayContour=True) 
    PSD = fao.powerSpectrumDensity()
    tcalc = (time.time() - t0) 
    print("Required time for total calculation (s)\t : {:f}".format(tcalc))
    #nmRMS= fao.wfeTot
    SR   = fao.SR 
    FWHM = fao.FWHM
    EE   = fao.EncE
    EEtrue = []
    
    eewidthIDiamInMas = np.array([75, 150, 300]) 
    for i  in range(len(eewidthIDiamInMas)) : 
        nntrue      = eewidthIDiamInMas[i]/fao.psInMas/2
        nn2         = int(nntrue)
        add         = (nntrue - nn2)*EE[nn2+1,:,0] + (nn2+1-nntrue)*EE[nn2,:,0]
        EEtrue.append(add)  
    
    if MOAO == 1 : 
        parser = ConfigParser()
        parser.optionxform = str
        parser.read(parfile)
        EvalRad     = []
        EvalAz      = []
        EvalRad.append(aopt)
        EvalAz.append(bopt)
        parser.set('PSF_DIRECTIONS','ScienceZenith',str(EvalRad))
        parser.set('PSF_DIRECTIONS','ScienceAzimuth',str(EvalAz))
        with open(parfile, 'w') as configfile:
                parser.write(configfile)
        fao = fourierModel(parfile,calcPSF=True,verbose=False,display=False,getErrorBreakDown=False,getFWHM=True,getEncircledEnergy=True,getEnsquaredEnergy=False, displayContour=False) 
        PSD = fao.powerSpectrumDensity()
        SR   = fao.SR 
        FWHMMOAO = fao.FWHM
        EE   = fao.EncE
        EEMOAOtrue = []
        eewidthIDiamInMas = np.array([80, 120,150, 180]) 
        for i  in range(len(eewidthIDiamInMas)) : 
            nntrue      = eewidthIDiamInMas[i]/fao.psInMas/2
            nn2         = int(nntrue)
            add         = (nntrue - nn2)*EE[nn2+1,:,0] + (nn2+1-nntrue)*EE[nn2,:,0]
            EEMOAOtrue.append(add)  
        tcalc = (time.time() - t0) 
        print('<FWHMx>        at %4.0fµm \t = \t %4.0f %s'%(fao.wvlSrc[0]*1e6,FWHMMOAO[0],'\tmas'))
        print('<FWHMy>        at %4.0fµm \t = \t %4.0f %s'%(fao.wvlSrc[0]*1e6,FWHMMOAO[1],'\tmas'))
        print('<EE> in 75 mas at %4.0fµm \t = \t %4.0f %s'%(fao.wvlSrc[0]*1e6,EEMOAOtrue[0],' \t%'))
        print('<EE> in 150mas at %4.0fµm \t = \t %4.0f %s'%(fao.wvlSrc[0]*1e6,EEMOAOtrue[1],' \t%'))
        print('<EE> in 300mas at %4.0fµm \t = \t %4.0f %s'%(fao.wvlSrc[0]*1e6,EEMOAOtrue[2],' \t%'))
    else :
        print('<FWHMx>        at %4.0fµm \t = \t %4.0f \t\u00B1 %4.0f%s'%(fao.wvlSrc[0]*1e6,np.mean(FWHM[0]),np.sqrt(np.var(FWHM[0])),'\tmas'))
        print('<FWHMy>        at %4.0fµm \t = \t %4.0f \t\u00B1 %4.0f%s'%(fao.wvlSrc[0]*1e6,np.mean(FWHM[1]),np.sqrt(np.var(FWHM[1])),'\tmas'))
        print('<EE> in 75 mas at %4.0fµm \t = \t %4.0f \t\u00B1 %4.0f%s'%(fao.wvlSrc[0]*1e6,np.mean(EEtrue[0]),np.sqrt(np.var(EEtrue[0])),' \t%'))
        print('<EE> in 150mas at %4.0fµm \t = \t %4.0f \t\u00B1 %4.0f%s'%(fao.wvlSrc[0]*1e6,np.mean(EEtrue[1]),np.sqrt(np.var(EEtrue[1])),' \t%'))
        print('<EE> in 300mas at %4.0fµm \t = \t %4.0f \t\u00B1 %4.0f%s'%(fao.wvlSrc[0]*1e6,np.mean(EEtrue[2]),np.sqrt(np.var(EEtrue[2])),' \t%'))

#PSFMOAO = fao.PSF 
#EEMOAO = fao.EE 
               
#plt.plot(EEMOAO[:,:,0,0])

#plt.imshow(np.log10(PSFMOAO[:,:,0,0]))

# #%% chnage the parameters
# wfe = np.zeros(nP)
# for k in range(nP):
#     print('Doing case {:d} over {:d}\n'.format(k+1,nP))
#     # update the parameter in the .ini file
#     parser.set('DM','DmPitchs',str([dm_pitch[k]]))
#     with open(parfile, 'w') as configfile:
#         parser.write(configfile)
#     # Instantiating the Fourier model object
#     # Note : calcPSF =false -> only instantiation of matrices
#     fao = fourierModel(parfile,calcPSF=False,verbose=False,display=False,getErrorBreakDown=False)
#     # Get the PSD in nm^2 (accounts for the pixel scale already)
#     PSD = fao.powerSpectrumDensity()
#     # wavefront error
#     wfe[k] = np.sqrt(PSD.sum())

# #%% Save results
# fits.writeto('wfe.fits',wfe,overwrite=True) 
   
# #%% PLot
#     plt.close('all')
# # to have nice latex fonts
# plt.rcParams.update({
#     "text.usetex": True,
#     "font.family": "serif",
#     "font.serif": ["Palatino"],
# })
# mpl.rcParams['font.size'] = 18

# plt.figure(figsize=(10,10))
# plt.plot(dm_pitch,wfe,'bs--',label='Total error')
# plt.plot(dm_pitch,np.sqrt(0.23*(dm_pitch/fao.atm.r0)**(5/3))*fao.wvlSrc[0]*1e9/2/np.pi,'r-',label='Pure fitting error')
# plt.ylabel('Wavefront error (nm)')
# plt.xlabel('DM actuators pitch (m)')
# plt.legend()
