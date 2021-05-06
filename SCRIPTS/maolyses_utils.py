#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu May  6 16:29:18 2021

@author: omartin
"""
import numpy as np

def FromProfileToAtmosphereClass(Cn2_eso,name):
    
    h    = Cn2_eso[:,1]
    v    = Cn2_eso[:,2]
    vDir = Cn2_eso[:,8]*180/np.pi

    if name == 'JQ1' : 
        Cn2Val      = Cn2_eso[:,4]/100
        r0z         = 0.234
        zenithAngle = 30
    if name == 'JQ2' : 
        Cn2Val      = Cn2_eso[:,5]/100
        r0z         = 0.178
        zenithAngle = 30
    if name == 'JQ3' : 
        Cn2Val      = Cn2_eso[:,6]/100
        r0z         = 0.139
        zenithAngle = 30
    if name == 'JQ4' : 
        Cn2Val      = Cn2_eso[:,7]/100
        r0z         = 0.097 
        zenithAngle = 45
    if name == 'Median' : 
        Cn2Val      = Cn2_eso[:,3]/100
        r0z         = 0.157
        zenithAngle = 30

    
    seeing = 0.987*0.5/r0z/4.85
    Cn2Weights  = []
    Cn2Heights  = []
    wSpeed      = []
    wDir        = []
    for ii in range(len(Cn2_eso[:,0])):
        Cn2Weights.append(Cn2Val[ii]) 
        Cn2Heights.append(h[ii]) 
        wSpeed.append(v[ii]) 
        wDir.append(vDir[ii])
        
    return zenithAngle, seeing, Cn2Weights, Cn2Heights, wSpeed, wDir