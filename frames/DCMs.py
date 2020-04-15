# -*- coding: utf-8 -*-
"""
Created on Thu Apr  2 16:02:48 2020

@author: Samuel
"""

def getNS(t,params):
    OM = t * params.p.om
    
    SN = getM3(OM)
    NS = SN.T
    
    return NS

def getSN(t,params):
    OM = t * params.p.om
    
    SN = getM3(OM)
    
    return SN

def getM3(th):
    import numpy as np
    
    M3 = np.array([
        [np.cos(th), np.sin(th), 0],
        [-np.sin(th), np.cos(th), 0],
        [0, 0, 1]
        ])
    
    return M3