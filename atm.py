# -*- coding: utf-8 -*-
"""
Created on Thu Apr  2 16:38:00 2020

@author: Samuel
"""

import numpy as np


def getWind(h):
    
    wvec = np.array([0,0,0])
    
    return wvec

def getRho_from_table(atmdat, h):
    rho = np.interp(h, atmdat['alt'], atmdat['density'])
    
    return rho

def getAtmdat(filename):
    atmdat = np.genfromtxt(filename, delimiter=',', names=True,
                           encoding='utf-8-sig')
    
    return atmdat