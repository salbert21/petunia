# -*- coding: utf-8 -*-
"""
Created on Wed Jun 10 13:37:05 2020

@author: Samuel Albert
"""

import numpy as np
import scipy.linalg as LA

from conversions import LLAEHV2RV, RV2LLAEHV
from frames import DCMs

class params:
    def __init__(self, dMode):
        self.dMode = dMode
        
    class p:
        pass


### CONSTANTS
omE = 2 * np.pi / (0.99726968 * 86400) # rad/s rotation rate of Earth
# radE = 6378.1363
radE = 6378
mu = 3.986e5 # km^3/s^2

### ATM MODE
# dMode = 'fun'
dMode = 'table'

params = params(dMode)
params.p.mu = mu
params.p.om = omE
params.p.rad = radE

lat = 40
lon = -160
alt = 300
fpa = -5
hda = 70
vmag = 10
t = 80

rvec, vvec = LLAEHV2RV(lat, lon, alt, fpa, hda, vmag, params, t)
latc, lonc, altc, fpac, hdac, vmagc = RV2LLAEHV(rvec, vvec, params, t)

# fpacheck = np.degrees(np.arcsin(np.dot(rvec,vvec) / (LA.norm(rvec)*LA.norm(vvec))))
# print(fpacheck)

# M3T = DCMs.getM3(np.radians(lon)).T

# uD_S = -np.cos(np.radians(lat)) * np.array([1,0,0]) + \
#        -np.sin(np.radians(lat)) * np.array([0,0,1])
# uD_S = M3T @ uD_S
# uE_S = M3T @ np.array([0,1,0])

# SN = DCMs.getSN(t,params)
# vvec_S = SN @ vvec

# vvecH_S = vvec_S - np.dot(vvec_S,uD_S) / (LA.norm(uD_S)**2) * uD_S

# hdacheck = np.degrees(np.arccos(np.dot(vvecH_S,uE_S) / (LA.norm(vvecH_S) * LA.norm(uE_S))))
# print(hdacheck)













