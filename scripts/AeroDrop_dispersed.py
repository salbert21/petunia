# -*- coding: utf-8 -*-
"""
Created on Mon Feb  8 15:00:25 2021

AeroDrop_dispersed.py:
    Simulates orbiter and probe using FNP(A/E)G guidance under dispersions.

@author: Samuel Albert
"""


import numpy as np
import sys
import matplotlib as mpl
import matplotlib.pyplot as plt

import constants
from conversions import LLAEHV2RV, RV2LLAEHV, VN2Vinf, Vinf2VN, getApses
from sim import Params
from guidance import dynFNPAGPhase1, dynFNPAGPhase2

plt.close('all')
mpl.rcParams["figure.autolayout"] = True
mpl.rcParams.update({'font.size': 15})

### CREATE params INPUT CLASS
params = Params()
params.p = constants.MARS
params.returnTimeVectors = True

### INPUT ATM TABLE - GET ATM TABLE FROM GRAM DATA FILE
params.dMode = 'table'
filename = '../data/dens_Mars_nom.txt'
atmdata = np.genfromtxt(filename, names=True)
params.atmdat = np.array([atmdata['Var_X'], atmdata['DENSAV']])
    
### VEHICLE PARAMS
params.m = 2920 # kg, roughly MSL mass
params.CD = 1.6 # roughly MSL CD
params.LD = 0.25
params.BC = 130

params.A = params.m / (params.BC * params.CD)
params.CL = params.CD * params.LD
params.Rn = np.sqrt(params.A / np.pi) / 2

### WIND-RELATIVE INITIAL STATE
params.lat = 18.38
params.lon = -77.58
params.alt = params.p.halt
params.efpaWR = -12
params.hdaWR = 0
params.vmagWR = 6

### GET OTHER STATE TYPES (and assume t0 = 0)
params.x0, params.vInfvec_N = LLAEHV2RV(params.lat, params.lon,
                                        params.alt, params.efpaWR,
                                        params.hdaWR, params.vmagWR,
                                        params, 0)
params.v0 = Vinf2VN(params.x0, params.vInfvec_N, params, 0)
_, _, _, params.efpa, params.hda, params.vmag = \
    RV2LLAEHV(params.x0, params.v0, params, 0)
    

### TEST FUNCTIONS
t = 0
ts = 154
sig0 = 0
sigd = 180
params.rtol = 1e-10
params.atol = 1e-10
params.hmin = 10
params.hmax = params.p.halt + 1e-7 + 10
params.tf = 3000
xxvec = dynFNPAGPhase1(np.block([params.x0, params.v0]),
                       t, ts, sig0, sigd, params)

raf, rpf = getApses(xxvec[:3,-1], xxvec[3:,-1], params)
print(raf - params.p.rad)

h = np.linalg.norm(xxvec[:3,:], axis = 0)
vmag = np.linalg.norm(xxvec[3:,:], axis = 0)

fig = plt.figure()
ax = fig.add_subplot(111)
ax.plot(vmag, h)
ax.grid()




















