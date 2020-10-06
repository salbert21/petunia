# -*- coding: utf-8 -*-
"""
Created on Mon Oct  5 10:56:19 2020

Mars_single.py:
    main script for single runs at Mars of aerocapture/entry scenarios for ISG

@author: Samuel Albert
"""

import numpy as np
import time
import matplotlib.pyplot as plt

import constants
import ODEs
from sim import Params, Outs, mainAD
from atm import getMarsGRAMDensTable, getMarsGRAMDensTableAll

tic = time.time()
plt.close('all')

# =============================================================================
# Create Params input class for Mars
# =============================================================================
params = Params()
params.p = constants.MARS
params.returnTimeVectors = True

# =============================================================================
# Generate atm table for Nmc profiles
# =============================================================================
filename = 'data/Mars_0.1_5000.txt'
Nmc = 1
densAll, densMean, h = getMarsGRAMDensTable(filename, Nmc)
# densAll, densMean, h = getMarsGRAMDensTableAll(filename, Nmc)
params.dMode = 'table'

# =============================================================================
# Load atm table for this run
# =============================================================================
i_trial = 0
params.atmdat = np.array([h,densAll[:,i_trial]])
# make sure atmdata is sorted by ascending altitude
params.atmdat = params.atmdat[:,params.atmdat[0,:].argsort()]


# =============================================================================
# Set vehicle params
# =============================================================================
params.m = 3000 # kg
params.CD = 1.59
params.LD = 0 # L/D = 0 --> CL= 0
params.BC = 120 # kg/m^2

# =============================================================================
# Wind-Relative initial state
# =============================================================================
params.inputType = 'wind-relative angles'
params.lat = 0
params.lon = 0
params.alt = params.p.halt
params.efpaWR = -10.6
params.hdaWR = 0
params.vmagWR = 6 # km/s

# =============================================================================
# Control state
# =============================================================================
params.bank = 0

# =============================================================================
# Time vector and exit conditions
# =============================================================================
tspan = (0,30000)

params.hmin = 20
params.hmax = params.p.halt + 10

event1 = lambda t, y: ODEs.below_min_alt(t, y, params)
event1.terminal = True

event2 = lambda t, y: ODEs.above_max_alt(t, y, params)
event2.terminal = True

events = (event1, event2)

# =============================================================================
# Run sim
# =============================================================================
outs = Outs()
outs = mainAD(params, tspan, events, outs)

# =============================================================================
# Plot results
# =============================================================================
print(outs.fpaf)
print(outs.engf)
print(outs.raf)

alt = np.linalg.norm(outs.rvec_N, axis=0) - params.p.rad #/1e3
vmag = np.linalg.norm(outs.vvec_N, axis=0)

fig = plt.figure(1)
ax = fig.add_subplot(111)
ax.plot(vmag, alt)
ax.set_xlabel('Inertial velocity, km/s')
ax.set_ylabel('Altitude, km')
ax.grid()


















