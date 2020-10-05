# -*- coding: utf-8 -*-
"""
Created on Mon Oct  5 12:37:27 2020

Mars_MC.py:
    main script for Monte Carlo trials of mars aerocapture/entry for ISG

@author: Samuel Albert
"""

import numpy as np
import time
import matplotlib.pyplot as plt
from scipy.stats import norm, uniform

import constants
import ODEs
from sim import Params, Outs, mainAD
from atm import getMarsGRAMDensTable

tic = time.time()
plt.close('all')

# =============================================================================
# Number of Monte Carlo trials
# =============================================================================
Nmc = 10

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
densAll, densMean, h = getMarsGRAMDensTable(filename, Nmc)
params.dMode = 'table'

# =============================================================================
# Set non-dispersed params
# =============================================================================
BCnom = 120 # kg/m^2
CDnom = 1.59
mnom = 3000
Anom = mnom / (CDnom * BCnom)

params.LD = 0 # L/D = 0 --> CL= 0

# wind-relative initial states
params.inputType = 'wind-relative angles'
params.lat = 0
params.lon = 0
params.alt = params.p.halt
params.hdaWR = 0

efpanom = -10.6 # deg
vmagnom = 6 # km/s

# control state
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
# Monte Carlo trials
# =============================================================================

for i_trial in range(Nmc):
    # generate input realizations
    params.efpaWR = norm.rvs(size = 1, loc = efpanom, scale = 0.2/3)
    params.vmagWR = norm.rvs(size = 1, loc = vmagnom, scale = 1e-2/3)
    params.CD = uniform.rvs(size = 1,
                            loc = CDnom - 0.1*CDnom, scale = 0.2*CDnom)
    params.m = uniform.rvs(size = 1, loc = mnom - 0.05*mnom, scale = 0.1*mnom)
    params.BC = params.m / (params.CD * Anom)
    
    # get atmosphere table for this trial
    params.atmdat = np.array([h,densAll[:,i_trial]])
    params.atmdat = params.atmdat[:,params.atmdat[0,:].argsort()]
    
    # run sim
    outs = Outs()
    outs = mainAD(params, tspan, events, outs)
    
    
    # print status
    print('Trial {}: fpaf = {}, engf = {}'.format(i_trial,
                                                  outs.fpaf, outs.engf))








