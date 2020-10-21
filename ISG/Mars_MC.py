# -*- coding: utf-8 -*-
"""
Created on Mon Oct  5 12:37:27 2020

Mars_MC.py:
    main script for Monte Carlo trials of mars aerocapture/entry for ISG

@author: Samuel Albert
"""

import numpy as np
import time
from datetime import datetime
import matplotlib.pyplot as plt
from scipy.stats import norm, uniform

import constants
import ODEs
from sim import Params, Outs, mainAD
from atm import getMarsGRAMDensTableAll

def getQoIParams(params):
    paramsQ = Params()
    paramsQ.CD = params.CD
    paramsQ.BC = params.BC
    paramsQ.m = params.m
    paramsQ.efpaWR = params.efpaWR
    paramsQ.vmagWR = params.vmagWR
    paramsQ.atmdat = params.atmdat
    
    return paramsQ

tic = time.time()
datestring = datetime.now().strftime('%m%d%H%M%S')
plt.close('all')

# =============================================================================
# Number of Monte Carlo trials
# =============================================================================
Nmc = 5000

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
densAll, densMean, h = getMarsGRAMDensTableAll(filename, Nmc)
params.dMode = 'table'

# =============================================================================
# Set non-dispersed params
# =============================================================================
params.LD = 0 # L/D = 0 --> CL= 0

# wind-relative initial states
params.inputType = 'wind-relative angles'
params.lat = 0
params.lon = 0
params.alt = params.p.halt
params.hdaWR = 0

# control state
params.bank = 0

# =============================================================================
# Define mean and variance (or bounds if uniform) for dispersed inputs
# =============================================================================
BCnom = 120 # kg/m^2
CDnom = 1.59
CD_LB = CDnom - 0.1 * CDnom
CD_UB = CDnom + 0.1 * CDnom
mnom = 3000
m_LB = mnom - 0.05 * mnom
m_UB = mnom + 0.05 * mnom
Anom = mnom / (CDnom * BCnom)

efpanom = -10.6 # deg
efpa_sig = 0.2/3
vmagnom = 6 # km/s
vmag_sig = 1e-2/3


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

paramsList = []
outsList = []
outname = './results/' + params.p.name + '_' + str(Nmc) + '_' + datestring

for i_trial in range(Nmc):
    # generate input realizations
    params.efpaWR = norm.rvs(size = 1, loc = efpanom, scale = efpa_sig)[0]
    params.vmagWR = norm.rvs(size = 1, loc = vmagnom, scale = vmag_sig)[0]
    params.CD = uniform.rvs(size = 1,
                            loc = CD_LB, scale = 2*(CD_UB - CDnom))[0]
    params.m = uniform.rvs(size = 1,
                           loc = m_LB, scale = 2*(m_UB - mnom))[0]
    params.BC = params.m / (params.CD * Anom)
    
    # get atmosphere table for this trial
    params.atmdat = np.array([h,densAll[:,i_trial]])
    params.atmdat = params.atmdat[:,params.atmdat[0,:].argsort()]
    
    # run sim
    print('\nTrial {}'.format(i_trial+1))
    outs = Outs()
    outsList.append(mainAD(params, tspan, events, outs))
    paramsList.append(getQoIParams(params))
    
    
    # every 50 trials, save all results to a file
    if not i_trial % 50:
        
        ## Save results to a file            
        np.savez(outname,
                 paramsList = paramsList,
                 outsList = outsList,
                 CD_LB = CD_LB,
                 CD_UB = CD_UB,
                 m_LB = m_LB,
                 m_UB = m_UB,
                 efpanom = efpanom,
                 efpa_sig = efpa_sig,
                 vmagnom = vmagnom,
                 vmag_sig = vmag_sig,
                 paramsf = params
                 )
    
    
    # # print status
    # print('Trial {}: fpaf = {}, engf = {}'.format(i_trial,
    #                                               outs.fpaf, outs.engf))
    
    
    
# save final values to file
## Save results to a file
    
np.savez(outname,
         paramsList = paramsList,
         outsList = outsList,
         CD_LB = CD_LB,
         CD_UB = CD_UB,
         m_LB = m_LB,
         m_UB = m_UB,
         efpanom = efpanom,
         efpa_sig = efpa_sig,
         vmagnom = vmagnom,
         vmag_sig = vmag_sig,
         paramsf = params
         )


toc = time.time()
print('\nAnalysis complete!')
print('{} trajectories simulated'.format(Nmc))
print('Total run time: {0:.0f} seconds'.format(toc-tic))






