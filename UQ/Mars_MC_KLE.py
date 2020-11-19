# -*- coding: utf-8 -*-
"""
Created on Mon Oct  5 12:37:27 2020

Mars_MC_KLE.py:
    main script for Monte Carlo trials of mars aerocapture/entry.
    Uses KLE model for density, for subsequent PCE analysis.

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
# from atm import getMarsGRAMDensTableAll
from UQ import getEproblem, getKLEdensfun

def getQoIParams(params):
    paramsQ = Params()
    paramsQ.CD = params.CD
    paramsQ.BC = params.BC
    paramsQ.m = params.m
    paramsQ.efpaWR = params.efpaWR
    paramsQ.vmagWR = params.vmagWR
    paramsQ.m_Y = params.m_Y
    paramsQ.CD_Y = params.CD_Y
    paramsQ.efpa_Y = params.efpa_Y
    paramsQ.vmag_Y = params.vmag_Y
    paramsQ.atm_Ys = params.atm_Ys
    
    return paramsQ

tic = time.time()
datestring = datetime.now().strftime('%m%d%H%M%S')
plt.close('all')

# =============================================================================
# Number of Monte Carlo trials
# =============================================================================
Nmc = 10000

# =============================================================================
# Create Params input class for Mars
# =============================================================================
params = Params()
params.p = constants.MARS
params.returnTimeVectors = True

# =============================================================================
# Generate eigenvalues, eigenvectors to use in KLE function
# =============================================================================
params.dMode = 'fun'

filename = '../data/Mars_0.1_50000.npz'
alpha = 0.95
pfact = 1 # use GRAM dispersions directly

evals, evecs, densSampMean, d, h = getEproblem(filename, alpha, pfact)
params.evals = evals
params.evecs = evecs
params.d = d
params.h = h


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
# Define mean and variance for dispersed inputs
# =============================================================================
# TODO: using all normal distributions for now. should update some to uniform.
BCmean = 129 # kg/m^2

# change to U[CD_LB, CD_UB] later
CDmean = 1.59
CD_LB = CDmean - 0.1 * CDmean
CD_UB = CDmean + 0.1 * CDmean
CDstd = (CD_UB - CD_LB) / np.sqrt(12)

mmean = 3000 # kg
m_LB = mmean - 0.05 * mmean
m_UB = mmean + 0.05 * mmean
mstd = (m_UB - m_LB) / np.sqrt(12)

Anom = mmean / (CDmean * BCmean)

efpamean = -10.6 # deg
efpastd = 0.2/3

vmagmean = 6 # km/s
vmagstd = 1e-2/3

# #### super high dispersions
# # change to U[CD_LB, CD_UB] later
# CDmean = 1.59
# CD_LB = CDmean - 0.3 * CDmean
# CD_UB = CDmean + 0.3 * CDmean
# CDstd = (CD_UB - CD_LB) / np.sqrt(12)

# mmean = 3000 # kg
# m_LB = mmean - 0.15 * mmean
# m_UB = mmean + 0.15 * mmean
# mstd = (m_UB - m_LB) / np.sqrt(12)

# Anom = mmean / (CDmean * BCmean)

# efpamean = -10.6 # deg
# efpastd = 0.2

# vmagmean = 6 # km/s
# vmagstd = 1e-2

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

m_YList = []
CD_YList = []
efpa_YList = []
vmag_YList = []
atm_YsList = []

outname = './results/' + params.p.name + '_' + str(Nmc) + '_' + datestring

for i_trial in range(Nmc):
    # generate density function for this trial
    params.dFun, params.atm_Ys = getKLEdensfun(evals, evecs, densSampMean, d, h)
    
    # generate input realizations
    params.efpa_Y = norm.rvs(size = 1).item(0)
    params.efpaWR = efpamean + efpastd * params.efpa_Y
    
    params.vmag_Y = norm.rvs(size = 1).item(0)
    params.vmagWR = vmagmean + vmagstd * params.vmag_Y
    
    params.CD_Y = norm.rvs(size = 1).item(0)
    params.CD = CDmean + CDstd * params.CD_Y
    
    params.m_Y = norm.rvs(size = 1).item(0)
    params.m = mmean + mstd * params.m_Y
    
    params.BC = params.m / (params.CD * Anom)
    params.A = Anom
    
    
    # run sim
    print('\nTrial {}'.format(i_trial+1))
    outs = Outs()
    outsList.append(mainAD(params, tspan, events, outs))
    paramsList.append(getQoIParams(params))
    m_YList.append(params.m_Y)
    CD_YList.append(params.CD_Y)
    efpa_YList.append(params.efpa_Y)
    vmag_YList.append(params.vmag_Y)
    atm_YsList.append(params.atm_Ys)
    
    
    # every 50 trials, save all results to a file
    if not i_trial % 50:
        
        ## Save results to a file            
        np.savez(outname,
                 paramsList = paramsList,
                 outsList = outsList,
                 m_YList = m_YList,
                 CD_YList = CD_YList,
                 vmag_YList = vmag_YList,
                 efpa_YList = efpa_YList,
                 atm_YsList = atm_YsList,
                 CDmean = CDmean,
                 CDstd = CDstd,
                 mmean = mmean,
                 mstd = mstd,
                 efpamean = efpamean,
                 efpastd = efpastd,
                 vmagmean = vmagmean,
                 vmagstd = vmagstd
                 )
    
    
    # # print status
    # print('Trial {}: fpaf = {}, engf = {}'.format(i_trial,
    #                                               outs.fpaf, outs.engf))




# save final results
np.savez(outname,
         paramsList = paramsList,
         outsList = outsList,
         m_YList = m_YList,
         CD_YList = CD_YList,
         vmag_YList = vmag_YList,
         efpa_YList = efpa_YList,
         atm_YsList = atm_YsList,
         CDmean = CDmean,
         CDstd = CDstd,
         mmean = mmean,
         mstd = mstd,
         efpamean = efpamean,
         efpastd = efpastd,
         vmagmean = vmagmean,
         vmagstd = vmagstd
         )


toc = time.time()
print('\nAnalysis complete!')
print('{} trajectories simulated'.format(Nmc))
print('Total run time: {0:.0f} seconds'.format(toc-tic))

























