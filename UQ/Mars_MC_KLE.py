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
from atm import getMarsGRAMDensTableAll
from UQ import getEproblem

# def getQoIParams(params):
#     paramsQ = Params()
#     paramsQ.CD = params.CD
#     paramsQ.BC = params.BC
#     paramsQ.m = params.m
#     paramsQ.efpaWR = params.efpaWR
#     paramsQ.vmagWR = params.vmagWR
#     paramsQ.atmdat = params.atmdat
    
#     return paramsQ

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
# Generate atm function using KLE
# =============================================================================
params.dMode = 'fun'

filename = '../data/Mars_0.1_5000.npz'
alpha = 0.99
pfact = 1 # use GRAM dispersions directly

evals, evecs, densSampMean, d, h = getEproblem(filename, alpha, pfact)




# params.dFun, Ys = getKLEdensfun(evals, evecs, sampMean, d, h)

























