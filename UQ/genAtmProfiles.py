# -*- coding: utf-8 -*-
"""
Created on Tues Feb 16 15:37:27 2020

genAtmProfiles.py:
    Generate density profiles with both KLE and GRAM for comparison.

@author: Samuel Albert
"""

import numpy as np
import matplotlib.pyplot as plt

import planetaryconstants as constants
from sim import Params
from UQ import getEproblem, getKLEdensfun


plt.close('all')

# =============================================================================
# Create Params input class for Mars
# =============================================================================
params = Params()
params.p = constants.MARS
params.returnTimeVectors = True

# =============================================================================
# Load MarsGRAM density data
# =============================================================================
filename = '../data/Mars_0.1_50000.npz'
densdata = np.load(filename)
densTot = densdata['densTot']
h = densdata['h']
densMean = densdata['densMean']
densCentered = densTot - densMean[:,None]

# =============================================================================
# Generate eigenvalues, eigenvectors to use in KLE function
# =============================================================================
alpha = 1
pfact = 1 # use GRAM dispersions directly

evals, evecs, densSampMean, d, h = getEproblem(filename, alpha, pfact)
params.evals = evals
params.evecs = evecs
params.d = d
params.h = h

d = 200

# =============================================================================
# Generate N profiles of density and plot
# =============================================================================
N = 100

fig1 = plt.figure()
fig2 = plt.figure()
ax1 = fig1.add_subplot(111)
ax2 = fig2.add_subplot(111)

rhoKLEList = []
rhoGRAMList = []

for i in range(N):
    # KLE density
    dFun, _ = getKLEdensfun(evals, evecs, densSampMean, d, h)
    rho = dFun(h)
    rhoKLEList.append(rho)
    rhoPert = rho / densMean - 1
    ax1.plot(rhoPert, h, 'r')
    
    
    # GRAM density
    rho = densTot[:,i]
    rhoGRAMList.append(rho)
    rhoPert = rho / densMean - 1
    ax2.plot(rhoPert, h, 'b')
    
    
ax1.grid()
ax2.grid()
    
    
    
    




































