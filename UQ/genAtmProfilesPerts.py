# -*- coding: utf-8 -*-
"""
Created on Tues Feb 16 15:37:27 2020

genAtmProfilesPerts.py:
    Generate density profiles with both KLE and GRAM for comparison.
    Uses KLE on density variations instead of actual density.

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
filename = '../data/Mars_0.1_50000_Perts.npz'
densdata = np.load(filename)
densPert = densdata['densTot']
h = densdata['h']
densMean = densdata['densMean']

# =============================================================================
# Generate eigenvalues, eigenvectors to use in KLE function
# =============================================================================
alpha = 0.99
pfact = 1 # use GRAM dispersions directly

evals, evecs, densSampMean, d, hKLE = getEproblem(filename, alpha, pfact)
params.evals = evals
params.evecs = evecs
params.d = d
params.h = h

# =============================================================================
# Generate N profiles of density and plot
# =============================================================================
N = 100

fig1 = plt.figure()
# fig2 = plt.figure()
ax1 = fig1.add_subplot(111)
# ax2 = fig2.add_subplot(111)

rhoPertKLEList = []
rhoPertGRAMList = []

for i in range(N):
    # KLE density
    dFun, _ = getKLEdensfun(evals, evecs, densSampMean, d, hKLE)
    rhoPert = dFun(hKLE)
    rhoPertKLEList.append(rhoPert)
    # rhoPert = rho / densMean - 1
    if i < 5:
        ax1.plot(rhoPert, hKLE, 'r', alpha = 0.5)
    
    
    # GRAM density
    rhoPert = densPert[:,i]
    rhoPertGRAMList.append(rhoPert)
    # rhoPert = rho / densMean - 1
    if i < 5:
        ax1.plot(rhoPert, h, 'b', alpha = 0.5)
    
    
ax1.grid()
# ax2.grid()

rhoPertKLE = np.asarray(rhoPertKLEList)
KLEMEAN = np.mean(rhoPertKLE, axis = 0)
KLESTD = np.std(rhoPertKLE, axis = 0)

rhoPertGRAM = np.asarray(rhoPertGRAMList)
GRAMMEAN = np.mean(rhoPertGRAM, axis = 0)
GRAMSTD = np.std(rhoPertGRAM, axis = 0)

ax1.plot(KLEMEAN - 3 * KLESTD, hKLE, 'r--', linewidth = 3)
ax1.plot(KLEMEAN + 3 * KLESTD, hKLE, 'r--', linewidth = 3, label = 'KLE')
ax1.plot(GRAMMEAN - 3 * GRAMSTD, h, 'b--', linewidth = 3)
ax1.plot(GRAMMEAN + 3 * GRAMSTD, h, 'b--', linewidth = 3, label = 'GRAM')
ax1.set_xlabel('rho/rhoMean - 1')
ax1.set_ylabel('h, km')
ax1.set_title('alpha = {:.2f}, d = {:d}, N = {:d}'.format(alpha, d, N))
ax1.legend()



alpha = 0.95
pfact = 1 # use GRAM dispersions directly

evals, evecs, densSampMean, d95, h = getEproblem(filename, alpha, pfact)

alpha = 0.99
evals, evecs, densSampMean, d99, h = getEproblem(filename, alpha, pfact)

fig, ax = plt.subplots(1,1)
ax.plot(evals, '.', alpha = 0.75)
ax.plot(range(evals.shape[0]), np.ones(evals.shape[0]) * evals[d95], '--',
        linewidth = 2, label = r'$\alpha$ = 0.95 cutoff')
ax.plot(range(evals.shape[0]), np.ones(evals.shape[0]) * evals[d99], '-.',
        linewidth = 2, label = r'$\alpha$ = 0.99 cutoff')
ax.set_yscale('log')
ax.grid()
ax.set_xlabel('index i')
ax.set_ylabel('eigenvalues')
ax.legend()
    
    
    




































