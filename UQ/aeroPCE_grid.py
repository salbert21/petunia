# -*- coding: utf-8 -*-
"""
aeroPCE_grid.py:
    generates PCE solve for aerocapture data for multiple settings
Created on Thu Oct 29 13:34:25 2020

@author: Samuel Albert
"""

import numpy as np
import scipy.linalg as LA
from scipy import stats
from spgl1 import spg_bpdn 
# import matplotlib.pyplot as plt
import time
import sys

import constants
from sim import Params, Outs
# from conversions import RV2LLAEHV

def solveSPGL1(Psi, u_samps, N, sigma):
    delta = sigma * LA.norm(u_samps[:N])
    
    C_SPG, r, g, i = spg_bpdn(Psi[0:N,:], u_samps[0:N],
                                  delta, iter_lim=100000, verbosity=2,
                                  opt_tol=1e-9, bp_tol=1e-9)

    mean_SPG = C_SPG[0]
    var_SPG = np.sum(C_SPG[1:]**2)
    
    return C_SPG, mean_SPG, var_SPG, i

tic = time.time()

# =============================================================================
# Settings
# =============================================================================
MCname = 'Mars_70000_1118142637.npz'
QOI = 'engf'
p = 3
filtered = True
sigmaList = [0.01, 0.001, 0.0001]
NfactList = np.array([0.2, 0.5, 0.9, 1.1, 4])


# =============================================================================
# Load data and set QoI choice
# =============================================================================
MCfilename = './results/' + MCname

data = np.load(MCfilename, allow_pickle = True)

outsList = data['outsList']
N = len(outsList)

# hard-coded expression for d (assumes exactly 4 dispersed inputs plus atm)
d = 4 + data['atm_YsList'].shape[1]

# # Choose QoI:
if QOI == 'fpaf':
    u_samps = np.asarray([out.fpaf for out in outsList])
elif QOI == 'engf':
    u_samps = np.asarray([out.engf for out in outsList])
elif QOI == 'qpeak':
    u_samps = np.asarray([out.qpeak for out in outsList])
else:
    sys.exit('QOI choice not recognized')



if filtered:
    # =========================================================================
    # FILTERED
    # =========================================================================
    # load Psi matrix
    Psifilename = './psifiles/Psi_p' + str(p) + '_FILT_' + MCname
    Psidata = np.load(Psifilename)
    Psi = Psidata['Psi']
    P = Psi.shape[1]
    
    # filter out u_samps associated with impact cases
    fpafList = np.asarray([out.fpaf for out in outsList])
    u_samps = u_samps[fpafList >= 0]
    # =========================================================================
    # =========================================================================
else:
    # =========================================================================
    # NOT FILTERED
    # =========================================================================
    # load Psi matrix
    Psifilename = './psifiles/Psi_p' + str(p) + '_' + MCname
    Psidata = np.load(Psifilename)
    Psi = Psidata['Psi']
    P = Psi.shape[1]
    # =========================================================================
    # =========================================================================
    
# =============================================================================
# Monte Carlo params ("truth" for comparison)
# =============================================================================
## NOTE that if we are filtered, the Monte Carlo "truth" stats are also filtered
ss = stats.describe(u_samps)
mean_MC = ss[2]
var_MC = ss[3]


# =============================================================================
# Loop through sigma and Nfact lists, report stats for each case
# =============================================================================
C_SPGList = []
mean_SPGList = []
meanErrList = []
varErrList = []
var_SPGList = []
runtimeList = []
nitersList = [] # note that this is not equal to the number of line search its!

for sigma in sigmaList:
    for Nfact in NfactList:
        tici = time.time()
        N_SPG = int(Nfact * P)
        # if we have maxed out size of N_SPG
        if N_SPG > len(u_samps):
            N_SPG = len(u_samps)
        
        C_SPG, mean_SPG, var_SPG, inf = solveSPGL1(Psi, u_samps, N_SPG, sigma)
        toci = time.time()
        C_SPGList.append(C_SPG)
        mean_SPGList.append(mean_SPG)
        meanErrList.append(mean_SPG - mean_MC)
        var_SPGList.append(var_SPG)
        varErrList.append(var_SPG - var_MC)
        runtimeList.append(toci - tici)
        nitersList.append(inf['niters'])
        
        # =====================================================================
        # Print
        # =====================================================================
        print('\n\ncompleted run for Psi file: ' + Psifilename)
        print('\nsigma = {}'.format(sigma))
        print('N_SPG/P = {}'.format(N_SPG/P))
        print('SPGL1 mean: {}'.format(mean_SPG))
        print('SPGL1 variance: {}'.format(var_SPG))
        print('Total Monte Carlo mean: {}'.format(mean_MC))
        print('Total Monte Carlo variance: {}'.format(var_MC))
        
        
        
# =============================================================================
# Save data to file for later analysis
# =============================================================================
C_SPGList = np.asarray(C_SPGList)
mean_SPGList = np.asarray(mean_SPGList)
meanErrList = np.asarray(meanErrList)
var_SPGList = np.asarray(var_SPGList)
varErrList = np.asarray(varErrList)
runtimeList = np.asarray(runtimeList)
nitersList = np.asarray(nitersList)

if filtered:
    outname = './results/' + QOI + 'res_p' + str(p) + '_FILT_' + MCname
else:
    outname = './results/' + QOI + 'res_p' + str(p) + '_' + MCname
    
np.savez(outname,
         C_SPGList = C_SPGList,
         mean_SPGList = mean_SPGList,
         meanErrList = meanErrList,
         var_SPGList = var_SPGList,
         varErrList = varErrList,
         runtimeList = runtimeList,
         nitersList = nitersList,
         mean_MC = mean_MC,
         var_MC = var_MC,
         P = P,
         d = d,
         N = N,
         QOI = QOI
         )


toc = time.time()
print('\nTotal run time: {0:.0f} seconds'.format(toc-tic))



















