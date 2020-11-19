# -*- coding: utf-8 -*-
"""
aeroPCE.py:
    generates PCE solve for aerocapture data
Created on Thu Oct 29 13:34:25 2020

@author: Samuel Albert
"""

import numpy as np
import scipy.linalg as LA
from scipy import stats
from spgl1 import spg_bpdn 
import matplotlib.pyplot as plt
import time

import constants
from sim import Params, Outs
from conversions import RV2LLAEHV

def solveSPGL1(Psi, u_samps, N, sigma):
    delta = sigma * LA.norm(u_samps[:N])
    
    C_SPG, r, g, i = spg_bpdn(Psi[0:N,:], u_samps[0:N],
                                  delta, iter_lim=100000, verbosity=2,
                                  opt_tol=1e-9, bp_tol=1e-9)

    mean_SPG = C_SPG[0]
    var_SPG = np.sum(C_SPG[1:]**2)
    
    return C_SPG, mean_SPG, var_SPG

tic = time.time()


# =============================================================================
# Load data and set QoI choice
# =============================================================================
# MCfilename = './results/Mars_60000_1022201133.npz'
# MCfilename = './results/Mars_10000_highvar_1116164345.npz'
MCfilename = './results/Mars_70000_1118142637.npz'
data = np.load(MCfilename, allow_pickle = True)

outsList = data['outsList']
N = len(outsList)

# # Choose QoI:
# fpaf:
u_samps = np.asarray([out.fpaf for out in outsList])

# engf:
# u_samps = np.asarray([out.engf for out in outsList])

# # qpeak:
# u_samps = np.asarray([out.qpeak for out in outsList])

# # efpa Y list (debug):
# u_samps = np.asarray(data['efpa_YList'])

# # vmag:
# ## NOTE: the below uses inertial velocity, whereas we dispersed atm. rel. vel.
# uHist = [np.linalg.norm(out.vvec_N, axis = 0) for out in outsList]

# # vmag:
# u_samps = np.asarray([out.vmagf for out in outsList])


# # =============================================================================
# # NOT FILTERED
# # =============================================================================
# # load Psi matrix
# # Psifilename = 'Psi_p3_Mars_60000_1022201133.npz'
# # Psifilename = 'Psi_p2_Mars_10000_highvar_1116164345.npz'
# Psifilename = './psifiles/Psi_p3_Mars_70000_1118142637.npz'

# Psidata = np.load(Psifilename)
# Psi = Psidata['Psi']
# P = Psi.shape[1]
# # =============================================================================
# # =============================================================================

# =============================================================================
# FILTERED
# =============================================================================
# load Psi matrix
# Psifilename = 'Psi_p2_FILT_Mars_60000_1022201133.npz'
# Psifilename = 'Psi_p2_FILT_Mars_10000_highvar_1116164345.npz'
Psifilename = './psifiles/Psi_p3_FILT_Mars_70000_1118142637.npz'

Psidata = np.load(Psifilename)
Psi = Psidata['Psi']
P = Psi.shape[1]

# filter out u_samps associated with impact cases
fpafList = np.asarray([out.fpaf for out in outsList])
u_samps = u_samps[fpafList >= 0]
# =============================================================================
# =============================================================================






# # =============================================================================
# # Least-Squares Solve
# # =============================================================================

# N_LST = len(u_samps)
# C_LST = LA.lstsq(Psi[0:N_LST, :], u_samps[0:N_LST])[0]
# mean_LST = C_LST[0]
# var_LST = np.sum(C_LST[1:]**2)

# =============================================================================
# SPGL1 Solve
# =============================================================================
N_SPG = int(0.2*P)
if N_SPG > len(u_samps):
    N_SPG = len(u_samps)
sigma = 0.0001

C_SPG, mean_SPG, var_SPG = solveSPGL1(Psi, u_samps, N_SPG, sigma)


# =============================================================================
# Monte Carlo params ("truth" for comparison)
# =============================================================================
ss = stats.describe(u_samps)
mean_MC = ss[2]
var_MC = ss[3]

# stats for just the data we fed into SPGL1
u_samps_short = u_samps[0:N_SPG]
ss_short = stats.describe(u_samps_short)
mean_MC_short = ss_short[2]
var_MC_short = ss_short[3]

# lenList = [len(u) for u in uHist]
# maxLen = max(lenList)

# =============================================================================
# Print
# =============================================================================
print('\n\ncompleted run for Psi file: ' + Psifilename)
print('\nsigma = {}'.format(sigma))
print('N_SPG/P = {}'.format(N_SPG/P))
print('SPGL1 mean: {}'.format(mean_SPG))
print('SPGL1 variance: {}'.format(var_SPG))
print('Total Monte Carlo mean: {}'.format(mean_MC))
print('Total Monte Carlo variance: {}'.format(var_MC))




toc = time.time()
print('\nTotal run time: {0:.0f} seconds'.format(toc-tic))



















