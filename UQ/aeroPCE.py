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

import constants
from sim import Params, Outs
from conversions import RV2LLAEHV

def solveSPGL1(Psi, u_samps, N, sigma):
    delta = sigma * LA.norm(u_samps[:N])
    
    C_SPG, r, g, i = spg_bpdn(Psi[0:N,:], u_samps[0:N],
                                  delta, iter_lim=120000, verbosity=2,
                                  opt_tol=1e-9, bp_tol=1e-9)

    mean_SPG = C_SPG[0]
    var_SPG = np.sum(C_SPG[1:]**2)
    
    return C_SPG, mean_SPG, var_SPG

# =============================================================================
# Load data and Psi matrix
# =============================================================================

MCfilename = './results/Mars_60000_1022201133.npz'
data = np.load(MCfilename, allow_pickle = True)
# paramsList = data['paramsList']

# Psifilename = 'Psi_p3_Mars_60000_1022201133.npz'
Psifilename = 'Psi_p2_Mars_60000_1022201133.npz'
Psidata = np.load(Psifilename)

Psi = Psidata['Psi']
P = Psi.shape[1]

# =============================================================================
# Set QoI choice
# =============================================================================
outsList = data['outsList']
N = len(outsList)

# vmag:
## NOTE: the below uses inertial velocity, whereas we dispersed atm. rel. vel.
uHist = [np.linalg.norm(out.vvec_N, axis = 0) for out in outsList]

# rmag:
# uHist = [np.linalg.norm(out.rvec_N, axis = 0) for out in outsList]

# # fpa:
# # make a dummy params just to use an function argument
# params = Params()
# params.p = constants.MARS

# uHist = []
# for out in outsList:
#     fpahist = np.empty(len(out.tvec))
#     fpahist[:] = np.NaN
#     for ind, ti in enumerate(out.tvec):
#         lat, lon, alt, fpaWR, hdaWR, vmagWR = RV2LLAEHV(out.rvec_N[:,ind],
#                                                          out.vvec_N[:,ind],
#                                                          params, out.tvec[ind])
#         fpahist[ind] = fpaWR
#     uHist.append(fpahist)
        
        
uLen = [len(uHisti) for uHisti in uHist]
M = max(uLen)

u_samps = np.empty([N, M])
u_samps[:] = np.NaN
for i, mm in enumerate(uLen):
    u_samps[i, 0:mm] = uHist[i]
    if mm < M: # if some empty columns left over for this row
        u_samps[i, mm:] = uHist[i][-1] # pad with final value


# # =============================================================================
# # Least-Squares Solve
# # =============================================================================

# N_LST = len(u_samps)
# C_LST = LA.lstsq(Psi[0:N_LST, :], u_samps[0:N_LST, :])[0]
# mean_LST = C_LST[0]
# var_LST = np.sum(C_LST[1])

# # =============================================================================
# # SPGL1 Solve
# # =============================================================================
# N_SPG = int(0.9*P)
# sigma = 0.05

# C_SPG, mean_SPG, var_SPG = solveSPGL1(Psi, u_samps, N_SPG, sigma)




lenList = [len(u) for u in uHist]
maxLen = max(lenList)
























