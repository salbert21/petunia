# -*- coding: utf-8 -*-
"""
analyze.py:
    analyze Monte Carlo results
Created on Mon Nov 16 17:44:50 2020

@author: Samuel Albert
"""

import numpy as np
# import scipy.linalg as LA
from scipy import stats
# from spgl1 import spg_bpdn 
import matplotlib.pyplot as plt
# import time

# import constants
# from sim import Params, Outs
# from conversions import RV2LLAEHV

plt.close('all')

# =============================================================================
# Load data
# =============================================================================

# MCfilename = './results/Mars_60000_1022201133.npz'
MCfilename = 'results/Mars_10000_highvar_1116164345.npz'
data = np.load(MCfilename, allow_pickle = True)

outsList = data['outsList']
N = len(outsList)

fpafList = np.asarray([out.fpaf for out in outsList])
engfList = np.asarray([out.engf for out in outsList])
vmagfList = np.asarray([out.vmagf for out in outsList])
rafList = np.asarray([out.raf for out in outsList])
qpeakList = np.asarray([out.qpeak for out in outsList])
gpeakList = np.asarray([out.gpeak for out in outsList])
QloadList = np.asarray([out.Qload for out in outsList])

# =============================================================================
# Trim data
# =============================================================================
noImpact = fpafList > 0
noEscape = engfList < 0
Captured = np.logical_and(noImpact, noEscape)

rafListCap = rafList[Captured]

# =============================================================================
# Histograms
# =============================================================================
bins = 200
plt.figure()
plt.grid()
plt.hist(fpafList, bins)
plt.title('final flight path angle')

plt.figure()
plt.grid()
plt.hist(engfList, bins)
plt.title('final energy')

plt.figure()
plt.grid()
plt.hist(vmagfList, bins)
plt.title('final velocity')

plt.figure()
plt.grid()
plt.hist(rafList, bins)
plt.title('final radius of apoapsis')

plt.figure()
plt.grid()
plt.hist(rafListCap, bins)
plt.title('final radius of apoapsis for captured cases')

plt.figure()
plt.grid()
plt.hist(qpeakList, bins)
plt.title('peak heat rate')

plt.figure()
plt.grid()
plt.hist(gpeakList, bins)
plt.title('peak deceleration')

plt.figure()
plt.grid()
plt.hist(QloadList, bins)
plt.title('total heat load')

























