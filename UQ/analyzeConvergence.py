# -*- coding: utf-8 -*-
"""
analyzeConvergence.py:
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
# MCfilename = 'results/Mars_10000_highvar_1116164345.npz'
MCfilename = './results/Mars_70000_1118142637.npz'
data = np.load(MCfilename, allow_pickle = True)

outsList = data['outsList']
Ntot = len(outsList)

fpafList = np.asarray([out.fpaf for out in outsList])
engfList = np.asarray([out.engf for out in outsList])
vmagfList = np.asarray([out.vmagf for out in outsList])
rafList = np.asarray([out.raf for out in outsList])
qpeakList = np.asarray([out.qpeak for out in outsList])
gpeakList = np.asarray([out.gpeak for out in outsList])
QloadList = np.asarray([out.Qload for out in outsList])

ss = stats.describe(fpafList)
fpaf_mean = ss[2]
fpaf_var = ss[3]

ss = stats.describe(engfList)
engf_mean = ss[2]
engf_var = ss[3]

# =============================================================================
# Get nt trendlines of mean and variance at different sample sizes
# =============================================================================
nt = 15 # number of trendlines to compute
# N values for P = 2:
NList = np.array([ 38.,  95., 171., 209., 760., 1000])

fpafMeanErrListList = []
fpafVarErrListList = []
engfMeanErrListList = []
engfVarErrListList = []

for Ni in NList.astype(int):
    fpafMeanErrList = []
    fpafVarErrList = []
    engfMeanErrList = []
    engfVarErrList = []
    
    for i in range(nt):
        inds = np.random.choice(Ntot, size = Ni, replace = False)
        ss = stats.describe(fpafList[inds])
        fpafMeanErrList.append(abs((ss[2] - fpaf_mean) / fpaf_mean) * 100)
        fpafVarErrList.append(abs((ss[3] - fpaf_var) / fpaf_var) * 100)

        ss = stats.describe(engfList[inds])
        engfMeanErrList.append(abs((ss[2] - engf_mean) / engf_mean) * 100)
        engfVarErrList.append(abs((ss[3] - engf_var) / engf_var) * 100)
        
    fpafMeanErrListList.append(fpafMeanErrList)
    fpafVarErrListList.append(fpafVarErrList)
    engfMeanErrListList.append(engfMeanErrList)
    engfVarErrListList.append(engfVarErrList)
    
fpafMeanErrListList = np.asarray(fpafMeanErrListList)
fpafVarErrListList = np.asarray(fpafVarErrListList)
engfMeanErrListList = np.asarray(engfMeanErrListList)
engfVarErrListList = np.asarray(engfVarErrListList)
    
fig1, ax1 = plt.subplots(1,1)
for i, y in enumerate(fpafMeanErrListList.T):
    ax1.plot(NList, y, '.--', label = 'i = {}'.format(i))

ax1.grid()
ax1.legend()
ax1.set_title('FPAF MEAN')

fig2, ax2 = plt.subplots(1,1)
for i, y in enumerate(fpafVarErrListList.T):
    ax2.plot(NList, y, '.--', label = 'i = {}'.format(i))

ax2.grid()
ax2.legend()
ax2.set_title('FPAF VARIANCE')

fig3, ax3 = plt.subplots(1,1)
for i, y in enumerate(engfMeanErrListList.T):
    ax3.plot(NList, y, '.--', label = 'i = {}'.format(i))

ax3.grid()
ax3.legend()
ax3.set_title('FPAF MEAN')

fig4, ax4 = plt.subplots(1,1)
for i, y in enumerate(engfVarErrListList.T):
    ax4.plot(NList, y, '.--', label = 'i = {}'.format(i))

ax4.grid()
ax4.legend()
ax4.set_title('FPAF VARIANCE')

np.savez('convergenceData.npz',
         fpafMeanErrListList = fpafMeanErrListList,
         fpafVarErrListList = fpafVarErrListList,
         engfMeanErrListList = engfMeanErrListList,
         engfVarErrListList = engfVarErrListList)


    
    

# # =============================================================================
# # Get mean and variance at different sample sizes
# # =============================================================================
# nt = 10000
# # NList = np.array([25, 38, 95, 171, 209, 266, 300, 400, 500, 665, 750, 1000,
# #                   1197, 1463, 2000,5000, 5320, 1e4,2e4,3e4,4e4,5e4,6e4,7e4])
# NList = np.array([ 38.,  95., 171., 209., 760., 1000, 10000, 50000])

# fpafMeanList = []
# fpafMeanErrList = []
# fpafVarList = []
# fpafVarErrList = []
# engfMeanList = []
# engfMeanErrList = []
# engfVarList = []
# engfVarErrList = []

# for Ni in NList.astype(int):
#     fpafMeansNi = []
#     fpafVarsNi = []
#     engfMeansNi = []
#     engfVarsNi = []
#     for i in range(nt):
#         inds = np.random.choice(Ntot, size = Ni, replace = False)
#         ss = stats.describe(fpafList[inds])
#         fpafMeansNi.append(ss[2])
#         fpafVarsNi.append(ss[3])
        
#         ss = stats.describe(engfList[inds])
#         engfMeansNi.append(ss[2])
#         engfVarsNi.append(ss[3])
        
#     fpafMeanNi = np.mean(fpafMeansNi)
#     fpafMeanList.append(fpafMeanNi)
#     fpafMeanErrList.append(abs((fpafMeanNi - fpaf_mean) / fpaf_mean))
    
#     fpafVarNi = np.mean(fpafVarsNi)
#     fpafVarList.append(fpafVarNi)
#     fpafVarErrList.append(abs((fpafVarNi - fpaf_var) / fpaf_var))
    
#     engfMeanNi = np.mean(engfMeansNi)
#     engfMeanList.append(engfMeanNi)
#     engfMeanErrList.append(abs((engfMeanNi - engf_mean) / engf_mean))
    
#     engfVarNi = np.mean(engfVarsNi)
#     engfVarList.append(engfVarNi)
#     engfVarErrList.append(abs((engfVarNi - engf_var) / engf_var))
    
# fpafMeanErrList = np.asarray(fpafMeanErrList)
# fpafVarErrList = np.asarray(fpafVarErrList)
# engfMeanErrList = np.asarray(engfMeanErrList)
# engfVarErrList = np.asarray(engfVarErrList)
    
# # =============================================================================
# # Plots
# # =============================================================================
# fig, ax = plt.subplots(1,1)
# ax.plot(NList, fpafMeanErrList * 100, '.--', label = 'fpaf')
# ax.plot(NList, engfMeanErrList * 100, '.--', label = 'engf')
# ax.grid()
# ax.legend()

# fig2, ax2 = plt.subplots(1,1)
# ax2.plot(NList, fpafVarErrList * 100, '.--', label = 'fpaf')
# ax2.plot(NList, engfVarErrList * 100, '.--', label = 'engf')
# ax2.grid()    
# ax2.legend()


# # =============================================================================
# # Trim data
# # =============================================================================
# noImpact = fpafList > 0
# noEscape = engfList < 0
# Captured = np.logical_and(noImpact, noEscape)

# rafListCap = rafList[Captured]

# # =============================================================================
# # Histograms
# # =============================================================================
# bins = 200
# plt.figure()
# plt.grid()
# plt.hist(fpafList, bins)
# plt.title('final flight path angle')

# plt.figure()
# plt.grid()
# plt.hist(engfList, bins)
# plt.title('final energy')

# plt.figure()
# plt.grid()
# plt.hist(vmagfList, bins)
# plt.title('final velocity')

# plt.figure()
# plt.grid()
# plt.hist(rafList, bins)
# plt.title('final radius of apoapsis')

# plt.figure()
# plt.grid()
# plt.hist(rafListCap, bins)
# plt.title('final radius of apoapsis for captured cases')

# plt.figure()
# plt.grid()
# plt.hist(qpeakList, bins)
# plt.title('peak heat rate')

# plt.figure()
# plt.grid()
# plt.hist(gpeakList, bins)
# plt.title('peak deceleration')

# plt.figure()
# plt.grid()
# plt.hist(QloadList, bins)
# plt.title('total heat load')

























