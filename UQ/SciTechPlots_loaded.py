# -*- coding: utf-8 -*-
"""
SciTechPlots_loaded.py:
    plots compiled PCE data for SciTech 2021 paper, from data files
Created on Tue Nov 17 16:26:56 2020

@author: Samuel Albert
"""

import numpy as np
from scipy import stats
import matplotlib.pyplot as plt


plt.close('all')

# =============================================================================
# Load all data
# =============================================================================
# set MC run of interest
MCname = 'Mars_70000_1118142637.npz'

# hard-code sigma and Nfact lists
sigmaList = [0.01, 0.001, 0.0001]
NfactList = np.array([0.2, 0.5, 0.9, 1.1, 4])

# fpaf unfiltered p = 2
PCEfilename = './results/fpafres_p2_' + MCname
# PCEfilename = './results/fpafres_multi_p2_' + MCname
data = np.load(PCEfilename, allow_pickle = True)
fpaf_meanErr_all_2 = data['meanErrList']
fpaf_varErr_all_2  = data['varErrList']
fpaf_runtime_all_2 = data['runtimeList']
fpaf_niters_all_2  = data['nitersList']
fpaf_meanMC_all    = data['mean_MC']
fpaf_varMC_all     = data['var_MC']
P_2 = data['P']
d = data['d']
N = data['N']

# fpaf filtered p = 2
PCEfilename = './results/fpafres_p2_FILT_' + MCname
# PCEfilename = './results/fpafres_multi_p2_FILT_' + MCname
data = np.load(PCEfilename, allow_pickle = True)
fpaf_meanErr_filt_2 = data['meanErrList']
fpaf_varErr_filt_2  = data['varErrList']
fpaf_runtime_filt_2 = data['runtimeList']
fpaf_niters_filt_2  = data['nitersList']
fpaf_meanMC_filt    = data['mean_MC']
fpaf_varMC_filt     = data['var_MC']

# fpaf unfiltered p = 3
PCEfilename = './results/fpafres_p3_' + MCname
data = np.load(PCEfilename, allow_pickle = True)
fpaf_meanErr_all_3 = data['meanErrList']
fpaf_varErr_all_3  = data['varErrList']
fpaf_runtime_all_3 = data['runtimeList']
fpaf_niters_all_3  = data['nitersList']
P_3 = data['P']

# fpaf filtered p = 3
PCEfilename = './results/fpafres_p3_FILT_' + MCname
data = np.load(PCEfilename, allow_pickle = True)
fpaf_meanErr_filt_3 = data['meanErrList']
fpaf_varErr_filt_3  = data['varErrList']
fpaf_runtime_filt_3 = data['runtimeList']
fpaf_niters_filt_3  = data['nitersList']

# engf unfiltered p = 2
PCEfilename = './results/engfres_p2_' + MCname
# PCEfilename = './results/engfres_multi_p2_' + MCname
data = np.load(PCEfilename, allow_pickle = True)
engf_meanErr_all_2 = data['meanErrList']
engf_varErr_all_2  = data['varErrList']
engf_runtime_all_2 = data['runtimeList']
engf_niters_all_2  = data['nitersList']
engf_meanMC_all    = data['mean_MC']
engf_varMC_all     = data['var_MC']
P_2 = data['P']
d = data['d']
N = data['N']

# engf filtered p = 2
PCEfilename = './results/engfres_p2_FILT_' + MCname
# PCEfilename = './results/engfres_multi_p2_FILT_' + MCname
data = np.load(PCEfilename, allow_pickle = True)
engf_meanErr_filt_2 = data['meanErrList']
engf_varErr_filt_2  = data['varErrList']
engf_runtime_filt_2 = data['runtimeList']
engf_niters_filt_2  = data['nitersList']
engf_meanMC_filt    = data['mean_MC']
engf_varMC_filt     = data['var_MC']

# engf unfiltered p = 3
PCEfilename = './results/engfres_p3_' + MCname
data = np.load(PCEfilename, allow_pickle = True)
engf_meanErr_all_3 = data['meanErrList']
engf_varErr_all_3  = data['varErrList']
engf_runtime_all_3 = data['runtimeList']
engf_niters_all_3  = data['nitersList']
P_3 = data['P']

# engf filtered p = 3
PCEfilename = './results/engfres_p3_FILT_' + MCname
data = np.load(PCEfilename, allow_pickle = True)
engf_meanErr_filt_3 = data['meanErrList']
engf_varErr_filt_3  = data['varErrList']
engf_runtime_filt_3 = data['runtimeList']
engf_niters_filt_3  = data['nitersList']

# =============================================================================
# Load Monte Carlo results for partial Ns (unfiltered) from file
# =============================================================================
N_2 = NfactList * P_2
N_3 = NfactList * P_3

data = np.load('MCtrendlines_p2_70000.npz', allow_pickle = True)
fpafMeanErrAvg_2 = data['fpafMeanErrAvg'][:-1]
fpafMeanErrBest_2 = data['fpafMeanErrBest'][:-1]
fpafMeanErrWorst_2 = data['fpafMeanErrWorst'][:-1]
fpafVarErrAvg_2 = data['fpafVarErrAvg'][:-1]
fpafVarErrBest_2 = data['fpafVarErrBest'][:-1]
fpafVarErrWorst_2 = data['fpafVarErrWorst'][:-1]
engfMeanErrAvg_2 = data['engfMeanErrAvg'][:-1]
engfMeanErrBest_2 = data['engfMeanErrBest'][:-1]
engfMeanErrWorst_2 = data['engfMeanErrWorst'][:-1]
engfVarErrAvg_2 = data['engfVarErrAvg'][:-1]
engfVarErrBest_2 = data['engfVarErrBest'][:-1]
engfVarErrWorst_2 = data['engfVarErrWorst'][:-1]

# # =============================================================================
# # Get Monte Carlo results for partial Ns (unfiltered), 100 sample sets each pt.
# # =============================================================================
# N_2 = NfactList * P_2
# N_3 = NfactList * P_3

# MCfilename = './results/' + MCname
# data = np.load(MCfilename, allow_pickle = True)
# outsList = data['outsList']
# Ntot = len(outsList)

# fpaf_MC_all = np.asarray([out.fpaf for out in outsList])
# engf_MC_all = np.asarray([out.engf for out in outsList])

# fpaf_meanErr_MC_2 = []
# fpaf_varErr_MC_2 = []
# engf_meanErr_MC_2 = []
# engf_varErr_MC_2 = []

# nt = 10 # number of Monte Carlo trendlines to compute

# for Ni in N_2.astype(int):
#     fpaf_MCmeanList = []
#     fpaf_MCvarList = []
#     engf_MCmeanList = []
#     engf_MCvarList = []
    
#     for i in range(nt):
#         # generate Ni random indices in range of full Monte Carlo
#         inds = np.random.choice(Ntot, size = Ni, replace = False)
#         ss = stats.describe(fpaf_MC_all[inds])
#         fpaf_MCmeanList.append(ss[2])
#         fpaf_MCvarList.append(ss[3])
        
#         ss = stats.describe(engf_MC_all[inds])
#         engf_MCmeanList.append(ss[2])
#         engf_MCvarList.append(ss[3])
        
#     fpaf_meanErr_MC_2.append(np.mean(fpaf_MCmeanList) - fpaf_meanMC_all)
#     fpaf_varErr_MC_2.append(np.mean(fpaf_MCvarList) - fpaf_varMC_all)
#     engf_meanErr_MC_2.append(np.mean(engf_MCmeanList) - engf_meanMC_all)
#     engf_varErr_MC_2.append(np.mean(engf_MCvarList) - engf_varMC_all)
    
# fpaf_meanErr_MC_3 = []
# fpaf_varErr_MC_3 = []
# engf_meanErr_MC_3 = []
# engf_varErr_MC_3 = []

# for Ni in N_3.astype(int):
#     fpaf_MCmeanList = []
#     fpaf_MCvarList = []
#     engf_MCmeanList = []
#     engf_MCvarList = []
    
#     for i in range(nt):
#         # generate Ni random indices in range of full Monte Carlo
#         inds = np.random.choice(Ntot, size = Ni, replace = False)
#         ss = stats.describe(fpaf_MC_all[inds])
#         fpaf_MCmeanList.append(ss[2])
#         fpaf_MCvarList.append(ss[3])
        
#         ss = stats.describe(engf_MC_all[inds])
#         engf_MCmeanList.append(ss[2])
#         engf_MCvarList.append(ss[3])
        
#     fpaf_meanErr_MC_3.append(np.mean(fpaf_MCmeanList) - fpaf_meanMC_all)
#     fpaf_varErr_MC_3.append(np.mean(fpaf_MCvarList) - fpaf_varMC_all)
#     engf_meanErr_MC_3.append(np.mean(engf_MCmeanList) - engf_meanMC_all)
#     engf_varErr_MC_3.append(np.mean(engf_MCvarList) - engf_varMC_all)

    
# fpaf_meanErr_MC_2 = np.asarray(fpaf_meanErr_MC_2)
# fpaf_varErr_MC_2 = np.asarray(fpaf_varErr_MC_2)
# engf_meanErr_MC_2 = np.asarray(engf_meanErr_MC_2)
# engf_varErr_MC_2 = np.asarray(engf_varErr_MC_2)
# fpaf_meanErr_MC_3 = np.asarray(fpaf_meanErr_MC_3)
# fpaf_varErr_MC_3 = np.asarray(fpaf_varErr_MC_3)
# engf_meanErr_MC_3 = np.asarray(engf_meanErr_MC_3)
# engf_varErr_MC_3 = np.asarray(engf_varErr_MC_3)


# =============================================================================
# Plots
# =============================================================================
## FPAF MEAN p = 2
fig1, ax1 = plt.subplots(1,1)
ax1.plot(N_2, abs(fpaf_meanErr_all_2[0:5] / fpaf_meanMC_all) * 100, ':',
         label = 'unfiltered, sigma = {}'.format(sigmaList[0]))
ax1.plot(N_2, abs(fpaf_meanErr_all_2[5:10] / fpaf_meanMC_all) * 100, ':',
         label = 'unfiltered, sigma = {}'.format(sigmaList[1]))
ax1.plot(N_2, abs(fpaf_meanErr_all_2[10:15] / fpaf_meanMC_all) * 100, ':',
         label = 'unfiltered, sigma = {}'.format(sigmaList[2]))
ax1.plot(N_2, abs(fpaf_meanErr_filt_2[0:5] / fpaf_meanMC_all) * 100, ':',
         label = 'filtered, sigma = {}'.format(sigmaList[0]))
ax1.plot(N_2, abs(fpaf_meanErr_filt_2[5:10] / fpaf_meanMC_all) * 100, ':',
         label = 'filtered, sigma = {}'.format(sigmaList[1]))
ax1.plot(N_2, abs(fpaf_meanErr_filt_2[10:15] / fpaf_meanMC_all) * 100, ':',
         label = 'filtered, sigma = {}'.format(sigmaList[2]))
# ax1.plot(N_2, abs(fpaf_meanErr_MC_2 / fpaf_meanMC_all) * 100, 'r--', linewidth = 2,
#          label = 'Monte Carlo (unfiltered)')
ax1.plot(N_2, fpafMeanErrAvg_2, 'b--', linewidth = 2,
         label = 'Monte Carlo - average case')
ax1.plot(N_2, fpafMeanErrBest_2, 'g--', linewidth = 2,
         label = 'Monte Carlo - best case')
ax1.plot(N_2, fpafMeanErrWorst_2, 'r--', linewidth = 2,
         label = 'Monte Carlo - worst case')
ax1.set_yscale('log')
ax1.grid()
ax1.legend()
plt.xlabel('Number of trials')
plt.ylabel('Percent error against 70,000 Monte Carlo case')
plt.title('FPAF Mean Error, p = 2')

## FPAF VAR p = 2
fig2, ax2 = plt.subplots(1,1)
ax2.plot(N_2, abs(fpaf_varErr_all_2[0:5] / fpaf_varMC_all) * 100, ':',
          label = 'unfiltered, sigma = {}'.format(sigmaList[0]))
ax2.plot(N_2, abs(fpaf_varErr_all_2[5:10] / fpaf_varMC_all) * 100, ':',
          label = 'unfiltered, sigma = {}'.format(sigmaList[1]))
ax2.plot(N_2, abs(fpaf_varErr_all_2[10:15] / fpaf_varMC_all) * 100, ':',
          label = 'unfiltered, sigma = {}'.format(sigmaList[2]))
ax2.plot(N_2, abs(fpaf_varErr_filt_2[0:5] / fpaf_varMC_all) * 100, ':',
          label = 'filtered, sigma = {}'.format(sigmaList[0]))
ax2.plot(N_2, abs(fpaf_varErr_filt_2[5:10] / fpaf_varMC_all) * 100, ':',
          label = 'filtered, sigma = {}'.format(sigmaList[1]))
ax2.plot(N_2, abs(fpaf_varErr_filt_2[10:15] / fpaf_varMC_all) * 100, ':',
          label = 'filtered, sigma = {}'.format(sigmaList[2]))
# ax2.plot(N_2, abs(fpaf_varErr_MC_2 / fpaf_varMC_all) * 100, 'r--', linewidth = 2,
#           label = 'Monte Carlo (unfiltered)')
ax2.plot(N_2, fpafVarErrAvg_2, 'b--', linewidth = 2,
         label = 'Monte Carlo - average case')
ax2.plot(N_2, fpafVarErrBest_2, 'g--', linewidth = 2,
         label = 'Monte Carlo - best case')
ax2.plot(N_2, fpafVarErrWorst_2, 'r--', linewidth = 2,
         label = 'Monte Carlo - worst case')
ax2.set_yscale('log')
ax2.grid()
ax2.legend()
plt.xlabel('Number of trials')
plt.ylabel('Percent error magnitude against 70,000 Monte Carlo case')
plt.title('FPAF Variance Error, p = 2')

## ENGF MEAN p = 2
fig3, ax3 = plt.subplots(1,1)
ax3.plot(N_2, abs(engf_meanErr_all_2[0:5] / engf_meanMC_all) * 100, ':',
          label = 'unfiltered, sigma = {}'.format(sigmaList[0]))
ax3.plot(N_2, abs(engf_meanErr_all_2[5:10] / engf_meanMC_all) * 100, ':',
          label = 'unfiltered, sigma = {}'.format(sigmaList[1]))
ax3.plot(N_2, abs(engf_meanErr_all_2[10:15] / engf_meanMC_all) * 100, ':',
          label = 'unfiltered, sigma = {}'.format(sigmaList[2]))
ax3.plot(N_2, abs(engf_meanErr_filt_2[0:5] / engf_meanMC_all) * 100, ':',
          label = 'filtered, sigma = {}'.format(sigmaList[0]))
ax3.plot(N_2, abs(engf_meanErr_filt_2[5:10] / engf_meanMC_all) * 100, ':',
          label = 'filtered, sigma = {}'.format(sigmaList[1]))
ax3.plot(N_2, abs(engf_meanErr_filt_2[10:15] / engf_meanMC_all) * 100, ':',
          label = 'filtered, sigma = {}'.format(sigmaList[2]))
# ax3.plot(N_2, abs(engf_meanErr_MC_2 / engf_meanMC_all) * 100, 'r--', linewidth = 2,
#           label = 'Monte Carlo (unfiltered)')
ax3.plot(N_2, engfMeanErrAvg_2, 'b--', linewidth = 2,
         label = 'Monte Carlo - average case')
ax3.plot(N_2, engfMeanErrBest_2, 'g--', linewidth = 2,
         label = 'Monte Carlo - best case')
ax3.plot(N_2, engfMeanErrWorst_2, 'r--', linewidth = 2,
         label = 'Monte Carlo - worst case')
ax3.set_yscale('log')
ax3.grid()
ax3.legend()
plt.xlabel('Number of trials')
plt.ylabel('Percent error magnitude against 70,000 Monte Carlo case')
plt.title('ENGF Mean Error, p = 2')

## ENGF VAR p = 2
fig4, ax4 = plt.subplots(1,1)
ax4.plot(N_2, abs(engf_varErr_all_2[0:5] / engf_varMC_all) * 100, ':',
          label = 'unfiltered, sigma = {}'.format(sigmaList[0]))
ax4.plot(N_2, abs(engf_varErr_all_2[5:10] / engf_varMC_all) * 100, ':',
          label = 'unfiltered, sigma = {}'.format(sigmaList[1]))
ax4.plot(N_2, abs(engf_varErr_all_2[10:15] / engf_varMC_all) * 100, ':',
          label = 'unfiltered, sigma = {}'.format(sigmaList[2]))
ax4.plot(N_2, abs(engf_varErr_filt_2[0:5] / engf_varMC_all) * 100, ':',
          label = 'filtered, sigma = {}'.format(sigmaList[0]))
ax4.plot(N_2, abs(engf_varErr_filt_2[5:10] / engf_varMC_all) * 100, ':',
          label = 'filtered, sigma = {}'.format(sigmaList[1]))
ax4.plot(N_2, abs(engf_varErr_filt_2[10:15] / engf_varMC_all) * 100, ':',
          label = 'filtered, sigma = {}'.format(sigmaList[2]))
# ax4.plot(N_2, abs(engf_varErr_MC_2 / engf_varMC_all) * 100, 'r--', linewidth = 2,
#           label = 'Monte Carlo (unfiltered)')
ax4.plot(N_2, engfVarErrAvg_2, 'b--', linewidth = 2,
         label = 'Monte Carlo - average case')
ax4.plot(N_2, engfVarErrBest_2, 'g--', linewidth = 2,
         label = 'Monte Carlo - best case')
ax4.plot(N_2, engfVarErrWorst_2, 'r--', linewidth = 2,
         label = 'Monte Carlo - worst case')
ax4.set_yscale('log')
ax4.grid()
ax4.legend()
plt.xlabel('Number of trials')
plt.ylabel('Percent error magnitude against 70,000 Monte Carlo case')
plt.title('ENGF Variance Error, p = 2')

# ## FPAF MEAN p = 3
# fig5, ax5 = plt.subplots(1,1)
# ax5.plot(N_3, abs(fpaf_meanErr_all_3[0:5] / fpaf_meanMC_all) * 100, ':',
#          label = 'unfiltered, sigma = {}'.format(sigmaList[0]))
# ax5.plot(N_3, abs(fpaf_meanErr_all_3[5:10] / fpaf_meanMC_all) * 100, ':',
#          label = 'unfiltered, sigma = {}'.format(sigmaList[1]))
# ax5.plot(N_3, abs(fpaf_meanErr_all_3[10:15] / fpaf_meanMC_all) * 100, ':',
#          label = 'unfiltered, sigma = {}'.format(sigmaList[2]))
# ax5.plot(N_3, abs(fpaf_meanErr_filt_3[0:5] / fpaf_meanMC_all) * 100, ':',
#          label = 'filtered, sigma = {}'.format(sigmaList[0]))
# ax5.plot(N_3, abs(fpaf_meanErr_filt_3[5:10] / fpaf_meanMC_all) * 100, ':',
#          label = 'filtered, sigma = {}'.format(sigmaList[1]))
# ax5.plot(N_3, abs(fpaf_meanErr_filt_3[10:15] / fpaf_meanMC_all) * 100, ':',
#          label = 'filtered, sigma = {}'.format(sigmaList[2]))
# ax5.plot(N_3, abs(fpaf_meanErr_MC_3 / fpaf_meanMC_all) * 100, 'r--', linewidth = 2,
#          label = 'Monte Carlo (unfiltered)')
# ax5.grid()
# ax5.legend()
# plt.xlabel('Number of trials')
# plt.ylabel('Percent error magnitude against 70,000 Monte Carlo case')
# plt.title('FPAF Mean Error, p = 3')

# ## FPAF VAR p = 3
# fig6, ax6 = plt.subplots(1,1)
# ax6.plot(N_3, abs(fpaf_varErr_all_3[0:5] / fpaf_varMC_all) * 100, ':',
#          label = 'unfiltered, sigma = {}'.format(sigmaList[0]))
# ax6.plot(N_3, abs(fpaf_varErr_all_3[5:10] / fpaf_varMC_all) * 100, ':',
#          label = 'unfiltered, sigma = {}'.format(sigmaList[1]))
# ax6.plot(N_3, abs(fpaf_varErr_all_3[10:15] / fpaf_varMC_all) * 100, ':',
#          label = 'unfiltered, sigma = {}'.format(sigmaList[2]))
# ax6.plot(N_3, abs(fpaf_varErr_filt_3[0:5] / fpaf_varMC_all) * 100, ':',
#          label = 'filtered, sigma = {}'.format(sigmaList[0]))
# ax6.plot(N_3, abs(fpaf_varErr_filt_3[5:10] / fpaf_varMC_all) * 100, ':',
#          label = 'filtered, sigma = {}'.format(sigmaList[1]))
# ax6.plot(N_3, abs(fpaf_varErr_filt_3[10:15] / fpaf_varMC_all) * 100, ':',
#          label = 'filtered, sigma = {}'.format(sigmaList[2]))
# ax6.plot(N_3, abs(fpaf_varErr_MC_3 / fpaf_varMC_all) * 100, 'r--', linewidth = 2,
#          label = 'Monte Carlo (unfiltered)')
# ax6.grid()
# ax6.legend()
# plt.xlabel('Number of trials')
# plt.ylabel('Percent error magnitude against 70,000 Monte Carlo case')
# plt.title('FPAF Variance Error, p = 3')

# ## ENGF MEAN p = 3
# fig7, ax7 = plt.subplots(1,1)
# ax7.plot(N_3, abs(engf_meanErr_all_3[0:5] / engf_meanMC_all) * 100, ':',
#          label = 'unfiltered, sigma = {}'.format(sigmaList[0]))
# ax7.plot(N_3, abs(engf_meanErr_all_3[5:10] / engf_meanMC_all) * 100, ':',
#          label = 'unfiltered, sigma = {}'.format(sigmaList[1]))
# ax7.plot(N_3, abs(engf_meanErr_all_3[10:15] / engf_meanMC_all) * 100, ':',
#          label = 'unfiltered, sigma = {}'.format(sigmaList[2]))
# ax7.plot(N_3, abs(engf_meanErr_filt_3[0:5] / engf_meanMC_all) * 100, ':',
#          label = 'filtered, sigma = {}'.format(sigmaList[0]))
# ax7.plot(N_3, abs(engf_meanErr_filt_3[5:10] / engf_meanMC_all) * 100, ':',
#          label = 'filtered, sigma = {}'.format(sigmaList[1]))
# ax7.plot(N_3, abs(engf_meanErr_filt_3[10:15] / engf_meanMC_all) * 100, ':',
#          label = 'filtered, sigma = {}'.format(sigmaList[2]))
# ax7.plot(N_3, abs(engf_meanErr_MC_3 / engf_meanMC_all) * 100, 'r--', linewidth = 2,
#          label = 'Monte Carlo (unfiltered)')
# ax7.grid()
# ax7.legend()
# plt.xlabel('Number of trials')
# plt.ylabel('Percent error magnitude against 70,000 Monte Carlo case')
# plt.title('ENGF Mean Error, p = 3')

# ## ENGF VAR p = 3
# fig8, ax8 = plt.subplots(1,1)
# ax8.plot(N_3, abs(engf_varErr_all_3[0:5] / engf_varMC_all) * 100, ':',
#          label = 'unfiltered, sigma = {}'.format(sigmaList[0]))
# ax8.plot(N_3, abs(engf_varErr_all_3[5:10] / engf_varMC_all) * 100, ':',
#          label = 'unfiltered, sigma = {}'.format(sigmaList[1]))
# ax8.plot(N_3, abs(engf_varErr_all_3[10:15] / engf_varMC_all) * 100, ':',
#          label = 'unfiltered, sigma = {}'.format(sigmaList[2]))
# ax8.plot(N_3, abs(engf_varErr_filt_3[0:5] / engf_varMC_all) * 100, ':',
#          label = 'filtered, sigma = {}'.format(sigmaList[0]))
# ax8.plot(N_3, abs(engf_varErr_filt_3[5:10] / engf_varMC_all) * 100, ':',
#          label = 'filtered, sigma = {}'.format(sigmaList[1]))
# ax8.plot(N_3, abs(engf_varErr_filt_3[10:15] / engf_varMC_all) * 100, ':',
#          label = 'filtered, sigma = {}'.format(sigmaList[2]))
# ax8.plot(N_3, abs(engf_varErr_MC_3 / engf_varMC_all) * 100, 'r--', linewidth = 2,
#          label = 'Monte Carlo (unfiltered)')
# ax8.grid()
# ax8.legend()
# plt.xlabel('Number of trials')
# plt.ylabel('Percent error magnitude against 70,000 Monte Carlo case')
# plt.title('ENGF Variance Error, p = 3')







# # =============================================================================
# # Plots
# # =============================================================================
# ## FPAF MEAN
# fig1, ax1 = plt.subplots(1,1)
# ax1.plot(N, abs(fpaf_meanErr_all_2_01), '.--',
#          label = 'Unfiltered, sigma = 0.01')
# ax1.plot(N, abs(fpaf_meanErr_all_2_001), '.--',
#          label = 'Unfiltered, sigma = 0.001')
# ax1.plot(N, abs(fpaf_meanErr_all_2_0001), '.--',
#          label = 'Unfiltered, sigma = 0.0001')
# ax1.plot(N, abs(fpaf_meanErr_filt_2_01), '.--',
#          label = 'Filtered, sigma = 0.01')
# ax1.plot(N, abs(fpaf_meanErr_filt_2_001), '.--',
#          label = 'Filtered, sigma = 0.001')
# ax1.plot(N, abs(fpaf_meanErr_filt_2_0001), '.--',
#          label = 'Filtered, sigma - 0.0001')
# ax1.plot(N, abs(fpaf_meanErr_MC), '.:',
#          label = 'Monte Carlo')
# ax1.grid()
# ax1.legend()
# plt.xlabel('Number of trials')
# plt.ylabel('Percent error magnitude against 60,000 Monte Carlo case')
# plt.title('FPAF Mean Error')


# ## FPAF VAR
# fig2, ax2 = plt.subplots(1,1)
# ax2.plot(N, abs(fpaf_varErr_all_2_01), '.--',
#          label = 'Unfiltered, sigma = 0.01')
# ax2.plot(N, abs(fpaf_varErr_all_2_001), '.--',
#          label = 'Unfiltered, sigma = 0.001')
# ax2.plot(N, abs(fpaf_varErr_all_2_0001), '.--',
#          label = 'Unfiltered, sigma = 0.0001')
# ax2.plot(N, abs(fpaf_varErr_filt_2_01), '.--',
#          label = 'Filtered, sigma = 0.01')
# ax2.plot(N, abs(fpaf_varErr_filt_2_001), '.--',
#          label = 'Filtered, sigma = 0.001')
# ax2.plot(N, abs(fpaf_varErr_filt_2_0001), '.--',
#          label = 'Filtered, sigma - 0.0001')
# ax2.plot(N, abs(fpaf_varErr_MC), '.:',
#          label = 'Monte Carlo')
# ax2.grid()
# ax2.legend()
# plt.xlabel('Number of trials')
# plt.ylabel('Percent error magnitude against 60,000 Monte Carlo case')
# plt.title('FPAF Variance Error')


# ## ENGF MEAN
# fig3, ax3 = plt.subplots(1,1)
# ax3.plot(N, abs(engf_meanErr_all_2_01), '.--',
#          label = 'Unfiltered, sigma = 0.01')
# ax3.plot(N, abs(engf_meanErr_all_2_001), '.--',
#          label = 'Unfiltered, sigma = 0.001')
# ax3.plot(N, abs(engf_meanErr_all_2_0001), '.--',
#          label = 'Unfiltered, sigma = 0.0001')
# ax3.plot(N, abs(engf_meanErr_filt_2_01), '.--',
#          label = 'Filtered, sigma = 0.01')
# ax3.plot(N, abs(engf_meanErr_filt_2_001), '.--',
#          label = 'Filtered, sigma = 0.001')
# ax3.plot(N, abs(engf_meanErr_filt_2_0001), '.--',
#          label = 'Filtered, sigma - 0.0001')
# ax3.plot(N, abs(engf_meanErr_MC), '.:',
#          label = 'Monte Carlo')
# ax3.grid()
# ax3.legend()
# plt.xlabel('Number of trials')
# plt.ylabel('Percent error magnitude against 60,000 Monte Carlo case')
# plt.title('ENGF Mean Error')


# ## ENGF VAR
# fig4, ax4 = plt.subplots(1,1)
# ax4.plot(N, abs(engf_varErr_all_2_01), '.--',
#          label = 'Unfiltered, sigma = 0.01')
# ax4.plot(N, abs(engf_varErr_all_2_001), '.--',
#          label = 'Unfiltered, sigma = 0.001')
# ax4.plot(N, abs(engf_varErr_all_2_0001), '.--',
#          label = 'Unfiltered, sigma = 0.0001')
# ax4.plot(N, abs(engf_varErr_filt_2_01), '.--',
#          label = 'Filtered, sigma = 0.01')
# ax4.plot(N, abs(engf_varErr_filt_2_001), '.--',
#          label = 'Filtered, sigma = 0.001')
# ax4.plot(N, abs(engf_varErr_filt_2_0001), '.--',
#          label = 'Filtered, sigma - 0.0001')
# ax4.plot(N, abs(engf_varErr_MC), '.:',
#          label = 'Monte Carlo')
# ax4.grid()
# ax4.legend()
# plt.xlabel('Number of trials')
# plt.ylabel('Percent error magnitude against 60,000 Monte Carlo case')
# plt.title('ENGF Variance Error')





















