# -*- coding: utf-8 -*-
"""
SciTechPlots.py:
    plots compiled PCE data for SciTech 2021 paper
Created on Tue Nov 17 16:26:56 2020

@author: Samuel Albert
"""

import numpy as np
from scipy import stats
import matplotlib.pyplot as plt
import matplotlib as mpl


plt.close('all')
mpl.rcParams["figure.autolayout"] = True
mpl.rcParams.update({'font.size': 16})

# =============================================================================
# Load MC data
# =============================================================================
MCname = 'Mars_60000_1022201133.npz'

# hard-code sigma and Nfact lists
sigmaList = [0.01, 0.001, 0.0001]
NfactList = np.array([0.2, 0.5, 0.9, 1.1, 4])

# =============================================================================
# Input compilation data (from Google Sheets document)
#   naming convention: QOI_mean/var_all/filt_p_sigma
# =============================================================================
# N = N_SPG / P (same for all lists below)
P_2 = 2415
N_2 = NfactList * P_2


## FPAF unfiltered
fpaf_mean_all_2_01 = np.array([8.40554478, 8.419652716, 8.404467019,
                               8.39921254, 8.406985406])
fpaf_var_all_2_01 = np.array([0.5013498414, 0.5339571899, 0.7265565189,
                              2.434263774, 0.8067869536])
fpaf_mean_all_2_001 = np.array([8.407821277, 8.419579338, 8.412704258,
                                8.399212223, 8.406985406])
fpaf_var_all_2_001 = np.array([0.5241022444, 0.5550029678, 1.225703281,
                               2.434261993, 0.8067869531])
fpaf_mean_all_2_0001 = np.array([8.407489335, 8.418432861, 8.414135301,
                                 8.399212392, 8.406985405])
fpaf_var_all_2_0001 = np.array([0.5244383605, 0.5582453235, 1.35261509,
                                2.434262256, 0.8067869609])

## FPAF filtered
fpaf_mean_filt_2_01 = np.array([8.40554478, 8.419652716, 8.419891812,
                                8.421890827, 8.425419559])
fpaf_var_filt_2_01 = np.array([0.5013498414, 0.5339571899, 0.5696962683,
                               0.5760721124, 0.6099267138])
fpaf_mean_filt_2_001 = np.array([8.407821277, 8.419579338, 8.421262549,
                                 8.431840116, 8.42541956])
fpaf_var_filt_2_001 = np.array([0.5241022444, 0.5550029678, 0.6669506266,
                                0.7238611884, 0.6099267156])
fpaf_mean_filt_2_0001 = np.array([8.407489335, 8.418432861, 8.420288777,
                                  8.431840141, 8.425419559])
fpaf_var_filt_2_0001 = np.array([0.5244383605, 0.5582453235, 0.718924618,
                                 0.7238611619, 0.609926714])

## ENGF unfiltered
engf_mean_all_2_01 = np.array([-1.657402306, -1.667863887, -1.676346795,
                               -1.6817111, -1.677584163])
engf_var_all_2_01 = np.array([1.605149465, 1.638323649, 1.680892757,
                              1.948518577, 1.720564389])
engf_mean_all_2_001 = np.array([-1.658105564, -1.669946084, -1.675073166,
                                -1.681711134, -1.677584163])
engf_var_all_2_001 = np.array([1.612665145, 1.645937618, 1.743228843,
                               1.948518133, 1.720564388])
engf_mean_all_2_0001 = np.array([-1.658195154, -1.669935346, -1.674192582,
                                 -1.681711077, -1.677584163])
engf_var_all_2_0001 = np.array([1.613157557, 1.646305012, 1.756624473,
                                1.948518366, 1.720564391])

## ENGF filtered
engf_mean_filt_2_01 = np.array([-1.657402306, -1.667863887, -1.669460124,
                                -1.669571773, -1.670278878])
engf_var_filt_2_01 = np.array([1.605149465, 1.638323649, 1.645238106,
                               1.645977621, 1.653430171])
engf_mean_filt_2_001 = np.array([-1.658105564, -1.669946084, -1.670405717,
                                 -1.669485697, -1.670659536])
engf_var_filt_2_001 = np.array([1.612665145, 1.645937618, 1.652228506,
                                1.655632229, 1.658348952])
engf_mean_filt_2_0001 = np.array([-1.658195154, -1.669935346, -1.670121347,
                                  -1.6694857, -1.670659536])
engf_var_filt_2_0001 = np.array([1.613157557, 1.646305012, 1.654903856,
                                 1.655632238, 1.658348953])

# =============================================================================
# Load Monte Carlo results for comparison from file
# =============================================================================
data = np.load('MCtrendlines_p2_60000.npz', allow_pickle = True)
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

# =============================================================================
# Load data for comparison
# =============================================================================
MCfilename = './results/Mars_60000_1022201133.npz'
data = np.load(MCfilename, allow_pickle = True)

outsList = data['outsList']

# unfiltered statistics
fpaf_MC_all = np.asarray([out.fpaf for out in outsList])
fpaf_stats_MC_all = stats.describe(fpaf_MC_all)
fpaf_mean_MC_all = fpaf_stats_MC_all[2]
fpaf_var_MC_all = fpaf_stats_MC_all[3]


engf_MC_all = np.asarray([out.engf for out in outsList])
engf_stats_MC_all = stats.describe(engf_MC_all)
engf_mean_MC_all = engf_stats_MC_all[2]
engf_var_MC_all = engf_stats_MC_all[3]

# filtered statistics
fpaf_MC_filt = fpaf_MC_all[fpaf_MC_all >= 0]
fpaf_stats_MC_filt = stats.describe(fpaf_MC_filt)
fpaf_mean_MC_filt = fpaf_stats_MC_filt[2]
fpaf_var_MC_filt = fpaf_stats_MC_filt[3]

engf_MC_filt = engf_MC_all[fpaf_MC_all >= 0]
engf_stats_MC_filt = stats.describe(engf_MC_filt)
engf_mean_MC_filt = engf_stats_MC_filt[2]
engf_var_MC_filt = engf_stats_MC_filt[3]

# =============================================================================
# Get error compared to full 60,000 MC "truth" value
# =============================================================================
## FPAF unfiltered
fpaf_meanErr_all_2_01 = fpaf_mean_all_2_01 - fpaf_mean_MC_all
fpaf_meanErr_all_2_001 = fpaf_mean_all_2_001 - fpaf_mean_MC_all
fpaf_meanErr_all_2_0001 = fpaf_mean_all_2_0001 - fpaf_mean_MC_all
fpaf_varErr_all_2_01 = fpaf_var_all_2_01 - fpaf_var_MC_all
fpaf_varErr_all_2_001 = fpaf_var_all_2_001 - fpaf_var_MC_all
fpaf_varErr_all_2_0001 = fpaf_var_all_2_0001 - fpaf_var_MC_all

## FPAF filtered
fpaf_meanErr_filt_2_01 = fpaf_mean_filt_2_01 - fpaf_mean_MC_filt
fpaf_meanErr_filt_2_001 = fpaf_mean_filt_2_001 - fpaf_mean_MC_filt
fpaf_meanErr_filt_2_0001 = fpaf_mean_filt_2_0001 - fpaf_mean_MC_filt
fpaf_varErr_filt_2_01 = fpaf_var_filt_2_01 - fpaf_var_MC_filt
fpaf_varErr_filt_2_001 = fpaf_var_filt_2_001 - fpaf_var_MC_filt
fpaf_varErr_filt_2_0001 = fpaf_var_filt_2_0001 - fpaf_var_MC_filt

## ENGF unfiltered
engf_meanErr_all_2_01 = engf_mean_all_2_01 - engf_mean_MC_all
engf_meanErr_all_2_001 = engf_mean_all_2_001 - engf_mean_MC_all
engf_meanErr_all_2_0001 = engf_mean_all_2_0001 - engf_mean_MC_all
engf_varErr_all_2_01 = engf_var_all_2_01 - engf_var_MC_all
engf_varErr_all_2_001 = engf_var_all_2_001 - engf_var_MC_all
engf_varErr_all_2_0001 = engf_var_all_2_0001 - engf_var_MC_all

## ENGF filtered
engf_meanErr_filt_2_01 = engf_mean_filt_2_01 - engf_mean_MC_filt
engf_meanErr_filt_2_001 = engf_mean_filt_2_001 - engf_mean_MC_filt
engf_meanErr_filt_2_0001 = engf_mean_filt_2_0001 - engf_mean_MC_filt
engf_varErr_filt_2_01 = engf_var_filt_2_01 - engf_var_MC_filt
engf_varErr_filt_2_001 = engf_var_filt_2_001 - engf_var_MC_filt
engf_varErr_filt_2_0001 = engf_var_filt_2_0001 - engf_var_MC_filt



# =============================================================================
# Plots
# =============================================================================
## FPAF MEAN
fig1, ax1 = plt.subplots(1,1, figsize=(8,8))
# ax1.plot(N_2, abs(fpaf_meanErr_all_2_01 / fpaf_mean_MC_all) * 100, 'C0:.',
#           label = 'Unfiltered, sigma = 0.01', linewidth = 2)
ax1.plot(N_2, abs(fpaf_meanErr_all_2_001 / fpaf_mean_MC_all) * 100, 'C1:.',
         label = 'Unfiltered, sigma = 0.001', linewidth = 2)
ax1.plot(N_2, abs(fpaf_meanErr_all_2_0001 / fpaf_mean_MC_all) * 100, 'C2:.',
         label = 'Unfiltered, sigma = 0.0001', linewidth = 2)
# ax1.plot(N_2, abs(fpaf_meanErr_filt_2_01 / fpaf_mean_MC_filt) * 100, 'C0-..',
#           label = 'Filtered, sigma = 0.01', linewidth = 2)
ax1.plot(N_2, abs(fpaf_meanErr_filt_2_001 / fpaf_mean_MC_filt) * 100, 'C1-..',
         label = 'Filtered, sigma = 0.001', linewidth = 2)
ax1.plot(N_2, abs(fpaf_meanErr_filt_2_0001 / fpaf_mean_MC_filt) * 100, 'C2-..',
         label = 'Filtered, sigma - 0.0001', linewidth = 2)
ax1.plot(N_2, fpafMeanErrAvg_2, 'C0--.', linewidth = 2,
         label = 'Monte Carlo - average case')
ax1.plot(N_2, fpafMeanErrBest_2, 'C1--.', linewidth = 2,
         label = 'Monte Carlo - best case')
# ax1.plot(N_2, fpafMeanErrWorst_2, 'C2--.', linewidth = 2,
#           label = 'Monte Carlo - worst case')
ax1.set_yscale('log')
ax1.grid()
ax1.legend()
plt.xlabel(r'Sample size $N$')
plt.ylabel('Percent error against 70,000 Monte Carlo case')
# plt.title('FPAF Mean Error, p = 2')

# extra x-axis ticks
axsub1 = ax1.twiny()
axsub1.set_xlim(ax1.get_xlim())
axsub1.set_xticks(N_2)
axsub1.set_xticklabels(NfactList)
axsub1.set_xlabel(r'Sample size : number of coefficients $N/(P+1)$')

## FPAF VAR
fig2, ax2 = plt.subplots(1,1, figsize=(8,8))
ax2.plot(N_2, abs(fpaf_varErr_all_2_01 / fpaf_var_MC_all) * 100, 'C0:.',
         label = 'Unfiltered, sigma = 0.01', linewidth = 2)
ax2.plot(N_2, abs(fpaf_varErr_all_2_001 / fpaf_var_MC_all) * 100, 'C1:.',
         label = 'Unfiltered, sigma = 0.001', linewidth = 2)
# ax2.plot(N_2, abs(fpaf_varErr_all_2_0001 / fpaf_var_MC_all) * 100, 'C2:.',
#           label = 'Unfiltered, sigma = 0.0001', linewidth = 2)
ax2.plot(N_2, abs(fpaf_varErr_filt_2_01 / fpaf_var_MC_filt) * 100, 'C0-..',
          label = 'Filtered, sigma = 0.01', linewidth = 2)
ax2.plot(N_2, abs(fpaf_varErr_filt_2_001 / fpaf_var_MC_filt) * 100, 'C1-..',
         label = 'Filtered, sigma = 0.001', linewidth = 2)
# ax2.plot(N_2, abs(fpaf_varErr_filt_2_0001 / fpaf_var_MC_filt) * 100, 'C2-..',
#          label = 'Filtered, sigma - 0.0001', linewidth = 2)
ax2.plot(N_2, fpafVarErrAvg_2, 'C0--.', linewidth = 2,
         label = 'Monte Carlo - average case')
ax2.plot(N_2, fpafVarErrBest_2, 'C1--.', linewidth = 2,
         label = 'Monte Carlo - best case')
# ax2.plot(N_2, fpafVarErrWorst_2, 'C2--.', linewidth = 2,
#           label = 'Monte Carlo - worst case')
ax2.set_yscale('log')
ax2.grid()
ax2.legend()
plt.xlabel(r'Sample size $N$')
plt.ylabel('Percent error against 70,000 Monte Carlo case')
# plt.title('FPAF Var Error, p = 2')

# extra x-axis ticks
axsub2 = ax2.twiny()
axsub2.set_xlim(ax2.get_xlim())
axsub2.set_xticks(N_2)
axsub2.set_xticklabels(NfactList)
axsub2.set_xlabel(r'Sample size : number of coefficients $N/(P+1)$')

## ENGF MEAN
fig3, ax3 = plt.subplots(1,1, figsize=(8,8))
# ax3.plot(N_2, abs(engf_meanErr_all_2_01 / engf_mean_MC_all) * 100, 'C0:.',
#           label = 'Unfiltered, sigma = 0.01', linewidth = 2)
ax3.plot(N_2, abs(engf_meanErr_all_2_001 / engf_mean_MC_all) * 100, 'C1:.',
         label = 'Unfiltered, sigma = 0.001', linewidth = 2)
ax3.plot(N_2, abs(engf_meanErr_all_2_0001 / engf_mean_MC_all) * 100, 'C2:.',
         label = 'Unfiltered, sigma = 0.0001', linewidth = 2)
ax3.plot(N_2, abs(engf_meanErr_filt_2_01 / engf_mean_MC_filt) * 100, 'C0-..',
          label = 'Filtered, sigma = 0.01', linewidth = 2)
ax3.plot(N_2, abs(engf_meanErr_filt_2_001 / engf_mean_MC_filt) * 100, 'C1-..',
         label = 'Filtered, sigma = 0.001', linewidth = 2)
# ax3.plot(N_2, abs(engf_meanErr_filt_2_0001 / engf_mean_MC_filt) * 100, 'C2-..',
#          label = 'Filtered, sigma - 0.0001', linewidth = 2)
ax3.plot(N_2, engfMeanErrAvg_2, 'C0--.', linewidth = 2,
         label = 'Monte Carlo - average case')
ax3.plot(N_2, engfMeanErrBest_2, 'C1--.', linewidth = 2,
         label = 'Monte Carlo - best case')
# ax3.plot(N_2, engfMeanErrWorst_2, 'C2--.', linewidth = 2,
#           label = 'Monte Carlo - worst case')
ax3.set_yscale('log')
ax3.grid()
ax3.legend()
plt.xlabel(r'Sample size $N$')
plt.ylabel('Percent error against 70,000 Monte Carlo case')
# plt.title('ENGF Mean Error, p = 2')

# extra x-axis ticks
axsub3 = ax3.twiny()
axsub3.set_xlim(ax3.get_xlim())
axsub3.set_xticks(N_2)
axsub3.set_xticklabels(NfactList)
axsub3.set_xlabel(r'Sample size : number of coefficients $N/(P+1)$')

## ENGF VAR
fig4, ax4 = plt.subplots(1,1, figsize=(8,8))
# ax4.plot(N_2, abs(engf_varErr_all_2_01 / engf_var_MC_all) * 100, 'C0:.',
#           label = 'Unfiltered, sigma = 0.01', linewidth = 2)
ax4.plot(N_2, abs(engf_varErr_all_2_001 / engf_var_MC_all) * 100, 'C1:.',
         label = 'Unfiltered, sigma = 0.001', linewidth = 2)
ax4.plot(N_2, abs(engf_varErr_all_2_0001 / engf_var_MC_all) * 100, 'C2:.',
         label = 'Unfiltered, sigma = 0.0001', linewidth = 2)
ax4.plot(N_2, abs(engf_varErr_filt_2_01 / engf_var_MC_filt) * 100, 'C0-..',
          label = 'Filtered, sigma = 0.01', linewidth = 2)
ax4.plot(N_2, abs(engf_varErr_filt_2_001 / engf_var_MC_filt) * 100, 'C1-..',
         label = 'Filtered, sigma = 0.001', linewidth = 2)
# ax4.plot(N_2, abs(engf_varErr_filt_2_0001 / engf_var_MC_filt) * 100, 'C2-..',
#          label = 'Filtered, sigma - 0.0001', linewidth = 2)
ax4.plot(N_2, engfVarErrAvg_2, 'C0--.', linewidth = 2,
         label = 'Monte Carlo - average case')
ax4.plot(N_2, engfVarErrBest_2, 'C1--.', linewidth = 2,
         label = 'Monte Carlo - best case')
# ax4.plot(N_2, engfVarErrWorst_2, 'C2--.', linewidth = 2,
#           label = 'Monte Carlo - worst case')
ax4.set_yscale('log')
ax4.grid()
ax4.legend()
plt.xlabel(r'Sample size $N$')
plt.ylabel('Percent error against 70,000 Monte Carlo case')
# plt.title('ENGF Var Error, p = 2')

# extra x-axis ticks
axsub4 = ax4.twiny()
axsub4.set_xlim(ax4.get_xlim())
axsub4.set_xticks(N_2)
axsub4.set_xticklabels(NfactList)
axsub4.set_xlabel(r'Sample size : number of coefficients $N/(P+1)$')



