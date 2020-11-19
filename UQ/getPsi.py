# -*- coding: utf-8 -*-
"""
getPsi.py:
    build and save Psi matrix for aerocapture PCE solve
Created on Thu Oct 29 13:35:26 2020

@author: Samuel Albert
"""


import numpy as np

from UQ import nD_poly_array, piset

# =============================================================================
# Build Psi and save into binary file
# =============================================================================

# MCfilename = './results/Mars_60000_1022201133.npz'
# MCfilename = './results/Mars_10000_highvar_1116164345.npz'
MCfileshort = 'Mars_70000_1118142637.npz'
MCfilename = './results/' + MCfileshort

# set order p
p = 3

filtered = False

if filtered:
    Psifilename = 'Psi_p' + str(p) + '_FILT_' + MCfileshort
else:
    Psifilename = 'Psi_p' + str(p) + '_' + MCfileshort

# load MC data
data = np.load(MCfilename, allow_pickle = True)

# paramsList = data['paramsList']
# outsList = data['outsList']
m_YList = data['m_YList']
CD_YList = data['CD_YList']
vmag_YList = data['vmag_YList']
efpa_YList = data['efpa_YList']
atm_YsList = data['atm_YsList']
# CDmean = data['CDmean'].item()
# CDstd = data['CDstd'].item()
# mmean = data['mmean'].item()
# mstd = data['mstd'].item()
# efpamean = data['efpamean'].item()
# efpastd = data['efpastd'].item()
# vmagmean = data['vmagmean'].item()
# vmagstd = data['vmagstd'].item()

if filtered:
    # =========================================================================
    # Optional: filter-out impact cases
    # =========================================================================
    outsList = data['outsList']
    fpafList = np.asarray([out.fpaf for out in outsList])
    
    m_YList = m_YList[fpafList >= 0]
    CD_YList = CD_YList[fpafList >= 0]
    vmag_YList = vmag_YList[fpafList >= 0]
    efpa_YList = efpa_YList[fpafList >= 0]
    atm_YsList = atm_YsList[fpafList >= 0]
    # =========================================================================
    # =========================================================================


y_samps = np.block([m_YList[:,None], CD_YList[:,None], vmag_YList[:,None],
                    efpa_YList[:,None], atm_YsList])

# set total order and dimensions
d = y_samps.shape[1] # # dispersed inputs + dimensions in KLE

# get multi-indices, compute max number of samples we will use
index_pc = nD_poly_array(d,p)
P = index_pc.shape[0] # actually P+1
Ndat = y_samps.shape[0]

# now build Psi row-by-row for all rows
Psi = np.empty([Ndat, P])
Psi[:] = np.NaN

for i in range(Ndat):
    Psi[i,:] = piset(y_samps[i,:], index_pc, distribution = 'Gaussian')
    print(i)
    
np.savez(Psifilename, Psi = Psi)





















