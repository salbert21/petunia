# -*- coding: utf-8 -*-
"""
reader_KLE.py:
    reader to extract data from .npz generated from Monte Carlo trials via KLE
Created on Wed Oct 21 12:13:13 2020

@author: Samuel Albert
"""


import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl

# from sim import Params, Outs
import constants
###

plt.close('all')
mpl.rcParams["figure.autolayout"] = True
mpl.rcParams.update({'font.size': 16})

# =============================================================================
# Import raw data
# =============================================================================

# filename = './results/Mars_5000_1021122652.npz'
# filename = './results/Mars_60000_1022201133.npz'
filename = './results/Mars_70000_1118142637.npz'

data = np.load(filename, allow_pickle=True)
paramsList = data['paramsList']
outsList = data['outsList']
m_YList = data['m_YList']
CD_YList = data['CD_YList']
vmag_YList = data['vmag_YList']
efpa_YList = data['efpa_YList']
atm_YsList = data['atm_YsList']
CDmean = data['CDmean'].item()
CDstd = data['CDstd'].item()
mmean = data['mmean'].item()
mstd = data['mstd'].item()
efpamean = data['efpamean'].item()
efpastd = data['efpastd'].item()
vmagmean = data['vmagmean'].item()
vmagstd = data['vmagstd'].item()


# =============================================================================
# Create list for each dispersed input variable
# =============================================================================
CDs = [param.CD for param in paramsList]
ms = [param.m for param in paramsList]
BCs = [param.BC for param in paramsList]
efpas = [param.efpaWR for param in paramsList]
vmags = [param.vmagWR for param in paramsList]


# =============================================================================
# Create list for each quantity of interest (QoI)
# =============================================================================
rafList = [out.raf for out in outsList] # apoapsis radius at tf, km
engfList = [out.engf for out in outsList] # specific energy at tf, km^2/s^2
fpafList = [out.fpaf for out in outsList] # flight path angle at tf, deg
rvecfList = [out.rvec_N for out in outsList] # final inertial pos. vec., km
vvecfList = [out.vvec_N for out in outsList] # final inertial vel. vec., km/s
qpeakList = [out.qpeak for out in outsList] # peak heat rate, W/cm^2
QloadList = [out.Qload for out in outsList] # total heat load, J/cm^2
gloadList = [out.gload for out in outsList] # peak deceleration, Earth g's
minAltList = [min(np.linalg.norm(out.rvec_N, axis=0) - constants.EARTH.rad)\
              for out in outsList] # minimum altitude, km
tofList = [out.tvec[-1] - out.tvec[0] for out in outsList] # time of flight, s


# =============================================================================
# Quick analytics
# =============================================================================
Nmc = len(ms)

# number of runs that impacted the surface
NCrash = sum(fpaf < 0 for fpaf in fpafList)

# number of runs that escaped (skip-out trajectories)
NEscape = sum(engf > 0 for engf in engfList)

NCapture = sum(engf < 0 and fpaf > 0 for engf, fpaf in zip(engfList, fpafList))

# can check that the sume of the above counts equals total MC count

print('{} cases crashed, {} cases escaped, and {} cases aerocaptured'\
      .format(NCrash, NEscape, NCapture))
    
# =============================================================================
# Plot to show bifurcation of data
# =============================================================================
fpafList = np.asarray(fpafList)
engfList = np.asarray(engfList)

fpafCaptured = fpafList[np.logical_and(engfList < 0, fpafList > 0)]
fpafEscaped = fpafList[engfList > 0]
fpafImpacted = fpafList[fpafList < 0]

engfCaptured = engfList[np.logical_and(engfList < 0, fpafList > 0)]
engfEscaped = engfList[engfList > 0]
engfImpacted = engfList[fpafList < 0]

fig, ax = plt.subplots(1,1)
ax.plot(engfCaptured, fpafCaptured, '.', label = 'Captured')
ax.plot(engfEscaped, fpafEscaped, 'C2.', label = 'Escaped')
ax.plot(engfImpacted, fpafImpacted, 'C1.', label = 'Impacted')

ax.grid()
ax.legend()
ax.set_xlabel(r'Final specific orbital energy, $km^2s^{-2}$')
ax.set_ylabel('Final flight path angle, deg')





















