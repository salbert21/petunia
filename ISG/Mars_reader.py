# -*- coding: utf-8 -*-
"""
Created on Mon Oct  5 13:46:25 2020

Mars_reader.py:
    script to load data from .npz file and get relevant QoIs

@author: Samuel Albert
"""

import numpy as np

from sim import Params, Outs
import constants


filename = './results/Mars_5000_1005170208.npz'

data = np.load(filename, allow_pickle=True)
paramsList = data['paramsList']
outsList = data['outsList']
CD_LB = data['CD_LB'].item()
CD_UB = data['CD_UB'].item()
m_LB = data['m_LB'].item()
m_UB = data['m_UB'].item()
efpanom = data['efpanom'].item()
efpa_sig = data['efpa_sig'].item()
vmagnom = data['vmagnom'].item()
vmag_sig = data['vmag_sig'].item()
paramsf = data['paramsf'].item()


# =============================================================================
# Create list for each dispersed input variable
# =============================================================================
CDs = [param.CD for param in paramsList]
ms = [param.m for param in paramsList]
BCs = [param.BC for param in paramsList]
efpas = [param.efpaWR for param in paramsList]
vmags = [param.vmagWR for param in paramsList]
atmprofiles = [param.atmdat for param in paramsList]


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
minAltList = [min(np.linalg.norm(out.rvec_N, axis=0) - paramsf.p.rad)\
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











