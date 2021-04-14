# -*- coding: utf-8 -*-
"""
marvel.py:
    contains script and functions for SHIELD+MARVEL concept
Created on Tue Apr 13 12:08:55 2021

@author: Samuel Albert
"""
import pdb

import numpy as np
import time
import matplotlib.pyplot as plt
import sys
import scipy.interpolate as interp
from scipy.integrate import solve_ivp
import copy

import planetaryconstants as constants
from conversions import LLAEHV2RV, Vinf2VN
import ODEs
from sim import Params, Outs, mainAD
from atm import getRho_from_table

tic = time.time()
plt.close('all')

# =============================================================================
# Main function
# =============================================================================
def marvelProbe(params):
    '''
    main function for MARVEL probes.
    1) back-propagate probes for params.tback seconds
    2) split probe into N separate probes
    3) apply delta-V of params.DV to each probe in hard-coded directions
    4) forward propagate to entry interface
    5) call mainAD for each probe from entry interface to landing
    '''
    
    # get inertial cartesian state from given planet-relative spherical state
    params.x0, params.vInfvec_N = LLAEHV2RV(params.lat, params.lon,
                                            params.alt, params.efpaWR,
                                            params.hdaWR, params.vmagWR,
                                            params, params.tspan[0])
    params.v0 = Vinf2VN(params.x0, params.vInfvec_N, params, params.tspan[0])
    
    # propagate vehicle backwards in time by params.tback seconds
    yy0 = np.block([params.x0, params.v0])
    
    solBack = solve_ivp(lambda t, y: ODEs.TBP(t, y, params),
                        (params.tback, 0), yy0.flatten(), rtol = params.rtol,
                        atol = params.atol)
    
    xx0vecCenter = solBack.y[:,-1]
    x0vecCenter = xx0vecCenter[:3]
    v0vecCenter = xx0vecCenter[3:]
    
    # split probe into 4 separate probes
    rDir = x0vecCenter / np.linalg.norm(x0vecCenter)
    vDir = v0vecCenter / np.linalg.norm(v0vecCenter)
    hDir = np.cross(rDir, vDir)
    hDir = hDir / np.linalg.norm(hDir)
    
    thDir = np.cross(hDir, rDir)
    thDir = thDir / np.linalg.norm(thDir)
    
    phi = np.radians(0)
    eta = np.radians(0)
    
    aDir = np.cos(phi) * thDir + np.sin(phi) * hDir
    bDir = -aDir
    cDir = np.cos(eta) * hDir - np.sin(eta) * thDir
    dDir = -cDir
    
    
    # pdb.set_trace()
    
    # print(np.degrees(np.arccos(np.dot(-rDir, vDir))))
    
    v0vecDown = v0vecCenter + params.DV * aDir
    v0vecUp = v0vecCenter + params.DV * bDir
    
    v0vecLeft = v0vecCenter + params.DV * cDir
    v0vecRight = v0vecCenter + params.DV * dDir
    
    
    xx0vecDown = np.hstack([x0vecCenter, v0vecDown])
    xx0vecUp = np.hstack([x0vecCenter, v0vecUp])
    xx0vecLeft = np.hstack([x0vecCenter, v0vecLeft])
    xx0vecRight = np.hstack([x0vecCenter, v0vecRight])
    
    # pdb.set_trace()
    
    # propagate each probe forward until atm interface
    eventEI = lambda t, y: ODEs.above_max_alt(t, y, params)
    eventEI.terminal = True
    
    paramsDown = copy.deepcopy(params)
    solDown = solve_ivp(lambda t, y: ODEs.TBP(t, y, params),
                        (0, 10 * params.tback), xx0vecDown, rtol = params.rtol,
                        atol = params.atol, events = (eventEI))
    paramsDown.x0 = solDown.y[:3,-1]
    paramsDown.v0 = solDown.y[3:,-1]
    paramsDown.inputType = 'inertial vectors'
    
    paramsUp = copy.deepcopy(params)
    solUp = solve_ivp(lambda t, y: ODEs.TBP(t, y, params),
                        (0, 10 * params.tback), xx0vecUp, rtol = params.rtol,
                        atol = params.atol, events = (eventEI))
    paramsUp.x0 = solUp.y[:3,-1]
    paramsUp.v0 = solUp.y[3:,-1]
    paramsUp.inputType = 'inertial vectors'
    
    paramsLeft = copy.deepcopy(params)
    solLeft = solve_ivp(lambda t, y: ODEs.TBP(t, y, params),
                        (0, 10 * params.tback), xx0vecLeft, rtol = params.rtol,
                        atol = params.atol, events = (eventEI))
    paramsLeft.x0 = solLeft.y[:3,-1]
    paramsLeft.v0 = solLeft.y[3:,-1]
    paramsLeft.inputType = 'inertial vectors'
    
    paramsRight = copy.deepcopy(params)
    solRight = solve_ivp(lambda t, y: ODEs.TBP(t, y, params),
                        (0, 10 * params.tback), xx0vecRight, rtol = params.rtol,
                        atol = params.atol, events = (eventEI))
    paramsRight.x0 = solRight.y[:3,-1]
    paramsRight.v0 = solRight.y[3:,-1]
    paramsRight.inputType = 'inertial vectors'
    
    
    # propagate all 4 probes to the surface
    outsDown = Outs()
    outsUp = Outs()
    outsLeft = Outs()
    outsRight = Outs()
    
    outsDown = mainAD(paramsDown, params.tspan, params.events, outsDown)
    outsUp = mainAD(paramsUp, params.tspan, params.events, outsUp)
    outsLeft = mainAD(paramsLeft, params.tspan, params.events, outsLeft)
    outsRight = mainAD(paramsRight, params.tspan, params.events, outsRight)
    
    # pdb.set_trace()
    
    print(outsLeft.lon0)
    print(outsLeft.lat0)
    
    
    
    
    return outsDown, outsUp, outsLeft, outsRight





# =============================================================================
# Set up probe params
# =============================================================================
### CREATE params INPUT CLASS
params = Params()
params.p = constants.MARS
params.returnTimeVectors = True
params.atmMod = 'nom'

### INPUT ATM TABLE - GET ATM TABLE FROM GRAM DATA FILE
params.dMode = 'table'
filename = '../data/dens_Mars_nom.txt'
atmdata = np.genfromtxt(filename, names=True)
params.atmdat = np.array([atmdata['Var_X'], atmdata['DENSAV']])

### VEHICLE PARAMS (NOT CHANGED DURING GRID SEARCH)
params.m = 2920 # kg, roughly MSL mass
params.CD = 1.6 # roughly MSL CD
params.BC = 10 # kg/m^2
params.LD = 0
params.bank = 0 # value doesn't matter since L/D = 0

### WIND-RELATIVE INITIAL STATE (COMPONENTS NOT CHANGED DURING GRID SEARCH)
params.inputType = 'wind-relative angles'
params.lat = 18.38
params.lon = -77.58
params.alt = params.p.halt
params.efpaWR = -12
params.hdaWR = 0
params.vmagWR = 6

### SET SEPARATION DELTA-V AND TIMING
params.DV = 1e-3 # km/s
params.tback = 1 * 60*60 # s

### TIME VECTOR AND EXIT CONDITIONS
params.tspan = (params.tback, 30000)

# EXIT CONDITIONS:
params.hmin = 10
params.hmax = params.p.halt + 1e-7

event1 = lambda t, y: ODEs.below_min_alt(t, y, params)
event1.terminal = True
event2 = lambda t, y: ODEs.above_max_alt(t, y, params)
event2.terminal = True
event2.direction = 1
events = (event1, event2)
params.events = events

### INTEGRATION TOLERANCE
params.rtol = 1e-11
params.atol = 1e-11



# =============================================================================
# Call main MARVEL probe function
# =============================================================================
# outsDown, outsUp, outsLeft, outsRight = marvelProbe(params)


# =============================================================================
# Plots
# =============================================================================
fig = plt.figure()
ax = fig.add_subplot(111)

# fig2 = plt.figure()
# ax2 = fig2.add_subplot(111)

# ax.plot(outsUp.lonf, outsUp.latf, 'o', color = 'C0', label = 'Up')
# ax.plot(outsDown.lonf, outsDown.latf, 'o', color = 'C1', label = 'Down')
# ax.plot(outsLeft.lonf, outsLeft.latf, 'o', color = 'C2', label = 'Left')
# ax.plot(outsRight.lonf, outsRight.latf, 'o', color = 'C3', label = 'Right')



# print('\napproximate up-down distance:')
# print(np.linalg.norm(outsDown.rvec_Nf - outsUp.rvec_Nf))

# print('approximate left-right distance:')
# print(np.linalg.norm(outsLeft.rvec_Nf - outsRight.rvec_Nf))


DVList = [1, 2, 3, 4, 5, 6]
tbackList = [0.5, 1, 1.5, 2, 2.5, 3]
# tbackList = [2, 4, 6, 8, 10]
# tbackList = [3]

# for i, DV in enumerate(DVList):
#     params.DV = DV * 1e-3

for i, tback in enumerate(tbackList):
    params.tback = tback * 3600
    
    outsDown, outsUp, outsLeft, outsRight = marvelProbe(params)
    
    ax.plot(outsUp.lonf, outsUp.latf, 'o', color = 'C0', label = 'Up')
    ax.plot(outsDown.lonf, outsDown.latf, 'o', color = 'C1', label = 'Down')
    ax.plot(outsLeft.lonf, outsLeft.latf, 'o', color = 'C2', label = 'Left')
    ax.plot(outsRight.lonf, outsRight.latf, 'o', color = 'C3', label = 'Right')
    
    # h = np.linalg.norm(outsLeft.rvec_N, axis = 0) - params.p.rad #/1e3
    # v = np.linalg.norm(outsLeft.vvec_N, axis = 0)
    # ax2.plot(v, h)
    
    if i == 0:
        ax.legend()
        
print('\napproximate up-down distance:')
print(np.linalg.norm(outsDown.rvec_Nf - outsUp.rvec_Nf))

print('approximate left-right distance:')
print(np.linalg.norm(outsLeft.rvec_Nf - outsRight.rvec_Nf))





ax.set_xlabel('longitude, deg')
ax.set_ylabel('latitude, deg')
# ax.legend()
ax.grid()
ax.axis('equal')


# ax2.grid()

