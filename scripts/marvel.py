# -*- coding: utf-8 -*-
"""
marvel.py:
    contains script and functions for SHIELD+MARVEL concept
Created on Tue Apr 13 12:08:55 2021

@author: Samuel Albert
"""

from scipy.integrate import solve_ivp
import matplotlib.pyplot as plt
import matplotlib as mpl
import numpy as np
import copy
import sys
import pdb

from conversions import LLAEHV2RV, Vinf2VN, greatCircleDistDeg
import planetaryconstants as constants
import ODEs
from sim import Params, Outs, mainAD

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
    
    v0vecDown = v0vecCenter + params.DV * aDir
    v0vecUp = v0vecCenter + params.DV * bDir
    
    v0vecLeft = v0vecCenter + params.DV * cDir
    v0vecRight = v0vecCenter + params.DV * dDir
    
    
    xx0vecDown = np.hstack([x0vecCenter, v0vecDown])
    xx0vecUp = np.hstack([x0vecCenter, v0vecUp])
    xx0vecLeft = np.hstack([x0vecCenter, v0vecLeft])
    xx0vecRight = np.hstack([x0vecCenter, v0vecRight])
    
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
    
    outsDown = mainAD(paramsDown, params.tspan, params.events, outsDown,
                      verbose = False)
    outsUp = mainAD(paramsUp, params.tspan, params.events, outsUp,
                    verbose = False)
    outsLeft = mainAD(paramsLeft, params.tspan, params.events, outsLeft,
                      verbose = False)
    outsRight = mainAD(paramsRight, params.tspan, params.events, outsRight,
                       verbose = False)
    
    
    # make sure that all probes enter atm at least 1 km away from orbiter
    dista = greatCircleDistDeg(params.lon, params.lat,
                               outsDown.lon0, outsDown.lat0, params)
    distb = greatCircleDistDeg(params.lon, params.lat,
                               outsUp.lon0, outsUp.lat0, params)
    distc = greatCircleDistDeg(params.lon, params.lat,
                               outsLeft.lon0, outsLeft.lat0, params)
    distd = greatCircleDistDeg(params.lon, params.lat,
                               outsRight.lon0, outsRight.lat0, params)
    
    if min(dista, distb, distc, distd) < 1:
        sys.exit('a probe entered the atmosphere within {0:.3f} km of the'\
                 ' orbtier'.format(min(dista, distb, distc, distd)))
    
    # find minimum distance between two landing locations
    distAB = greatCircleDistDeg(outsDown.lonf, outsDown.latf,
                                outsUp.lonf, outsUp.latf, params)
    distAC = greatCircleDistDeg(outsDown.lonf, outsDown.latf,
                                outsLeft.lonf, outsLeft.latf, params)
    distAD = greatCircleDistDeg(outsDown.lonf, outsDown.latf,
                                outsRight.lonf, outsRight.latf, params)
    distBC = greatCircleDistDeg(outsUp.lonf, outsUp.latf,
                                outsLeft.lonf, outsLeft.latf, params)
    distBD = greatCircleDistDeg(outsUp.lonf, outsUp.latf,
                                outsRight.lonf, outsRight.latf, params)
    distCD = greatCircleDistDeg(outsLeft.lonf, outsLeft.latf,
                                outsRight.lonf, outsRight.latf, params)
    minDist = min(distAB, distAC, distAD, distBC, distBD, distCD)
    maxDist = max(distAB, distAC, distAD, distBC, distBD, distCD)
    
    return outsDown, outsUp, outsLeft, outsRight, (minDist, maxDist)


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

### SET DEFAULT SEPARATION DELTA-V AND TIMING
DVnom = 1e-3 # km/s
tbacknom = 1 * 60*60 # s

### TIME VECTOR AND EXIT CONDITIONS
params.tspan = (tbacknom, 30000)

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
# Vary Separation Delta-V
# =============================================================================
print('\n\nVARYING SEPARATION DELTA-V:')
fig1 = plt.figure()
ax1 = fig1.add_subplot(111)
ax1.set_xlabel('longitude, deg')
ax1.set_ylabel('latitude, deg')
ax1.set_title('Landing locations for varying $\Delta$V, separation at'\
              ' E-{0:.0f} minutes'.format(tbacknom/60))
ax1.grid()
cmap = mpl.cm.get_cmap('viridis')

DVList = [1, 2, 3, 4, 5, 6]

cList = np.linspace(1, 0, len(DVList))
cList = cmap(cList)
legendElements = [None] * len(cList)

for i, DV in enumerate(DVList):
    params.DV = DV / 1e3
    params.tback = tbacknom
    
    outsDown, outsUp, outsLeft, outsRight, dists = marvelProbe(params)
    
    ax1.plot(outsUp.lonf, outsUp.latf, 'o', color = cList[i])
    ax1.plot(outsDown.lonf, outsDown.latf, 'o', color = cList[i])
    ax1.plot(outsLeft.lonf, outsLeft.latf, 's', color = cList[i])
    ax1.plot(outsRight.lonf, outsRight.latf, 's', color = cList[i])
    
    legendElements[i] = mpl.patches.Patch(facecolor = cList[i],
                                          label = r'$\Delta$V = {0:.0f} m/s'\
                                              .format(DV))
    
    print('\nDV = {0:.2f} m/s, separation at E-{1:.0f} minutes:'\
          .format(params.DV*1e3, params.tback/60))
    print('Separation (min, max) of ({0:.3f}, {1:.3f}) km'\
          .format(dists[0], dists[1]))

cirEl = mpl.lines.Line2D([0], [0], marker = 'o', color = 'w',
                         label = '+/- along-track $\Delta$V', markersize = 10,
                         markerfacecolor = 'k')
legendElements.append(cirEl)
sqEl = mpl.lines.Line2D([0], [0], marker = 's', color = 'w',
                         label = '+/- cross-track $\Delta$V', markersize = 10,
                         markerfacecolor = 'k')
legendElements.append(sqEl)
ax1.legend(handles = legendElements)

# =============================================================================
# Vary Separation Timing
# =============================================================================
print('\n\nVARYING SEPARATION TIMING:')
fig2 = plt.figure()
ax2 = fig2.add_subplot(111)
ax2.set_xlabel('longitude, deg')
ax2.set_ylabel('latitude, deg')
ax2.set_title('Landing locations for varying separation time, $\Delta$V = '\
              ' {0:.0f} m/s'.format(DVnom*1e3))
ax2.grid()

tbackList = [0.5, 1, 1.5, 2, 2.5, 3]

cList = np.linspace(1, 0, len(tbackList))
cList = cmap(cList)
legendElements = [None] * len(cList)

for i, tback in enumerate(tbackList):
    params.DV = DVnom
    params.tback = tback * 60 * 60
    
    outsDown, outsUp, outsLeft, outsRight, dists = marvelProbe(params)
    
    ax2.plot(outsUp.lonf, outsUp.latf, 'o', color = cList[i])
    ax2.plot(outsDown.lonf, outsDown.latf, 'o', color = cList[i])
    ax2.plot(outsLeft.lonf, outsLeft.latf, 's', color = cList[i])
    ax2.plot(outsRight.lonf, outsRight.latf, 's', color = cList[i])
    
    legendElements[i] = mpl.patches.Patch(facecolor = cList[i],
                                          label = 'E-{0:.0f} min.'\
                                              .format(tback * 60))
    
    print('\nDV = {0:.2f} m/s, separation at E-{1:.0f} minutes:'\
          .format(params.DV*1e3, params.tback/60))
    print('Separation (min, max) of ({0:.3f}, {1:.3f}) km'\
          .format(dists[0], dists[1]))

legendElements.append(cirEl)
legendElements.append(sqEl)
ax2.legend(handles = legendElements)





