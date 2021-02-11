# -*- coding: utf-8 -*-
"""
Created on Mon Feb  8 15:00:25 2021

AeroDrop_dispersed.py:
    Simulates orbiter and probe using FNP(A/E)G guidance under dispersions.

@author: Samuel Albert
"""

from conversions import LLAEHV2RV, RV2LLAEHV, Vinf2VN, getApsesSphPR
from sim import Params
from guidance import updateFNPAG, dynFNPAGPhase1Sph
import ODEs
import planetaryconstants as constants

import numpy as np
from scipy.integrate import solve_ivp
import sys
import matplotlib as mpl
import matplotlib.pyplot as plt
import time

tic = time.time()

plt.close('all')
mpl.rcParams["figure.autolayout"] = True
mpl.rcParams.update({'font.size': 15})

# =============================================================================
# Main FNPAG Function
# =============================================================================
def doFNPAG(mode, paramsTrue, paramsNom, t0, xx0vec, sig0, sigd, ts,
            verbose = True, plotsOn = True, updatesOn = True):
    '''
    main function for FNPAG with variable mode (supports modes 1 and 2).
    paramsTrue holds real values used in propagation, paramsNom holds nominal
        values used in prediction for guidance updates.
    Uses spherical planet-relative EOMs.
    '''
    
    # =========================================================================
    # Set up events
    # =========================================================================
    event1 = lambda t, y: ODEs.above_max_alt_sph(t, y, paramsTrue)
    event1.terminal = True
    event1.direction = 1
    event2 = lambda t, y: ODEs.below_min_alt_sph(t, y, paramsTrue)
    event2.terminal = True
    event3 = lambda t, y: ODEs.switchEvent(t, y, ts)
    event3.terminal = True
    
    # =========================================================================
    # Main integration loop
    # =========================================================================
    # # make sure correct bank angles are set
    # paramsTrue.bank = sig0
    # paramsNom.bank = sig0
    
    # initialize
    tvec = np.arange(t0, paramsNom.tf + paramsTrue.dtGdn, paramsTrue.dtGdn)
    tvecEval = np.empty(1)
    tvecEval[0] = t0
    
    xxvec = np.empty((len(xx0vec), 1))
    xxvec[:] = np.NaN
    xxvec[:,0] = xx0vec
    
    tsList = []
    sigdList = []
    
    t = t0
    ind = 0
    
    xx0veci = xx0vec
    
    ## PHASE 1 ##
    phase = 1
    
    for ind, t in enumerate(tvec):
        # update guidance
        if updatesOn:
            ts = updateFNPAG(xx0veci, t, ts, sig0, sigd, phase, mode,
                             paramsNom, sphericalEOMs = True)
            # update event with new ts value
            event3 = lambda t, y: ODEs.switchEvent(t, y, ts)
            event3.terminal = True
            
        if verbose:
            print('PHASE 1: updating guidance at time {0:.3f}    ts = {1:.3f}'\
                  .format(t, ts))
        tsList.append(ts)
        
        ### TROUBLESHOOTING: stop for profiler
        if t > 2:
            sys.exit('stopped execution for profiler')
        
        # propagate until next guidance update or switching time
        tspan = (t, t + paramsTrue.dtGdn)
        soli = solve_ivp(lambda t, y: ODEs.sphericalEntryEOMs(t, y,
                                                              np.radians(sig0),
                                                              paramsTrue),
                         tspan, xx0veci,
                         rtol = paramsTrue.rtol, atol = paramsTrue.atol,
                         events = (event1, event2, event3))
        
        # append results
        xxvec = np.append(xxvec, soli.y, axis = 1)
        tvecEval = np.append(tvecEval, soli.t)
        xx0veci = xxvec[:,-1]
        
        # either proceed to next guidance update or enter phase 2
        if soli.status == 0:
            # reached next guidance update
            continue
        elif len(soli.t_events[2]) > 0:
            # reached switching time
            break
        
        else:
            sys.exit('propagator never reached next guidance update or'\
                     ' switching time in phase 1')
        
    tvecP1 = tvecEval * 1
    
    ## create new tvec for phase 2
    tvec = np.arange(tvecEval[-1], paramsNom.tf + paramsTrue.dtGdn, paramsTrue.dtGdn)
    
    ## PHASE 2 ##
    phase = 2
    # paramsTrue.bank = sigd
    # paramsNom.bank = sigd
    print()
    
    for ind, t in enumerate(tvec):
        # update guidance
        if updatesOn:
            sigd = updateFNPAG(xxvec[:,-1], t, ts, sig0, sigd, phase, mode,
                               paramsNom)
            # paramsTrue.bank = sigd
            # paramsNom.bank = sigd
            
        
        # make sure sigd is in [-180, 180] deg range
        sigd = (sigd + 180) % (360) - 180
        
        if verbose:
            print('PHASE 2: updating guidance at time {0:.3f}'\
                  '    sigd = {1:.3f} deg'.format(t, sigd))
        sigdList.append(sigd)
        
        # propagate until next guidance update or final state
        tspan = (t, t + paramsTrue.dtGdn)
        soli = solve_ivp(lambda t, y: ODEs.sphericalEntryEOMs(t, y,
                                                              np.radians(sigd),
                                                              paramsTrue),
                         tspan, xx0veci,
                         rtol = paramsTrue.rtol, atol = paramsTrue.atol,
                         events = (event1, event2))
        
        # append results
        xxvec = np.append(xxvec, soli.y, axis = 1)
        tvecEval = np.append(tvecEval, soli.t)
        xx0veci = xxvec[:,-1]
        
        if soli.status == 1:
            # reached terminal state
            break
        
    # trim off first array elements
    xxvec = xxvec[:,1:]
    tvecEval = tvecEval[1:]
    
    tvecP2 = tvecEval * 1
    
    if plotsOn:
        # fig = plt.figure()
        # ax = fig.add_subplot(111)
        ax.plot(xxvec[3,:], xxvec[0,:], '--', label = 'true trajectory')
        ax.legend()
    
    # =========================================================================
    # Compute apoapsis error and DV magnitudes
    # =========================================================================
    raf, rpf = getApsesSphPR(xxvec[:,-1], paramsTrue)
    raErr = raf - paramsTrue.raStar
    
    mu = paramsTrue.p.mu * 1e9
    
    DV1 = np.sqrt(2 * mu)\
        * abs((np.sqrt(1/paramsTrue.raStar - 1\
                       / (paramsTrue.raStar + paramsTrue.rpStar))\
               - np.sqrt(1/paramsTrue.raStar - 1 / (raf + rpf))))
    DV2 = np.sqrt(2 * mu)\
        * abs((np.sqrt(1/paramsTrue.rpStar - 1\
                       / (paramsTrue.raStar + paramsTrue.rpStar))\
               - np.sqrt(1/paramsTrue.rpStar - 1 / (raf + paramsTrue.rpStar))))
    DV = DV1 + DV2
    
    print('\n\nfinal apoapsis error: {0:.3e} m'.format(raErr))
    print('delta-V for periapsis raise: {0:.3f} m/s'.format(DV1))
    print('delta-V for apoapsis correction: {0:.3f} m/s'.format(DV2))
    print('total delta-V: {0:.3f} m/s'.format(DV))
    
    return xxvec, tvecEval, raf, rpf, raErr, DV,\
        tsList, sigdList, tvecP1, tvecP2
    

# =============================================================================
# Script
# =============================================================================
### CREATE params INPUT CLASS FOR NOMINAL VALUES
paramsNom = Params()
paramsNom.p = constants.MARS
paramsNom.returnTimeVectors = True

### INPUT ATM TABLE - GET ATM TABLE FROM GRAM DATA FILE
paramsNom.dMode = 'table'
filename = '../data/dens_Mars_nom.txt'
atmdata = np.genfromtxt(filename, names=True)
paramsNom.atmdat = np.array([atmdata['Var_X'], atmdata['DENSAV']])
    
### VEHICLE PARAMS
paramsNom.m = 2920 # kg, roughly MSL mass
paramsNom.CD = 1.6 # roughly MSL CD
paramsNom.LD = 0.25
paramsNom.BC = 130

paramsNom.A = paramsNom.m / (paramsNom.BC * paramsNom.CD)
paramsNom.CL = paramsNom.CD * paramsNom.LD
paramsNom.Rn = np.sqrt(paramsNom.A / np.pi) / 2

### WIND-RELATIVE INITIAL STATE
paramsNom.lat = 18.38
paramsNom.lon = -77.58
paramsNom.alt = paramsNom.p.halt
paramsNom.efpaWR = -12
# paramsNom.efpaWR = -11
paramsNom.hdaWR = 0
paramsNom.vmagWR = 6

xx0vec = np.array([(paramsNom.alt + paramsNom.p.rad) * 1e3,
                   np.radians(paramsNom.lon),
                   np.radians(paramsNom.lat),
                   paramsNom.vmagWR * 1e3,
                   np.radians(paramsNom.efpaWR),
                   np.radians(paramsNom.hdaWR + 90)])
    
### NOMINAL SIMULATION PARAMS ###
t0 = 0
# xx0vec = np.block([paramsNom.x0, paramsNom.v0])
sig0 = 15
sigd = 150
ts = 157

# search brackets for Brent's Method
paramsNom.sig1 = 15
paramsNom.sig2 = 180
paramsNom.ts1 = 100
paramsNom.ts2 = 200

# other settings and constants
paramsNom.rtol = 1e-10
paramsNom.atol = 1e-10
paramsNom.errtol1 = 0
paramsNom.errtol2 = 0
paramsNom.dtGdn = 1 # s
paramsNom.hmin = 10
paramsNom.hmax = paramsNom.p.halt + 1e7
paramsNom.tf = 6000
paramsNom.raStar = (250 + paramsNom.p.rad) * 1e3
paramsNom.rpStar = (250 + paramsNom.p.rad) * 1e3

### SET TRUE MC SIMULATION PARAMS ###
paramsTrue = paramsNom

# add dispersions here as desired

# =============================================================================
# Demonstrate dynamics
# =============================================================================
xxvecs1 = dynFNPAGPhase1Sph(xx0vec, 0, ts, sig0, sigd, paramsNom)
raf, rpf = getApsesSphPR(xxvecs1[:,-1], paramsNom)
print(raf/1e3 - paramsNom.p.rad)

h = np.linalg.norm(xxvecs1[:3,:], axis = 0)
vmag = np.linalg.norm(xxvecs1[3:,:], axis = 0)

fig = plt.figure()
ax = fig.add_subplot(111)
# ax.plot(vmag, h, label = 'predicted trajectory')
ax.grid()

# =============================================================================
# Call FNPAG
# =============================================================================
mode = 1
xxvec, tvecEval, raf, rpf, raErr, DV,\
        tsList, sigdList, tvecP1, tvecP2 = doFNPAG(mode, paramsTrue, paramsNom,
                                                    t0, xx0vec, sig0, sigd, ts)



    

toc = time.time()

print('Total time elapsed: {0:.2f} s'.format(toc-tic))




















