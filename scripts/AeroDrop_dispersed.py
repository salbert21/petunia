# -*- coding: utf-8 -*-
"""
Created on Mon Feb  8 15:00:25 2021

AeroDrop_dispersed.py:
    Simulates orbiter and probe using FNP(A/E)G guidance under dispersions.

@author: Samuel Albert
"""

from conversions import getApsesSphPR
from sim import Params
from guidance import updateFNPAG
import ODEs
import planetaryconstants as constants
from atm import getMarsGRAMDensTableAll

import numpy as np
from scipy.integrate import solve_ivp
from scipy.stats import norm, uniform
from datetime import datetime
import sys
import matplotlib as mpl
import matplotlib.pyplot as plt
import time
import copy

tic = time.time()
datestring = datetime.now().strftime('%m%d%H%M%S')
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
                
    xswitchvec = xxvec[:,-1] # save state at switching time
        
    tvecP1 = tvecEval * 1
    
    ## create new tvec for phase 2
    tvec = np.arange(tvecEval[-1], paramsNom.tf + paramsTrue.dtGdn, paramsTrue.dtGdn)
    
    ## PHASE 2 ##
    phase = 2
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
        # uses global variable ax for plotting!
        ax.plot(xxvec[3,:]/1e3, xxvec[0,:]/1e3 - paramsNom.p.rad)
        ax.plot(xswitchvec[3]/1e3, xswitchvec[0]/1e3 - paramsNom.p.rad, 'o',
                markersize = 6, color=plt.gca().lines[-1].get_color())
    
    # =========================================================================
    # Compute apoapsis error and DV magnitudes
    # =========================================================================
    raf, rpf, engf = getApsesSphPR(xxvec[:,-1], paramsTrue, returnEng = True)
    raErr = raf - paramsTrue.raStar
    
    mu = paramsTrue.p.mu * 1e9
    
    # if hyperbolic orbit, assign NaN to both DV values
    if raf < 0:
        DV1 = np.NaN
        DV2 = np.NaN
    else:
        DV1 = np.sqrt(2 * mu)\
            * abs((np.sqrt(1/raf - 1\
                           / (raf + paramsTrue.rpStar))\
                   - np.sqrt(1/raf - 1 / (raf + rpf))))
        DV2 = np.sqrt(2 * mu)\
            * abs((np.sqrt(1/paramsTrue.rpStar - 1\
                           / (paramsTrue.raStar + paramsTrue.rpStar))\
                   - np.sqrt(1/paramsTrue.rpStar - 1\
                             / (raf + paramsTrue.rpStar))))
                
    DV = DV1 + DV2 # m/s
    
    print('\n\nfinal apoapsis error: {0:.3e} m'.format(raErr))
    print('delta-V for periapsis raise: {0:.3f} m/s'.format(DV1))
    print('delta-V for apoapsis correction: {0:.3f} m/s'.format(DV2))
    print('total delta-V: {0:.3f} m/s'.format(DV))
    
    return xxvec, tvecEval, raf, rpf, engf, raErr, DV,\
        tsList, sigdList, tvecP1, tvecP2, xswitchvec
    

# =============================================================================
# Orbiter Setup
# =============================================================================
### CREATE params INPUT CLASS FOR NOMINAL VALUES
paramsNom_O = Params()
paramsNom_O.p = constants.MARS
paramsNom_O.returnTimeVectors = True

### INPUT ATM TABLE - GET ATM TABLE FROM GRAM DATA FILE
paramsNom_O.dMode = 'table'
filename = '../data/dens_Mars_nom.txt'
atmdata = np.genfromtxt(filename, names=True)
paramsNom_O.atmdat = np.array([atmdata['Var_X'], atmdata['DENSAV']])
    
### VEHICLE PARAMS
paramsNom_O.m = 2920 # kg, roughly MSL mass
paramsNom_O.CD = 1.6 # roughly MSL CD
paramsNom_O.LD = 0.25
paramsNom_O.BC = 130

paramsNom_O.A = paramsNom_O.m / (paramsNom_O.BC * paramsNom_O.CD)
paramsNom_O.CL = paramsNom_O.CD * paramsNom_O.LD
paramsNom_O.Rn = np.sqrt(paramsNom_O.A / np.pi) / 2

### WIND-RELATIVE INITIAL STATE
paramsNom_O.lat = 18.38
paramsNom_O.lon = -77.58
paramsNom_O.alt = paramsNom_O.p.halt
paramsNom_O.efpaWR = -12
paramsNom_O.hdaWR = 0
paramsNom_O.vmagWR = 6

xx0vec = np.array([(paramsNom_O.alt + paramsNom_O.p.rad) * 1e3,
                   np.radians(paramsNom_O.lon),
                   np.radians(paramsNom_O.lat),
                   paramsNom_O.vmagWR * 1e3,
                   np.radians(paramsNom_O.efpaWR),
                   np.radians(paramsNom_O.hdaWR + 90)])
    
### NOMINAL SIMULATION PARAMS ###
t0 = 0
sig0 = 15
sigd = 150
ts = 158.043
# ts = 160

# search brackets for Brent's Method
paramsNom_O.sig1 = 0
paramsNom_O.sig2 = 180
paramsNom_O.ts1 = 100
paramsNom_O.ts2 = 300

# other settings and constants
paramsNom_O.rtol = 1e-10
paramsNom_O.atol = 1e-10
paramsNom_O.errtol1 = 0
paramsNom_O.errtol2 = 0
paramsNom_O.dtGdn = 1 # s
paramsNom_O.hmin = 10
paramsNom_O.hmax = paramsNom_O.p.halt + 1e-7
paramsNom_O.tf = 6000
paramsNom_O.raStar = (250 + paramsNom_O.p.rad) * 1e3
paramsNom_O.rpStar = (250 + paramsNom_O.p.rad) * 1e3

# =============================================================================
# Demonstrate dynamics for nominal scenario
# =============================================================================
mode = 1
xxvec, tvecEval, raf, rpf, engf, raErr, DV,\
        tsList, sigdList, tvecP1, tvecP2, xswitchvec = doFNPAG(mode,
                                                               paramsNom_O,
                                                               paramsNom_O, t0,
                                                               xx0vec, sig0,
                                                               sigd, ts,
                                                               plotsOn = False,
                                                               updatesOn = False)

fig = plt.figure()
ax = fig.add_subplot(111)
ax.plot(xxvec[3,:]/1e3, xxvec[0,:]/1e3 - paramsNom_O.p.rad, 'k', linewidth = 4,
        label = 'nominal trajectory')
ax.plot(xswitchvec[3]/1e3, xswitchvec[0]/1e3 - paramsNom_O.p.rad, 'ko',
        markersize = 10)
ax.grid()
ax.legend()
ax.set_xlabel('planet-relative velocity, km/s')
ax.set_ylabel('spherical altitude, km')
ax.set_title('Orbiter trajectories')

# =============================================================================
# Orbiter Dispersions Setup
# =============================================================================
# set number of MC runs
Nmc = 2

# copy, DO NOT ASSIGN, paramsTrue_O for dispersions. leave paramsNom_O alone.
paramsTrue_O = copy.deepcopy(paramsNom_O)

# load full 5000-case Monte Carlo atmsophere dataset from MarsGRAM-2010
filename = '../data/Mars_0.1_5000.txt'
densAll, densMean, h = getMarsGRAMDensTableAll(filename, Nmc)

# define dispersions for other parameters
BCMean_O = paramsNom_O.BC
# BCLB_O = 0.9 * paramsNom_O.BC # kg/m^2
# BCRng_O = 0.2 * paramsNom_O.BC # kg/m^2
BCLB_O = 0.95 * paramsNom_O.BC # kg/m^2
BCRng_O = 0.1 * paramsNom_O.BC # kg/m^2

L_DMean_O = paramsNom_O.LD
# L_DLB_O = 0.9 * paramsNom_O.LD
# L_DRng_O = 0.2 * paramsNom_O.LD
L_DLB_O = 0.95 * paramsNom_O.LD
L_DRng_O = 0.1 * paramsNom_O.LD

gam0Mean = np.radians(paramsNom_O.efpaWR)
gam0STD = np.radians(0.2) / 3

v0Mean = paramsNom_O.vmagWR * 1e3
v0STD = 10/3 # m/s

# initialize lists for actual dispersed parameters and results
xxvecList_O = []
tvecList_O = []
sigdvecList_O = []
tsvecList_O = []

raErrList_O = []
DVList_O = []
tsfList_O = []
sigdfList_O = []
rafList_O = []
rpfList_O = []
engfList_O = []

BCList_O = []
L_DList_O = []
gam0List = []
v0List = []

# output filename
outname = '../results/AeroDrop_dispersed_' + str(Nmc) + '_' + datestring

# =============================================================================
#  Monte Carlo loop
# =============================================================================
for i in range(Nmc):
    # initialize random variables
    BC_O = uniform.rvs(size = 1, loc = BCLB_O, scale = BCRng_O)[0]
    L_D_O = uniform.rvs(size = 1, loc = L_DLB_O, scale = L_DRng_O)[0]
    gam0 = norm.rvs(size = 1, loc = gam0Mean, scale = gam0STD)[0]
    v0 = norm.rvs(size = 1, loc = v0Mean, scale = v0STD)[0]
    
    paramsTrue_O.BC = BC_O
    paramsTrue_O.LD = L_D_O
    xx0vec[3] = v0
    xx0vec[4] = gam0
    
    # get GRAM atm profile
    paramsTrue_O.atmdat = np.array([h,densAll[:,i]])
    paramsTrue_O.atmdat = paramsTrue_O.atmdat[:,paramsTrue_O.atmdat[0,:].argsort()]
    
    # append inputs to lists (skip GRAM profiles, can always get those later)
    BCList_O.append(BC_O)
    L_DList_O.append(L_D_O)
    gam0List.append(gam0)
    v0List.append(v0)
    
    # run FNPAG simulation
    xxvec, tvecEval, raf, rpf, engf, raErr, DV, tsList, sigdList, tvecP1,\
        tvecP2, xswwitchvec = doFNPAG(mode, paramsTrue_O, paramsNom_O, t0,
                              xx0vec, sig0, sigd, ts)
    # xxvec, tvecEval, raf, rpf, engf, raErr, DV, tsList, sigdList, tvecP1,\
    #     tvecP2, xswwitchvec = doFNPAG(mode, paramsTrue_O, paramsTrue_O, t0,
    #                           xx0vec, sig0, sigd, ts)
    
    # append outputs to lists
    xxvecList_O.append(xxvec)
    tvecList_O.append(tvecEval)
    sigdvecList_O.append(sigdList)
    tsvecList_O.append(tsList)
    raErrList_O.append(raErr)
    DVList_O.append(DV)
    tsfList_O.append(tsList[-1])
    sigdfList_O.append(sigdList[-1])
    rafList_O.append(raf)
    rpfList_O.append(rpf)
    engfList_O.append(engf)
    
    # save data to file for analysis
    print('saving file...')
    xxvecArr_O = np.empty(i+1, object)
    xxvecArr_O[:] = xxvecList_O
    tvecArr_O = np.empty(i+1, object)
    tvecArr_O[:] = tvecArr_O
    sigdvecArr_O = np.empty(i+1, object)
    sigdvecArr_O[:] = sigdvecList_O
    tsvecArr_O = np.empty(i+1, object) 
    tsvecArr_O[:] = tsvecList_O
    np.savez(outname,
             xxvecArr_O = xxvecArr_O,
             tvecArr_O = tvecArr_O,
             sigdvecArr_O = sigdvecArr_O,
             tsvecArr_O = tsvecArr_O,
             BCList_O = BCList_O,
             L_DList_O = L_DList_O,
             gam0List = gam0List,
             v0List = v0List,
             raErrList_O = raErrList_O,
             DVList_O = DVList_O,
             tsfList_O = tsfList_O,
             sigdfList_O = sigdfList_O,
             rafList_O = rafList_O,
             rpfList_O = rpfList_O,
             engfList_O = engfList_O)
    
    toc = time.time()
    print('\nRun {0:d} complete, {1:.2f} s elapsed'.format(i+1, toc-tic))






toc = time.time()

print('Total time elapsed: {0:.2f} s'.format(toc-tic))




















