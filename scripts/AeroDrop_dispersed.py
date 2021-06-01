# -*- coding: utf-8 -*-
"""
Created on Mon Feb  8 15:00:25 2021

AeroDrop_dispersed.py:
    Simulates orbiter and probe using FNP(A/E)G guidance under dispersions.
    
    NOTE: the way this is set up now, NPC has perfect knowledge of dispersed
    state, including at entry, but has only nominal values for BC, L/D, and atm

@author: Samuel Albert
"""

from conversions import getApsesSphPR
from sim import Params
from guidance import updateFNPAG, updateFNPEG, engEvent
import ODEs
import planetaryconstants as constants
from atm import getRho_from_table

import numpy as np
from scipy.integrate import solve_ivp
from scipy.stats import norm, uniform
import scipy.interpolate as interp
from random import randint
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
# SETTINGS
# =============================================================================
# number of Monte Carlo runs
Nmc = 100

# updates on during nominal profiles?
updatesNominal = False

# verbose during nominal profiles?
verboseNominal = False

# udpates on during Monte Carlo trials?
updatesMC = True

# verbose during Monte Carlo trials?
verboseMC = False

# dummy data file?
dummy = False

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
    print('\nRUNNING FNPAG...')
    
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
    if verbose:
        print()
    
    for ind, t in enumerate(tvec):
        # update guidance
        if updatesOn:
            sigd = updateFNPAG(xxvec[:,-1], t, ts, sig0, sigd, phase, mode,
                               paramsNom)
            
        
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
    
    print('final apoapsis error: {0:.3e} m'.format(raErr))
    print('delta-V for periapsis raise: {0:.3f} m/s'.format(DV1))
    print('delta-V for apoapsis correction: {0:.3f} m/s'.format(DV2))
    print('total delta-V: {0:.3f} m/s'.format(DV))
    
    return xxvec, tvecEval, raf, rpf, engf, raErr, DV,\
        tsList, sigdList, tvecP1, tvecP2, xswitchvec


# =============================================================================
# Main FNPEG Function
# =============================================================================

def doFNPEG(paramsTrue, paramsNom, t0, xx0vec, sig0, e0, verbose = True,
            plotsOn = True, updatesOn = True):
    '''
    main function for FNPEG. paramsTrue holds real values used in propagation,
        paramsNom holds nominal values used in prediction for guidance updates.
    Uses augmented spherical planet-relative EOMs.
    '''
    print('\nRUNNING FNPEG...')
    
    # =========================================================================
    # Set up events
    # =========================================================================
    event1 = lambda t, y: ODEs.above_max_alt_sph(t, y, paramsTrue)
    event1.terminal = True
    event1.direction = 1
    event2 = lambda t, y: ODEs.below_min_alt_sph(t, y, paramsTrue)
    event2.terminal = True
    
    event3 = lambda t, y: engEvent(t, y, paramsTrue)
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
    
    evec = np.empty(1)
    evec[:] = np.NaN
    sigvec = np.empty(1)
    sigvec[:] = np.NaN
    
    sig0List = []
    
    t = t0
    ind = 0
    xx0veci = xx0vec * 1
    
    for ind, t in enumerate(tvec):
        
        # update guidance
        if updatesOn:
            sig0 = updateFNPEG(xx0veci, t, sig0, e0, paramsNom)
            
        if verbose:
            print('FNPEG: updating guidance at time {0:.3f}  '\
                  '  sig0 = {1:.3f} deg'.format(t, sig0))
        sig0List.append(sig0)
        
        # propagate until next guidance update or switching time
        tspan = (t, t + paramsTrue.dtGdn)
        soli = solve_ivp(lambda t, y: ODEs.sphericalEntryEOMsAug(t, y,
                                                                 sig0, e0,
                                                                 paramsTrue),
                         tspan, xx0veci,
                         rtol = paramsTrue.rtol, atol = paramsTrue.atol,
                         events = (event1, event2, event3))
        
        # get bank angle history
        eveci = paramsTrue.p.mu * 1e9 / soli.y[0,:] - soli.y[3,:]**2 / 2
        sigveci = sig0 + (eveci - e0) /\
            (paramsTrue.ef - e0) * (paramsTrue.sigf - sig0)
        
        # append results
        xxvec = np.append(xxvec, soli.y, axis = 1)
        tvecEval = np.append(tvecEval, soli.t)
        evec = np.append(evec, eveci)
        sigvec = np.append(sigvec, sigveci)
        xx0veci = xxvec[:,-1]
        
        if soli.status == 1:
            break
        
    # trim off first array elements
    xxvec = xxvec[:,1:]
    tvecEval = tvecEval[1:]
    evec = evec[1:]
    sigvec = sigvec[1:]
    
    if plotsOn:
        # uses global variable ax2 for plotting!
        ax2.plot(xxvec[3,:]/1e3, xxvec[0,:]/1e3 - paramsNom.p.rad)
        
    # =========================================================================
    # Compute target errors
    # =========================================================================
    xxfvec = xxvec[:,-1]
    sfErr = np.degrees(xxfvec[6])
    hfErr = xxfvec[0] - paramsTrue.rf
    vfErr = xxfvec[3] - paramsTrue.vf
    
    print('Final range error: {0:.6e} deg'.format(sfErr))
    print('Final altitude error: {0:.6f} m'.format(hfErr))
    print('Final velocity error: {0:.6f} m/s'.format(vfErr))
    
    return xxvec, tvecEval, sfErr, hfErr, vfErr, evec, sigvec, sig0List


# =============================================================================
# Main Ballistic Probe Function
# =============================================================================
def doBallisticProbe(params, t0, xx0vec, plotsOn = True):
    '''
    main function for passive ballistic probe (no control or guidance).
    Uses augmented spherical planet-relative EOMs.
    '''
    print('\nRUNNING PASSIVE PROBE...')
    # set dummy sig0 and e0 values
    sig0 = 90
    e0 = 0
    
    # =========================================================================
    # Set up events
    # =========================================================================
    event1 = lambda t, y: ODEs.above_max_alt_sph(t, y, params)
    event1.terminal = True
    event1.direction = 1
    event2 = lambda t, y: ODEs.below_min_alt_sph(t, y, params)
    event2.terminal = True
    
    # =========================================================================
    # Main integration loop
    # =========================================================================
    tspan = (t0, params.tf)
    sol = solve_ivp(lambda t, y: ODEs.sphericalEntryEOMsAug(t, y, sig0, e0,
                                                            params),
                    tspan, xx0vec, rtol = params.rtol, atol = params.atol,
                    events = (event1, event2))
    
    xxvec = sol.y
    tvec = sol.t
    
    if plotsOn:
        # uses global variable ax3 for plotting!
        ax3.plot(xxvec[3,:]/1e3, xxvec[0,:]/1e3 - params.p.rad)
    
    # =========================================================================
    # Compute target errors
    # =========================================================================
    xxfvec = xxvec[:,-1]
    sfErr = np.degrees(xxfvec[6])
    hfErr = xxfvec[0] - params.rf
    vfErr = xxfvec[3] - params.vf
    
    print('Final range error: {0:.6e} deg'.format(sfErr))
    print('Final altitude error: {0:.6f} m'.format(hfErr))
    print('Final velocity error: {0:.6f} m/s'.format(vfErr))
    
    return xxvec, tvec, sfErr, hfErr, vfErr
    
        

# =============================================================================
# Generic Setup
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
# BC is set individually below

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
    
# other settings and constants
paramsNom_O.rtol = 1e-10
paramsNom_O.atol = 1e-10
paramsNom_O.errtol1 = 0
paramsNom_O.errtol2 = 0
paramsNom_O.dtGdn = 1 # s
paramsNom_O.hmin = 3
paramsNom_O.hmax = paramsNom_O.p.halt + 1e-7
paramsNom_O.tf = 6000

# copy generic params for probe and ballistic probe
paramsNom_P = copy.deepcopy(paramsNom_O)

# =============================================================================
# Orbiter Setup
# =============================================================================
# vehicle
paramsNom_O.BC = 130
paramsNom_O.A = paramsNom_O.m / (paramsNom_O.BC * paramsNom_O.CD)
paramsNom_O.CL = paramsNom_O.CD * paramsNom_O.LD
paramsNom_O.Rn = np.sqrt(paramsNom_O.A / np.pi) / 2

# search brackets for Brent's Method
paramsNom_O.sig1 = 0
paramsNom_O.sig2 = 180
paramsNom_O.ts1 = 1
paramsNom_O.ts2 = 300

# target state
paramsNom_O.raStar = (250 + paramsNom_O.p.rad) * 1e3
paramsNom_O.rpStar = (250 + paramsNom_O.p.rad) * 1e3

# nominal simulation values
t0 = 0
sig0_O = 0
sigd = 150
ts = 152.64287589323638

# =============================================================================
# Lifting Probe Setup
# =============================================================================
# vehicle
paramsNom_P.BC = 35
paramsNom_P.A = paramsNom_P.m / (paramsNom_P.BC * paramsNom_P.CD)
paramsNom_P.CL = paramsNom_P.CD * paramsNom_P.LD
paramsNom_P.Rn = np.sqrt(paramsNom_P.A / np.pi) / 2

# allow full range of bank angle magnitude:
paramsNom_P.sig01 = 0
paramsNom_P.sig02 = 180

# target states and constraints
paramsNom_P.sigf = 60 # fixed final constraint, chosen arbitrarily
paramsNom_P.sf = 0 # always target 0 range-to-go
paramsNom_P.rf = (15 + paramsNom_P.p.rad) * 1e3 # 15 km altitude
paramsNom_P.vf = 300 # from nominal trajectory
# comput target e from target r and v values
paramsNom_P.ef = paramsNom_P.p.mu * 1e9 / paramsNom_P.rf - paramsNom_P.vf**2 / 2
# paramsNom_P.ef = 2e7

# augment initial state to include range
s0 = 0.20628354331374213 # from nominal trajectory
# s0 = 0
xx0vecAug = np.append(xx0vec, s0)

# other sim params
# sig0 = np.radians(20) # initial guess
# sig0_P = 118.52940036607743
sig0_P = 139.31511427799256
r0 = (paramsNom_P.alt + paramsNom_P.p.rad) * 1e3
v0 = paramsNom_P.vmagWR * 1e3
e0 = paramsNom_P.p.mu * 1e9 / r0 - v0**2 / 2

# =============================================================================
# Ballistic Probe Setup
# =============================================================================
paramsNom_PBC = copy.deepcopy(paramsNom_P)
paramsNom_PBC.LD = 0
paramsNom_PBC.CL = paramsNom_PBC.CD * paramsNom_PBC.LD

# change minimum altitude to equal target altitude of lifting probe
paramsNom_PBC.hmin = paramsNom_P.rf/1e3 - paramsNom_P.p.rad

# separate target params for velocity and range
paramsNom_PBC.vf = 353.09785084111587
s0_PBC = 0.21635403165302095
xx0vecAug_PBC = np.append(xx0vec, s0_PBC)

# =============================================================================
# Demonstrate dynamics for nominal scenarios
# =============================================================================
mode = 1
xxvec1, tvecEval, raf, rpf, engf, raErr, DV, tsList, sigdList, tvecP1, tvecP2,\
    xswitchvec = doFNPAG(mode, paramsNom_O, paramsNom_O, t0, xx0vec, sig0_O,
                         sigd, ts, plotsOn = False, updatesOn = updatesNominal,
                         verbose = verboseNominal)

fig = plt.figure()
ax = fig.add_subplot(111)
ax.plot(xxvec1[3,:]/1e3, xxvec1[0,:]/1e3 - paramsNom_O.p.rad, 'k', linewidth = 4,
        label = 'nominal trajectory')
ax.plot(xswitchvec[3]/1e3, xswitchvec[0]/1e3 - paramsNom_O.p.rad, 'ko',
        markersize = 10)
ax.grid()
ax.legend()
ax.set_xlabel('planet-relative velocity, km/s')
ax.set_ylabel('spherical altitude, km')
ax.set_title('Orbiter trajectories')

xxvec2, tvecEval, sfErr, hfErr, vfErr, evec, sigvec,\
    sig0List = doFNPEG(paramsNom_P, paramsNom_P, t0, xx0vecAug, sig0_P, e0,
                       plotsOn = False, updatesOn = updatesNominal,
                       verbose = verboseNominal)

fig2 = plt.figure()
ax2 = fig2.add_subplot(111)
ax2.plot(xxvec2[3,:]/1e3, xxvec2[0,:]/1e3 - paramsNom_P.p.rad, 'k', linewidth = 4,
        label = 'nominal trajectory')
ax2.grid()
ax2.legend()
ax2.set_xlabel('planet-relative velocity, km/s')
ax2.set_ylabel('spherical altitude, km')
ax2.set_title('Lifting probe trajectories')

xxvec3, tvecEval, sfErr, hfErr, vfErr = doBallisticProbe(paramsNom_PBC, t0,
                                                        xx0vecAug_PBC,
                                                        plotsOn = False)
fig3 = plt.figure()
ax3 = fig3.add_subplot(111)
ax3.plot(xxvec3[3,:]/1e3, xxvec3[0,:]/1e3 - paramsNom_P.p.rad, 'k', linewidth = 4,
        label = 'nominal trajectory')
ax3.grid()
ax3.legend()
ax3.set_xlabel('planet-relative velocity, km/s')
ax3.set_ylabel('spherical altitude, km')
ax3.set_title('Ballistic probe trajectories')

# =============================================================================
# Compute target FNPEG params from nominal FNPEG run (comment out during MC)
# =============================================================================
print('\nCOMPUTING NOMINAL PARAMETERS...')
r = xxvec2[0,:]
vmag = xxvec2[3,:] / 1e3

h = (r - paramsNom_PBC.p.rad * 1e3) / 1e3

# get density and speed of sound at each altitude step
rho = []
assfun = interp.interp1d(paramsNom_PBC.p.sound[0,:], paramsNom_PBC.p.sound[1,:])
ass = assfun(h)
for hi in h:
    rho.append(getRho_from_table(paramsNom_PBC.atmdat, hi))
rho = np.asarray(rho)
    
# compute mach number and dynamic pressure at each step
Mvec = vmag * 1e3 / ass
Qinc = 1/2 * rho * (vmag*1e3)**2

# range
lon0 = xxvec2[1,0]
lonf = xxvec2[1,-1]
lat0 = xxvec2[2,0]
latf = xxvec2[2,-1]
dlon = abs(lonf - lon0)

dsig = np.arccos(np.sin(lat0) * np.sin(latf)\
                 + np.cos(lat0) * np.cos(latf) * np.cos(dlon))
    
# print
print('Lifting Probe final values:')
print('Final altitude: {:.3f} km'.format(h[-1]))
print('Final velocity {:.3f} m/s'.format(vmag[-1]*1e3))
print('Final Mach number: {:.3f}'.format(Mvec[-1]))
print('Final dynamic pressure: {:.3f} Pa'.format(Qinc[-1]))
print('Range traversed: {:.3f} rad\n'.format(dsig))

# =============================================================================
# Compute target passive probe params from nominal run (comment out during MC)
# =============================================================================
r = xxvec3[0,:]
vmag = xxvec3[3,:] / 1e3

h = (r - paramsNom_PBC.p.rad * 1e3) / 1e3

# get density and speed of sound at each altitude step
rho = []
assfun = interp.interp1d(paramsNom_PBC.p.sound[0,:], paramsNom_PBC.p.sound[1,:])
ass = assfun(h)
for hi in h:
    rho.append(getRho_from_table(paramsNom_PBC.atmdat, hi))
rho = np.asarray(rho)

# compute mach number and dynamic pressure at each step
Mvec = vmag * 1e3 / ass
Qinc = 1/2 * rho * (vmag*1e3)**2

# range
lon0 = xxvec3[1,0]
lonf = xxvec3[1,-1]
lat0 = xxvec3[2,0]
latf = xxvec3[2,-1]
dlon = abs(lonf - lon0)

dsig2 = np.arccos(np.sin(lat0) * np.sin(latf)\
                 + np.cos(lat0) * np.cos(latf) * np.cos(dlon))
    
# print
print('\nPassive Probe final values:')
print('Final altitude: {:.3f} km'.format(h[-1]))
print('Final velocity {:.3f} m/s'.format(vmag[-1]*1e3))
print('Final Mach number: {:.3f}'.format(Mvec[-1]))
print('Final dynamic pressure: {:.3f} Pa'.format(Qinc[-1]))
print('Range traversed: {:.3f} rad\n'.format(dsig2))

# =============================================================================
# Orbiter Dispersions Setup
# =============================================================================
# copy, DO NOT ASSIGN, paramsTrue for dispersions. leave paramsNom alone.
paramsTrue_O = copy.deepcopy(paramsNom_O)

# load full 5000-case Monte Carlo atmsophere dataset from MarsGRAM-2010
filename = '../data/Mars_0.1_50000.npz'
data = np.load(filename)
densAll = data['densTot']
h = data['h']
NATM = densAll.shape[1]
atmindices = [randint(0, NATM) for p in range(Nmc)] # all Nmc atm indices
del data
# filename = '../data/Mars_0.1_5000.txt'
# densAll, densMean, h = getMarsGRAMDensTableAll(filename, Nmc)

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

atmindList = [] # atm indices for only trials that actually get run

# =============================================================================
# Lifting Probe Dispersions Setup
# =============================================================================
paramsTrue_P = copy.deepcopy(paramsNom_P)

# define vehicle dispersions
BCMean_P = paramsNom_P.BC
BCLB_P = 0.95 * paramsNom_P.BC
BCRng_P = 0.1 * paramsNom_P.BC

L_DMean_P = paramsNom_P.LD
L_DLB_P = 0.95 * paramsNom_P.LD
L_DRng_P = 0.1 * paramsNom_P.LD

# initialize lists for actual dispersed parameters and results
xxvecList_P = []
tvecList_P = []
sig0ListList = []
sigvecList = []

sfErrList_P = []
hfErrList_P = []
vfErrList_P = []

BCList_P = []
L_DList_P = []

# =============================================================================
# Passive Probe Dispersions Setup
# =============================================================================
paramsTrue_PBC = copy.deepcopy(paramsNom_PBC)

# vehicle dispersions same as lifting probe

# initialize lists
xxvecList_PBC = []
tvecList_PBC = []

sfErrList_PBC = []
hfErrList_PBC = []
vfErrList_PBC = []


# output filename
if dummy:
    outname = '../results/DUMMY_AeroDrop_dispersed_'\
        + str(Nmc) + '_' + datestring
else:
    outname = '../results/AeroDrop_dispersed_' + str(Nmc) + '_' + datestring

# =============================================================================
#  Monte Carlo loop
# =============================================================================
# for i in range(Nmc):
runi = 1
for i in atmindices:
    # initialize random variables
    BC_O = uniform.rvs(size = 1, loc = BCLB_O, scale = BCRng_O)[0]
    L_D_O = uniform.rvs(size = 1, loc = L_DLB_O, scale = L_DRng_O)[0]
    BC_P = uniform.rvs(size = 1, loc = BCLB_P, scale = BCRng_P)[0]
    L_D_P = uniform.rvs(size = 1, loc = L_DLB_P, scale = L_DRng_P)[0]
    gam0 = norm.rvs(size = 1, loc = gam0Mean, scale = gam0STD)[0]
    v0 = norm.rvs(size = 1, loc = v0Mean, scale = v0STD)[0]
    
    # orbiter true params
    paramsTrue_O.BC = BC_O
    paramsTrue_O.LD = L_D_O
    xx0vec[3] = v0
    xx0vec[4] = gam0
    
    # lifting probe true params
    paramsTrue_P.BC = BC_P
    paramsTrue_P.LD = L_D_P
    xx0vecAug = np.append(xx0vec, s0)
    
    # passive probe true params
    paramsTrue_PBC.BC = BC_P
    xx0vecAug_PBC = np.append(xx0vec, s0_PBC)
    
    # get GRAM atm profile
    paramsTrue_O.atmdat = np.array([h,densAll[:,i]])
    paramsTrue_O.atmdat = paramsTrue_O.atmdat[:,paramsTrue_O.atmdat[0,:].argsort()]
    paramsTrue_P.atmdat = paramsTrue_O.atmdat
    paramsTrue_PBC.atmdat = paramsTrue_O.atmdat
    
    # append inputs to lists (skip GRAM profiles, can always get those later)
    BCList_O.append(BC_O)
    L_DList_O.append(L_D_O)
    BCList_P.append(BC_P)
    L_DList_P.append(L_D_P)
    gam0List.append(gam0)
    v0List.append(v0)
    atmindList.append(i)
    
    # run FNPAG simulation
    xxvec, tvecEval, raf, rpf, engf, raErr, DV, tsList, sigdList, tvecP1,\
        tvecP2, xswwitchvec = doFNPAG(mode, paramsTrue_O, paramsNom_O, t0,
                              xx0vec, sig0_O, sigd, ts, updatesOn = updatesMC,
                              verbose = verboseMC)
    
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
    
    # run FNPEG simulation
    xxvec, tvecEval, sfErr, hfErr, vfErr, evec, sigvec,\
        sig0List = doFNPEG(paramsTrue_P, paramsNom_P, t0, xx0vecAug,
                           sig0_P, e0, updatesOn = updatesMC,
                           verbose = verboseMC)
    
    # append outputs to lists
    xxvecList_P.append(xxvec)
    tvecList_P.append(tvecEval)
    sig0ListList.append(sig0List)
    sigvecList.append(sigvec)
    sfErrList_P.append(sfErr)
    hfErrList_P.append(hfErr)
    vfErrList_P.append(vfErr)
    
    # run passive probe simulation
    xxvec, tvecEval, sfErr, hfErr,\
        vfErr = doBallisticProbe(paramsTrue_PBC, t0, xx0vecAug_PBC)
    xxvecList_PBC.append(xxvec)
    tvecList_PBC.append(tvecEval)
    sfErrList_PBC.append(sfErr)
    hfErrList_PBC.append(hfErr)
    vfErrList_PBC.append(vfErr)
    
    
    
    # save data to file for analysis
    print('saving file...')
    xxvecArr_O = np.empty(runi, object)
    xxvecArr_O[:] = xxvecList_O
    tvecArr_O = np.empty(runi, object)
    tvecArr_O[:] = tvecList_O
    sigdvecArr_O = np.empty(runi, object)
    sigdvecArr_O[:] = sigdvecList_O
    tsvecArr_O = np.empty(runi, object) 
    tsvecArr_O[:] = tsvecList_O
    xxvecArr_P = np.empty(runi, object)
    xxvecArr_P[:] = xxvecList_P
    tvecArr_P = np.empty(runi, object)
    tvecArr_P[:] = tvecList_P
    sigvecArr = np.empty(runi, object)
    sigvecArr[:] = sigvecList
    sig0ListArr = np.empty(runi, object)
    sig0ListArr[:] = sig0ListList
    xxvecArr_PBC = np.empty(runi, object)
    xxvecArr_PBC[:] = xxvecList_PBC
    tvecArr_PBC = np.empty(runi, object)
    tvecArr_PBC[:] = tvecList_PBC
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
                engfList_O = engfList_O,
                xxvecArr_P = xxvecArr_P,
                tvecArr_P = tvecArr_P,
                sigvecArr = sigvecArr,
                sig0ListArr = sig0ListArr,
                sfErrList_P = sfErrList_P,
                hfErrList_P = hfErrList_P,
                vfErrList_P = vfErrList_P,
                xxvecArr_PBC = xxvecArr_PBC,
                tvecArr_PBC = tvecArr_PBC,
                sfErrList_PBC = sfErrList_PBC,
                hfErrList_PBC = hfErrList_PBC,
                vfErrList_PBC = vfErrList_PBC,
                atmindList = atmindList)
    
    toc = time.time()
    print('\nRun {0:d} complete, {1:.2f} s elapsed, atmID: {2:d}'\
          .format(runi, toc-tic, i))
    runi += 1






toc = time.time()

print('Total time elapsed: {0:.2f} s'.format(toc-tic))




















