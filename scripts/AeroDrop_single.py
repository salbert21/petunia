# -*- coding: utf-8 -*-
"""
Created on Wed Jun 17 15:25:31 2020

@author: Samuel Albert
"""

# =============================================================================
# IMPORTS AND CLASS DEFINITIONS
# =============================================================================

import numpy as np
import time
import matplotlib.pyplot as plt
import sys
from scipy.integrate import trapz

import constants
import ODEs
from sim import simRun
from atm import getMCdens, getRho_from_table
from conversions import LLAEHV2RV, RV2LLAEHV, VN2Vinf

# #### TEMP DEBUGGING CODE ####
# import warnings
# warnings.simplefilter('error')
# ####

class params:
       
    class p:
        pass 
    
class Outs:
    
    class s0:
        pass
    class sf:
        pass
    class v:
        pass

# =============================================================================
# MAIN FUNCTION DEFINITION
# =============================================================================

def main(params, tspan, events, outs):
    # get CD from BC
    params.CD = params.m / (params.BC * params.A)
    
    ### CONVERT TO INERTIAL INITIAL STATE
    params.x0, params.v0 = LLAEHV2RV(params.lat, params.lon, params.alt,
                           params.efpa, params.hda, params.vmag, params,
                           tspan[0])
    
    
    ### CALL SIMRUN
    sol = simRun(params, tspan, events, verbose=True)
    rvec_N = sol.y[0:3,:]
    vvec_N = sol.y[3:6,:] 
    
    ### GET OUTPUT PARAMS OF INTEREST
    # final inertial state
    rfvec_N = rvec_N[:,-1]
    vfvec_N = vvec_N[:,-1]
    tf = sol.t[-1]
    
    # convert final state params (these are inertial)
    lat, lon, alt, fpa, hda, vmag = RV2LLAEHV(rfvec_N, vfvec_N, params, tf)
    
    ## Compute peak heat rate, add to output
    # get airspeed at each time, vInf
    vInfvec_N = []
    for i in range(vvec_N.shape[1]):
        vInfvec_N.append(VN2Vinf(rvec_N[:,i], vvec_N[:,i], params, sol.t[i]))
    
    vInfvec_N = np.asarray(vInfvec_N).T
    vInf = np.linalg.norm(vInfvec_N, axis=0)
    
    # get density at each time
    r = np.linalg.norm(rvec_N, axis=0)
    h = r - params.p.rad
    if params.dMode == 'fun':
        rho = params.dFun(h)
    elif params.dMode == 'table':
        rho = getRho_from_table(params.atmdat, h)
    else:
        sys.exit('atm mode not recognized')
        
    # calculate S-G heat rate without coefficient, SG
    SG = []
    for i in range(len(rho)):
        SG.append(np.sqrt(rho[i] / params.Rn) * vInf[i]**3)
    
    SG = np.asarray(SG)
    SGpeak = SG.max()
    
    # calculate peak heat rate
    qpeak = params.p.k * SGpeak * 1e5 # puts q in W/cm^2 units
    q = params.p.k * SG * 1e5
    
    # now integrate numerically to get heat load
    Qload = trapz(q, sol.t) # J / cm^2
    
    ## Compute max g, add to output
    # run through dynamics again to get forces output
    Fgvec_N = np.empty((3, len(sol.t)))
    Fgvec_N[:] = np.NaN
    FLvec_N = np.empty((3, len(sol.t)))
    FLvec_N[:] = np.NaN
    FDvec_N = np.empty((3, len(sol.t)))
    FDvec_N[:] = np.NaN
    
    for i, (ti, xi, vi) in enumerate(zip(sol.t, rvec_N.T, vvec_N.T)):
        yyi = np.block([xi, vi])
        Fgvec_N[:,i], FLvec_N[:,i], FDvec_N[:,i] = ODEs.dynamics(ti, yyi,
                                                  params, returnForces=True)
    # compute g-load at each time
    gload = (np.linalg.norm(FLvec_N + FDvec_N - Fgvec_N, axis=0) / params.m)\
        / (constants.G0 / 1e3) # divide g0 by 1000 to get it in km/s^2 units
        
    gpeak = gload.max()
    
    ## Compute apoapsis at final time, add to output
    vf = np.linalg.norm(vfvec_N)
    rf = np.linalg.norm(rfvec_N)
    engf = vf**2 / 2 - params.p.mu / rf
    hfvec_N = np.cross(rfvec_N, vfvec_N)
    hf = np.linalg.norm(hfvec_N)
    
    af = - params.p.mu / (2 * engf)
    eccf = np.sqrt(1 + 2 * engf * hf**2 / params.p.mu**2)
    
    raf = af * (1 + eccf)
    haf = raf - params.p.rad
    
    ### ASSIGN OUTPUTS TO outs CLASS   
    # initial state
    outs.s0.lat = params.lat
    outs.s0.lon = params.lon
    outs.s0.alt = params.alt
    outs.s0.efpa = params.efpa
    outs.s0.hda = params.hda
    outs.s0.vmag = params.vmag
    outs.s0.rvec_N = params.x0
    outs.s0.vvec_N = params.v0
    
    outs.tspan = tspan
    
    # vehicle
    outs.v.m = params.m
    outs.v.A = params.A
    outs.v.CL = params.CL
    outs.v.CD = params.CD
    outs.v.BC = params.BC
    outs.v.Rn = params.Rn
    outs.v.bank = params.bank
    
    
    # final state
    outs.sf.lat = lat
    outs.sf.lon = lon
    outs.sf.alt = alt
    outs.sf.fpa = fpa
    outs.sf.hda = hda
    outs.sf.vmag = vmag
    outs.sf.rvec_N = rfvec_N
    outs.sf.vvec_N = vfvec_N
    outs.sf.raf = raf
    outs.sf.haf = haf
    outs.sf.t = tf
    
    # peak values and total loads
    outs.SGpeak = SGpeak
    outs.qpeak = qpeak
    outs.Qload = Qload
    outs.gpeak = gpeak

    
    
    # may not want these always on since they have values at each time step
    outs.rvec_N = rvec_N
    outs.vvec_N = vvec_N
    outs.tvec = sol.t
    outs.q = q
    outs.gload = gload

    return outs
    
    
                    
    

# =============================================================================
# MAIN - GRID SEARCH
# =============================================================================
tic = time.time()
plt.close('all')

### CREATE params INPUT CLASS
params = params()
params.p = constants.EARTH

### INPUT ATM TABLE - GET ATM TABLE FROM EARTHGRAM DATA FILE
params.dMode = 'table'
filename = '../data/dat_raw_Earth_nom.txt'
atmdata = np.genfromtxt(filename, names=True)
atmdata.sort(order='Hgtkm') # put in ascending altitude order
params.atmdat = np.array([atmdata['Hgtkm'], atmdata['DensMean']])

# ## get Nmc atmosphere profiles
# Nmc = 1
# i_trial = 0
# densPert, densMean, h = getMCdens(filename, Nmc)
# # at some point would be good to build this as a pandas df instead of np array
# # rhoTable = np.array([h,densPert[:,i_trial]])
# # load nominal density:
# rhoTable = np.array([h,densPert[:,i_trial]])
# params.atmdat = rhoTable

### VEHICLE PARAMS (NOT CHANGED DURING GRID SEARCH)
params.m = 2000 # kg 
params.Rn = 0.25 # m, effective nose radius # TODO: CHANGE THIS VALUE
# params.A = 15 # m^2
params.A = 30
# params.CL = 0.1212
params.CL = 0.6

### INITIAL STATE (COMPONENTS NOT CHANGED DURING GRID SEARCH)
params.lat = 40
params.lon = 100
params.alt = 100 - 1e-7
params.hda = 0
params.vmag = 11

### CONTROL STATE
params.bank = 60 # deg

### TIME VECTOR AND EXIT CONDITIONS
# should always stop on an exit condition
tspan = (0,1000) # don't make too long or results get choppy!

# exit conditions:
params.hmin = 30
params.hmax = 100

event1 = lambda t, y: ODEs.below_min_alt(t, y, params)
event1.terminal = True

event2 = lambda t, y: ODEs.above_max_alt(t, y, params)
event2.terminal = True

events = (event1, event2)


### SINGLE RUN
params.efpa = -5
params.BC = 80
outs = Outs()
outs = main(params, tspan, events, outs)



### PLOTS
alt = np.linalg.norm(outs.rvec_N, axis=0) - params.p.rad #/1e3
vmag = np.linalg.norm(outs.vvec_N, axis=0)

fig = plt.figure(2)
ax = fig.add_subplot(111)
ax.plot(vmag, alt)
ax.set_xlabel('Inertial velocity, km/s')
ax.set_ylabel('Altitude, km')
ax.grid()

fig = plt.figure()
ax = fig.add_subplot(111)
ax.plot(outs.tvec, outs.q)
ax.set_xlabel('time')
ax.set_ylabel('heat rate')
ax.grid()

fig = plt.figure()
ax = fig.add_subplot(111)
ax.plot(outs.tvec, vmag)
ax.set_xlabel('time')
ax.set_ylabel('vmag')
ax.grid()

toc = time.time()
print('Time elapsed: %.2f s' % (toc-tic))








