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

import constants
import ODEs
from sim import simRun
from atm import getMCdens, getRho_from_table
from conversions import LLAEHV2RV, RV2LLAEHV, VN2Vinf

class params:
       
    class p:
        pass
    
class outs:
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
    sol = simRun(params, tspan, events)
    rvec_N = sol.y[0:3,:]
    vvec_N = sol.y[3:6,:] 
    
    ### GET OUTPUT PARAMS OF INTEREST
    # final inertial state
    rfvec_N = rvec_N[:,-1]
    vfvec_N = vvec_N[:,-1]
    tf = sol.t[-1]
    
    # convert final state params (these are inertial)
    lat, lon, alt, fpa, hda, vmag = RV2LLAEHV(rfvec_N, vfvec_N, params, tf)
    
    ## TODO - compute peak heat rate, add to output
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
    qpeak = params.p.k * SGpeak / 1e5 # puts q in W/cm^2 units (I think)
    
    ## TODO - compute max g, add to output
    
    ## TODO - compute apoapsis, add to output
    
    
    ### ASSIGN OUTPUTS TO outs CLASS   
    # final state
    # outs = outs()
    outs.lat = lat
    outs.lon = lon
    outs.alt = alt
    outs.fpa = fpa
    outs.hda = hda
    outs.vmag = vmag
    outs.rfvec_N = rfvec_N
    outs.vfvec_N = vfvec_N
    outs.t = tf
    outs.SGpeak = SGpeak
    outs.qpeak = qpeak
    
    # may not want these always on
    outs.rvec_N = rvec_N
    outs.vvec_N = vvec_N
    
    
    
    
    
    
    return outs
    
    
                    
    

# =============================================================================
# MAIN - GRID SEARCH
# =============================================================================
tic = time.time()
# plt.close('all')

### CREATE params INPUT CLASS
params = params()
params.p = constants.EARTH

# ### INPUT ATM TABLE - GET ATM TABLE FROM BINARY EARTHGRAM DATA FILE
# params.dMode = 'table'
# filename = '../data/rawOutput.txt'
# # get Nmc atmosphere profiles
# Nmc = 1
# i_trial = 0
# densPert, densMean, h = getMCdens(filename, Nmc)
# # at some point would be good to build this as a pandas df instead of np array
# rhoTable = np.array([h,densPert[:,i_trial]])
# params.atmdat = rhoTable

### GET ATM TABLE FROM OLD EARTH GRAM .CSV FILE
params.dMode = 'table'
filename = '../data/atm_earth_gram2016.csv'
atmdata_raw = np.genfromtxt(filename, delimiter=',', names=True,
                            encoding='utf-8-sig')
# at some point would be good to build this as a pandas df instead of np array
rhoTable = np.array([atmdata_raw['alt']/1e3,atmdata_raw['density']])
params.atmdat = rhoTable

### VEHICLE PARAMS (NOT CHANGED DURING GRID SEARCH)
params.m = 2000 # kg 
params.Rn = 0.25 # m, effective nose radius # TODO: CHANGE THIS VALUE
# params.A = 15 # m^2
params.A = 30
# params.CL = 0.1212
params.CL = 0.6

# ### INITIAL STATE (COMPONENTS NOT CHANGED DURING GRID SEARCH)
# params.lat = 0
# params.lon = 0
# params.alt = 100.0
# params.hda = 0
# params.vmag = 11

### CONTROL STATE
params.bank = 60 # deg

### TIME VECTOR AND EXIT CONDITIONS
# should always stop on an exit condition
tspan = np.linspace(0,1800,10000) # don't make too long or results get choppy!

# exit conditions:
params.hmin = 30
params.hmax = 100

event1 = lambda t, y: ODEs.below_min_alt(t, y, params)
event1.terminal = True

event2 = lambda t, y: ODEs.above_max_alt(t, y, params)
event2.terminal = True

events = (event1, event2)

# ### GRID SEARCH
# # currently a single value each for testing the main function
# params.efpa = -3.5
params.BC = 51.282051282051285

# =============================================================================
# TEMPORARY DEBUG - overwrite initial params with inertial position from asim
# =============================================================================

r0vec_N = np.array([6478100,           0,           0]) / 1e3
# v0vec_N = np.array([-671.533934883426,            472.3899576546,          10979.4827826405]) / 1e3
v0vec_N = np.array([-0.67153393,  0.47238996, 10.97948278])

params.lat, params.lon, params.alt, params.efpa, params.hda, params.vmag =\
    RV2LLAEHV(r0vec_N, v0vec_N, params, 0)

outs = outs() # blank instance of output class
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

toc = time.time()
print('Time elapsed: %.2f s' % (toc-tic))








