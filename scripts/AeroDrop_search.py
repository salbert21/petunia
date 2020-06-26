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

import constants
import ODEs
from sim import simRun
from atm import getMCdens
from conversions import LLAEHV2RV, RV2LLAEHV

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
    
    # convert final state params
    lat, lon, alt, fpa, hda, vmag = RV2LLAEHV(rfvec_N, vfvec_N, params, tf)
    
    ## TODO - compute peak heat rate, add to output
    
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
    outs.rvec_N = rfvec_N
    outs.vvec_N = vfvec_N
    outs.t = tf
    
    
    
    
    
    
    return outs
    
    
                    
    

# =============================================================================
# MAIN - GRID SEARCH
# =============================================================================
tic = time.time()

### CREATE params INPUT CLASS
params = params()
params.p = constants.EARTH

### INPUT ATM TABLE - GET ATM TABLE FROM BINARY EARTHGRAM DATA FILE
params.dMode = 'table'
filename = '../data/rawOutput.txt'
# get Nmc atmosphere profiles
Nmc = 1
i_trial = 0
densPert, densMean, h = getMCdens(filename, Nmc)
# at some point would be good to build this as a pandas df instead of np array
rhoTable = np.array([h,densPert[:,i_trial]])
params.atmdat = rhoTable

### VEHICLE PARAMS (NOT CHANGED DURING GRID SEARCH)
params.m = 2000 # kg 
params.A = 15 # m^2
params.CL = 0.1212

### INITIAL STATE (COMPONENTS NOT CHANGED DURING GRID SEARCH)
params.lat = 0
params.lon = 0
params.alt = 100.0
params.hda = 0
params.vmag = 11

### CONTROL STATE
params.bank = 0 # deg

### TIME VECTOR AND EXIT CONDITIONS
# should always stop on an exit condition
tspan = np.linspace(0,5000,10000) # don't make too long or results get choppy!

# exit conditions:
params.hmin = 10
params.hmax = 125

event1 = lambda t, y: ODEs.below_min_alt(t, y, params)
event1.terminal = True

event2 = lambda t, y: ODEs.above_max_alt(t, y, params)
event2.terminal = True

events = (event1, event2)

### GRID SEARCH
# currently a single value each for testing the main function
params.efpa = -7
params.BC = 130

outs = outs() # blank instance of output class
outs = main(params, tspan, events, outs)










