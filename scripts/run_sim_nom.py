# -*- coding: utf-8 -*-
"""
Created on Thu May 28 15:29:37 2020

@author: Samuel Albert
"""

import numpy as np
import matplotlib.pyplot as plt
import time

import ODEs
from sim import simRun
from atm import getMCdens
import constants

class params:
    def __init__(self, dMode):
        self.dMode = dMode
    
# plt.close('all')

tic = time.time()


dMode = 'table'
params = params(dMode)
params.p = constants.EARTH

### GET ATM TABLE FROM OLD EARTH GRAM .CSV FILE
filename = '../data/atm_earth_gram2016.csv'
atmdata_raw = np.genfromtxt(filename, delimiter=',', names=True,
                            encoding='utf-8-sig')
# at some point would be good to build this as a pandas df instead of np array
rhoTable = np.array([atmdata_raw['alt']/1e3,atmdata_raw['density']])
params.atmdat = rhoTable


# ### GET ATM TABLE FROM BINARY EARTHGRAM DATA FILE
# filename = '../data/rawOutput.txt'
# # get Nmc atmosphere profiles
# Nmc = 1
# i_trial = 0
# densPert, densMean, h = getMCdens(filename, Nmc)
# # at some point would be good to build this as a pandas df instead of np array
# rhoTable = np.array([h,densPert[:,i_trial]])
# params.atmdat = rhoTable

### VEHICLE PARAMS
params.m = 2000 # kg 
params.A = 30 # m^2
params.CL = 0.6
params.BC = 51.282051282051285
# get CD from BC
# params.CD = params.m / (params.BC * params.A)
params.CD = 1.3

### INITIAL STATE
# params.x0 = np.array([6478.100, 0, 0]) # inertial frame and coordinates
# params.v0 = np.array([-1.25145758, 0.47238996, 7.9015096])

params.x0 = np.array([6478100,           0,           0]) / 1e3
# v0vec_N = np.array([-671.533934883426,            472.3899576546,          10979.4827826405]) / 1e3
params.v0 = np.array([-0.67153393,  0.47238996, 10.97948278])

### CONTROL STATE
params.bank = 60 # full-lift-up, deg

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

### CALL SIMRUN
sol = simRun(params, tspan, events)
rvec_N = sol.y[0:3,:]
vvec_N = sol.y[3:6,:] 

### PLOTS
alt = np.linalg.norm(rvec_N, axis=0) - params.p.rad #/1e3
vmag = np.linalg.norm(vvec_N, axis=0)

fig = plt.figure(2)
ax = fig.add_subplot(111)
ax.plot(vmag, alt)
ax.set_xlabel('Inertial velocity, km/s')
ax.set_ylabel('Altitude, km')
ax.grid()

toc = time.time()
print('Time elapsed: %.2f s' % (toc-tic))