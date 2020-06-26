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

class params:
    def __init__(self, dMode):
        self.dMode = dMode
    
plt.close('all')

tic = time.time()


### CONSTANTS
omE = 2 * np.pi / (0.99726968 * 86400) # rad/s rotation rate of Earth
radE = 6378.1363
mu = 3.986e5 # km^3/s^2

### ATM MODE
# dMode = 'fun'
dMode = 'table'

params = params(dMode)
params.p.mu = mu
params.p.om = omE
params.p.rad = radE

# ### GET ATM TABLE FROM OLD EARTH GRAM .CSV FILE
# filename = '../data/atm_earth_gram2016.csv'
# atmdata_raw = np.genfromtxt(filename, delimiter=',', names=True,
#                             encoding='utf-8-sig')
# # at some point would be good to build this as a pandas df instead of np array
# rhoTable = np.array([atmdata_raw['alt']/1e3,atmdata_raw['density']])
# params.atmdat = rhoTable


### GET ATM TABLE FROM BINARY EARTHGRAM DATA FILE
filename = '../data/rawOutput.txt'
# get Nmc atmosphere profiles
Nmc = 1
i_trial = 0
densPert, densMean, h = getMCdens(filename, Nmc)
# at some point would be good to build this as a pandas df instead of np array
rhoTable = np.array([h,densPert[:,i_trial]])
params.atmdat = rhoTable

### VEHICLE PARAMS
params.m = 2000 # kg 
params.A = 15 # m^2
params.CL = 0.1212
params.BC = 110
# get CD from BC
params.CD = params.m / (params.BC * params.A)

### INITIAL STATE
params.x0 = np.array([6478.100, 0, 0]) # inertial frame and coordinates
params.v0 = np.array([-1.25145758, 0.47238996, 7.9015096])

### CONTROL STATE
params.bank = 0 # full-lift-up, deg

### TIME VECTOR AND EXIT CONDITIONS
# should always stop on an exit condition
tspan = np.linspace(0,10000) # don't make too long or results get choppy!

# exit conditions:
params.hmin = 10
params.hmax = 125

event1 = lambda t, y: ODEs.below_min_alt(t, y, params)
event1.terminal = True

event2 = lambda t, y: ODEs.above_max_alt(t, y, params)
event2.terminal = True

events = (event1, event2)

### CALL SIMRUN
rvec_N, vvec_N = simRun(params, tspan, events)

### PLOTS
alt = np.linalg.norm(rvec_N, axis=0) - radE #/1e3
vmag = np.linalg.norm(vvec_N, axis=0)

fig = plt.figure(2)
ax = fig.add_subplot(111)
ax.plot(vmag, alt)
ax.set_xlabel('Inertial velocity, km/s')
ax.set_ylabel('Altitude, km')
ax.grid()

toc = time.time()
print('Time elapsed: %.2f s' % (toc-tic))