# -*- coding: utf-8 -*-
"""
Created on Wed Jun 17 15:25:31 2020

@author: Samuel Albert
"""

import numpy as np
import time
import matplotlib.pyplot as plt
import sys

import constants
import ODEs
from sim import Params, Outs, mainAD

tic = time.time()
plt.close('all')

### CREATE params INPUT CLASS
params = Params()
params.p = constants.NEPTUNE
params.returnTimeVectors = True
params.atmMod = '20% low'

### INPUT ATM TABLE - GET ATM TABLE FROM EARTHGRAM DATA FILE
params.dMode = 'table'
filename = '../data/dens_Neptune_nom.txt'
atmdata = np.genfromtxt(filename, names=True)
atmdata.sort(order='Var_X') # put in ascending altitude order
params.atmdat = np.array([atmdata['Var_X'], atmdata['DENSAV']])

# alter altitude if requested
if params.atmMod == 'nom':
    pass
elif params.atmMod == '20% high':
    params.atmdat[1,:] = 1.2 * params.atmdat[1,:]
elif params.atmMod == '20% low':
    params.atmdat[1,:] = 0.8 * params.atmdat[1,:]
else:
    sys.exit('atmMod not recognized')

### VEHICLE PARAMS (NOT CHANGED DURING GRID SEARCH)
params.m = 2920 # kg, roughly MSL mass
params.CD = params.m / (115 * np.pi * (4.5/2)**2) # roughly MSL CD

params.LD = 0.25

### WIND-RELATIVE INITIAL STATE (COMPONENTS NOT CHANGED DURING GRID SEARCH)
params.inputType = 'wind-relative angles'
params.lat = 22.0
params.lon = 48.0
params.alt = params.p.halt
params.hdaWR = 0
params.vmagWR = 27

### CONTROL STATE
params.bank = 0 # deg

### TIME VECTOR AND EXIT CONDITIONS
# should always stop on an exit condition
tspan = (0,30000) # don't make too long or results get choppy!

# exit conditions:
params.hmin = 125
params.hmax = params.p.halt + 1e-7 + 10

event1 = lambda t, y: ODEs.below_min_alt(t, y, params)
event1.terminal = True

event2 = lambda t, y: ODEs.above_max_alt(t, y, params)
event2.terminal = True

events = (event1, event2)

### SINGLE RUN
params.efpaWR = -14.5
params.BC = 200



outs = Outs()
outs = mainAD(params, tspan, events, outs)
print(outs.fpaf)
print(outs.engf)



### PLOTS
alt = np.linalg.norm(outs.rvec_N, axis=0) - params.p.rad #/1e3
vmag = np.linalg.norm(outs.vvec_N, axis=0)

fig = plt.figure(1)
ax = fig.add_subplot(111)
ax.plot(vmag, alt)
ax.set_xlabel('Inertial velocity, km/s')
ax.set_ylabel('Altitude, km')
ax.grid()

fig = plt.figure(2)
ax = fig.add_subplot(111)
ax.plot(outs.tvec, outs.q)
ax.set_xlabel('time')
ax.set_ylabel('heat rate')
ax.grid()

fig = plt.figure(3)
ax = fig.add_subplot(111)
ax.plot(outs.tvec, vmag)
ax.set_xlabel('time')
ax.set_ylabel('vmag')
ax.grid()

fig = plt.figure(4)
ax = fig.add_subplot(111)
ax.plot(outs.tvec, alt)
ax.set_xlabel('time')
ax.set_ylabel('altitude, km')
ax.grid()

toc = time.time()
print('Time elapsed: %.2f s' % (toc-tic))








