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
params.p = constants.MARS
params.returnTimeVectors = True
params.atmMod = 'nom'

### INPUT ATM TABLE - GET ATM TABLE FROM GRAM DATA FILE
params.dMode = 'table'
filename = '../data/dens_Mars_nom.txt'
# filename = '../data/dat_raw_Earth_nom.txt'
atmdata = np.genfromtxt(filename, names=True)
# atmdata.sort(order='Var_X') # put in ascending altitude order
params.atmdat = np.array([atmdata['Var_X'], atmdata['DENSAV']])
# atmdata.sort(order='Hgtkm') # put in ascending altitude order
# params.atmdat = np.array([atmdata['Hgtkm'], atmdata['DensMean']])

# alter density if requested
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
params.CD = 1.6 # roughly MSL CD

params.LD = 0.25

### WIND-RELATIVE INITIAL STATE (COMPONENTS NOT CHANGED DURING GRID SEARCH)
params.inputType = 'wind-relative angles'
params.lat = 18.38
params.lon = -77.58
params.alt = params.p.halt
params.hdaWR = 0
params.vmagWR = 6

### CONTROL STATE
params.bank = 30 # deg

### TIME VECTOR AND EXIT CONDITIONS
# should always stop on an exit condition
tspan = (0,30000)

# exit conditions:
params.hmin = 10
params.hmax = params.p.halt + 1e-7 + 10

event1 = lambda t, y: ODEs.below_min_alt(t, y, params)
event1.terminal = True

event2 = lambda t, y: ODEs.above_max_alt(t, y, params)
event2.terminal = True

events = (event1, event2)

### SINGLE RUN
params.efpaWR = -12
params.BC = 130

### SIMULATE ORBITER
outsOrbiter = Outs()
outsOrbiter = mainAD(params, tspan, events, outsOrbiter)
print('\nORBITER:')
print('Final energy: {:.3f} km^2/s^2'.format(outsOrbiter.engf))
print('Apoapsis altitude: {:.3f} km'.format(outsOrbiter.haf))
print('Bank angle set to {:.1f} deg\n'.format(params.bank))

### SIMULATE PROBE
params.LD = 0
outsProbe = Outs()
outsProbe = mainAD(params, tspan, events, outsProbe)
print('\nPROBE:')
print('Final energy: {:.3f} km^2/s^2'.format(outsProbe.engf))
print('Apoapsis altitude: {:.3f} km'.format(outsProbe.haf))


### PLOTS
altOrbiter = np.linalg.norm(outsOrbiter.rvec_N, axis=0) - params.p.rad #/1e3
vmagOrbiter = np.linalg.norm(outsOrbiter.vvec_N, axis=0)

altProbe = np.linalg.norm(outsProbe.rvec_N, axis=0) - params.p.rad #/1e3
vmagProbe = np.linalg.norm(outsProbe.vvec_N, axis=0)

fig = plt.figure(1)
ax = fig.add_subplot(111)
ax.plot(vmagOrbiter, altOrbiter, label = 'Orbiter')
ax.plot(vmagProbe, altProbe, label = 'Probe')
ax.set_xlabel('Inertial velocity, km/s')
ax.set_ylabel('Altitude, km')
ax.grid()
ax.legend()

fig = plt.figure(2)
ax = fig.add_subplot(111)
ax.plot(outsOrbiter.tvec, outsOrbiter.q, label = 'Orbiter')
ax.plot(outsProbe.tvec, outsProbe.q, label = 'Probe')
ax.set_xlabel('time')
ax.set_ylabel('heat rate')
ax.grid()
ax.legend()

fig = plt.figure(3)
ax = fig.add_subplot(111)
ax.plot(outsOrbiter.tvec, vmagOrbiter, label = 'Orbiter')
ax.plot(outsProbe.tvec, vmagProbe, label = 'Probe')
ax.set_xlabel('time')
ax.set_ylabel('vmag')
ax.grid()
ax.legend()

fig = plt.figure(4)
ax = fig.add_subplot(111)
ax.plot(outsOrbiter.tvec, altOrbiter, label = 'Orbiter')
ax.plot(outsProbe.tvec, altProbe, label = 'Probe')
ax.set_xlabel('time')
ax.set_ylabel('altitude, km')
ax.grid()
ax.legend()

toc = time.time()
print('Time elapsed: %.2f s' % (toc-tic))








