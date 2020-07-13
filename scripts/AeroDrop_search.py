# -*- coding: utf-8 -*-
"""
Created on Wed Jun 17 15:25:31 2020

@author: Samuel Albert
"""

import numpy as np
import time
import matplotlib.pyplot as plt
from datetime import datetime
import sys

import constants
import ODEs
from sim import Params, Outs, mainAD

tic = time.time()
plt.close('all')

### CREATE params INPUT CLASS
params = Params()
params.p = constants.VENUS
params.returnTimeVectors = False
params.atmMod = 'nom'

### INPUT ATM TABLE - GET ATM TABLE FROM EARTHGRAM DATA FILE
params.dMode = 'table'
filename = '../data/dens_Venus_nom.txt'
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

params.LD = 0 #0.25

### WIND-RELATIVE INITIAL STATE (COMPONENTS NOT CHANGED DURING GRID SEARCH)
params.inputType = 'wind-relative angles'
params.lat = -31.3
params.lon = 317.0
params.alt = params.p.halt
params.hdaWR = 0
params.vmagWR = 11.5

### CONTROL STATE
params.bank = 0 # deg

### TIME VECTOR AND EXIT CONDITIONS
# should always stop on an exit condition
tspan = (0,3000) # don't make too long or results get choppy!

# exit conditions:
params.hmin = 30
params.hmax = params.p.halt + 1e-7

event1 = lambda t, y: ODEs.below_min_alt(t, y, params)
event1.terminal = True

event2 = lambda t, y: ODEs.above_max_alt(t, y, params)
event2.terminal = True

events = (event1, event2)

### GRID SEARCH
efpaList = np.arange(-3.5, -7.3, -0.02)
BCList = np.arange(10, 200, 2.5)
outsList = []

for params.efpaWR in efpaList:
    for params.BC in BCList:
        outs = Outs() # create new blank outs class instance
        outsList.append(mainAD(params, tspan, events, outs))

## Save results to a file
outname = './../results/sweeps/' + params.p.name + '_' + str(params.vmagWR) + '_'\
        + params.atmMod + '_' + str(params.LD) + '_' + str(params.bank) + '_'\
        + datetime.now().strftime('%m%d%H%M%S')
    
np.savez(outname,
         params = np.array([params]),
         outsList = outsList,
         efpaList = efpaList,
         BCList = BCList
         )

toc = time.time()
print('%d trajectories simulated' %(len(efpaList) * len(BCList)))
print('Time elapsed: %.2f s' % (toc-tic))








