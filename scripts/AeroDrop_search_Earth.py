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

import planetaryconstants as constants
import ODEs
from sim import Params, Outs, mainAD

tic = time.time()
plt.close('all')

### CREATE params INPUT CLASS
params = Params()
params.p = constants.EARTH
params.returnTimeVectors = False

# =============================================================================
# comment/uncomment the below blocks of code for desired scenario at Mars
# =============================================================================
#### Make all BC and efpa ranges the same ####
BCList = np.arange(10, 200, 2.5)
efpaList = np.arange(-3.7, -7.4, -0.02)
# efpaList = np.arange(-4.4, -8.2, -0.02) # for 12 km/s entry



# ### Lift-up, nominal atmosphere
# params.atmMod = 'nom'
# params.LD = 0.25
# params.bank = 0 # deg
# modestr = 'lift-up, nominal atmosphere'
# ###

# ### Lift-up, 3sigHigh atmosphere
# params.atmMod = '3sigHigh'
# params.LD = 0.25
# params.bank = 0 # deg
# modestr = 'Lift-up, 3sigHigh atmosphere'
# ###

# ### Lift-up, 3sigLow atmosphere
# params.atmMod = '3sigLow'
# params.LD = 0.25
# params.bank = 0 # deg
# modestr = 'Lift-up, 3sigLow atmosphere'
# ###

# ### Lift-down, nominal atmosphere
# params.atmMod = 'nom'
# params.LD = 0.25
# params.bank = 180 # deg
# modestr = 'Lift-down, nominal atmosphere'
# ###

# ### Lift-down, 3sigHigh atmosphere
# params.atmMod = '3sigHigh'
# params.LD = 0.25
# params.bank = 180 # deg
# modestr = 'Lift-down, 3sigHigh atmosphere'
# ###

# ### Lift-down, 3sigLow atmosphere
# params.atmMod = '3sigLow'
# params.LD = 0.25
# params.bank = 180 # deg
# modestr = 'Lift-down, 3sigLow atmosphere'
# ###

# ### Ballistic, nominal atmosphere
# params.atmMod = 'nom'
# params.LD = 0
# params.bank = 0
# modestr = 'Ballistic, nominal atmosphere'
# ###

# ### Ballistic, 3sigHigh atmosphere
# params.atmMod = '3sigHigh'
# params.LD = 0
# params.bank = 0
# modestr = 'Ballistic, 3sigHigh atmosphere'
# ###

### Ballistic, 3sigLow atmosphere
params.atmMod = '3sigLow'
params.LD = 0
params.bank = 0
modestr = 'Ballistic, 3sigLow atmosphere'
###

# =============================================================================
# =============================================================================
### INPUT ATM TABLE - GET ATM TABLE FROM EARTHGRAM DATA FILE
params.dMode = 'table'
filename = '../data/dat_raw_Earth_nom.txt'
atmdata = np.genfromtxt(filename, names=True)
atmdata.sort(order='Hgtkm') # put in ascending altitude order
params.atmdat = np.array([atmdata['Hgtkm'], atmdata['DensMean']])

# alter altitude if requested
if params.atmMod == 'nom':
    pass
elif params.atmMod == '20% high':
    params.atmdat[1,:] = 1.2 * params.atmdat[1,:]
elif params.atmMod == '20% low':
    params.atmdat[1,:] = 0.8 * params.atmdat[1,:]
elif params.atmMod == '3sigHigh':
    params.atmdat[1,:] = params.atmdat[1,:] + \
        3 * (atmdata['DensMean'] * atmdata['SDden']/100)
elif params.atmMod == '3sigLow':
    params.atmdat[1,:] = params.atmdat[1,:] - \
        3 * (atmdata['DensMean'] * atmdata['SDden']/100)
else:
    sys.exit('atmMod not recognized')

### WIND-RELATIVE INITIAL STATE (COMPONENTS NOT CHANGED DURING GRID SEARCH)
params.inputType = 'wind-relative angles'
params.lat = 40
params.lon = 100
params.alt = params.p.halt
params.hdaWR = 0
params.vmagWR = 11

### ASSUME Rn = 1 m
params.Rn = 1

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
outsList = []

for params.efpaWR in efpaList:
    for params.BC in BCList:
        outs = Outs() # create new blank outs class instance
        outsList.append(mainAD(params, tspan, events, outs))
        print(modestr)

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








