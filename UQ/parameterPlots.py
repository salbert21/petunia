# -*- coding: utf-8 -*-
"""
Created on Thu Oct  8 17:51:35 2020

parameterPlots.py:
    runs a few representative EDL scenarios, then plots a range of parameters

@author: Samuel Albert
"""

import numpy as np
import time
import matplotlib.pyplot as plt
import matplotlib as mpl

import constants
import ODEs
from sim import Params, Outs, mainAD
from atm import getMarsGRAMDensTable
from conversions import RV2LLAEHV, VN2Vinf

tic = time.time()
plt.close('all')
mpl.rcParams["figure.autolayout"] = True
mpl.rcParams.update({'font.size': 16})

outsList = []


# =============================================================================
# Create Params input class for Mars
# =============================================================================
params = Params()
params.p = constants.MARS
params.returnTimeVectors = True

# =============================================================================
# Generate atm table for Nmc profiles
# =============================================================================
filename = './../data/Mars_0.1_5000.txt'
Nmc = 1
densAll, densMean, h = getMarsGRAMDensTable(filename, Nmc)
# densAll, densMean, h = getMarsGRAMDensTableAll(filename, Nmc)
params.dMode = 'table'

# =============================================================================
# Load atm table for this run
# =============================================================================
i_trial = 0
params.atmdat = np.array([h,densAll[:,i_trial]])
# make sure atmdata is sorted by ascending altitude
params.atmdat = params.atmdat[:,params.atmdat[0,:].argsort()]

# =============================================================================
# Set vehicle params
# =============================================================================
params.m = 3000 # kg
params.CD = 1.59
params.LD = 0 # L/D = 0 --> CL= 0
params.BC = 120 # kg/m^2

# =============================================================================
# Wind-Relative initial state
# =============================================================================
params.inputType = 'wind-relative angles'
params.lat = 0
params.lon = 0
params.alt = params.p.halt
params.hdaWR = 0
params.vmagWR = 6 # km/s

# =============================================================================
# Control state
# =============================================================================
params.bank = 0

# =============================================================================
# Time vector and exit conditions
# =============================================================================
tspan = (0,30000)

params.hmin = 20
params.hmax = params.p.halt + 10

event1 = lambda t, y: ODEs.below_min_alt(t, y, params)
event1.terminal = True

event2 = lambda t, y: ODEs.above_max_alt(t, y, params)
event2.terminal = True

events = (event1, event2)


# =============================================================================
# Run sim for a few values of EFPA
# =============================================================================
efpaList = np.around(np.arange(-10.3, -11.5, -0.025), 7)
for efpai in efpaList:
    params.efpaWR = efpai

    outs = Outs()
    outs = mainAD(params, tspan, events, outs)
    
    # get EFPA at each time step
    outs.latvec = np.empty(len(outs.tvec))
    outs.latvec[:] = np.NaN
    outs.lonvec = np.empty(len(outs.tvec))
    outs.lonvec[:] = np.NaN
    outs.altvec = np.empty(len(outs.tvec))
    outs.altvec[:] = np.NaN
    outs.efpaWRvec = np.empty(len(outs.tvec))
    outs.efpaWRvec[:] = np.NaN
    outs.hdaWRvec = np.empty(len(outs.tvec))
    outs.hdaWRvec[:] = np.NaN
    outs.vmagWRvec = np.empty(len(outs.tvec))
    outs.vmagWRvec[:] = np.NaN
    outs.engvec = np.empty(len(outs.tvec))
    outs.engvec[:] = np.NaN
    
    for ind, ti in enumerate(outs.tvec):
        outs.vInfvec_N = VN2Vinf(outs.rvec_N[:,ind], outs.vvec_N[:,ind],
                                   params, outs.tvec[ind])
        outs.latvec[ind], outs.lonvec[ind], outs.altvec[ind],\
            outs.efpaWRvec[ind], outs.hdaWRvec[ind], outs.vmagWRvec[ind] = \
            RV2LLAEHV(outs.rvec_N[:,ind], outs.vvec_N[:,ind],
                      params, outs.tvec[ind])
            
        outs.engvec[ind] = np.linalg.norm(outs.vvec_N[:,ind])**2 / 2\
            - params.p.mu / np.linalg.norm(outs.rvec_N[:,ind])
            
            
    outs.engInt = np.trapz(outs.engvec, outs.tvec)
    
    # get final slope of vmagWR
    outs.slf = (outs.vmagWRvec[-1] - outs.vmagWRvec[-6])\
        / (outs.tvec[-1] - outs.tvec[-6])
    

        
    
    outsList.append(outs)

# =============================================================================
# Plot results
# =============================================================================

NUM_COLORS = len(efpaList)
LEGENDS_ON = False

yline = 5

cm = plt.get_cmap('gist_rainbow')

fig1 = plt.figure(1)
ax1 = fig1.add_subplot()
ax1.set_prop_cycle(color = [cm(1.*i/NUM_COLORS) for i in range(NUM_COLORS)])

fig2 = plt.figure(2)
ax2 = fig2.add_subplot()
ax2.set_prop_cycle(color = [cm(1.*i/NUM_COLORS) for i in range(NUM_COLORS)])

fig3 = plt.figure(3)
ax3 = fig3.add_subplot()
ax3.set_prop_cycle(color = [cm(1.*i/NUM_COLORS) for i in range(NUM_COLORS)])

fig4 = plt.figure(4)
ax4 = fig4.add_subplot()
ax4.set_prop_cycle(color = [cm(1.*i/NUM_COLORS) for i in range(NUM_COLORS)])

fig5 = plt.figure(5)
ax5 = fig5.add_subplot()
ax5.set_prop_cycle(color = [cm(1.*i/NUM_COLORS) for i in range(NUM_COLORS)])

fig6 = plt.figure(6)
ax6 = fig6.add_subplot()
ax6.set_prop_cycle(color = [cm(1.*i/NUM_COLORS) for i in range(NUM_COLORS)])

fig7 = plt.figure(7)
ax7 = fig7.add_subplot()
ax7.set_prop_cycle(color = [cm(1.*i/NUM_COLORS) for i in range(NUM_COLORS)])

fig8 = plt.figure(8)
ax8 = fig8.add_subplot()
ax8.set_prop_cycle(color = [cm(1.*i/NUM_COLORS) for i in range(NUM_COLORS)])

fig9 = plt.figure(9)
ax9 = fig9.add_subplot()
ax9.set_prop_cycle(color = [cm(1.*i/NUM_COLORS) for i in range(NUM_COLORS)])

for ind, outs in enumerate(outsList):
    
    ## ALT vs. VMAG
    ax1.plot(outs.vmagWRvec, outs.altvec,
             label = 'EFPA = {}'.format(efpaList[ind]))
    
    ## FPA vs. TIME
    ax2.plot(outs.tvec, outs.efpaWRvec,
             label = 'EFPA = {}'.format(efpaList[ind]))
    
    ## ALT vs. TIME
    ax3.plot(outs.tvec, outs.altvec,
             label = 'EFPA = {}'.format(efpaList[ind]))
    
    ## VMAG vs. TIME
    ax4.plot(outs.tvec, outs.vmagWRvec,
             label = 'EFPA = {}'.format(efpaList[ind]))
    
    ## q vs. TIME
    ax5.plot(outs.tvec, outs.q,
             label = 'EFPA = {}'.format(efpaList[ind]))
    
    ## gload vs. TIME
    ax6.plot(outs.tvec, outs.gload,
             label = 'EFPA = {}'.format(efpaList[ind]))
    
    ## eng vs. TIME
    ax7.plot(outs.tvec, outs.engvec,
             label = 'EFPA = {}'.format(efpaList[ind]))
    
    ## engInt number line
    ax8.plot(outs.engInt, yline, 'o',
             label = 'EFPA = {}'.format(efpaList[ind]))
    
    ax9.plot(outs.slf, '.')

    
    
## ALT vs. VMAG
ax1.set_xlabel('Inertial velocity, km/s')
ax1.set_ylabel('Altitude, km')
ax1.grid()

## FPA vs. TIME
ax2.set_xlabel('Time, s')
ax2.set_ylabel('Flight path angle, deg')
ax2.grid()

## ALT vs. TIME
ax3.set_xlabel('Time, s')
ax3.set_ylabel('Altitude, km')
ax3.grid()

## VMAG vs. TIME
ax4.set_xlabel('Time, s')
ax4.set_ylabel('Airspeed, km/s')
ax4.grid()

## q vs. TIME
ax5.set_xlabel('Time, s')
ax5.set_ylabel('Heating rate, W/cm^2')
ax5.grid()

## gload vs. TIME
ax6.set_xlabel('Time, s')
ax6.set_ylabel('G-load, Earth g\'s')
ax6.grid()

## eng vs. TIME
ax7.set_xlabel('Time, s')
ax7.set_ylabel(r'Specific orbital energy, $km^2 s^{-2}$')
ax7.grid()

## engInt vs. TIME
ax8.set_xlabel('Integral of specific mechanic energy, km^2/s')
ax8.grid()

## final vmag slope vs. time
ax9.set_ylabel('final airspeed slope')
ax9.grid()



if LEGENDS_ON:
    ax1.legend()
    ax2.legend()
    ax3.legend()
    ax4.legend()
    ax5.legend()
    ax6.legend()
    ax7.legend()
    ax8.legend()



# fig1.savefig('alt_vs_vmag.png')
# fig2.savefig('fpa_vs_time.png')
# fig3.savefig('alt_vs_time.png')
# fig4.savefig('v_vs_time.png')
# fig5.savefig('q_vs_time.png')
# fig6.savefig('g_vs_time.png')
# fig7.savefig('eng_vs_time.png')












