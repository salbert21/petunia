# -*- coding: utf-8 -*-
"""
Created on Tue Feb 23 16:12:25 2021

AeroDrop_dispersed_analysis.py:
    loads and analyzes compiled data generated by AeroDrop_dispersed.py

@author: Samuel Albert
"""

import planetaryconstants as constants
from sim import Params


import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
from scipy import stats

plt.close('all')
mpl.rcParams["figure.autolayout"] = True
mpl.rcParams.update({'font.size': 15})

# =============================================================================
# Load Mars constants
# =============================================================================
params = Params()
params.p = constants.MARS


# =============================================================================
# Load data
# =============================================================================
filename = '../results/AeroDrop_dispersed_500_0302193435.npz'
data = np.load(filename, allow_pickle = True)

xxvecArr_O = data['xxvecArr_O']
tvecArr_O = data['tvecArr_O']
sigdvecArr_O = data['sigdvecArr_O']
tsvecArr_O = data['tsvecArr_O']
BCList_O = data['BCList_O']
L_DList_O = data['L_DList_O']
gam0List = data['gam0List']
v0List = data['v0List']
raErrList_O = data['raErrList_O']
DVList_O = data['DVList_O']
tsfList_O = data['tsfList_O']
sigdfList_O = data['sigdfList_O']
rafList_O = data['rafList_O']
rpfList_O = data['rpfList_O']
engfList_O = data['engfList_O']
xxvecArr_P = data['xxvecArr_P']
tvecArr_P = data['tvecArr_P']
sigvecArr = data['sigvecArr']
sig0ListArr = data['sig0ListArr']
sfErrList_P = data['sfErrList_P']
hfErrList_P = data['hfErrList_P']
vfErrList_P = data['vfErrList_P']
xxvecArr_PBC = data['xxvecArr_PBC']
tvecArr_PBC = data['tvecArr_PBC']
sfErrList_PBC = data['sfErrList_PBC']
hfErrList_PBC = data['hfErrList_PBC']
vfErrList_PBC = data['vfErrList_PBC']
atmindList = data['atmindList']

del data



# =============================================================================
# Print statistics
# =============================================================================
print('\n\nORBITER STATISTICS:')
print('\nTotal Delta-V Cost (m/s):')
print(stats.describe(DVList_O))
print('Standard deviation: {0:.5f} m/s'.format(stats.tstd(DVList_O)))

print('\nApoapsis Radius Error (km):')
print(stats.describe(raErrList_O/1e3))
print('Standard deviation: {0:.5f} km'.format(stats.tstd(raErrList_O/1e3)))

print('\n\nGUIDED PROBE STATISTICS:')
print('\nRange Error (km):')
print(stats.describe(params.p.rad * np.radians(sfErrList_P)))
print('Standard deviation: {:0.5f} km'.\
      format(stats.tstd(params.p.rad * np.radians(sfErrList_P))))

print('\nAltitude Error (m):')
print(stats.describe(hfErrList_P))
print('Standard deviation: {0:.5f} m'.format(stats.tstd(hfErrList_P)))

print('\nVelocity Error (m/s):')
print(stats.describe(vfErrList_P))
print('Standard deviation: {0:.5f} m/s'.format(stats.tstd(vfErrList_P)))

print('\n\nPASSIVE PROBE STATISTICS:')
print('\nRange Error (km):')
print(stats.describe(params.p.rad * np.radians(sfErrList_PBC)))
print('Standard deviation: {0:.5f} km'.\
      format(stats.tstd(params.p.rad * np.radians(sfErrList_PBC))))

print('\nVelocity Error (m/s):')
print(stats.describe(vfErrList_PBC))
print('Standard deviation: {0:.5f} m/s'.format(stats.tstd(vfErrList_PBC)))



# =============================================================================
# Plots
# =============================================================================
Nbin = 'auto'
rwidth = 1

# Orbiter DV histogram
fig1 = plt.figure()
ax1 = fig1.add_subplot(111)
ax1.hist(DVList_O, bins = Nbin, rwidth = rwidth)
ax1.grid()
ax1.set_title('Total Delta-V Cost for Orbiter')
ax1.set_xlabel('Total Delta-V Cost, m/s')
ax1.set_ylabel('Occurrences')

# Orbiter apoapsis error histogram
fig3 = plt.figure()
ax3 = fig3.add_subplot(111)
ax3.hist(raErrList_O/1e3, bins = Nbin, rwidth = rwidth)
ax3.grid()
ax3.set_title('Apoapsis Radius Error for Orbiter')
ax3.set_xlabel('Apoapsis targeting error, km')
ax3.set_ylabel('Occurrnces')

# Lifting probe range error histogram
fig4 = plt.figure()
ax4 = fig4.add_subplot(111)
ax4.hist(params.p.rad * np.radians(sfErrList_P), bins = Nbin, rwidth = rwidth)
ax4.grid()
ax4.set_title('Range Error for Guided Probe')
ax4.set_xlabel('Range targeting error, km')
ax4.set_ylabel('Occurrences')

# Lifting probe altitude error histogram
fig6 = plt.figure()
ax6 = fig6.add_subplot(111)
ax6.hist(hfErrList_P, bins = Nbin, rwidth = rwidth)
ax6.grid()
ax6.set_title('Altitude Error for Guided Probe')
ax6.set_xlabel('Altitude targeting error, m')
ax6.set_ylabel('Occurrences')

# Lifting probe velocity error histogram
fig7 = plt.figure()
ax7 = fig7.add_subplot(111)
ax7.hist(vfErrList_P, bins = Nbin, rwidth = rwidth)
ax7.grid()
ax7.set_title('Velocity Error for Guided Probe')
ax7.set_xlabel('Velocity targeting error, m/s')
ax7.set_ylabel('Occurrences')

# Passive probe range error histogram
fig8 = plt.figure()
ax8 = fig8.add_subplot(111)
ax8.hist(params.p.rad * np.radians(sfErrList_PBC), bins = Nbin, rwidth = rwidth)
ax8.grid()
ax8.set_title('Range Error for Passive Probe')
ax8.set_xlabel('Range targeting error, km')
ax8.set_ylabel('Occurrences')

# NOTE: don't need this one because passive probe terminates on altitude, so
#           altitude error is always 0
# # Passive probe altitude error histogram
# fig9 = plt.figure()
# ax9 = fig9.add_subplot(111)
# ax9.hist(hfErrList_PBC, bins = Nbin)
# ax9.grid()
# ax9.set_title('Altitude Error for Passive Probe')
# ax9.set_xlabel('Altitude targeting error, m')
# ax9.set_ylabel('Occurrences')

# Passive probe velocity error histogram
fig10 = plt.figure()
ax10 = fig10.add_subplot(111)
ax10.hist(vfErrList_PBC, bins = Nbin, rwidth = rwidth)
ax10.grid()
ax10.set_title('Velocity Error for Passive Probe')
ax10.set_xlabel('Velocity targeting error, m/s')
ax10.set_ylabel('Occurrences')


# DV vs. EFPA
fig2 = plt.figure()
ax2 = fig2.add_subplot(111)
ax2.plot(np.degrees(gam0List), DVList_O, 'o')
ax2.grid()
ax2.set_title('Total Delta-V vs. EFPA for Orbiter')
ax2.set_xlabel('Planet-relative entry flight path angle, deg')
ax2.set_ylabel('Total Delta-V Cost, m/s')

# Orbiter apoapsis error vs EFPA
fig5 = plt.figure()
ax5 = fig5.add_subplot(111)
ax5.plot(np.degrees(gam0List), raErrList_O/1e3, 'o')
ax5.grid()
ax5.set_title('Apoapsis Radius Error vs. EFPA for Orbiter')
ax5.set_xlabel('Planet-relative entry flight path angle, deg')
ax5.set_ylabel('Apoapsis targeting error, km')

# Orbiter apoapsis error vs BC
fig11 = plt.figure()
ax11 = fig11.add_subplot(111)
ax11.plot(BCList_O, raErrList_O/1e3, 'o')
ax11.grid()
ax11.set_title('Apoapsis Radius Error vs. Ballistic Coefficient for Orbiter')
ax11.set_xlabel('Ballistic coefficient, kg/m^2')
ax11.set_ylabel('APoapsis targeting error, km')

# DV vs. apoapsis error
fig12 = plt.figure()
ax12 = fig12.add_subplot(111)
ax12.plot(raErrList_O/1e3, DVList_O, 'o')
ax12.grid()
ax12.set_title('Total delta-V cost vs. apoapsis radius error')
ax12.set_xlabel('Apoapsis targeting error, km')
ax12.set_ylabel('Total delta-V cost, m/s')

# DV vs. final periapsis
fig13 = plt.figure()
ax13 = fig13.add_subplot(111)
ax13.plot(rpfList_O/1e3 - params.p.rad, DVList_O, 'o')
ax13.grid()
ax13.set_title('Total delta-V cost vs. periapsis altitude')
ax13.set_xlabel('Periapsis altitude upon atm exit, km')
ax13.set_ylabel('Total delta-V cost, m/s')

























