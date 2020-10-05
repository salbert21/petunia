# -*- coding: utf-8 -*-
"""
Created on Thu Apr 16 11:45:54 2020

Tests new functionality for atm module interacting with Earth-GRAM2016 data
for many Monte Carlo runs.

@author: Samuel
"""

import numpy as np
import matplotlib.pyplot as plt
from atm import getMCdens

plt.close('all')
plt.rcParams.update({'font.size': 15})

filename = '../data/rawOutput.txt'

# get Nmc atmosphere profiles
Nmc = 10
densPert, densMean, h = getMCdens(filename, Nmc)

# compare with exponential model
rho0 = 1.225
H = 7.257 # km
densExp = rho0 * np.exp(-h/H)
MeanErr = (densExp - densMean) / densExp

PertErr = (densExp[:,None] - densPert) / densExp[:,None]

fig = plt.figure()
ax = fig.add_subplot(111)
ax.plot(densMean, h, label='GRAM Mean')
ax.plot(densExp, h, label='Exponential')
ax.set_xscale('log')
ax.set_xlabel('Density, kg/m^2')
ax.set_ylabel('Altitude, km')
ax.legend()
ax.grid()

fig = plt.figure()
ax = fig.add_subplot(111)
ax.plot(MeanErr, h, label='GRAM mean vs. exponential')
# ax.plot(PertErr, h, label='perturbed GRAM vs exponential')
# ax.legend()
ax.grid()
ax.set_xlabel('Relative density error normalized by exponential model')
ax.set_ylabel('Altitude, km')


# # test functionality of interpolating to get density at specific altitudes
# hlist = np.array([125, 80.7, 100.3])
# rholist = getRho_from_EarthGRAM(atmList[0], hlist)

# # plot some altitude profiles, include only data above 50 km
# fig = plt.figure(1)
# ax = fig.add_subplot(111)
# ax.grid()
# for i in range(Nmc):
#     dpert = (densMean - atmList[i]['DensPert']) / densMean
#     ax.plot(dpert.where(atmList[i]['Hgtkm'] > 50),
#             atmList[i].where(atmList[i]['Hgtkm'] > 50)['Hgtkm'],'b')
#     # print('Plotted MC profile %d\n' %i)