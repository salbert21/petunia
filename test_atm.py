# -*- coding: utf-8 -*-
"""
Created on Thu Apr 16 11:45:54 2020

Tests new functionality for atm module interacting with Earth-GRAM2016 data
for many Monte Carlo runs.

@author: Samuel
"""

import numpy as np
import matplotlib.pyplot as plt
from atm import getMCAtmdat, getMCdens, getRho_from_EarthGRAM

plt.close('all')

filename = 'data/rawOutput.txt'

# get 3 atmosphere profiles
Nmc = 10
densPert, densMean, h = getMCdens(filename, Nmc)


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