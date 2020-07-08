# -*- coding: utf-8 -*-
"""
Created on Tue Jul  7 17:05:31 2020

@author: Samuel Albert
"""

import numpy as np
import matplotlib.pyplot as plt

from sim import Params, Outs
import constants

plt.close('all')

## Load file archive and get data
filename = './../data/sweeps/FAKE_Earth_11_0_0_0708163447.npz'
data = np.load(filename, allow_pickle=True)
params = data['params'][0] # array of 1
outsList = data['outsList']
efpaList = data['efpaList']
BCList = data['BCList']

# recreate BClist and EFPAlist from outsList
# BClist = [outs.BC for outs in outsList]
# EFPAlist = [outs.efpa0 for outs in outsList]

## Create mesh grid for contour plots, reshape result arrays
BCgrid, EFPAgrid = np.meshgrid(BCList, efpaList)

ggrid = np.reshape([out.gpeak for out in outsList], BCgrid.shape)
hafgrid = np.reshape([out.haf for out in outsList], BCgrid.shape)
qgrid = np.reshape([out.qpeak for out in outsList], BCgrid.shape)
QLgrid = np.reshape([out.Qload for out in outsList], BCgrid.shape)
fpafgrid = np.reshape([out.fpaf for out in outsList], BCgrid.shape)

##  Data pruning
hafgrid[hafgrid < 0] = np.inf # set haf to inf for hyperbolic cases
hafgrid[fpafgrid < 0] = np.nan # ignore haf values for landed cases

# find line between landing and aerocapture
landline_BC = BCList
landline_EFPA = []
for rind in range(fpafgrid.shape[1]):
    ind = (next(ind for ind, val in enumerate(fpafgrid[:,rind])\
                              if val < 0))
    landline_EFPA.append(efpaList[ind])

## Make countour plots
fig = plt.figure()
ax = fig.add_subplot(111)
haflevels = [1e3, 10e3, 30e3, 100e3]
cp = ax.contour(EFPAgrid.T, BCgrid.T, hafgrid.T, haflevels)
ax.clabel(cp, inline = True, fontsize = 10, fmt = '%1.0f')

ax.plot(landline_EFPA, landline_BC)