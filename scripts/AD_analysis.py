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
filename = './../results/sweeps/Titan_6_nom_0_0_0716043817.npz'
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
engfgrid = np.reshape([out.engf for out in outsList], BCgrid.shape)

##  Data pruning
hafgrid[hafgrid < 0] = np.inf # set haf to inf for hyperbolic cases
hafgrid[fpafgrid < 0] = np.nan # ignore haf values for landed cases

# find line between landing and aerocapture
landline_BC = BCList
landline_EFPA = []
for rind in range(fpafgrid.shape[1]):
    ind = next(ind for ind, val in enumerate(fpafgrid[:,rind])\
               if val < 0)
    landline_EFPA.append(efpaList[ind])
    
# find line between aerocapture and escape
spaceline_BC = []
spaceline_EFPA = []
for rind in range(engfgrid.shape[1]):
    ind = [ind for ind, val in enumerate(engfgrid[:,rind]) if val < 0]
    if len(ind) == 0:
        break
    else:
        spaceline_BC.append(BCList[rind])
        spaceline_EFPA.append(efpaList[ind[0]])



## Make countour plots
fig = plt.figure()
ax = fig.add_subplot(111)
ax.grid()
haflevels = [5e3, 10e3, 30e3, 120e3]
cp = ax.contour(EFPAgrid.T, BCgrid.T, hafgrid.T, haflevels, colors = 'blue')
ax.clabel(cp, inline = True, colors = 'blue', fontsize = 10, fmt = '%1.0f')

ax.plot(landline_EFPA, landline_BC, color = 'black')
ax.plot(spaceline_EFPA, spaceline_BC, color = 'magenta')
ax.fill_between(landline_EFPA, landline_BC, 10,
                facecolor = 'blue', alpha = 0.3)

edgy = np.minimum(spaceline_EFPA, max(landline_EFPA))
ax.fill_betweenx(landline_BC, landline_EFPA, edgy,
                 facecolor = 'green', alpha = 0.3)
