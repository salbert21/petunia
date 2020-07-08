# -*- coding: utf-8 -*-
"""
Created on Tue Jul  7 17:05:31 2020

@author: Samuel Albert
"""

import numpy as np
import matplotlib.pyplot as plt

class Params:
       
    class p:
        pass 
    
class Outs:
    pass

plt.close('all')

## Load file archive and get data
filename = './../data/sweeps/Earth_11_0.0_0_0707195056.npz'
data = np.load(filename, allow_pickle=True)
params = data['params']
outsList = data['outsList']

# recreate BClist and EFPAlist from outsList
# BClist = [outs.BC for outs in outsList]
# EFPAlist = [outs.efpa0 for outs in outsList]

## Create mesh grid for contour plots, reshape result arrays
BCgrid, EFPAgrid = np.meshgrid(outsList[0].BCList, outsList[0].efpaList)

ggrid = np.reshape([out.gpeak for out in outsList], BCgrid.shape)
hafgrid = np.reshape([out.haf for out in outsList], BCgrid.shape)
qgrid = np.reshape([out.qpeak for out in outsList], BCgrid.shape)
QLgrid = np.reshape([out.Qload for out in outsList], BCgrid.shape)
fpafgrid = np.reshape([out.fpaf for out in outsList], BCgrid.shape)

# TODO: need to add some data pruning (i.e. get rid of negative haf values)

## Make countour plots
fig = plt.figure()
ax = fig.add_subplot(111)
cp = ax.contour(BCgrid, EFPAgrid, hafgrid)
ax.clabel(cp, inline = True, fontsize = 10, fmt = '%1.0f')