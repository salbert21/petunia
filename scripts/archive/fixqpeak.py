# -*- coding: utf-8 -*-
"""
Created on Wed Jul 22 14:34:24 2020

@author: Samuel Albert
"""

import numpy as np

from sim import Params, Outs
import constants

# fixing qpeak values that got screwed up for Earth

## Load file archive and get data
filename = './../results/sweeps/Earth_12_nom_0_0_0710033611.npz'
data = np.load(filename, allow_pickle=True)
params = data['params'][0] # array of 1
outsList = data['outsList']
efpaList = data['efpaList']
BCList = data['BCList']

for i, outi in enumerate(outsList):
    outi.qpeak = params.p.k * outi.SGpeak * 1e5
    
np.savez(filename,
          params = np.array([params]),
          outsList = outsList,
          efpaList = efpaList,
          BCList = BCList
          )