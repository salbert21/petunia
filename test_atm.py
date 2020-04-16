# -*- coding: utf-8 -*-
"""
Created on Thu Apr 16 11:45:54 2020

Tests new functionality for atm module interacting with Earth-GRAM2016 data
for many Monte Carlo runs.

@author: Samuel
"""

import numpy as np
from atm import getMCAtmdat, getRho_from_EarthGRAM

filename = 'data/rawOutput.txt'

atmList = getMCAtmdat(filename, 3)

hlist = np.array([125, 80.7, 100.3])
rholist = getRho_from_EarthGRAM(atmList[0], hlist)