# -*- coding: utf-8 -*-
"""
Created on Tue Jul  7 17:05:31 2020

@author: Samuel Albert
"""

import numpy as np

class Params:
       
    class p:
        pass 
    
class Outs:
    
    class s0:
        pass
    class sf:
        pass
    class v:
        pass

# load file archive and get data
filename = './../data/sweeps/FAKEEarth_11_0.0_0_0707174059.npz'
data = np.load(filename, allow_pickle=True)
params = data['params']
outsList = data['outsList']

# recreate BClist and EFPAlist from outsList
# BClist = [outs.BC for outs in outsList]
# EFPAlist = [outs.efpa0 for outs in outsList]

