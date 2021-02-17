# -*- coding: utf-8 -*-
"""
dens_txt2npz_Mars_Perts.py:
    takes MarsGRAM .txt file input, loads density and saves data in .npz file.
    Saves density perturbations to file instead of actual density.
    Run this script just once per MarsGRAM batch output.
    
Created on Tues Feb 16 15:37:27 2020

@author: Samuel Albert
"""


import numpy as np

from atm import getMarsGRAMDensTableAll

filename = '../data/Mars_0.1_50000.txt'

# get density from GRAM output file
Nmc = 10657 # note: there are actually only this many profiles in the 50,000 case GRAM data file...
densTot, densMean, h = getMarsGRAMDensTableAll(filename, Nmc)

# get density perturbations
for i in range(Nmc):
    densTot[:,i] = densTot[:,i] / densMean - 1

# save in binary to .npz file
outname = '../data/Mars_0.1_50000_Perts.npz'

np.savez(outname,
         densTot = densTot,
         densMean = densMean,
         h = h)

