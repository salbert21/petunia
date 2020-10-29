# -*- coding: utf-8 -*-
"""
dens_txt2npz_Mars.py:
    takes MarsGRAM .txt file input, loads density and saves data in .npz file.
    Run this script just once per MarsGRAM batch output.
Created on Tue Oct 20 14:27:57 2020

@author: Samuel Albert
"""


import numpy as np

from atm import getMarsGRAMDensTableAll

filename = '../data/Mars_0.1_50000.txt'

# get density from GRAM output file
Nmc = 5000
densTot, densMean, h = getMarsGRAMDensTableAll(filename, Nmc)

# save in binary to .npz file
outname = '../data/Mars_0.1_50000.npz'

np.savez(outname,
         densTot = densTot,
         densMean = densMean,
         h = h)

