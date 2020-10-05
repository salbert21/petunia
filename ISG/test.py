# -*- coding: utf-8 -*-
"""
Created on Mon Oct  5 11:00:13 2020

test.py: testing monte carlo functionality for atmposheric functions

@author: Samuel Albert
"""

from atm import getMarsGRAMDensTable


## Test Mars atmosphere density file load
filename = 'data/Mars_0.1_5000.txt'

# get Nmc atmosphere profiles
Nmc = 10
densPert, densMean, h = getMarsGRAMDensTable(filename, Nmc)