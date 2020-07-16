# -*- coding: utf-8 -*-
"""
Created on Wed Jun 17 15:33:58 2020

@author: Samuel Albert
"""

import numpy as np

class p:
    def __init__(self, name):
        self.name = name

# =============================================================================
# PLANETARY CONSTANTS (from Vallado, in km, kg, s, rad)
# =============================================================================
EARTH = p('Earth')
EARTH.rad = 6378.1363 # equatorial radius, km
EARTH.mu = 3.986004415e5 # gravitational parameter, km^3/s^2
EARTH.om = 2 * np.pi / (0.99726968 * 86400) # rotation rate, rad/s
EARTH.k = 1.748e-04 # S-G coefficient (see spreadsheet), kg^0.5/m
EARTH.halt = 125 # atmospheric interface altitude, km

MARS = p('Mars')
MARS.rad = 3397.2 # equatorial radius, km
MARS.mu = 4.305e4 # gravitational parameter, km^3/s^2
MARS.om = 2 * np.pi / (1.02595675 * 86400) # rotation rate, rad/s
MARS.k = 1.904e-04 # S-G coefficient (see spreadsheet), kg^0.5/m
MARS.halt = 125 # atmospheric interface altitude, km

VENUS = p('Venus')
VENUS.rad = 6052.0 # equatorial radius, km
VENUS.mu = 3.257e5 # gravitational parameter, km^3/s^2
VENUS.om = 2 * np.pi / (-243 * 86400) # rotation rate, rad/s
VENUS.k = 1.897E-04 # S-G coefficient (see spreadsheet), kg^0.5/m
VENUS.halt = 135 # atmospheric interface altitude, km

TITAN = p('Titan')
TITAN.rad = 2574.73 # mean radius, km
TITAN.mu = 8978.1382 # gravitational parameter, km^3/s^2
TITAN.om = 2 * np.pi / (15.945 * 86400) # rotation rate, rad/s
TITAN.k = 1.758E-04
TITAN.halt = 800

G0 = 9.80665 # Earth standard gravity, m/s^2
