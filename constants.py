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

MARS = p('Mars')
MARS.rad = 3397.2 # equatorial radius, km
MARS.mu = 4.305e4 # gravitational parameter, km^3/s^2
MARS.om = 2 * np.pi / (1.02595675 * 86400) # rotation rate, rad/s
