# -*- coding: utf-8 -*-
"""
Created on Tue Feb  9 15:10:34 2021

@author: Samuel Albert
"""

from conversions import LLAEHV2RV, sph2cart, RV2LLAEHV, cart2sph
from sim import Params
import planetaryconstants as constants

import numpy as np

params = Params()
params.p = constants.MARS
# params.p.om = 0

# create spherical state
r = 4000 # km
lon = 15 # deg
lat = 100 # deg
v = 6 # km/s
fpa = -7 # deg
hda = 30 # deg

t = 0

# params.p.om = 0

# LLAEHV2RV
alt = r - params.p.rad
rvec1, vvec1 = LLAEHV2RV(lat, lon, alt, fpa, hda, v, params, t)
lat1, lon1, alt1, fpa1, hda1, vmag1 = RV2LLAEHV(rvec1, vvec1, params, t)

# sph2cart
xsphvec = np.array([r, np.radians(lon), np.radians(lat), v, np.radians(fpa),
                    np.radians(hda + 90)])
rvec2, vvec2 = sph2cart(xsphvec)
xsphvec2 = cart2sph(rvec2, vvec2)

print(rvec1)
print(rvec2)
print()
print(vvec1)
print(vvec2)

