# -*- coding: utf-8 -*-
"""
Created on Tue Feb  9 13:19:06 2021

@author: Samuel Albert
"""
from conversions import LLAEHV2RV, RV2LLAEHV, Vinf2VN, getApses, getApsesSphPR, sph2cart, cart2sph, VN2Vinf
from sim import Params
from guidance import updateFNPAG, dynFNPAGPhase1, dynFNPAGPhase1Sph, dynFNPAGPhase2, dynFNPAGPhase2Sph, dynFNPEGSph
import ODEs
import planetaryconstants as constants

import numpy as np
from scipy.integrate import solve_ivp
import sys
import matplotlib as mpl
import matplotlib.pyplot as plt
import time

plt.close('all')
mpl.rcParams["figure.autolayout"] = True
mpl.rcParams.update({'font.size': 15})

# =============================================================================
# Script
# =============================================================================
### CREATE params INPUT CLASS FOR NOMINAL VALUES
paramsNom = Params()
paramsNom.p = constants.MARS
# paramsNom.p.J2 = 0 # OVERRIDE - assume spherical planet for apples-to-apples comparison
# paramsNom.p.om = 0
paramsNom.returnTimeVectors = True

### INPUT ATM TABLE - GET ATM TABLE FROM GRAM DATA FILE
paramsNom.dMode = 'table'
filename = '../data/dens_Mars_nom.txt'
atmdata = np.genfromtxt(filename, names=True)
paramsNom.atmdat = np.array([atmdata['Var_X'], atmdata['DENSAV']])
    
### VEHICLE PARAMS
paramsNom.m = 2920 # kg, roughly MSL mass
paramsNom.CD = 1.6 # roughly MSL CD
# paramsNom.LD = 0.25
paramsNom.LD = 0
# paramsNom.BC = 130
paramsNom.BC = 35

paramsNom.A = paramsNom.m / (paramsNom.BC * paramsNom.CD)
paramsNom.CL = paramsNom.CD * paramsNom.LD
paramsNom.Rn = np.sqrt(paramsNom.A / np.pi) / 2

### WIND-RELATIVE INITIAL STATE
paramsNom.lat = 18.38
paramsNom.lon = -77.58
paramsNom.alt = paramsNom.p.halt
paramsNom.efpaWR = -12
paramsNom.hdaWR = 0
paramsNom.vmagWR = 6

### GET OTHER STATE TYPES (and assume t0 = 0)
paramsNom.x0, paramsNom.vInfvec_N = LLAEHV2RV(paramsNom.lat, paramsNom.lon,
                                        paramsNom.alt, paramsNom.efpaWR,
                                        paramsNom.hdaWR, paramsNom.vmagWR,
                                        paramsNom, 0)
paramsNom.v0 = Vinf2VN(paramsNom.x0, paramsNom.vInfvec_N, paramsNom, 0)
_, _, _, paramsNom.efpa, paramsNom.hda, paramsNom.vmag = \
    RV2LLAEHV(paramsNom.x0, paramsNom.v0, paramsNom, 0)

# xsphvec = np.array([paramsNom.alt + paramsNom.p.rad,
#                                       np.radians(paramsNom.lon), 
#                                       np.radians(paramsNom.lat),
#                                       paramsNom.vmagWR,
#                                       np.radians(paramsNom.efpaWR),
#                                       np.radians(paramsNom.hdaWR)])

# paramsNom.x0, paramsNom.v0 = sph2cart(xsphvec)

### NOMINAL SIMULATION PARAMS ###
t0 = 0
xx0vec1 = np.block([paramsNom.x0, paramsNom.v0])
# xx0vec1 = np.block([paramsNom.x0, paramsNom.vInfvec_N])
sig0 = 15
sigd = 150
ts = 160
ts1 = 137.9
ts2 = 140

ts = 0

paramsNom.rtol = 1e-10
paramsNom.atol = 1e-10
paramsNom.errtol1 = 0
paramsNom.errtol2 = 0
paramsNom.dtGdn = 1 # s
paramsNom.hmin = 10
paramsNom.hmax = paramsNom.p.halt + 1e-7
paramsNom.tf = 5000
paramsNom.raStar = 250 + paramsNom.p.rad
paramsNom.rpStar = 250 + paramsNom.p.rad

### SET TRUE MC SIMULATION PARAMS ###
paramsTrue = paramsNom


# =============================================================================
# Demonstrate dynamics
# =============================================================================
# cartesian:
tic1 = time.time()
xxvecs1, tvec1 = dynFNPAGPhase1(xx0vec1, 0, ts, sig0, sigd, paramsNom, returnTime = True)
toc1 = time.time()
# xxvecs1, tvec1 = dynFNPAGPhase2(xx0vec1, 0, sigd, paramsNom, returnTime = True)
raf, rpf = getApses(xxvecs1[:3,-1], xxvecs1[3:,-1], paramsNom)
print(raf - paramsNom.p.rad)

# now convert vmag to planet-relative for plotting
rvecs1 = xxvecs1[:3,:]
vvecs1 = xxvecs1[3:,:]
vInfvecs1 = []
for rvec, vvec, t in zip(rvecs1.T,vvecs1.T, tvec1):
    vInfvecs1.append(VN2Vinf(rvec, vvec, paramsNom, t))
    
vInfvecs1 = np.asarray(vInfvecs1)
vmagWR = np.linalg.norm(vInfvecs1, axis=1)



h = np.linalg.norm(xxvecs1[:3,:], axis = 0)
vmag = np.linalg.norm(xxvecs1[3:,:], axis = 0)

fig = plt.figure()
ax = fig.add_subplot(111)
ax.plot(vmagWR, h, label = 'cartesian EOMs')
ax.grid()

# spherical:
## FIRST convert spherical to cartesian using my method
# # x0, vInfvec_N = LLAEHV2RV(paramsNom.lat, paramsNom.lon,
# #                                         paramsNom.alt, paramsNom.efpaWR,
# #                                         paramsNom.hdaWR, paramsNom.vmagWR,
# #                                         paramsNom, 0)
# x0 = paramsNom.x0 * 1e3
# vInfvec_N = paramsNom.vInfvec_N * 1e3
# ## NEXT convert *back* to spherical, but now using Lu's definitions
# xsphvec2 = cart2sph(x0, vInfvec_N)
# xx0vec2 = xsphvec2


xx0vec2 = np.array([(paramsNom.alt + paramsNom.p.rad) * 1e3,
                    np.radians(paramsNom.lon), np.radians(paramsNom.lat),
                    paramsNom.vmagWR * 1e3, np.radians(paramsNom.efpaWR),
                    np.radians(paramsNom.hdaWR + 90)])
tic2 = time.time()
xxvecs2 = dynFNPAGPhase1Sph(xx0vec2, 0, ts, sig0, sigd, paramsNom)
toc2 = time.time()
# xxvecs2 = dynFNPAGPhase2Sph(xx0vec2, 0, sigd, paramsNom)
raf, rpf = getApsesSphPR(xxvecs2[:,-1], paramsNom)
print(raf/1e3 - paramsNom.p.rad)

h = np.linalg.norm(xxvecs2[:3,:], axis = 0) / 1e3
vmag = np.linalg.norm(xxvecs2[3:,:], axis = 0) / 1e3

ax.plot(vmag, h, label = 'spherical EOMs')
ax.legend()

# spherical aug
xx0vec3 = np.append(xx0vec2, 0)
tic3 = time.time()
e0 = 0
paramsNom.ef = 100
paramsNom.sigf = sigd
event1 = lambda t, y: ODEs.above_max_alt_sph(t, y, paramsNom)
event1.terminal = True
event1.direction = 1
event2 = lambda t, y: ODEs.below_min_alt_sph(t, y, paramsNom)
event2.terminal = True
tspan = (t0, paramsNom.tf)
sol = solve_ivp(lambda t, y: ODEs.sphericalEntryEOMsAug(t, y, sig0, e0,
                                                        paramsNom),
                tspan, xx0vec3, rtol = paramsNom.rtol, atol = paramsNom.atol,
                events = (event1, event2))
xxvecs3 = sol.y
toc3 = time.time()
raf, rpf = getApsesSphPR(xxvecs3[:,-1], paramsNom)
print(raf/1e3 - paramsNom.p.rad)

h = np.linalg.norm(xxvecs3[:3,:], axis = 0) / 1e3
vmag = np.linalg.norm(xxvecs3[3:,:], axis = 0) / 1e3

ax.plot(vmag, h, label = 'spherical EOMs Aug')
ax.legend()




print('cartesian EOMs took {0:.5f} s'.format(toc1-tic1))
print('spherical EOMs took {0:.3f} s'.format(toc2-tic2))
print('augmented spherical EOMs took {0:.3f} s'.format(toc3-tic3))

print('NOTE: J2 = {0:.3e}, planetary rotation rate Om = {1:.3e}'\
      .format(paramsNom.p.J2, paramsNom.p.om))


lon0 = xxvecs3[1,0]
lonf = xxvecs3[1,-1]
lat0 = xxvecs3[2,0]
latf = xxvecs3[2,-1]
dlon = abs(lonf - lon0)

dsig = np.arccos(np.sin(lat0) * np.sin(latf)\
                 + np.cos(lat0) * np.cos(latf) * np.cos(dlon))
drange = paramsNom.p.rad * dsig

print(-dsig)
print(xxvecs3[6,-1])

# NOTE: if there is any nonplanar motion, dsig will not match the computed
#   range xxvecs[6,-1]. dsig is the great circle arc, whereas the computed
#   range is the surface projection of the actual trajectory.



    
    
