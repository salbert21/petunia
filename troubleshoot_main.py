# -*- coding: utf-8 -*-
"""
Created on Fri Apr 10 14:24:42 2020

@author: Samuel
"""


import matplotlib.pyplot as plt
from mpl_toolkits import mplot3d

import numpy as np
import ODEs
from scipy.integrate import solve_ivp
from atm import getAtmdat

class params:
    def __init__(self, mu, m):
        self.mu = mu
        self.m = m
        
    class p:
        pass
    
plt.close('all')
        

omE = 2 * np.pi / (0.99726968 * 86400) # rad/s rotation rate of Earth
radE = 6378.1363 * 1e3
mu = 3.986e5 * 1e9 # m^3/s^2

## TROUBLESHOOTING
radE = 6378100
omE = 7.2921066e-05
mu = 398600000000000
##

m = 3380

params1 = params(mu, m) # mu in m^3/s^2 !
params1.p.om = omE
params1.p.rad = radE

# params1.CL = 0.363905772129617
# params1.CD = 1.45562308851847
# params1.A = 15.6893896627335
# params1.bank = np.radians(0)

params1.CL = 0.363905772129617
params1.CD = 1.45562308851847
params1.A = 23.2202967008456
params1.bank = np.radians(0)

params1.atmdat = getAtmdat('data/atm_earth_gram2016.csv')

params1.hmin = 30 * 1e3
params1.hmax = 125 * 1e3

tspan = np.linspace(0,300,5000) # integrate for 1 day, 5,000 time steps


# r0vec_N = np.array([-6402,-1809,1065]) * 1e3
# v0vec_N = np.array([0.999,-6.471,-4.302]) * 1e3

r0vec_N = np.array([6478100, 0, 0])
v0vec_N = np.array([-1251.45758126121, 472.3899576546, 7901.50959768473])


y0 = np.block([r0vec_N, v0vec_N])

dydt = ODEs.dynamics(0,y0,params1)


# event1 = lambda t, y: ODEs.below_min_alt(t, y, params1)
# event1.terminal = True

# event2 = lambda t, y: ODEs.above_max_alt(t, y, params1)
# event2.terminal = True


# sol = solve_ivp(lambda t, y: ODEs.dynamics(t,y,params1), 
#                 [tspan[0], tspan[-1]], y0.flatten(),
#                 rtol=1e-9,atol=1e-9, t_eval=tspan,
#                 events=(event1,event2)) 
#                 # normally set tols to 1e-12
# print(sol.message)

# rvec_N = sol.y[0:3,:] / 1e3 # convert to km
# vvec_N = sol.y[3:6,:] / 1e3 # convert to km/s

# alt = np.linalg.norm(rvec_N, axis=0) - radE/1e3
# vmag = np.linalg.norm(vvec_N, axis=0)

# fig = plt.figure(2)
# ax = fig.add_subplot(111)
# ax.plot(vmag, alt)

# fig = plt.figure(11)
# ax = fig.add_subplot(111)
# ax.plot(rvec_N[0,:], rvec_N[1,:])
# ax.set_xlabel('x inertial position (km)')
# ax.set_ylabel('y inertial position (km)')

# fig = plt.figure(12)
# ax = fig.add_subplot(111)
# ax.plot(rvec_N[1,:], rvec_N[2,:])
# ax.set_xlabel('y inertial position (km)')
# ax.set_ylabel('z inertial position (km)')

# fig = plt.figure(13)
# ax = fig.add_subplot(111)
# ax.plot(rvec_N[0,:], rvec_N[2,:])
# ax.set_xlabel('x inertial position (km)')
# ax.set_ylabel('z inertial position (km)')



# fig = plt.figure(1)
# ax = plt.axes(projection="3d")

# ax.plot3D(rvec_N[0,:],rvec_N[1,:],rvec_N[2,:], label='petunia solution')
