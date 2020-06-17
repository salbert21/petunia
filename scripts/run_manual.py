# -*- coding: utf-8 -*-
"""
Created on Thu Apr  2 14:36:28 2020

@author: Samuel
"""

import matplotlib.pyplot as plt
from mpl_toolkits import mplot3d

import numpy as np
import ODEs
from scipy.integrate import solve_ivp
from atm import getMCdens
import constants

class params:
    def __init__(self, dMode):
        self.dMode = dMode
    
plt.close('all')
        

dMode = 'table'
params1 = params(dMode)
params1.p = constants.EARTH

params1.m = 2000

params1.CL = 0.29
params1.CD = 1.2121
params1.A = 15
params1.bank = np.radians(0)

### GET ATM TABLE FROM BINARY EARTHGRAM DATA FILE
filename = '../data/rawOutput.txt'
# get Nmc atmosphere profiles
Nmc = 1
i_trial = 0
densPert, densMean, h = getMCdens(filename, Nmc)
# at some point would be good to build this as a pandas df instead of np array
rhoTable = np.array([h,densPert[:,i_trial]])
params.atmdat = rhoTable

params1.dMode = 'table'

params1.hmin = 30
params1.hmax = 125

tspan = np.linspace(0,300,5000) # integrate, 5,000 time steps


# r0vec_N = np.array([-6402,-1809,1065]) * 1e3
# v0vec_N = np.array([0.999,-6.471,-4.302]) * 1e3

r0vec_N = np.array([6478.100, 0, 0])
# v0vec_N = np.array([-1251.45758126121, 472.3899576546, 7901.50959768473])
v0vec_N = np.array([-1.25145758, 0.47238996, 7.9015096])


y0 = np.block([r0vec_N, v0vec_N])


event1 = lambda t, y: ODEs.below_min_alt(t, y, params1)
event1.terminal = True

event2 = lambda t, y: ODEs.above_max_alt(t, y, params1)
event2.terminal = True


sol = solve_ivp(lambda t, y: ODEs.dynamics(t,y,params1), 
                [tspan[0], tspan[-1]], y0.flatten(),
                rtol=1e-9,atol=1e-9, t_eval=tspan,
                events=(event1,event2)) 
                # normally set tols to 1e-12
print(sol.message)

rvec_N = sol.y[0:3,:] #/ 1e3 # convert to km
vvec_N = sol.y[3:6,:] #/ 1e3 # convert to km/s

alt = np.linalg.norm(rvec_N, axis=0) - params1.p.rad #/1e3
vmag = np.linalg.norm(vvec_N, axis=0)

fig = plt.figure(2)
ax = fig.add_subplot(111)
ax.plot(vmag, alt)
ax.set_xlabel('Inertial velocity, km/s')
ax.set_ylabel('Altitude, km')
ax.grid()

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
