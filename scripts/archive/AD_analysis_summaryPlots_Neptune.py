# -*- coding: utf-8 -*-
"""
Created on Tue Jul  7 17:05:31 2020

@author: Samuel Albert
"""

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.lines as mlines

from sim import Params, Outs
import constants

plt.close('all')

landcol = 'black'
spacecol = 'C4'
hafcol = 'C0'
gcol = 'C2'
qcol = 'C1'
QLcol = 'C3'

hafstyle = 'dashed'
gstyle = 'dotted'
qstyle = 'dashdot'
QLstyle = 'dotted'

# =============================================================================
# First file: full-lift-down
# =============================================================================

### Start with uncertainty bars on landline and spaceline
## Load file archive and get data
filename = './../results/sweeps/Neptune_27_20% low_0.25_180_0720052245.npz'
data = np.load(filename, allow_pickle=True)
params = data['params'][0] # array of 1
outsList = data['outsList']
efpaList = data['efpaList']
BCList = data['BCList']

## Create mesh grid for contour plots, reshape result arrays
BCgrid, EFPAgrid = np.meshgrid(BCList, efpaList)
fpafgrid = np.reshape([out.fpaf for out in outsList], BCgrid.shape)
engfgrid = np.reshape([out.engf for out in outsList], BCgrid.shape)

# find line between landing and aerocapture
landline_BC_low = BCList
landline_EFPA_low = []
for rind in range(fpafgrid.shape[1]):
    ind = next(ind for ind, val in enumerate(fpafgrid[:,rind])\
               if val < 0)
    landline_EFPA_low.append(efpaList[ind])
    
# find line between aerocapture and escape
spaceline_BC_low = []
spaceline_EFPA_low = []
for rind in range(engfgrid.shape[1]):
    ind = [ind for ind, val in enumerate(engfgrid[:,rind]) if val < 0]
    if len(ind) == 0:
        break
    else:
        spaceline_BC_low.append(BCList[rind])
        spaceline_EFPA_low.append(efpaList[ind[0]])
        
## Make countour plots
# fig = plt.figure(figsize = (6, 8.5))
# ax = fig.add_subplot(311)

fig, (ax1, ax2, ax3) = plt.subplots(3, 1, sharex = True, sharey = True)
fig.set_size_inches(6, 9)

ax1.grid()
ax1.set_title('full-lift-down', fontsize = 10)

# plot landline
ax1.plot(landline_EFPA_low, landline_BC_low, color = landcol,
        linewidth = 1, alpha = 0.3)
# plot spaceline
ax1.plot(spaceline_EFPA_low, spaceline_BC_low, color = spacecol,
        linewidth = 1, alpha = 0.3)

## Load file archive and get data
filename = './../results/sweeps/Neptune_27_20% high_0.25_180_0720071717.npz'
data = np.load(filename, allow_pickle=True)
params = data['params'][0] # array of 1
outsList = data['outsList']
efpaList = data['efpaList']
BCList = data['BCList']

## Create mesh grid for contour plots, reshape result arrays
BCgrid, EFPAgrid = np.meshgrid(BCList, efpaList)
fpafgrid = np.reshape([out.fpaf for out in outsList], BCgrid.shape)
engfgrid = np.reshape([out.engf for out in outsList], BCgrid.shape)

# find line between landing and aerocapture
landline_BC_high = BCList
landline_EFPA_high = []
for rind in range(fpafgrid.shape[1]):
    ind = next(ind for ind, val in enumerate(fpafgrid[:,rind])\
               if val < 0)
    landline_EFPA_high.append(efpaList[ind])
    
# find line between aerocapture and escape
spaceline_BC_high = []
spaceline_EFPA_high = []
for rind in range(engfgrid.shape[1]):
    ind = [ind for ind, val in enumerate(engfgrid[:,rind]) if val < 0]
    if len(ind) == 0:
        break
    else:
        spaceline_BC_high.append(BCList[rind])
        spaceline_EFPA_high.append(efpaList[ind[0]])
        
## Make countour plots

# plot landline
ax1.plot(landline_EFPA_high, landline_BC_high, color = landcol,
        linewidth = 1, alpha = 0.3)
# plot spaceline
ax1.plot(spaceline_EFPA_high, spaceline_BC_high, color = spacecol,
        linewidth = 1, alpha = 0.3)

# fill in uncertainty regions
ax1.fill_betweenx(landline_BC_low, landline_EFPA_low, landline_EFPA_high,
                facecolor = landcol, alpha = 0.3)
ax1.fill_betweenx(spaceline_BC_low, spaceline_EFPA_low, spaceline_EFPA_high,
                 facecolor = spacecol, alpha = 0.3)

### Now nominal atmosphere

## Load file archive and get data
filename = './../results/sweeps/Neptune_27_nom_0.25_180_0720062742.npz'
data = np.load(filename, allow_pickle=True)
params = data['params'][0] # array of 1
outsList = data['outsList']
efpaList = data['efpaList']
BCList = data['BCList']

## Create mesh grid for contour plots, reshape result arrays
BCgrid, EFPAgrid = np.meshgrid(BCList, efpaList)

ggrid = np.reshape([out.gpeak for out in outsList], BCgrid.shape)
hafgrid = np.reshape([out.haf for out in outsList], BCgrid.shape)
qgrid = np.reshape([out.qpeak for out in outsList], BCgrid.shape)
QLgrid = np.reshape([out.Qload for out in outsList], BCgrid.shape)
fpafgrid = np.reshape([out.fpaf for out in outsList], BCgrid.shape)
engfgrid = np.reshape([out.engf for out in outsList], BCgrid.shape)

##  Data pruning
hafgrid[hafgrid < 0] = np.inf # set haf to inf for hyperbolic cases
hafgrid[fpafgrid < 0] = np.nan # ignore haf values for landed cases

# find line between landing and aerocapture
landline_BC = BCList
landline_EFPA = []
for rind in range(fpafgrid.shape[1]):
    ind = next(ind for ind, val in enumerate(fpafgrid[:,rind])\
               if val < 0)
    landline_EFPA.append(efpaList[ind])
    
# find line between aerocapture and escape
spaceline_BC = []
spaceline_EFPA = []
for rind in range(engfgrid.shape[1]):
    ind = [ind for ind, val in enumerate(engfgrid[:,rind]) if val < 0]
    if len(ind) == 0:
        break
    else:
        spaceline_BC.append(BCList[rind])
        spaceline_EFPA.append(efpaList[ind[0]])



# plot landline
ax1.plot(landline_EFPA, landline_BC, color = landcol, linewidth = 3,
        label = 'cutoff b/w orbiters & probes')
# plot spaceline
ax1.plot(spaceline_EFPA, spaceline_BC, color = spacecol, linewidth = 3,
        label = 'cutoff b/w aerocapture & escape')

# contours:
haflevels = []
cp = ax1.contour(EFPAgrid.T, BCgrid.T, hafgrid.T, haflevels,
                colors = hafcol, linestyles = hafstyle)
ax1.clabel(cp, inline = True, colors = hafcol, fontsize = 9, fmt = '%1.0f')

glevels = [2, 35, 45, 55, 65, 75, 85]
cp = ax1.contour(EFPAgrid.T, BCgrid.T, ggrid.T, glevels,
                colors = gcol, linestyles = gstyle)
ax1.clabel(cp, inline = True, colors = gcol, fontsize = 9, fmt = '%1.0f')

qlevels = [200, 600, 1000]
cp = ax1.contour(EFPAgrid.T, BCgrid.T, qgrid.T, qlevels,
                colors = qcol, linestyles = qstyle)
ax1.clabel(cp, inline = True, colors = qcol, fontsize = 9, fmt = '%1.0f')

QLlevels = [2e4, 6e4]
cp = ax1.contour(EFPAgrid.T, BCgrid.T, QLgrid.T, QLlevels,
                colors = QLcol, linestyles = QLstyle)
ax1.clabel(cp, inline = True, colors = QLcol, fontsize = 9, fmt = '%1.0f')

# ## shrink plot to make room for legend
# box = ax1.get_position()
# ax1.set_position([box.x0, box.y0 ,
#                  box.width, box.height * 0.9])

land_proxy = mlines.Line2D([], [], color = landcol, linewidth = 3,
                           label = r'orbiters/probes cutoff, '
                                r'$\mathregular{\pm20\%\ \rho}$')
space_proxy = mlines.Line2D([], [], color = spacecol, linewidth = 3,
                            label = r'orbiters/escape cutoff, '
                                r'$\mathregular{\pm20\%\ \rho}$')
haf_proxy = mlines.Line2D([], [], color = hafcol, linestyle = hafstyle,
                          label = 'apoapsis altitude, km')
g_proxy = mlines.Line2D([], [], color = gcol, linestyle = gstyle,
                        label = 'peak g load, Earth g\'s')
q_proxy = mlines.Line2D([], [], color = qcol, linestyle = qstyle,
                        label = 'peak heat rate, '
                            '$\mathregular{W/cm^2}$')
QL_proxy = mlines.Line2D([], [], color = QLcol, linestyle = QLstyle,
                         label = 'heat load, '
                             '$\mathregular{J/cm^2}$')

fig.legend(handles=[land_proxy, space_proxy, haf_proxy,
                    g_proxy, q_proxy, QL_proxy],
           loc='lower center', bbox_to_anchor = (0.5, 0.91),
           ncol=2, borderaxespad=0., prop = {'size': 9})


# =============================================================================
# Second file: ballistic
# =============================================================================

### Start with uncertainty bars on landline and spaceline
## Load file archive and get data
filename = './../results/sweeps/Neptune_27_20% low_0_0_0720052600.npz'
data = np.load(filename, allow_pickle=True)
params = data['params'][0] # array of 1
outsList = data['outsList']
efpaList = data['efpaList']
BCList = data['BCList']

## Create mesh grid for contour plots, reshape result arrays
BCgrid, EFPAgrid = np.meshgrid(BCList, efpaList)
fpafgrid = np.reshape([out.fpaf for out in outsList], BCgrid.shape)
engfgrid = np.reshape([out.engf for out in outsList], BCgrid.shape)

# find line between landing and aerocapture
landline_BC_low = BCList
landline_EFPA_low = []
for rind in range(fpafgrid.shape[1]):
    ind = next(ind for ind, val in enumerate(fpafgrid[:,rind])\
               if val < 0)
    landline_EFPA_low.append(efpaList[ind])
    
# find line between aerocapture and escape
spaceline_BC_low = []
spaceline_EFPA_low = []
for rind in range(engfgrid.shape[1]):
    ind = [ind for ind, val in enumerate(engfgrid[:,rind]) if val < 0]
    if len(ind) == 0:
        break
    else:
        spaceline_BC_low.append(BCList[rind])
        spaceline_EFPA_low.append(efpaList[ind[0]])
        
## Make countour plots
ax2.set_title('ballistic', fontsize = 10)
ax2.grid()

# plot landline
ax2.plot(landline_EFPA_low, landline_BC_low, color = landcol,
        linewidth = 1, alpha = 0.3)
# plot spaceline
ax2.plot(spaceline_EFPA_low, spaceline_BC_low, color = spacecol,
        linewidth = 1, alpha = 0.3)

## Load file archive and get data
filename = './../results/sweeps/Neptune_27_20% high_0_0_0720071518.npz'
data = np.load(filename, allow_pickle=True)
params = data['params'][0] # array of 1
outsList = data['outsList']
efpaList = data['efpaList']
BCList = data['BCList']

## Create mesh grid for contour plots, reshape result arrays
BCgrid, EFPAgrid = np.meshgrid(BCList, efpaList)
fpafgrid = np.reshape([out.fpaf for out in outsList], BCgrid.shape)
engfgrid = np.reshape([out.engf for out in outsList], BCgrid.shape)

# find line between landing and aerocapture
landline_BC_high = BCList
landline_EFPA_high = []
for rind in range(fpafgrid.shape[1]):
    ind = next(ind for ind, val in enumerate(fpafgrid[:,rind])\
               if val < 0)
    landline_EFPA_high.append(efpaList[ind])
    
# find line between aerocapture and escape
spaceline_BC_high = []
spaceline_EFPA_high = []
for rind in range(engfgrid.shape[1]):
    ind = [ind for ind, val in enumerate(engfgrid[:,rind]) if val < 0]
    if len(ind) == 0:
        break
    else:
        spaceline_BC_high.append(BCList[rind])
        spaceline_EFPA_high.append(efpaList[ind[0]])
        
## Make countour plots

# plot landline
ax2.plot(landline_EFPA_high, landline_BC_high, color = landcol,
        linewidth = 1, alpha = 0.3)
# plot spaceline
ax2.plot(spaceline_EFPA_high, spaceline_BC_high, color = spacecol,
        linewidth = 1, alpha = 0.3)

# fill in uncertainty regions
ax2.fill_betweenx(landline_BC_low, landline_EFPA_low, landline_EFPA_high,
                facecolor = landcol, alpha = 0.3)
ax2.fill_betweenx(spaceline_BC_low, spaceline_EFPA_low, spaceline_EFPA_high,
                 facecolor = spacecol, alpha = 0.3)

## Load file archive and get data
filename = './../results/sweeps/Neptune_27_nom_0_0_0720062828.npz'
data = np.load(filename, allow_pickle=True)
params = data['params'][0] # array of 1
outsList = data['outsList']
efpaList = data['efpaList']
BCList = data['BCList']

## Create mesh grid for contour plots, reshape result arrays
BCgrid, EFPAgrid = np.meshgrid(BCList, efpaList)

ggrid = np.reshape([out.gpeak for out in outsList], BCgrid.shape)
hafgrid = np.reshape([out.haf for out in outsList], BCgrid.shape)
qgrid = np.reshape([out.qpeak for out in outsList], BCgrid.shape)
QLgrid = np.reshape([out.Qload for out in outsList], BCgrid.shape)
fpafgrid = np.reshape([out.fpaf for out in outsList], BCgrid.shape)
engfgrid = np.reshape([out.engf for out in outsList], BCgrid.shape)

##  Data pruning
hafgrid[hafgrid < 0] = np.inf # set haf to inf for hyperbolic cases
hafgrid[fpafgrid < 0] = np.nan # ignore haf values for landed cases

# find line between landing and aerocapture
landline_BC = BCList
landline_EFPA = []
for rind in range(fpafgrid.shape[1]):
    ind = next(ind for ind, val in enumerate(fpafgrid[:,rind])\
               if val < 0)
    landline_EFPA.append(efpaList[ind])
    
# find line between aerocapture and escape
spaceline_BC = []
spaceline_EFPA = []
for rind in range(engfgrid.shape[1]):
    ind = [ind for ind, val in enumerate(engfgrid[:,rind]) if val < 0]
    if len(ind) == 0:
        break
    else:
        spaceline_BC.append(BCList[rind])
        spaceline_EFPA.append(efpaList[ind[0]])

## Make countour plots

# plot landline
ax2.plot(landline_EFPA, landline_BC, color = landcol, linewidth = 3)
# plot spaceline
ax2.plot(spaceline_EFPA, spaceline_BC, color = spacecol, linewidth = 3)

## Contours:

haflevels = [5e4]
cp = ax2.contour(EFPAgrid.T, BCgrid.T, hafgrid.T, haflevels,
                colors = hafcol, linestyles = hafstyle)
ax2.clabel(cp, inline = True, colors = hafcol, fontsize = 9, fmt = '%1.0f')

glevels = [2, 4, 15, 25, 35, 45]
cp = ax2.contour(EFPAgrid.T, BCgrid.T, ggrid.T, glevels,
                colors = gcol, linestyles = gstyle)
ax2.clabel(cp, inline = True, colors = gcol, fontsize = 9, fmt = '%1.0f')

qlevels = [200, 600, 1000]
cp = ax2.contour(EFPAgrid.T, BCgrid.T, qgrid.T, qlevels,
                colors = qcol, linestyles = qstyle)
ax2.clabel(cp, inline = True, colors = qcol, fontsize = 9, fmt = '%1.0f')

QLlevels = [2e4, 6e4, 1e5]
cp = ax2.contour(EFPAgrid.T, BCgrid.T, QLgrid.T, QLlevels,
                colors = QLcol, linestyles = QLstyle)
ax2.clabel(cp, inline = True, colors = QLcol, fontsize = 9, fmt = '%1.0f')


# =============================================================================
# Third file: full-lift-up
# =============================================================================

### Start with uncertainty bars on landline and spaceline
## Load file archive and get data
filename = './../results/sweeps/Neptune_27_20% low_0.25_0_0720080534.npz'
data = np.load(filename, allow_pickle=True)
params = data['params'][0] # array of 1
outsList = data['outsList']
efpaList = data['efpaList']
BCList = data['BCList']

## Create mesh grid for contour plots, reshape result arrays
BCgrid, EFPAgrid = np.meshgrid(BCList, efpaList)
fpafgrid = np.reshape([out.fpaf for out in outsList], BCgrid.shape)
engfgrid = np.reshape([out.engf for out in outsList], BCgrid.shape)

# find line between landing and aerocapture
landline_BC_low = BCList
landline_EFPA_low = []
for rind in range(fpafgrid.shape[1]):
    ind = next(ind for ind, val in enumerate(fpafgrid[:,rind])\
               if val < 0)
    landline_EFPA_low.append(efpaList[ind])
    
# find line between aerocapture and escape
spaceline_BC_low = []
spaceline_EFPA_low = []
for rind in range(engfgrid.shape[1]):
    ind = [ind for ind, val in enumerate(engfgrid[:,rind]) if val < 0]
    if len(ind) == 0:
        break
    else:
        spaceline_BC_low.append(BCList[rind])
        spaceline_EFPA_low.append(efpaList[ind[0]])
        
## Make countour plots
ax3.grid()
ax3.set_title('full-lift-up', fontsize = 10)

# plot landline
ax3.plot(landline_EFPA_low, landline_BC_low, color = landcol,
        linewidth = 1, alpha = 0.3)
# plot spaceline
ax3.plot(spaceline_EFPA_low, spaceline_BC_low, color = spacecol,
        linewidth = 1, alpha = 0.3)

## Load file archive and get data
filename = './../results/sweeps/Neptune_27_20% high_0.25_0_0720094323.npz'
data = np.load(filename, allow_pickle=True)
params = data['params'][0] # array of 1
outsList = data['outsList']
efpaList = data['efpaList']
BCList = data['BCList']

## Create mesh grid for contour plots, reshape result arrays
BCgrid, EFPAgrid = np.meshgrid(BCList, efpaList)
fpafgrid = np.reshape([out.fpaf for out in outsList], BCgrid.shape)
engfgrid = np.reshape([out.engf for out in outsList], BCgrid.shape)

# find line between landing and aerocapture
landline_BC_high = BCList
landline_EFPA_high = []
for rind in range(fpafgrid.shape[1]):
    ind = next(ind for ind, val in enumerate(fpafgrid[:,rind])\
               if val < 0)
    landline_EFPA_high.append(efpaList[ind])
    
# find line between aerocapture and escape
spaceline_BC_high = []
spaceline_EFPA_high = []
for rind in range(engfgrid.shape[1]):
    ind = [ind for ind, val in enumerate(engfgrid[:,rind]) if val < 0]
    if len(ind) == 0:
        break
    else:
        spaceline_BC_high.append(BCList[rind])
        spaceline_EFPA_high.append(efpaList[ind[0]])
        
## Make countour plots

# plot landline
ax3.plot(landline_EFPA_high, landline_BC_high, color = landcol,
        linewidth = 1, alpha = 0.3)
# plot spaceline
ax3.plot(spaceline_EFPA_high, spaceline_BC_high, color = spacecol,
        linewidth = 1, alpha = 0.3)

# fill in uncertainty regions
ax3.fill_betweenx(landline_BC_low, landline_EFPA_low, landline_EFPA_high,
                facecolor = landcol, alpha = 0.3)
ax3.fill_betweenx(spaceline_BC_low, spaceline_EFPA_low, spaceline_EFPA_high,
                 facecolor = spacecol, alpha = 0.3)


## Load file archive and get data
filename = './../results/sweeps/Neptune_27_nom_0.25_0_0720085952.npz'
data = np.load(filename, allow_pickle=True)
params = data['params'][0] # array of 1
outsList = data['outsList']
efpaList = data['efpaList']
BCList = data['BCList']

## Create mesh grid for contour plots, reshape result arrays
BCgrid, EFPAgrid = np.meshgrid(BCList, efpaList)

ggrid = np.reshape([out.gpeak for out in outsList], BCgrid.shape)
hafgrid = np.reshape([out.haf for out in outsList], BCgrid.shape)
qgrid = np.reshape([out.qpeak for out in outsList], BCgrid.shape)
QLgrid = np.reshape([out.Qload for out in outsList], BCgrid.shape)
fpafgrid = np.reshape([out.fpaf for out in outsList], BCgrid.shape)
engfgrid = np.reshape([out.engf for out in outsList], BCgrid.shape)

##  Data pruning
hafgrid[hafgrid < 0] = np.inf # set haf to inf for hyperbolic cases
hafgrid[fpafgrid < 0] = np.nan # ignore haf values for landed cases

# find line between landing and aerocapture
landline_BC = BCList
landline_EFPA = []
for rind in range(fpafgrid.shape[1]):
    ind = next(ind for ind, val in enumerate(fpafgrid[:,rind])\
               if val < 0)
    landline_EFPA.append(efpaList[ind])
    
# find line between aerocapture and escape
spaceline_BC = []
spaceline_EFPA = []
for rind in range(engfgrid.shape[1]):
    ind = [ind for ind, val in enumerate(engfgrid[:,rind]) if val < 0]
    if len(ind) == 0:
        break
    else:
        spaceline_BC.append(BCList[rind])
        spaceline_EFPA.append(efpaList[ind[0]])

## Make countour plots

# plot landline
ax3.plot(landline_EFPA, landline_BC, color = landcol, linewidth = 3)
# plot spaceline
ax3.plot(spaceline_EFPA, spaceline_BC, color = spacecol, linewidth = 3)

## contours:
    
haflevels = [2e3, 1e4, 2e4, 5e4, 15e4]
cp = ax3.contour(EFPAgrid.T, BCgrid.T, hafgrid.T, haflevels,
                colors = hafcol, linestyles = hafstyle)
ax3.clabel(cp, inline = True, colors = hafcol, fontsize = 9, fmt = '%1.0f')

glevels = [2, 4, 10, 15, 25]
cp = ax3.contour(EFPAgrid.T, BCgrid.T, ggrid.T, glevels,
                colors = gcol, linestyles = gstyle)
ax3.clabel(cp, inline = True, colors = gcol, fontsize = 9, fmt = '%1.0f')

qlevels = [200, 600, 1000]
cp = ax3.contour(EFPAgrid.T, BCgrid.T, qgrid.T, qlevels,
                colors = qcol, linestyles = qstyle)
ax3.clabel(cp, inline = True, colors = qcol, fontsize = 9, fmt = '%1.0f')

QLlevels = [2e4, 6e4, 1e5]
cp = ax3.contour(EFPAgrid.T, BCgrid.T, QLgrid.T, QLlevels,
                colors = QLcol, linestyles = QLstyle)
ax3.clabel(cp, inline = True, colors = QLcol, fontsize = 9, fmt = '%1.0f')


# =============================================================================
# Legend, Plot cleanup
# =============================================================================
plt.xlabel('Entry Flight Path Angle, deg', fontsize = 9)
plt.ylabel('Ballistic Coefficient, $\mathregular{kg/m^2}$', fontsize = 9)
# plt.tight_layout()


plt.subplots_adjust(left = 0.11,
                bottom = 0.05,
                right = 0.9,
                top = 0.88,
                wspace = 0.2,
                hspace = 0.17)





















