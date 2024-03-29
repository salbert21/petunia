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
filename = './../results/sweeps/Earth_11_20% low_0.25_180_0709192205.npz'
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

fig, (ax2, ax3) = plt.subplots(2, 1, sharex = True, sharey = True)
fig.set_size_inches(12, 12)

ax1.grid()
ax1.set_title('full-lift-down', fontsize = 10)

# plot landline
ax1.plot(landline_EFPA_low, landline_BC_low, color = landcol,
        linewidth = 1, alpha = 0.3)
# plot spaceline
ax1.plot(spaceline_EFPA_low, spaceline_BC_low, color = spacecol,
        linewidth = 1, alpha = 0.3)

## Load file archive and get data
filename = './../results/sweeps/Earth_11_20% high_0.25_180_0709184117.npz'
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
filename = './../results/sweeps/Earth_11_nom_0.25_180_0709182147.npz'
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
haflevels = [25e3, 1e5, 3e5]
cp = ax1.contour(EFPAgrid.T, BCgrid.T, hafgrid.T, haflevels,
                colors = hafcol, linestyles = hafstyle)
ax1.clabel(cp, inline = True, colors = hafcol, fontsize = 9, fmt = '%1.0f')

glevels = [25, 35, 40, 45]
cp = ax1.contour(EFPAgrid.T, BCgrid.T, ggrid.T, glevels,
                colors = gcol, linestyles = gstyle)
ax1.clabel(cp, inline = True, colors = gcol, fontsize = 9, fmt = '%1.0f')

qlevels = [50, 150, 250, 350]
cp = ax1.contour(EFPAgrid.T, BCgrid.T, qgrid.T, qlevels,
                colors = qcol, linestyles = qstyle)
ax1.clabel(cp, inline = True, colors = qcol, fontsize = 9, fmt = '%1.0f')

QLlevels = [5e3, 10e3, 15e3]
cp = ax1.contour(EFPAgrid.T, BCgrid.T, QLgrid.T, QLlevels,
                colors = QLcol, linestyles = QLstyle)
ax1.clabel(cp, inline = True, colors = QLcol, fontsize = 9, fmt = '%1.0f')

# ## shrink plot to make room for legend
# box = ax1.get_position()
# ax1.set_position([box.x0, box.y0 ,
#                  box.width, box.height * 0.9])




# =============================================================================
# Second file: ballistic
# =============================================================================

### Start with uncertainty bars on landline and spaceline
## Load file archive and get data
filename = './../results/sweeps/Earth_11_20% low_0_0_0709180245.npz'
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
ax2.set_title('ballistic', fontsize = 20)
ax2.grid()

# plot landline
ax2.plot(landline_EFPA_low, landline_BC_low, color = landcol,
        linewidth = 1, alpha = 0.3)
# plot spaceline
ax2.plot(spaceline_EFPA_low, spaceline_BC_low, color = spacecol,
        linewidth = 1, alpha = 0.3)

## Load file archive and get data
filename = './../results/sweeps/Earth_11_20% high_0_0_0709183314.npz'
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
filename = './../results/sweeps/Earth_11_nom_0_0_0709181558.npz'
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

# ## Contours:
    
haflevels = [5e3, 25e3, 1e5, 3e5]
cp = ax2.contour(EFPAgrid.T, BCgrid.T, hafgrid.T, haflevels,
                colors = hafcol, linestyles = hafstyle)
ax2.clabel(cp, inline = True, colors = hafcol, fontsize = 15, fmt = '%1.0f')

glevels = [1, 3, 5, 10, 15, 20, 25]
cp = ax2.contour(EFPAgrid.T, BCgrid.T, ggrid.T, glevels,
                colors = gcol, linestyles = gstyle)
ax2.clabel(cp, inline = True, colors = gcol, fontsize = 15, fmt = '%1.0f')

qlevels = [50, 150, 250, 350]
cp = ax2.contour(EFPAgrid.T, BCgrid.T, qgrid.T, qlevels,
                colors = qcol, linestyles = qstyle)
ax2.clabel(cp, inline = True, colors = qcol, fontsize = 15, fmt = '%1.0f')

QLlevels = [5e3, 10e3, 15e3, 20e3]
cp = ax2.contour(EFPAgrid.T, BCgrid.T, QLgrid.T, QLlevels,
                colors = QLcol, linestyles = QLstyle)
ax2.clabel(cp, inline = True, colors = QLcol, fontsize = 15, fmt = '%1.0f')

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

ax2.set_xlim((-7.3800000000000034, -3.7))
ax2.set_ylim((10.0, 197.5))

leg = plt.legend(handles=[land_proxy, space_proxy,
                            haf_proxy, g_proxy, q_proxy, QL_proxy],
            loc='lower center', bbox_to_anchor = (0.5, 0.8),
            ncol=2, borderaxespad=0., prop = {'size': 15})
leg.set_draggable(True)



# =============================================================================
# Third file: full-lift-up
# =============================================================================

### Start with uncertainty bars on landline and spaceline
## Load file archive and get data
filename = './../results/sweeps/Earth_11_20% low_0.25_0_0709162350.npz'
data = np.load(filename, allow_pickle=True)
params = data['params'][0] # array of 1
outsList = data['outsList']
efpaList = data['efpaList']
BCList = data['BCList']

## Create mesh grid for contour plots, reshape result arrays
BCgrid2, EFPAgrid = np.meshgrid(BCList, efpaList)
fpafgrid2 = np.reshape([out.fpaf for out in outsList], BCgrid2.shape)
engfgrid2 = np.reshape([out.engf for out in outsList], BCgrid2.shape)

# find line between landing and aerocapture
landline_BC_low = BCList
landline_EFPA_low = []
for rind in range(fpafgrid2.shape[1]):
    ind = next(ind for ind, val in enumerate(fpafgrid2[:,rind])\
                if val < 0)
    landline_EFPA_low.append(efpaList[ind])
    
# find line between aerocapture and escape
spaceline_BC_low = []
spaceline_EFPA_low = []
for rind in range(engfgrid2.shape[1]):
    ind = [ind for ind, val in enumerate(engfgrid2[:,rind]) if val < 0]
    if len(ind) == 0:
        break
    else:
        spaceline_BC_low.append(BCList[rind])
        spaceline_EFPA_low.append(efpaList[ind[0]])
        
## Make countour plots
ax3.grid()
ax3.set_title('full-lift-up', fontsize = 20)

# plot landline
ax3.plot(landline_EFPA_low, landline_BC_low, color = landcol,
        linewidth = 1, alpha = 0.3)
# plot spaceline
ax3.plot(spaceline_EFPA_low, spaceline_BC_low, color = spacecol,
        linewidth = 1, alpha = 0.3)

## Load file archive and get data
filename = './../results/sweeps/Earth_11_20% high_0.25_0_0709165632.npz'
data = np.load(filename, allow_pickle=True)
params = data['params'][0] # array of 1
outsList = data['outsList']
efpaList = data['efpaList']
BCList = data['BCList']

## Create mesh grid for contour plots, reshape result arrays
BCgrid2, EFPAgrid = np.meshgrid(BCList, efpaList)
fpafgrid2 = np.reshape([out.fpaf for out in outsList], BCgrid2.shape)
engfgrid2 = np.reshape([out.engf for out in outsList], BCgrid2.shape)

# find line between landing and aerocapture
landline_BC_high = BCList
landline_EFPA_high = []
for rind in range(fpafgrid2.shape[1]):
    ind = next(ind for ind, val in enumerate(fpafgrid2[:,rind])\
                if val < 0)
    landline_EFPA_high.append(efpaList[ind])
    
# find line between aerocapture and escape
spaceline_BC_high = []
spaceline_EFPA_high = []
for rind in range(engfgrid2.shape[1]):
    ind = [ind for ind, val in enumerate(engfgrid2[:,rind]) if val < 0]
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
filename = './../results/sweeps/Earth_11_nom_0.25_0_0709164051.npz'
data = np.load(filename, allow_pickle=True)
params = data['params'][0] # array of 1
outsList = data['outsList']
efpaList = data['efpaList']
BCList = data['BCList']

## Create mesh grid for contour plots, reshape result arrays
BCgrid2, EFPAgrid = np.meshgrid(BCList, efpaList)

ggrid2 = np.reshape([out.gpeak for out in outsList], BCgrid2.shape)
hafgrid2 = np.reshape([out.haf for out in outsList], BCgrid2.shape)
qgrid2 = np.reshape([out.qpeak for out in outsList], BCgrid2.shape)
QLgrid2 = np.reshape([out.Qload for out in outsList], BCgrid2.shape)
fpafgrid2 = np.reshape([out.fpaf for out in outsList], BCgrid2.shape)
engfgrid2 = np.reshape([out.engf for out in outsList], BCgrid2.shape)

##  Data pruning
hafgrid2[hafgrid2 < 0] = np.inf # set haf to inf for hyperbolic cases
hafgrid2[fpafgrid2 < 0] = np.nan # ignore haf values for landed cases

# find line between landing and aerocapture
landline_BC = BCList
landline_EFPA = []
for rind in range(fpafgrid2.shape[1]):
    ind = next(ind for ind, val in enumerate(fpafgrid2[:,rind])\
                if val < 0)
    landline_EFPA.append(efpaList[ind])
    
# find line between aerocapture and escape
spaceline_BC = []
spaceline_EFPA = []
for rind in range(engfgrid2.shape[1]):
    ind = [ind for ind, val in enumerate(engfgrid2[:,rind]) if val < 0]
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

# contours:
haflevels = [150, 500, 3e3, 1e4, 3e4, 3e5]
cp = ax3.contour(EFPAgrid.T, BCgrid2.T, hafgrid2.T, haflevels,
                colors = hafcol, linestyles = hafstyle)
ax3.clabel(cp, inline = True, colors = hafcol, fontsize = 15, fmt = '%1.0f')

glevels = [1, 3, 5, 8, 12, 15]
cp = ax3.contour(EFPAgrid.T, BCgrid2.T, ggrid2.T, glevels,
                colors = gcol, linestyles = gstyle)
ax3.clabel(cp, inline = True, colors = gcol, fontsize = 15, fmt = '%1.0f')

qlevels = [50, 150, 250]
cp = ax3.contour(EFPAgrid.T, BCgrid2.T, qgrid2.T, qlevels,
                colors = qcol, linestyles = qstyle)
ax3.clabel(cp, inline = True, colors = qcol, fontsize = 15, fmt = '%1.0f')

QLlevels = [5e3, 10e3, 15e3]
cp = ax3.contour(EFPAgrid.T, BCgrid2.T, QLgrid2.T, QLlevels,
                colors = QLcol, linestyles = QLstyle)
ax3.clabel(cp, inline = True, colors = QLcol, fontsize = 15, fmt = '%1.0f')


# =============================================================================
# Legend, Plot cleanup
# =============================================================================
plt.xlabel('Entry Flight Path Angle, deg', fontsize = 15)
plt.ylabel('Ballistic Coefficient, $\mathregular{kg/m^2}$', fontsize = 15)
# fig.legend(handles=[land_proxy, space_proxy, haf_proxy,
#                      g_proxy, q_proxy, QL_proxy], loc='upper right',
#             ncol=2, prop = {'size': 15})
# plt.tight_layout()


# plt.subplots_adjust(left = 0.11,
#                 bottom = 0.05,
#                 right = 0.9,
#                 top = 0.88,
#                 wspace = 0.2,
#                 hspace = 0.17)

plt.subplots_adjust(top=0.87,
                    bottom=0.045,
                    left=0.045,
                    right=0.985,
                    hspace=0.2,
                    wspace=0.2)





















