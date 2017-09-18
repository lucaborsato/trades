#!/usr/bin/env python
# -*- coding: utf-8 -*-

from __future__ import division # no more "zero" integer division bugs!:P
# import argparse
import numpy as np # array
# import h5py
# import astropy.time as atime
import os
import sys
#import emcee.interruptible_pool as emcee_pool
# import gc

import matplotlib as mpl
#from matplotlib import use as mpluse
mpl.use("Agg", warn=False)
#mpluse("Qt4Agg")
#mpluse("TkAgg")
import matplotlib.pyplot as plt
plt.rc('font',**{'family':'serif','serif':['Computer Modern Roman']})
plt.rc('text', usetex=True)

#import matplotlib.mlab as mlab
import matplotlib.cm as cm
from  matplotlib import colors as mc
#from mpl_toolkits.axes_grid1 import make_axes_locatable
from matplotlib.ticker import ScalarFormatter

script_path = os.path.realpath(__file__)
module_path = os.path.abspath(os.path.join(os.path.dirname(script_path), '../'))
sys.path.append(module_path)

import ancillary as anc
import constants as cst
# from pytrades_lib import pytrades
#import pytrades_lib.pytrades as pytrades

#import emcee.interruptible_pool as emcee_pool
#import multiprocessing as mp

# ==============================================================================

def compute_K():
  
  Ms_sun = 1.030
  
  dM_ear = 5.
  Mi_nep = 1.
  Mi_ear = Mi_nep * cst.Mneps*cst.Msear
  Mf_jup = 3.
  Mf_ear = Mf_jup * cst.Mjear
  Mgrid_ear = np.arange(Mi_ear, Mf_ear+dM_ear*0.5, dM_ear)
  Mgrid_jup = Mgrid_ear * cst.Mejup
  n_M = np.shape(Mgrid_ear)[0]
  
  dP_day = 1.
  Pi_day = 40.
  Pf_day = 500.
  Pgrid_day = np.arange(Pi_day, Pf_day+dP_day*0.5, dP_day)
  n_P = np.shape(Pgrid_day)[0]
  
  
  mesh_Mj = np.repeat(Mgrid_jup, n_P).reshape((n_M, n_P))
  mesh_Pd = np.tile(Pgrid_day, n_M).reshape((n_M, n_P))
  mesh_rv = anc.compute_Kms(Ms_sun, mesh_Mj, 88.55, mesh_Pd, 0.)
  
  rvmin = np.min(np.min(mesh_rv))
  rvmax = np.max(np.max(mesh_rv))
  print(' | RV min = %.6f m/s max = %.6f m/s' %(rvmin, rvmax))
  
  fig = plt.figure(figsize=(12,12))
  ax = plt.subplot2grid((1, 1), (0, 0))
  ax.set_xlabel('$P [\\textrm{days}]$')
  ax.set_ylabel('$M [M_\\textrm{Jup}]$')
  
  cnrv = ax.contourf(mesh_Pd, mesh_Mj, mesh_rv, 
                    levels=[1., 2., 3., 5., 25., 50., rvmax],
                    #cmap=greys_r_t,
                    #cmap=cm.Greys,
                    cmap=cm.viridis,
                    norm=mc.LogNorm(),
                    #vmin=rvmin, vmax=rvmax,
                    zorder=6
                    )
  print(' | mesh/contourf rv done')
  ax.tick_params(direction='in')
  
  # 1 Msat
  print('1 Msat = %.6f Mjup' %(cst.Msats*cst.Msjup))
  ax.axhline(cst.Msats*cst.Msjup, color='black', ls='-', lw=1.5, zorder=7)
  ax.text(Pgrid_day[1], cst.Msats*cst.Msjup+0.01, r'1 Msat == %.3f Mjup' %(cst.Msats*cst.Msjup), fontsize=10, zorder=7)
  
  # colorbar
  cbar = plt.colorbar(cnrv, ax=ax, orientation='vertical', pad=0.01, format='%5.2f')
  cbar.ax.set_ylabel('$K_\\textrm{RV}$ m/s', size=8)
  cbar.ax.tick_params(labelsize=8, direction='in')
  
  
    
  out_file = os.path.join('/home/borsato/Dropbox/CHEOPS/CHEOPS_Themes/Explore/ToMWarmj/wasp-8/plots_grid/rv_grid', 'rv_grid.png')
  fig.savefig('%s' %(out_file), bbox_inches='tight', dpi=300)
  plt.close(fig)
  print(' Saved plot to %s' %(out_file))
  
  return

# ==============================================================================

# ==============================================================================
# RUN MAIN

if __name__ == "__main__":
  compute_K()

# ==============================================================================

