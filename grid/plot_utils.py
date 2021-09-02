#!/usr/bin/env python
# -*- coding: utf-8 -*-

import numpy as np # array
import matplotlib.pyplot as plt
import copy
import os
import sys

script_path = os.path.realpath(__file__)
module_path = os.path.abspath(os.path.join(os.path.dirname(script_path),
                                           '../pytrades'))
sys.path.append(module_path)

import constants as cst
# ==============================================================================
# ==============================================================================
# READ AND PLOT FUNCTIONS

# TTV data
class TTVdata:
  def __init__(self, TTV_file, old_type=False):
    if old_type:
      self.sim_num, _, _, \
        self.Mpert, self.Ppert, _, _, \
          self.ATTV, self.PTTV, _, \
            self.dof, self.rchisqr = np.genfromtxt(TTV_file, unpack=True)
      self.ATTV_MC, self.PTTV_MC = self.ATTV, self.PTTV
    else:
      self.sim_num, _, _, \
        self.Mpert, self.Ppert, _, _, \
          self.ATTV, self.PTTV, _, self.dof, self.rchisqr, \
            self.ATTV_MC, self.PTTV_MC = np.genfromtxt(TTV_file, unpack=True)

# RV data
class RVdata:
  def __init__(self, RV_file):
    self.sim_num, self.Mpert, self.Ppert, self.Krv_mps = np.genfromtxt(RV_file, unpack=True)

# Plot/s
def get_limits(v, scale):
  vmin = np.min(v)
  vmax = np.max(v)
  dv   = vmax - vmin
  lim_min = vmin - dv * scale
  lim_max = vmax + dv * scale
  return [lim_min, lim_max]

def scatter_plot_TTV(TTV):
  
  fig = plt.figure()
  
  vir_map = copy.copy(plt.cm.get_cmap("viridis_r"))
  
  plt.scatter(TTV.Ppert, TTV.Mpert, c=TTV.ATTV*cst.day2min,
              vmin=0.0,
              vmax=10.0,
              cmap=vir_map,
             )
  plt.colorbar(label=r'A$_\mathrm{TTV}$ (min)')
  plt.xlabel(r"$P_\mathrm{perturber}$ (days)")
  plt.ylabel(r"$M_\mathrm{perturber}$ ($M_\oplus$)")
  plt.show()
  return fig

# ======================================================================

def grid_plot(TTV, RV, grid_on=True, rchisq=False):
    
  n_P = np.shape(np.unique(TTV.Ppert))[0]
  n_M = np.shape(np.unique(TTV.Mpert))[0]
  
  mesh_P   = TTV.Ppert.reshape((n_M, n_P))
  mesh_M   = TTV.Mpert.reshape((n_M, n_P))
  mesh_TTV = TTV.ATTV.reshape((n_M, n_P))*cst.day2min
  mesh_RV  = RV.Krv_mps.reshape((n_M, n_P))
  mesh_c2r = TTV.rchisqr.reshape((n_M, n_P))
  
  fig, ax = plt.subplots()

  # TTV amplitude in min
  TTVlevels_min = [np.float(i) for i in range(0,11)] + [np.max(TTV.ATTV)*cst.day2min]
  vir_map = copy.copy(plt.cm.get_cmap("viridis_r"))
  vir_map.set_under("white")
  vir_map.set_over("black")
  
  nTTVlev = np.shape(TTVlevels_min)[0]
  TTV_colors = vir_map(np.linspace(0.2,1.0,num=nTTVlev,endpoint=True))
  
  TTV_ax = ax.contourf(mesh_P, mesh_M, mesh_TTV,
                       levels=TTVlevels_min,
                       colors=TTV_colors,
                       origin=None,
                       extend='max',
                       zorder=5
                      )
  
  TTV_bar = fig.colorbar(TTV_ax)
  TTV_bar.set_label(r'A$_\mathrm{TTV}$ (min)')
  
  # RV amplitude in m/s
#   RVlevels_mps = [1.0, 2.0, 3.0, 4.0, 5.0, 6.0] # m/s
  # RVlevels_mps = [float(i) for i in range(1, 11)] + [int(0.5*(10.0+np.max(RV.Krv_mps)))]  # m/s
  RVlevels_mps = [float(i) for i in range(2, 12, 2)] + [int(0.5*(10.0+np.max(RV.Krv_mps)))]  # m/s
#   print(RVlevels_mps)
  # nRVlev       = len(RVlevels_mps)
  gval         = 0.5
  # RV_colors    = [(gval, gval, gval, i) for i in np.linspace(0.0, 1.0, num=nRVlev, endpoint=True)]
  
#   RV_ax = ax.contourf(mesh_P, mesh_M, mesh_RV,
#                       levels=RVlevels_mps,
#                       colors=RV_colors,
#                       origin=None,
#                       extend='max',
#                       zorder=6
#                      )
  
#   RV_bar = fig.colorbar(RV_ax)
#   RV_bar.set_label(r'$K_\mathrm{RV}$ (m/s)')
  
  RV_lines = ax.contour(mesh_P, mesh_M, mesh_RV,
    levels=RVlevels_mps,
    #                         colors=RV_colors,
    colors='white',
    origin=None,
    extend='max',
    linewidths=0.5,
    antialiased=True,
    zorder=7
  )
  clabels = ax.clabel(RV_lines, inline=True, fontsize=plt.rcParams['font.size']-4, colors='white', fmt='%.0f m/s')
  for clabel in clabels:
    clabel.set_zorder(8)
  
  if(rchisq):
    
#     # \chi^_r as lines contour
    c2rlevels = [0.0, 1.0, 2.0, 3.0, 5.0]
    nc2rlev   = len(c2rlevels)  

    gval         = 0.3
    c2r_colors    = [(gval, gval, gval, i) for i in np.linspace(0.0, 0.7, num=nc2rlev, endpoint=True)]
    c2r = ax.contourf(mesh_P, mesh_M, mesh_c2r,
                       levels=c2rlevels,
                       colors=c2r_colors,
                       origin=None,
                       extend='max',
                       zorder=6
                      )
#     cbaxes  = fig.add_axes([0.0, 1.0, 0.8, 0.03]) # [left, bottom, width, height]
    c2r_bar = fig.colorbar(c2r, orientation="horizontal",
#                            cax=cbaxes
                           anchor=(0.0, 0.0)
    )
    c2r_bar.set_label(r'$\chi^2_r$')
    
  
  if(grid_on): ax.plot(TTV.Ppert, TTV.Mpert, color='black', marker='o', ms=1.5, mec='None', ls='', alpha=0.5, zorder=8)
  
  
  scale = 0.01
  ax.set_xlim(get_limits(TTV.Ppert, scale))
  ax.set_ylim(get_limits(TTV.Mpert, scale))
  
  ax.set_xlabel(r"$P_\mathrm{perturber}$ (days)")
  ax.set_ylabel(r"$M_\mathrm{perturber}$ ($M_\oplus$)")
  plt.tight_layout()
  plt.show()
  
  return fig

# ==============================================================================
# ==============================================================================