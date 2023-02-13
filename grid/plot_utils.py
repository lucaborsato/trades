#!/usr/bin/env python
# -*- coding: utf-8 -*-

import numpy as np  # array
import matplotlib.pyplot as plt
import copy
import os
import sys
from scipy import ndimage

# CUSTOM MODULES
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
                self.ATTV_MC, self.PTTV_MC = np.genfromtxt(
                    TTV_file, unpack=True)

# RV data


class RVdata:
    def __init__(self, RV_file):
        self.sim_num, self.Mpert, self.Ppert, self.Krv_mps = np.genfromtxt(
            RV_file, unpack=True)

# Plot/s


def get_limits(v, scale):
    vmin = np.min(v)
    vmax = np.max(v)
    dv = vmax - vmin
    lim_min = vmin - dv * scale
    lim_max = vmax + dv * scale
    return [lim_min, lim_max]


def scatter_plot_TTV(TTV):

    fig = plt.figure()

    vir_map = copy.copy(plt.cm.get_cmap("viridis_r"))

    plt.scatter(TTV.Ppert, TTV.Mpert, c=TTV.ATTV * cst.day2min,
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


def grid_plot(
    TTV,
    RV=None,
    grid_on=True,
    rchisq=False,
    Mlim=None,
    Plim=None,
    ATTV_MC=False,
    maxATTVmin=None,
    deltaATTVmin=1.0,
    rvplot=True,
    deltaRVmps=2.0,
    apply_filter=True,
    kwargs_filter={
        "filter": "percentile",
        "size": 5,
        "percentile_value": 50,
    },
    Pval=None,
    Mval=None,
):

    n_P = np.shape(np.unique(TTV.Ppert))[0]
    n_M = np.shape(np.unique(TTV.Mpert))[0]

    mesh_P = TTV.Ppert.reshape((n_M, n_P))
    mesh_M = TTV.Mpert.reshape((n_M, n_P))
    if ATTV_MC:
        mesh_TTV = TTV.ATTV_MC.reshape((n_M, n_P)) * cst.day2min
    else:
        mesh_TTV = TTV.ATTV.reshape((n_M, n_P)) * cst.day2min
    mesh_c2r = TTV.rchisqr.reshape((n_M, n_P))

    if apply_filter:
        if kwargs_filter["filter"].lower() in ["percentile", "perc"]:
            # mesh_TTV = ndimage.median_filter(mesh_TTV, size=size_filter)
            # mesh_c2r = ndimage.median_filter(mesh_c2r, size=size_filter)
            size_filter = kwargs_filter["size"]
            percentile_val = kwargs_filter["percentile_value"]
            mesh_TTV = ndimage.percentile_filter(mesh_TTV, percentile=percentile_val, size=size_filter)
            mesh_c2r = ndimage.percentile_filter(mesh_c2r, percentile=percentile_val, size=size_filter)
        elif kwargs_filter["filter"].lower() == "gaussian":
            sigma=kwargs_filter["sigma"]
            mesh_TTV = ndimage.gaussian_filter(mesh_TTV, sigma=sigma)
            mesh_c2r = ndimage.gaussian_filter(mesh_c2r, sigma=sigma)
        elif kwargs_filter["filter"].lower() == "rank":
            size = kwargs_filter["size"]
            rank = kwargs_filter["rank"]
            mesh_TTV = ndimage.rank_filter(mesh_TTV, rank=rank, size=size)
            mesh_c2r = ndimage.rank_filter(mesh_c2r, rank=rank, size=size)
        elif kwargs_filter["filter"].lower() == "uniform":
            size = kwargs_filter["size"]
            mesh_TTV = ndimage.rank_filter(mesh_TTV, size=size)
            mesh_c2r = ndimage.rank_filter(mesh_c2r, size=size)
        elif kwargs_filter["filter"].lower() == "max":
            threshold_ATTV = kwargs_filter["threshold_ATTVmin"]
            threshold_c2r = kwargs_filter["threshold_c2r"]
            sel = mesh_TTV > threshold_ATTV
            mesh_TTV[sel] = threshold_ATTV
            sel = mesh_c2r > threshold_c2r
            mesh_c2r[sel] = threshold_c2r

    fig, ax = plt.subplots()

    # TTV amplitude in min
    # TTVlevels_min = [float(i) for i in range(
    #     0, 11)] + [np.max(TTV.ATTV)*cst.day2min]
    if maxATTVmin is None:
        if ATTV_MC:
            maxATTVmin = np.max(TTV.ATTV_MC) * cst.day2min
        else:
            maxATTVmin = np.max(TTV.ATTV) * cst.day2min
    if maxATTVmin < deltaATTVmin:
        TTVlevels_min = [0.0, maxATTVmin]
    else:
        TTVlevels_min = np.arange(0.0, maxATTVmin + deltaATTVmin, deltaATTVmin)

    vir_map = copy.copy(plt.cm.get_cmap("viridis_r"))
    vir_map.set_under("white")
    vir_map.set_over("black")

    nTTVlev = np.shape(TTVlevels_min)[0]
    TTV_colors = vir_map(np.linspace(0.2, 1.0, num=nTTVlev, endpoint=True))

    TTV_ax = ax.contourf(
        mesh_P, mesh_M, mesh_TTV,
        levels=TTVlevels_min,
        colors=TTV_colors,
        origin='lower',
        extend='max',
        zorder=5
    )

    TTV_bar = fig.colorbar(TTV_ax)
    TTV_bar.set_label(r'A$_\mathrm{TTV}$ (min)')

    if rvplot and RV is not None:
        # RV amplitude in m/s
        mesh_RV = RV.Krv_mps.reshape((n_M, n_P))
        # RVlevels_mps = [1.0, 2.0, 3.0, 4.0, 5.0, 6.0] # m/s
        # RVlevels_mps = [float(i) for i in range(1, 11)] + [int(0.5*(10.0+np.max(RV.Krv_mps)))]  # m/s
        # RVlevels_mps = [float(i) for i in range(2, 12, 2)] + \
        #     [int(0.5*(10.0+np.max(RV.Krv_mps)))]  # m/s
        maxRVlevel = int(0.5 * (10.0 + np.max(RV.Krv_mps)))
        if maxRVlevel <= 2.0:
            RVlevels_mps = [maxRVlevel]
        elif maxRVlevel > 12.0:
            RVlevels_mps = [float(i) for i in range(
                int(deltaRVmps), 12, int(deltaRVmps))] + [maxRVlevel]
        else:
            RVlevels_mps = np.arange(
                deltaRVmps, maxRVlevel + deltaRVmps, deltaRVmps)
        # print(RVlevels_mps)
        # nRVlev       = len(RVlevels_mps)
        gval = 0.5
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
        clabels = ax.clabel(RV_lines, inline=True,
                            fontsize=plt.rcParams['font.size'] - 4, colors='white', fmt='%.0f m/s')
        for clabel in clabels:
            clabel.set_zorder(8)

    if(rchisq):

        #     # \chi^_r as lines contour
        c2rlevels = [0.0, 1.0, 2.0, 3.0, 5.0]
        nc2rlev = len(c2rlevels)

        gval = 0.3
        c2r_colors = [(gval, gval, gval, i)
                      for i in np.linspace(0.0, 0.7, num=nc2rlev, endpoint=True)]
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

    if(grid_on):
        ax.plot(TTV.Ppert, TTV.Mpert, color='black', marker='o',
                ms=1.5, mec='None', ls='', alpha=0.5, zorder=8)

    if Pval is not None:
        ax.axvline(
            Pval,
            color='C1',
            ls='-',
            lw=0.8,
            alpha=0.7,
            zorder=9
        )
    if Mval is not None:
        ax.axhline(
            Mval,
            color='C1',
            ls='-',
            lw=0.8,
            alpha=0.7,
            zorder=9
        )

    scale = 0.01
    if Mlim is None or len(Mlim) <= 1:
        ax.set_ylim(get_limits(TTV.Mpert, scale))
    else:
        ax.set_ylim(Mlim)
    if Plim is None or len(Plim) <= 1:
        ax.set_xlim(get_limits(TTV.Ppert, scale))
    else:
        ax.set_xlim(Plim)

    ax.set_xlabel(r"$P_\mathrm{perturber}$ (days)")
    ax.set_ylabel(r"$M_\mathrm{perturber}$ ($M_\oplus$)")
    plt.tight_layout()
    plt.show()

    return fig

# ==============================================================================
# ==============================================================================
