#!/usr/bin/env python
# -*- coding: utf-8 -*-

# no more "zero" integer division bugs!:P
# import argparse
import warnings
# warnings.simplefilter("ignore", np.RankWarning)
warnings.filterwarnings('ignore')
import os
import sys
import numpy as np  # array

from . import ancillary as anc

import h5py
from scipy.stats import norm as scipy_norm

import pygtc

import matplotlib as mpl

# mpl.use("Agg")
import matplotlib.pyplot as plt

# matplotlib rc params
anc.set_rcParams()

import matplotlib.cm as cm
from matplotlib.ticker import FormatStrFormatter


def set_xaxis(
    ax,
    label_size,
    label_separation,
    label_pad,
    ticklabel_size,
    kel_label,
    ticks_formatter,
    tick_fmt="%.4f",
):
    ax.get_xaxis().set_visible(True)
    ax.xaxis.set_tick_params(labelsize=ticklabel_size)
    ax.xaxis.set_label_coords(0.5, label_separation)
    # ax.ticklabel_format(style='plain', axis='both', useOffset=False)
    plt.setp(ax.xaxis.get_majorticklabels(), rotation=70)
    # ax.xaxis.labelpad = label_pad
    ax.set_xlabel(kel_label, fontsize=label_size, rotation=70.0)
    # tick_step = (ticks_formatter[1] - ticks_formatter[0]) / ticks_formatter[2]
    # ax.xaxis.set_ticks(np.arange(ticks_formatter[0], ticks_formatter[1], tick_step))
    ax.xaxis.set_ticks(np.linspace(ticks_formatter[0], ticks_formatter[1], ticks_formatter[2], endpoint=True))
    tick_formatter = FormatStrFormatter(tick_fmt)
    ax.xaxis.set_major_formatter(tick_formatter)
    return


def set_yaxis(
    ax,
    label_size,
    label_separation,
    label_pad,
    ticklabel_size,
    kel_label,
    ticks_formatter,
    tick_fmt="%.4f",
):
    ax.get_yaxis().set_visible(True)
    ax.yaxis.set_tick_params(labelsize=ticklabel_size)
    ax.yaxis.set_label_coords(label_separation, 0.5)
    # ax.ticklabel_format(style='plain', axis='both', useOffset=False)
    # ax.yaxis.labelpad = label_pad
    ax.set_ylabel(kel_label, fontsize=label_size, rotation=0.0)
    # tick_step = (ticks_formatter[1] - ticks_formatter[0]) / ticks_formatter[2]
    # ax.yaxis.set_ticks(np.arange(ticks_formatter[0], ticks_formatter[1], tick_step))
    ax.yaxis.set_ticks(np.linspace(ticks_formatter[0], ticks_formatter[1], ticks_formatter[2], endpoint=True))
    tick_formatter = FormatStrFormatter(tick_fmt)
    ax.yaxis.set_major_formatter(tick_formatter)
    return

def plot_triangle(cli, log_folder, plot_folder, fitting_posterior, fitting_names, fitting_minmax):
    
    os.makedirs(log_folder, exist_ok=True)
    log_file = os.path.join(log_folder, "log_correlation_plot_fitted.txt")
    olog = open(log_file, 'w')

    anc.print_both("", output=olog)
    anc.print_both(" ================== ", output=olog)
    anc.print_both(" FITTED - CORRELATION PLOTS", output=olog)
    anc.print_both(" ================== ", output=olog)
    anc.print_both("", output=olog)

    nfit = len(fitting_names)

    label_size=12
    ticklabel_size=6
    label_separation=-0.5
    label_pad=12
    label_separation -= (0.1 * (nfit - 2))
    minl, maxl, dl = 2, label_size, 0.6
    label_size = max(minl, maxl - dl*nfit)

    mintl, maxtl, dtl = 2, ticklabel_size, 0.3
    ticklabel_size = max(mintl, maxtl - dtl*nfit)
    # computes mass conversion factor
    m_factor, m_unit, r_factor, r_unit = anc.mass_radius_type_factor(mtype=cli.m_type)
    # set label and legend names
    kel_plot_labels = anc.keplerian_legend(fitting_names, cli.m_type)
    k = anc.get_auto_bins(fitting_posterior)

    anc.print_both("Set overplot ...", output=olog)
    overp_leg, read_fit, overp_fit = None, None, None
    overplot = anc.set_overplot(cli.overplot)
    if overplot is not None:
        overp_leg = overplot.split("sim_")[1]
        try:
            ## OPEN summary_parameters.hdf5 FILE
            with h5py.File(
                os.path.join(cli.full_path, "summary_parameters.hdf5"), "r"
            ) as s_h5f:
                # take only the selected sample
                read_fit = s_h5f["parameters/{:s}/fitted/parameters".format(overplot)][...]

            # fitted parameters has always Mp/Ms in Msun/Mstar, so it is needed to rescale it properly
            overp_fit = read_fit.copy()
        except:
            # reset to None in case cannot find parameters
            overplot = None

    plot_posterior = fitting_posterior.copy()
    plot_minmax    = fitting_minmax.copy()
    for ifit, name in enumerate(fitting_names):
        if "Ms" in name:
            plot_posterior[:,ifit] *= m_factor
            plot_minmax[ifit, :] *= m_factor
            if overp_fit is not None:
                overp_fit[ifit] *= m_factor
        elif (name[0] in ["w"]) or (name[0:2] == "mA" or ("lambda" in name)):
            plot_posterior[:,ifit] %= 360.0
            plot_minmax[ifit, :] %= 360.0
            if overp_fit is not None:
                overp_fit[ifit] %= 360.0

    if cli.corner_type.lower() in ["custom", "both"]:

        anc.print_both("Custom Correlation Plot ...", output=olog)
        sys.stdout.flush()

        fig = plt.figure()
        fig.subplots_adjust(hspace=0.05, wspace=0.05)

        for ix in range(0, nfit, 1):
            xname = fitting_names[ix]
            x_data = plot_posterior[:, ix]

            x_min, x_max = anc.compute_limits(x_data, 0.05)
            if x_min == x_max:
                x_min = plot_minmax[ix, 0]
                x_max = plot_minmax[ix, 1]

            for iy in range(nfit - 1, -1, -1):
                yname = fitting_names[iy]
                if "Ms" in yname:
                    convy = m_factor
                else:
                    convy = 1.0
                y_data = plot_posterior[:, iy]
                y_min, y_max = anc.compute_limits(y_data, 0.05)
                if y_min == y_max:
                    y_min = plot_minmax[iy, 0]
                    y_max = plot_minmax[iy, 1]

                if iy > ix:  # correlation plot
                    anc.print_both(
                        "{:s} vs {:s}".format(xname, yname),
                        output=olog
                    )

                    ax = plt.subplot2grid((nfit + 1, nfit), (iy, ix))

                    new_k = k
                    hist2d_counts_2, xedges_2, yedges_2 = np.histogram2d(
                        x_data,
                        y_data,
                        bins=new_k,
                        range=[[x_data.min(), x_data.max()], [y_data.min(), y_data.max()]],
                        density=True,
                    )
                    xedges = xedges_2
                    yedges = yedges_2

                    x_bins = [
                        0.5 * (xedges_2[i] + xedges_2[i + 1]) for i in range(0, new_k)
                    ]
                    y_bins = [
                        0.5 * (yedges_2[i] + yedges_2[i + 1]) for i in range(0, new_k)
                    ]

                    nl = 3
                    ax.contour(
                        x_bins,
                        y_bins,
                        hist2d_counts_2.T,
                        nl,
                        cmap=cm.viridis,
                        linestyles="solid",
                        linewidths=0.5,
                    )

                    if overplot is not None:
                        # plot selected overplot sample
                        ax.axvline(overp_fit[ix], color="C0", ls="--", lw=1.1, alpha=0.5)
                        ax.axhline(overp_fit[iy], color="C0", ls="--", lw=1.1, alpha=0.5)

                    ax.get_xaxis().set_visible(False)
                    ax.get_yaxis().set_visible(False)
                    if iy == nfit - 1:
                        set_xaxis(
                            ax,
                            label_size,
                            label_separation-0.1,
                            label_pad,
                            ticklabel_size,
                            kel_plot_labels[ix],
                            [xedges[0], xedges[-1], 3],
                        )
                    if ix == 0:
                        set_yaxis(
                            ax,
                            label_size,
                            label_separation-0.1,
                            label_pad,
                            ticklabel_size,
                            kel_plot_labels[iy],
                            [yedges[0], yedges[-1], 3],
                        )

                    ax.set_ylim([y_min, y_max])
                    ax.set_xlim([x_min, x_max])
                    # plt.draw()

                elif iy == ix:  # distribution plot
                    
                    anc.print_both(
                        "{:s} distribution".format(xname),
                        output=olog
                    )
                    ax = plt.subplot2grid((nfit + 1, nfit), (ix, ix))
                    if ix == nfit - 1:
                        hist_orientation = "horizontal"
                    else:
                        hist_orientation = "vertical"

                    idx = np.argsort(x_data)

                    # HISTOGRAM
                    hist_counts, edges, patces = ax.hist(
                        x_data,
                        bins=k,
                        range=[x_data.min(), x_data.max()],
                        histtype="stepfilled",
                        color="darkgrey",
                        # edgecolor='lightgray',
                        edgecolor="None",
                        align="mid",
                        orientation=hist_orientation,
                        density=True,
                        stacked=True,
                    )

                    if overplot is not None:
                        # plot selected overplot sample
                        ax.axhline(
                            overp_fit[ix], color="C0", ls="--", lw=1.1, alpha=0.5
                        )
                    if ix == nfit - 1:
                        ax.set_ylim([y_min, y_max])
                    else:
                        ax.set_xlim([x_min, x_max])

                    ax.get_xaxis().set_visible(False)
                    ax.get_yaxis().set_visible(False)
                    ax.set_title(kel_plot_labels[ix], fontsize=label_size)

                    # plt.draw()
                sys.stdout.flush()

        triangle_file = os.path.join(plot_folder, "emcee_triangle_fitted.png")
        anc.print_both("Saving file {:s}".format(triangle_file.replace(".png", "")), output=olog)
        fig.savefig(triangle_file, bbox_inches="tight", dpi=300)
        triangle_file = os.path.join(plot_folder, "emcee_triangle_fitted.pdf")
        fig.savefig(triangle_file, bbox_inches="tight", dpi=96)
        plt.close(fig)

    if cli.corner_type.lower() in ["pygtc", "both"]:

        anc.print_both("Pygtc Correlation Plot ...", output=olog)
        sys.stdout.flush()

        GTC = pygtc.plotGTC(
            chains=plot_posterior,
            paramNames=kel_plot_labels,
            nContourLevels=3,
            nBins=k,
            truths=overp_fit,
            truthLabels=(overp_leg),
            figureSize=plt.rcParams["figure.figsize"][0],
            mathTextFontSet=plt.rcParams["mathtext.fontset"],
            customLabelFont={"family": plt.rcParams["font.family"], "size": label_size},
            customTickFont={"family": plt.rcParams["font.family"], "size": ticklabel_size},
            customLegendFont={"family": plt.rcParams["font.family"], "size": label_size},
            legendMarker='All',
            labelRotation=(True, False),
        )
        axs = GTC.axes
        for ax in axs:
            ax.tick_params(
                direction='inout',
                pad=4,
                size=3,
                labelsize=ticklabel_size
            )
            lb = ax.get_xlabel()
            if lb != "":
                ax.xaxis.set_label_coords(0.5, label_separation)
                ax.set_xlabel(lb, fontsize=label_size, rotation=45.0)
            
            lb = ax.get_ylabel()
            if lb != "":
                ax.yaxis.set_label_coords(label_separation, 0.5)
                ax.set_ylabel(lb, fontsize=label_size, rotation=45.0)

            for axis in ['top', 'bottom', 'left', 'right']:
                    ax.spines[axis].set_linewidth(0.6)

        triangle_file = os.path.join(plot_folder, "emcee_triangle_fitted_pygtc.png")
        anc.print_both("Saving file {:s}".format(triangle_file.replace(".png", "")), output=olog)
        GTC.savefig(triangle_file, bbox_inches="tight", dpi=300)
        triangle_file = os.path.join(plot_folder, "emcee_triangle_fitted_pygtc.pdf")
        GTC.savefig(triangle_file, bbox_inches="tight", dpi=96)
        plt.close(GTC)

    olog.close()

    return

def plot_triangle_from_file(cli):

    # set emcee folder
    emcee_folder = cli.full_path
    log_folder = os.path.join(emcee_folder, "logs")
    os.makedirs(log_folder, exist_ok=True)
    # plot_folder = prepare_plot_folder(working_path)
    emcee_plots = anc.prepare_emcee_plot_folder(cli.full_path)

    # computes mass conversion factor
    m_factor, m_unit, r_factor, r_unit = anc.mass_radius_type_factor(mtype=cli.m_type)

    # set emcee and trades folder
    emcee_file, emcee_best, folder_best = anc.get_emcee_file_and_best(
        emcee_folder, cli.temp_status
    )

    # get data from the hdf5 file
    (
        fitting_names,
        parameter_boundaries,
        chains,
        acceptance_fraction,
        autocor_time,
        lnprobability,
        ln_err_const,
        completed_steps,
    ) = anc.get_data(emcee_file, cli.temp_status)

    nfit, nwalkers, nruns, nburnin, nruns_sel = anc.get_emcee_parameters(
        chains, cli.temp_status, cli.nburnin, completed_steps
    )

    chains_T_full, parameter_boundaries = anc.select_transpose_convert_chains(
        nfit,
        nwalkers,
        nburnin,
        nruns,
        nruns_sel,
        m_factor,
        fitting_names,
        parameter_boundaries,
        chains,
    )

    chains_T, fitting_posterior, lnprob_burnin, thin_steps = anc.thin_the_chains(
        cli.use_thin,
        nburnin,
        nruns,
        nruns_sel,
        autocor_time,
        chains_T_full,
        lnprobability,
        burnin_done=False,
    )

    fitting_posterior = anc.fix_lambda(fitting_posterior, fitting_names)


    plot_triangle(cli, log_folder, emcee_plots, fitting_posterior, fitting_names, parameter_boundaries)

    return

# =============================================================================
def main():

    cli = anc.get_args()

    plot_triangle_from_file(cli)

    return

# =============================================================================
if __name__ == "__main__":
    main()
