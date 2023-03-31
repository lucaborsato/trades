#!/usr/bin/env python
# -*- coding: utf-8 -*-

# no more "zero" integer division bugs!:P
# import argparse
import os
import sys
import numpy as np  # array
import h5py
import ancillary as anc
# from scipy.stats import norm as scipy_norm
import pygtc
import matplotlib.pyplot as plt

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


## PLOT PHYSICAL POSTERIOR CORRELATION PLOT
def plot_triangle(cli, log_folder, plot_folder, physical_posterior, physical_names):

    os.makedirs(log_folder, exist_ok=True)
    log_file = os.path.join(log_folder, "log_correlation_plot_physical.txt")
    olog = open(log_file, 'w')

    anc.print_both("", output=olog)
    anc.print_both(" ================== ", output=olog)
    anc.print_both(" PHYSICAL - CORRELATION PLOTS", output=olog)
    anc.print_both(" ================== ", output=olog)
    anc.print_both("", output=olog)

    nphy = len(physical_names)

    label_size=12
    ticklabel_size=6
    label_separation=-0.5
    label_pad=12
    label_separation -= (0.1 * (nphy - 2))
    minl, maxl, dl = 2, label_size, 0.6
    label_size = max(minl, maxl - dl*nphy)

    mintl, maxtl, dtl = 2, ticklabel_size, 0.3
    ticklabel_size = max(mintl, maxtl - dtl*nphy)
    # computes mass conversion factor
    # m_factor, m_unit, r_factor, r_unit = anc.mass_radius_type_factor(mtype=cli.m_type)
    # set label and legend names
    physical_labels = anc.derived_labels(physical_names, cli.m_type)
    k = anc.get_auto_bins(physical_posterior)

    anc.print_both("Set overplot ...", output=olog)
    overp_leg, read_phy, overp_phy = None, None, None
    overplot = anc.set_overplot(cli.overplot)
    if overplot is not None:
        overp_leg = overplot.split("sim_")[1]
        try:   
            ## OPEN summary_parameters.hdf5 FILE
            with h5py.File(os.path.join(cli.full_path, "summary_parameters.hdf5"), "r") as s_h5f:
                # take only the selected sample
                read_phy = s_h5f["parameters/{:s}/physical/parameters".format(overplot)][...]

            overp_phy = read_phy.copy()
        except:
            # reset to None in case cannot find parameters
            overplot = None

    plot_posterior = physical_posterior.copy()

    # for iphy, name in enumerate(physical_names):
    #     if (name[0] in ["w"]) or (name[0:2] == "mA" or ("lambda" in name)):
    #         plot_posterior[:,iphy] %= 360.0
    #         if overp_phy is not None:
    #             overp_phy[iphy] %= 360.0

    if cli.corner_type.lower() in ["custom", "both"]:

        anc.print_both("Custom Correlation Plot ...", output=olog)
        sys.stdout.flush()

        fig = plt.figure(figsize=(6, 6))
        fig.subplots_adjust(hspace=0.05, wspace=0.05)

        for ix in range(0, nphy):
            xname = physical_names[ix]
            x_data = plot_posterior[:, ix]
            x_min, x_max = anc.compute_limits(x_data, 0.05)
            if x_min == x_max:
                x_min = x_min - 1.0
                x_max = x_max + 1.0

            for iy in range(0, nphy):
                yname = physical_names[iy]
                y_data = plot_posterior[:, iy]
                y_min, y_max = anc.compute_limits(y_data, 0.05)
                if y_min == y_max:
                    y_min = y_min - 1.0
                    y_max = y_max + 1.0

                if iy > ix:  # correlation plot
                    anc.print_both(
                        "{:s} vs {:s}".format(xname, yname),
                        output=olog
                    )

                    ax = plt.subplot2grid((nphy + 1, nphy), (iy, ix))

                    # hist2d_counts, xedges, yedges, image2d = ax.hist2d(
                    #     x_data,
                    #     y_data,
                    #     bins=k,
                    #     range=[[x_data.min(), x_data.max()], [y_data.min(), y_data.max()]],
                    #     cmap=cm.gray_r,
                    #     density=False,
                    # )

                    # new_k = int(k/3)
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
                        ax.axhline(overp_phy[ix], color="C0", ls="--", lw=1.1, alpha=0.5)
                        ax.axhline(overp_phy[iy], color="C0", ls="--", lw=1.1, alpha=0.5)

                    ax.get_xaxis().set_visible(False)
                    ax.get_yaxis().set_visible(False)
                    if iy == nphy - 1:
                        set_xaxis(
                            ax,
                            label_size,
                            label_separation,
                            label_pad,
                            ticklabel_size,
                            physical_labels[ix],
                            [xedges[0], xedges[-1], 3],
                        )
                    if ix == 0:
                        set_yaxis(
                            ax,
                            label_size,
                            label_separation,
                            label_pad,
                            ticklabel_size,
                            physical_labels[iy],
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

                    ax = plt.subplot2grid((nphy + 1, nphy), (ix, ix))
                    if ix == nphy - 1:
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
                            overp_phy[ix], color="C0", ls="--", lw=1.1, alpha=0.5
                        )
                    if ix == nphy - 1:
                        ax.set_ylim([y_min, y_max])
                    else:
                        ax.set_xlim([x_min, x_max])

                    ax.get_xaxis().set_visible(False)
                    ax.get_yaxis().set_visible(False)
                    ax.set_title(physical_labels[ix], fontsize=label_size)

                    # plt.draw()
                sys.stdout.flush()

        correlation_file = os.path.join(plot_folder, "emcee_triangle_physical.png")
        anc.print_both("Saving file {:s}".format(correlation_file.replace(".png", "")), output=olog)
        fig.savefig(correlation_file, bbox_inches="tight", dpi=300)
        correlation_file = os.path.join(plot_folder, "emcee_triangle_physical.pdf")
        fig.savefig(correlation_file, bbox_inches="tight", dpi=96)
        plt.close(fig)


    if cli.corner_type.lower() in ["pygtc", "both"]:

        anc.print_both("Pygtc Correlation Plot ...", output=olog)
        sys.stdout.flush()

        GTC = pygtc.plotGTC(
            chains=plot_posterior,
            paramNames=physical_labels,
            nContourLevels=3,
            nBins=k,
            truths=overp_phy,
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

        emcee_fig_file = os.path.join(plot_folder, "emcee_triangle_physical_pygtc.png")
        anc.print_both("Saving file {:s}".format(emcee_fig_file.replace(".png", "")), output=olog)
        GTC.savefig(emcee_fig_file, bbox_inches="tight", dpi=300)
        emcee_fig_file = os.path.join(plot_folder, "emcee_triangle_physical_pygtc.pdf")
        GTC.savefig(emcee_fig_file, bbox_inches="tight", dpi=96)
        plt.close(GTC)

    olog.close()

    return

def plot_triangle_from_file(cli):

    # set emcee folder
    emcee_folder = cli.full_path
    log_folder = os.path.join(emcee_folder, "logs")
    emcee_plots = anc.prepare_emcee_plot_folder(cli.full_path)

    # read physical posterior file
    posterior_file = os.path.join(cli.full_path, "posterior.hdf5")
    with h5py.File(posterior_file, "r") as h5f:
        physical_names = anc.decode_list(h5f["physical_names"][...])
        physical_posterior_in = np.array(h5f["posterior_physical"], dtype=np.float64)

    physical_posterior = anc.derived_posterior_check(physical_names, physical_posterior_in)

    plot_triangle(cli, log_folder, emcee_plots, physical_posterior, physical_names)

    return

# =============================================================================
def main():

    cli = anc.get_args()

    plot_triangle_from_file(cli)

    return

# =============================================================================
#
# Here run the proper function
#
if __name__ == "__main__":
    main()
