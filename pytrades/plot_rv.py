#!/usr/bin/env python
# -*- coding: utf-8 -*-

# no more "zero" integer division bugs!:P
import argparse
import os  #  os: operating system
# import time  # time: execution time
import glob  # glob: globbing file...loading multiple files as *.pippa
# import sys  # sys: system
import numpy as np  # array

import h5py

# import matplotlib as mpl
# mpl.use("Agg")

import matplotlib.pyplot as plt

# custom modules
import constants as cst
import ancillary as anc

anc.set_rcParams()

# ==============================================================================

filled_markers = anc.filled_markers
size_markers = anc.size_markers

CLI_RV = anc.CLI_RV
# ==============================================================================


def get_sim_file(cli):

    main_folder = cli.full_path
    # sim_file = os.path.join(main_folder, '%s_%s_simRV.dat' %(cli.idsim, cli.lmflag))
    sim_file = os.path.join(main_folder, "%d_%d_simRV.dat" % (cli.idsim, cli.lmflag))

    return sim_file

# ==============================================================================
def get_jitter(cli):

    main_folder = cli.full_path
    finalpar_pattern = os.path.join(main_folder, "%d_%d_final*par.dat" % (cli.idsim, cli.lmflag))
    finalpar_file = glob.glob(finalpar_pattern)[0]

    with open(finalpar_file, 'r') as of:
        lines = of.readlines()
    
    jitters = []
    for line in lines:
        if (line.strip()[0] != "#") and ("l2j_" in line):
            l2j = float(line.strip().split()[1])
            jitters.append(np.power(2.0, l2j))
    return jitters

# ==============================================================================
class rv_data:
    def __init__(self, t, RVo, eRVo, rvs, RVs, rv_gamma, rv_trend, jitter, RVset, tscale):
        self.RVset = RVset
        self.t = t
        self.time_plot = self.t - tscale
        self.RVo = RVo
        self.rvo = RVo - rv_gamma
        self.eRVo = eRVo
        self.rvs = rvs
        self.RVs = RVs - rv_gamma
        self.rv_gamma = rv_gamma
        self.rv_trend = rv_trend
        self.jitter = jitter
        self.dRVos = RVo - RVs
        self.nRV = np.shape(self.RVo)[0]
        self.eRV_j = np.sqrt(eRVo**2 + jitter**2)
        self.wres = self.dRVos/self.eRV_j

    def add_jitter(self, jitter):

        self.eRV_j = np.sqrt(self.eRVo**2 + jitter**2)
        self.wres = self.dRVos/self.eRV_j


# ==============================================================================


def get_sim_data(file_in, jitters, tscale=2440000.5):

    # OLD
    # JD 0 RVobs 1 eRVobs 2 rv_sim 3 RV_sim 4 gamma 5 e_gamma 6 RVsetID 7 RV_stat 8
    # itime, iRVo, ieRVo = 0,1,2
    # irvs, iRVs, igamma, itrend, isetid = 3,4,5,6,7

    # NEW
    # JD 0 RVobs 1 eRVobs 2 rv_sim 3 RV_sim 4 rv_gamma 5 rv_trend 6 jitter 7 RVsetID 8 RV_stat 9
    itime, iRVo, ieRVo = 0,1,2
    irvs, iRVs, igamma, itrend, ijitter, isetid = 3,4,5,6,7,8

    sim_in = np.genfromtxt(file_in)
    rvsetid = np.unique(sim_in[:, isetid].astype(int))
    nset = np.shape(rvsetid)[0]

    print("rvsetid: ", rvsetid)
    print("nset = ", nset)
    rv_os = []  # observations and simulations data
    for i in range(0, nset):
        print("create rv set for i = ", i)
        rvsel = sim_in[:, isetid].astype(int) == rvsetid[i]
        print("nrvsel = ", np.sum(rvsel))
        orv = rv_data(
            sim_in[rvsel, itime],
            sim_in[rvsel, iRVo],
            sim_in[rvsel, ieRVo],
            sim_in[rvsel, irvs],
            sim_in[rvsel, iRVs],
            sim_in[rvsel, igamma],
            sim_in[rvsel, itrend],
            sim_in[rvsel, ijitter],
            rvsetid[i],
            tscale
        )
        orv.add_jitter(jitters[i])
        print("added jitter = {:.3f} m/s".format(jitters[i]))
        rv_os.append(orv)

    return rv_os


# ==============================================================================


def get_sim_model(cli):

    # orb_file = os.path.join(cli.full_path, '%s_%s_rotorbit.dat' %(cli.idsim, cli.lmflag))

    # data = np.genfromtxt(orb_file, usecols=(0,1,-1)) # time, lte, rv
    # tmod = data[:,0]+data[:,1]-cli.tscale
    # rvmod = data[:,2]

    tscale = cli.tscale

    mod_file = h5py.File(
        os.path.join(cli.full_path, "{:s}_models.hdf5".format(cli.sim_type)), "r"
    )
    # gr = "{:04d}".format(int(cli.idsim))
    gr = cli.sim_type
    tepoch = mod_file[gr].attrs["tepoch"]
    tmod = mod_file["{:s}/time_rv_mod".format(gr)][...]
    tmod += tepoch - tscale
    rvmod = mod_file["{:s}/rv_mod".format(gr)][...]

    return tmod, rvmod


# ==============================================================================


class rv_sample:
    def __init__(self, tepoch, time_rv_mod, rv_mod, tscale):
        self.tepoch = tepoch
        self.time_rv_mod = time_rv_mod + tepoch
        self.tscale = tscale
        self.time_plot = time_rv_mod + tepoch - tscale
        self.rv_mod = rv_mod


def read_samples(cli, samples_file=None):

    if samples_file is not None:

        sh5 = h5py.File(os.path.abspath(samples_file), "r")
        # n_samples = len(list(sh5.keys()))
        samples = []
        for igr in list(sh5.keys()):  # loop on groups, one for each sample
            tepoch = sh5[igr].attrs["tepoch"]
            trv = sh5["%s/time_rv_mod" % (igr)][...]
            rv = sh5["%s/rv_mod" % (igr)][...]
            samples.append(rv_sample(tepoch, trv, rv, cli.tscale))
        sh5.close()

    else:

        samples = None

    return samples


# ==============================================================================


def set_axis_default(ax, ticklabel_size=10, aspect="equal", labeldown=True):

    ax.ticklabel_format(useOffset=False)
    ax.tick_params(direction="inout", labelsize=ticklabel_size, labelbottom=labeldown)
    ax.set_aspect(aspect)

    return


def set_symmetric_lim(ax, val):

    ax.set_xlim(-val, val)
    ax.set_ylim(-val, val)

    return


def axtitle(ax, labtitle="", fontsize=8):

    ax.text(
        0.5,
        1.02,
        labtitle,
        horizontalalignment="center",
        verticalalignment="bottom",
        fontsize=fontsize,
        transform=ax.transAxes,
    )

    return


# ==============================================================================
set_colors = anc.set_colors
# ==============================================================================

def plot_rv(cli, figsize=(5,5), samples=None, save_plot=True, show_plot=False):

    tscale = cli.tscale
    sim_file = get_sim_file(cli)
    jitters = get_jitter(cli)
    rv_os = get_sim_data(sim_file, jitters, tscale)
    try:
        tmod, rvmod = get_sim_model(cli)
        print("rvmod = [{:.4f}, {:.4f}]".format(np.min(rvmod), np.max(rvmod)))
    except:
        tmod, rvmod = None, None
        print("tmod and rvmod set to None")

    if tscale == 0.0:
        xlabel = "BJD$_\mathrm{TDB}$"
    else:
        xlabel = "BJD$_\mathrm{{TDB}} - {:.3f}$".format(tscale)

    axs = []

    fig = plt.figure(figsize=figsize)
    fig.subplots_adjust(hspace=0.05, wspace=0.25)

    nrows = 3
    ncols = 1

    lfont = plt.rcParams["font.size"]
    tfont = plt.rcParams["xtick.labelsize"]
    dlfont = 2

    print("Plotting ...")

    # ==========================================
    # RV plot
    ax = plt.subplot2grid((nrows, ncols), (0, 0), rowspan=2)
    set_axis_default(ax, ticklabel_size=tfont, aspect="auto", labeldown=False)
    ax.set_ylabel("RV (m/s)", fontsize=lfont)
    axtitle(ax, labtitle="Radial Velocities", fontsize=lfont)

    ax.axhline(0.0, color="black", ls="-", lw=0.7, zorder=2)

    lobs_a = []
    lsim_a = []

    x = []
    y = []
    px = 0.05
    py = 0.05
    
    nset = len(rv_os)
    label_data = ["RV dataset \#{:d}".format(i + 1) for i in range(0, nset)]
    if cli.labels is not None:
        nlabels = len(list(cli.labels))
        label_data[0:nlabels] = [ll for ll in list(cli.labels)]
    n_idsrc = len(label_data)
    ocolors = [""]*n_idsrc

    ## ocolors = set_colors(nset, colormap=cli.color_map)
    # set default colors:
    ocolors = set_colors(n_idsrc, colormap="nipy_spectral")
    # check color_map option
    if isinstance(cli.color_map, str):
        print("color_map is a string")
        if cli.color_map in anc.all_color_maps:
            colormap = cli.color_map
        else:
            print("color_map not found in all_color_maps")
            colormap = "nipy_spectral"
        full_colors = anc.set_colors(n_idsrc, colormap=colormap)
        # print("full_colors = ", full_colors)
        ocolors = full_colors
    elif isinstance(cli.color_map, dict):
        print("color_map is a dictionary")
        for ilab, lab in enumerate(label_data):
            if lab in cli.color_map.keys():
                ocolors[ilab] = cli.color_map[lab]

    zsim = 9
    zobs = 8
    zmod = 7
    zsmp = 6
    zleg = zsim + 1

    cfsm = plt.get_cmap("gray")
    gval = 0.6
    dg = 0.1


    for i in range(0, nset):
        mso = size_markers[i]
        mss = mso - 0.5
        # print(len(rv_os), type(rv_os))
        # print(rv_os)
        rv = rv_os[i]
        # print(type(rv))
        # print(rv)
        chi2 = np.sum(rv.wres*rv.wres)
        print("RVs(set {}) contribution to the chi_square = {:.3f}".format(i, chi2))

        sort_obs = np.argsort(rv.time_plot)
        x.append(np.tile(rv.time_plot[sort_obs], 3))  # for obs limits
        y.append(
            np.concatenate(
                (
                    rv.rvo[sort_obs] - rv.eRVo[sort_obs],
                    rv.rvo[sort_obs] + rv.eRVo[sort_obs],
                    rv.rvs[sort_obs] + rv.rv_trend,
                )
            )
        )  # for obs limits
        # data
        lobs = ax.errorbar(
            rv.time_plot[sort_obs],
            rv.rvo[sort_obs],
            yerr=rv.eRVo[sort_obs],
            color=ocolors[i],
            ecolor=ocolors[i],
            fmt=filled_markers[i],
            ms=mso,
            mec="black",
            mew=0.3,
            ls="",
            elinewidth=0.6,
            capsize=0,
            zorder=zobs,
            label=label_data[i],
        )
        lobs_a.append(lobs)
        # trades
        (lsim,) = ax.plot(
            rv.time_plot[sort_obs],
            rv.rvs[sort_obs]+rv.rv_trend[sort_obs],
            # marker=filled_markers[i], ms=mss, mfc='None',
            color="C0",
            marker=filled_markers[i],
            ms=anc.size_default * 0.5,
            mfc="black",
            # mec=ocolors[i],
            mec="black",
            mew=0.4,
            ls="",
            zorder=zsim,
            # label="simulations",
        )
        lsim_a.append(lsim)

    # prepare limits
    x = np.concatenate(x)
    dx = np.max(x) - np.min(x)
    minx = np.min(x) - px * dx
    maxx = np.max(x) + px * dx
    y = np.concatenate(y)
    dy = np.max(y) - np.min(y)
    miny = np.min(y) - py * dy
    maxy = np.max(y) + py * dy
    # xlims = ax.get_xlim()

    if tmod is not None:
        # plot model and samples
        sort_mod = np.argsort(tmod)
        (lmod,) = ax.plot(
            tmod[sort_mod],
            rvmod[sort_mod],
            color="black",
            marker="None",
            ls="-",
            lw=0.3,
            alpha=1,
            zorder=zmod,
            label="RV model",
        )
    
    xx = x
    yy = y
    
    if samples is not None:
        print("samples plot ... ")
        # # ==== plot all the samples o-c
        # lsmp_a = []
        # # nsmp = len(samples)
        # for smp in samples:
        #     sort_smp = np.argsort(smp.time_plot)
        #     xx = np.concatenate((xx, smp.time_plot[sort_smp]))
        #     yy = np.concatenate((yy, smp.rv_mod[sort_smp]))
        #     (lsmp,) = ax.plot(
        #         smp.time_plot[sort_smp],
        #         smp.rv_mod[sort_smp],
        #         color="gray",
        #         marker="None",
        #         ls="-",
        #         lw=0.3,
        #         alpha=0.8,
        #         zorder=1,
        #         label="RV samples",
        #     )
        #     lsmp_a.append(lsmp)

        # # ==== plot 1/2/3 CI
        if tmod is not None:
            tsort = np.sort(tmod)
        else:
            tsort == np.sort(ss[0].time_plot)
        n_sim = len(tsort)
        print("n_sim (RVs) = {:d}".format(n_sim))
        rv_smp = []
        print("creating rv_smp list ...")
        for ss in samples:
            sort_smp = np.argsort(ss.time_plot)
            rv_smp.append(ss.rv_mod[sort_smp])
        print("rv_smp list to array ... ")
        rv_smp = np.array(rv_smp).T
        print("shape(rv_smp) = ", np.shape(rv_smp))
        c1, c2, c3 = 0.6827, 0.9544, 0.9974
        hc1, hc2, hc3 = c1*0.5, c2*0.5, c3*0.5

        hdi1 = np.percentile(rv_smp, [50 - (100*hc1), 50 + (100*hc1)], axis=1).T
        hdi2 = np.percentile(rv_smp, [50 - (100*hc2), 50 + (100*hc2)], axis=1).T
        hdi3 = np.percentile(rv_smp, [50 - (100*hc3), 50 + (100*hc3)], axis=1).T

        print("shape(hdi1) = ", np.shape(hdi1))
        print("shape(hdi2) = ", np.shape(hdi2))
        print("shape(hdi3) = ", np.shape(hdi3))
        
        print("plot HDI ... ")
        yy = np.concatenate((yy, hdi3[:, 0], hdi3[:, 1]))
        xx = np.concatenate((xx, tsort, tsort))

        ax.fill_between(
            tsort,
            hdi1[:, 0],
            hdi1[:, 1],
            color=cfsm(gval),
            alpha=1.0,
            lw=0.0,
            zorder=zsmp,
        )
        ax.fill_between(
            tsort,
            hdi2[:, 0],
            hdi2[:, 1],
            color=cfsm(gval+dg),
            alpha=1.0,
            lw=0.0,
            zorder=zsmp-1,
        )
        ax.fill_between(
            tsort,
            hdi3[:, 0],
            hdi3[:, 1],
            color=cfsm(gval+(dg*2)),
            alpha=1.0,
            lw=0.0,
            zorder=zsmp-2,
        )

    if cli.limits == "sam":
        dx = np.max(xx) - np.min(xx)
        minx = np.min(xx) - px * dx
        maxx = np.max(xx) + px * dx
        dy = np.max(yy) - np.min(yy)
        miny = np.min(yy) - py * dy
        maxy = np.max(yy) + py * dy
    else:  # if obs & samples not None adjust y limits
        yyy = yy[np.logical_and(xx >= minx, xx <= maxx)]
        dy = np.max(yyy) - np.min(yyy)
        miny = np.min(yyy) - py * dy
        maxy = np.max(yyy) + py * dy

    # ax.set_xlim(xlims)
    ax.set_xlim(minx, maxx)
    ax.set_ylim(miny, maxy)

    # if samples is None and tmod is None:
    #     lhand = lobs_a #+ [lsim_a[0]]
    # elif samples is None and tmod is not None:
    #     lhand = lobs_a + [lmod] # [lsim_a[0], lmod]
    # elif samples is not None and tmod is None:
    #     lhand = lobs_a + [lsmp_a[0]]#[lsim_a[0], lsmp_a[0]]
    # else:
    #     lhand = lobs_a + [lmod, lsmp_a[0]] #[lsim_a[0], lmod, lsmp_a[0]]
    lhand = lobs_a

    # legend in a column right of the boxes
    if "out" in cli.legend:
        ax.legend(handles=lhand,
            loc="center left",
            bbox_to_anchor=(1., 0.5),
            ncol=1,
            fontsize=tfont-dlfont,
        ).set_zorder(zleg)
    else:
        ax.legend(handles=lhand,
            loc="best",
            # bbox_to_anchor=(1., 0.5),
            # ncol=1,
            fontsize=tfont-dlfont,
        ).set_zorder(zleg)

    axs.append(ax)
    
    # ===============================================
    # residuals plot
    ax = plt.subplot2grid((nrows, ncols), (2, 0), rowspan=1)
    set_axis_default(ax, ticklabel_size=tfont, aspect="auto")
    ax.set_ylabel("res (m/s)", fontsize=lfont)
    ax.set_xlabel(xlabel, fontsize=lfont)

    ax.axhline(0.0, color="black", ls="-", lw=0.7, zorder=1)


    for i in range(0, nset):
        mso = size_markers[i]
        # mss = mso - 0.5
        rv = rv_os[i]
        sort_obs = np.argsort(rv.time_plot)
        ax.errorbar(
            rv.time_plot[sort_obs],
            rv.dRVos[sort_obs],
            yerr=rv.eRVo[sort_obs],
            color=ocolors[i],
            fmt=filled_markers[i],
            ms=mso,
            mec="black",
            mew=0.3,
            ls="",
            ecolor=ocolors[i],
            elinewidth=0.6,
            capsize=0,
            zorder=zobs,
        )
        # error only, with jitter
        cjb = 'black'
        ax.errorbar(
            rv.time_plot[sort_obs],
            rv.dRVos[sort_obs],
            yerr=rv.eRV_j[sort_obs],
            color=cjb,
            marker="None",
            # mew=0.4,
            ls="",
            ecolor=cjb,
            elinewidth=0.8,
            capsize=0,
            zorder=zobs-1,
        )

    ax.set_xlim(minx, maxx)

    axs.append(ax)

    fig.align_ylabels(axs)

    if save_plot:
        print("full_path", cli.full_path)
        # folder_out = os.path.join(os.path.dirname(cli.full_path), 'plots')
        folder_out = os.path.join(cli.full_path, "plots")
        print("folder_out", folder_out)
        if not os.path.isdir(folder_out):
            os.makedirs(folder_out)
        fname = os.path.basename(sim_file)
        plt_file = os.path.join(folder_out, os.path.splitext(fname)[0])
        fig.savefig("{:s}.png".format(plt_file), bbox_inches="tight",
        )
        print("Saved plot into:")
        print("{:s}".format("{:s}.png".format(plt_file)))
        fig.savefig("{:s}.pdf".format(plt_file), bbox_inches="tight",
            # dpi=72
        )
        print("{:s}".format("{:s}.pdf".format(plt_file)))

    if show_plot:
        plt.show()

    # plt.close(fig)

    return fig

# ==============================================================================

# class CLI_RV:
#     def __init__(
#         self,
#         full_path=os.path.abspath("."),
#         idsim=1,
#         lmflag=0,
#         tscale=2440000.5,
#         labels="None",
#         samples_file=None,
#         limits="obs"
#     ):
#         self.full_path = os.path.abspath(full_path)
#         self.idsim=int(idsim)
#         self.lmflag=int(lmflag)
#         self.tscale=tscale
#         if str(labels).lower() != "none":
#             self.labels=labels.split()
#         else:
#             self.labels=None
#         self.samples_file=samples_file
#         self.limits=limits
        
#         return

# ==============================================================================


def get_args():

    parser = argparse.ArgumentParser()

    parser.add_argument(
        "-p",
        "--p",
        "-path",
        "--path",
        "-path-folder",
        "--path-folder",
        action="store",
        dest="full_path",
        required=True,
        help="Folder path",
    )
    parser.add_argument(
        "-s",
        "--s",
        "-sim-id",
        "--sim-id",
        action="store",
        dest="idsim",
        required=True,
        help="Simulation ID number",
    )
    parser.add_argument(
        "-lm",
        "--lm",
        "-lm-flag",
        "--lm-flag",
        action="store",
        dest="lmflag",
        default="0",
        help="LM flag: 0 = not used LM (default), 1 = used LM",
    )
    parser.add_argument(
        "-ts",
        "--ts" "-t-scale",
        "--t-scale",
        action="store",
        dest="tscale",
        default=None,
        help="Value to be subtract to the x-axis (time) in the plots. Default is None that means 2440000.5 will be subtract.",
    )
    parser.add_argument(
        "-labels",
        "--labels",
        action="store",
        dest="labels",
        default="None",
        help='Sequence of labels for RV data set, avoid spaces within 1 label, e.g.: "HARPS FIES". Put labels sorted as in the obsRV.dat file: "HARPS FIES" means 1 -> HARPS, 2 -> FIES. Default is None, meaning that it will plot "RV dataset #x" where x is an increasing number. ',
    )
    parser.add_argument(
        "-samples-file",
        "--samples-file",
        action="store",
        dest="samples_file",
        default="None",
        help="HDF5 file with T0 and RV from emcee samples to overplot on O-Cs.",
    )
    parser.add_argument(
        "-limits",
        "--limits",
        action="store",
        dest="limits",
        default="obs",
        help="Axis limits based on observations (obs) or synthetic samples (sam). Default obs.",
    )

    cli = parser.parse_args()

    cli.full_path = os.path.abspath(cli.full_path)

    cli.idsim = int(cli.idsim)

    cli.lmflag = int(cli.lmflag)

    if str(cli.tscale).lower() != "none":
        try:
            cli.tscale = float(cli.tscale)
        except:
            cli.tscale = 2440000.5
    else:
        cli.tscale = 2440000.5

    if str(cli.labels).lower() != "none":
        try:
            cli.labels = cli.labels.split()
        except:
            cli.labels = None
    else:
        cli.labels = None

    if str(cli.samples_file).lower() == "none":
        cli.samples_file = None
    else:
        if os.path.isfile(os.path.abspath(cli.samples_file)):
            cli.samples_file = os.path.abspath(cli.samples_file)
        else:
            cli.samples_file = None

    if cli.limits.lower()[:3] == "sam":
        cli.limits = "sam"
    else:
        cli.limits = "obs"

    return cli


# ==============================================================================


def main():

    cli = get_args()
    samples = read_samples(cli, samples_file=cli.samples_file)
    fig = plot_rv(cli, samples=samples, save_plot=True, show_plot=False)
    plt.close(fig)

    return


# ==============================================================================
# ==============================================================================

if __name__ == "__main__":
    main()
# ==============================================================================
