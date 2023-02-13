#!/usr/bin/env python
# -*- coding: utf-8 -*-

# no more "zero" integer division bugs!:P
import argparse
import os  #  os: operating system

# import time # time: execution time
import glob

# import sys # sys: system
import numpy as np  # array

import h5py

# import matplotlib as mpl
# mpl.use('Agg')
import matplotlib.pyplot as plt

# matplotlib rc params
# plt.rcParams['text.usetex'] = True
plt.rcParams["text.usetex"] = False
# plt.rcParams['font.family']       = 'sans-serif'
plt.rcParams["font.family"] = "serif"
plt.rcParams["font.serif"] = ["Computer Modern Roman", "Palatino", "DejaVu Serif"]
plt.rcParams["mathtext.fontset"] = "cm"
plt.rcParams["figure.figsize"] = [5, 5]
plt.rcParams["figure.facecolor"] = "white"
plt.rcParams["savefig.facecolor"] = "white"
plt.rcParams["figure.dpi"] = 200
plt.rcParams["savefig.dpi"] = 300
plt.rcParams["font.size"] = 12
plt.rcParams["xtick.labelsize"] = plt.rcParams["font.size"] - 2
plt.rcParams["ytick.labelsize"] = plt.rcParams["xtick.labelsize"]

# custom modules
import constants as cst
import ancillary as anc

CLI_OC = anc.CLI_OC

# ==============================================================================
class sim_data:
    def __init__(self):
        self.body_id = 0
        self.nTTs, self.ncols = 0, 0

        self.epo = []
        self.TTo = []
        self.eTTo = []
        self.TTs = []
        self.dTTos = []
        self.TTstat = []
        self.Teph = []
        self.Peph = []
        self.eTPeph = []
        self.TTlin = []
        self.oc_o = []
        self.oc_s = []

        self.T41o = []
        self.eT41o = []
        self.T41s = []
        self.dT41os = []
        self.T41stat = []

        self.period = []
        self.sma = []
        self.ecc = []
        self.inc = []
        self.meana = []
        self.argp = []
        self.truea = []
        self.longn = []

    def update_sim_data(self, body_id, sim_in, kep_ele=False):
        self.body_id = body_id
        self.nTTs, self.ncols = np.shape(sim_in)

        self.epo = sim_in[:, 0].astype(int)
        self.TTo = sim_in[:, 1]
        self.eTTo = sim_in[:, 2]
        self.TTs = sim_in[:, 3]
        self.dTTos = sim_in[:, 4]
        self.TTstat = sim_in[:, 5]

        self.epo, self.Teph, self.Peph, self.eTPeph = anc.compute_lin_ephem(
            self.TTo,
            eT0=self.eTTo,
            epoin=self.epo,
            # modefit='odr'
            modefit="wls",
        )
        self.TTlin = self.Teph + self.epo * self.Peph
        self.oc_o = self.TTo - self.TTlin
        self.oc_s = self.TTs - self.TTlin

        if kep_ele is False:
            if self.ncols == 11:  # T0 + T41
                self.T41o = sim_in[:, 6] / 1440.0  # duration saved in min
                self.eT41o = sim_in[:, 7] / 1440.0
                self.T41s = sim_in[:, 8] / 1440.0
                self.dT41os = sim_in[:, 9] / 1440.0
                self.T41stat = sim_in[:, 10]
        elif kep_ele is True:
            if self.ncols == 14:  # T0 + kep elements
                self.period = sim_in[:, 6]
                self.sma = sim_in[:, 7]
                self.ecc = sim_in[:, 8]
                self.inc = sim_in[:, 9]
                self.meana = sim_in[:, 10]
                self.argp = sim_in[:, 11]
                self.truea = sim_in[:, 12]
                self.longn = sim_in[:, 13]
            elif self.ncols == 19:  # T0 + T41 + kep elements
                self.T41o = sim_in[:, 6] / 1440.0  # duration saved in min
                self.eT41o = sim_in[:, 7] / 1440.0
                self.T41s = sim_in[:, 8] / 1440.0
                self.dT41os = sim_in[:, 9] / 1440.0
                self.T41stat = sim_in[:, 10]
                self.period = sim_in[:, 11]
                self.sma = sim_in[:, 12]
                self.ecc = sim_in[:, 13]
                self.inc = sim_in[:, 14]
                self.meana = sim_in[:, 15]
                self.argp = sim_in[:, 16]
                self.truea = sim_in[:, 17]
                self.longn = sim_in[:, 18]

        return

    def plot_kep_elem(self, cli, show_plot=False):

        nkep = 8

        lfont = plt.rcParams["font.size"]
        tfont = plt.rcParams["xtick.labelsize"]

        x = self.TTo - cli.tscale
        if cli.tscale > 0.0:
            xlabel = "BJD$_\mathrm{{TDB}} - {0:.3f}$".format(cli.tscale)
        else:
            xlabel = "BJD$_\mathrm{TDB}$"

        kep = np.column_stack(
            (
                self.period,
                self.sma,
                self.ecc,
                self.inc,
                self.meana,
                self.argp,
                self.truea,
                self.longn,
            )
        )
        kep_names = [
            "$P\ (d)$",
            "$a\ (\mathrm{{au}})$",
            "$e$",
            "$i\ (^\circ)$",
            "$\mathcal{{M}}\ (^\circ)$",
            "$\omega\ (^\circ)$",
            "$f\ (^\circ)$",
            "$\Omega\ (^\circ)$",
        ]

        fig, axs = plt.subplots(nkep, 1, sharex=True, figsize=(6, 6))

        for i_k, kname in enumerate(kep_names):
            axs[i_k].plot(
                x, kep[:, i_k], color="C0", marker="o", ms=3, mec="None", ls=""
            )
            axs[i_k].set_ylabel(kname, fontsize=lfont)
            axs[i_k].ticklabel_format(useOffset=False)
            axs[i_k].tick_params(direction="inout", labelsize=tfont)

        axs[-1].set_xlabel(xlabel, fontsize=lfont)
        axtitle(axs[0], labtitle="planet {:s} Keplerian Elements".format(self.planet))

        fig.align_ylabels(axs)

        file_in = self.file_in
        folder_out = os.path.join(os.path.dirname(file_in), "plots")
        if not os.path.isdir(folder_out):
            os.makedirs(folder_out)
        fname = os.path.basename(file_in)
        plt_file = os.path.join(
            folder_out, "{}_KepElements".format(os.path.splitext(fname)[0])
        )
        fig.savefig("%s.png" % (plt_file), bbox_inches="tight", dpi=200)
        print("Saved plot into:")
        print("%s" % ("%s.png" % (plt_file)))
        fig.savefig("%s.pdf" % (plt_file), bbox_inches="tight", dpi=72)
        print("%s" % ("%s.pdf" % (plt_file)))
        if show_plot:
            plt.show()
        plt.close(fig)

        return


# ==============================================================================


class oc_sample:
    def __init__(self, idnum, TTs_s, T41s_s):
        self.idplanet = int(idnum)
        self.TTs = np.array(TTs_s)
        self.nTTs = np.shape(self.TTs)[0]
        self.epo = np.arange(0, self.nTTs, 1)
        self.TTlin = np.zeros((self.nTTs))
        self.oc = np.zeros((self.nTTs))
        self.T41s = np.array(T41s_s) / 1440.0  # computed in min in trades

    def update(self, Teph, Peph):
        self.epo = anc.calculate_epoch(self.TTs, Teph, Peph)
        self.TTlin = Teph + self.epo * Peph
        self.oc = self.TTs - self.TTlin


# ==============================================================================


def get_sim_data(file_in, kep_ele=False):

    # epo 0 TTo 1 eTTo 2 TTs 3 dTTos 4 TTstat 5 T41o 6 eT41o 7 T41s 8 dT41os 9 T41stat 10
    sim_in = np.genfromtxt(file_in)
    body_id = int(file_in.split("NB")[1].split("_simT0")[0])

    sim = sim_data()
    sim.file_in = file_in
    sim.update_sim_data(body_id, sim_in, kep_ele=kep_ele)

    return sim


# ==============================================================================


def get_simT0_file_list(cli):

    # file_pattern = os.path.join(cli.full_path, '%d_%d_NB*_simT0.dat' \
    #                                             %(cli.idsim, cli.lmflag)
    #                             )
    file_dict = {}

    pttrn = "{:d}_{:d}_NB*_simT0.dat".format(cli.idsim, cli.lmflag)
    file_pattern = os.path.join(cli.full_path, pttrn)
    file_list = np.sort(glob.glob(file_pattern))
    if len(file_list) == 0:
        print("Cannot find in path {:s}".format(cli.full_path))
        print("files with pattern {:s}".format(pttrn))
    else:
        for f in file_list:
            id_body = int(f.split("NB")[1].split("_")[0])
            file_dict[id_body] = f
    return file_dict


# ==============================================================================


def read_samples(samples_file=None):

    if samples_file is not None:

        read_file = os.path.abspath(samples_file)
        print("Reading samples file {}".format(read_file))
        sh5 = h5py.File(read_file, "r")
        # n_samples = len(list(sh5.keys()))
        samples = []
        for igr in list(sh5.keys()):  # loop on groups, one for each sample
            TTfield = list(sh5[igr].keys())
            for kk in TTfield:
                if "TTs_" in kk:
                    idpl = int(kk.split("_")[1])
                    TTs = sh5["%s/%s" % (igr, kk)][...]
                    try:
                        T41s = sh5["%s/%s" % (igr, "T41s_%d" % (idpl))][...]
                    except:
                        T41s = np.zeros(np.shape(TTs))
                    samples.append(oc_sample(idpl, TTs, T41s))

        sh5.close()

    else:

        samples = None

    return samples


# ==============================================================================


def set_unit(cli, Aoc_d):

    u_in = cli.ocunit

    try:
        if u_in.lower() in "s sec second seconds".split():
            ocu = [86400.0, "s"]
        elif u_in.lower() in "m min minute minutes".split():
            ocu = [1440.0, "min"]
        elif u_in.lower() in "h hour hours".split():
            ocu = [24.0, "hours"]
        elif u_in.lower() == "auto":
            ocu = anc.set_automatic_unit_time(Aoc_d)
        else:
            ocu = [1.0, "d"]
    except:
        ocu = [1.0, "days"]

    return ocu


# ==============================================================================


def set_axis_default(
    ax, ticklabel_size=plt.rcParams["xtick.labelsize"], aspect="equal", labeldown=True
):

    ax.ticklabel_format(useOffset=False)
    ax.tick_params(direction="inout", labelsize=ticklabel_size, labelbottom=labeldown)
    ax.set_aspect(aspect)

    return


def set_symmetric_lim(ax, val):

    ax.set_xlim(-val, val)
    ax.set_ylim(-val, val)

    return


def axtitle(ax, labtitle="", fontsize=plt.rcParams["xtick.labelsize"]):

    ax.text(
        0.5,
        1.02,
        labtitle,
        horizontalalignment="center",
        # verticalalignment='center',
        verticalalignment="bottom",
        fontsize=fontsize,
        transform=ax.transAxes,
    )

    return


# ==============================================================================


def plot_oc_T41(cli, file_in, planet_name=None, samples=None, figsize=(5,5), save_plot=True, show_plot=False):

    sim = get_sim_data(file_in, kep_ele=cli.kep_ele)

    if cli.tscale is None:
        tscale = 2440000.5
    else:
        tscale = float(cli.tscale)
        # if tscale <= 0.0:
        #     tscale = 0.0

    Aoc_d = 0.5 * (np.max(sim.oc_o) - np.min(sim.oc_o))
    ocu = set_unit(cli, Aoc_d)
    sim.oc_unit = ocu

    letters = "a b c d e f g h i j k l m n o p q r s t u v w x y z".split()
    ibd = sim.body_id - 1
    if planet_name is None:
        planet = letters[ibd]
    else:
        planet = planet_name
    sim.planet = planet

    print("Planet {}".format(sim.planet))
    print("Observed  A_TTV = (MAX OC - MIN OC)/2 = {} {}".format(Aoc_d*ocu[0], ocu[1]))

    model_file = os.path.join(
        cli.full_path,
        "{:s}_models.hdf5".format(cli.sim_type)
    )
    oc_models = read_samples(samples_file=model_file)
    oc_model = [oo for oo in oc_models if oo.idplanet == sim.body_id][0]
    oc_model.update(sim.Teph, sim.Peph)

    Aoc_model = 0.5 * (np.max(oc_model.oc) - np.min(oc_model.oc)) * ocu[0]
    print("Simulated A_TTV = (MAX OC - MIN OC)/2 = {} {}".format(Aoc_model, ocu[1]))

    if samples is not None:
        samples_plt = [ss for ss in samples if ss.idplanet == sim.body_id]
        nsmp = len(samples_plt)
    else:
        nsmp = 0

    lfont = plt.rcParams["font.size"]
    tfont = plt.rcParams["xtick.labelsize"]
    dlfont = 4

    if tscale > 0.0:
        xlabel = "BJD$_\mathrm{{TDB}} - {0:.3f}$".format(tscale)
    else:
        xlabel = "BJD$_\mathrm{TDB}$"

    fig = plt.figure(figsize=(6, 6))
    fig.subplots_adjust(hspace=0.05, wspace=0.25)
    
    axs = []

    nrows = 3
    if sim.ncols in [11, 19]:
        ncols = 2
    else:
        ncols = 1

    malpha = 0.88
    # obs marker
    mso = 3
    mewo = 0.7
    ewo = 0.9
    cfo = "white"
    ceo = "C1"
    # sim marker
    mss = 3
    mews = 0.5
    cfs = "C0"
    ces = "white"
    # model line
    cfm = "black"
    # samples line
    cfsm = "gray"

    px = 0.05
    py = 0.05

    # =========================================
    # O-C
    ax = plt.subplot2grid((nrows, ncols), (0, 0), rowspan=2)
    set_axis_default(ax, ticklabel_size=tfont, aspect="auto", labeldown=False)
    ax.set_ylabel("O-C (%s)" % (ocu[1]), fontsize=lfont)
    axtitle(ax, labtitle="planet %s: transit times" % (planet), fontsize=lfont)

    ax.axhline(0.0, color="black", ls="-", lw=0.7)

    x = sim.TTo - tscale
    dx = np.max(x) - np.min(x)
    minx = np.min(x) - px * dx
    maxx = np.max(x) + px * dx
    # needed for obs limits
    xxx = np.concatenate((x, x, x))
    y = np.concatenate(
        (
            sim.oc_o * ocu[0] - sim.eTTo * ocu[0],
            sim.oc_o * ocu[0] + sim.eTTo * ocu[0],
            sim.oc_s * ocu[0],
        )
    )
    dy = np.max(y) - np.min(y)
    miny = np.min(y) - py * dy
    maxy = np.max(y) + py * dy

    # data
    lobs = ax.errorbar(
        x,
        sim.oc_o * ocu[0],
        yerr=sim.eTTo * ocu[0],
        color=ceo,
        ecolor=ceo,
        fmt="o",
        ms=mso,
        mfc=cfo,
        mec=ceo,
        mew=mewo,
        ls="",
        elinewidth=ewo,
        capsize=0,
        zorder=5,
        label="observations",
    )
    # trades
    (lsim,) = ax.plot(
        x,
        sim.oc_s * ocu[0],
        marker="o",
        color=cfs,
        ms=mss,
        mec=ces,
        mew=mews,
        ls="",
        zorder=6,
        alpha=malpha,
        label="simulations",
    )
    # model
    xx = oc_model.TTlin - tscale
    xxx = np.concatenate((xxx, xx))
    yy = np.concatenate((y, oc_model.oc * ocu[0]))
    (lmod, ) = ax.plot(
        xx,
        oc_model.oc * ocu[0],
        color=cfm,
        marker="None",
        ls="-",
        lw=0.6,
        zorder=5,
        alpha=1.0,
        label="model",
    )

    # samples
    if nsmp > 0:
        lsmp = []
        # yy = y
        for ss in samples_plt:
            ss.update(sim.Teph, sim.Peph)
            xx = ss.TTlin - tscale
            xxx = np.concatenate((xxx, xx))
            yy = np.concatenate((yy, ss.oc * ocu[0]))
            (ll,) = ax.plot(
                xx,
                ss.oc * ocu[0],
                color=cfsm,
                marker="None",
                ls="-",
                lw=0.4,
                zorder=4,
                alpha=0.4,
                label="samples",
            )
            lsmp.append(ll)

        ax.legend(handles=[lobs, lsim, lmod, lsmp[0]], loc="best", fontsize=lfont - dlfont)
    else:
        ax.legend(loc="best", fontsize=lfont - dlfont)

    if cli.limits == "sam":
        dx = np.max(xxx) - np.min(xxx)
        minx = np.min(xxx) - px * dx
        maxx = np.max(xxx) + px * dx
        dy = np.max(yy) - np.min(yy)
        miny = np.min(yy) - py * dy
        maxy = np.max(yy) + py * dy
    else:  # if obs & nsmp > 0 adjust y limits
        yyy = yy[np.logical_and(xxx >= minx, xxx <= maxx)]
        dy = np.max(yyy) - np.min(yyy)
        miny = np.min(yyy) - py * dy
        maxy = np.max(yyy) + py * dy

    ax.set_xlim(minx, maxx)
    # xlims = ax.get_xlim()
    ax.set_ylim(miny, maxy)

    axs.append(ax)

    # ==============================================
    # residuals TTo-TTs
    ax = plt.subplot2grid((nrows, ncols), (2, 0), rowspan=1)
    set_axis_default(ax, ticklabel_size=tfont, aspect="auto")
    ax.set_ylabel("res (%s)" % (ocu[1]), fontsize=lfont)
    ax.set_xlabel(xlabel, fontsize=lfont)

    ax.axhline(0.0, color="black", ls="-", lw=0.7)

    # data - trades
    ax.errorbar(
        x,
        sim.dTTos * ocu[0],
        yerr=sim.eTTo * ocu[0],
        color=ceo,
        ecolor=ceo,
        fmt="o",
        ms=mso,
        mfc=cfo,
        mec=ceo,
        mew=mewo,
        ls="",
        elinewidth=ewo,
        capsize=0,
        zorder=6,
    )
    rms_res = np.percentile(np.abs(sim.dTTos * ocu[0]), 68.27)
    print("rms Obs-Sim = 68.27th percentile of |residuals| = {} {}".format(rms_res, ocu[1]))

    ax.set_xlim(minx, maxx)

    axs.append(ax)

    if sim.ncols in [11, 19]:
        # ==============================================
        # T41 duration

        y = np.concatenate(
            (
                sim.T41o * ocu[0] - sim.eT41o * ocu[0],
                sim.T41o * ocu[0] + sim.eT41o * ocu[0],
                sim.T41s * ocu[0],
            )
        )
        dy = np.max(y) - np.min(y)
        miny = np.min(y) - py * dy
        maxy = np.max(y) + py * dy

        ax = plt.subplot2grid((nrows, ncols), (0, 1), rowspan=2)
        set_axis_default(ax, ticklabel_size=tfont, aspect="auto", labeldown=False)
        ax.set_ylabel("T$_{14}$ (%s)" % (ocu[1]), fontsize=lfont)
        axtitle(
            ax,
            labtitle="planet %s: durations as T$_4$-T$_1$" % (planet),
            fontsize=lfont,
        )

        # ax.axhline(np.median(data[:,6]), color='black', ls='-',lw=0.7)

        # data
        ax.errorbar(
            x,
            sim.T41o * ocu[0],
            yerr=sim.eT41o * ocu[0],
            color=ceo,
            ecolor=ceo,
            fmt="o",
            ms=mso,
            mfc=cfo,
            mec=ceo,
            mew=mewo,
            ls="",
            elinewidth=ewo,
            capsize=0,
            zorder=5,
            label="observations",
        )
        # trades
        ax.plot(
            x,
            sim.T41s * ocu[0],
            marker="o",
            color=cfs,
            ms=mss,
            mec=ces,
            mew=mews,
            ls="",
            zorder=6,
            alpha=malpha,
            label="simulations",
        )
        # model
        yy = np.concatenate((y, oc_model.T41s * ocu[0]))
        (lmod, ) = ax.plot(
            xx,
            oc_model.T41s * ocu[0],
            color=cfm,
            marker="None",
            ls="-",
            lw=0.6,
            zorder=5,
            alpha=1.0,
            label="model",
        )

        # samples
        if nsmp > 0:
            # yy = y
            for ss in samples_plt:
                # xx = ss.TTlin - tscale
                yy = np.concatenate((yy, ss.T41s * ocu[0]))
                # print np.min(ss.T41s*ocu[0]),np.max(ss.T41s*ocu[0])
                ax.plot(
                    xx,
                    ss.T41s * ocu[0],
                    color=cfsm,
                    marker="None",
                    ls="-",
                    lw=0.4,
                    zorder=4,
                    alpha=0.4,
                    label="samples",
                )

        if cli.limits == "sam":
            dx = np.max(xxx) - np.min(xxx)
            minx = np.min(xxx) - px * dx
            maxx = np.max(xxx) + px * dx
            dy = np.max(yy) - np.min(yy)
            miny = np.min(yy) - py * dy
            maxy = np.max(yy) + py * dy
        else:  # if obs & nsmp > 0 adjust y limits
            yyy = yy[np.logical_and(xxx >= minx, xxx <= maxx)]
            dy = np.max(yyy) - np.min(yyy)
            miny = np.min(yyy) - py * dy
            maxy = np.max(yyy) + py * dy

        ax.set_xlim(minx, maxx)
        # ax.set_xlim(xlims)
        ax.set_ylim(miny, maxy)

        axs.append(ax)

        # ==============================================
        # residuals T41o - T41s
        ax = plt.subplot2grid((nrows, ncols), (2, 1), rowspan=1)
        set_axis_default(ax, ticklabel_size=tfont, aspect="auto")
        ax.set_ylabel("res (min)", fontsize=lfont)
        ax.set_xlabel(xlabel, fontsize=lfont)

        ax.axhline(0.0, color="black", ls="-", lw=0.7)

        # data - trades
        ax.errorbar(
            x,
            sim.dT41os * ocu[0],
            yerr=sim.eT41o * ocu[0],
            color=ceo,
            ecolor=ceo,
            fmt="o",
            ms=mso,
            mfc=cfo,
            mec=ceo,
            mew=mewo,
            ls="",
            elinewidth=ewo,
            capsize=0,
            zorder=6,
        )

        ax.set_xlim(minx, maxx)
        # ax.set_xlim(xlims)

        axs.append(ax)

    fig.align_ylabels(axs)

    if save_plot == True:
        folder_out = os.path.join(os.path.dirname(file_in), "plots")
        if not os.path.isdir(folder_out):
            os.makedirs(folder_out)
        fname = os.path.basename(file_in)
        plt_file = os.path.join(folder_out, os.path.splitext(fname)[0])
        fig.savefig("%s.png" % (plt_file), bbox_inches="tight",
            # dpi=200
        )
        print("Saved plot into:")
        print("%s" % ("%s.png" % (plt_file)))
        fig.savefig("%s.pdf" % (plt_file), bbox_inches="tight",
            # dpi=72
        )
        print("%s" % ("%s.pdf" % (plt_file)))
    if show_plot:
        plt.show()
    # plt.close(fig)

    return fig


# ==============================================================================


# class CLI_OC:
#     def __init__(
#         self,
#         full_path=os.path.abspath("."),
#         idsim=1,
#         lmflag=0,
#         tscale=2440000.5,
#         ocunit="d",
#         samples_file=None,
#         limits="obs",
#         kep_ele=False
#     ):
#         self.full_path = os.path.abspath(full_path)
#         self.idsim = int(idsim)
#         self.lmflag = int(lmflag)
#         self.tscale = tscale
#         self.ocunit = ocunit
#         self.samples_file = samples_file
#         self.limits = limits
#         self.kep_ele = kep_ele

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
        "-u",
        "--u",
        "-unit",
        "--unit",
        action="store",
        dest="ocunit",
        default="d",
        help="Unit of the O-C/T41 plot. Possible (lower or upper case): d, day, days or h, hour, hours or m, min, minute, minutes or s, sec, second, seconds. Default is d.",
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
    parser.add_argument(
        "-kep",
        "--kep",
        "-kep-ele",
        "--kep-ele",
        action="store",
        dest="kep_ele",
        default=False,
        help="Keplerian elements are present in the simulated transit file? Default False",
    )

    cli = parser.parse_args()

    cli.full_path = os.path.abspath(cli.full_path)

    cli.idsim = int(cli.idsim)

    cli.lmflag = int(cli.lmflag)

    if cli.samples_file.lower() == "none":
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

    if str(cli.kep_ele).lower() in ["t", "tr", "tru", "true", "1"]:
        cli.kep_ele = True
    else:
        cli.kep_ele = False

    return cli


# ==============================================================================


def main():

    cli = get_args()
    samples = read_samples(samples_file=cli.samples_file)
    file_list = get_simT0_file_list(cli)
    # sims = []
    for idbody, f in file_list.items():
        pl_name = "{}".format(idbody)
        fig = plot_oc_T41(cli, f, planet_name =pl_name, samples=samples, save_plot=True, show_plot=False)
        plt.close(fig)

    return


# ==============================================================================
# ==============================================================================

if __name__ == "__main__":
    main()
# ==============================================================================
