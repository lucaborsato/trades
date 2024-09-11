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

# custom modules
import constants as cst
import ancillary as anc
import pytrades

CLI_OC = anc.CLI_OC
filled_markers = anc.filled_markers
size_markers = anc.size_markers

anc.set_rcParams()

# ==============================================================================
def set_observation_sources(sources_id, idsource_name=None):

    u_id = np.unique(sources_id)
    # print("idsource_name = {}".format(idsource_name))
    if len(u_id) == 1:
        idname = {1: "observations"}
    else:
        idname = {i: "obs.{:d}".format(i) for i in u_id}
        colors = ["C1"]*len(idname)
    
    # print("idname = {}".format(idname))
    if idsource_name is not None:
        for kin, vin in idsource_name.items():
            if kin in u_id:
                idname[kin] = vin
    # print("idname = {}".format(idname))

    # sources = np.array([""] * len(sources_id))
    sources = [""] * len(sources_id)
    # for kid, vid in idname.items():
    #     sources[sources_id == kid] = vid
    for i_id, s_id in enumerate(sources_id):
        sources[i_id] = idname[s_id]
    sources = np.array(sources)

    # print("sources = {}".format(sources))
    return sources, idname

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
        self.sources_id = []
        self.sources = []
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

    def update_sim_data(self, body_id, sim_in, idsource_name=None, kep_ele=False):
        self.body_id = body_id
        self.sim_in = sim_in
        self.nTTs, self.ncols_all = np.shape(sim_in)
        # remove source id column
        self.ncols = self.ncols_all - 1

        self.epo = sim_in[:, 0].astype(int)
        self.TTo = sim_in[:, 1]
        self.eTTo = sim_in[:, 2]
        self.TTs = sim_in[:, 3]
        self.dTTos = sim_in[:, 4]
        self.TTstat = sim_in[:, 5]
        self.sources_id = sim_in[:, -1].astype(int)
        self.sources, self.idname = set_observation_sources(self.sources_id, idsource_name=idsource_name)

        Tref, Pref, chi2 = pytrades.linear_fit(self.epo, self.TTo, ey=self.eTTo)
        self.Teph = Tref[0]
        self.Peph = Pref[0]
        self.eTPeph = [Tref[1], Pref[1]]
        self.chi2 = chi2

        self.TTlin = self.Teph + self.epo * self.Peph
        self.oc_o = self.TTo - self.TTlin
        self.oc_s = self.TTs - self.TTlin

        ncols_tra = 6
        ncols_kep = 8
        ncols_dur = 5

        self.ncols_tra = ncols_tra
        self.ncols_kep = ncols_kep
        self.ncols_dur = ncols_dur

        print("Kep_elem {} and ncols = {}".format(kep_ele, self.ncols))

        # only transits
        if self.ncols == ncols_tra:
            print("No Keplerian elements or Transit Durations in the file.")
        elif self.ncols == (ncols_tra + ncols_kep):
            if kep_ele:
                self.period = sim_in[:, 6]
                self.sma = sim_in[:, 7]
                self.ecc = sim_in[:, 8]
                self.inc = sim_in[:, 9]
                self.meana = sim_in[:, 10]
                self.argp = sim_in[:, 11]
                self.truea = sim_in[:, 12]
                self.longn = sim_in[:, 13]
        elif self.ncols == (ncols_tra + ncols_dur):
            if kep_ele:
                print("No Keplerian elements in the file, but found Transit durations")
            else:
                print("Found Transit durations in the file")
            self.T41o = sim_in[:, -6]  # duration saved in min
            self.eT41o = sim_in[:, -5]
            self.T41s = sim_in[:, -4]
            self.dT41os = sim_in[:, -3]
            self.T41stat = sim_in[:, -2]
        elif self.ncols == (ncols_tra + ncols_kep + ncols_dur):
            print("Found Keplerian elements and Transit durations in the file")
            if kep_ele:
                self.period = sim_in[:, 6]
                self.sma = sim_in[:, 7]
                self.ecc = sim_in[:, 8]
                self.inc = sim_in[:, 9]
                self.meana = sim_in[:, 10]
                self.argp = sim_in[:, 11]
                self.truea = sim_in[:, 12]
                self.longn = sim_in[:, 13]
                self.T41o = sim_in[:, -6]  # duration saved in min
                self.eT41o = sim_in[:, -5]
                self.T41s = sim_in[:, -4]
                self.dT41os = sim_in[:, -3]
                self.T41stat = sim_in[:, -2]
            else:
                self.T41o = sim_in[:, -6] # duration saved in min
                self.eT41o = sim_in[:, -5]
                self.T41s = sim_in[:, -4]
                self.dT41os = sim_in[:, -3]
                self.T41stat = sim_in[:, -2]

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
            "$P\ (\mathrm{{d}})$",
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
                x, kep[:, i_k], color="C0", marker="o", ms=2, mec="None", ls=""
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
        fig.savefig("%s.png" % (plt_file), bbox_inches="tight", dpi=300)
        print("Saved plot into:")
        print("%s" % ("%s.png" % (plt_file)))
        fig.savefig("%s.pdf" % (plt_file), bbox_inches="tight", dpi=300)
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
        self.T41s = np.array(T41s_s)  # computed in min in trades

    def update(self, Teph, Peph):
        self.epo = anc.calculate_epoch(self.TTs, Teph, Peph)
        self.TTlin = Teph + self.epo * Peph
        self.oc = self.TTs - self.TTlin

    def select_sub_sample(self, sel):
        self.TTs = self.TTs[sel]
        self.epo = self.epo[sel]
        self.TTlin = self.TTlin[sel]
        self.oc = self.oc[sel]
        self.T41s = self.T41s[sel]
        self.nTTs = len(self.TTs)


# ==============================================================================


def get_sim_data(file_in, idsource_name=None, kep_ele=False):
    # epo 0 TTo 1 eTTo 2 TTs 3 dTTos 4 TTstat 5 T41o 6 eT41o 7 T41s 8 dT41os 9 T41stat 10
    # or
    # epoT0obs 0 T0obs 1 eT0obs 2 T0_sim 3 T0obs-T0_sim 4 T0_stat 5 
    # period_d 6 sma_au 7 ecc 8 inc_deg 9 meana_deg 10 argp_deg 11 truea_deg 12 longn_deg 13 
    # Dur_obs 14 eDur_obs 15 Dur_sim 16 Dur_obs-Dur_sim 17 Dur_stat 18
    sim_in = np.genfromtxt(file_in)
    body_id = int(file_in.split("NB")[1].split("_simT0")[0])

    sim = sim_data()
    sim.file_in = file_in
    sim.update_sim_data(body_id, sim_in, idsource_name=idsource_name, kep_ele=kep_ele)

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


def set_unit_base(u_in, Aoc_d):
    try:
        if u_in.lower() in "s sec second seconds".split():
            ocu = [86400.0, "sec"]
        elif u_in.lower() in "m min minute minutes".split():
            ocu = [1440.0, "min"]
        elif u_in.lower() in "h hour hours".split():
            ocu = [24.0, "hours"]
        elif u_in.lower() == "auto":
            ocu = anc.set_automatic_unit_time(Aoc_d)
        else:
            ocu = [1.0, "days"]
    except:
        ocu = [1.0, "days"]

    return ocu


def set_unit(cli, Aoc_d):
    u_in = cli.ocunit

    ocu = set_unit_base(u_in, Aoc_d)

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


def plot_oc_T41(
    cli,
    file_in,
    planet_name=None,
    samples=None,
    figsize=(5, 5),
    save_plot=True,
    show_plot=False,
):
    sim = get_sim_data(file_in, idsource_name=cli.idsource_name, kep_ele=cli.kep_ele)
    u_src = np.unique(sim.sources)
    # n_src = len(u_src)
    # ocolors = anc.set_colors(n_src, colormap=cli.color_map)
    n_idsrc = np.max(list(cli.idsource_name.keys()))
    
    if isinstance(cli.color_map, dict):
        # full_colors = cli.color_map
        ocolors = [""]*len(cli.color_map)
        for k, v in cli.color_map.items():
            ocolors[k-1] = v
    elif isinstance(cli.color_map, str):
        if cli.color_map in anc.all_color_maps:
            colormap = cli.color_map
        else:
            colormap = "nipy_spectral"
        full_colors = anc.set_colors(n_idsrc, colormap=colormap)
        ocolors = [full_colors[k-1] for k in cli.idsource_name.keys()]

    if cli.tscale is None:
        tscale = 2440000.5
    else:
        tscale = float(cli.tscale)
        # if tscale <= 0.0:
        #     tscale = 0.0

    Aoc_d = 0.5 * (np.max(sim.oc_o) - np.min(sim.oc_o))
    ocu = set_unit(cli, Aoc_d)
    sim.oc_unit = ocu

    # letters = "a b c d e f g h i j k l m n o p q r s t u v w x y z".split()
    letters = anc.letters
    ibd = sim.body_id - 1
    if planet_name is None:
        planet = letters[ibd]
    else:
        planet = planet_name
    sim.planet = planet

    print("\nPlanet {}".format(sim.planet))
    print("Read {}".format(file_in))
    # print("sim.sources_id = {}".format(sim.sources_id))
    # print("sim.sources    = {}".format(sim.sources))
    print("idsource_name = {}".format(cli.idsource_name))
    print("OBS id -> name = {}".format(sim.idname))
    print(
        "Observed  A_TTV = (MAX OC - MIN OC)/2 = {} {}".format(Aoc_d * ocu[0], ocu[1])
    )
    print("Tref = {} +/- {}".format(sim.Teph, sim.eTPeph[0]))
    print("Pref = {} +/- {}".format(sim.Peph, sim.eTPeph[1]))

    model_file = os.path.join(cli.full_path, "{:s}_models.hdf5".format(cli.sim_type))
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
    dlfont = 2

    if tscale > 0.0:
        xlabel = "BJD$_\mathrm{{TDB}} - {0:.3f}$".format(tscale)
    else:
        xlabel = "BJD$_\mathrm{TDB}$"

    fig = plt.figure(figsize=figsize)
    fig.subplots_adjust(hspace=0.05, wspace=0.25)

    axs = []

    nrows = 3
    # if sim.ncols in [11, 19]:
    #     ncols = 2
    # else:
    #     ncols = 1
    if (sim.ncols == sim.ncols_tra or sim.ncols == (sim.ncols_tra+sim.ncols_kep)):
        ncols = 1
    else:
        ncols = 2

    malpha = 0.88
    
    # obs marker
    # mso = 10
    mewo = 0.6
    ewo = 1.1
    # cfo = "white"
    # ceo = "C1"
    zobs = 8
    ceo = "black"

    # sim marker
    # mss = 2.5
    # mss = mso * 0.35
    
    mss = anc.size_default * 0.5
    mews = 0.7
    cfs = "black"
    ces = "black"
    zsim = 9
    
    # model line
    cfm = "black" #"dimgray"
    zmod = 7
    
    # samples line
    # cfsm = "silver"
    cfsm = plt.get_cmap("gray")
    gval = 0.6
    dg = 0.1
    zsmp = 6

    zleg = zsim + 1

    px = 0.05
    py = 0.05

    leg_ncol = 1
    if "out" in cli.legend:
        leg_loc = "center left"
        leg_bbox = (1.0, 0.5)
    else:
        leg_loc = "best"
        leg_bbox = (0.025, 0.025, 0.925, 0.925)

    # =========================================
    # O-C
    ax = plt.subplot2grid((nrows, ncols), (0, 0), rowspan=2)
    set_axis_default(ax, ticklabel_size=tfont, aspect="auto", labeldown=False)
    ax.set_ylabel("O-C ({:s})".format(ocu[1]), fontsize=lfont)
    axtitle(ax, labtitle="planet {:s}: transit times".format(planet), fontsize=lfont)

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
    lobs = []
    # for isrc, src in enumerate(u_src):
    for ksrc, src in cli.idsource_name.items():
        isrc = ksrc - 1
        sel = src == sim.sources
        if np.sum(sel) > 0:
            eco = ocolors[isrc]
            lobs.append(ax.errorbar(
                x[sel],
                sim.oc_o[sel] * ocu[0],
                yerr=sim.eTTo[sel] * ocu[0],
                color=ocolors[isrc],
                ecolor=eco,
                fmt=filled_markers[isrc],
                ms=size_markers[isrc],
                mfc=ocolors[isrc],
                mec=ceo,
                mew=mewo,
                ls="",
                elinewidth=ewo,
                capsize=0,
                alpha=1.0,
                zorder=zobs,
                label=src,
            ))
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
        zorder=zsim,
        alpha=malpha,
        label="simulations",
    )
    # model
    xx = oc_model.TTlin - tscale
    # xx = oc_model.TTs - tscale
    xxx = np.concatenate((xxx, xx))
    yy = np.concatenate((y, oc_model.oc * ocu[0]))
    (lmod,) = ax.plot(
        xx,
        oc_model.oc * ocu[0],
        color=cfm,
        marker="None",
        ls="-",
        lw=0.6,
        zorder=zmod,
        alpha=1.0,
        label="model",
    )

    # samples
    if nsmp > 0:
        print("samples plot ... ")
        # # ==== plot all the samples o-c
        # lsmp = []
        # for ss in samples_plt:
        #     ss.update(sim.Teph, sim.Peph)
        #     xx = ss.TTlin - tscale
        #     xxx = np.concatenate((xxx, xx))
        #     yy = np.concatenate((yy, ss.oc * ocu[0]))
        #     (ll,) = ax.plot(
        #         xx,
        #         ss.oc * ocu[0],
        #         color=cfsm(gval),
        #         marker="None",
        #         ls="-",
        #         lw=0.4,
        #         zorder=zsmp,
        #         alpha=0.4,
        #         label="samples",
        #     )
        #     lsmp.append(ll)

        # # ==== plot 1/2/3 CI        
        n_sim = len(xx)
        print("n_sim (TTs) = {:d}".format(n_sim))
        oc_smp = []
        print("creating oc_smp list ...")
        for ss in samples_plt:
            ss.update(sim.Teph, sim.Peph)
            oc_smp.append(ss.oc)
        print("oc_smp list to array ... ")
        oc_smp = np.array(oc_smp).T *ocu[0]
        print("shape(oc_smp) = ", np.shape(oc_smp))
        c1, c2, c3 = 0.6827, 0.9544, 0.9974
        hc1, hc2, hc3 = c1*0.5, c2*0.5, c3*0.5

        # TESTING HDI
        # print("computing HDI @ 68.27-95.44-99.73 % ... ")
        # # hdi1, hdi2, hdi3 = [], [], []
        # hdi1, hdi2, hdi3 = np.zeros((n_sim, 2)), np.zeros((n_sim, 2)), np.zeros((n_sim, 2))
        # for i in range(n_sim):
        #     ocy = oc_smp[i, :]
        #     hhh = anc.hpd(ocy, cred=0.6827)
        #     hdi1[i, :] = hhh
        #     hhh = anc.hpd(ocy, cred=0.9544)
        #     hdi2[i, :] = hhh
        #     hhh = anc.hpd(ocy, cred=0.9974)
        #     hdi3[i, :] = hhh
        # TODO: TEST PERCENTILE
        # print("computing CI @ 68.27-95.44-99.73 % ... ")
        hdi1 = np.percentile(oc_smp, [50 - (100*hc1), 50 + (100*hc1)], axis=1).T
        hdi2 = np.percentile(oc_smp, [50 - (100*hc2), 50 + (100*hc2)], axis=1).T
        hdi3 = np.percentile(oc_smp, [50 - (100*hc3), 50 + (100*hc3)], axis=1).T

        print("shape(hdi1) = ", np.shape(hdi1))
        print("shape(hdi2) = ", np.shape(hdi2))
        print("shape(hdi3) = ", np.shape(hdi3))
        
        print("plot HDI ... ")
        yy = np.concatenate((yy, hdi3[:, 0], hdi3[:, 1]))
        xxx = np.concatenate((xxx, xx, xx))

        ax.fill_between(
            xx,
            hdi1[:, 0],
            hdi1[:, 1],
            color=cfsm(gval),
            alpha=1.0,
            lw=0.0,
            zorder=zsmp,
        )
        ax.fill_between(
            xx,
            hdi2[:, 0],
            hdi2[:, 1],
            color=cfsm(gval+dg),
            alpha=1.0,
            lw=0.0,
            zorder=zsmp-1,
        )
        ax.fill_between(
            xx,
            hdi3[:, 0],
            hdi3[:, 1],
            color=cfsm(gval+(dg*2)),
            alpha=1.0,
            lw=0.0,
            zorder=zsmp-2,
        )
    
    # =================================================================
    ax.legend(
        handles= lobs + [lsim], #, lmod, lsmp[0]],
        loc=leg_loc, 
        bbox_to_anchor=leg_bbox,
        ncols=leg_ncol,
        fontsize=lfont - dlfont,
    ).set_zorder(zleg)
    # else:
    #     ax.legend(
    #         loc=leg_loc, 
    #         bbox_to_anchor=leg_bbox,
    #         ncols=leg_ncol,
    #         fontsize=lfont - dlfont
    #     )

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
    # for isrc, src in enumerate(u_src):
    for ksrc, src in cli.idsource_name.items():
        isrc = ksrc - 1
        sel = src == sim.sources
        if np.sum(sel) > 0:
            eco = ocolors[isrc]
            ax.errorbar(
                x[sel],
                sim.dTTos[sel] * ocu[0],
                yerr=sim.eTTo[sel] * ocu[0],
                color=ocolors[isrc],
                ecolor=eco,
                fmt=filled_markers[isrc],
                ms=size_markers[isrc],
                mfc=ocolors[isrc],
                mec=ceo,
                mew=mewo,
                ls="",
                elinewidth=ewo,
                capsize=0,
                zorder=zobs,
            )
            rms_res = np.percentile(np.abs(sim.dTTos[sel] * ocu[0]), 68.27)
            print(
                "rms Obs-Sim {} = 68.27th percentile of |residuals| = {} {}".format(
                    src, rms_res, ocu[1]
                )
            )
    rms_res = np.percentile(np.abs(sim.dTTos * ocu[0]), 68.27)
    print(
        "rms Obs-Sim FULL = 68.27th percentile of |residuals| = {} {}".format(
            rms_res, ocu[1]
        )
    )
    wres = sim.dTTos/sim.eTTo
    chi2 = np.sum(wres*wres)
    print("T0s contribution to the chi_square = {:.3f}".format(chi2))

    ax.set_xlim(minx, maxx)

    axs.append(ax)

    n_cols = ncols == 2
    dur_fit = sim.T41o is not None and len(sim.T41o) > 0
    if n_cols and dur_fit:
        # ==============================================
        # T41 duration
        y = np.concatenate(
            (
                sim.T41o - sim.eT41o,
                sim.T41o + sim.eT41o,
                sim.T41s,
            )
        )
        dy = np.max(y) - np.min(y)
        miny = np.min(y) - py * dy
        maxy = np.max(y) + py * dy

        ax = plt.subplot2grid((nrows, ncols), (0, 1), rowspan=2)
        set_axis_default(ax, ticklabel_size=tfont, aspect="auto", labeldown=False)
        ax.set_ylabel("T$_{{14}}$ (min)", fontsize=lfont)
        axtitle(
            ax,
            labtitle="planet {:s}: durations as T$_4$-T$_1$".format(planet),
            fontsize=lfont,
        )

        # ax.axhline(np.median(data[:,6]), color='black', ls='-',lw=0.7)

        # data
        for ksrc, src in cli.idsource_name.items():
            isrc = ksrc - 1
            sel = src == sim.sources
            if np.sum(sel) > 0:
                eco = ocolors[isrc]
                ax.errorbar(
                    x[sel],
                    sim.T41o[sel],
                    yerr=sim.eT41o[sel],
                    color=eco,
                    ecolor=eco,
                    fmt="o",
                    ms=size_markers[isrc],
                    mfc=cfo,
                    mec=ceo,
                    mew=mewo,
                    ls="",
                    elinewidth=ewo,
                    capsize=0,
                    alpha=1.0,
                    zorder=zobs,
                    label=src,
                )

        # trades
        ax.plot(
            x,
            sim.T41s,
            marker="o",
            color=cfs,
            ms=mss,
            mec=ces,
            mew=mews,
            ls="",
            zorder=zsim,
            alpha=malpha,
            label="simulations",
        )
        # model
        yy = np.concatenate((y, oc_model.T41s))
        (lmod,) = ax.plot(
            xx,
            oc_model.T41s,
            color=cfm,
            marker="None",
            ls="-",
            lw=0.6,
            zorder=zmod,
            alpha=1.0,
            label="model",
        )

        # samples
        if nsmp > 0:
            # yy = y
            for ss in samples_plt:
                # xx = ss.TTlin - tscale
                yy = np.concatenate((yy, ss.T41s))
                # print np.min(ss.T41s*ocu[0]),np.max(ss.T41s*ocu[0])
                ax.plot(
                    xx,
                    ss.T41s,
                    color=cfsm,
                    marker="None",
                    ls="-",
                    lw=0.4,
                    zorder=zsmp,
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
        for ksrc, src in cli.idsource_name.items():
            isrc = ksrc - 1
            sel = src == sim.sources
            if np.sum(sel) > 0:
                eco = ocolors[isrc]
                ax.errorbar(
                    x[sel],
                    sim.dT41os[sel],
                    yerr=sim.eT41o[sel],
                    color=eco,
                    ecolor=eco,
                    fmt="o",
                    ms=size_markers[isrc],
                    mfc=cfo,
                    mec=ceo,
                    mew=mewo,
                    ls="",
                    elinewidth=ewo,
                    capsize=0,
                    zorder=zobs,
                )

        ax.set_xlim(minx, maxx)
        # ax.set_xlim(xlims)
        wres = sim.dT41os/sim.eT41o
        chi2 = np.sum(wres*wres)
        print("T14 contribution to the chi_square = {:.3f}".format(chi2))

        axs.append(ax)

    fig.align_ylabels(axs)
    plt.tight_layout()

    if save_plot == True:
        folder_out = os.path.join(os.path.dirname(file_in), "plots")
        if not os.path.isdir(folder_out):
            os.makedirs(folder_out)
        fname = os.path.basename(file_in)
        plt_file = os.path.join(folder_out, os.path.splitext(fname)[0])
        for ext in ["png", "pdf"]:
            fig.savefig(
                "{:s}.{:s}".format(plt_file, ext),
                bbox_inches="tight",
                dpi=300
            )   
            print("Saved plot into:")
            print("{:s}.{:s}".format(plt_file, ext))

    if show_plot:
        plt.show()
    # plt.close(fig)

    if cli.kep_ele:
        sim.plot_kep_elem(cli, show_plot=show_plot)
        
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
        fig = plot_oc_T41(
            cli,
            f,
            planet_name=pl_name,
            samples=samples,
            save_plot=True,
            show_plot=False,
        )
        plt.close(fig)

    return


# ==============================================================================
# ==============================================================================

if __name__ == "__main__":
    main()
# ==============================================================================
