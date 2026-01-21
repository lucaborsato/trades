#!/usr/bin/env python
# -*- coding: utf-8 -*-

# ==============================================================================
# IMPORT

# no more "zero" integer division bugs!:P
import argparse
import numpy as np  # array
import os

# import sys
# import glob
# import astropy.time as atime
# ==============================================================================
import matplotlib as mpl

mpl.use("Agg", warn=False)
# mpluse("Qt4Agg")
# mpluse("TkAgg")
import matplotlib.pyplot as plt

plt.rcParams["text.usetex"] = True
plt.rcParams["font.family"] = "serif"
# ==============================================================================
import gls

# ==============================================================================
# custom modules
# script_path = os.path.realpath(__file__)
# module_path = os.path.abspath(os.path.join(os.path.dirname(script_path),
#                                            '../pytrades'))
# if(module_path not in sys.path):
#   print(module_path)
#   sys.path.append(module_path)

import ancillary as anc
import constants as cst

# from pytrades_lib import f90trades
# ==============================================================================


def get_args():

    parser = argparse.ArgumentParser(description="CHEOPS TTV LEVEL ANALYSIS")

    parser.add_argument(
        "-f",
        "--f",
        "-folder-path",
        "--folder-path",
        action="store",
        dest="folder_path",
        required=True,
        help="The path (absolute or relative) with simulation files for TRADES.",
    )

    parser.add_argument(
        "-s",
        "--s",
        "-sim-id",
        "--sim-id",
        action="store",
        dest="sim_id",
        default="1",
        help="Number ID of the simulation. Default is 1.",
    )

    parser.add_argument(
        "-b",
        "--b",
        "-body-id",
        "--body-id",
        action="store",
        dest="body_id",
        default="2",
        help="Number ID of the body: star is 1, planets from 2 to NB. Default is 2.",
    )

    parser.add_argument(
        "-ts",
        "--ts",
        "-time-start",
        "--time-start",
        action="store",
        dest="t_start",
        default=0.0,
        help="Enter starting time of selection for the plot in same unit of the simulation. Default is 0.0",
    )

    parser.add_argument(
        "-te",
        "--te",
        "-time-end",
        "--time-end",
        action="store",
        dest="t_end",
        default=365.25,
        help="Enter ending time of selection for the plot in same unit of the simulation. Default is 365.25",
    )

    parser.add_argument(
        "-t",
        "--t",
        "-t-ref",
        "--t-ref",
        action="store",
        dest="tref",
        default=None,
        help="Transit time of reference for/from linear ephemeris. Default is None, and it will be computed from a linear ephemeris with guess the first transit time.",
    )

    parser.add_argument(
        "-p",
        "--p",
        "-p-ref",
        "--p-ref",
        action="store",
        dest="pref",
        default=None,
        help="Period (in days) of reference for/from linear ephemeris. Default is None, and it will computed from the median of the differences between each transit time.",
    )

    parser.add_argument(
        "-n",
        "--n",
        "-name",
        "--name",
        action="store",
        dest="name",
        default="planet X",
        help="Name of simulation appearing in the title plot",
    )

    parser.add_argument(
        "-u",
        "--u",
        "-unit",
        "--unit",
        action="store",
        dest="ocunit",
        default="min",
        help="Unit for O-C plot: sec, min, or hour. Default min.",
    )

    cli = parser.parse_args()

    cli.folder_path = os.path.abspath(cli.folder_path)

    try:
        cli.sim_id = int(cli.sim_id)
        if cli.sim_id <= 0:
            cli.sim_id = 1
    except:
        cli.sim_id = 1

    try:
        cli.body_id = int(cli.body_id)
        if cli.body_id <= 1:
            cli.body_id = 2
    except:
        cli.body_id = 2

    try:
        cli.t_start = float(cli.t_start)
    except:
        cli.t_start = 0.0

    try:
        cli.t_end = float(cli.t_end)
    except:
        cli.t_end = 365.25

    try:
        cli.tref = float(cli.tref)
    except:
        cli.tref = None

    try:
        cli.pref = float(cli.pref)
    except:
        cli.pref = None

    return cli


# ==============================================================================


def set_unit(cli):

    u_in = cli.ocunit

    try:
        if u_in.lower() in "s sec second seconds".split():
            ocu = [86400.0, "sec"]
        elif u_in.lower() in "m min minute minutes".split():
            ocu = [1440.0, "min"]
        elif u_in.lower() in "h hour hours".split():
            ocu = [24.0, "hours"]
        else:
            ocu = [1.0, "days"]
    except:
        ocu = [1.0, "days"]

    return ocu


# ==============================================================================


def printlog(line, log=None):

    print("%s" % (line))
    if log is not None:
        log.write("%s\n" % (line))

    return


# ==============================================================================


def read_file(full_file):

    ttra_lte = np.genfromtxt(full_file, usecols=(0, 1))
    ttra = ttra_lte[:, 0] + ttra_lte[:, 1]
    nttra = np.shape(ttra)[0]
    if nttra == 0:
        ttra = None

    return ttra, nttra


# ==============================================================================


class Transit:
    def __init__(self, full_file, sim_id, body_id, tref, pref, tstart, tend):
        self.full_file = full_file
        self.dirname = os.path.dirname(full_file)
        self.filename = os.path.basename(full_file)
        self.sim_id = sim_id
        self.body_id = body_id
        self.ttra, self.nttra = read_file(self.full_file)
        self.tref = tref
        self.pref = pref
        self.ttra_sel = self.ttra[
            np.logical_and(self.ttra >= tstart, self.ttra <= tend)
        ]
        self.nttra_sel = np.shape(self.ttra_sel)[0]


# ==============================================================================


def main():

    # get folder name (and sim id?)
    cli = get_args()

    # determines number of files with transit times
    full_file = os.path.join(
        cli.folder_path, "%d_0_NB%d_tra.dat" % (cli.sim_id, cli.body_id)
    )
    target = os.path.basename(cli.folder_path)
    target = target.replace("_", "\\_")

    olog = open(
        os.path.join(
            os.path.dirname(full_file), "log_oc_%d_NB%d.txt" % (cli.sim_id, cli.body_id)
        ),
        "w",
    )

    printlog("\n%s\nREAD FILE: %s" % (target, full_file), olog)

    ocu = set_unit(cli)

    # cheops launch date
    # t_launch = atime.Time('2018-06-30T06:00:00.0', format='isot', scale='tdb')
    # t_launch = atime.Time('2019-01-01T07:00:00.0', format='isot', scale='tdb')
    # t_launch = t_launch.jd

    # cheops end nominal mission
    # t_end = t_launch + 365.25*3.5
    # t_end = t_launch + 365.25*4. # it could start in the first 6 months of 2019
    # cheops_dur = t_end - t_launch
    # printlog('TIME SELECTION BJD: %.6f - %.6f' %(t_launch, t_end), olog)

    # t_launch = 0.0
    # t_end = 365.25*3.5
    t_launch = cli.t_start
    t_end = cli.t_end
    cheops_dur = t_end - t_launch
    printlog("TIME SELECTION BJD: %.6f - %.6f" % (t_launch, t_end), olog)

    # find transit times in file
    transits = Transit(
        full_file, cli.sim_id, cli.body_id, cli.tref, cli.pref, t_launch, t_end
    )

    Tin = transits.ttra_sel[transits.nttra_sel // 2]
    Pin = np.median(np.diff(transits.ttra_sel))
    epo_sel = np.rint( (transits.ttra_sel - Tin) / Pin)
    Tref, Pref, _ = pytrades.linear_fit(epo_sel, transits.ttra_sel)
    transits.tref = Tref[0]
    transits.pref = Pref[0]

    printlog("LINEAR EPHEMERIS: %.6f + N x %.6f" % (transits.tref, transits.pref), olog)
    printlog(
        "FOUND %d TRANSIT DURING CHEOPS ON A TOTAL OF %d"
        % (transits.nttra_sel, transits.nttra),
        olog,
    )

    tlin = transits.tref + epo_sel * transits.pref
    oc_d = transits.ttra_sel - tlin
    eoc_d = np.zeros((transits.nttra_sel)) + 1.0 / 86400.0

    tt = transits.ttra_sel - transits.tref
    gls_o = gls.Gls(
        (tt, oc_d, eoc_d),
        Pbeg=3.0 * transits.pref,
        Pend=2.0 * cheops_dur,
        ofac=10,
        verbose=False,
    )
    gls_o.info()
    amp_ttv = gls_o.hpstat["amp"]
    p_ttv = 1.0 / gls_o.hpstat["fbest"]
    printlog("amp_ttv = %.6f d P_ttv = %.6f d" % (amp_ttv, p_ttv), olog)

    oc_file = os.path.join(
        cli.folder_path, "oc_%d_NB%d.dat" % (cli.sim_id, cli.body_id)
    )

    np.savetxt(
        oc_file,
        np.column_stack((epo_sel, transits.ttra_sel, oc_d)),
        fmt="%5d %23.16e %23.16e",
        header="Tlin = %23.16e + N x %23.16e\namp_TTV = %.12f min P_TTV = %.12f d\nepoch TT_d OC_d"
        % (transits.tref, transits.pref, amp_ttv * 1440.0, p_ttv),
    )

    plot_folder = os.path.join(os.path.dirname(full_file), "plots")
    if not os.path.isdir(plot_folder):
        os.makedirs(plot_folder)
    plot_file = os.path.join(plot_folder, "oc_%d_NB%d" % (cli.sim_id, cli.body_id))

    printlog("SAVING O-C PLOT %s.png\n" % (plot_file), olog)

    tsize = 9
    tsize2 = tsize - 3
    lsize = 10
    msize = 4

    fig = plt.figure(figsize=(6, 6))

    plt.axhline(0.0, color="black", ls="-", lw=0.33, alpha=0.77)
    plt.plot(
        tt,
        oc_d * ocu[0],
        color="lightgray",
        marker="o",
        ms=msize,
        mfc="C0",
        ls="-",
        lw=0.45,
    )

    if cli.name is None:
        plt.title(
            "%s (O-C) body %d" % (target, cli.body_id), loc="left", fontsize=tsize
        )
    else:
        plt.title("%s (O-C)" % (cli.name,), loc="left", fontsize=tsize)

    title1 = "Tlin = {:.6f} + N x {:.6f}".format(transits.tref, transits.pref)
    # plt.title(title1,
    # loc='center', fontsize=tsize2
    # )
    title2 = (
        r"$\sigma_\textrm{{TTV}}={0:.3f}$ {1:s} $P_\textrm{{TTV}}={2:.4f}$ d".format(
            amp_ttv * ocu[0], ocu[1], p_ttv
        )
    )
    # plt.title(title2,
    # loc='right', fontsize=tsize2
    # )
    plt.title("%s\n%s" % (title1, title2), loc="right", fontsize=tsize2)

    plt.xlabel("BJD - %.4f" % (transits.tref), fontsize=lsize)
    plt.ylabel("O-C ({})".format(ocu[1]), fontsize=lsize)

    fig.savefig("%s.png" % (plot_file), bbox_inches="tight", dpi=300)
    plt.close(fig)

    olog.close()

    return


# ==============================================================================
# ==============================================================================
if __name__ == "__main__":
    main()
