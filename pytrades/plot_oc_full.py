#!/usr/bin/env python
# -*- coding: utf-8 -*-

# ==============================================================================
# IMPORT

# no more "zero" integer division bugs!:P
import argparse
import numpy as np  # array
import os
# import sys
import glob
import gc

import matplotlib as mpl

mpl.use("Agg")
import matplotlib.pyplot as plt

# matplotlib rc params
plt.rcParams['text.usetex']       = False
# plt.rcParams['font.family']       = 'sans-serif'
plt.rcParams['font.family']       = 'serif'
plt.rcParams['font.serif']        = ['Computer Modern Roman', 'Palatino', 'DejaVu Serif']
plt.rcParams['mathtext.fontset']  = 'cm'
plt.rcParams['figure.figsize']    = [5, 5]
plt.rcParams["figure.facecolor"]  = 'white'
plt.rcParams["savefig.facecolor"] = 'white'
plt.rcParams["figure.dpi"]        = 200
plt.rcParams["savefig.dpi"]       = 300
plt.rcParams["font.size"]         = 12
plt.rcParams["xtick.labelsize"]   = plt.rcParams["font.size"] - 2
plt.rcParams["ytick.labelsize"]   = plt.rcParams["xtick.labelsize"]
# custom modules
import ancillary as anc
import gls

# ==============================================================================


def get_data_from_trafile(tra_file):

    tra_name = os.path.basename(tra_file)
    idlmbody = tra_name.split("_")
    sim_id = idlmbody[0]
    lm_flag = idlmbody[1]
    body_id = int(idlmbody[2].split("NB")[1])
    xx = np.genfromtxt(tra_file, usecols=(0, 1))
    TTs = np.sort(xx[:, 0] + xx[:, 1])

    return sim_id, lm_flag, body_id, TTs


# ==============================================================================


class Transit:
    def __init__(self, full_file):
        self.full_file = full_file
        self.dirname = os.path.dirname(full_file)
        self.filename = os.path.basename(full_file)
        self.sim_id, self.lm_flag, self.body_id, self.TTs = get_data_from_trafile(
            full_file
        )
        self.nTTs = np.shape(self.TTs)[0]
        # self.epo = np.arange(0,self.nTTs)
        # self.Tref = self.TTs[0]
        # self.Pref = self.TTs[1]-self.TTs[0]
        self.epo, self.Tref, self.Pref, self.TPerr = anc.compute_lin_ephem(self.TTs)
        self.linTTs = self.Tref + self.epo * self.Pref
        self.ocd = self.TTs - self.linTTs
        self.amp_ttv = 0.0
        self.p_ttv = 0.0

    def compute_TTV_gls(self):
        dur = self.TTs[-1] - self.TTs[0]
        e1sec = np.ones(np.shape(self.TTs)) / 86400.0
        ogls = gls.Gls(
            (self.TTs, self.ocd, e1sec),
            Pbeg=3.0 * self.Pref,
            Pend=2.0 * dur,
            ofac=10,
            verbose=False,
        )
        self.amp_ttv = ogls.hpstat["amp"]
        self.p_ttv = 1.0 / ogls.hpstat["fbest"]
        del ogls
        gc.collect()


# ==============================================================================


def get_tra_file_list(folder_path, sim_id, lm_flag):

    pattern = "{0:s}_{1:s}_NB*_tra.dat".format(str(sim_id), lm_flag)
    tra_file_list = np.sort(glob.glob(os.path.join(folder_path, pattern)))

    return tra_file_list


# ==============================================================================


def get_args():

    parser = argparse.ArgumentParser(description="PLOT OC")

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
        "-l",
        "--l",
        "-lm-flag",
        "--lm-flag",
        action="store",
        dest="lm_flag",
        default="0",
        help="LM flag: 0 or 1. Default is 0.",
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

    cli = parser.parse_args()
    cli.folder_path = os.path.abspath(cli.folder_path)

    try:
        cli.sim_id = int(cli.sim_id)
        if cli.sim_id <= 0:
            cli.sim_id = 1
    except:
        cli.sim_id = 1

    if str(cli.lm_flag) == "1":
        cli.lm_flag = "1"
    else:
        cli.lm_flag = "0"

    return cli


# ==============================================================================


def main():

    cli = get_args()

    tra_file_list = get_tra_file_list(cli.folder_path, cli.sim_id, cli.lm_flag)

    n_files = len(tra_file_list)
    viridis = plt.get_cmap("viridis")
    colors = np.linspace(0.2, 1.0, endpoint=False, num=n_files)


    plot_folder = os.path.join(cli.folder_path, "plots")
    if not os.path.isdir(plot_folder):
        os.makedirs(plot_folder)
    plot_file = os.path.join(
        plot_folder, "{0:d}_{1:s}_oc_full.png".format(cli.sim_id, cli.lm_flag)
    )

    # fig = plt.figure(figsize=(12,12))
    # fig = plt.figure(figsize=(6,6))
    fig = plt.figure()
    ax = plt.subplot2grid((1,1), (0,0), fig=fig)

    # plt.tick_params(direction="in")

    ax.axhline(0.0, color="black", ls="-", lw=0.33, alpha=0.77)

    tra_obj = []
    cnt = 0
    for i_f, ff in enumerate(tra_file_list):
        print("READ FILE {}".format(ff))
        oo = Transit(ff)
        oo.compute_TTV_gls()
        tra_obj.append(oo)

        print("body n. {0:d}".format(oo.body_id))
        # print 'ephem: %20.8f + N x %20.8f'.format(oo.Tref, oo.Pref)

        # single planet O-C
        figx = plt.figure()
        axx = plt.subplot2grid((1,1), (0,0), fig=figx)
        axx.axhline(0.0, color="black", ls="-", lw=0.33, alpha=0.77)
        
        color = viridis(colors[i_f])

        # for both plots
        xscale = 0.0  # oo.TTs[0]
        if cli.tscale is not None:
            #xscale = oo.TTs - float(cli.tscale)
            xscale = float(cli.tscale)
            xlabel = "BJD - {0:.3f}".format(xscale)
        else:
            xlabel = "BJD"
        ylabel = "O-C (min)"
        axx.set_xlabel(xlabel)
        if cnt == 0:
            ax.set_xlabel(xlabel)
            ax.set_ylabel(ylabel)
        axx.set_ylabel(ylabel)
        print("plotting O-C")
        cnt += 1
        labTeph = "TTlin $= {0:17.6f} + N \\times {1:13.6f}$".format(oo.Tref, oo.Pref)
        print(labTeph)
        labPTTV = "$P_\mathrm{{TTV}} = {0:13.6f}$ days".format(oo.p_ttv)
        labATTV = "$A_\mathrm{{TTV}} = {0:10.3f}$ min".format(oo.amp_ttv * 1440.0)
        print("{} , {}".format(labPTTV, labATTV))
        ax.plot(
            oo.TTs - xscale,
            oo.ocd * 1440.0,
            color=color,
            marker="o",
            ms=4,
            mec="white",
            mew=0.3,
            ls="-",
            lw=0.45,
            alpha=0.77,
            label="body #{}: {} ; {} , {}".format(
                oo.body_id, labTeph, labPTTV, labATTV
            ),
        )
        axx.plot(
            oo.TTs - xscale,
            oo.ocd * 1440.0,
            color=color,
            marker="o",
            ms=4,
            mec="white",
            mew=0.3,
            ls="-",
            lw=0.45,
            alpha=1.0,
            label="body #{}: {} ; {} , {}".format(
                oo.body_id, labTeph, labPTTV, labATTV
            ),
        )
        axx.legend(loc="lower center", bbox_to_anchor=(0.5, 1.0), fontsize=6, ncol=1)
        
        plot_filex = plot_file.replace("full", "NB{}".format(oo.body_id))
        plt.tight_layout()
        print("Saving file: {}".format(plot_filex))
        figx.savefig(plot_filex, bbox_inches='tight')
        plt.close(figx)

    print("Saving file: {}".format(plot_file))
    ax.legend(loc="lower center", bbox_to_anchor=(0.5, 1.0), fontsize=6, ncol=1)

    # fig.savefig(plot_file, dpi=300)
    plt.tight_layout()
    fig.savefig(plot_file, bbox_inches='tight')
    plt.close(fig)
    print("DONE")

    return


# ==============================================================================
# ==============================================================================

if __name__ == "__main__":
    main()

# ==============================================================================

