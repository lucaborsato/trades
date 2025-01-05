#!/usr/bin/env python
# -*- coding: utf-8 -*-

# no more "zero" integer division bugs!:P
import argparse
import os
import numpy as np  # array
import glob
import pandas as pd
import sys

# import h5py

import matplotlib as mpl

# mpl.use("Agg")
import matplotlib.pyplot as plt

# local constants module
from . import constants as cst
from . import ancillary as anc

anc.set_rcParams()
# ======================================================================


class PSOEvolution:
    def __init__(self, pso_path):

        self.pso_path = pso_path

        all_pattern = os.path.join(pso_path, "*_allpso_*par.out")
        all_file = glob.glob(all_pattern)[0]
        print("Reading:\n{}".format(all_file))
        self.all_pso = pd.DataFrame(np.genfromtxt(all_file, names=True))
        self.all_names = self.all_pso.columns.values.tolist()

        self.last_iter = self.all_pso["iter"].max()

        best_pattern = os.path.join(pso_path, "*_bestpso_*par.out")
        best_file = glob.glob(best_pattern)[0]
        print("Reading:\n{}".format(best_file))
        self.best_pso = pd.DataFrame(np.genfromtxt(best_file, names=True))
        self.best_names = self.best_pso.columns.values.tolist()

    def convert_masses(self, Mstar=1.0, mass_unit="e"):

        if mass_unit.lower() in ["j", "ju", "jup"]:
            u = cst.Msjup
        elif mass_unit.lower() in ["n", "ne", "nep"]:
            u = cst.Msnep
        else:
            u = cst.Msear

        for i_name, name in enumerate(self.all_names):

            if "Ms" in name:
                new_name = name.replace("Ms", "").replace("m", "mass")
                self.all_pso[new_name] = self.all_pso[name] * Mstar * u
                if name == self.best_names[i_name]:
                    self.best_pso[new_name] = self.best_pso[name] * Mstar * u

        return

    def get_eccentricity(self):

        for i_name, name in enumerate(self.all_names):
            if "ecosw" in name:
                name2 = self.all_names[i_name + 1]
                square_all = self.all_pso[name] ** 2 + self.all_pso[name2] ** 2
                square_best = self.best_pso[name] ** 2 + self.best_pso[name2] ** 2
                if name[0:2] == "ec":
                    new_name = name.replace("ecosw", "ecc")
                    self.all_pso[new_name] = np.sqrt(square_all)
                    self.best_pso[new_name] = np.sqrt(square_best)
                else:
                    new_name = name.replace("secosw", "ecc")
                    self.all_pso[new_name] = square_all
                    self.best_pso[new_name] = square_best

        return

    def ax_scatter_plot(
        self,
        fig,
        ax,
        name1,
        name2,
        scale_name="pso_fitness",
        colorbar=True,
        thin_iter=1,
        percentile_max=50,
    ):

        max_iter = self.last_iter
        pso_iter = self.all_pso["iter"]
        if thin_iter > 1:
            plot_iter = np.arange(1, max_iter, thin_iter)
            if plot_iter[-1] < max_iter:
                plot_iter = np.concatenate((plot_iter, [max_iter]))
            sel_iter = np.in1d(pso_iter, plot_iter)
        else:
            sel_iter = np.ones((len(pso_iter))).astype(bool)

        x = self.all_pso[name1][sel_iter]
        y = self.all_pso[name2][sel_iter]
        # c = self.all_pso[scale_name][sel_iter]
        scale_full = self.all_pso[scale_name].copy()
        scale_perc = np.percentile(scale_full, percentile_max)
        scale_full[scale_full < scale_perc] = scale_perc
        c = scale_full[sel_iter]

        alpha = 1.0

        sc = ax.scatter(x, y, s=2, c=c, alpha=alpha)

        xbest = self.best_pso[name1]
        ybest = self.best_pso[name2]
        ax.plot(xbest, ybest, color="C0", marker="None", ls="-", lw=1.5)

        msize = 5
        ax.plot(
            xbest.iloc[0],
            ybest.iloc[0],
            color="black",
            marker="s",
            ms=msize,
            mec="white",
            mew=0.4,
            ls="",
            lw=1.5,
        )
        ax.plot(
            xbest.iloc[-1],
            ybest.iloc[-1],
            color="black",
            marker="^",
            ms=msize,
            mec="white",
            mew=0.4,
            ls="",
            lw=1.5,
        )

        ax.set_xlabel(name1)
        ax.set_ylabel(name2)
        if colorbar:
            # Adding the colorbar
            # left, bottom, width, height
            cbaxes = fig.add_axes([0.92, 0.05, 0.03, 0.85])
            plt.colorbar(
                sc,
                ax=ax,
                # location='right',
                cax=cbaxes,
                label=scale_name,
            )

        return

    def scatter_plot(
        self,
        name1,
        name2,
        scale_name="pso_fitness",
        colorbar=True,
        thin_iter=1,
        percentile_max=50,
    ):

        fig = plt.figure()
        nrow, ncol = 1, 1
        irow, icol = 0, 0
        ax = plt.subplot2grid((nrow, ncol), (irow, icol))
        self.ax_scatter_plot(
            fig,
            ax,
            name1,
            name2,
            scale_name=scale_name,
            colorbar=colorbar,
            thin_iter=thin_iter,
            percentile_max=percentile_max,
        )
        plt.show()
        plt.close(fig)

        return

    def threepanel_scatter_plot(
        self,
        three_pairs=[["P2", "P3"], ["mass2", "mass3"], ["ecc2", "ecc3"]],
        scale_name="pso_fitness",
        thin_iter=1,
        percentile_max=50,
    ):

        # c = self.all_pso[scale_name]

        best_fitness = self.best_pso["pso_fitness"].max()
        best_bic = self.best_pso["bic"].min()

        fig = plt.figure(figsize=(8, 5))

        title = "PSO FITNESS = {} BIC = {}".format(best_fitness, best_bic)
        print(title)
        plt.suptitle(title, fontsize=plt.rcParams["font.size"] - 4)

        # nrow, ncol = 1, 3

        # irow, icol = 0, 0
        # ax1 = plt.subplot2grid((nrow, ncol), (irow, icol))
        dl = 0.06
        l, b, w, h = dl + 0.0, 0.05, 0.22, 0.85
        ax1 = fig.add_axes([l, b, w, h])
        self.ax_scatter_plot(
            fig,
            ax1,
            three_pairs[0][0],
            three_pairs[0][1],
            scale_name=scale_name,
            colorbar=False,
            thin_iter=thin_iter,
            percentile_max=percentile_max,
        )

        # irow, icol = 0, 1
        # ax2 = plt.subplot2grid((nrow, ncol), (irow, icol))
        l += w + dl
        ax2 = fig.add_axes([l, b, w, h])
        self.ax_scatter_plot(
            fig,
            ax2,
            three_pairs[1][0],
            three_pairs[1][1],
            scale_name=scale_name,
            colorbar=False,
            thin_iter=thin_iter,
            percentile_max=percentile_max,
        )

        # irow, icol = 0, 2
        # ax3 = plt.subplot2grid((nrow, ncol), (irow, icol))
        l += w + dl
        ax3 = fig.add_axes([l, b, w, h])
        self.ax_scatter_plot(
            fig,
            ax3,
            three_pairs[2][0],
            three_pairs[2][1],
            scale_name=scale_name,
            colorbar=True,
            thin_iter=thin_iter,
            percentile_max=percentile_max,
        )

        # plt.tight_layout()
        plt.show()
        plt.close(fig)

        return

    def adhoc_scatter_plot(
        self, scale_name="pso_fitness", thin_iter=1, percentile_max=50
    ):

        self.threepanel_scatter_plot(
            scale_name=scale_name, thin_iter=thin_iter, percentile_max=percentile_max
        )

        return

    def histogram_plot(self, n_planets, label_type="mass"):

        fig = plt.figure()

        k = "pso_fitness"
        fitness = self.all_pso[k]
        nall = len(fitness)
        p50 = np.percentile(fitness, 50, interpolation="midpoint")
        sel_k = self.all_pso[k] >= p50

        viridis = plt.get_cmap("viridis")
        inferno = plt.get_cmap("inferno")
        nc = n_planets
        seq = np.linspace(0.0, 1.0, endpoint=True, num=nc)
        colors = viridis(seq)
        # colors_r = colors[::-1]
        colors_r = inferno(seq[::-1])

        for i_body in range(0, nc):
            k = "{}{}".format(label_type, i_body + 2)
            color = colors[i_body]
            color_r = colors_r[i_body]
            plt.hist(self.all_pso[k], bins=100, color=color, histtype="bar", label=k)
            plt.hist(
                self.all_pso[k][sel_k],
                bins=100,
                color=color_r,
                histtype="step",
                label="{} sel lnL".format(k),
                alpha=0.8,
            )

        plt.xlabel(label_type)

        plt.legend(loc="best", fontsize=8)
        plt.show()
        plt.close(fig)

        return

    def evolution_plot_one_par(
        self, label, true_pars_dictionary={}, thin_iter=1, percentile_max=50
    ):

        k = label

        # thin the data by iter
        max_iter = self.last_iter
        iters = self.all_pso["iter"]
        if thin_iter > 1:
            plot_iter = np.arange(1, max_iter, thin_iter)
            if plot_iter[-1] < max_iter:
                plot_iter = np.concatenate((plot_iter, [max_iter]))
            sel_iter = np.in1d(iters, plot_iter)
        else:
            sel_iter = np.ones((len(iters))).astype(bool)

        # re-scale by fitness
        kf = "pso_fitness"
        scale_full = self.all_pso[kf].copy()
        scale_perc = np.percentile(scale_full, percentile_max)
        scale_full[scale_full < scale_perc] = scale_perc

        x = iters[sel_iter]
        y = self.all_pso[k][sel_iter]
        c = scale_full[sel_iter]

        fig = plt.figure(figsize=(6, 4))

        sc = plt.scatter(x, y, c=c, s=2, alpha=0.2)
        plt.colorbar(sc, label=kf)
        plt.plot(
            self.best_pso["iter"],
            self.best_pso[k],
            color="C3",
            marker="None",
            ls="-",
            lw=1.3,
        )
        if len(true_pars_dictionary) > 0 and k in true_pars_dictionary.keys():
            plt.axhline(true_pars_dictionary[k], color="C1", ls="--", lw=1.5)
        plt.xlabel("N iteration")
        plt.ylabel(k)

        return fig


# ======================================================================

# read command line (cli) arguments
def get_args():
    parser = argparse.ArgumentParser(description="PSO PLOT")
    parser.add_argument(
        "-p",
        "--path",
        action="store",
        dest="full_path",
        required=True,
        help="The path (absolute or relative) with simulation files.",
    )
    parser.add_argument(
        "-m",
        "--mtype",
        "--mass-type",
        action="store",
        dest="m_type",
        default="e",
        help="Mass type: j = Jupiter, e = Earth, s = Sun. Default is Earth = e.",
    )

    cli = parser.parse_args()
    cli.full_path = os.path.abspath(cli.full_path)
    cli.m_type = str(cli.m_type).lower()

    return cli


# ======================================================================


def main():

    print()
    print(" ================== ")
    print(" PSO PLOTS")
    print(" ================== ")
    print()

    # read cli arguments
    cli = get_args()

    # set pso_file
    pso_file = os.path.join(cli.full_path, "pso_run.hdf5")
    (
        population,
        population_fitness,
        pso_parameters,
        pso_fitness,
        pso_best_evolution,
        parameters_minmax,
        parameter_names,
        pop_shape,
    ) = anc.get_pso_data(pso_file)
    nfit = pop_shape[0]
    npop = pop_shape[1]
    niter = pop_shape[2]

    nwalkers = npop

    mass_star = 1.0
    m_factor = 1.0
    # m_factor, m_unit = anc.mass_type_factor(
    #     Ms=mass_star, mtype=cli.m_type, mscale=False
    # )
    m_factor, m_unit, r_factor, r_unit = anc.mass_radius_type_factor(mtype=cli.m_type)

    iteration = np.arange(0, niter) + 1
    if isinstance(parameters_minmax, type(population_fitness)):
        parameters_minmax_bck = parameters_minmax.copy()

    # set label and legend names
    kel_legends = anc.keplerian_legend(parameter_names, cli.m_type)

    anc.print_memory_usage(population)
    anc.print_memory_usage(population_fitness)
    fit_min = np.min(population_fitness)
    fit_max = fit_min + 100
    print(
        "population fitness: min = {} (set max = min + 100 = {}".format(
            fit_min, fit_max
        )
    )

    pso_plots = os.path.join(cli.full_path, "plots")
    if not os.path.isdir(pso_plots):
        os.makedirs(pso_plots)

    # parameter_names and parameters_minmax in pso_run.hdf5
    for ii in range(0, nfit):
        print("parameter: {:s}".format(parameter_names[ii]))
        if parameter_names[ii][0] == "m" and parameter_names[ii][1] != "A":
            pop_plt = population[ii, :, :] * m_factor
            y_min = parameters_minmax[ii, 0] * m_factor
            y_max = parameters_minmax[ii, 1] * m_factor
        else:
            pop_plt = population[ii, :, :]
            y_min = parameters_minmax[ii, 0]
            y_max = parameters_minmax[ii, 1]

        print("boundaries: [{:.6f}, {:.6f}]".format(y_min, y_max))
        print(
            "    minmax: [{:.6f}, {:.6f}]".format(
                np.min(population[ii, :, :]), np.max(population[ii, :, :])
            )
        )
        pso_fig_file = os.path.join(
            pso_plots, "pso_evolution_{:s}.png".format(parameter_names[ii])
        )
        print(" {:s}".format(pso_fig_file), end=" ")

        fig = plt.figure(figsize=(5, 5))
        for jj in range(0, npop):
            plt.plot(
                iteration,
                pop_plt[jj, :],
                marker="o",
                mfc=(0.8, 0.8, 0.8, 0.8),
                mec="None",
                ls="",
                ms=0.5,
                zorder=6,
            )
            # # SCATTER plot
            # sc = plt.scatter(iteration, pop_plt[jj,:],
            #   s=1,
            #   c=population_fitness[jj,:],
            #   marker='.',
            #   vmin=fit_min, vmax=fit_max,
            #   alpha=0.45
            # )

        plt.ylim(y_min, y_max)

        if parameter_names[ii][0] == "m" and parameter_names[ii][1] != "A":
            par_plt = pso_best_evolution[ii, :] * m_factor
        else:
            par_plt = pso_best_evolution[ii, :]

        # plt.plot(iteration, par_plt,
        #   color='C0',
        #   marker='o',
        #   mfc='C0',
        #   mec='white',
        #   mew=0.1,
        #   ls='-',
        #   ms=3,
        #   alpha=0.45
        # )

        # sel_min = np.argmin(population_fitness, axis=0)
        par_fit = np.min(population_fitness, axis=0)
        sc = plt.scatter(
            iteration,
            par_plt,
            s=2,
            c=par_fit,
            marker="o",
            vmin=fit_min,
            vmax=fit_max,
            alpha=1.0,
            zorder=7,
        )

        plt.colorbar(sc)

        plt.xlabel("$N_\mathrm{iteration}$")
        plt.ylabel(kel_legends[ii])
        plt.draw()
        fig.savefig(pso_fig_file, bbox_inches="tight")
        print(" done")

    return


if __name__ == "__main__":
    main()

