import os
import numpy as np  # array
import pandas as pd
import h5py

# import h5py

import matplotlib as mpl

# mpl.use("Agg")
import matplotlib.pyplot as plt

# # matplotlib rc params
# plt.rcParams["text.usetex"] = False
# # plt.rcParams['font.family']       = 'sans-serif'
# plt.rcParams["font.family"] = "serif"
# plt.rcParams["font.serif"] = ["Computer Modern Roman", "Palatino", "DejaVu Serif"]
# plt.rcParams["mathtext.fontset"] = "cm"
# plt.rcParams["figure.figsize"] = [5, 5]
# plt.rcParams["figure.facecolor"] = "white"
# plt.rcParams["savefig.facecolor"] = "white"
# plt.rcParams["figure.dpi"] = 200
# plt.rcParams["savefig.dpi"] = 300
# plt.rcParams["font.size"] = 12
# plt.rcParams["xtick.labelsize"] = plt.rcParams["font.size"] - 2
# plt.rcParams["ytick.labelsize"] = plt.rcParams["xtick.labelsize"]


# local constants module
import constants as cst
import ancillary as anc

anc.set_rcParams()

# ======================================================================
class DEEvolution:
    def __init__(self, run_folder):
        de_file = os.path.join(run_folder, "de_run.hdf5")

        h5f = h5py.File(
            de_file,
            "r",
            # libver='latest', swmr=True
        )

        self.npop_de = h5f["population"].attrs["npop"]
        self.ngen_de = h5f["population"].attrs["ngen"]
        self.iter_de = h5f["population"].attrs["iter_de"]
        self.iter_global = h5f["population"].attrs["iter_global"]
        self.nfit = h5f["population"].attrs["nfit"]
        self.ndata = h5f["population"].attrs["ndata"]
        print("Read de_run.hdf5 file with:")
        print("npop_de = {}".format(self.npop_de))
        print("ngen_de = {}".format(self.ngen_de))
        print("iter_de = {}".format(self.iter_de))
        print("iter_global = {}".format(self.iter_global))
        print("nfit = {}".format(self.nfit))
        print("ndata = {}".format(self.ndata))
        print()

        if self.iter_de == 0:
            print("The file has been only initialised, wait a bit more ...")
            return
        if self.iter_de + 1 < self.ngen_de:
            self.ngen_de = self.iter_de + 1
            print("updated ngen_de = {}".format(self.ngen_de))

        self.de_pop = np.array(h5f["population"][: self.ngen_de, :, :])
        self.de_fit = h5f["population_fitness"][: self.ngen_de, :]

        self.de_pop_best = h5f["best_population"][: self.ngen_de, :]
        self.de_fit_best = h5f["best_fitness"][: self.ngen_de]

        self.de_par = h5f["de_parameters"][...]
        if np.all(self.de_par == 0.0):
            loc_best = np.argmax(self.de_fit_best[: self.ngen_de])
            self.de_par = self.de_pop_best[loc_best, :].copy()
        self.de_bounds = h5f["parameters_minmax"][...]
        self.par_names = anc.decode_list(h5f["parameter_names"])

        h5f.close()

        de_pop_flat = self.de_pop.copy().reshape(
            (self.ngen_de * self.npop_de, self.nfit)
        )
        iters = np.repeat(np.arange(self.ngen_de), self.npop_de)
        d = {}
        d["iter"] = iters
        for i_p, p in enumerate(self.par_names):
            d[p] = de_pop_flat[:, i_p]
        d["de_fitness"] = self.de_fit.copy().reshape((self.ngen_de * self.npop_de))
        # for k, v in d.items():
        #     print(k, len(v))

        self.de_pop_flat = pd.DataFrame(d)

        d = {}
        d["iter"] = np.arange(self.ngen_de)
        for i_p, p in enumerate(self.par_names):
            d[p] = self.de_pop_best[:, i_p].copy()
        d["de_fitness"] = self.de_fit_best.copy()
        # print()
        # for k, v in d.items():
        #     print(k, len(v))
        self.de_best_flat = pd.DataFrame(d)

    def convert_masses(self, Mstar=1.0, mass_unit="e"):

        if mass_unit.lower() in ["j", "ju", "jup"]:
            u = cst.Msjup
        elif mass_unit.lower() in ["n", "ne", "nep"]:
            u = cst.Msnep
        else:
            u = cst.Msear

        for i_name, name in enumerate(self.par_names):

            if "Ms" in name:
                new_name = name.replace("Ms", "").replace("m", "mass")
                self.de_pop_flat[new_name] = self.de_pop_flat[name] * Mstar * u
                self.de_best_flat[new_name] = self.de_best_flat[name] * Mstar * u

        return

    def convert_jitter(self):

        for i_name, name in enumerate(self.par_names):

            if "l2j" in name:
                new_name = name.replace("l2j", "jitter")
                self.de_pop_flat[new_name] = np.power(2, self.de_pop_flat[name])
                self.de_best_flat[new_name] = np.power(2, self.de_best_flat[name])

        return

    def recenter_angles(self):

        for i_name, name in enumerate(self.par_names):

            if (
                name[0] == "w"
                or name[0:2] == "mA"
                or name[0:2] == "lN"
                or "lambda" in name
            ):
                new_name = "rec_{}".format(name)
                self.de_pop_flat[new_name] = anc.recenter_angle_distribution(
                    self.de_pop_flat[name] % 360.0
                )
                self.de_best_flat[new_name] = anc.recenter_angle_distribution(
                    self.de_best_flat[name] % 360.0
                )

        return

    def get_eccentricity(self):

        for i_name, name in enumerate(self.par_names):
            if "ecosw" in name:
                name2 = self.par_names[i_name + 1]
                square_all = self.de_pop_flat[name] ** 2 + self.de_pop_flat[name2] ** 2
                angle_all = anc.recenter_angle_distribution(
                    (
                        np.arctan2(self.de_pop_flat[name2], self.de_pop_flat[name])
                        * cst.rad2deg
                    )
                    % 360,
                )
                square_best = (
                    self.de_best_flat[name] ** 2 + self.de_best_flat[name2] ** 2
                )
                angle_best = anc.recenter_angle_distribution(
                    (
                        np.arctan2(self.de_best_flat[name2], self.de_best_flat[name])
                        * cst.rad2deg
                    )
                    % 360,
                )
                if name[0:2] == "ec":
                    new_name = name.replace("ecosw", "ecc")
                    self.de_pop_flat[new_name] = np.sqrt(square_all)
                    self.de_best_flat[new_name] = np.sqrt(square_best)
                else:
                    new_name = name.replace("secosw", "ecc")
                    self.de_pop_flat[new_name] = square_all
                    self.de_best_flat[new_name] = square_best
                new_name = "w{}".format(name.split("w")[1])
                self.de_pop_flat[new_name] = angle_all
                self.de_best_flat[new_name] = angle_best

        return

    def ax_scatter_plot(
        self,
        fig,
        ax,
        name1,
        name2,
        scale_name="de_fitness",
        colorbar=True,
        thin_iter=1,
        percentile_max=50,
    ):

        max_iter = self.ngen_de
        de_iter = self.de_pop_flat["iter"]
        if thin_iter > 1:
            plot_iter = np.arange(1, max_iter, thin_iter)
            if plot_iter[-1] < max_iter:
                plot_iter = np.concatenate((plot_iter, [max_iter]))
            sel_iter = np.in1d(de_iter, plot_iter)
        else:
            sel_iter = np.ones((len(de_iter))).astype(bool)

        x = self.de_pop_flat[name1][sel_iter]
        y = self.de_pop_flat[name2][sel_iter]
        scale_full = self.de_pop_flat[scale_name].copy()
        scale_perc = np.percentile(scale_full, percentile_max)
        scale_full[scale_full < scale_perc] = scale_perc
        c = scale_full[sel_iter]

        alpha = 1.0

        sc = ax.scatter(x, y, s=2, c=c, alpha=alpha)

        xbest = self.de_best_flat[name1]
        ybest = self.de_best_flat[name2]
        ax.plot(xbest, ybest, color="C0", marker="None", ls="-", lw=1.5)

        msize = 5
        ax.plot(
            xbest.iloc[0],
            ybest.iloc[0],
            color="black",
            marker="s",
            ms=msize,
            mec="white",
            mew=0.8,
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
        scale_name="de_fitness",
        colorbar=True,
        thin_iter=1,
        percentile_max=50,
        show_plot=False,
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
        if show_plot:
            plt.show()
        # plt.close(fig)

        return fig

    def threepanel_scatter_plot(
        self,
        three_pairs=[["P2", "P3"], ["mass2", "mass3"], ["ecc2", "ecc3"]],
        scale_name="de_fitness",
        thin_iter=1,
        percentile_max=50,
        show_plot=False,
    ):

        # c = self.all_pso[scale_name]

        best_fitness = self.de_best_flat["de_fitness"].max()
        best_bic = -2 * best_fitness + self.nfit * np.log(self.ndata)

        fig = plt.figure(figsize=(8, 5))

        title = "DE FITNESS = {} BIC = {}".format(best_fitness, best_bic)
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
        if show_plot:
            plt.show()
        # plt.close(fig)

        return fig

    def adhoc_scatter_plot(
        self, scale_name="de_fitness", thin_iter=1, percentile_max=50, show_plot=False
    ):

        fig = self.threepanel_scatter_plot(
            scale_name=scale_name,
            thin_iter=thin_iter,
            percentile_max=percentile_max,
            show_plot=show_plot,
        )

        return fig

    def histogram_plot(self, n_planets, label_type="mass", show_plot=False):

        fig = plt.figure()

        k = "de_fitness"
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
            plt.hist(
                self.de_pop_flat[k], bins=100, color=color, histtype="bar", label=k
            )
            plt.hist(
                self.de_pop_flat[k][sel_k],
                bins=100,
                color=color_r,
                histtype="step",
                label="{} sel lnL".format(k),
                alpha=0.8,
            )

        plt.xlabel(label_type)

        plt.legend(loc="best", fontsize=8)
        if show_plot:
            plt.show()
        # plt.close(fig)

        return fig

    def evolution_plot_one_par(
        self,
        label,
        true_pars_dictionary={},
        thin_iter=1,
        percentile_max=50,
        show_plot=False,
    ):

        k = label

        # thin the data by iter
        max_iter = self.ngen_de
        de_iter = self.de_pop_flat["iter"]
        if thin_iter > 1:
            plot_iter = np.arange(1, max_iter, thin_iter)
            if plot_iter[-1] < max_iter:
                plot_iter = np.concatenate((plot_iter, [max_iter]))
            sel_iter = np.in1d(de_iter, plot_iter)
        else:
            sel_iter = np.ones((len(de_iter))).astype(bool)

        # re-scale by fitness
        kf = "de_fitness"
        scale_full = self.de_pop_flat[kf].copy()
        scale_perc = np.percentile(scale_full, percentile_max)
        scale_full[scale_full < scale_perc] = scale_perc

        x = de_iter[sel_iter]
        y = self.de_pop_flat[k][sel_iter]
        c = scale_full[sel_iter]

        fig = plt.figure(figsize=(6, 4))
        plt.title(
            "DE best {} = {:.8f}".format(k, self.de_best_flat[k].iloc[-1]),
            loc="center",
            # fontdict={"fontsize": 12}
        )
        sc = plt.scatter(x, y, c=c, s=2, alpha=0.2)
        cbar = plt.colorbar(sc, label=kf)
        cbar.solids.set(alpha=1)
        plt.plot(
            self.de_best_flat["iter"],
            self.de_best_flat[k],
            color="C3",
            marker="None",
            ls="-",
            lw=1.3,
        )
        if len(true_pars_dictionary) > 0 and k in true_pars_dictionary.keys():
            plt.axhline(true_pars_dictionary[k], color="C1", ls="--", lw=1.5)
        plt.xlabel("N iteration")
        plt.ylabel(k)
        if show_plot:
            plt.show()

        return fig
