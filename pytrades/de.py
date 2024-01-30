import os
import numpy as np  # array
import pandas as pd
import h5py
import time
import sys

# import h5py

import matplotlib as mpl

# mpl.use("Agg")
import matplotlib.pyplot as plt


# local constants module
import constants as cst
import ancillary as anc

from multiprocessing import Pool
from pyde import DiffEvol

anc.set_rcParams()

# =============================================================================
def de_get_best_loc(de_fit_best, iter_de, de_maximize=True):

    loc_best = (
        np.argmax(de_fit_best[: iter_de + 1])
        if de_maximize
        else np.argmin(de_fit_best[: iter_de + 1])
    )

    return loc_best

def de_get_best_parameters(de_fit_best, de_pop_best, iter_de, de_maximize=True):

    loc_best = de_get_best_loc(de_fit_best, iter_de, de_maximize=de_maximize)

    de_parameters = de_pop_best[loc_best, :].copy()

    return de_parameters

def de_hdf5_save_one_dataset(de_file, data_name, data, data_type, hdf5_mode="a"):

    de_hdf5 = h5py.File(
        de_file,
        hdf5_mode,
        # libver='latest'
    )
    # de_hdf5.swmr_mode = True # ERRORS WHEN ADDING/DELETING/UPDATING DATASETS AND ATTRS!!
    if data_name in de_hdf5:
        del de_hdf5[data_name]
    de_hdf5.create_dataset(data_name, data=data, dtype=data_type, compression="gzip")
    de_hdf5.close()

    return


def de_hdf5_update_attr(de_file, data_name, attrs_name, attrs_value):

    de_hdf5 = h5py.File(
        de_file,
        "r+",
        # libver='latest'
    )
    # de_hdf5.swmr_mode = True
    de_hdf5[data_name].attrs[attrs_name] = attrs_value
    de_hdf5.close()

    return


def de_save_evolution(
    de_file,
    npop_de,
    ngen_de,
    iter_de,
    iter_global,
    nfit,
    ndata,
    de_pop,
    de_fit,
    de_pop_best,
    de_fit_best,
    de_bounds,
    parameter_names,
    de_maximize=True,
):

    # de_file = os.path.join(de_path, 'de_run.hdf5')
    de_hdf5_save_one_dataset(de_file, "population", de_pop, "f8")
    de_hdf5_update_attr(de_file, "population", "npop", npop_de)
    de_hdf5_update_attr(de_file, "population", "ngen", ngen_de)
    de_hdf5_update_attr(de_file, "population", "iter_de", iter_de)
    de_hdf5_update_attr(de_file, "population", "iter_global", iter_global)
    de_hdf5_update_attr(de_file, "population", "nfit", nfit)
    de_hdf5_update_attr(de_file, "population", "ndata", ndata)

    de_hdf5_save_one_dataset(de_file, "population_fitness", de_fit, "f8")

    de_hdf5_save_one_dataset(de_file, "best_population", de_pop_best, "f8")

    de_hdf5_save_one_dataset(de_file, "best_fitness", de_fit_best, "f8")

    de_hdf5_save_one_dataset(de_file, "parameters_minmax", de_bounds, "f8")

    # best_loc = np.argmax(de_fit_best[:iter_de])
    # de_hdf5_save_one_dataset(de_file, 'de_parameters', de_pop_best[best_loc,:].copy(), 'f8')
    de_parameters = de_get_best_parameters(
        de_fit_best, de_pop_best, iter_de, de_maximize=de_maximize
    )
    de_hdf5_save_one_dataset(de_file, "de_parameters", de_parameters, "f8")

    de_hdf5_save_one_dataset(de_file, "parameter_names", parameter_names, "S10")

    return

def de_load_parameters(de_file):

    with h5py.File(
        de_file,
        mode="r",
        # libver='latest',
        swmr=True,
    ) as de_hdf5:
        de_parameters = de_hdf5["de_parameters"][...]
    # de_hdf5.close()

    return de_parameters

# ======================================================================
def run_de(conf, sim, fitness_func, working_folder, of_run=None):

    de_path = working_folder
    de_file = os.path.join(de_path, "de_run.hdf5")
    de_file_check = os.path.exists(de_file) and os.path.isfile(de_file)

    de_run      = conf.de_type.lower() == "run"
    de_resume   = conf.de_type.lower() == "resume"
    de_to_emcee = conf.de_type.lower() == "to_emcee"
    de_type     = conf.de_type.lower()

    npop_de = conf.npop_de
    ngen_de = conf.ngen_de
    nfit = sim.nfit
    ndata = sim.ndata

    nthreads = conf.nthreads

    de_save = conf.nsave_de
    de_bounds = sim.fitting_minmax
    de_f = conf.de_f
    de_c = conf.de_c
    de_maximize = conf.de_maximize
    fit_type = -1 if de_maximize else 1


    anc.print_both(
        "pyDE de_type = {} ==> de_run = {}, de_resume = {}, de_to_emcee = {}".format(
        de_type, de_run, de_resume, de_to_emcee
        ),
        output=of_run,
    )

    anc.print_both(
        "pyDE npop = {} ngen = {} nsave = {}".format(npop_de, ngen_de, de_save),
        output=of_run,
    )

    # prepare arrays to store pyDE evolution
    de_pop, de_fit = (
        np.zeros((ngen_de, npop_de, nfit)),
        np.zeros((ngen_de, npop_de)) - 1.0e-30,  # set to very small non-zero values
    )
    de_pop_best, de_fit_best = (
        np.zeros((ngen_de, nfit)),
        np.zeros((ngen_de)) - 1.0e-30,
    )
    last_iter_de = 0

    if de_to_emcee or de_resume:
        if de_file_check:
            anc.print_both(
                "pyDE file already exists --> loading it ...", output=of_run,
            )
            de_obj = DEEvolution(de_path, of_run=of_run)
            de_pop[: de_obj.ngen_de, :, :] = de_obj.de_pop.copy()
            de_fit[: de_obj.ngen_de, :] = de_obj.de_fit.copy()
            de_pop_best[: de_obj.ngen_de, :] = de_obj.de_pop_best.copy()
            de_fit_best[: de_obj.ngen_de] = de_obj.de_fit_best.copy()
            last_iter_de = de_obj.iter_de
            loc_best = de_get_best_loc(de_fit_best, last_iter_de, de_maximize=de_maximize)
            de_lnP = de_fit_best[loc_best]
            de_par = de_obj.de_par
            if de_resume:
                de_run = True
        else:
            de_run = True
            anc.print_both(
                "pyDE file {} does not exist, cannot load pyDE run, starting from scratch.".format(
                    de_file
                ),
                output=of_run,
            )
    # end if de_to_emcee

    if de_run:

        if nthreads > 1:
            pool = Pool(nthreads)
        else:
            pool = None

        # create pyDE object
        de_evol = DiffEvol(
            fitness_func,
            de_bounds,
            npop_de,
            f=de_f,
            c=de_c,
            seed=conf.seed,
            maximize=de_maximize,
            vectorize=False,
            pool=pool,
        )

        if de_resume and de_file_check:
            de_evol._population = de_obj.de_pop[-1, :, :].copy()
            de_evol._fitness = de_obj.de_fit[-1, :].copy()
            de_evol._minidx = (
                np.argmax(de_evol._fitness)
                if de_maximize
                else np.argmin(de_evol._fitness)
            )
        else: # in case de_run = True only (or de_to_emcee = True, but file not existing)
            if de_file_check:
                anc.print_both(
                    "pyDE file {} exists: deleting it!".format(de_file), output=of_run
                )
                os.remove(de_file)

        start_iter_de = time.time()
        de_start = start_iter_de
        anc.print_both(" DE - START ", output=of_run)
        sys.stdout.flush()
        if last_iter_de + 1 < ngen_de:
            anc.print_both(
                " last_iter_de + 1 = {} < ngen_de = {}".format(
                    last_iter_de + 1, ngen_de
                ), output=of_run,
            )
            sys.stdout.flush()

            for iter_de_temp, res_de in enumerate(de_evol(ngen_de - last_iter_de)):
                # current iter_de_temp starts from 0 and ends to ngen_de-last_iter_de
                iter_de = last_iter_de + iter_de_temp
                # iter_de goes from last_iter_de to ngen_de
                mod_save = (iter_de + 1) % de_save
                check_save = mod_save == 0
                check_iter = iter_de + 1 >= ngen_de

                de_pop[iter_de, :, :] = de_evol.population.copy()
                de_fit[iter_de, :] = fit_type * de_evol._fitness.copy()
                de_pop_best[iter_de, :] = de_evol.minimum_location.copy()
                de_fit_best[iter_de] = fit_type * de_evol.minimum_value

                # print(iter_de, mod_save, check_save, check_iter, fit_type * de_evol.minimum_value)

                if iter_de > 0:
                    if check_save or check_iter:
                        anc.print_both(" ============= ", output=of_run)
                        anc.print_both(
                            " pyDE - iter = {} / {} ==> saving".format(
                                iter_de + 1, ngen_de
                            ),
                            output=of_run,
                        )
                        anc.print_both(
                            " last best fitness = {}".format(de_fit_best[iter_de]),
                            output=of_run,
                        )
                        sys.stdout.flush()

                        # SAVE DE SIMULATION IN de_run.hdf5 FILE
                        de_save_evolution(
                            de_file,
                            npop_de,
                            ngen_de,
                            iter_de,
                            1,
                            nfit,
                            ndata,
                            de_pop,
                            de_fit,
                            de_pop_best,
                            de_fit_best,
                            de_bounds,
                            sim.fitting_names,
                            de_maximize=de_maximize,
                        )
                        anc.print_both(
                            " Updated DE hdf5 file: {}".format(de_file),
                            output=of_run,
                        )
                        elapsed_de = time.time() - start_iter_de
                        (
                            elapsed_de_d,
                            elapsed_de_h,
                            elapsed_de_m,
                            elapsed_de_s,
                        ) = anc.computation_time(elapsed_de)
                        anc.print_both(
                            " pyDE: {:d} steps in {:2d} day {:02d} hour {:02d} min {:.2f} sec".format(
                                de_save,
                                int(elapsed_de_d),
                                int(elapsed_de_h),
                                int(elapsed_de_m),
                                elapsed_de_s,
                            ),
                            output=of_run,
                        )
                        sys.stdout.flush()
                        start_iter_de = time.time()
            # end for iter_de_temp, res_de

            # FORCE
            # SAVE DE SIMULATION IN de_run.hdf5 FILE
            de_save_evolution(
                de_file,
                npop_de,
                ngen_de,
                iter_de,
                1,
                nfit,
                ndata,
                de_pop,
                de_fit,
                de_pop_best,
                de_fit_best,
                de_bounds,
                sim.fitting_names,
                de_maximize=de_maximize,
            )
            anc.print_both(
                " Updated DE hdf5 file: {}".format(de_file), output=of_run
            )
        else:
            iter_de = ngen - 1
        anc.print_both("", output=of_run)
        anc.print_both(" completed DE", output=of_run)
        anc.print_both("", output=of_run)

        elapsed = time.time() - de_start
        elapsed_d, elapsed_h, elapsed_m, elapsed_s = anc.computation_time(elapsed)
        anc.print_both(" ", output=of_run)
        anc.print_both(
            " DE FINISHED in {0:2d} day {1:02d} hour {2:02d} min {3:.2f} sec - bye bye".format(
                int(elapsed_d), int(elapsed_h), int(elapsed_m), elapsed_s
            ),
            output=of_run,
        )
        sys.stdout.flush()

        if nthreads > 1:
            # pool.join()
            pool.close()
            pool.terminate()

        loc_best = de_get_best_loc(de_fit_best, iter_de, de_maximize=de_maximize)
        de_par = de_pop_best[loc_best, :].copy()
        de_lnP = de_fit_best[loc_best]
    
    return de_par, de_lnP, de_pop_best, de_fit_best, de_pop, de_fit


# ======================================================================
class DEEvolution:
    def __init__(self, run_folder, of_run=None):
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
        anc.print_both("Read de_run.hdf5 file with:", output=of_run)
        anc.print_both("npop_de = {}".format(self.npop_de), output=of_run)
        anc.print_both("ngen_de = {}".format(self.ngen_de), output=of_run)
        anc.print_both("iter_de = {}".format(self.iter_de), output=of_run)
        anc.print_both("iter_global = {}".format(self.iter_global), output=of_run)
        anc.print_both("nfit = {}".format(self.nfit), output=of_run)
        anc.print_both("ndata = {}".format(self.ndata), output=of_run)
        anc.print_both("", output=of_run)

        if self.iter_de == 0:
            anc.print_both("The file has been only initialised, wait a bit more ...", output=of_run)
            return
        if self.iter_de + 1 < self.ngen_de:
            self.ngen_de = self.iter_de + 1
            anc.print_both("updated ngen_de = {}".format(self.ngen_de), output=of_run)

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
        self.all_names = self.par_names.copy()

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

        self.de_pop_flat = pd.DataFrame(d)

        d = {}
        d["iter"] = np.arange(self.ngen_de)
        for i_p, p in enumerate(self.par_names):
            d[p] = self.de_pop_best[:, i_p].copy()
        d["de_fitness"] = self.de_fit_best.copy()
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
                self.all_names.append(new_name)

        return

    def convert_jitter(self):

        for i_name, name in enumerate(self.par_names):

            if "l2j" in name:
                new_name = name.replace("l2j", "jitter")
                self.de_pop_flat[new_name] = np.power(2, self.de_pop_flat[name])
                self.de_best_flat[new_name] = np.power(2, self.de_best_flat[name])
                self.all_names.append(new_name)

        return

    def recenter_angles(self):

        ll = "lambda"
        nl = len(ll)
        for i_name, name in enumerate(self.all_names):
            print(i_name, name)
            if (
                name[0] == "w"
                or name[0:2] == "mA"
                or name[0:2] == "lN"
                or name[0:nl] == ll
            ):
                new_name = "rec_{}".format(name)
                print(" ==> ", name, new_name)
                self.de_pop_flat[new_name], new_type = anc.recenter_angle_distribution(
                    self.de_pop_flat[name] % 360.0, type_out=True
                )
                # self.de_best_flat[new_name] = anc.recenter_angle_distribution(
                #     self.de_best_flat[name] % 360.0
                # )
                if new_type == "scale": 
                    self.de_best_flat[new_name] = anc.get_arctan_angle(self.de_best_flat[name])
                else:
                    self.de_best_flat[new_name] = self.de_best_flat[name] % 360.0
                self.all_names.append(new_name)

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
                    self.all_names.append(new_name)
                else:
                    new_name = name.replace("secosw", "ecc")
                    self.de_pop_flat[new_name] = square_all
                    self.de_best_flat[new_name] = square_best
                    self.all_names.append(new_name)
                new_name = "w{}".format(name.split("w")[1])
                self.de_pop_flat[new_name] = angle_all
                self.de_best_flat[new_name] = angle_best
                self.all_names.append(new_name)

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
