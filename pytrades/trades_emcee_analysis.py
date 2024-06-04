#!/usr/bin/env python
# -*- coding: utf-8 -*-

# from pytrades_lib import pytrades
import pytrades
import ancillary as anc

# import constants as cst  # local constants module
from gelman_rubin import compute_gr
from geweke import compute_geweke
from convergence import compute_convergence, full_statistics, log_probability_trace
from chains_summary_plot import plot_chains
from fitted_correlation_plot import plot_triangle as correlation_fitted
from physical_correlation_plot import plot_triangle as correlation_physical
import plot_oc as poc
import plot_rv as prv
import sys

# import argparse
# import time
import os
import numpy as np  # array
import h5py
import emcee

import matplotlib.pyplot as plt

# import numba
# numba.set_num_threads(1)
# numba.config.THREADING_LAYER = "tbb"

# numba.config.DISABLE_JIT = 1

# ==============================================================================
# custom modules
script_path = os.path.realpath(__file__)
module_path = os.path.abspath(os.path.join(os.path.dirname(script_path), "../"))
sys.path.append(module_path)

# ==============================================================================
anc.set_rcParams()
# ==============================================================================



def get_emcee_run(cli):

    # read emcee data
    emcee_file, _, _ = anc.get_emcee_file_and_best(cli.full_path, cli.temp_status)

    if cli.emcee_old_save:
        (
            fitting_names_original,
            _,
            chains,
            _,
            autocor_time,
            lnprobability,
            _,
            completed_steps,
        ) = anc.get_data(emcee_file, cli.temp_status)

        anc.print_both("temp_status     = {}".format(cli.temp_status))
        _, _, nruns, nburnin, nruns_sel = anc.get_emcee_parameters(
            chains, cli.temp_status, cli.nburnin, completed_steps
        )
        anc.print_both("completed_steps = {}".format(completed_steps))
        anc.print_both("nruns           = {}".format(nruns))
        anc.print_both("nburnin         = {}".format(nburnin))

        (
            chains_posterior,
            fitting_posterior,
            lnprobability_posterior,
            lnprob_posterior,
            chains_full_thinned,
            lnprobability_full_thinned,
        ) = anc.thin_the_chains(
            cli.use_thin,
            nburnin,
            nruns,
            chains,
            lnprobability,
            burnin_done=False,
            full_chains_thinned=True,
        )

    else:
        fitting_names_original = None
        backend = emcee.backends.HDFBackend(emcee_file)
        completed_steps = backend.iteration  # --> nruns
        nruns = completed_steps
        anc.print_both("... loading chains ...")
        sys.stdout.flush()
        chains = backend.get_chain()  # nruns x nwalkers x nfit
        anc.print_both("... loading posterior chains ...")
        sys.stdout.flush()
        chains_posterior = backend.get_chain(
            discard=cli.nburnin, flat=False, thin=cli.use_thin
        )
        anc.print_both("... loading thinned chains ...")
        sys.stdout.flush()
        chains_full_thinned = backend.get_chain(thin=cli.use_thin)
        anc.print_both("... loading thinned flat posteriors ...")
        sys.stdout.flush()
        fitting_posterior = backend.get_chain(
            discard=cli.nburnin, flat=True, thin=cli.use_thin
        )
        anc.print_both("... loading log_prob ...")
        sys.stdout.flush()
        lnprobability = backend.get_log_prob()  # nruns x nwalkers
        anc.print_both("... loading thinned log_prob ...")
        sys.stdout.flush()
        lnprobability_full_thinned = backend.get_log_prob(thin=cli.use_thin)
        anc.print_both("... loading thinned posterior log_prob ...")
        sys.stdout.flush()
        lnprobability_posterior = backend.get_log_prob(
            discard=cli.nburnin, flat=False, thin=cli.use_thin
        )
        anc.print_both("... loading thinned flat posterior log_prob ...")
        sys.stdout.flush()
        lnprob_posterior = backend.get_log_prob(
            discard=cli.nburnin, flat=True, thin=cli.use_thin
        )
        sys.stdout.flush()

    anc.print_both("... done")
    sys.stdout.flush()

    return (
        chains,
        chains_full_thinned,
        chains_posterior,
        fitting_posterior,
        lnprobability,
        lnprobability_full_thinned,
        lnprobability_posterior,
        lnprob_posterior,
        emcee_file,
        fitting_names_original,
    )


# ================================================================================================


class AnalysisTRADES:
    def __init__(self, cli):
        self.cli = cli
        np.random.seed(cli.seed)

        # self.log_folder = os.path.join(cli.full_path, 'logs')
        # os.makedirs(self.log_folder, exist_ok=True)
        # self.olog = open(os.path.join(self.log_folder, "log_confidence_intervals.txt"), 'w')

        self.summary_file = os.path.join(self.cli.full_path, "summary_parameters.hdf5")

        if cli.save_parameters:
            self.summary_h5f = h5py.File(self.summary_file, "w")
        else:
            self.summary_h5f = None

        # # init trades
        anc.print_both("Init trades ... ")
        
        sim = pytrades.TRADESfolder(cli.full_path, seed=cli.seed, m_type=cli.m_type)
        sim.init_trades()
        sys.stdout.flush()

        self.sim = sim
        self.star = sim.MR_star
        self.fitting_names = sim.fitting_names
        anc.print_both("Fitting parameters: {}".format(sim.fitting_names))

        self.thin_steps = cli.use_thin
        # get data from the hdf5 file
        anc.print_both("\nGet data from hdf5 file ...")
        (
            self.chains,
            self.chains_full_thinned,
            self.chains_posterior,
            self.fitting_posterior,
            self.lnprobability,
            self.lnprobability_full_thinned,
            self.lnprobability_posterior,
            self.lnprob_posterior,
            self.emcee_file,
            self.fitting_names_original,
        ) = get_emcee_run(cli)
        self.npost, _ = np.shape(self.fitting_posterior)

        sys.stdout.flush()
        # check if emcee run with old parameterization in (icos/sinlN) or (cosicos/sinlN)
        if self.fitting_names_original is not None:
            (
                self.fitting_names,
                self.fitting_posterior,
            ) = anc.update_parameterisation_posterior(
                self.fitting_names_original,
                self.fitting_posterior,
            )
            (_, self.chains_posterior,) = anc.update_parameterisation_chains(
                self.fitting_names_original,
                self.chains_posterior,
            )
            (_, self.chains,) = anc.update_parameterisation_chains(
                self.fitting_names_original,
                self.chains,
            )
            (_, self.chains_full_thinned,) = anc.update_parameterisation_chains(
                self.fitting_names_original,
                self.chains_full_thinned,
            )

        sys.stdout.flush()

        self.fitting_type = [None]*len(self.fitting_names)
        # # fix angles posterior
        for ifit, fitn in enumerate(self.fitting_names):
            if (
                fitn[0] == "w"
                or fitn[0:2] == "mA"
                or fitn[0:2] == "lN"
                or "lambda" in fitn
            ):
                anc.print_both("Fixing parameter {} with idx {}".format(fitn, ifit))
                # print("shape of chains_posterior = {}".format(np.shape(self.chains_posterior)))
                self.fitting_posterior[:, ifit], i_type = anc.recenter_angle_distribution(
                    self.fitting_posterior[:, ifit] % 360.0, debug=True, type_out=True
                )
                # self.chains_posterior[:, :, ifit] = anc.recenter_angle_distribution(
                #     self.chains_posterior[:, :, ifit] % 360.0, debug=True, type_out=False
                # )
                self.fitting_type[ifit] = i_type
                if "mod" in i_type:
                    self.chains[:, :, ifit] %= 360.0
                    self.chains_full_thinned[:, :, ifit] %= 360.0
                    self.chains_posterior[:, :, ifit] %= 360.0
                else:
                    self.chains[:, :, ifit] = anc.get_arctan_angle(self.chains[:, :, ifit])
                    self.chains_full_thinned[:, :, ifit] = anc.get_arctan_angle(self.chains_full_thinned[:, :, ifit])
                    self.chains_posterior[:, :, ifit] = anc.get_arctan_angle(self.chains_posterior[:, :, ifit])
        anc.print_both("####################################################################")
        for fitn, fitt in zip(self.fitting_names, self.fitting_type):
            print(fitn, fitt)
        anc.print_both("####################################################################")

        sys.stdout.flush()
        anc.print_both("Determine conversion factors ...")
        anc.print_both(
            "sim.m_factor = {} sim.mass_unit = {}".format(sim.m_factor, sim.mass_unit)
        )
        (
            self.mass_conv_factor,
            self.mass_unit,
            self.radius_conv_factor,
            self.radius_unit,
        ) = (sim.m_factor, sim.mass_unit, sim.r_factor, sim.radius)
        anc.print_both(
            "mass_conv_factor = {} mass_unit = {}".format(
                self.mass_conv_factor, self.mass_unit
            )
        )

        # full chains
        nrt, nw, _ = np.shape(self.chains_full_thinned)
        sys.stdout.flush()
        # STELLAR MASS
        self.mass_conv_factor_post = [1.0]*3
        if "m1" in sim.fitting_names:
            idx_m1 = list(sim.fitting_names).index("m1")
            # full chains
            self.stellar_mass_normal_chains = self.chains_full_thinned[:,:,idx_m1]
            
        else:
            # full chains
            self.stellar_mass_normal_chains = np.random.normal(
                loc=self.star[0, 0], scale=self.star[0, 1], size=(nrt, nw)
            )
        
        # posterior chains
        self.stellar_mass_normal_post_chains = self.stellar_mass_normal_chains[cli.nburnin:, :]
        # posterior flat
        self.stellar_mass_posterior = self.stellar_mass_normal_post_chains.reshape(self.npost)
        
        self.mass_conv_factor_post[0] = (
            self.mass_conv_factor / self.star[0, 0]
            ) * self.stellar_mass_normal_chains
        self.mass_conv_factor_post[1] = (
            self.mass_conv_factor / self.star[0, 0]
            ) * self.stellar_mass_normal_post_chains
        self.mass_conv_factor_post[2] = (
            self.mass_conv_factor / self.star[0, 0]
            ) * self.stellar_mass_posterior

        # STELLAR RADIUS
        self.radius_conv_factor_post = [1.0]*3
        if "R1" in sim.fitting_names:
            idx_R1 = list(sim.fitting_names).index("R1")
            # full chains
            self.stellar_radius_normal_chains = self.chains_full_thinned[:,:,idx_R1]
        else:
            # full chains
            self.stellar_radius_normal_chains = np.random.normal(
                loc=self.star[1, 0], scale=self.star[1, 1], size=(nrt, nw)
            )

        # posterior chains
        self.stellar_radius_normal_post_chains = self.stellar_radius_normal_chains[cli.nburnin:, :]
        # posterior 
        self.stellar_radius_posterior = self.stellar_radius_normal_post_chains.reshape(self.npost)

        self.radius_conv_factor_post[0] = (
            self.radius_conv_factor / self.star[1, 0]
        ) * self.stellar_radius_normal_chains
        self.radius_conv_factor_post[1] = (
            self.radius_conv_factor / self.star[1, 0]
        ) * self.stellar_radius_normal_post_chains
        self.radius_conv_factor_post[2] = (
            self.radius_conv_factor / self.star[1, 0]
        ) * self.stellar_radius_posterior

        sys.stdout.flush()
        anc.print_both("Get physical parameters/posterior ... ")
        (
            self.physical_names,
            self.physical_parameters,
            self.physical_posterior,
            self.phys_type,
            self.physical_chains,
            self.physical_chains_posterior,
        ) = anc.compute_physical_parameters(
            sim.n_bodies,
            sim.fitting_parameters,
            sim.fitting_names,
            sim.system_parameters,
            mass_conv_factor=self.mass_conv_factor,
            radius_conv_factor=self.radius_conv_factor,
            posterior_fit=self.fitting_posterior,
            chains_full_fit=self.chains_full_thinned,
            chains_posterior_fit=self.chains_posterior,
            mass_post_conv_factor=self.mass_conv_factor_post,
            radius_post_conv_factor=self.radius_conv_factor_post,
        )

        sys.stdout.flush()
        anc.print_both("Get units ...")
        self.fitting_units = anc.get_units(self.fitting_names, self.mass_unit)
        self.physical_units = anc.get_units(self.physical_names, self.mass_unit)

        anc.print_both("Preparation of the analysis of TRADES simulation done.")

    # ---------------------------------
    def save_posterior(self):

        self.posterior_file = os.path.join(self.cli.full_path, "posterior.hdf5")
        with h5py.File(self.posterior_file, "w") as p_h5f:
            p_h5f.create_dataset(
                "posterior", data=self.fitting_posterior, dtype=np.float64
            )
            p_h5f.create_dataset(
                "loglikelihood", data=self.lnprob_posterior, dtype=np.float64
            )
            p_h5f["posterior"].attrs["nfit"] = self.sim.nfit
            p_h5f["posterior"].attrs["nposterior"] = self.npost

            p_h5f.create_dataset(
                "fitting_names", data=anc.encode_list(self.fitting_names), dtype="S15"
            )

            p_h5f.create_dataset(
                "posterior_physical", data=self.physical_posterior, dtype=np.float64
            )
            p_h5f.create_dataset(
                "physical_names", data=anc.encode_list(self.physical_names), dtype="S15"
            )

            p_h5f.create_dataset(
                "stellar_mass", data=self.stellar_mass_posterior, dtype=np.float64
            )
            p_h5f.create_dataset(
                "stellar_radius", data=self.stellar_radius_posterior, dtype=np.float64
            )

        anc.print_both(" Saved posterior file: {}".format(self.posterior_file))

    # ---------------------------------
    def save_observables_from_samples(self, samples_fit_par, smp_h5, olog=None):

        anc.print_both("=" * 40, olog)
        anc.print_both("= BEGIN = Getting Transit Times and RV from samples ...", olog)

        # n_samples, nfit = np.shape(samples_fit_par)
        n_samples, _ = np.shape(samples_fit_par)

        for ismp in range(0, n_samples):

            smp_name = "{0:04d}".format(ismp)

            self.sim.save_models_from_parameters(
                samples_fit_par[ismp],
                smp_h5,
                smp_name,
            )

        anc.print_both("= END = Getting Transit Times and RV from samples ...", olog)
        anc.print_both("=" * 40, olog)
        sys.stdout.flush()

        return

    def get_and_save_samples(self):

        if self.cli.n_samples > 0:
            anc.print_both(
                " Selecting {:d} samples from the posterior ...".format(
                    self.cli.n_samples
                )
            )
            sys.stdout.flush()
            samples_fit_par = anc.take_n_samples(
                self.fitting_posterior,
                lnprob=self.lnprob_posterior,
                post_ci=self.ci_fitted[0:2, :],
                n_samples=self.cli.n_samples,
            )
            anc.print_both(" Running TRADES and computing the T0s and RVs ...")
            samples_file = os.path.join(self.cli.full_path, "samples_ttra_rv.hdf5")
            anc.print_both(" Saving into {:s}".format(samples_file))
            with h5py.File(samples_file, "w") as smp_h5:
                self.save_observables_from_samples(
                    samples_fit_par,
                    smp_h5,
                )
            anc.print_both(" ... done")
            sys.stdout.flush()

        return

    # ---------------------------------
    def run_and_save_from_parameters(
        self,
        fit_par,
        phy_par,
        id_sim,
        sim_name,
        par_description,
        ci_fit,
        ci_phy,
    ):

        full_sim_name = "{:04d}_sim_{:s}".format(id_sim, sim_name)
        out_folder = os.path.join(self.cli.full_path, full_sim_name, "")
        anc.print_both(out_folder)
        os.makedirs(out_folder, exist_ok=True)

        # compute sigma fit/phy!!
        anc.print_both("#############################")
        anc.print_both("CHECK SCALE/MOD FITTING ANGULAR PARAMETER")
        fit_temp = np.copy(fit_par)
        for ifit, itype in enumerate(self.fitting_type):
            if itype is not None:
                if itype == "mod":
                    fit_temp[ifit] = fit_temp[ifit]%360.0
                else:
                    # cc = np.cos(fit_temp[ifit] * cst.deg2rad)
                    # ss = np.sin(fit_temp[ifit] * cst.deg2rad)
                    # fit_temp[ifit] = np.arctan2(ss, cc) * cst.rad2deg
                    fit_temp[ifit] = anc.get_arctan_angle(fit_par[ifit])
                anc.print_both("{0:s} = {1:.4f} ==> {2:.4f}".format(
                    self.fitting_names[ifit], fit_par[ifit], fit_temp[ifit]
                ))
        anc.print_both("#############################")

        sigma_hdi_fit = anc.hdi_to_sigma(fit_temp, ci_fit)
        mad_fit, rms_fit = anc.posterior_to_rms_mad(self.fitting_posterior, fit_temp)
        sigma_fit = np.column_stack((rms_fit, mad_fit, sigma_hdi_fit.T)).T

        sigma_hdi_phy = anc.hdi_to_sigma(phy_par, ci_phy)
        mad_phy, rms_phy = anc.posterior_to_rms_mad(self.physical_posterior, phy_par)
        sigma_phy = np.column_stack((rms_phy, mad_phy, sigma_hdi_phy.T)).T

        (
            chi_square,
            reduced_chi_square,
            lgllhd,
            lnprior,
            ln_const,
            bic,
            check,
        ) = self.sim.run_and_write_summary_files(
            full_sim_name,
            fit_par,
            phy_par,
            id_sim=id_sim,
            sigma_fit=sigma_fit,
            sigma_phy=sigma_phy,
            par_description=par_description,
            return_stats=True,
        )

        s_h5f = self.summary_h5f
        if s_h5f is not None:
            gr = s_h5f.create_group("parameters/{:s}".format(full_sim_name))
            gr.attrs["info"] = par_description
            gr.attrs["chi_square"] = chi_square
            gr.attrs["reduced_chi_square"] = reduced_chi_square
            gr.attrs["lgllhd"] = lgllhd
            gr.attrs["lnprior"] = lnprior
            gr.attrs["ln_const"] = ln_const
            gr.attrs["bic"] = bic
            gr.attrs["check"] = check

            gr.create_dataset(
                "fitted/parameters",
                data=fit_par,
                dtype=np.float64,
                compression="gzip",
            )
            gr.create_dataset(
                "fitted/names",
                data=anc.encode_list(self.fitting_names),
                dtype="S15",
                compression="gzip",
            )
            gr.create_dataset(
                "fitted/units",
                data=anc.encode_list(self.fitting_units),
                dtype="S15",
                compression="gzip",
            )
            gr.create_dataset(
                "fitted/sigma",
                data=sigma_fit.T,
                dtype=np.float64,
                compression="gzip",
            )
            gr["fitted/sigma"].attrs["percentiles"] = anc.percentile_val

            gr.create_dataset(
                "physical/parameters",
                data=phy_par,
                dtype=np.float64,
                compression="gzip",
            )
            gr.create_dataset(
                "physical/names",
                data=anc.encode_list(self.physical_names),
                dtype="S15",
                compression="gzip",
            )
            gr.create_dataset(
                "physical/units",
                data=anc.encode_list(self.physical_units),
                dtype="S15",
                compression="gzip",
            )
            gr.create_dataset(
                "physical/sigma",
                data=sigma_phy.T,
                dtype=np.float64,
                compression="gzip",
            )
            gr["physical/sigma"].attrs["percentiles"] = anc.percentile_val

        return

    def save_original_parameters(self):
        # ==============================================================================
        # ORIGINAL FITTING PARAMETERS ID == 0
        # ==============================================================================
        id_sim = 0
        sim_name = "initial"
        par_desc = "initial parameters"
        anc.print_both("\nPARAMETERS: {} ==> {}".format(par_desc, id_sim))
        sys.stdout.flush()

        self.run_and_save_from_parameters(
            self.sim.fitting_parameters,
            self.sim.physical_parameters,
            id_sim,
            sim_name,
            par_desc,
            self.ci_fitted,
            self.ci_physical,
        )

        return

    def save_pso_parameters(self):

        # ==============================================================================
        # IF PSO BEST FILE AVAILABLE PARAMETERS -> id 0004
        # ==============================================================================
        pso_file = os.path.join(  # or pso_pattern
            self.cli.full_path,
            # "*_bestpso_*par.out"
            "pso_run.hdf5",
        )
        pso = os.path.exists(pso_file) and os.path.isfile(pso_file)
        if pso:
            id_sim = 4
            sim_name = "pso"
            par_desc = "best-fit parameters from PSO"
            anc.print_both("\nPARAMETERS: {} ==> {}".format(par_desc, id_sim))
            sys.stdout.flush()
            with h5py.File(pso_file, "r") as fpso:
                pso_fit = fpso["pso_parameters"][...].copy()
            _, pso_phy, _, _ = anc.compute_physical_parameters(
                self.sim.n_bodies,
                pso_fit,
                self.fitting_names,
                self.sim.system_parameters,
                mass_conv_factor=self.mass_conv_factor,
                radius_conv_factor=self.radius_conv_factor,
            )
            pso_phy = anc.search_scale_parameter(pso_phy, self.phys_type)
            self.run_and_save_from_parameters(
                pso_fit,
                pso_phy,
                id_sim,
                sim_name,
                par_desc,
                self.ci_fitted,
                self.ci_physical,
            )
        return

    def save_de_parameters(self):
        # ==============================================================================
        # IF DE BEST FILE AVAILABLE PARAMETERS -> id 0006
        # ==============================================================================
        de_file = os.path.join(self.cli.full_path, "de_run.hdf5")
        de = os.path.exists(de_file) and os.path.isfile(de_file)
        if de:
            id_sim = 6
            sim_name = "de"
            par_desc = "best-fit parameters from DE"
            anc.print_both("\nPARAMETERS: {} ==> {}".format(par_desc, id_sim))
            sys.stdout.flush()
            with h5py.File(de_file, "r") as fde:
                de_fit = fde["de_parameters"][...].copy()
            _, de_phy, _, _, _, _ = anc.compute_physical_parameters(
                self.sim.n_bodies,
                de_fit,
                self.fitting_names,
                self.sim.system_parameters,
                mass_conv_factor=self.mass_conv_factor,
                radius_conv_factor=self.radius_conv_factor,
            )
            de_phy = anc.search_scale_parameter(de_phy, self.phys_type)
            self.run_and_save_from_parameters(
                de_fit,
                de_phy,
                id_sim,
                sim_name,
                par_desc,
                self.ci_fitted,
                self.ci_physical,
            )
        return

    def save_median_parameters(self):
        # ==============================================================================
        # MEDIAN of fitted and phyiscal posteriorior ID == 1050
        # ==============================================================================
        id_sim = 1050
        sim_name = "median"
        par_desc = "median of posterior and median of physical posterior"
        anc.print_both("\nPARAMETERS: {} ==> {}".format(par_desc, id_sim))
        median_fit = np.median(self.fitting_posterior, axis=0)
        median_phys = np.median(self.physical_posterior, axis=0)
        sys.stdout.flush()

        self.run_and_save_from_parameters(
            median_fit,
            median_phys,
            id_sim,
            sim_name,
            par_desc,
            self.ci_fitted,
            self.ci_physical,
        )

        # ==============================================================================
        id_sim = 1051
        sim_name = "median_to_physical"
        par_desc = "median of posterior and converted to physical"
        anc.print_both("\nPARAMETERS: {} ==> {}".format(par_desc, id_sim))
        _, median_phys, _, _, _, _ = anc.compute_physical_parameters(
            self.sim.n_bodies,
            median_fit,
            self.fitting_names,
            self.sim.system_parameters,
            mass_conv_factor=self.mass_conv_factor,
            radius_conv_factor=self.radius_conv_factor,
        )
        median_phys = anc.search_scale_parameter(median_phys, self.phys_type)
        sys.stdout.flush()
        self.run_and_save_from_parameters(
            median_fit,
            median_phys,
            id_sim,
            sim_name,
            par_desc,
            self.ci_fitted,
            self.ci_physical,
        )

        return

    def save_map_parameters(self):
        # ==============================================================================
        # MAP == maximum a posteriori probability ID == 2050
        # ==============================================================================
        id_sim = 2050
        sim_name = "map"
        par_desc = "map of posterior"
        anc.print_both("\nPARAMETERS: {} ==> {}".format(par_desc, id_sim))

        map_idx = np.argmax(self.lnprob_posterior)
        map_fit = self.fitting_posterior[map_idx, :].copy()
        map_phy = self.physical_posterior[map_idx, :].copy()

        # _, map_phy, _, _, _, _ = anc.compute_physical_parameters(
        #     self.sim.n_bodies,
        #     map_fit,
        #     self.fitting_names,
        #     self.sim.system_parameters,
        #     mass_conv_factor=self.mass_conv_factor,
        #     radius_conv_factor=self.radius_conv_factor,
        # )
        # map_phy = anc.search_scale_parameter(map_phy, self.phys_type)
        sys.stdout.flush()
        self.run_and_save_from_parameters(
            map_fit,
            map_phy,
            id_sim,
            sim_name,
            par_desc,
            self.ci_fitted,
            self.ci_physical,
        )

        # ==============================================================================
        id_sim = 2051
        sim_name = "map_to_physical"
        par_desc = "map of posterior converted to physical"
        anc.print_both("\nPARAMETERS: {} ==> {}".format(par_desc, id_sim))
        _, map_phy, _, _, _, _ = anc.compute_physical_parameters(
            self.sim.n_bodies,
            map_fit,
            self.fitting_names,
            self.sim.system_parameters,
            mass_conv_factor=self.mass_conv_factor,
            radius_conv_factor=self.radius_conv_factor,
        )
        map_phy = anc.search_scale_parameter(map_phy, self.phys_type)
        sys.stdout.flush()
        self.run_and_save_from_parameters(
            map_fit,
            map_phy,
            id_sim,
            sim_name,
            par_desc,
            self.ci_fitted,
            self.ci_physical,
        )

        return

    def save_mode_parameters(self):
        # ==============================================================================
        # MODE-LIKE == mode parameters as 3 x median - 2 x mean ID == 3050
        # ==============================================================================
        id_sim = 3050
        sim_name = "mode"
        par_desc = "mode of fitted and physical posterior"
        anc.print_both("\nPARAMETERS: {} ==> {}".format(par_desc, id_sim))

        median = np.percentile(
            self.fitting_posterior,
            50.0,
            axis=0,
            #
        )
        mean = np.mean(self.fitting_posterior, axis=0)
        mode_fit = 3 * median - 2 * mean

        if np.any(np.isnan(mode_fit)):
            anc.print_both(
                "MODE parameters: Some values are Nan, skip the mode parameters"
            )
        else:
            neg_values = False
            neg_index = None

            for p_n in self.fitting_names:
                if ("Ms" in p_n) or ("P" in p_n):
                    v = mode_fit[list(self.fitting_names).index(p_n)]
                    if v <= 0.0:
                        neg_values = True
                        neg_index = self.fitting_names.index(p_n)
                elif (
                    ("w" in p_n) or ("mA" in p_n) or ("lambda" in p_n) or ("lN" in p_n)
                ):
                    mode_fit[list(self.fitting_names).index(p_n)] %= 360.0

            if neg_values:
                anc.print_both(
                    "MODE parameters: Parameter {:s} <= 0.0!! skip the mode parameters".format(
                        self.fitting_names[neg_index]
                    )
                )
            else:

                median = np.percentile(
                    self.physical_posterior,
                    50.0,
                    axis=0,
                    #
                )
                mean = np.mean(self.physical_posterior, axis=0)
                mode_phy = 3 * median - 2 * mean

                for p_n in self.physical_names:
                    if (
                        ("w" in p_n)
                        or ("mA" in p_n)
                        or ("lambda" in p_n)
                        or ("lN" in p_n)
                    ):
                        mode_phy[list(self.physical_names).index(p_n)] %= 360.0

                self.run_and_save_from_parameters(
                    mode_fit,
                    mode_phy,
                    id_sim,
                    sim_name,
                    par_desc,
                    self.ci_fitted,
                    self.ci_physical,
                )

        return

    def save_from_file_parameters(self):
        # ==============================================================================
        # SELECT AD HOC PARAMETERS FROM FILE ID == 777
        # ==============================================================================
        if self.cli.from_file is not None:
            id_sim = 777
            sim_name = "from_file"
            par_desc = "parameters from file"
            anc.print_both("\nPARAMETERS: {} ==> {}".format(par_desc, id_sim))
            sys.stdout.flush()
            _, adhoc_fit = anc.read_fitted_file(self.cli.from_file)

            _, adhoc_phy, _, _, _, _ = anc.compute_physical_parameters(
                self.sim.n_bodies,
                adhoc_fit,
                self.fitting_names,
                self.sim.system_parameters,
                mass_conv_factor=self.mass_conv_factor,
                radius_conv_factor=self.radius_conv_factor,
            )
            adhoc_phy = anc.search_scale_parameter(adhoc_phy, self.phys_type)

            self.run_and_save_from_parameters(
                adhoc_fit,
                adhoc_phy,
                id_sim,
                sim_name,
                par_desc,
                self.ci_fitted,
                self.ci_physical,
            )

        return

    def save_map_hdi_parameters(self):
        # ==============================================================================
        # MAP_HDI == maximum log-probability=logL + logPrior of posterior within HDI ID == 668
        # ==============================================================================
        id_sim = 668
        sim_name = "map_hdi"
        par_desc = "map_hdi of posterior"
        anc.print_both("\nPARAMETERS: {} ==> {}".format(par_desc, id_sim))

        try:
            map_hdi_fit, _, map_hdi_idx = anc.select_maxlglhd_with_hdi(
                self.fitting_posterior,
                self.ci_fitted.T,
                self.lnprob_posterior,
                return_idmax=True,
            )
            map_hdi_phy = self.physical_posterior[map_hdi_idx, :]
            sys.stdout.flush()
            self.run_and_save_from_parameters(
                map_hdi_fit,
                map_hdi_phy,
                id_sim,
                sim_name,
                par_desc,
                self.ci_fitted,
                self.ci_physical,
            )

            # ==============================================================================
            id_sim = 669
            sim_name = "map_hdi_to_physical"
            par_desc = "map of posterior within HDI converted to physical"
            anc.print_both("\nPARAMETERS: {} ==> {}".format(par_desc, id_sim))
            _, map_hdi_phy, _, _, _, _ = anc.compute_physical_parameters(
                self.sim.n_bodies,
                map_hdi_fit,
                self.fitting_names,
                self.sim.system_parameters,
                mass_conv_factor=self.mass_conv_factor,
                radius_conv_factor=self.radius_conv_factor,
            )
            map_hdi_phy = anc.search_scale_parameter(map_hdi_phy, self.phys_type)
            sys.stdout.flush()
            self.run_and_save_from_parameters(
                map_hdi_fit,
                map_hdi_phy,
                id_sim,
                sim_name,
                par_desc,
                self.ci_fitted,
                self.ci_physical,
            )

        except:
            par_desc += " !! MAP WITHIN HDI NOT FOUND --> USING MAP"
            full_sim_name = "{:04d}_sim_{:s}".format(id_sim, sim_name)
            out_folder = os.path.join(os.path.join(self.cli.full_path, full_sim_name), "")
            os.makedirs(out_folder, exist_ok=True)
            warn_file = os.path.join(out_folder, "WARNING.txt")
            anc.print_both("WRITING WARNING FILE: {:s}".format(warn_file))
            with open(os.path.join(out_folder, "WARNING.txt"), "w") as warn_o:
                anc.print_both("*******\nWARNING: {}\n*******".format(par_desc), output=warn_o)
            anc.print_both("")
            anc.print_both("# " + "=" * 40)
            anc.print_both("")
        
        return

    # ---------------------------------
    def summary_parameters(self):

        anc.print_both("Computing HDI/CI")
        # hdi_ci, mode_parameters = anc.compute_hdi_full(fitting_posterior, mode_output=True)
        hdi_ci = anc.compute_hdi_full(self.fitting_posterior)
        ci_fitted = hdi_ci.T
        anc.print_both(
            " shape: hdi_ci = {} ci_fitted = {}".format(
                np.shape(hdi_ci), np.shape(ci_fitted)
            )
        )
        self.ci_fitted = ci_fitted
        # hdi_ci: nfit x nci
        # ci_fitted: nci x nfit
        # nci -> -1sigma(0) +1sigma(1) -2sigma(2) +2sigma(3) -3sigma(4) +3sigma(5)

        hdi_ci_phys = anc.compute_hdi_full(self.physical_posterior)
        self.ci_physical = hdi_ci_phys.T
        anc.print_both(
            " shape: hdi_ci_phys = {} ci_physical = {}".format(
                np.shape(hdi_ci_phys), np.shape(self.ci_physical)
            )
        )

        if self.cli.save_parameters:
            anc.print_both("Creates and populates summary_parameters.hdf5 file\n")
            sys.stdout.flush()

            s_h5f = self.summary_h5f
            if s_h5f is not None:
                s_h5f.create_dataset(
                    "confidence_intervals/fitted/ci",
                    data=ci_fitted.T,
                    dtype=np.float64,
                    compression="gzip",
                )
                s_h5f.create_dataset(
                    "confidence_intervals/fitted/names",
                    data=anc.encode_list(self.fitting_names),
                    dtype="S15",
                    compression="gzip",
                )
                s_h5f.create_dataset(
                    "confidence_intervals/fitted/units",
                    data=anc.encode_list(self.fitting_units),
                    dtype="S15",
                    compression="gzip",
                )

                s_h5f["confidence_intervals/fitted"].attrs["nfit"] = self.sim.nfit
                s_h5f["confidence_intervals/fitted"].attrs["nfree"] = self.sim.nfree
                s_h5f["confidence_intervals/fitted"].attrs["ndata"] = self.sim.ndata
                s_h5f["confidence_intervals/fitted"].attrs["dof"] = self.sim.dof

                s_h5f.create_dataset(
                    "confidence_intervals/physical/ci",
                    data=self.ci_physical.T,
                    dtype=np.float64,
                    compression="gzip",
                )
                s_h5f.create_dataset(
                    "confidence_intervals/physical/names",
                    data=anc.encode_list(self.physical_names),
                    dtype="S15",
                    compression="gzip",
                )
                s_h5f.create_dataset(
                    "confidence_intervals/physical/units",
                    data=anc.encode_list(self.physical_units),
                    dtype="S15",
                    compression="gzip",
                )

                self.save_original_parameters()
                self.save_pso_parameters()
                self.save_de_parameters()
                self.save_median_parameters()
                self.save_map_parameters()
                self.save_map_hdi_parameters()
                self.save_mode_parameters()
                self.save_from_file_parameters()

        self.get_and_save_samples()

        if self.summary_h5f is not None:
            self.summary_h5f.close()

        sys.stdout.flush()

        return


# ==============================================================================


def run_analysis(cli):

    anc.print_both("")
    anc.print_both(" TRADES analysis")
    anc.print_both("")

    analysis = AnalysisTRADES(cli)
    conf_run = anc.ConfigurationRun(cli.yaml_file)
    # os.environ["OMP_NUM_THREADS"] = str(cli.nthreads)
    sys.stdout.flush()

    # analysis.fit_to_physical()
    anc.print_both("ndata  = {}".format(analysis.sim.ndata))
    anc.print_both(" nfit  = {}".format(analysis.sim.nfit))
    anc.print_both("  dof  = {}".format(analysis.sim.dof))
    anc.print_both("star = {}".format(analysis.sim.MR_star))

    anc.print_both(
        "shape of chains_posterior = {}".format(np.shape(analysis.chains_posterior))
    )
    anc.print_both(
        "shape of posterior fit = {}".format(np.shape(analysis.fitting_posterior))
    )
    anc.print_both(
        "shape of posterior phy = {}".format(np.shape(analysis.physical_posterior))
    )
    anc.print_both("npost = {}".format(analysis.npost))

    logs_folder = os.path.join(cli.full_path, "logs")
    os.makedirs(logs_folder, exist_ok=True)
    plots_folder = os.path.join(cli.full_path, "plots")
    os.makedirs(plots_folder, exist_ok=True)

    # ===== Save Posterior ===== #
    if cli.save_posterior:
        analysis.save_posterior()

    # ===== CI/HDI and sample parameters ===== #
    analysis.summary_parameters()

    anc.print_both("Completed analysis of TRADES simulation.")

    # ===== LOG-PROB TRACE plot ===== #
    log_probability_trace(
        analysis.lnprobability_full_thinned, 
        analysis.lnprob_posterior, 
        plots_folder, 
        n_burn=cli.nburnin,
        n_thin=conf_run.thin_by,
        show_plot=False, 
        figsize=(6, 6)
    )

    # ===== CONVERGENCE ===== #
    if cli.gelman_rubin or cli.geweke or cli.chain:
        anc.print_both("\nComputing converging stats ...")
        par_file = os.path.join(cli.full_path, "summary_parameters.hdf5")
        stats_file = os.path.join(logs_folder, "convergence_stats.logs")
        overplot = anc.set_overplot(cli.overplot)
        with open(stats_file, 'w') as olog:
            with h5py.File(par_file, "r") as s_h5f:
                l = "fitted"
                if overplot is not None:
                    sim_id_str = "{}".format(overplot)
                    overp_par = s_h5f["parameters/{:s}/{:s}/parameters".format(sim_id_str, l)][
                        ...
                    ]
                else:
                    overp_par = analysis.fitting_posterior[np.argmax(analysis.lnprob_posterior), :]
                exp_acf_fit, exp_steps_fit = full_statistics(
                    analysis.chains_full_thinned,
                    analysis.fitting_posterior,
                    analysis.fitting_names,
                    overp_par,
                    analysis.lnprob_posterior,
                    plots_folder,
                    olog=olog,
                    ilast=0,
                    n_burn=cli.nburnin,
                    n_thin=conf_run.thin_by,
                    show_plot=False,
                    figsize=(6, 6),
                )
                
                l = "physical"
                if overplot is not None:
                    sim_id_str = "{}".format(overplot)
                    overp_par = s_h5f["parameters/{:s}/{:s}/parameters".format(sim_id_str, l)][
                        ...
                    ]
                else:
                    overp_par = analysis.physical_posterior[np.argmax(analysis.lnprob_posterior), :]
                exp_acf_phy, exp_steps_phy = full_statistics(
                    analysis.physical_chains,
                    analysis.physical_posterior,
                    analysis.physical_names,
                    overp_par,
                    analysis.lnprob_posterior,
                    plots_folder,
                    olog=olog,
                    ilast=analysis.sim.nfit,
                    n_burn=cli.nburnin,
                    n_thin=conf_run.thin_by,
                    show_plot=False,
                    figsize=(6, 6),
                )
            anc.print_both("", output=olog)
            anc.print_both(
                "All expected steps   for each parameter needed to reach full convergence:\n{}".format(
                    exp_steps_fit
                ),
                output=olog,
            )
            anc.print_both(
                "All expected ACF len for each parameter needed to reach full convergence:\n{}".format(
                    exp_acf_fit
                ),
                output=olog,
            )
            imax_acf = np.argmax(exp_acf_fit)
            anc.print_both(
                "MAX ACF = {} ==> needed chains of {} steps\n".format(
                    exp_acf_fit[imax_acf], exp_steps_fit[imax_acf]
                ),
                output=olog,
            )

    # ===== Fitted Correlation Plot ===== #
    if cli.correlation_fitted:
        anc.print_both("\nPlotting Fitted Correlation Diagram ... ")
        correlation_fitted(
            cli,
            logs_folder,
            plots_folder,
            analysis.fitting_posterior,
            analysis.fitting_names,
            analysis.sim.fitting_minmax,
        )

    # ===== Physical Correlation Plot ===== #
    if cli.correlation_physical:
        anc.print_both("\nPlotting Physical Correlation Diagram ... ")
        correlation_physical(
            cli,
            logs_folder,
            plots_folder,
            analysis.physical_posterior,
            analysis.physical_names,
        )

    # ===== OC ===== #
    if cli.plot_oc:
        samples_oc = None
        if len(cli.ocs) > 0:
            for coc in cli.ocs:
                if coc.tscale is None or str(coc.tscale).lower() == "none":
                    coc.tscale = analysis.sim.tepoch
                print(coc.samples_file)
                if samples_oc is None and coc.samples_file is not None:
                    samples_oc = poc.read_samples(samples_file=coc.samples_file)
                anc.print_both("\n ================ ")
                anc.print_both("OC of {}".format(os.path.basename(coc.full_path)))
                file_list = poc.get_simT0_file_list(coc)
                for idbody, flist in file_list.items():
                    if (cli.idplanet_name is None) or (len(cli.idplanet_name) == 0):
                        pl_name = None
                    else:
                        pl_name = cli.idplanet_name[idbody]
                    fig = poc.plot_oc_T41(
                        coc,
                        flist,
                        planet_name=pl_name,
                        samples=samples_oc,
                        save_plot=True,
                        show_plot=False,
                    )
                    plt.close(fig)
                    anc.print_both(flist)
    anc.print_both("")

    # ===== RV ===== #
    if cli.plot_rv:
        samples_rv = None
        if len(cli.rvs) > 0:
            for crv in cli.rvs:
                if crv.tscale is None or str(crv.tscale).lower() == "none":
                    crv.tscale = analysis.sim.tepoch
                if samples_rv is None and crv.samples_file is not None:
                    samples_rv = prv.read_samples(crv, samples_file=crv.samples_file)
                anc.print_both("\n ================ ")
                anc.print_both("RV of {}".format(os.path.basename(crv.full_path)))
                try:
                    fig = prv.plot_rv(
                        crv, samples=samples_rv, save_plot=True, show_plot=False
                    )
                    plt.close(fig)
                except:
                    print("Issues plotting RV for {}".format(crv.full_path))
    anc.print_both("")

    anc.print_both("")
    analysis.sim.reset()

    return


# ==============================================================================
def analysis_from_file(yml_file):

    cli = anc.ConfigurationAnalysis(yml_file)
    run_analysis(cli)

    return


# ==============================================================================
# main
def main():

    yml_file = anc.get_input_file()
    analysis_from_file(yml_file)

    return


# ==============================================================================
# ==============================================================================


if __name__ == "__main__":
    main()
