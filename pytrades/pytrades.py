# pytrades module, it imports the fortran-python library and creates a trades object
import numpy as np
import os
import h5py
import sys

# from pytrades_lib import f90trades as f90trades
from pytrades_lib import f90trades
import ancillary as anc

# ============================


class TRADES:
    def __init__(self, input_path, sub_folder="", nthreads=1, seed=42, m_type="e"):

        self.work_path = os.path.join(input_path, "")
        self.sub_folder = sub_folder
        self.full_path = os.path.join(self.work_path, self.sub_folder, "")
        self.nthreads = nthreads
        self.seed = seed
        self.m_type = m_type

        return

    def init_trades_from_path(self, work_path, sub_path):

        print("Initialising TRADES by f90 modules ...")
        sys.stdout.flush()

        self.reset()

        w_path = os.path.join(work_path, "")

        # INITIALISE TRADES WITH SUBROUTINE WITHIN TRADES_LIB -> PARAMETER NAMES, MINMAX, INTEGRATION ARGS, READ DATA ...
        f90trades.initialize_trades(
            w_path, sub_path, self.nthreads
        )

        print("... done")
        sys.stdout.flush()

        self.full_path = f90trades.get_path().decode("utf-8").strip()
        self.full_path = os.path.join(os.path.abspath(self.full_path), "")

        self.n_bodies = f90trades.n_bodies  # NUMBER OF TOTAL BODIES OF THE SYSTEM
        self.n_planets = self.n_bodies - 1  # NUMBER OF PLANETS IN THE SYSTEM
        self.ndata = f90trades.ndata  # TOTAL NUMBER OF DATA AVAILABLE
        self.nkel = f90trades.nkel # NUMBER OF FITTING KEPLERIAN ELEMENTS
        self.nfit = f90trades.nfit  # NUMBER OF PARAMETERS TO FIT
        self.nfree = f90trades.nfree  # NUMBER OF FREE PARAMETERS (ie nrvset)
        self.dof = f90trades.dof  # NUMBER OF DEGREES OF FREEDOM = NDATA - NFIT
        self.inv_dof = f90trades.inv_dof

        self.tepoch = f90trades.tepoch

        # fitting parameters
        str_len = f90trades.str_len
        temp_names = f90trades.get_parameter_names(self.nfit, str_len)
        # names
        self.fitting_names = anc.convert_fortran_charray2python_strararray(temp_names)
        # fitting values
        self.fitting_parameters = f90trades.fitting_parameters.copy()
        # boundaries
        self.fitting_minmax = f90trades.parameters_minmax.copy()

        # all system_parameters
        self.system_parameters = f90trades.system_parameters.copy()

        # Stellar Mass and Radius
        self.MR_star = f90trades.mr_star.copy()
        # Mass conversion factor for output
        # m_factor_0, self.mass_unit = anc.mass_type_factor(
        #     Ms=1.0, mtype=self.m_type, mscale=False
        # )
        m_factor_0, self.mass_unit, r_factor_0, self.radius_unit = anc.mass_radius_type_factor(mtype=self.m_type)
        self.m_factor = m_factor_0 * self.MR_star[0, 0]
        self.m_factor_post = self.m_factor # temporary
        self.r_factor = r_factor_0 * self.MR_star[1, 0]
        self.r_factor_post = self.r_factor # temporary

        # units of the fitted parameters
        self.fitting_units = anc.get_units(self.fitting_names, self.mass_unit)

        # get information on how to decompose fitted parameters in physical/derived ones
        # (
        #     _,
        #     _,
        #     self.bodies_file,
        #     self.id_fit,
        #     self.id_all,
        #     self.nfit_list,
        #     self.cols_list,
        #     self.case,
        # ) = anc.get_fitted(self.work_path)

        # set kepler elements, divided by names ...
        self.mass = np.zeros((self.n_bodies))
        (
            self.radius,
            self.period,
            self.sma,
            self.ecc,
            self.argp,
            self.meanA,
            self.inc,
            self.longN,
        ) = (
            np.zeros_like(self.mass),
            np.zeros_like(self.mass),
            np.zeros_like(self.mass),
            np.zeros_like(self.mass),
            np.zeros_like(self.mass),
            np.zeros_like(self.mass),
            np.zeros_like(self.mass),
            np.zeros_like(self.mass),
        )
        (
            self.mass,
            self.radius,
            self.period,
            self.sma,
            self.ecc,
            self.argp,
            self.meanA,
            self.inc,
            self.longN,
            self.checkpar,
        ) = f90trades.convert_trades_par_to_kepelem(
            self.system_parameters, self.fitting_parameters, self.n_bodies
        )
        # and in array, without the start!
        # M_Msun R_Rsun P_d a_AU e w_deg mA_deg inc_deg lN_deg
        self.kep_elem = np.zeros((self.n_planets, 9))
        self.kep_elem[:, 0] = self.mass[1:]
        self.kep_elem[:, 1] = self.radius[1:]
        self.kep_elem[:, 2] = self.period[1:]
        self.kep_elem[:, 3] = self.sma[1:]
        self.kep_elem[:, 4] = self.ecc[1:]
        self.kep_elem[:, 5] = self.argp[1:]
        self.kep_elem[:, 6] = self.meanA[1:]
        self.kep_elem[:, 7] = self.inc[1:]
        self.kep_elem[:, 8] = self.longN[1:]

        # computes the derived/physical parameters and their names
        self.physical_names, self.physical_parameters, self.physical_posterior = anc.compute_physical_parameters(
            self.n_bodies, 
            self.fitting_parameters, 
            self.fitting_names,
            self.system_parameters,
            mass_conv_factor=self.m_factor,
            radius_conv_factor=self.r_factor
        )
        # units of the derived/physical parameters
        self.physical_units = anc.get_units(self.physical_names, self.mass_unit)
        self.nphy = len(self.physical_parameters)

        # RADIAL VELOCITIES SET
        self.n_rv = f90trades.nrv
        self.n_set_rv = f90trades.nrvset  # number of jitter parameters

        # TRANSITS SET
        self.n_t0 = f90trades.nt0
        self.n_t0_sum = f90trades.ntts
        self.n_set_t0 = 0
        for i in range(0, self.n_planets):
            if self.n_t0[i] > 0:
                self.n_set_t0 += 1

        self.pephem = f90trades.pephem
        self.tephem = f90trades.tephem
        self.wrttime = f90trades.wrttime

        return

    def init_trades(self):

        self.init_trades_from_path(self.work_path, self.sub_folder)

        return

    # def set_args(
    #     self,
    # ):

    #     return

    def path_change(self, new_path):

        f90trades.path_change(os.path.join(new_path, ""))

        return

    def write_time_change(self, new_wrttime):

        f90trades.wrttime = new_wrttime
        self.wrttime = f90trades.wrttime

        return

    def update_system_fitting_parameters_from_keplerian_input(
        self,
        mass,
        radius,
        period,
        ecc,
        argp,
        meanA,
        inc,
        longN,
    ):

        n_all_par = len(self.system_parameters)
        all_par = np.zeros((n_all_par))
        fit_par = np.zeros((self.nfit))
        (
            all_par,
            fit_par,
        ) = f90trades.update_parameters_from_keplerian(
            mass,
            radius,
            period,
            ecc,
            argp,
            meanA,
            inc,
            longN,
        )
        self.system_parameters = all_par.copy()
        self.fitting_parameters = fit_par.copy()

        return

    def update_system_fitting_parameters(self):

        self.update_system_fitting_parameters_from_keplerian_input(
            self.mass,
            self.radius,
            self.period,
            self.ecc,
            self.argp,
            self.meanA,
            self.inc,
            self.longN,
        )

        return

    def run_and_get_stats_from_parameters(self, fit_pars):

        (
            chi_square,
            reduced_chi_square,
            lgllhd,
            lnprior,
            ln_const,
            bic,
            check,
        ) = f90trades.fortran_fitness_function(np.asarray(fit_pars, dtype=np.float64))

        return chi_square, reduced_chi_square, lgllhd, lnprior, ln_const, bic, check

    def run_and_get_stats(self):

        (
            chi_square,
            reduced_chi_square,
            lgllhd,
            lnprior,
            ln_const,
            bic,
            check,
        ) = self.run_and_get_stats_from_parameters(self.fitting_parameters)

        return chi_square, reduced_chi_square, lgllhd, lnprior, ln_const, bic, check

    def get_all_observables_from_parameters(self, fit_par):

        # run TRADES to get all TTs and RVs during the integration time
        nt0_full, nrv_nmax = f90trades.get_ntt_nrv(fit_par)
        # anc.print_both("full nRV = {:d} and nT0s = {}".format(nrv_nmax, nt0_full))
        (
            ttra_full,
            dur_full,
            id_ttra_full,
            stats_ttra,
            time_rv_nmax,
            rv_nmax,
            stats_rv,
        ) = f90trades.fit_par_to_ttra_rv(fit_par, nt0_full, nrv_nmax)
        
        return ttra_full, dur_full, id_ttra_full, stats_ttra, time_rv_nmax, rv_nmax, stats_rv

    def save_models_from_parameters(self, fit_par, smp_h5, smp_name):

        (
            ttra_full,
            dur_full,
            id_ttra_full,
            stats_ttra,
            time_rv_nmax,
            rv_nmax,
            stats_rv,
        ) = self.get_all_observables_from_parameters(fit_par)


        gr = smp_h5.create_group(smp_name)
        gr.create_dataset("fitting_parameters", data=fit_par, dtype=np.float64, compression="gzip")
        gr["fitting_parameters"].attrs["fitting_names"] = self.fitting_names
        gr.attrs["tepoch"] = self.tepoch

        stats_rv = np.array(stats_rv).astype(bool)
        time_rvx, rvx = time_rv_nmax[stats_rv], rv_nmax[stats_rv]
        id_rv = np.argsort(time_rvx)
        time_rv, rv = time_rvx[id_rv], rvx[id_rv]

        # anc.print_both("time_rv min = {:.5f} max = {:.5f}".format(np.min(time_rv), np.max(time_rv)))
        gr.create_dataset("time_rv_mod", data=time_rv,
                        dtype=np.float64, compression="gzip")
        gr.create_dataset("rv_mod", data=rv, dtype=np.float64, compression="gzip")

        stats_ttra = np.array(stats_ttra).astype(bool)
        for ipl in range(2, self.n_bodies + 1):
            sel_tra = np.logical_and(
                np.logical_and(id_ttra_full == ipl,
                            stats_ttra), ttra_full > -9.0e10
            )
            ttra = np.array(ttra_full)[sel_tra]
            dur = np.array(dur_full)[sel_tra]
            idx_tra = np.argsort(ttra)
            ttra = ttra[idx_tra]
            dur = dur[idx_tra]
            gr.create_dataset(
                "TTs_{:d}".format(ipl), data=ttra, dtype=np.float64, compression="gzip"
            )
            gr.create_dataset(
                "T41s_{:d}".format(ipl), data=dur, dtype=np.float64, compression="gzip"
            )
            # if np.sum(sel_tra) > 0:
            #     anc.print_both("T0_{:02d} min = {:.5f} max = {:.5f}".format(ipl, np.min(ttra), np.max(ttra)))

        return

    def run_and_write_summary_files(self, sim_name, fit_pars, phy_pars, id_sim=1, sigma_fit = None, sigma_phy=None, par_description="", return_stats=False):

        print("RUNNING sim {}".format(sim_name))
        # print("with parameters:")
        # print("{}".format(fit_pars))

        out_folder = os.path.join(self.full_path, sim_name, "")
        os.makedirs(out_folder, exist_ok=True)
        f90trades.path_change(out_folder)

        (
            chi_square,
            reduced_chi_square,
            lgllhd,
            lnprior,
            ln_const,
            bic,
            check,
        ) = f90trades.write_summary_files(id_sim, fit_pars)

        # smp_name = "{:04d}_{:s}".format(id_sim, sim_name)
        smp_name = sim_name
        smp_file = os.path.join(out_folder, "{}_models.hdf5".format(smp_name))
        with h5py.File(smp_file, "w") as smp_h5:

            self.save_models_from_parameters(
                fit_pars, smp_h5, smp_name
            )

            gr = smp_h5[smp_name]
            gr.attrs["chi_square"] = chi_square
            gr.attrs["reduced_chi_square"] = reduced_chi_square
            gr.attrs["lgllhd"] = lgllhd
            gr.attrs["lnprior"] = lnprior
            gr.attrs["ln_const"] = ln_const
            gr.attrs["bic"] = bic

        mass = np.zeros_like(self.mass)
        radius, period, sma, ecc, argp, meanA, inc, longN = (
            np.zeros_like(mass),
            np.zeros_like(mass),
            np.zeros_like(mass),
            np.zeros_like(mass),
            np.zeros_like(mass),
            np.zeros_like(mass),
            np.zeros_like(mass),
            np.zeros_like(mass),
        )

        (
            mass,
            radius,
            period,
            sma,
            ecc,
            argp,
            meanA,
            inc,
            longN,
            checkpar,
        ) = f90trades.convert_trades_par_to_kepelem(
            self.system_parameters, fit_pars, self.n_bodies
        )
        # M_Msun R_Rsun P_d a_AU e w_deg mA_deg inc_deg lN_deg
        kep_elem = np.zeros((self.n_planets, 9))
        kep_elem[:, 0] = mass[1:]
        kep_elem[:, 1] = radius[1:]
        kep_elem[:, 2] = period[1:]
        kep_elem[:, 3] = sma[1:]
        kep_elem[:, 4] = ecc[1:]
        kep_elem[:, 5] = argp[1:]
        kep_elem[:, 6] = meanA[1:]
        kep_elem[:, 7] = inc[1:]
        kep_elem[:, 8] = longN[1:]


        top_header, header = anc.get_header(anc.percentile_val)
        out_file = os.path.join(out_folder, "parameters_summary.txt")
        with open(out_file, "w") as out:

            if sigma_fit is None:
                sigma_print = np.zeros((8, self.nfit))
            else:
                sigma_print = sigma_fit
            anc.print_both("", out)
            anc.print_both("# " + "=" * 40, out)
            anc.print_both("\n# SUMMARY: {:s} ==> {:s}".format(sim_name, par_description), out)
            anc.print_both("# FITTED PARAMETERS", out)
            anc.print_parameters(
                top_header,
                header,
                self.fitting_names,
                self.fitting_units,
                fit_pars,
                sigma_parameters=sigma_print,
                output=out,
            )
            
            if sigma_phy is None:
                sigma_print = np.zeros((8, self.nphy))
            else:
                sigma_print = sigma_phy
            anc.print_both("# DERIVED PARAMETERS", out)
            anc.print_parameters(
                top_header,
                header,
                self.physical_names,
                self.physical_units,
                phy_pars,
                sigma_parameters=sigma_print,
                output=out,
            )
            if not bool(check):
                anc.print_both(
                    "WRITING WARNING FILE: {:s}".format(
                        os.path.join(out_folder, "WARNING.txt")
                    ),
                    out,
                )
                warn_o = open(os.path.join(out_folder, "WARNING.txt"), "w")
                warn_o.write(
                    "*******\nWARNING: FITTED PARAMETERS COULD NOT BE PHYSICAL!\nWARNING: BE VERY CAREFUL WITH THIS PARAMETER SET!\n*******"
                )
                warn_o.close()

            anc.print_both("", out)
            anc.print_both("# " + "=" * 40, out)
            anc.print_both("", out)

        f90trades.path_change(self.full_path)

        if return_stats:
            return chi_square, reduced_chi_square, lgllhd, lnprior, ln_const, bic, check

        return


    

    def reset(self):

        f90trades.deallocate_variables()

        return
