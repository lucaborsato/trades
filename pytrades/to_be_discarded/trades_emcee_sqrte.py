#!/usr/bin/env python
# -*- coding: utf-8 -*-

# no more "zero" integer division bugs!:P
import argparse
import os
import numpy as np  # array
import h5py
import sys
import time

# import glob

from scipy import optimize as sopt

# import multiprocessing as mp
from multiprocessing import Pool

# from schwimmbad import JoblibPool as Pool

import emcee

from constants import Mjups, Msear
import ancillary as anc
from pytrades_lib import f90trades

# =============================================================================


def get_args():
    parser = argparse.ArgumentParser(description="TRADES+EMCEE")

    # PATH FOLDER: full_path
    parser.add_argument(
        "-p",
        "--path",
        action="store",
        dest="full_path",
        required=True,
        help="The path (absolute or relative) with simulation files for TRADES.",
    )

    # SUBFOLDER TO SAVE ALL SIMULATION DATA
    parser.add_argument(
        "-s",
        "--sub-folder",
        "--sb",
        action="store",
        dest="sub_folder",
        default="emcee_run",
        help="Sub-folder name, without full path. Default = emcee_run",
    )

    # NUMBER OF CPU TO USE WITH EMCE !!
    parser.add_argument(
        "-c",
        "--cpu",
        "--nthreads",
        action="store",
        dest="nthreads",
        default=1,
        help="Number of threads to use. default nthreads = 1.",
    )

    # NUMBER OF WALKERS TO USE WITH EMCEE
    parser.add_argument(
        "-nw",
        "--nwalkers",
        "-np",
        "--npop",
        action="store",
        dest="nwalkers",
        default=1,
        help="Number of walkers (or number of chains) to use. default nwalkers = nfit*2",
    )

    # NUMBER OF STEPS/RUNS TO DO FOR EACH WALKER OF EMCEE
    parser.add_argument(
        "-nr",
        "--nruns",
        "-ns",
        "--nsteps",
        action="store",
        dest="nruns",
        default=10000,
        help="Number of runs/steps to use for each chain. default nruns = 10000.",
    )

    # NUMBER OF STEPS/RUNS TO SAVE TEMPORARY EMCEE SIMULATION
    parser.add_argument(
        "--isave",
        "--iter-save",
        "--iterations-save",
        action="store",
        dest="nsave",
        default="False",
        help="Number of iterations to do for save temporary chain. No intermediate save.",
    )

    # COMPUTE OR NOT CONSTANT TO ADD TO THE LGLIKELIHOOD: - (1/2dof) * SUM( ln 2pi * sigma_obs^2 )
    parser.add_argument(
        "-l",
        "--ln-err",
        "--ln-err-const",
        action="store",
        dest="ln_flag",
        default=True,
        help="Computes or not constant to add to the lglikelihood: - (1/2dof) * SUM( ln 2pi * sigma_obs^2 ). Default = True. Set to False if not use this value.",
    )

    parser.add_argument(
        "-ds",
        "--d-sigma",
        "--delta-sigma",
        action="store",
        dest="delta_sigma",
        default="1.e-4",
        type=str,
        help="Value of the sigma to compute initial walkers from initial solution. Default=1.e-4",
    )

    parser.add_argument(
        "-e",
        "--emcee-previous",
        action="store",
        dest="emcee_previous",
        default="None",
        type=str,
        help='Provide an existing "emcee_summary.hdf5" file if you wanto to start from the last step of that simulation. Default is None, create new initial walkers.',
    )

    parser.add_argument(
        "-progress",
        "--progress",
        "-emcee-progress",
        "--emcee-progress",
        action="store",
        dest="emcee_progress",
        default=True,
        help="Show the progress of emcee or not. Default True",
    )

    parser.add_argument(
        "--trades-fin",
        "--trades-final-previous",
        "--trades-final",
        action="store",
        dest="trades_previous",
        default="None",
        help="Define file from a previous TRADES simulation. File name and structure should be of type X_Y_finalNpar.dat. The parameters from this file will be the new original parameters",
    )

    parser.add_argument(
        "-seed",
        "--seed",
        action="store",
        dest="seed",
        default="None",
        help="Seed for random number generator. Default is None.",
    )

    cli = parser.parse_args()

    cli.full_path = os.path.join(os.path.abspath(cli.full_path), "")
    cli.sub_folder = os.path.join(os.path.relpath(cli.sub_folder), "")

    cli.nthreads = int(cli.nthreads)
    cli.nwalkers = int(cli.nwalkers)
    cli.nruns = int(cli.nruns)

    cli.ln_flag = anc.set_bool_argument(cli.ln_flag)

    cli.emcee_previous = anc.set_adhoc_file(cli.emcee_previous)
    cli.trades_previous = anc.set_adhoc_file(cli.trades_previous)

    cli.emcee_progress = anc.set_bool_argument(cli.emcee_progress)

    cli.seed = anc.set_int_or_none(cli.seed)
    # try:
    #     cli.seed = int(cli.seed)
    #     if cli.seed <= 0:
    #         cli.seed = None
    # except:
    #     cli.seed = None

    return cli


# =============================================================================

# INITIALISE FOLDER AND LOG FILE
init_folder = anc.init_folder

# ==============================================================================


def get_emcee_arguments(cli, nfit):

    # NUMBER OF WALKERS
    if cli.nwalkers < nfit * 2:
        nwalkers = nfit * 2
    else:
        nwalkers = cli.nwalkers
    if (nwalkers % 2) != 0:
        nwalkers += 1

    # NUMBER OF STEPS/RUNS FOR EACH WALKER
    if cli.nruns < 1:
        nruns = 10000
    else:
        nruns = cli.nruns

    try:
        nsave = int(cli.nsave)
        if nsave <= 0 or nsave >= nruns:
            nsave = False
    except:
        nsave = False

    npost = 0

    return nwalkers, nruns, nsave, npost


# =============================================================================

compute_proper_sigma = anc.compute_proper_sigma
compute_initial_walkers = anc.compute_initial_walkers

# =============================================================================

# =============================================================================
# =============================================================================
# MAIN SCRIPT - NOT IN FUNCTION DUE TO ISSUE WITH PARALLEL AND PICKLE FUNCTION OF LNPROB..
# =============================================================================
# =============================================================================

# def main():

# MAIN -- TRADES + EMCEE
# READ COMMAND LINE ARGUMENTS
cli = get_args()

# STARTING TIME
start = time.time()

# RENAME
working_path = cli.full_path
nthreads = cli.nthreads
np.random.seed(cli.seed)

# INITIALISE TRADES WITH SUBROUTINE WITHIN TRADES_LIB -> PARAMETER NAMES, MINMAX, INTEGRATION ARGS, READ DATA ...
pytrades.initialize_trades(working_path, cli.sub_folder, nthreads)

# RETRIEVE DATA AND VARIABLES FROM TRADES_LIB MODULE

n_bodies = pytrades.n_bodies  # NUMBER OF TOTAL BODIES OF THE SYSTEM
n_planets = n_bodies - 1  # NUMBER OF PLANETS IN THE SYSTEM
ndata = pytrades.ndata  # TOTAL NUMBER OF DATA AVAILABLE
nfit = pytrades.nfit  # NUMBER OF PARAMETERS TO FIT
nfree = pytrades.nfree  # NUMBER OF FREE PARAMETERS (ie nrvset)
dof = pytrades.dof  # NUMBER OF DEGREES OF FREEDOM = NDATA - NFIT
global inv_dof
inv_dof = pytrades.inv_dof

# READ THE NAMES OF THE PARAMETERS FROM THE TRADES_LIB AND CONVERT IT TO PYTHON STRINGS
str_len = pytrades.str_len
temp_names = pytrades.get_parameter_names(nfit, str_len)
trades_names = anc.convert_fortran_charray2python_strararray(temp_names)
# OLD (UNTIL TRADES v2.17.1)
# parameter_names = anc.trades_names_to_emcee(trades_names)
# NEW (TRADES v2.18.0)
parameter_names = trades_names

if cli.trades_previous is not None:
    temp_names, trades_parameters = anc.read_fitted_file(cli.trades_previous)
    if nfit != np.shape(trades_parameters)[0]:
        anc.print_both(
            " NUMBER OF PARAMETERS (%d) IN TRADES-PREVIOUS FILE DOES NOT"
            "MATCH THE CURRENT CONFIGURATION nfit=%d\nSTOP"
            % (np.shape(trades_parameters)[0], nfit)
        )
        sys.exit()
    del temp_names
else:
    # INITIAL PARAMETER SET (NEEDED ONLY TO HAVE THE PROPER ARRAY/VECTOR)
    trades_parameters = pytrades.fitting_parameters

# save initial_fitting parameters into array
original_fit_parameters = trades_parameters.copy()
# OLD v2.17.1
# fitting_parameters = anc.e_to_sqrte_fitting(trades_parameters, trades_names) # it should not change sqrt(e) as its from v2.18.0
# trades_minmax = pytrades.parameters_minmax  # PARAMETER BOUNDARIES
# parameters_minmax = anc.e_to_sqrte_boundaries(trades_minmax, trades_names)
# NEW v2.18.0
fitting_parameters = trades_parameters.copy()
parameters_minmax = pytrades.parameters_minmax.copy()

# RADIAL VELOCITIES SET
n_rv = pytrades.nrv
n_set_rv = pytrades.nrvset  # number of jitter parameters

# TRANSITS SET
n_t0 = pytrades.nt0
n_t0_sum = pytrades.ntts
n_set_t0 = 0
for i in range(0, n_bodies - 1):
    if n_t0[i] > 0:
        n_set_t0 += 1

# compute global constant for the loglhd
# global ln_err_const
# ln_err_const = pytrades.ln_err_const # with no RV jitters!

# SET EMCEE PARAMETERS:
nwalkers, nruns, nsave, _ = get_emcee_arguments(cli, nfit)

# uses the ancillary.get_fitted(full_path)
# to obtain info for the conversion from fitted to physical parameters
# nfit, NB, bodies_file, id_fit, id_all, nfit_list, cols_list, case =  anc.get_fitted(working_path)
_, _, _, id_fit, _, _, cols_list, case = anc.get_fitted(working_path)
# stellar mass in solar unit to earth, priors m in Earth masses
m_factor = pytrades.mr_star[0, 0] * Msear

# 2022-05-09 reads the priors from fortran and compute lnpriors in fortran
#
# LOGPROBABILITY FUNCTION NEEDED BY EMCEE
#
def lnprob(fitting_parameters):

    loglhd = 0.0
    check = 1
    loglhd, logprior, check = pytrades.fortran_logprob(
        np.asarray(fitting_parameters, dtype=float)
    )
    # not needed to add ln_err_const, because it is included in the loglhd now
    # loglhd = loglhd + ln_err_const # ln_err_const: global variable
    if check == 0:
        return -np.inf

    return loglhd + logprior

def minimize_func(fitting_parameters):

    lnL = lnprob(fitting_parameters)
    to_min = -2.0 * lnL
    if np.isinf(to_min):
        to_min = cst.huge

    return to_min

def residuals_func(fitting_parameters):

    to_min = minimize_func(fitting_parameters)
    resw = np.zeros(ndata) + to_min/ndata

    return resw


# INITIALISE SCRIPT FOLDER/LOG FILE
working_folder, _, of_run = init_folder(working_path, cli.sub_folder)

anc.print_both("", of_run)
anc.print_both(" ======== ", of_run)
anc.print_both(" pyTRADES", of_run)
anc.print_both(" ======== ", of_run)
anc.print_both("", of_run)
anc.print_both(" WORKING PATH = {:s}".format(working_path), of_run)
anc.print_both(" NUMBER OF THREADS = {:d}".format(nthreads), of_run)
anc.print_both(
    " dof = ndata({:d}) - nfit({:d}) - nfree({:d}) = {:d}".format(
        ndata, nfit, nfree, dof
    ),
    of_run,
)
anc.print_both(" Total N_RV = {:d} for {:d} set(s)".format(n_rv, n_set_rv), of_run)
anc.print_both(
    " Total N_T0 = {:d} for {:d} out of {:d} planet(s)".format(
        n_t0_sum, n_set_t0, n_planets
    ),
    of_run,
)
anc.print_both(" seed = {:s}".format(str(cli.seed)), of_run)

# INITIALISES THE WALKERS
if cli.emcee_previous is not None:
    anc.print_both(
        " Use a previous emcee simulation: %s" % (cli.emcee_previous), of_run
    )
    last_p0, old_nwalkers, last_done = anc.get_last_emcee_iteration(
        cli.emcee_previous, nwalkers
    )
    if not last_done:
        anc.print_both(
            "**STOP: USING A DIFFERENT NUMBER OF WALKERS (%d) W.R.T. PREVIOUS EMCEE SIMULATION (%d)."
            % (nwalkers, old_nwalkers),
            of_run,
        )
        sys.exit()
    p0 = last_p0

else:

    if cli.trades_previous is not None:
        anc.print_both(
            "\n ******\n INITIAL FITTING PARAMETERS FROM PREVIOUS"
            " TRADES-EMCEE SIM IN FILE:\n %s\n ******\n" % (cli.trades_previous),
            of_run,
        )

    anc.print_both(" ORIGINAL PARAMETER VALUES -> 0000", of_run)
    (
        chi_square_0000,
        reduced_chi_square_0000,
        lgllhd_0000,
        lnprior_0000,
        ln_const_0000,
        bic_0000,
        check_0000,
    ) = pytrades.write_summary_files(0, original_fit_parameters)
    anc.print_both(" ", of_run)
    anc.print_both(" TESTING LNPROB_SQ ...", of_run)

    # lgllhd_zero = lnprob(trades_parameters)
    # lgllhd_sq_zero = lnprob_sq(fitting_parameters, parameter_names)
    # lgllhd_sq_zero = lnprob_sq(fitting_parameters, parameter_names)
    lgllhd_sq_zero = lnprob(fitting_parameters)
    lgllhd_zero = lgllhd_sq_zero

    anc.print_both(" ", of_run)
    anc.print_both(
        " %15s %23s %23s %15s %23s"
        % ("trades_names", "original_trades", "trades_par", "emcee_names", "emcee_par"),
        of_run,
    )
    for ifit in range(0, nfit):
        anc.print_both(
            " %15s %23.16e %23.16e %15s %23.16e"
            % (
                trades_names[ifit],
                original_fit_parameters[ifit],
                trades_parameters[ifit],
                parameter_names[ifit],
                fitting_parameters[ifit],
            ),
            of_run,
        )
    anc.print_both(" ", of_run)
    anc.print_both(
        " %15s %23.16e %23.16e %15s %23.16e"
        % ("lnprob", lgllhd_0000, lgllhd_zero, "lnprob_sq", lgllhd_sq_zero),
        of_run,
    )
    anc.print_both(" ", of_run)

    if pytrades.lmon == 1:

        # TESTING A LOCAL FIT BEFORE EMCEE
        anc.print_both("TESTING A LOCAL FIT BEFORE EMCEE", of_run)

        # minimize with Nelder-Mead
        # res_opt = sopt.minimize(
        #     minimize_func,
        #     fitting_parameters,
        #     method="Nelder-Mead",
        #     options={"disp": True, "adaptive": False},
        # )
        # fitting_parameters = res_opt.x.copy()
        # least_squares
        res_opt = sopt.least_squares(
            residuals_func,
            fitting_parameters,
            method='lm',
            verbose=1
        )
        fitting_parameters = res_opt.x.copy()

        # lnL_opt = lnprob_sq(fitting_parameters, parameter_names)
        lnL_opt = lnprob(fitting_parameters)
        for i, p in enumerate(parameter_names):
            anc.print_both("{} {}".format(p, fitting_parameters[i]), of_run)

        # OLD v2.17.1
        # fit_trades_opt = anc.sqrte_to_e_fitting(fitting_parameters, parameter_names)
        # NEW v2.18.0
        fit_trades_opt = fitting_parameters.copy()
        (
            chi_square_opt,
            reduced_chi_square_opt,
            lgllhd_opt,
            lnc_opt,
            bic_opt,
            check_opt,
        ) = pytrades.write_summary_files(1, fit_trades_opt)
        anc.print_both("", of_run)
        anc.print_both("lnL = {}".format(lnL_opt), of_run)
        anc.print_both("lLH = {}".format(lgllhd_opt), of_run)
        anc.print_both("lnC = {}".format(lnc_opt), of_run)
        anc.print_both("rc2 = {}".format(reduced_chi_square_opt), of_run)
        anc.print_both("bic = {}".format(bic_opt), of_run)
        anc.print_both("", of_run)

    p0 = compute_initial_walkers(
        # lnprob_sq,
        lnprob,
        nfit,
        nwalkers,
        fitting_parameters,
        parameters_minmax,
        cli.delta_sigma,
        parameter_names,
        args=(),
        of_run=of_run,
    )

anc.print_both(" emcee chain: nwalkers = %d nruns = %d" % (nwalkers, nruns), of_run)
anc.print_both(" sampler ... ", of_run)

if nthreads > 1:
    threads_pool = Pool(nthreads)
else:
    threads_pool = None

sampler = emcee.EnsembleSampler(
    nwalkers,
    nfit,
    lnprob,
    pool=threads_pool,
    moves=[
        (emcee.moves.DEMove(), 0.8),
        (emcee.moves.DESnookerMove(), 0.2),
    ],
)

anc.print_both(" ready to go", of_run)
anc.print_both(" with nsave = {}".format(nsave), of_run)
sys.stdout.flush()

if nsave != False:
    # save temporary sampling during emcee every nruns*10%
    if os.path.exists(
        os.path.join(working_folder, "emcee_summary.hdf5")
    ) and os.path.isfile(os.path.join(working_folder, "emcee_summary.hdf5")):
        os.remove(os.path.join(working_folder, "emcee_summary.hdf5"))

    f_hdf5 = h5py.File(os.path.join(working_folder, "emcee_summary.hdf5"), "a")
    f_hdf5.create_dataset(
        "parameter_names", data=anc.encode_list(parameter_names), dtype="S10"
    )
    f_hdf5.create_dataset("boundaries", data=parameters_minmax, dtype=float)
    f_hdf5.create_dataset("chains", (nwalkers, nruns, nfit), dtype=float)
    f_hdf5["chains"].attrs["nwalkers"] = nwalkers
    f_hdf5["chains"].attrs["nruns"] = nruns
    f_hdf5["chains"].attrs["nfit"] = nfit
    f_hdf5["chains"].attrs["nfree"] = nfree
    f_hdf5.create_dataset("lnprobability", (nwalkers, nruns), dtype=float)
    # f_hdf5['lnprobability'].attrs['ln_err_const'] = ln_err_const
    f_hdf5.create_dataset(
        "acceptance_fraction", data=np.zeros((nfit)), dtype=float
    )
    f_hdf5.create_dataset("autocor_time", data=np.zeros((nfit)), dtype=float)
    f_hdf5.close()

    pos = p0
    niter_save = int(nruns / nsave)
    anc.print_both(" Running emcee with temporary saving", of_run)
    sys.stdout.flush()
    for i in range(0, niter_save):
        anc.print_both("", of_run)
        anc.print_both(" iter: {0:6d} ".format(i + 1), of_run)
        aaa = i * nsave
        bbb = aaa + nsave
        pos = sampler.run_mcmc(
            pos,
            nsave,
            progress=cli.emcee_progress,
            tune=True,  # default False
            skip_initial_state_check=False,  # default False
        )
        anc.print_both("completed {0:d} steps of {1:d}".format(bbb, nruns), of_run)

        f_hdf5 = h5py.File(os.path.join(working_folder, "emcee_summary.hdf5"), "a")
        temp_dset = f_hdf5["chains"]
        temp_dset[:, aaa:bbb, :] = sampler.chain[:, aaa:bbb, :]
        temp_dset.attrs["completed_steps"] = bbb
        temp_lnprob = f_hdf5["lnprobability"]
        try:
            temp_lnprob[:, aaa:bbb] = sampler.lnprobability[:, aaa:bbb]
        except:
            temp_lnprob[:, aaa:bbb] = sampler.lnprobability.T[:, aaa:bbb]
        mean_acceptance_fraction = np.mean(sampler.acceptance_fraction)
        acor_time = anc.compute_acor_time(sampler, steps_done=bbb)
        temp_acor = f_hdf5["autocor_time"]
        temp_acor[...] = acor_time
        f_hdf5.close()
        sys.stdout.flush()

    anc.print_both("", of_run)
    anc.print_both(
        "...done with saving temporary total shape = {}".format(
            str(np.shape(sampler.chain))
        ),
        of_run,
    )
    anc.print_both("", of_run)
    sys.stdout.flush()

else:
    # GOOD COMPLETE SINGLE RUNNING OF EMCEE, WITHOUT REMOVING THE BURN-IN
    anc.print_both(" Running full emcee ...", of_run)
    sys.stdout.flush()
    sampler.run_mcmc(
        p0,
        nruns,
        progress=cli.emcee_progress,
        tune=True,  # default False
        skip_initial_state_check=False,  # default False
    )
    anc.print_both("done", of_run)
    anc.print_both("", of_run)
    sys.stdout.flush()

    mean_acceptance_fraction = np.mean(sampler.acceptance_fraction)
    acor_time = anc.compute_acor_time(sampler)
    lnprobability = sampler.lnprobability

    # save chains with original shape as hdf5 file
    f_hdf5 = h5py.File(os.path.join(working_folder, "emcee_summary.hdf5"), "w")
    f_hdf5.create_dataset("chains", data=sampler.chain, dtype=float)
    f_hdf5["chains"].attrs["nwalkers"] = nwalkers
    f_hdf5["chains"].attrs["nruns"] = nruns
    f_hdf5["chains"].attrs["nfit"] = nfit
    f_hdf5["chains"].attrs["nfree"] = nfree
    f_hdf5["chains"].attrs["completed_steps"] = nruns
    f_hdf5.create_dataset(
        "parameter_names", data=anc.encode_list(parameter_names), dtype="S10"
    )
    f_hdf5.create_dataset("boundaries", data=parameters_minmax, dtype=float)
    f_hdf5.create_dataset(
        "acceptance_fraction", data=sampler.acceptance_fraction, dtype=float
    )
    f_hdf5["acceptance_fraction"].attrs[
        "mean_acceptance_fraction"
    ] = mean_acceptance_fraction
    f_hdf5.create_dataset("autocor_time", data=acor_time, dtype=float)
    f_hdf5.create_dataset("lnprobability", data=lnprobability, dtype=float)
    # f_hdf5['lnprobability'].attrs['ln_err_const'] = ln_err_const
    f_hdf5.close()

anc.print_both(
    " Mean_acceptance_fraction should be between [0.25-0.5] = {0:.6f}".format(
        mean_acceptance_fraction
    ),
    of_run,
)
anc.print_both("", of_run)

if threads_pool is not None:
    # close the pool of threads
    threads_pool.close()
    # threads_pool.terminate()
    threads_pool.join()

anc.print_both("COMPLETED EMCEE", of_run)

elapsed = time.time() - start
elapsed_d, elapsed_h, elapsed_m, elapsed_s = anc.computation_time(elapsed)

anc.print_both("", of_run)
anc.print_both(
    " pyTRADES: EMCEE FINISHED in {0:2d} day {1:02d} hour {2:02d} min {3:.2f} sec - bye bye".format(
        int(elapsed_d), int(elapsed_h), int(elapsed_m), elapsed_s
    ),
    of_run,
)

anc.print_both("", of_run)
of_run.close()
pytrades.deallocate_variables()

# return

# ==============================================================================
# ==============================================================================

# if __name__ == "__main__":
#   main()
