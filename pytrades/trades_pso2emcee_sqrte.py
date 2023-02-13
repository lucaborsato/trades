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

# import multiprocessing as mp
from multiprocessing import Pool

# from schwimmbad import JoblibPool as Pool

import emcee

from constants import Msear
import ancillary as anc
from pytrades_lib import pytrades

# =============================================================================


def get_args():
    parser = argparse.ArgumentParser(description="TRADES+PSO+EMCEE")

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

    # NUMBER OF CPU TO USE WITH EMCE AND PSO!!
    parser.add_argument(
        "-c",
        "--cpu",
        "--nthreads",
        action="store",
        dest="nthreads",
        default=1,
        help="Number of threads to use. default nthreads = 1. For PSO with openMP you have to set export OMP_NUM_THREADS=CPUNUMBER before to run this script",
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
        help="Number of iterations to do for save temporary chain. default each 0.1 of nruns.",
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
        "-progress",
        "--progress",
        "-emcee-progress",
        "--emcee-progress",
        action="store",
        dest="emcee_progress",
        default=True,
        help="Show the progress of emcee or not. Default True",
    )

    # HOW TO RUN PSO OR SKIP
    parser.add_argument(
        "-r",
        "--pso",
        "--pso-type",
        action="store",
        dest="pso_type",
        default="run",
        help='Define PSO run type: "skip" = do not run PSO, only emcee with random walkers; "run" = run PSO normally; "exists" = do not run PSO, but read previous PSO (pso_run.hdf5) file and start emcee from its population',
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
    cli.emcee_progress = anc.set_bool_argument(cli.emcee_progress)

    if cli.pso_type not in ["skip", "run", "exists"]:
        cli.pso_type = "run"

    try:
        cli.seed = int(cli.seed)
        if cli.seed <= 0:
            cli.seed = None
    except:
        cli.seed = None

    return cli


# =============================================================================
# INITIALISE FOLDER AND LOG FILE
init_folder = anc.init_folder

# =============================================================================


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

pso_to_emcee = anc.pso_to_emcee

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

# INITIAL PARAMETER SET (NEEDED ONLY TO HAVE THE PROPER ARRAY/VECTOR)
trades_parameters = pytrades.fitting_parameters.copy()
# OLD v2.17.1
# fitting_parameters = anc.e_to_sqrte_fitting(trades_parameters, trades_names)
# trades_minmax = pytrades.parameters_minmax  # PARAMETER BOUNDARIES
# parameters_minmax = anc.e_to_sqrte_boundaries(trades_minmax, trades_names)
# NEW v2.18.0
fitting_parameters = trades_parameters.copy()
trades_minmax = pytrades.parameters_minmax.copy()
parameters_minmax = trades_minmax

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
    # loglhd, check = pytrades.fortran_loglikelihood(
    #     np.asarray(fitting_parameters, dtype=np.float64)
    # )
    loglhd, logprior, check = pytrades.fortran_logprob(
        np.asarray(fitting_parameters, dtype=np.float64)
    )
    # not needed to add ln_err_const, because it is included in the loglhd now
    # loglhd = loglhd + ln_err_const # ln_err_const: global variable
    if check == 0:
        return -np.inf

    return loglhd + logprior


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

# INITIALISE PSO ARGUMENTS FROM pso.opt FILE
pytrades.init_pso(1, working_path)  # read PSO options
# PSO VARIABLES
np_pso = pytrades.np_pso
nit_pso = pytrades.nit_pso
n_global = pytrades.n_global
# n_global = 1
anc.print_both(
    " PSO n_global = {} npop = {} ngen = {}".format(n_global, np_pso, nit_pso), of_run
)
anc.print_both(" PSO run type: {:s}".format(cli.pso_type))

# RUN PSO+EMCEE n_global TIMES
for iter_global in range(0, n_global):

    # commented 2019-05-21
    # threads_pool = mp.Pool(1)

    # CREATES PROPER WORKING PATH AND NAME
    i_global = iter_global + 1
    pso_path = os.path.join(
        os.path.join(working_folder, "{0:04d}_pso2emcee".format(i_global)), ""
    )
    pytrades.path_change(pso_path)

    anc.print_both(
        "\n\n GLOBAL RUN {0:04d} INTO PATH: {1:s}\n".format(i_global, pso_path), of_run
    )

    if cli.pso_type == "run":
        # RUN PSO
        anc.print_both(" RUN PSO ...", of_run)

        pso_start = time.time()
        if not os.path.exists(pso_path):
            os.makedirs(pso_path)
        # copy files
        anc.copy_simulation_files(working_path, pso_path)

        # CALL RUN_PSO SUBROUTINE FROM TRADES_LIB: RUNS PSO AND COMPUTES THE BEST SOLUTION, SAVING ALL THE POPULATION EVOLUTION
        pso_parameters = trades_parameters.copy()
        pso_fitness = 0.0
        pso_parameters, pso_fitness = pytrades.pyrun_pso(nfit, i_global)
        anc.print_both(" completed run_pso", of_run)

        pso_best_evolution = np.asarray(
            pytrades.pso_best_evolution[...], dtype=np.float64
        )
        anc.print_both(" pso_best_evolution retrieved", of_run)

        anc.print_both(" last pso_best_evolution", of_run)
        last_pso_fitness = pso_best_evolution[-1, -1].astype(np.float64)
        anc.print_both(" fitness = {}".format(last_pso_fitness), of_run)

        # SAVE PSO SIMULATION IN pso_run.hdf5 FILE
        print(
            " Creating pso hdf5 file: {}".format(os.path.join(pso_path, "pso_run.hdf5"))
        )
        pso_hdf5 = h5py.File(os.path.join(pso_path, "pso_run.hdf5"), "w")
        pso_hdf5.create_dataset(
            "population", data=pytrades.population, dtype=np.float64
        )
        pso_hdf5.create_dataset(
            "population_fitness", data=pytrades.population_fitness, dtype=np.float64
        )
        pso_hdf5.create_dataset("pso_parameters", data=pso_parameters, dtype=np.float64)
        pso_hdf5.create_dataset(
            "pso_fitness", data=np.array(pso_fitness), dtype=np.float64
        )
        pso_hdf5.create_dataset(
            "pso_best_evolution", data=pso_best_evolution, dtype=np.float64
        )
        pso_hdf5.create_dataset(
            "parameters_minmax", data=trades_minmax, dtype=np.float64
        )
        pso_hdf5.create_dataset(
            "parameter_names", data=anc.encode_list(trades_names), dtype="S10"
        )
        pso_hdf5["population"].attrs["npop"] = np_pso
        pso_hdf5["population"].attrs["niter"] = nit_pso
        pso_hdf5["population"].attrs["iter_global"] = iter_global + 1
        pso_hdf5["population"].attrs["nfit"] = nfit
        pso_hdf5.close()

        anc.print_both(" ", of_run)
        # fitness_iter, lgllhd_iter, check_iter = pytrades.write_summary_files(i_global, pso_parameters)
        elapsed = time.time() - pso_start
        elapsed_d, elapsed_h, elapsed_m, elapsed_s = anc.computation_time(elapsed)
        anc.print_both(" ", of_run)
        anc.print_both(
            " PSO FINISHED in {0:2d} day {1:02d} hour {2:02d} min {3:.2f} sec - bye bye".format(
                int(elapsed_d), int(elapsed_h), int(elapsed_m), elapsed_s
            ),
            of_run,
        )

        # TRADES/PSO USES ECOSW/ESINW --> HERE EMCEE USES SQRTECOSW/SQRTESINW
        # emcee_parameters = anc.e_to_sqrte_fitting(pso_parameters, trades_names)
        emcee_parameters = pso_parameters.copy()
        if float(cli.delta_sigma) <= 0.0:
            p0 = pso_to_emcee(
                nfit,
                nwalkers,
                pytrades.population,
                parameter_names,
                #   sqrt_par=True
                sqrt_par=False,
            )
        else:
            # p0, pso_fitness_p0 = pso_to_emcee(nfit, nwalkers, population, population_fitness, pso_parameters, pso_fitness, pso_best_evolution)
            p0 = compute_initial_walkers(
                # lnprob_sq,
                lnprob,
                nfit,
                nwalkers,
                emcee_parameters,
                parameters_minmax,
                parameter_names,
                cli.delta_sigma,
                of_run,
            )

    elif cli.pso_type == "exists":
        # READ PREVIOUS PSO_RUN.HDF5 FILE AND INITIALISE POPULATION FOR EMCEE
        anc.print_both(
            " READ PREVIOUS PSO_RUN.HDF5 FILE AND INITIALISE POPULATION FOR EMCEE",
            of_run,
        )

        (
            population,
            population_fitness,
            pso_parameters,
            pso_fitness,
            pso_best_evolution,
            parameters_minmax,
            parameter_names,
            pop_shape,
        ) = anc.get_pso_data(os.path.join(pso_path, "pso_run.hdf5"))
        # population, population_fitness, pso_parameters, pso_fitness, pso_best_evolution, parameters_minmax, parameter_names, pop_shape = get_pso_data(os.path.join(pso_path, 'pso_run.hdf5'))

        # fitness_iter, lgllhd_iter, check_iter = pytrades.write_summary_files(i_global, pso_parameters)
        _, _, _ = pytrades.write_summary_files(i_global, pso_parameters)

        pso_fitness = float(pso_fitness)
        anc.print_both(
            " read pso_run.hdf5 file with best pso_fitness = {:.7f}".format(
                pso_fitness
            ),
            of_run,
        )

        # TRADES/PSO USES ECOSW/ESINW --> HERE EMCEE USES SQRTECOSW/SQRTESINW
        # emcee_parameters = anc.e_to_sqrte_fitting(pso_parameters, trades_names)
        emcee_parameters = pso_parameters.copy()
        if float(cli.delta_sigma) <= 0.0:
            p0 = pso_to_emcee(
                nfit,
                nwalkers,
                population,
                trades_names,
                # sqrt_par=True
                sqrt_par=False,
            )
        else:
            # p0, pso_fitness_p0 = pso_to_emcee(nfit, nwalkers, population, population_fitness, pso_parameters, pso_fitness, pso_best_evolution)
            p0 = compute_initial_walkers(
                # lnprob_sq,
                lnprob,
                nfit,
                nwalkers,
                emcee_parameters,
                parameters_minmax,
                parameter_names,
                cli.delta_sigma,
                of_run,
            )

    elif cli.pso_type == "skip":
        # DO NOT RUN PSO, ONLY EMCEE
        anc.print_both(" DO NOT RUN PSO, ONLY EMCEE", of_run)

        # p0 = [parameters_minmax[:,0] + np.random.random(nfit)*delta_parameters for i in range(0, nwalkers)]
        p0 = compute_initial_walkers(
            # lnprob_sq,
            lnprob,
            nfit,
            nwalkers,
            fitting_parameters,
            parameters_minmax,
            parameter_names,
            cli.delta_sigma,
            of_run,
        )
    # end if cli.pso_type

    anc.print_both(
        " emcee chain: nwalkers = {} nruns = {}".format(nwalkers, nruns), of_run
    )
    anc.print_both(" sampler ... ", of_run)

    # close the pool of threads
    # commented 2019-05-21
    # threads_pool.close()
    # threads_pool.join()

    if nthreads > 1:
        threads_pool = Pool(nthreads)
    else:
        threads_pool = None

    sampler = emcee.EnsembleSampler(
        nwalkers,
        nfit,
        # lnprob_sq,
        lnprob,
        pool=threads_pool,
        # args=[parameter_names],
        moves=[
            (emcee.moves.DEMove(), 0.8),
            (emcee.moves.DESnookerMove(), 0.2),
        ],
    )  # needed to use sqrt(e) in emcee instead of e (in fortran)

    anc.print_both(" ready to go", of_run)
    anc.print_both(" with nsave = {}".format(nsave), of_run)
    sys.stdout.flush()

    if nsave != False:
        # save temporary sampling during emcee every nruns*10%
        if os.path.exists(
            os.path.join(pso_path, "emcee_summary.hdf5")
        ) and os.path.isfile(os.path.join(pso_path, "emcee_summary.hdf5")):
            os.remove(os.path.join(pso_path, "emcee_summary.hdf5"))

        f_hdf5 = h5py.File(os.path.join(pso_path, "emcee_summary.hdf5"), "a")
        f_hdf5.create_dataset(
            "parameter_names", data=anc.encode_list(parameter_names), dtype="S10"
        )
        f_hdf5.create_dataset("boundaries", data=parameters_minmax, dtype=np.float64)
        f_hdf5.create_dataset("chains", (nwalkers, nruns, nfit), dtype=np.float64)
        f_hdf5["chains"].attrs["nwalkers"] = nwalkers
        f_hdf5["chains"].attrs["nruns"] = nruns
        f_hdf5["chains"].attrs["nfit"] = nfit
        f_hdf5["chains"].attrs["nfree"] = nfree
        f_hdf5.create_dataset("lnprobability", (nwalkers, nruns), dtype=np.float64)
        # f_hdf5['lnprobability'].attrs['ln_err_const'] = ln_err_const
        f_hdf5.create_dataset(
            "acceptance_fraction", data=np.zeros((nfit)), dtype=np.float64
        )
        f_hdf5.create_dataset("autocor_time", data=np.zeros((nfit)), dtype=np.float64)
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

            f_hdf5 = h5py.File(os.path.join(pso_path, "emcee_summary.hdf5"), "a")
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
            "...done with saving temporary total shape = {:s}".format(
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
        f_hdf5 = h5py.File(os.path.join(pso_path, "emcee_summary.hdf5"), "w")
        f_hdf5.create_dataset("chains", data=sampler.chain, dtype=np.float64)
        f_hdf5["chains"].attrs["nwalkers"] = nwalkers
        f_hdf5["chains"].attrs["nruns"] = nruns
        f_hdf5["chains"].attrs["nfit"] = nfit
        f_hdf5["chains"].attrs["nfree"] = nfree
        f_hdf5["chains"].attrs["completed_steps"] = nruns
        f_hdf5.create_dataset(
            "parameter_names", data=anc.encode_list(parameter_names), dtype="S10"
        )
        f_hdf5.create_dataset("boundaries", data=parameters_minmax, dtype=np.float64)
        f_hdf5.create_dataset(
            "acceptance_fraction", data=sampler.acceptance_fraction, dtype=np.float64
        )
        f_hdf5["acceptance_fraction"].attrs[
            "mean_acceptance_fraction"
        ] = mean_acceptance_fraction
        f_hdf5.create_dataset("autocor_time", data=acor_time, dtype=np.float64)
        f_hdf5.create_dataset("lnprobability", data=lnprobability, dtype=np.float64)
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

# if __name__ == "__main__":
#   main()
