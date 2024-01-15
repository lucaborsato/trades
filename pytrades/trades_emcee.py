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

import emcee
# try:
#     from pytransit.utils.de import DiffEvol
# except:
#     print("PyDE not within PyTransit")
#     try:
#         from pyde.de import DiffEvol
#     except:
#         print("PyDE not installed, errors could occur.")
from de import DiffEvol

from constants import Mjups, Msear
# from pytrades_lib import f90trades
import pytrades
import ancillary as anc
import de_plot as dep


# import numba
# numba.set_num_threads(1)
# numba.config.THREADING_LAYER = "tbb"
# numba.config.THREADING_LAYER = "'workqueue'"
# # numba.config.DISABLE_JIT = 1

# =============================================================================


# =============================================================================

init_folder = anc.init_folder
compute_proper_sigma = anc.compute_proper_sigma
compute_initial_walkers = anc.compute_initial_walkers

# =============================================================================

# =============================================================================
# =============================================================================
# MAIN SCRIPT - NOT IN FUNCTION DUE TO ISSUE WITH PARALLEL AND PICKLE FUNCTION OF LNPROB..
# =============================================================================
# =============================================================================

# MAIN -- TRADES + EMCEE
# READ COMMAND LINE ARGUMENTS
# cli = get_args()

yml_file = anc.get_input_file()
cli = anc.ConfigurationRun(yml_file)

os.environ["OMP_NUM_THREADS"] = "1"

# STARTING TIME
start = time.time()

# RENAME
working_path = cli.full_path
nthreads = cli.nthreads
np.random.seed(cli.seed)
# os.environ["OMP_NUM_THREADS"] = str(nthreads)

# INITIALISE TRADES WITH SUBROUTINE WITHIN TRADES_LIB -> PARAMETER NAMES, MINMAX, INTEGRATION ARGS, READ DATA ...
# pytrades.initialize_trades(working_path, cli.sub_folder, nthreads)
anc.print_both("Init trades ... ")
sim = pytrades.TRADESfolder(
    working_path,
    sub_folder=cli.sub_folder, 
    nthreads=cli.nthreads,
    seed=cli.seed,
    m_type=cli.m_type
)
sim.init_trades()
sys.stdout.flush()


def lnprob(fitting_parameters):

    (
        # chi_square,
        # reduced_chi_square,
        # lgllhd,
        # lnprior,
        # ln_const,
        # bic,
        # check,
        _,
        _,
        lgllhd,
        lnprior,
        _,
        _,
        check,
    ) = sim.run_and_get_stats_from_parameters(fitting_parameters)

    if check == 0:
        return -np.inf
        # return -2.0e30

    return lgllhd+lnprior

def lnprob_de(fitting_parameters):
    lnP = lnprob(fitting_parameters)
    if np.isinf(lnP):
        lnP = -2.0e30
    return lnP

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
sys.stdout.flush()

anc.print_both("")
anc.print_both("fitting parameters:")
anc.print_both(" ".join([n for n in sim.fitting_names]))
anc.print_both("")

# check pyde args
anc.print_both("\n Check pyDE configuration ... ", output=of_run)
de_run = cli.de_type.lower() == "run"
de_resume = cli.de_type.lower() == "resume"
de_to_emcee = cli.de_type.lower() == "to_emcee"

# SET INITIAL PARAMETERS FOR EMCEE

# PYDE
if de_run or de_resume or de_to_emcee:

    de_path = working_folder
    de_file = os.path.join(de_path, "de_run.hdf5")

    npop_de = cli.npop_de
    ngen_de = cli.ngen_de
    nfit    = sim.nfit
    ndata   = sim.ndata

    de_save = cli.nsave_de
    de_bounds = sim.fitting_minmax
    de_f = cli.de_f
    de_c = cli.de_c
    de_maximize = cli.de_maximize
    fit_type = -1 if de_maximize else 1
    anc.print_both(
        " pyDE de_type = {} ==> de_run = {}, de_resume = {}, de_to_emcee = {}".format(
            cli.de_type.lower(), de_run, de_resume, de_to_emcee
        ), output=of_run
    )

    if de_to_emcee:

        anc.print_both(
            " pyDE RUN already exist --> loading best-fit parameters and send them to emcee ...", output=of_run
        )
        fitting_parameters = anc.de_load_parameters(de_file)

    else:
        anc.print_both(" pyDE npop = {} ngen = {} nsave = {}".format(npop_de, ngen_de, de_save), output=of_run)

        # prepare arrays to store pyDE evolution
        de_pop, de_fit = (
            np.zeros((ngen_de, npop_de, nfit)),
            np.zeros((ngen_de, npop_de))-1.0e-30, # set to very small non-zero values
        )
        de_pop_best, de_fit_best = np.zeros((ngen_de, nfit)), np.zeros((ngen_de))-1.0e-30

        with Pool(nthreads) as threads_pool:
            # create pyDE object
            de_evol = DiffEvol(
                lnprob_de,
                de_bounds,
                npop_de,
                f=de_f,
                c=de_c,
                seed=cli.seed,
                maximize=de_maximize,
                vectorize=False,
                pool=threads_pool,
            )

            last_iter_de = 0
            if de_run:
                if os.path.exists(de_file) and os.path.isfile(de_file):
                    anc.print_both("File {} exists: deleting it!".format(de_file), output=of_run)
                    os.remove(de_file)
            elif de_resume:
                if os.path.exists(de_file) and os.path.isfile(de_file):
                    de_old = dep.DEEvolution(de_path)
                    # update de_pop, de_fit, de_pop_best, de_fit_best
                    de_pop[: de_old.ngen_de, :, :] = de_old.de_pop.copy()
                    de_fit[: de_old.ngen_de, :] = de_old.de_fit.copy()
                    de_pop_best[: de_old.ngen_de, :] = de_old.de_pop_best.copy()
                    de_fit_best[: de_old.ngen_de] = de_old.de_fit_best.copy()
                    # last_iter_de = de_old.ngen_de
                    last_iter_de = de_old.iter_de
                    de_evol._population = de_old.de_pop[-1, :, :].copy()
                    de_evol._fitness = de_old.de_fit[-1, :].copy()
                    de_evol._minidx = (
                        np.argmax(de_evol._fitness)
                        if de_maximize
                        else np.argmin(de_evol._fitness)
                    )
                else:
                    anc.print_both(" File {} does not exist, cannot resume pyDE, starting from scratch.".format(de_file), output=of_run)

            start_iter_de = time.time()
            de_start = start_iter_de
            anc.print_both(" DE - START ", output=of_run)
            sys.stdout.flush()
            if last_iter_de + 1 < ngen_de:
                anc.print_both(" last_iter_de + 1 = {} < ngen_de = {}".format(last_iter_de+1, ngen_de), output=of_run)
                
                # print("iter_de mod_save  check_save check_iter best_fitness")
                
                sys.stdout.flush()
                
                for iter_de_temp, res_de in enumerate(de_evol(ngen_de-last_iter_de)):
                    # current iter_de_temp starts from 0 and ends to ngen_de-last_iter_de
                    iter_de = last_iter_de + iter_de_temp
                    # iter_de goes from last_iter_de to ngen_de
                    mod_save = (iter_de + 1) % de_save
                    check_save = mod_save == 0
                    check_iter = iter_de + 1 >= ngen_de
                    # anc.print_both("iter_de = {}".format(iter_de), output=of_run)
                    # sys.stdout.flush()

                    de_pop[iter_de, :, :]   = de_evol.population.copy()
                    de_fit[iter_de, :]      = fit_type * de_evol._fitness.copy()
                    de_pop_best[iter_de, :] = de_evol.minimum_location.copy()
                    de_fit_best[iter_de]    = fit_type * de_evol.minimum_value

                    # print(iter_de, mod_save, check_save, check_iter, fit_type * de_evol.minimum_value)
                    
                    if iter_de > 0:
                        if (check_save or check_iter):
                            anc.print_both(" ============= ", output=of_run)
                            anc.print_both(" pyDE - iter = {} / {} ==> saving".format(iter_de+ 1, ngen_de), output=of_run)
                            anc.print_both(
                                " last best fitness = {}".format(de_fit_best[iter_de]), output=of_run
                            )
                            sys.stdout.flush()

                            # sys.exit()

                            # SAVE DE SIMULATION IN de_run.hdf5 FILE
                            anc.de_save_evolution(
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
                            anc.print_both(" Updated DE hdf5 file: {}".format(de_file), output=of_run)
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
                                ), output=of_run
                                
                            )
                            sys.stdout.flush()
                            start_iter_de = time.time()
                # end for iter_de_temp, res_de

                # FORCE
                # SAVE DE SIMULATION IN de_run.hdf5 FILE
                anc.de_save_evolution(
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
                anc.print_both(" Updated DE hdf5 file: {}".format(de_file), output=of_run)
            else:
                iter_de = ngen -1
        anc.print_both("", output=of_run)
        anc.print_both(" completed DE", output=of_run)
        anc.print_both("", output=of_run)

        elapsed = time.time() - de_start
        elapsed_d, elapsed_h, elapsed_m, elapsed_s = anc.computation_time(elapsed)
        anc.print_both(" ", output=of_run)
        anc.print_both(
            " DE FINISHED in {0:2d} day {1:02d} hour {2:02d} min {3:.2f} sec - bye bye".format(
                int(elapsed_d), int(elapsed_h), int(elapsed_m), elapsed_s
            ), output=of_run
        )
        sys.stdout.flush()
        fitting_parameters = anc.de_get_best_parameters(
            de_fit_best, de_pop_best, iter_de, de_maximize=de_maximize
        )
        # END DE

elif cli.trades_previous is not None:
    temp_names, fitting_parameters = anc.read_fitted_file(cli.trades_previous)
    if sim.nfit != np.shape(fitting_parameters)[0]:
        anc.print_both(
            " NUMBER OF PARAMETERS (%d) IN TRADES-PREVIOUS FILE DOES NOT"
            " MATCH THE CURRENT CONFIGURATION nfit=%d\nSTOP"
            % (np.shape(fitting_parameters)[0], sim.nfit), output=of_run
        )
        sys.exit()
    del temp_names
    anc.print_both("INPUT AD-HOC INITIAL FITTING PARAMETERS", output=of_run)
    anc.print_both("{:20s} {:23s}".format("name", "parameter"), output=of_run)
    for n, p in zip(sim.fitting_names, fitting_parameters):
        anc.print_both("{:20s} {:23.16e}".format(n, p), output=of_run)
    anc.print_both("", output=of_run)

else:
    # INITIAL PARAMETER SET (NEEDED ONLY TO HAVE THE PROPER ARRAY/VECTOR)
    anc.print_both(
            """
                pyDE type is not run or resume or to_emcee.
                trades_previous file not provided.
                Initial parameters directly to emcee.
            """, output=of_run
        )
    fitting_parameters = sim.fitting_parameters.copy()

sys.stdout.flush()
# numba.set_num_threads(1)

# save initial_fitting parameters into array
initial_parameters = fitting_parameters.copy()
(
    chi_square,
    reduced_chi_square,
    lgllhd,
    lnprior,
    ln_const,
    bic,
    check,
) = sim.run_and_get_stats_from_parameters(initial_parameters)
anc.print_both("Stats for the initial set of parameters:", output=of_run)
anc.print_both("chi_square = {}".format(chi_square), output=of_run)
anc.print_both("reduced_chi_square = {}".format(reduced_chi_square), output=of_run)
anc.print_both("lgllhd = {}".format(lgllhd), output=of_run)
anc.print_both("lnprior = {}".format(lnprior), output=of_run)
anc.print_both("ln_const = {}".format(ln_const), output=of_run)
anc.print_both("bic = {}".format( bic), output=of_run)
anc.print_both("check = {}".format(check), output=of_run)
anc.print_both("", output=of_run)
sys.stdout.flush()
# SET EMCEE PARAMETERS:
# nwalkers, nruns, nsave, _ = get_emcee_arguments(cli, nfit)

if cli.nruns > 0:

    anc.print_both(" emcee chain: nwalkers = %d nruns = %d" % (cli.nwalkers, cli.nruns), output=of_run)
    if cli.nwalkers < 2 * sim.nfit:
        cli.nwalkers = 2 * sim.nfit
        anc.print_both("WARNING: nwalkers < 2 x nfit ==> setting nalkers = 2 x nfit  = {}".format(
            cli.nwalkers
        ), output=of_run)
    sys.stdout.flush()

    # anc.print_both(" sampler ... ")

    backend_filename = os.path.join(working_folder, "emcee_summary.hdf5")
    backend = emcee.backends.HDFBackend(backend_filename, compression='gzip')
    if os.path.exists(backend_filename):
        completed_steps = backend.iteration
    else:
        completed_steps = 0

    if cli.emcee_restart:
        anc.print_both("keyword emcee_restart is True ==> resetting file content!")
        backend.reset(cli.nwalkers, sim.nfit)

    anc.print_both("Initial size of backend: {}".format(completed_steps), output=of_run)
    if  completed_steps > 0:
        # p0 = None
        p0 = backend.get_last_sample()
        anc.print_both("continue emcee analysis ...", output=of_run)
    else:
        anc.print_both("new run of emcee analysis ...", output=of_run)
        if cli.pre_optimise:
            def neg_lnprob(p):
                lnP = lnprob_de(p)
                return -lnP
            opt_res = sopt.minimize(
                neg_lnprob, 
                fitting_parameters, 
                method="Nelder-Mead", 
                bounds=[(bd[0], bd[1]) for bd in sim.fitting_minmax]
            )
            anc.print_both(
                "Optimised fitting parameters with Minimize(Nelder-Mead)", output=of_run
            )
            fitting_parameters = opt_res.x
            anc.print_both(
                "Opt Message: {}".format(opt_res.message), output=of_run
            )
            anc.print_both(
                "Opt Success: {}".format(opt_res.success), output=of_run
            )
        p0 = compute_initial_walkers(
            # lnprob_sq,
            lnprob,
            sim.nfit,
            cli.nwalkers,
            fitting_parameters,
            sim.fitting_minmax,
            cli.delta_sigma,
            sim.fitting_names,
            args=(),
            # of_run=None,
            of_run=of_run,
        )
    sys.stdout.flush()

    anc.print_both("emcee moves: {}".format(cli.emcee_move), output=of_run)

    # os.environ["OMP_NUM_THREADS"] = "1"
    with Pool(nthreads) as threads_pool:
        sampler = emcee.EnsembleSampler(
            cli.nwalkers,
            sim.nfit,
            lnprob,
            pool=threads_pool,
            # moves=[
            #     (emcee.moves.DEMove(), 0.8),
            #     (emcee.moves.DESnookerMove(), 0.2),
            # ],
            moves=cli.emcee_move,
            backend=backend,
        )

        anc.print_both("ready to go", output=of_run)
        sys.stdout.flush()

        sampler.run_mcmc(
            p0, 
            cli.nruns,
            progress=cli.emcee_progress,
            thin_by=cli.thin_by,
            tune=True,  # default False
            skip_initial_state_check=False,  # default False
        )
    sys.stdout.flush()
    anc.print_both("", output=of_run)
    anc.print_both("COMPLETED EMCEE", output=of_run)

elapsed = time.time() - start
elapsed_d, elapsed_h, elapsed_m, elapsed_s = anc.computation_time(elapsed)

anc.print_both("", output=of_run)
anc.print_both(
    " pyTRADES: PyDE/EMCEE FINISHED in {0:2d} day {1:02d} hour {2:02d} min {3:.2f} sec - bye bye".format(
        int(elapsed_d), int(elapsed_h), int(elapsed_m), elapsed_s
    ), output=of_run
)

anc.print_both("", output=of_run)
of_run.close()
sim.reset()

# return

# ==============================================================================
# ==============================================================================

# if __name__ == "__main__":
#   main()
