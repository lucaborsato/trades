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

from pytrades.constants import Mjups, Msear, huge

# from pytrades_lib import f90trades
from pytrades import pytrades
from pytrades import ancillary as anc
from pytrades.pyde import DiffEvol
from pytrades import de


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
    m_type=cli.m_type,
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

    return lgllhd + lnprior


def lnprob_real(fitting_parameters):
    lnP = lnprob(fitting_parameters)
    if np.isinf(lnP):
        lnP = -huge
    return lnP


def minimize_func(fitting_parameters):
    lnL = lnprob(fitting_parameters)
    to_min = -2.0 * lnL
    if np.isinf(to_min):
        to_min = huge

    return to_min

def neg_lnprob(p):
    lnP = lnprob_real(p)
    return -lnP

def residuals_func(fitting_parameters):
    to_min = minimize_func(fitting_parameters)
    resw = np.zeros(ndata) + to_min / ndata

    return resw


# INITIALISE SCRIPT FOLDER/LOG FILE
working_folder, _, of_run = init_folder(working_path, cli.sub_folder)
sys.stdout.flush()

anc.print_both("")
# anc.print_both("fitting parameters:")
# anc.print_both(" ".join([n for n in sim.fitting_names]))
anc.print_both(
    "{:>20s} {:>13s} {:>13s} {:>13s}".format("name", "parameter", "min", "max")
)
for n, p, bd in zip(sim.fitting_names, sim.fitting_parameters, sim.fitting_minmax):
    anc.print_both("{:20s} {:13.6f} {:13.6f} {:13.6f}".format(n, p, bd[0], bd[1]))
anc.print_both("")

initial_parameters = sim.fitting_parameters.copy()

(
    chi_square,
    reduced_chi_square,
    lgllhd,
    lnprior,
    ln_const,
    bic,
    check,
) = sim.run_and_get_stats_from_parameters(initial_parameters)
anc.print_both("\nStats for the initial set of parameters:", output=of_run)
anc.print_both("chi_square = {}".format(chi_square), output=of_run)
anc.print_both("reduced_chi_square = {}".format(reduced_chi_square), output=of_run)
anc.print_both("lgllhd = {}".format(lgllhd), output=of_run)
anc.print_both("lnprior = {}".format(lnprior), output=of_run)
anc.print_both("ln_const = {}".format(ln_const), output=of_run)
anc.print_both("bic = {}".format(bic), output=of_run)
anc.print_both("check = {}".format(check), output=of_run)
anc.print_both("", output=of_run)
sys.stdout.flush()

fitting_parameters = sim.fitting_parameters.copy()

# check pyde args
anc.print_both("\n Check pyDE configuration ... ", output=of_run)
de_run = cli.de_type.lower() == "run"
de_resume = cli.de_type.lower() == "resume"
de_to_emcee = cli.de_type.lower() == "to_emcee"

# SET INITIAL PARAMETERS FOR EMCEE

# PYDE
if de_run or de_resume or de_to_emcee:
    (
        de_par, de_lnP, de_pop_best, de_fit_best, de_pop, de_fit
    ) = de.run_de(
        cli, sim, lnprob, working_folder, of_run=of_run
    )
    fitting_parameters = de_par.copy()


elif cli.trades_previous is not None:
    temp_names, fitting_parameters = anc.read_fitted_file(cli.trades_previous)
    if sim.nfit != np.shape(fitting_parameters)[0]:
        anc.print_both(
            " NUMBER OF PARAMETERS (%d) IN TRADES-PREVIOUS FILE DOES NOT"
            " MATCH THE CURRENT CONFIGURATION nfit=%d\nSTOP"
            % (np.shape(fitting_parameters)[0], sim.nfit),
            output=of_run,
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
        """,
        output=of_run,
    )
    fitting_parameters = initial_parameters.copy()

sys.stdout.flush()

for ifit, fitn in enumerate(sim.fitting_names):
    if fitn[0] == "w" or fitn[0:2] == "mA" or fitn[0:2] == "lN" or "lambda" in fitn:
        p = fitting_parameters[ifit]
        if (p < 0.0) or (p > 360.0):
            fitting_parameters[ifit] = p % 360.0

# save initial_fitting parameters into array
# initial_parameters = fitting_parameters.copy()

(
    chi_square,
    reduced_chi_square,
    lgllhd,
    lnprior,
    ln_const,
    bic,
    check,
) = sim.run_and_get_stats_from_parameters(fitting_parameters)
anc.print_both("\nStats for the fitting set of parameters:", output=of_run)
anc.print_both("chi_square = {}".format(chi_square), output=of_run)
anc.print_both("reduced_chi_square = {}".format(reduced_chi_square), output=of_run)
anc.print_both("lgllhd = {}".format(lgllhd), output=of_run)
anc.print_both("lnprior = {}".format(lnprior), output=of_run)
anc.print_both("ln_const = {}".format(ln_const), output=of_run)
anc.print_both("bic = {}".format(bic), output=of_run)
anc.print_both("check = {}".format(check), output=of_run)
anc.print_both("", output=of_run)
sys.stdout.flush()

# SET EMCEE PARAMETERS:
# nwalkers, nruns, nsave, _ = get_emcee_arguments(cli, nfit)

anc.print_both(
    "{:>20s} {:>13s} {:>13s} {:>13s}".format("name", "parameter", "min", "max")
)
for n, p, bd in zip(sim.fitting_names, fitting_parameters, sim.fitting_minmax):
    anc.print_both("{:20s} {:13.6f} {:13.6f} {:13.6f}".format(n, p, bd[0], bd[1]))
anc.print_both("")

if cli.nruns > 0:
    anc.print_both(
        "emcee chain: nwalkers = %d nruns = %d" % (cli.nwalkers, cli.nruns),
        output=of_run,
    )
    if cli.nwalkers < 2 * sim.nfit:
        cli.nwalkers = 2 * sim.nfit
        anc.print_both(
            "WARNING: nwalkers < 2 x nfit ==> setting nalkers = 2 x nfit  = {}".format(
                cli.nwalkers
            ),
            output=of_run,
        )
    sys.stdout.flush()

    # anc.print_both(" sampler ... ")

    backend_filename = os.path.join(working_folder, "emcee_summary.hdf5")
    backend = emcee.backends.HDFBackend(backend_filename, compression="gzip")
    if os.path.exists(backend_filename):
        completed_steps = backend.iteration
    else:
        completed_steps = 0

    if cli.emcee_restart:
        anc.print_both("keyword emcee_restart is True ==> resetting file content!")
        backend.reset(cli.nwalkers, sim.nfit)

    anc.print_both("Initial size of backend: {}".format(completed_steps), output=of_run)
    if completed_steps > 0:
        # p0 = None
        p0 = backend.get_last_sample()
        anc.print_both("continue emcee analysis ...", output=of_run)
    else:
        anc.print_both("new run of emcee analysis ...", output=of_run)
        if cli.pre_optimise:
            anc.print_both(
                "\n-- Optimise fitting parameters with Minimize(Nelder-Mead) ...",
                output=of_run,
            )

            opt_res = sopt.minimize(
                neg_lnprob,
                fitting_parameters,
                method="Nelder-Mead",
                bounds=[(bd[0], bd[1]) for bd in sim.fitting_minmax],
            )
            anc.print_both("Opt Message: {}".format(opt_res.message), output=of_run)
            anc.print_both("Opt Success: {}".format(opt_res.success), output=of_run)
            if opt_res.success:
                fitting_parameters = opt_res.x
                anc.print_both("\nStats for the fitting set of parameters:", output=of_run)
                anc.print_both("chi_square = {}".format(chi_square), output=of_run)
                anc.print_both("reduced_chi_square = {}".format(reduced_chi_square), output=of_run)
                anc.print_both("lgllhd = {}".format(lgllhd), output=of_run)
                anc.print_both("lnprior = {}".format(lnprior), output=of_run)
                anc.print_both("ln_const = {}".format(ln_const), output=of_run)
                anc.print_both("bic = {}".format(bic), output=of_run)
                anc.print_both("check = {}".format(check), output=of_run)
                anc.print_both("", output=of_run)
                anc.print_both(
                    "{:>20s} {:>13s} {:>13s} {:>13s}".format("name", "parameter", "min", "max")
                )
                for n, p, bd in zip(sim.fitting_names, fitting_parameters, sim.fitting_minmax):
                    anc.print_both("{:20s} {:13.6f} {:13.6f} {:13.6f}".format(n, p, bd[0], bd[1]))
                anc.print_both("")
                sys.stdout.flush()
            else:
                anc.print_both("Not changing fitting parameters")
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

        pka = {
            'ncols': 75,
            'dynamic_ncols': False,
            'position': 0
        }
        sampler.run_mcmc(
            p0,
            cli.nruns,
            thin_by=cli.thin_by,
            tune=True,  # default False
            skip_initial_state_check=False,  # default False
            progress=cli.emcee_progress,
            progress_kwargs=pka,
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
    ),
    output=of_run,
)

anc.print_both("", output=of_run)
of_run.close()
sim.reset()

# return

# ==============================================================================
# ==============================================================================

# if __name__ == "__main__":
#   main()
