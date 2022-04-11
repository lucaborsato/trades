#!/usr/bin/env python
# -*- coding: utf-8 -*-

# no more "zero" integer division bugs!:P
import argparse
import os
import numpy as np  # array
import h5py
import sys
import time
import glob

# import multiprocessing as mp
from multiprocessing import Pool

# from schwimmbad import JoblibPool as Pool

import ultranest
from ultranest import stepsampler

from constants import Mjups, Msear
import ancillary as anc
from pytrades_lib import pytrades

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

    # ULTRANEST

    # NUMBER OF WALKERS TO USE WITH EMCEE
    parser.add_argument(
        "--nl",
        "--n-live-points",
        "--n-live",
        "--n-points",
        "--npop",
        action="store",
        dest="n_live_points",
        default=1,
        help="Number of live points (or number of chains) to use. default n_live_points = nfit x 2",
    )

    # resume flag: 'resume', 'resume-similar', 'overwrite' or 'subfolder'
    parser.add_argument(
        "-r",
        "--resume-flag",
        action="store",
        dest="resume_flag",
        default='resume',
        help="Ultranest resume flag: 'resume', 'resume-similar', 'overwrite' or 'subfolder'",
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

    cli.n_live_points = int(cli.n_live_points)

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

# ==============================================================================


# =============================================================================
# =============================================================================
# MAIN SCRIPT - NOT IN FUNCTION DUE TO ISSUE WITH PARALLEL AND PICKLE FUNCTION OF LNPROB..
# =============================================================================
# =============================================================================

# def main():

# MAIN -- TRADES + ULTRANEST
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
# parameter_names = anc.trades_names_to_emcee(trades_names)
parameter_names = trades_names.copy()

# INITIAL PARAMETER SET (NEEDED ONLY TO HAVE THE PROPER ARRAY/VECTOR)
trades_parameters = pytrades.fitting_parameters

# save initial_fitting parameters into array
original_fit_parameters = trades_parameters.copy()
# fitting_parameters = anc.e_to_sqrte_fitting(trades_parameters, trades_names)
fitting_parameters = trades_parameters.copy()

trades_minmax = pytrades.parameters_minmax  # PARAMETER BOUNDARIES
# parameters_minmax = anc.e_to_sqrte_boundaries(trades_minmax, trades_names)
parameters_minmax = trades_minmax.copy()

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

# uses the ancillary.get_fitted(full_path)
# to obtain info for the conversion from fitted to physical parameters
# nfit, NB, bodies_file, id_fit, id_all, nfit_list, cols_list, case =  anc.get_fitted(working_path)
_, _, _, id_fit, _, _, cols_list, case = anc.get_fitted(working_path)
# stellar mass in solar unit to earth, priors m in Earth masses
m_factor = pytrades.mr_star[0, 0] * Msear
# read priors file: priors.in
priors = anc.read_priors(working_path)
kep_elem = anc.all_parameters_to_kep_elem(
    pytrades.system_parameters, n_bodies
)
# lnL_priors(p, priors, names_par, kep_elem, id_fit, case_list, cols_list, m_factor)
ln_prior = anc.lnL_priors(
    fitting_parameters,
    priors,
    parameter_names,
    kep_elem,
    id_fit,
    case,
    cols_list,
    m_factor,
)

#
# FUNCTION NEEDED BY ULTRANEST
#
def my_prior_transform(cube):
    params = np.array(cube.copy())

    dminmax = parameters_minmax[:,1] - parameters_minmax[:,0]
    params = parameters_minmax[:,0] + np.array(cube) * dminmax

    return params

def lnprob(fitting_parameters):
    loglhd = 0.0
    # check = 1
    loglhd, _ = pytrades.fortran_loglikelihood(
        np.array(fitting_parameters, dtype=np.float64)
    )
    # not neede to add ln_err_const, because it is included in the loglhd now
    # loglhd = loglhd + ln_err_const # ln_err_const: global variable
    # if check == 0:
    #     loglhd = -np.inf
    return loglhd


# def lnprob_sq(fitting_parameters, names_par):
#     fitting_trades = anc.sqrte_to_e_fitting(fitting_parameters, names_par)
#     loglhd = lnprob(fitting_trades)
def lnprob_sq(fitting_parameters, parameter_names):
    loglhd = lnprob(fitting_parameters)
    ln_prior = anc.lnL_priors(
        fitting_parameters,
        priors,
        parameter_names,
        kep_elem,
        id_fit,
        case,
        cols_list,
        m_factor,
    )
    loglhd += ln_prior

    return loglhd

# def lnL_ultranest(cube_parameters):
#     fitting_parameters = my_prior_transform(cube_parameters)
def lnL_ultranest(fitting_parameters):
    loglhd = lnprob_sq(fitting_parameters, parameter_names)

    return loglhd

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


anc.print_both(" ORIGINAL PARAMETER VALUES -> 0000", of_run)
chi_square_0000, reduced_chi_square_0000, \
    lgllhd_0000, lnc_0000, bic_000, check_0000 = pytrades.write_summary_files(
    0, original_fit_parameters
)

anc.print_both(" ", of_run)
anc.print_both(" TESTING LNPROB_SQ ...", of_run)

lgllhd_zero = lnprob(trades_parameters)
lgllhd_sq_zero = lnprob_sq(fitting_parameters, parameter_names)

anc.print_both(" ", of_run)
anc.print_both(
    " %15s %23s %23s %15s %23s"
    % ("trades_names", "original_trades", "trades_par", "python_names", "python_par"),
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
sys.stdout.flush()

anc.print_both("ULTRANEST", of_run)
resume_flag = ['resume', 'resume-similar', 'overwrite' or 'subfolder']
if cli.resume_flag not in resume_flag:
    anc.print_both("ERROR: ultranest resume flag not available. please select among: {}".format(resume_flag))
    sys.exit()

wrapped_pars = anc.check_wrapped_parameters(parameter_names)

nfit_min = 64
if cli.n_live_points < nfit_min:
    n_live_points = nfit_min
    anc.print_both("set n_live_points = {}".format(n_live_points))
else:
    n_live_points = cli.n_live_points


sampler = ultranest.ReactiveNestedSampler(
    parameter_names,
    lnL_ultranest,
    my_prior_transform,
    log_dir=working_folder, # folder where to store files
    resume=cli.resume_flag, # whether to resume from there (otherwise start from scratch)
    # vectorized=True, # NOT WORKING WITH TRADES ... 
    storage_backend='hdf5',
    wrapped_params=wrapped_pars,
    draw_multiple=True,
)

nsteps = 2 * nfit
# create step sampler:
sampler.stepsampler = stepsampler.RegionSliceSampler(
    nsteps=nsteps,
    adaptive_nsteps='move-distance'
)

sys.stdout.flush()
result = sampler.run(
    min_num_live_points=n_live_points,
    dlogz=0.5, # desired accuracy on logz
    min_ess=n_live_points*2, # number of effective samples
    # update_interval_iter_fraction=0.4, # how often to update region !!NOT HERE?!!
    max_num_improvement_loops=3, # how many times to go back and improve
)

elapsed = time.time() - start
elapsed_d, elapsed_h, elapsed_m, elapsed_s = anc.computation_time(elapsed)

anc.print_both("COMPLETED ULTRANEST", of_run)
sys.stdout.flush()

sampler.print_results()
sys.stdout.flush()

sampler.plot_run()
sampler.plot_trace()
sampler.plot_corner()

anc.print_both("", of_run)
anc.print_both(
    " pyTRADES: ULTRANEST FINISHED in {0:2d} day {1:02d} hour {2:02d} min {3:.2f} sec - bye bye".format(
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
