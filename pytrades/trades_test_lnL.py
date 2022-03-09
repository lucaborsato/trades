#!/usr/bin/env python
# -*- coding: utf-8 -*-

# no more "zero" integer division bugs!:P
import argparse
import os
import numpy as np  # array
import sys
import time

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
    # cli.nruns = int(cli.nruns)

    # cli.emcee_previous = anc.set_adhoc_file(cli.emcee_previous)
    # cli.trades_previous = anc.set_adhoc_file(cli.trades_previous)

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
parameter_names = anc.trades_names_to_emcee(trades_names)

# INITIAL PARAMETER SET (NEEDED ONLY TO HAVE THE PROPER ARRAY/VECTOR)
trades_parameters = pytrades.fitting_parameters

# save initial_fitting parameters into array
original_fit_parameters = trades_parameters.copy()
fitting_parameters = anc.e_to_sqrte_fitting(trades_parameters, trades_names)

trades_minmax = pytrades.parameters_minmax  # PARAMETER BOUNDARIES
parameters_minmax = anc.e_to_sqrte_boundaries(trades_minmax, trades_names)

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


def lnprob_sq(fitting_parameters, names_par):

    fitting_trades = anc.sqrte_to_e_fitting(fitting_parameters, names_par)
    loglhd = lnprob(fitting_trades)
    # if np.isfinite(loglhd):
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


n_sim = 10
elapsed_times = np.zeros((n_sim))
elapsed_start, elapsed_end = 0.0, 0.0
anc.print_both(" CHECKING TIMING FOR n_sim = {}".format(n_sim), of_run)
for i_sim in range(0,n_sim):
    elapsed_start = time.time()
    _ = lnL_ultranest(fitting_parameters)
    elapsed_end = time.time()
    elapsed_times[i_sim] = elapsed_end - elapsed_start
    anc.print_both("Elapsed time for sim number {} : {:6.2f}".format(i_sim, elapsed_times[i_sim]),
        of_run
    )
mean_ela = np.mean(elapsed_times)
std_ela = np.std(elapsed_times, ddof=1)/np.sqrt(n_sim)
anc.print_both("Elapsed time for {} simulations: {:6.2f}".format(n_sim, np.sum(elapsed_times)),
    of_run
)
anc.print_both("Mean Elapsed time for 1 simulation: {:6.2f} +\- {:6.2f}".format(mean_ela, std_ela),
    of_run
)

anc.print_both("", of_run)
of_run.close()
pytrades.deallocate_variables()

# return

# ==============================================================================
# ==============================================================================

# if __name__ == "__main__":
#   main()
