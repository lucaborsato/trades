#!/usr/bin/env python
# -*- coding: utf-8 -*-

# no more "zero" integer division bugs!:P
import argparse
from logging.handlers import DEFAULT_SOAP_LOGGING_PORT
import os
from random import seed
import numpy as np  # array
import h5py
import sys
import time

# import glob

# import multiprocessing as mp
from multiprocessing import Pool

# from schwimmbad import JoblibPool as Pool

from pyde.de import DiffEvol
import emcee
import yaml

from constants import Msear
import ancillary as anc
import de_plot as dep
from pytrades_lib import f90trades

# =============================================================================


def get_args():
    parser = argparse.ArgumentParser(description="TRADES+PyDE+EMCEE")

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

    parser.add_argument(
        "-seed",
        "--seed",
        action="store",
        dest="seed",
        default="None",
        help="Seed for random number generator. Default is None.",
    )

    parser.add_argument(
        "-de",
        "--de",
        "-de-type",
        "--de-type",
        action="store",
        dest="de_type",
        default="run",
        help="""
            DE run type:
            "run" (Default) it run a DE and overwrites files;
            "to_emcee" reads previous de_run.hdf5 file and load the best-fit parameters and provides them to emcee;
            "resume" reads old hdf5 file with DE evolution and start from last iteration
        """,
    )

    cli = parser.parse_args()

    cli.full_path = os.path.join(os.path.abspath(cli.full_path), "")
    cli.sub_folder = os.path.join(os.path.relpath(cli.sub_folder), "")

    cli.nthreads = int(cli.nthreads)

    cli.nwalkers = int(cli.nwalkers)
    cli.nruns = int(cli.nruns)

    cli.ln_flag = anc.set_bool_argument(cli.ln_flag)

    cli.emcee_progress = anc.set_bool_argument(cli.emcee_progress)

    cli.de_type = str(cli.de_type).lower()

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

# compute_proper_sigma = anc.compute_proper_sigma
compute_initial_walkers = anc.compute_initial_walkers

# =============================================================================


def load_de_configuration(cli):

    de_file = os.path.join(cli.full_path, "de.yml")

    de_conf_default = {
        "npop": 84,
        "ngen": 4200,
        "n_global": 1,
        "save": 100,
        "f": 0.5,
        "c": 0.5,
        "maximize": True,
    }

    de_yml = de_conf_default.copy()
    de_labels = de_conf_default.keys()
    with open(de_file, "r") as stream:
        try:
            de_yml_read = yaml.safe_load(stream)
            for k, v in de_yml_read.items():
                if k in de_labels:
                    de_yml[k] = v
        except yaml.YAMLError as exc:
            print("de.yml not present, using defaults.")

    return de_yml


# =============================================================================

def hdf5_save_one_dataset(de_file, data_name, data, data_type, hdf5_mode="a"):

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


def hdf5_update_attr(de_file, data_name, attrs_name, attrs_value):

    de_hdf5 = h5py.File(
        de_file,
        "r+",
        # libver='latest'
    )
    # de_hdf5.swmr_mode = True
    de_hdf5[data_name].attrs[attrs_name] = attrs_value
    de_hdf5.close()

    return


def save_de_evolution(
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
    de_yml,
):

    # de_file = os.path.join(de_path, 'de_run.hdf5')
    hdf5_save_one_dataset(de_file, "population", de_pop, "f8")
    hdf5_update_attr(de_file, "population", "npop", npop_de)
    hdf5_update_attr(de_file, "population", "ngen", ngen_de)
    hdf5_update_attr(de_file, "population", "iter_de", iter_de)
    hdf5_update_attr(de_file, "population", "iter_global", iter_global + 1)
    hdf5_update_attr(de_file, "population", "nfit", nfit)
    hdf5_update_attr(de_file, "population", "ndata", ndata)

    hdf5_save_one_dataset(de_file, "population_fitness", de_fit, "f8")

    hdf5_save_one_dataset(de_file, "best_population", de_pop_best, "f8")

    hdf5_save_one_dataset(de_file, "best_fitness", de_fit_best, "f8")

    hdf5_save_one_dataset(de_file, "parameters_minmax", de_bounds, "f8")

    # best_loc = np.argmax(de_fit_best[:iter_de])
    # hdf5_save_one_dataset(de_file, 'de_parameters', de_pop_best[best_loc,:].copy(), 'f8')
    de_parameters = get_best_de_parameters(de_fit_best, de_pop_best, iter_de, de_yml)
    hdf5_save_one_dataset(de_file, "de_parameters", de_parameters, "f8")

    hdf5_save_one_dataset(de_file, "parameter_names", parameter_names, "S10")

    return


# =============================================================================


def get_best_de_parameters(de_fit_best, de_pop_best, iter_de, de_yml):

    de_maximize = de_yml["maximize"]
    loc_best = (
        np.argmax(de_fit_best[: iter_de + 1])
        if de_maximize
        else np.argmin(de_fit_best[: iter_de + 1])
    )

    de_parameters = de_pop_best[loc_best, :].copy()

    return de_parameters


# =============================================================================
def load_de_parameters(de_file):

    de_hdf5 = h5py.File(
        de_file,
        "r",
        # libver='latest', swmr=True
    )
    de_parameters = de_hdf5["de_parameters"][...]
    de_hdf5.close()

    return de_parameters


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

def lnprob(fitting_parameters):

    # loglhd, logprior = 0.0, 0.0
    # check = 1
    loglhd, logprior, check = pytrades.fortran_logprob(
        # np.asarray(fitting_parameters, dtype=float)
        fitting_parameters
    )
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


# INITIALISE DE ARGUMENTS FROM de.yml FILE
de_yml = load_de_configuration(cli)
anc.print_both("\n DE configuration: ", of_run)
for k, v in de_yml.items():
    anc.print_both(" {:10s} = {}".format(k, v))

de_run = cli.de_type.lower() == "run"
de_resume = cli.de_type.lower() == "resume"
de_to_emcee = cli.de_type.lower() == "to_emcee"
npop_de = de_yml["npop"]
ngen_de = de_yml["ngen"]
n_global = de_yml["n_global"]
de_save = de_yml["save"]
de_bounds = parameters_minmax.copy()
de_f = de_yml["f"]
de_c = de_yml["c"]
de_maximize = de_yml["maximize"]
fit_type = -1 if de_maximize else 1
anc.print_both(
    " DE n_global = {} npop = {} ngen = {}".format(n_global, npop_de, ngen_de), of_run
)
anc.print_both(
    " DE de_type = {} ==> de_run = {}, de_resume = {}, de_to_emcee = {}".format(
        cli.de_type.lower(), de_run, de_resume, de_to_emcee
    )
)

# RUN DE+EMCEE n_global TIMES
for iter_global in range(0, n_global):

    # CREATES PROPER WORKING PATH AND NAME
    i_global = iter_global + 1
    de_path = os.path.join(
        os.path.join(working_folder, "{0:04d}_de2emcee".format(i_global)), ""
    )
    pytrades.path_change(de_path)

    anc.print_both(
        "\n\n GLOBAL RUN {0:04d} INTO PATH: {1:s}\n".format(i_global, de_path), of_run
    )

    if nthreads > 1:
        threads_pool = Pool(nthreads)
        anc.print_both(" using Pool with {} threads".format(nthreads), of_run)
    else:
        threads_pool = None
        anc.print_both(" serial mode", of_run)

    if de_run:

        # ==============================================================================
        # START DE
        anc.print_both(" RUN DE ...", of_run)

        de_start = time.time()
        if not os.path.exists(de_path):
            os.makedirs(de_path)
        # copy files
        anc.copy_simulation_files(working_path, de_path)

        de_file = os.path.join(de_path, "de_run.hdf5")
        if os.path.exists(de_file) and os.path.isfile(de_file):
            anc.print_both("File {} exists: deleting it!".format(de_file), of_run)
            os.remove(de_file)

        de_parameters = trades_parameters.copy()
        de_fitness = 0.0

        de_pop, de_fit = (
            np.zeros((ngen_de, npop_de, nfit)),
            np.zeros((ngen_de, npop_de)),
        )
        de_pop_best, de_fit_best = np.zeros((ngen_de, nfit)), np.zeros((ngen_de))

        de_evol = DiffEvol(
            # lnprob_sq,
            lnprob,
            de_bounds,
            npop_de,
            f=de_f,
            c=de_c,
            seed=cli.seed,
            maximize=de_maximize,
            vectorize=False,
            pool=threads_pool,
            # args=(parameter_names,)
        )

        start_iter_de = time.time()
        anc.print_both(" DE - START ")

        sys.stdout.flush()
        for iter_de, res_de in enumerate(de_evol(ngen_de)):
            
            de_pop[iter_de, :, :] = de_evol.population.copy()
            de_fit[iter_de, :] = fit_type * de_evol._fitness.copy()
            de_pop_best[iter_de, :] = de_evol.minimum_location.copy()
            de_fit_best[iter_de] = fit_type * de_evol.minimum_value

            if iter_de > 0:
                if ((iter_de + 1) % de_save == 0) or (iter_de + 1 >= ngen_de):
                    anc.print_both(" DE - iter = {} ==> saving".format(iter_de))
                    anc.print_both(
                        " last best fitness = {}".format(de_fit_best[iter_de]), of_run
                    )

                    # SAVE DE SIMULATION IN de_run.hdf5 FILE
                    save_de_evolution(
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
                        de_yml,
                    )
                    print(" Updated DE hdf5 file: {}".format(de_file))
                    elapsed_de = time.time() - start_iter_de
                    (
                        elapsed_de_d,
                        elapsed_de_h,
                        elapsed_de_m,
                        elapsed_de_s,
                    ) = anc.computation_time(elapsed_de)
                    anc.print_both(
                        " DE: {:d} steps in {:2d} day {:02d} hour {:02d} min {:.2f} sec".format(
                            de_save,
                            int(elapsed_de_d),
                            int(elapsed_de_h),
                            int(elapsed_de_m),
                            elapsed_de_s,
                        ),
                        of_run,
                    )
                    sys.stdout.flush()
                    start_iter_de = time.time()

        print()
        anc.print_both(" completed DE", of_run)
        anc.print_both(" ", of_run)

        elapsed = time.time() - de_start
        elapsed_d, elapsed_h, elapsed_m, elapsed_s = anc.computation_time(elapsed)
        anc.print_both(" ", of_run)
        anc.print_both(
            " DE FINISHED in {0:2d} day {1:02d} hour {2:02d} min {3:.2f} sec - bye bye".format(
                int(elapsed_d), int(elapsed_h), int(elapsed_m), elapsed_s
            ),
            of_run,
        )
        sys.stdout.flush()
        # best_loc = np.argmax(de_fit_best)
        # de_parameters = de_pop_best[best_loc,:].copy()
        de_parameters = get_best_de_parameters(
            de_fit_best, de_pop_best, ngen_de-1, de_yml
        )
        # END DE
        # ==============================================================================
        # sys.exit()
    elif (de_resume) and (os.path.isdir(de_path)):
        # ==============================================================================
        # START DE

        anc.print_both(
            "DE RUN already exist --> loading and resume DE ...",
            of_run,
        )
        anc.print_both(" RESUME DE ...", of_run)

        de_start = time.time()

        de_old = dep.DEEvolution(de_path)
        de_file = os.path.join(de_path, "de_run.hdf5")

        # CALL RUN_PSO SUBROUTINE FROM TRADES_LIB: RUNS PSO AND COMPUTES THE BEST SOLUTION, SAVING ALL THE POPULATION EVOLUTION
        de_parameters = trades_parameters.copy()

        if de_old.ngen_de == ngen_de:
            de_parameters = load_de_parameters(de_file)
        else:
            de_fitness = 0.0

            de_pop, de_fit = (
                np.zeros((ngen_de, npop_de, nfit)),
                np.zeros((ngen_de, npop_de)),
            )
            de_pop_best, de_fit_best = np.zeros((ngen_de, nfit)), np.zeros((ngen_de))

            # update de_pop, de_fit, de_pop_best, de_fit_best
            de_pop[: de_old.ngen_de, :, :] = de_old.de_pop.copy()
            de_fit[: de_old.ngen_de, :] = de_old.de_fit.copy()
            de_pop_best[: de_old.ngen_de, :] = de_old.de_pop_best.copy()
            de_fit_best[: de_old.ngen_de] = de_old.de_fit_best.copy()
            last_iter_de = de_old.iter_de

            de_evol = DiffEvol(
                lnprob,
                de_bounds,
                npop_de,
                f=de_f,
                c=de_c,
                seed=cli.seed,
                maximize=de_maximize,
                vectorize=False,
                pool=threads_pool,
                # args=(parameter_names,)
            )

            de_evol._population = de_old.de_pop[-1, :, :].copy()
            de_evol._fitness = de_old.de_fit[-1, :].copy()
            de_evol._minidx = (
                np.argmax(de_evol._fitness)
                if de_maximize
                else np.argmin(de_evol._fitness)
            )
            if last_iter_de + 1 < ngen_de:
                start_iter_de = time.time()
                anc.print_both(
                    " DE - START: CONTINUE FROM ITER = {}".format(last_iter_de)
                )
                for iter_de_temp, res_de in enumerate(de_evol(ngen_de - last_iter_de)):
                    # current iter_de_temp starts from 0 and ends to ngen_de-last_iter_de
                    iter_de = last_iter_de + iter_de_temp
                    # iter_de goes from last_iter_de to ngen_de

                    de_pop[iter_de, :, :] = de_evol.population.copy()
                    de_fit[iter_de, :] = fit_type * de_evol._fitness.copy()
                    de_pop_best[iter_de, :] = de_evol.minimum_location.copy()
                    de_fit_best[iter_de] = fit_type * de_evol.minimum_value

                    if iter_de > 0:
                        if ((iter_de + 1) % de_save == 0) or (iter_de + 1 >= ngen_de):
                            anc.print_both(" DE - iter = {} ==> saving".format(iter_de))
                            anc.print_both(
                                " last best fitness = {}".format(de_fit_best[iter_de]),
                                of_run,
                            )

                            # SAVE DE SIMULATION IN de_run.hdf5 FILE
                            save_de_evolution(
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
                                de_yml,
                            )
                            print(" Updated DE hdf5 file: {}".format(de_file))
                            elapsed_de = time.time() - start_iter_de
                            (
                                elapsed_de_d,
                                elapsed_de_h,
                                elapsed_de_m,
                                elapsed_de_s,
                            ) = anc.computation_time(elapsed_de)
                            anc.print_both(
                                " DE: {:d} steps in {:2d} day {:02d} hour {:02d} min {:.2f} sec".format(
                                    de_save,
                                    int(elapsed_de_d),
                                    int(elapsed_de_h),
                                    int(elapsed_de_m),
                                    elapsed_de_s,
                                ),
                                of_run,
                            )
                            sys.stdout.flush()
                            start_iter_de = time.time()

                anc.print_both(" completed DE", of_run)
                anc.print_both(" ", of_run)

                elapsed = time.time() - de_start
                elapsed_d, elapsed_h, elapsed_m, elapsed_s = anc.computation_time(
                    elapsed
                )
                anc.print_both(" ", of_run)
                anc.print_both(
                    " DE FINISHED in {0:2d} day {1:02d} hour {2:02d} min {3:.2f} sec - bye bye".format(
                        int(elapsed_d), int(elapsed_h), int(elapsed_m), elapsed_s
                    ),
                    of_run,
                )
                sys.stdout.flush()
            else:
                iter_de = last_iter_de
            # best_loc = np.argmax(de_fit_best)
            # de_parameters = de_pop_best[best_loc,:].copy()
            de_parameters = get_best_de_parameters(
                de_fit_best, de_pop_best, iter_de, de_yml
            )
        # END DE
        # ==============================================================================

    elif (de_to_emcee) and (os.path.isdir(de_path)):
        anc.print_both(
            "DE RUN already exist --> loading best-fit parameters and send them to emcee ...",
            of_run,
        )
        de_parameters = load_de_parameters(de_path)

    else:
        anc.print_both(
            "de_type is not run or resume or to_emcee. Initial parameters directly to emcee",
            of_run,
        )
        de_parameters = trades_parameters.copy()


    sys.stdout.flush()

    # TRADES/PSO USES ECOSW/ESINW --> HERE EMCEE USES SQRTECOSW/SQRTESINW
    # emcee_parameters = anc.e_to_sqrte_fitting(pso_parameters, trades_names)
    emcee_parameters = de_parameters.copy()

    if float(cli.delta_sigma) <= 0.0:
        delta_sigma = 1.0e-4
    delta_sigma = cli.delta_sigma
    p0 = compute_initial_walkers(
        # lnprob_sq,
        lnprob,
        nfit,
        nwalkers,
        emcee_parameters,
        parameters_minmax,
        delta_sigma,
        parameter_names,
        args=(),
        of_run=of_run,
    )


    anc.print_both(
        " emcee chain: nwalkers = {} nruns = {}".format(nwalkers, nruns), of_run
    )
    anc.print_both(" sampler ... ", of_run)

    # if(nthreads > 1):
    #     threads_pool = Pool(nthreads)
    # else:
    #     threads_pool = None

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
    )

    anc.print_both(" ready to go", of_run)
    anc.print_both(" with nsave = {}".format(nsave), of_run)
    sys.stdout.flush()

    if nsave != False:
        # save temporary sampling during emcee every nruns*10%
        if os.path.exists(
            os.path.join(de_path, "emcee_summary.hdf5")
        ) and os.path.isfile(os.path.join(de_path, "emcee_summary.hdf5")):
            os.remove(os.path.join(de_path, "emcee_summary.hdf5"))

        f_hdf5 = h5py.File(os.path.join(de_path, "emcee_summary.hdf5"), "a")
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

            f_hdf5 = h5py.File(os.path.join(de_path, "emcee_summary.hdf5"), "a")
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
        f_hdf5 = h5py.File(os.path.join(de_path, "emcee_summary.hdf5"), "w")
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

# if __name__ == "__main__":
#   main()
