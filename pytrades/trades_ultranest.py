#!/usr/bin/env python
# -*- coding: utf-8 -*-

# no more "zero" integer division bugs!:P
import argparse
import os
import numpy as np  # array

# import h5py
import sys
import time

# import glob

# import multiprocessing as mp
# from multiprocessing import Pool


import ultranest
from ultranest import stepsampler

from constants import Mjups, Msear
import ancillary as anc
# from pytrades_lib import f90trades
import pytrades

# =============================================================================

def mpi_print(l, rank, log=None):

    if rank == 0:
        anc.print_both(l, output=log)

    return


try:
    from mpi4py import MPI
    comm = MPI.COMM_WORLD
    rank = comm.Get_rank()
    nthreads = comm.Get_size()
except:
    comm = None
    rank = 0
    nthreads = 1
    
# from mpi4py import MPI
# comm = MPI.COMM_WORLD
# rank = comm.Get_rank()
# nthreads = comm.Get_size()
print("Hi! I am thread {} out of {}".format(rank, nthreads))
# =============================================================================

# sys.exit()

init_folder = anc.init_folder
compute_proper_sigma = anc.compute_proper_sigma
compute_initial_walkers = anc.compute_initial_walkers

# ==============================================================================


# =============================================================================
# =============================================================================
# MAIN SCRIPT - NOT IN FUNCTION DUE TO ISSUE WITH PARALLEL AND PICKLE FUNCTION OF LNPROB..
# =============================================================================
# =============================================================================

# MAIN -- TRADES + ULTRANEST
# READ COMMAND LINE ARGUMENTS
yml_file = anc.get_input_file()
cli = anc.ConfigurationRun(yml_file)

os.environ["OMP_NUM_THREADS"] = "1"
# os.environ["HDF5_USE_FILE_LOCKING"] = "FALSE"


# STARTING TIME
start = time.time()

# RENAME
working_path = cli.full_path
# nthreads = cli.nthreads
np.random.seed(cli.seed)

# INITIALISE TRADES WITH SUBROUTINE WITHIN TRADES_LIB -> PARAMETER NAMES, MINMAX, INTEGRATION ARGS, READ DATA ...
# pytrades.initialize_trades(working_path, cli.sub_folder, nthreads)
mpi_print("Init trades ... ", rank)
if rank == 0:
    sim = pytrades.TRADES(
        working_path,
        sub_folder=cli.sub_folder, 
        nthreads=cli.nthreads,
        seed=cli.seed,
        m_type=cli.m_type
    )
    sim.init_trades()
sys.stdout.flush()


#
# FUNCTION NEEDED BY ULTRANEST
#
def my_prior_transform(cube):

    params = np.array(cube.copy())
    bounds = sim.fitting_minmax
    db = bounds[:, 1] - bounds[:, 0]
    params = bounds[:, 0] + np.array(cube) * db

    return params


def lnL_ultranest(fitting_parameters):

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
        lnL,
        lnp,
        _,
        _,
        check,
    ) = sim.run_and_get_stats_from_parameters(fitting_parameters)

    return lnL + lnp


# INITIALISE SCRIPT FOLDER/LOG FILE
working_folder, _, of_run = init_folder(working_path, cli.sub_folder)
sys.stdout.flush()



mpi_print("ULTRANEST", rank, log=of_run)
resume_flag = ["resume", "resume-similar", "overwrite", "subfolder"]
if cli.resume_flag not in resume_flag:
    mpi_print(
        "ERROR: ultranest resume flag not available. please select among: {}".format(
            resume_flag
        ), rank, log=of_run
    )
    sys.exit()

wrapped_pars = anc.check_wrapped_parameters(sim.fitting_names)

nfit_min = 2 * sim.nfit
if cli.live_points < nfit_min:
    n_live_points = nfit_min
    anc.print_both("set n_live_points = {}".format(n_live_points))
else:
    n_live_points = cli.live_points


sampler = ultranest.ReactiveNestedSampler(
    sim.fitting_names,
    lnL_ultranest,
    my_prior_transform,
    log_dir=working_folder,  # folder where to store files
    resume=cli.resume_flag,  # whether to resume from there (otherwise start from scratch)
    # vectorized=True, # NOT WORKING WITH TRADES ...
    storage_backend="hdf5",
    wrapped_params=wrapped_pars,
    draw_multiple=True,
)

nsteps = 2 * sim.nfit
# create step sampler:
sampler.stepsampler = stepsampler.RegionSliceSampler(
    nsteps=nsteps, adaptive_nsteps="move-distance"
)

sys.stdout.flush()
result = sampler.run(
    min_num_live_points = n_live_points,
    dlogz               = cli.dlogz,  # desired accuracy on logz
    # min_ess=n_live_points,  # number of effective samples
    # update_interval_iter_fraction=0.4, # how often to update region !!NOT HERE?!!
    # max_num_improvement_loops=3,  # how many times to go back and improve
)

elapsed = time.time() - start
elapsed_d, elapsed_h, elapsed_m, elapsed_s = anc.computation_time(elapsed)

mpi_print("COMPLETED ULTRANEST", rank, log=of_run)
sys.stdout.flush()

mpi_print("", rank, log=of_run)
mpi_print(
    " pyTRADES: ULTRANEST FINISHED in {0:2d} day {1:02d} hour {2:02d} min {3:.2f} sec - bye bye".format(
        int(elapsed_d), int(elapsed_h), int(elapsed_m), elapsed_s
    ),
    rank,
    log=of_run,
)
mpi_print("", rank, log=of_run)
sys.stdout.flush()

sampler.print_results()
sys.stdout.flush()

sampler.plot_run()
sampler.plot_trace()
sampler.plot_corner()

mpi_print("", rank, log=of_run)
of_run.close()
sim.reset()

# return

# ==============================================================================
# ==============================================================================

# if __name__ == "__main__":
#   main()
