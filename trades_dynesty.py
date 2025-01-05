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

import pickle

# import multiprocessing as mp
# from multiprocessing import Pool

import matplotlib.pyplot as plt

import dynesty
from dynesty.pool import Pool
from dynesty import plotting as dyplot

from pytrades.constants import Mjups, Msear
from pytrades import ancillary as anc
# from pytrades_lib import f90trades
from pytrades import pytrades

# =============================================================================
# =============================================================================

anc.set_rcParams()

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
sim = pytrades.TRADESfolder(
    working_path,
    sub_folder=cli.sub_folder, 
    nthreads=cli.nthreads,
    seed=cli.seed,
    m_type=cli.m_type
)
sim.init_trades()
sys.stdout.flush()


#
# FUNCTION NEEDED BY DYNESTY
#
def my_prior_transform(cube):

    params = np.array(cube.copy())
    bounds = sim.fitting_minmax
    db = bounds[:, 1] - bounds[:, 0]
    params = bounds[:, 0] + np.array(cube) * db

    return params


def lnL_dynesty(fitting_parameters):

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

dynesty_output_file = os.path.join(
    working_folder,
    "dynesty_results.pickle"
)
dyne_output_exists = os.path.exists(dynesty_output_file)
anc.print_both("Dynesty: just plots = {}".format(cli.dynesty_just_plots), output=of_run)
anc.print_both("Dynesty: output_file exists? {}".format(dyne_output_exists), output=of_run)

if cli.dynesty_just_plots and dyne_output_exists:
    with open(dynesty_output_file, "rb") as handle:
        results = pickle.load(handle)

else:
    wrapped_pars = anc.check_wrapped_parameters(sim.fitting_names)
    if np.any(wrapped_pars):
        wrapped_pars = np.arange(0,sim.nfit,1)[wrapped_pars]
    else:
        wrapped_pars = None

    nlive_min = 2 * sim.nfit
    if cli.dynesty_live_points < nlive_min:
        n_live_points = nlive_min
        anc.print_both("set n_live_points = {}".format(n_live_points), output=of_run)
    else:
        n_live_points = cli.dynesty_live_points

    dynesty_save_file = os.path.join(
        working_folder,
        "dynesty.save"
    )
    save_file_exists = os.path.exists(dynesty_save_file)
    resume = cli.dynesty_restore and save_file_exits
    anc.print_both("live_points = {}\ncheckout save file = {}\nresume = {}\nbound = {}\nsampler = {}\ndlogz = {}\npfrac = {}\n".format(
        n_live_points,
        dynesty_save_file,
        resume,
        cli.dynesty_bound,
        cli.dynesty_sample,
        cli.dynesty_dlogz,
        cli.dynesty_pfrac
        ), output=of_run
    )


    # PARALLEL
    if cli.nthreads > 1:


        with Pool(cli.nthreads, lnL_dynesty, my_prior_transform) as pool:
        # with Pool(cli.nthreads) as pool:

            if resume:
                sampler = dynesty.DynamicNestedSampler.restore(
                    dynesty_save_file,
                    pool=pool
                )
            else:
                sampler = dynesty.DynamicNestedSampler(
                    pool.loglike,
                    pool.prior_transform,
                    # lnL_dynesty,
                    # my_prior_transform,
                    sim.nfit,
                    nlive=n_live_points,
                    bound=cli.dynesty_bound,
                    sample=cli.dynesty_sample,
                    periodic=wrapped_pars,
                    pool=pool,
                    queue_size=cli.nthreads,
                    use_pool={'prior_transform': False}
                )
            sampler.run_nested(
                dlogz_init=cli.dynesty_dlogz,
                resume=resume,
                checkpoint_file=dynesty_save_file,
                wt_kwargs={'pfrac': cli.dynesty_pfrac},
            )

    # SERIAL
    else:
        if resume:
                sampler = dynesty.DynamicNestedSampler.restore(
                    dynesty_save_file
                )
        else:
            sampler = dynesty.DynamicNestedSampler(
                lnL_dynesty,
                my_prior_transform,
                sim.nfit,
                nlive=n_live_points,
                bound=cli.dynesty_bound,
                sample=cli.dynesty_sample,
                periodic=wrapped_pars,
            )

        sampler.run_nested(
            dlogz_init=cli.dynesty_dlogz,
            resume=resume,
            checkpoint_file=dynesty_save_file,
            wt_kwargs={'pfrac': cli.dynesty_pfrac},
        )


    sys.stdout.flush()


    elapsed = time.time() - start
    elapsed_d, elapsed_h, elapsed_m, elapsed_s = anc.computation_time(elapsed)

    sys.stdout.flush()

    results = sampler.results

    anc.print_both(
        " pyTRADES: DYNESTY FINISHED in {0:2d} day {1:02d} hour {2:02d} min {3:.2f} sec".format(
            int(elapsed_d), int(elapsed_h), int(elapsed_m), elapsed_s
        ),
        output=of_run,
    )
    sys.stdout.flush()

    #taken from dynesty/dynesty/results.py  but without the nlive point causing an error
    res = ("niter: {:d}\n"
            "ncall: {:d}\n"
            "eff(%): {:6.3f}\n"
            "logz: {:6.3f} +/- {:6.3f}"
            .format(results.niter, sum(results.ncall),
                    results.eff, results.logz[-1], results.logzerr[-1]))

    anc.print_both("", output=of_run)
    anc.print_both('Summary\n=======\n'+res, output=of_run)
    anc.print_both("", output=of_run)
    sys.stdout.flush()

    # with open('filename.pickle', 'wb') as handle:
    #     pickle.dump(a, handle, protocol=pickle.HIGHEST_PROTOCOL)

    # with open('filename.pickle', 'rb') as handle:
    #     b = pickle.load(handle)

    with open(dynesty_output_file, "wb") as handle:
        pickle.dump(results, handle, protocol=pickle.HIGHEST_PROTOCOL)

anc.print_both("Plotting ... ", output=of_run)
plot_folder = os.path.join(
    working_folder, "plots"
)
os.makedirs(plot_folder, exist_ok=True)

# analytic evidence solution
anc.print_both("... runplot: analytic evidence solution ... ", output=of_run)
fig, axes = dyplot.runplot(
    results, 
    color='C0',
    logplot=True, 
)
for ax in axes:
    ax.xaxis.label.set_fontsize(18)
    ax.yaxis.label.set_fontsize(18)
    ax.tick_params(axis='both', labelsize=14)
fig.tight_layout()
fig.align_ylabels(axes)
fig_file = os.path.join(plot_folder, "01_dynesty_runplot.png")
fig.savefig(fig_file, dpi=300, bbox_inches="tight")
plt.close(fig)

# traceplot
anc.print_both("... traceplot ... ", output=of_run)
fig, axes = dyplot.traceplot(
    results, 
    show_titles=True, 
    trace_cmap='viridis',
    title_kwargs={'fontsize': 8, 'y': 0.5, 'x': 1.05, "loc":"left"},
    labels=sim.fitting_names,
    quantiles=None,
    fig=plt.subplots(sim.nfit, 2, figsize=(14, 20))
)
# fig.tight_layout()
fig.align_ylabels(axes)
fig_file = os.path.join(plot_folder, "02_dynesty_traceplot.png")
fig.savefig(fig_file, dpi=300, bbox_inches="tight")
plt.close(fig)

# cornerplot
anc.print_both("... cornerplot ... ", output=of_run)
fig, axes = dyplot.cornerplot(
    results, 
    color="C0",
    show_titles=True, 
    title_kwargs={'y': 1.05, "x": 0.05, "loc":"left"},
    labels=sim.fitting_names,
    quantiles=None,
    fig=plt.subplots(sim.nfit, sim.nfit, figsize=(14, 14))
)

fig.align_ylabels(axes)
for axs in axes:
    for ax in axs:
        ax.xaxis.labelpad = 50
        ax.yaxis.labelpad = 50

# fig.tight_layout()
fig_file = os.path.join(plot_folder, "03_dynesty_cornerplot.png")
fig.savefig(fig_file, dpi=300, bbox_inches="tight")
plt.close(fig)

of_run.close()
sim.reset()

# return

# ==============================================================================
# ==============================================================================

# if __name__ == "__main__":
#   main()
