#!/usr/bin/env python
# -*- coding: utf-8 -*-

# no more "zero" integer division bugs!:P
# import argparse
import os
import sys
import numpy as np  # array
import time

# import h5py
# import random
# import constants as cst # local constants module
# from scipy.stats import norm as scipy_norm
import ancillary as anc

import matplotlib as mpl

import matplotlib.pyplot as plt

anc.set_rcParams()


def compute_gr(cli, log_folder, plot_folder, chains, fitting_names):

    os.makedirs(log_folder, exist_ok=True)
    log_file = os.path.join(log_folder, "log_GelmanRubin.txt")
    olog = open(log_file, 'w')

    anc.print_both("", output=olog)
    anc.print_both(" ======================== ", output=olog)
    anc.print_both(" GELMAN-RUBIN PLOTS", output=olog)
    anc.print_both(" ======================== ", output=olog)
    anc.print_both("", output=olog)

    # if cli.temp_status:
    #     n_steps = completed_steps
    # else:
    #     n_steps = nruns
    n_steps, _, nfit = np.shape(chains)
    
    gr_steps = int(cli.gr_steps)
    if gr_steps == 0:
        gr_steps = 10

    # steps = np.linspace(start=0, stop=n_steps, num=gr_steps, endpoint=True, dtype=int)
    # steps[0] = 10
    # gr_steps = len(steps)

    # step = int(n_steps/gr_steps)
    # steps = np.arange(step, n_steps+step, step)
    step = int((n_steps -10)/gr_steps)
    steps = np.rint(np.linspace(step, n_steps, endpoint=True, num=gr_steps)).astype(int)

    anc.print_both('number of GR steps: {}'.format(gr_steps), output=olog)

    gr_Rc_2 = np.ones((gr_steps, nfit)) + 99.0

    anc.print_both("Plotting GR ...", output=olog)
    sys.stdout.flush()

    # for ifit in range(0, nfit):
    #     pnames = str(fitting_names[ifit])
    for ifit, pnames in enumerate(fitting_names):
        anc.print_both("Parameter: {:13s}".format(pnames), output=olog)
        fig = plt.figure()
        ax = plt.subplot2grid((1, 1), (0, 0))

        # for istep in range(0, gr_steps):
        #     gr_Rc_2[istep, ifit] = anc.GelmanRubin(chains[: steps[istep], :, ifit])
        for istep, vstep in enumerate(steps):
            gr_Rc_2[istep, ifit] = anc.GelmanRubin(chains[: vstep, :, ifit])

        ax.axhline(1.01, color="gray")
        ax.plot(steps, gr_Rc_2[:, ifit], "-", color="k", lw=1.3, label="LBo 2")
        ax.set_ylim(0.95, 2.3)
        ax.set_xlabel("steps ({:s})".format(pnames.strip()))
        ax.tick_params(axis='x', labelrotation = 45)
        ax.legend(loc="center left", fontsize=9, bbox_to_anchor=(1, 0.5))
        plot_file = os.path.join(
                plot_folder,
                "Gelman-Rubin_{:03d}_{:s}.png".format(ifit + 1, pnames),
            )
        fig.savefig(
            plot_file,
            bbox_inches="tight",
        )
        plt.close(fig)
        anc.print_both(
            "Saved plot {:s}".format(
                plot_file, output=olog
            )
        )

    olog.close()

    return

def compute_gr_from_file(cli):

    emcee_folder= cli.full_path
    log_folder = os.path.join(emcee_folder, "logs")

    emcee_plots = anc.prepare_emcee_plot_folder(cli.full_path)

    # computes mass conversion factor
    m_factor, m_unit, r_factor, r_unit = anc.mass_radius_type_factor(mtype=cli.m_type)

    emcee_file, _, _ = anc.get_emcee_file_and_best(emcee_folder, cli.temp_status)

    # get data from the hdf5 file
    # parameter_names_emcee, parameter_boundaries, chains, acceptance_fraction, autocor_time, lnprobability, ln_err_const, completed_steps = anc.get_data(emcee_file, cli.temp_status)
    (
        parameter_names_emcee,
        parameter_boundaries,
        chains,
        _,
        _,
        _,
        _,
        completed_steps,
    ) = anc.get_data(emcee_file, cli.temp_status)

    nfit, nwalkers, nruns, nburnin, nruns_sel = anc.get_emcee_parameters(
        chains, cli.temp_status, cli.nburnin, completed_steps
    )

    chains_T, parameter_boundaries = anc.select_transpose_convert_chains(
        nfit,
        nwalkers,
        nburnin,
        nruns,
        nruns_sel,
        m_factor,
        parameter_names_emcee,
        parameter_boundaries,
        chains,
    )

    compute_gr(cli, log_folder, emcee_plots, chains_T, parameter_names_emcee)

    return

def main():

    # read cli arguments
    cli = anc.get_args()

    compute_gr_from_file(cli)

    return


if __name__ == "__main__":
    main()
