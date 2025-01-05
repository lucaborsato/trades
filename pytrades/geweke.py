#!/usr/bin/env python
# -*- coding: utf-8 -*-

# no more "zero" integer division bugs!:P
# import argparse
import os
import sys
import numpy as np  # array
# import time

# import h5py
# import random
# import constants as cst # local constants module
# from scipy.stats import norm as scipy_norm
from . import ancillary as anc

import matplotlib.pyplot as plt

anc.set_rcParams()

def compute_geweke(cli, log_folder, plot_folder, chains, fitting_names):

    os.makedirs(log_folder, exist_ok=True)
    log_file = os.path.join(log_folder, "log_Geweke.txt")
    olog = open(log_file, 'w')

    anc.print_both("", output=olog)
    anc.print_both(" ======================== ", output=olog)
    anc.print_both(" GEWEKE PLOTS", output=olog)
    anc.print_both(" ======================== ", output=olog)
    anc.print_both("", output=olog)

    n_steps, nwalkers, nfit = np.shape(chains)
    gk_steps = int(cli.gk_steps)
    if gk_steps == 0:
        gk_steps = 10
    
    cols = 1 + int(np.rint(nwalkers / 40.0))

    anc.print_both("Plotting Geweke ...", output=olog)
    sys.stdout.flush()

    # for ifit in range(0, nfit):
    #     pnames = str(fitting_names[ifit])
    for ifit, pnames in enumerate(fitting_names):
        anc.print_both("Parameter: {:13s}".format(pnames), output=olog)
        fig = plt.figure()

        lower_interval, z_score = anc.geweke_test(
            chains[:, :, ifit], start_frac=0.01, n_sel_steps=gk_steps
        )

        ax = plt.subplot2grid((1, 1), (0, 0))
        for i_c in range(0, nwalkers):
            ax.plot(
                lower_interval,
                z_score[:, i_c],
                ".-",
                label="walker {:d}".format(i_c + 1),
                alpha=0.8,
            )

        ax.axhline(2.0, color="lightgray")
        ax.axhline(-2.0, color="lightgray")
        ax.set_xlabel("steps ({:s})".format(pnames))
        ax.tick_params(axis='x', labelrotation = 45)
        plot_file = os.path.join(
            plot_folder, "Geweke_{:03d}_{:s}.png".format(ifit + 1, pnames)
        )
        fig.savefig(
            plot_file,
            bbox_inches="tight",
            dpi=200,
        )
        plt.close(fig)
        anc.print_both(
            "Saved plot {:s}".format(plot_file), output=olog
        )

    olog.close()


    return


def compute_geweke_from_file(cli):

    # plot_folder = prepare_plot_folder(working_path)
    emcee_folder= cli.full_path
    log_folder = os.path.join(emcee_folder, "logs")

    emcee_file, _, _ = anc.get_emcee_file_and_best(emcee_folder, cli.temp_status)
    emcee_plots = anc.prepare_emcee_plot_folder(cli.full_path)

    # computes mass conversion factor
    m_factor, m_unit, r_factor, r_unit = anc.mass_radius_type_factor(mtype=cli.m_type)

    # get data from the hdf5 file
    # fitting_names, parameter_boundaries, chains, acceptance_fraction, autocor_time, lnprobability, ln_err_const, completed_steps = anc.get_data(emcee_file, cli.temp_status)
    (
        fitting_names,
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
        fitting_names,
        parameter_boundaries,
        chains,
    )

    compute_geweke(cli, log_folder, emcee_plots, chains_T, fitting_names)

    return

def main():

    # read cli arguments
    cli = anc.get_args()

    compute_geweke_from_file(cli)

    return



if __name__ == "__main__":
    main()
