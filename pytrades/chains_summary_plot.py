#!/usr/bin/env python
# -*- coding: utf-8 -*-

# no more "zero" integer division bugs!:P
# import argparse
import os
import sys
import numpy as np  # array
import h5py

# import random
import constants as cst  # local constants module
from scipy.stats import norm as scipy_norm
import ancillary as anc

import matplotlib as mpl

import matplotlib.pyplot as plt

anc.set_rcParams()


def plot_chains(
    cli,
    log_folder,
    plot_folder,
    fitting_posterior_in,
    chains_in,
    lnprob_posterior,
    lnprobability,
    fitting_names,
    thin_steps,
    fitting_minmax,
    show_plots = False
):

    os.makedirs(log_folder, exist_ok=True)
    log_file = os.path.join(log_folder, "log_chains.txt")
    olog = open(log_file, "w")

    anc.print_both("", output=olog)
    anc.print_both(" ======================== ", output=olog)
    anc.print_both(" CHAIN PLOTS", output=olog)
    anc.print_both(" ======================== ", output=olog)
    anc.print_both("", output=olog)

    fitting_posterior = fitting_posterior_in.copy()
    chains = chains_in.copy()

    # set label and legend names
    kel_labels = anc.keplerian_legend(fitting_names, cli.m_type)
    k = anc.get_auto_bins(fitting_posterior)
    m_factor, m_unit, r_factor, r_unit = anc.mass_radius_type_factor(mtype=cli.m_type)

    sys.stdout.flush()
    ## OPEN summary_parameters.hdf5 FILE
    anc.print_both("Reading summary_parameters.hdf5 file ...", output=olog)
    with h5py.File(
        os.path.join(cli.full_path, "summary_parameters.hdf5"), "r"
    ) as s_h5f:
        # sample_parameters
        ci_fitted = s_h5f["confidence_intervals/fitted/ci"][...]

        anc.print_both("Set overplot ...", output=olog)
        overplot = anc.set_overplot(cli.overplot)
        if overplot is not None:
            sim_id_str = "{}".format(overplot)
            overp_par = s_h5f["parameters/{:s}/fitted/parameters".format(sim_id_str)][
                ...
            ]
    sys.stdout.flush()

    # nruns_full, _, _ = np.shape(chains_full)
    # nruns_thinned, _, _ = np.shape(chains_full_thinned)
    # if cli.use_thin or cli.use_thin > 0:
    #     nburnin_plt = np.rint(cli.nburnin / thin_steps).astype(int)
    #     nend = np.rint(nruns_thinned / thin_steps).astype(int)
    # else:
    #     nburnin_plt = cli.nburnin
    #     nend = nruns_full
    nruns, _, _ = np.shape(chains) # it could be full or full_thinned, with full steps or completed_steps
    if cli.use_thin > 1:
        nburnin_plt = np.rint(cli.nburnin / thin_steps).astype(int)
        # nend = np.rint(nruns / thin_steps).astype(int)
    else:
        nburnin_plt = cli.nburnin
        # nend = nruns
    nend = nruns # it is always right? It depends on the chains_in, that could be full or thinned
    anc.print_both("thinning? {}".format(cli.use_thin), output=olog)
    anc.print_both("original nburnin = {}".format(cli.nburnin), output=olog)
    anc.print_both("current  nburnin = {}".format(nburnin_plt), output=olog)
    anc.print_both("original   nruns = {}".format(nruns), output=olog)
    anc.print_both("current    nruns = {}".format(nend), output=olog)
    sys.stdout.flush()

    fitting_boundaries = fitting_minmax.copy()

    anc.print_both("Plotting chains ...", output=olog)
    sys.stdout.flush()

    for i, pnames in enumerate(fitting_names):
        if "Ms" in pnames:
            conv_plot = m_factor
        else:
            conv_plot = 1.0

        if (pnames[0] in ["w"]) or (pnames[0:2] == "mA") or ("lambda" in pnames):
            overp_par[i] %= 360.0
            ci_fitted[i, :] %= 360.0
            fitting_boundaries[i, :] %= 360.0
            fitting_posterior[:, i] %= 360.0
            chains[:, :, i] %= 360.0
            # if("lambda" in pnames):
            #     print("{} min = {} max = {}".format(
            #         pnames, 
            #         np.min(fitting_posterior[:, i]),
            #         np.max(fitting_posterior[:, i]),
            #     ))

        chain_file = os.path.join(
            plot_folder,
            "Chain_{:03d}_{:s}.png".format(i + 1, pnames.strip()),
        )
        anc.print_both("{:s} ==> {:s}".format(pnames, chain_file), output=olog)
        sys.stdout.flush()

        fig, (axChain, axHist) = plt.subplots(nrows=1, ncols=2, figsize=(6, 6))

        # ===================
        # POSTERIOR HISTOGRAM
        (counts, bins_val, patches) = axHist.hist(
            fitting_posterior[:, i] * conv_plot,
            bins=k,
            range=(
                fitting_posterior[:, i].min() * conv_plot,
                fitting_posterior[:, i].max() * conv_plot,
            ),
            orientation="horizontal",
            density=True,
            stacked=True,
            histtype="stepfilled",
            color="darkgrey",
            edgecolor="lightgray",
            align="mid",
        )

        xpdf = scipy_norm.pdf(
            fitting_posterior[:, i] * conv_plot,
            loc=fitting_posterior[:, i].mean() * conv_plot,
            scale=fitting_posterior[:, i].std() * conv_plot,
        )
        idx = np.argsort(fitting_posterior[:, i] * conv_plot)
        axHist.plot(
            xpdf[idx],
            fitting_posterior[idx, i] * conv_plot,
            color="black",
            marker="None",
            ls="-.",
            lw=1.5,
            label="pdf",
        )

        # ===================
        # CHAINS PLOT
        # # chains with the burn-in
        # if cli.use_thin > 1:
        #     axChain.plot(chains_full_thinned[:, :, i] * conv_plot, "-", alpha=0.3)
        # else:
        #     axChain.plot(chains_full[:, :, i] * conv_plot, "-", alpha=0.3)
        axChain.plot(chains[:, :, i] * conv_plot, "-", alpha=0.3)

        axChain.axvspan(0, nburnin_plt, color="gray", alpha=0.45)
        axChain.axvline(nburnin_plt, color="gray", ls="-", lw=1.5)

        if overplot is not None:
            axChain.axhline(
                overp_par[i] * conv_plot,
                marker="None",
                c="black",
                ls="--",
                lw=2.5,
                alpha=0.6,
                label="overplot {}".format(overplot.split("sim_")[1]),
            )

        # plot ci
        axChain.axhline(
            ci_fitted[i, 0] * conv_plot,
            marker="None",
            c="forestgreen",
            ls="-",
            lw=2.1,
            alpha=1.0,
            label="CI/HDI 16th (%.5f)" % (ci_fitted[i, 0] * conv_plot),
        )
        axChain.axhline(
            ci_fitted[i, 1] * conv_plot,
            marker="None",
            c="forestgreen",
            ls="-",
            lw=2.1,
            alpha=1.0,
            label="CI/HDI 84th (%.5f)" % (ci_fitted[i, 1] * conv_plot),
        )

        axChain.ticklabel_format(useOffset=False)
        xlabel = "$N_\mathrm{steps}$"
        if cli.use_thin > 1:
            xlabel = "$N_\mathrm{steps} \\times %d$" % (thin_steps)
        axChain.set_xlabel(xlabel)
        axChain.set_xlim([0, nend])
        axChain.set_ylabel(kel_labels[i])

        y_min = fitting_posterior[:, i].min() * conv_plot
        y_max = fitting_posterior[:, i].max() * conv_plot

        axChain.set_ylim([y_min, y_max])
        stitle = "Full chain {:s}:=[{:.3f} , {:.3f}]".format(
            kel_labels[i],
            fitting_boundaries[i, 0] * conv_plot,
            fitting_boundaries[i, 1] * conv_plot,
        )
        axChain.set_title(
            stitle,
            fontsize=plt.rcParams["xtick.labelsize"] - 2,
        )
        axChain.tick_params(axis='x', labelrotation = 45)

        axHist.ticklabel_format(useOffset=False)
        axHist.tick_params(direction="inout", labelleft=False)
        axHist.set_ylim([y_min, y_max])

        if overplot is not None:
            axHist.axhline(
                overp_par[i] * conv_plot,
                marker="None",
                c="black",
                ls="--",
                lw=2.5,
                alpha=0.8,
                label="overplot {}".format(overplot.split("sim_")[1]),
            )

        # plot ci
        axHist.axhline(
            ci_fitted[i, 0] * conv_plot,
            marker="None",
            c="forestgreen",
            ls="-",
            lw=2.1,
            alpha=1.0,
            label="CI/HDI 16th (%.5f)" % (ci_fitted[i, 0] * conv_plot),
        )
        axHist.axhline(
            ci_fitted[i, 1] * conv_plot,
            marker="None",
            c="forestgreen",
            ls="-",
            lw=2.1,
            alpha=1.0,
            label="CI/HDI 84th (%.5f)" % (ci_fitted[i, 1] * conv_plot),
        )

        axHist.set_title(
            "Distribution of posterior chain",
            fontsize=plt.rcParams["xtick.labelsize"] - 2,
        )
        axHist.legend(loc="center left", fontsize=9, bbox_to_anchor=(1, 0.5))

        if show_plots:
            plt.show(fig)
        fig.savefig(chain_file, bbox_inches="tight")
        plt.close(fig)

    fig = plt.figure(figsize=(6, 6))

    # lnprob

    lnprob_file = os.path.join(plot_folder, "emcee_lnprobability.png")
    anc.print_both("lnprob ==> {:s}".format(lnprob_file), output=olog)
    xlabel = "$N_\mathrm{steps}$"
    xlabel_post = "{:s} + {:d}".format(xlabel, nburnin_plt)
    if cli.use_thin > 1:
        xlabel = "$N_\mathrm{{steps}} \\times {:d}$".format(thin_steps)
        xlabel_post = "$N_\mathrm{{steps}} + {:d} \\times {:d}$".format(nburnin_plt, thin_steps)

    nrows, ncols = 3, 1

    axs = []

    ax = plt.subplot2grid((nrows, ncols), (0, 0), rowspan=1)
    ax.plot(lnprobability, "-", alpha=0.3)
    ax.axvline(nburnin_plt, color="gray", ls="-", lw=1.5)
    ax.set_ylabel("lnprob")
    ax.set_xlabel(xlabel)
    axs.append(ax)

    ax = plt.subplot2grid((nrows, ncols), (1, 0), rowspan=2)
    ax.plot(lnprob_posterior, "-", alpha=0.3)

    min_lnp = np.min(lnprob_posterior, axis=0).min()
    max_lnp = np.max(lnprob_posterior, axis=0).max()
    y_min, y_max = anc.compute_limits(np.asarray([min_lnp, max_lnp]), 0.05)
    ax.set_ylim((y_min, y_max))
    ax.set_ylabel("lnprob")
    ax.set_xlabel(xlabel_post)
    axs.append(ax)

    fig.align_ylabels(axs)
    plt.tight_layout()

    if show_plots:
        plt.show(fig)
    fig.savefig(
        lnprob_file,
        bbox_inches="tight",
        # dpi=150,
    )
    plt.close(fig)

    olog.close()

    return



## WARNING ##
## TO BE UPDATED ! ##
def plot_chains_from_file(cli):

    # set emcee folder
    emcee_folder = cli.full_path
    log_folder = os.path.join(emcee_folder, "logs")
    emcee_plots = os.path.join(cli.full_path, "plots")
    os.makedirs(emcee_plots, exist_ok=True)

    # computes mass conversion factor
    m_factor, m_unit, r_factor, r_unit = anc.mass_radius_type_factor(mtype=cli.m_type)
    # and best folder
    emcee_file, emcee_best, folder_best = anc.get_emcee_file_and_best(
        emcee_folder, cli.temp_status
    )

    (
        fitting_names,
        fitting_minmax,
        chains,
        acceptance_fraction,
        autocor_time,
        lnprobability,
        ln_err_const,
        completed_steps,
    ) = anc.get_data(emcee_file, cli.temp_status)

    nfit, nwalkers, nruns, nburnin, nruns_sel = anc.get_emcee_parameters(
        chains, cli.temp_status, cli.nburnin, completed_steps
    )
    chains_T_full, fitting_minmax = anc.select_transpose_convert_chains(
        nfit,
        nwalkers,
        nburnin,
        nruns,
        nruns_sel,
        m_factor,
        fitting_names,
        fitting_minmax,
        chains,
    )

    if cli.use_thin or cli.use_thin > 0:
        (
            chains_T,
            fitting_posterior,
            lnprob_posterior,
            thin_steps,
            chains_T_full_thinned,
        ) = anc.thin_the_chains(
            cli.use_thin,
            nburnin,
            nruns,
            nruns_sel,
            autocor_time,
            chains_T_full,
            lnprobability,
            burnin_done=False,
            full_chains_thinned=True,
        )
        nburnin_plt = np.rint(nburnin / thin_steps).astype(int)
        nend = np.rint(nruns / thin_steps).astype(int)

    else:
        (chains_T, fitting_posterior, lnprob_posterior, thin_steps,) = anc.thin_the_chains(
            cli.use_thin,
            nburnin,
            nruns,
            nruns_sel,
            autocor_time,
            chains_T_full,
            lnprobability,
            burnin_done=False,
            full_chains_thinned=False,
        )
        nburnin_plt = nburnin
        nend = nruns

    # lambda fix
    for name in fitting_names:
        if "lambda" in name:
            i_fit = list(fitting_names).index(name)
            chains_T_full[:, :, i_fit] %= 360.0
            chains_T[:, :, i_fit] %= 360.0
            fitting_posterior[:, i_fit] %= 360.0
            if cli.use_thin or cli.use_thin > 0:
                chains_T_full_thinned[:, :, i_fit] %= 360.0

    plot_chains(
        cli,
        log_folder,
        emcee_plots,
        fitting_posterior,
        chains_T_full,
        chains_T_full_thinned,
        lnprob_posterior,
        fitting_names,
        thin_steps,
        fitting_minmax,
    )


    return


def main():
    # read cli arguments
    cli = anc.get_args()
    plot_chains_from_file(cli)

    return


if __name__ == "__main__":
    main()
