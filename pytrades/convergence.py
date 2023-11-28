#!/usr/bin/env python
# -*- coding: utf-8 -*-

# no more "zero" integer division bugs!:P
# import argparse
import os
import sys
import numpy as np  # array
import time

import emcee

# import h5py
# import random
# import constants as cst # local constants module
# from scipy.stats import norm as scipy_norm
import ancillary as anc

import matplotlib as mpl

import matplotlib.pyplot as plt

anc.set_rcParams()


def compute_convergence(
    chains, fit_names, log_folder, plots_folder, n_cv=10, n_thin=1, figsize=(5, 5)
):
    os.makedirs(log_folder, exist_ok=True)
    log_file = os.path.join(log_folder, "log_convergence.txt")
    with open(log_file, "w") as olog:
        anc.print_both("", output=olog)
        anc.print_both(" ======================== ", output=olog)
        anc.print_both(" CONVERGENCE PLOTS", output=olog)
        anc.print_both(" ======================== ", output=olog)
        anc.print_both("", output=olog)

        n_steps, n_walkers, n_fit = np.shape(chains)

        if n_cv < 1:
            n_cv = 10

        step = int((n_steps - 10) / n_cv)
        steps = np.rint(np.linspace(step, n_steps, endpoint=True, num=n_cv)).astype(int)

        expected_acf = np.zeros((n_fit)).astype(int)
        expected_steps = np.zeros((n_fit)).astype(int)

        for ifit, iname in enumerate(fit_names):
            anc.print_both("\nParameter {}".format(iname), output=olog)

            fig, axs = plt.subplots(nrows=3, ncols=1, sharex=False, figsize=figsize)
            fig.suptitle(iname)

            ax = axs[0]
            # anc.print_both("Gelman-Rubin", output=olog)
            gr = np.zeros((n_cv)) + 100
            for istep in range(n_cv):
                gr[istep] = anc.GelmanRubin(chains[: steps[istep], :, ifit])
            ax.plot(
                steps,
                gr,
                color="black",
                marker="o",
                ms=4,
                mec="white",
                mew=0.3,
                ls="-",
                lw=0.7,
                zorder=5,
            )
            ax.axhline(1.01, color="gray", ls="--", lw=0.7, zorder=4)
            ax.set_ylim(0.95, 1.2)
            ax.set_ylabel("G-R ($\^{R}$)")
            ax.set_xlabel("steps $\\times {}$".format(n_thin))

            ax = axs[1]
            # anc.print_both("Geweke", output=olog)
            lower_interval, z_score = anc.geweke_test(
                chains[:, :, ifit], start_frac=0.01, n_sel_steps=n_cv
            )
            for i_c in range(0, n_walkers):
                ax.plot(
                    lower_interval,
                    z_score[:, i_c],
                    marker="o",
                    ms=2,
                    mec="None",
                    ls="-",
                    lw=0.4,
                    label="walker {:d}".format(i_c + 1),
                    alpha=0.6,
                )
            # ax.legend(loc='best', fontsize=3)
            ax.axhline(
                +2.0,
                color="lightgray",
                ls="-",
                lw=0.7,
            )
            ax.axhline(
                -2.0,
                color="lightgray",
                ls="-",
                lw=0.7,
            )
            ax.set_ylabel("Geweke")
            ax.set_xlabel("steps $\\times {}$".format(n_thin))
            ax.set_ylim(-3, +3)

            ax = axs[2]
            # anc.print_both("ACF", output=olog)
            tolerance = 50
            integrated_ACF = emcee.autocorr.integrated_time(
                chains[:, :, ifit], tol=tolerance, quiet=True
            )
            acf_len = int(np.nanmax(integrated_ACF))
            n_expected = acf_len * tolerance
            anc.print_both(
                "ACF {}x{} expected chain long as n = {}x{} (current {} steps)".format(
                    acf_len, n_thin, n_expected, n_thin, n_steps
                ),
                output=olog,
            )
            expected_steps[ifit] = n_expected
            expected_acf[ifit] = acf_len

            n_acf = acf_len * tolerance

            n_acf = 10
            acf_steps = np.rint(
                np.linspace(acf_len // 2, n_steps, endpoint=True, num=n_acf)
            ).astype(int)
            tau_est = np.zeros((n_acf))
            for i_acf, n_s in enumerate(acf_steps):
                acf_mean = np.zeros((n_s))
                for iw in range(0, n_walkers):
                    acf = emcee.autocorr.function_1d(chains[:n_s, iw, ifit])
                    acf_mean += acf
                acf_mean /= n_walkers

                c = 5
                taus = 2.0 * np.cumsum(acf_mean) - 1.0
                window = emcee.autocorr.auto_window(taus, c)
                tau_est[i_acf] = taus[window]

            ax.plot(
                acf_steps,
                tau_est,
                color="C0",
                marker="o",
                ms=2,
                mec="None",
                ls="-",
                lw=0.5,
                label="$\\tau$",
                zorder=6,
            )
            ax.axhline(
                acf_len,
                color="black",
                ls="-",
                lw=0.7,
                label="ACF = ${}\\times{}$".format(acf_len, n_thin),
                zorder=5,
            )
            ax.plot(
                acf_steps,
                acf_steps / tolerance,
                color="gray",
                marker="None",
                ls="--",
                lw=0.5,
                label="$\\tau = N/({}\\times{})$".format(tolerance, n_thin),
                zorder=4,
            )
            ax.legend(loc="best", fontsize=4)
            ax.set_ylabel("ACF")
            ax.set_xlabel("steps $\\times {}$".format(n_thin))

            plt.tight_layout()
            fig.align_ylabels(axs)
            out_file = os.path.join(
                plots_folder, "{:03d}_{}_convergence.png".format(ifit, iname)
            )
            fig.savefig(out_file, dpi=300, bbox_inches="tight")

            plt.close(fig)

        anc.print_both("", output=olog)
        anc.print_both(
            "All expected steps   for each parameter needed to reach full convergence:\n{}".format(
                expected_steps
            ),
            output=olog,
        )
        anc.print_both(
            "All expected ACF len for each parameter needed to reach full convergence:\n{}".format(
                expected_acf
            ),
            output=olog,
        )
        imax_acf = np.argmax(expected_acf)
        anc.print_both(
            "MAX ACF = {} ==> needed chains of {} steps\n".format(
                expected_acf[imax_acf], expected_steps[imax_acf]
            ),
            output=olog,
        )
    # close olog

    return


def full_statistics(
    chains,
    # post,
    flat_post,
    names,
    pars,
    lnp_post,
    output_folder,
    olog=None,
    ilast=0,
    n_burn=0,
    n_thin=1,
    show_plot=False,
    figsize=(8, 8),
):
    # 68.27% (15.87th-84.13th) ==> alpha = 1. - 0.6827 = 0.3173
    # 95.44% ( 2.28th-97.72th) ==> alpha = 1. - 0.9544 = 0.0456
    # 99.74% ( 0.13th-99.87th) ==> alpha = 1. - 0.9974 = 0.0026
    cred1 = 0.6827
    scred1 = "{:.2f}".format(100 * cred1)
    cred2 = 0.9544
    scred2 = "{:.2f}".format(100 * cred2)
    cred3 = 0.9974
    scred3 = "{:.2f}".format(100 * cred3)

    lsize = 10
    tsize = lsize - 3

    n_steps, n_walkers, n_par = np.shape(chains)

    n_gr = 10
    if n_steps <= n_gr:
        n_gr = n_steps
        step = 1
    else:
        step = int((n_steps - 10) / n_gr)
    steps = np.rint(np.linspace(step, n_steps, endpoint=True, num=n_gr)).astype(int)

    expected_acf = np.zeros((n_par)).astype(int)
    expected_steps = np.zeros((n_par)).astype(int)

    for ipar, pname in enumerate(names):
        p = pars[ipar]
        if (
            pname[0] == "w"
            or pname[0:2] == "mA"
            or pname[0:2] == "lN"
            or "lambda" in pname
        ):
            pmod = p%360.0
            patg = anc.get_arctan_angle(pmod)
            pmin, pmax = flat_post[:, ipar].min(), flat_post[:, ipar].max()
            if np.logical_and(patg >= pmin, patg <= pmax):
                p = patg
            else:
                p = pmod
        
        hdi1 = anc.hpd(flat_post[:, ipar], cred=cred1)
        hdi2 = anc.hpd(flat_post[:, ipar], cred=cred2)
        hdi3 = anc.hpd(flat_post[:, ipar], cred=cred3)
        err = np.array(hdi1) - p
        med = np.median(flat_post[:, ipar])
        warn = ""
        if err[0] > 0 or err[1] < 0:
            warn = "!!WARNING MAP OUT OF HDI{}%!!".format(cred1)

        fig = plt.figure(figsize=figsize, layout="constrained")
        spec = fig.add_gridspec(3, 3)
        axs = []

        # TOP-LEFT ==== Gelman-Rubin
        topl = fig.add_subplot(spec[0, 0])
        topl.tick_params(axis="x", labelrotation=0, labelsize=tsize)
        topl.tick_params(axis="y", labelrotation=45, labelsize=tsize)
        gr = np.zeros((n_gr)) + 100
        for istep in range(n_gr):
            gr[istep] = anc.GelmanRubin(chains[: steps[istep], :, ipar])
        topl.plot(
            steps,
            gr,
            color="black",
            marker="o",
            ms=4,
            mec="white",
            mew=0.3,
            ls="-",
            lw=0.7,
            zorder=5,
        )
        topl.axhline(1.01, color="gray", ls="--", lw=0.7, zorder=4)
        ylim0 = topl.get_ylim()
        topl.set_ylim(max(0.95, ylim0[0]), min(1.2, ylim0[1]))
        # topl.set_ylabel("G-R ($\^{R}$)", fontsize=lsize)
        topl.set_xlabel("steps $\\times {}$".format(n_thin), fontsize=lsize)

        axs.append(topl)

        # TOP-CENTRE ==== Geweke
        topc = fig.add_subplot(spec[0, 1])
        topc.tick_params(axis="x", labelrotation=0, labelsize=tsize)
        topc.tick_params(axis="y", labelrotation=45, labelsize=tsize)
        lower_interval, z_score = anc.geweke_test(
            chains[:, :, ipar], start_frac=0.01, n_sel_steps=n_gr
        )
        for i_c in range(0, n_walkers):
            topc.plot(
                lower_interval,
                z_score[:, i_c],
                marker="o",
                ms=2,
                mec="None",
                ls="-",
                lw=0.4,
                # label="walker {:d}".format(i_c + 1),
                alpha=0.6,
            )
        topc.axhline(
            +2.0,
            color="lightgray",
            ls="-",
            lw=0.7,
        )
        topc.axhline(
            -2.0,
            color="lightgray",
            ls="-",
            lw=0.7,
        )
        topc.set_ylim(-3, +3)
        # topc.set_ylabel("Geweke", fontsize=lsize)
        topc.set_xlabel("steps $\\times {}$".format(n_thin), fontsize=lsize)

        axs.append(topc)

        # TOP-RIGHT ==== ACF
        topr = fig.add_subplot(spec[0, 2])
        topr.tick_params(axis="x", labelrotation=0, labelsize=tsize)
        topr.tick_params(axis="y", labelrotation=45, labelsize=tsize)
        tolerance = 50
        integrated_ACF = emcee.autocorr.integrated_time(
            chains[:, :, ipar], tol=tolerance, quiet=True
        )
        acf_len = int(np.nanmax(integrated_ACF))
        n_expected = acf_len * tolerance
        anc.print_both(
            "ACF {}x{} expected chain long as n = {}x{} (current {} steps)".format(
                acf_len, n_thin, n_expected, n_thin, n_steps
            ),
            output=olog,
        )
        expected_steps[ipar] = n_expected
        expected_acf[ipar] = acf_len

        n_acf = acf_len * tolerance  # ??

        n_acf = 10
        acf_steps = np.rint(
            np.linspace(acf_len // 2, n_steps, endpoint=True, num=n_acf)
        ).astype(int)
        tau_est = np.zeros((n_acf))
        for i_acf, n_s in enumerate(acf_steps):
            acf_mean = np.zeros((n_s))
            for iw in range(0, n_walkers):
                acf = emcee.autocorr.function_1d(chains[:n_s, iw, ipar])
                acf_mean += acf
            acf_mean /= n_walkers

            c = 5
            taus = 2.0 * np.cumsum(acf_mean) - 1.0
            window = emcee.autocorr.auto_window(taus, c)
            tau_est[i_acf] = taus[window]

        topr.plot(
            acf_steps,
            tau_est,
            color="C0",
            marker="o",
            ms=2,
            mec="None",
            ls="-",
            lw=0.5,
            label="$\\tau$",
            zorder=6,
        )
        topr.axhline(
            acf_len,
            color="black",
            ls="-",
            lw=0.7,
            label="ACF = ${}\\times{}$".format(acf_len, n_thin),
            zorder=5,
        )
        topr.plot(
            acf_steps,
            acf_steps / tolerance,
            color="gray",
            marker="None",
            ls="--",
            lw=0.5,
            label="$\\tau = N/({}\\times{})$".format(tolerance, n_thin),
            zorder=4,
        )
        topr.legend(loc="best", fontsize=tsize - 2)
        # topr.set_ylabel("ACF", fontsize=lsize)
        topr.set_xlabel("steps $\\times {}$".format(n_thin), fontsize=lsize)

        axs.append(topr)

        # MIDLEFT ==== trace full chains
        midl = fig.add_subplot(spec[1, 0])
        midl.tick_params(axis="x", labelrotation=0, labelsize=tsize)
        midl.tick_params(axis="y", labelrotation=45, labelsize=tsize)
        midl.plot(chains[:, :, ipar], ls="-", lw=0.2, alpha=0.3)
        midl.axvline(n_burn, color="gray", ls="-", lw=1.3, alpha=0.7)
        midl.axhline(p, color="C1", ls="-", lw=1.4, alpha=0.7)
        # midl.set_ylabel("{} (full)".format(pname), fontsize=lsize)
        midl.set_xlabel("steps $\\times {}$".format(n_thin), fontsize=lsize)

        axs.append(midl)

        # MIDCENTRE ==== trace posterior chains
        midc = fig.add_subplot(spec[1, 1])
        midc.tick_params(axis="x", labelrotation=0, labelsize=tsize)
        midc.tick_params(axis="y", labelrotation=45, labelsize=tsize)
        midc.plot(chains[:, :, ipar], ls="-", lw=0.2, alpha=0.3)
        midc.axvspan(
            0, n_burn, facecolor="gray", edgecolor="None", ls="-", lw=1.3, alpha=0.5
        )
        midc.axvline(n_burn, color="gray", ls="-", lw=1.3, alpha=0.7)
        midc.axhline(p, color="C1", ls="-", lw=1.4, alpha=0.7)
        y = flat_post[:, ipar]
        dy = np.ptp(y)
        midc.set_ylim([y.min() - 0.03 * dy, y.max() + 0.03 * dy])
        midc.set_xlabel("steps $\\times {}$".format(n_thin), fontsize=lsize)
        axs.append(midc)

        # MIDRIGHT ==== posterior distribution
        midr = fig.add_subplot(spec[1, 2])
        midr.tick_params(axis="x", labelbottom=False)
        midr.tick_params(axis="y", labelleft=False)
        midr.hist(
            flat_post[:, ipar],
            bins=33,
            color="black",
            density=False,
            orientation="horizontal",
            zorder=3,
        )
        midr.axhline(
            p, color="C1", ls="-", lw=1.3, alpha=1.0, label="MAP", zorder=5
        )
        midr.axhline(
            hdi1[0],
            color="C2",
            ls="--",
            lw=0.95,
            alpha=1.0,
            label="HDI{}%".format(scred1),
            zorder=4,
        )
        midr.axhline(hdi1[1], color="C2", ls="--", lw=0.95, alpha=1.0, zorder=4)
        midr.axhline(
            hdi2[0],
            color="C3",
            ls="--",
            lw=0.50,
            alpha=1.0,
            label="HDI{}%".format(scred2),
            zorder=4,
        )
        midr.axhline(hdi2[1], color="C3", ls="--", lw=0.50, alpha=1.0, zorder=5)
        midr.axhline(
            hdi3[0],
            color="C4",
            ls="--",
            lw=0.50,
            alpha=1.0,
            label="HDI{}%".format(scred3),
            zorder=4,
        )
        midr.axhline(hdi3[1], color="C4", ls="--", lw=0.50, alpha=1.0, zorder=6)
        midr.axhline(
            med, color="C0", ls="--", lw=1.0, alpha=1.0, label="MEDIAN", zorder=6
        )
        # midr.legend(loc='best', fontsize=tsize-3)

        axs.append(midr)

        # BOTTOM ==== lnP = f(par)
        bot = fig.add_subplot(spec[2, :])
        bot.tick_params(axis="x", labelrotation=0, labelsize=tsize)
        bot.tick_params(axis="y", labelrotation=45, labelsize=tsize)
        # bot.set_title(pname, fontsize=lsize+1)
        bot.plot(
            flat_post[:, ipar],
            lnp_post,
            color="black",
            marker="o",
            ms=1,
            mec="None",
            ls="",
            alpha=0.33,
            zorder=2,
        )
        bot.axvline(
            p, color="C1", ls="-", lw=1.3, alpha=1.0, label="MAP", zorder=5
        )
        bot.axvline(
            hdi1[0],
            color="C2",
            ls="--",
            lw=0.95,
            alpha=1.0,
            label="HDI{}%".format(scred1),
            zorder=4,
        )
        bot.axvline(hdi1[1], color="C2", ls="--", lw=0.95, alpha=1.0, zorder=4)
        bot.axvline(
            hdi2[0],
            color="C3",
            ls="--",
            lw=0.50,
            alpha=1.0,
            label="HDI{}%".format(scred2),
            zorder=4,
        )
        bot.axvline(hdi2[1], color="C3", ls="--", lw=0.50, alpha=1.0, zorder=5)
        bot.axvline(
            hdi3[0],
            color="C4",
            ls="--",
            lw=0.50,
            alpha=1.0,
            label="HDI{}%".format(scred3),
            zorder=4,
        )
        bot.axvline(hdi3[1], color="C4", ls="--", lw=0.50, alpha=1.0, zorder=6)
        bot.axvline(
            med, color="C0", ls="--", lw=1.0, alpha=1.0, label="MEDIAN", zorder=6
        )
        bot.legend(
            # loc='center left',
            # bbox_to_anchor=(1.01, 0.5),
            loc="best",
            fontsize=tsize - 3,
        )

        bot.set_ylabel("$\ln\mathcal{P}$", fontsize=lsize)
        bot.set_xlabel("{} (posterior)".format(pname), fontsize=lsize)

        axs.append(bot)

        plt.tight_layout()
        # fig.align_ylabels(axs)
        # save figure
        output_file = os.path.join(
            output_folder, "{:03d}_{}.png".format(ipar + ilast, pname)
        )
        anc.print_both("Saving {}".format(output_file), output=olog)
        fig.savefig(output_file, dpi=300, bbox_inches="tight")
        if show_plot:
            plt.show()
        plt.close(fig)

    return expected_acf, expected_steps


def log_probability_trace(
    log_prob, lnprob_posterior, plot_folder, n_burn=0, n_thin=1, show_plot=False, figsize=(8, 8), olog=None
):
    lsize = 10
    tsize = lsize - 3

    map_lnprob = np.max(lnprob_posterior)

    fig = plt.figure(figsize=figsize, layout="constrained")
    spec = fig.add_gridspec(3, 1)
    axs = []

    top = fig.add_subplot(spec[0, 0])
    top.tick_params(axis="x", labelrotation=45, labelsize=tsize, labelbottom=False)
    top.tick_params(axis="y", labelrotation=45, labelsize=tsize)
    top.plot(log_prob, ls="-", lw=0.2, alpha=0.3)
    top.axvline(n_burn, color="gray", ls="-", lw=1.3, alpha=0.7)
    top.axhline(map_lnprob, color="C1", ls="-", lw=1.2, alpha=0.7)
    top.set_ylabel("$\ln\mathcal{P}$ (full)")
    axs.append(top)

    mid = fig.add_subplot(spec[1, 0])
    mid.tick_params(axis="x", labelrotation=45, labelsize=tsize)
    mid.tick_params(axis="y", labelrotation=45, labelsize=tsize)
    mid.plot(log_prob, ls="-", lw=0.2, alpha=0.3)
    mid.axvline(n_burn, color="gray", ls="-", lw=1.3, alpha=0.7)
    mid.axhline(map_lnprob, color="C1", ls="-", lw=1.2, alpha=0.7)
    mid.set_ylabel("$\ln\mathcal{P}$ (post.)")
    mid.set_xlabel("steps $\\times {}$".format(n_thin))

    dlnP = np.ptp(lnprob_posterior)
    mid.set_ylim(
        [lnprob_posterior.min() - 0.03 * dlnP, lnprob_posterior.max() + 0.03 * dlnP]
    )
    axs.append(mid)

    bot = fig.add_subplot(spec[2, 0])
    bot.tick_params(axis="x", labelrotation=45, labelsize=tsize)
    bot.tick_params(axis="y", labelrotation=45, labelsize=tsize, labelleft=False)
    bot.hist(
        lnprob_posterior,
        bins=33,
        color="black",
        density=False,
        orientation="vertical",
        zorder=3,
    )
    bot.axvline(map_lnprob, color="C1", ls="-", lw=1.2, alpha=0.7)
    bot.set_xlabel("$\ln\mathcal{P}$ (post.)")
    axs.append(bot)

    plt.tight_layout()
    fig.align_ylabels(axs)
    output_file = os.path.join(plot_folder, "lnprob_trace.png")
    anc.print_both("\nSaving {}".format(output_file), output=olog)
    fig.savefig(output_file, dpi=300, bbox_inches="tight")
    if show_plot:
        plt.show()
    plt.close(fig)

    return
