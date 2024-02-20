# pytrades module, it imports the fortran-python library and creates a trades object
import numpy as np
import os
import h5py
import sys

from pytrades_lib import f90trades
import ancillary as anc
import constants as cst

from scipy.interpolate import interp1d

import pytransit
import numba

# os.environ["OMP_NUM_THREADS"] = "1"

# numba.set_num_threads(1)
# numba.config.THREADING_LAYER = "tbb"

# numba.config.DISABLE_JIT = 1


from plot_oc import set_unit_base

import matplotlib.pyplot as plt

# =============================================================================
anc.set_rcParams()

# =============================================================================
set_fitness = f90trades.set_fitness
convert_trades_par_to_kepelem = f90trades.convert_trades_par_to_kepelem
compute_ln_priors = f90trades.compute_ln_priors
check_boundaries = f90trades.check_boundaries
set_one_fit_par_boundaries = f90trades.set_one_fit_par_boundaries
reset_all_fit_boundaries = f90trades.reset_all_fit_boundaries
orbits_to_elements = f90trades.orbits_to_elements
# =============================================================================

def args_init(n_body, duration_check, t_epoch=None, t_start=None, t_int=None):

    f90trades.args_init(
        n_body,
        duration_check,
    )

    if t_epoch is not None:
        f90trades.set_time_epoch(t_epoch)
    if t_start is not None:
        f90trades.set_time_start(t_start)
    if t_int is not None:
        f90trades.set_time_int(t_int)

    return

# =============================================================================
# get_priors = f90trades.get_priors
def get_priors():

    n_priors = (f90trades.n_priors).item()
    names, values = f90trades.get_priors(n_priors)
    names = [n.decode("utf-8") for n in names]
    return names, values

# =============================================================================
def deallocate_rv_dataset():

    f90trades.deallocate_rv_dataset()

    return


# =============================================================================
def set_rv_dataset(t_rv, rv_obs, erv_obs, rv_setid=None, n_rvset=1):

    if rv_setid is not None:
        n = len(np.unique(rv_setid))
        if n != n_rvset:
            n_rvset = n
    else:
        rv_setid = np.ones((len(t_rv)))
    f90trades.set_rv_dataset(t_rv, rv_obs, erv_obs, rv_setid.astype(int), n_rvset)

    return


# =============================================================================
def set_t0_dataset(body_id, epo_b, t0_b, et0_b):

    f90trades.set_t0_dataset(body_id, epo_b, t0_b, et0_b)

    return


# =============================================================================
def deallocate_t0_dataset(body_id):

    f90trades.deallocate_t0_dataset(body_id)

    return


# =============================================================================
def kelements_to_rv_and_t0s(
    t_start,
    t_epoch,
    t_int,
    M_msun,
    R_rsun,
    P_day,
    ecc,
    argp_deg,
    mA_deg,
    inc_deg,
    lN_deg,
    transit_flag,
    # n_rv,
    # n_T0s,
):

    n_kep = 8  # number of keplerian elements in output for each T0s
    n_rv = f90trades.nrv
    n_T0s = f90trades.ntts
    (
        rv_sim,
        body_T0_sim,
        epo_sim,
        t0_sim,
        t14_sim,
        kel_sim,
    ) = f90trades.kelements_to_rv_and_t0s(
        t_start,
        t_epoch,
        t_int,
        M_msun,
        R_rsun,
        P_day,
        ecc,
        argp_deg,
        mA_deg,
        inc_deg,
        lN_deg,
        transit_flag,
        n_rv,
        n_T0s,
        n_kep,
    )

    return rv_sim, body_T0_sim, epo_sim, t0_sim, t14_sim, kel_sim


# =============================================================================
def kelements_to_rv(
    t_start,
    t_epoch,
    t_int,
    M_msun,
    R_rsun,
    P_day,
    ecc,
    argp_deg,
    mA_deg,
    inc_deg,
    lN_deg,
    # n_rv,
):

    n_rv = f90trades.nrv
    rv_sim = f90trades.kelements_to_rv(
        t_start,
        t_epoch,
        t_int,
        M_msun,
        R_rsun,
        P_day,
        ecc,
        argp_deg,
        mA_deg,
        inc_deg,
        lN_deg,
        n_rv,
    )

    return rv_sim


# =============================================================================
def kelements_to_t0s(
    t_start,
    t_epoch,
    t_int,
    M_msun,
    R_rsun,
    P_day,
    ecc,
    argp_deg,
    mA_deg,
    inc_deg,
    lN_deg,
    transit_flag,
    # n_T0s,
):

    n_kep = 8  # number of keplerian elements in output for each T0s
    n_T0s = f90trades.ntts
    body_T0_sim, epo_sim, t0_sim, t14_sim, kel_sim = f90trades.kelements_to_t0s(
        t_start,
        t_epoch,
        t_int,
        M_msun,
        R_rsun,
        P_day,
        ecc,
        argp_deg,
        mA_deg,
        inc_deg,
        lN_deg,
        transit_flag,
        n_T0s,
        n_kep,
    )

    return body_T0_sim, epo_sim, t0_sim, t14_sim, kel_sim


# =============================================================================
def kelements_to_orbits(
    steps, M_msun, R_rsun, P_day, ecc, argp_deg, mA_deg, inc_deg, lN_deg
):

    nb_dim = len(M_msun) * 6

    orbits = f90trades.kelements_to_orbits(
        nb_dim,
        steps,
        M_msun,
        R_rsun,
        P_day,
        ecc,
        argp_deg,
        mA_deg,
        inc_deg,
        lN_deg,
    )

    return orbits


# =============================================================================


def kelements_to_orbits_full(
    t_epoch,
    t_start,
    t_int,
    M_msun,
    R_rsun,
    P_day,
    ecc,
    argp_deg,
    mA_deg,
    inc_deg,
    lN_deg,
    specific_times=None,
    step_size=None,
    n_steps_smaller_orbits=10.0,
):

    t_end = t_start + t_int
    dt_start_epoch = t_start - t_epoch
    dt_end_epoch = t_end - t_epoch

    n_body = len(M_msun)
    nb_dim = n_body * 6

    if step_size is None:
        if n_steps_smaller_orbits is None:
            n_steps_smaller_orbits = 10.0
        step_size = np.min(P_day[1:]) / n_steps_smaller_orbits
    # print("step_size = {:.5f}".format(step_size))

    if (specific_times is not None) and (len(specific_times) > 0):
        t_add = specific_times
    else:
        t_add = np.array([])

    if np.abs(dt_start_epoch) > 0.0:
        # backward integration
        steps_backward = np.arange(
            t_epoch, t_start, step_size * np.sign(dt_start_epoch)
        )
        if steps_backward[-1] > t_start:
            steps_backward = np.concatenate((steps_backward, [t_start]))
        # add specific times to compute the orbits, such as observed RV times
        t_in = t_add[t_add <= t_epoch]
        steps_backward = np.sort(np.concatenate((steps_backward, t_in)))[
            ::-1
        ]  # reverse the sort

        # n_steps_backward = len(steps_backward)
        orbits_backward = kelements_to_orbits(
            steps_backward,
            M_msun,
            R_rsun,
            P_day,
            ecc,
            argp_deg,
            mA_deg,
            inc_deg,
            lN_deg,
        )
    else:
        steps_backward = np.array([])
        orbits_backward = np.empty((0, nb_dim))

    if np.abs(dt_end_epoch) > 0.0:
        # forward integration
        steps_forward = np.arange(t_epoch, t_end, step_size * np.sign(dt_end_epoch))
        if steps_forward[-1] < t_end:
            steps_forward = np.concatenate((steps_forward, [t_end]))
        # add specific times to compute the orbits, such as observed RV times
        t_in = t_add[t_add >= t_epoch]
        steps_forward = np.sort(np.concatenate((steps_forward, t_in)))

        # n_steps_forward = len(steps_forward)
        orbits_forward = kelements_to_orbits(
            steps_forward, M_msun, R_rsun, P_day, ecc, argp_deg, mA_deg, inc_deg, lN_deg
        )
    else:
        steps_forward = np.array([])
        orbits_forward = np.empty((0, nb_dim))

    time_steps = np.concatenate((steps_backward, steps_forward))
    orbits = np.concatenate((orbits_backward, orbits_forward))

    return time_steps, orbits


# =============================================================================


def orbits_to_rvs(M_msun, orbits):

    rvs = f90trades.orbits_to_rvs(M_msun, orbits)

    return rvs


# =============================================================================


def orbits_to_transits(
    n_all_transits, time_steps, M_msun, R_rsun, orbits, transiting_body
):

    transits, durations, kep_elem, body_flag = f90trades.orbits_to_transits(
        n_all_transits, time_steps, M_msun, R_rsun, orbits, transiting_body
    )
    # remove elements with body_flag < 2
    sel = body_flag > 1
    transits, durations, kep_elem, body_flag = (
        transits[sel],
        durations[sel],
        kep_elem[sel, :],
        body_flag[sel],
    )

    return transits, durations, kep_elem, body_flag


# =============================================================================


def linear_fit(x, y, ey=None):

    if ey is None:
        m, err_m, q, err_q = f90trades.linear_fit_no_errors(x, y)
        res = y - (q + m * x)
    else:
        m, err_m, q, err_q = f90trades.linear_fit_errors(x, y, ey)
        res = (y - (q + m * x)) / ey
    chi2 = np.dot(res, res)

    return (q, err_q), (m, err_m), chi2


# =============================================================================
# PHOTO-DYNAMICAL TRADES
# FUNCTION AND CLASS TO LOAD/STORE DATA (PHOTOMETRY, TRANSIT TIMES, AND RADIAL VELOCIITES)
# INTEGRATES THE ORBITS FROM MASSES, RADIUS AND ORBITAL PARAMETERS
# AND COMPUTE PHOTOMETRY, TRANSIT TIMES AND RADIAL VELOCITIES IF THEY ARE PROVIDE
# DOES NOT IMPLEMENT LOG-LIKELIHOOD FUNCTION
# =============================================================================
def orbital_parameters_to_transits(
    t_start, t_epoch, t_int, mass, radius, period, ecc, w, ma, inc, long, t_rv_obs
):

    time_steps, orbits = kelements_to_orbits_full(
        t_start,
        t_epoch,
        t_int,
        mass,
        radius,
        period,
        ecc,
        w,
        ma,
        inc,
        long,
        specific_times=t_rv_obs,  # RV TIMES
    )

    sel_t_rv = np.isin(time_steps, t_rv_obs)
    rv = orbits_to_rvs(mass, orbits[sel_t_rv, :])
    rv_sim = {"time": time_steps[sel_t_rv], "rv": rv}

    transiting_body = 1  # all planets
    # determine the number of transits ... it has to be done in advance
    n_transits = (t_int / period[1:]).astype(int)
    n_body = len(mass)
    n_all_transits = np.sum(n_transits) + (
        n_body - 1
    )  # star has no transits by definition
    transits, durations, kep_elem, body_flag = orbits_to_transits(
        n_all_transits, time_steps, mass, radius, orbits, transiting_body
    )
    # kep_elem == period 0, sma 1, ecc 2, inc 3, meana 4, argp 5, truea 6, longn 7
    return time_steps, orbits, transits, durations, kep_elem, body_flag, rv_sim


def set_transit_parameters(radius, transits, body_flag, kep_elem):

    rp_rs = np.zeros(len(transits)) + radius[body_flag - 1] / radius[0]
    per = kep_elem[:, 0]
    aRs = kep_elem[:, 1] / (radius[0] * cst.RsunAU)
    inc = kep_elem[:, 3] * cst.deg2rad
    ecc = kep_elem[:, 2]
    w = kep_elem[:, 5] * cst.deg2rad

    return rp_rs, per, aRs, inc, ecc, w


def get_simulate_flux(
    tm, 
    photometry, 
    transits, 
    durations, 
    rp_rs, ld_quad, per, aRs, inc, ecc, w, 
    time_key="time"
):
    if "full" in time_key:
        n_over_key = "n_oversample_full"
        t_exp_key = "t_exp_d_full"
    else:
        n_over_key = "n_oversample"
        t_exp_key = "t_exp_d"
    sim_photometry = {}
    for k, vis in photometry.items():
        t = vis[time_key]

        tra_in_t = np.logical_and(transits >= t.min(), transits <= t.max())
        n_tra = np.sum(tra_in_t)
        tra_dur_in_t = np.logical_and(
            transits - 0.5 * durations * cst.min2day >= t.min(),
            transits + 0.5 * durations * cst.min2day <= t.max(),
        )
        n_dur = np.sum(tra_dur_in_t)
        n = max(n_tra, n_dur)
        if n > 0:
            # print("n(t) = ", len(t), " ld_quad = ", ld_quad, "n = ", n)
            tm.set_data(t, nsamples=vis[n_over_key], exptimes=vis[t_exp_key])
            flux = tm.evaluate(
                k=rp_rs[tra_in_t],
                # ldc=np.tile(ld_quad, [n,1,1]),
                ldc=np.array([ld_quad] * n),
                t0=transits[tra_in_t],
                p=per[tra_in_t],
                a=aRs[tra_in_t],
                i=inc[tra_in_t],
                e=ecc[tra_in_t],
                w=w[tra_in_t],
            )
            f2d = np.atleast_2d(flux - 1.0)
            flux_ = np.sum(f2d, axis=0) + 1.0
        else:
            flux_ = np.ones((len(t)))
        sim_photometry[k] = flux_

    return sim_photometry


def add_photometry_trend(
    obs_photometry, 
    sim_photometry, 
    photometry_coeff_trend, 
    ancillary_coeff, 
    time_key="time"
):

    out_trend = {}
    out_photometry = {}
    for name_phot, phot in obs_photometry.items():
        sim_phot = sim_photometry[name_phot]
        out_photometry[name_phot] = {}
        out_trend[name_phot] = {}
        for k, vis in phot.items():
            t = vis[time_key]
            nt = len(t)
            flux = sim_phot[k].copy()
            pho_trend = photometry_coeff_trend[name_phot][k]
            tscale = t - t.min()
            ftrend = np.ones((nt))
            if pho_trend is not None:
                ftrend = np.zeros((nt))
                for o, c in enumerate(pho_trend):
                    ftrend += c * (tscale**o)
            anc_phot = vis["ancillary_interp"]
            acoeff = ancillary_coeff[name_phot][k]
            if (anc_phot is not None) and (acoeff is not None):
                if (len(anc_phot) > 0) and (len(acoeff) > 0):
                    fanc = 0.0
                    for i_c, kk in enumerate(anc_phot.keys()):
                        fanc += acoeff[i_c] * anc_phot[kk](t)
                    ftrend += fanc
            flux *= ftrend

            out_photometry[name_phot][k] = flux
            out_trend[name_phot][k] = ftrend

    return out_photometry, out_trend


def plot_photometry(
    photometry,
    sim_photometry,
    mod_photometry=None,
    trend_photometry=None,
    figsize=(3, 3),
    show_plot=True,
    output_folder=None,
    return_rms=False,
):

    rms_photometry = {}

    for name_phot, phot in photometry.items():
        sim_phot = sim_photometry[name_phot]
        if mod_photometry is not None:
            mod_phot = mod_photometry[name_phot]
        if trend_photometry is not None:
            trend_phot = trend_photometry[name_phot]

        rms_vis = {}
        for i_k, vis in phot.items():
            t = vis["time"]
            f = vis["flux"]
            ef = vis["flux_err"]
            flux = sim_phot[i_k]
            res = f - flux

            lsize = plt.rcParams["font.size"] - 2

            fig = plt.figure(figsize=figsize)
            title = "{} - id {}".format(name_phot, i_k)
            # plt.title(title, fontsize=lsize+1)
            axs = []
            nrows, ncols = 3, 1

            irow, icol = 0, 0
            rows = 2
            ax = plt.subplot2grid((nrows, ncols), (irow, icol), rowspan=rows)
            ax.set_title(title, fontsize=lsize + 1)
            ax.ticklabel_format(useOffset=False)
            ax.get_xaxis().set_ticks([])
            ax.tick_params(axis="both", labelsize=plt.rcParams["xtick.labelsize"] - 4)

            ax.errorbar(
                t,
                f,
                yerr=ef,
                color="black",
                marker="o",
                ms=1.8,
                mew=0.4,
                mec="white",
                ls="",
                elinewidth=0.7,
                ecolor="gray",
                capsize=0,
                zorder=4,
                label="obs",
            )
            if mod_photometry is not None:
                m = mod_phot[i_k]
                ax.plot(
                    t, m, color="C1", marker="o", ms=0.6, ls="", zorder=6, label="mod"
                )
            if trend_photometry is not None:
                m = trend_phot[i_k]
                ax.plot(
                    t, m, color="C2", marker="o", ms=0.4, ls="", zorder=5, label="trend"
                )
            ax.plot(
                t,
                flux,
                color="C0",
                marker="o",
                ms=1.0,
                ls="",
                zorder=7,
                label="mod+trend",
            )
            ax.set_ylabel("flux", fontsize=lsize)
            ax.legend(
                bbox_to_anchor=(1.02, 0.5),
                loc="center left",
                fontsize=5,
            )
            axs.append(ax)

            irow, icol = 2, 0
            rows = 1
            ax = plt.subplot2grid((nrows, ncols), (irow, icol), rowspan=rows)
            ax.ticklabel_format(useOffset=False)
            ax.tick_params(axis="both", labelsize=plt.rcParams["xtick.labelsize"] - 4)

            ax.axhline(0.0, color="black", ls="-", lw=0.8, zorder=2)
            ax.errorbar(
                t,
                res * 1.0e6,
                yerr=ef * 1.0e6,
                color="black",
                marker="o",
                ms=1.8,
                mew=0.4,
                mec="white",
                ls="",
                elinewidth=0.7,
                ecolor="gray",
                capsize=0,
                zorder=4,
            )
            ax.set_ylabel("res. (ppm)", fontsize=lsize)
            ax.set_xlabel("time", fontsize=lsize)
            axs.append(ax)

            # plt.tight_layout()
            plt.subplots_adjust(hspace=0.05)
            fig.align_ylabels(axs)
            if output_folder is not None:
                output_file = os.path.join(
                    output_folder, "{}_lc{:02d}.png".format(name_phot, i_k)
                )
                fig.savefig(output_file, dpi=300, bbox_inches="tight")
            if show_plot:
                plt.show()
            plt.close(fig)
            rms_vis[i_k] = np.std(res)
        rms_photometry[name_phot] = rms_vis

    if return_rms:
        return rms_photometry
    else:
        return None


def phototrades_plot_oc(
    transits_obs,
    Tref,
    Pref,
    transits_sim,
    pl_names,
    figsize=(8, 4),
    show_plot=True,
    output_folder=None,
):

    fig = plt.figure(figsize=figsize)

    npl = len(pl_names)

    cmap = plt.get_cmap("nipy_spectral")
    cmm = [0.1, 0.9]
    dmap = np.ptp(cmm) / npl
    colval = cmm[0]

    ms = 3

    nrows, ncols = 3 * npl, 1

    tl_size = plt.rcParams["xtick.labelsize"] - 4
    ll_size = plt.rcParams["xtick.labelsize"] - 2
    le_size = tl_size

    xlims = []
    axs = []

    icol = 0
    irow = 0
    for bd, obs in transits_obs.items():
        Tr, Pr = Tref[bd][0], Pref[bd][0]

        epo_obs = obs["epoch"]
        tra_obs = obs["T0s"]
        etra_obs = obs["err_T0s"]
        tlin = Tr + epo_obs * Pr
        oc_obs = tra_obs - tlin
        Aoc_d = np.ptp(oc_obs)
        u = set_unit_base("auto", Aoc_d)

        tra_sim = transits_sim[bd]
        oc_sim = tra_sim - tlin

        res = oc_obs - oc_sim

        # plot O-C
        ax = plt.subplot2grid((nrows, ncols), (irow, icol), rowspan=2)
        ax.ticklabel_format(useOffset=False)
        # ax.tick_params(axis="both", labelsize=tl_size, labelbottom=False)

        ax.axhline(0.0, color="black", ls="-", lw=0.8, zorder=2)
        ax.errorbar(
            tlin,
            oc_obs * u[0],
            yerr=etra_obs * u[0],
            color=cmap(colval),
            marker="o",
            ms=ms,
            mfc=cmap(colval),
            mec="white",
            mew=0.5,
            ls="",
            ecolor=cmap(colval),
            elinewidth=0.7,
            capsize=0,
            zorder=3,
            label="obs planet {}".format(pl_names[bd]),
        )
        ax.plot(
            tlin,
            oc_sim * u[0],
            color=cmap(colval),
            marker="o",
            ms=ms + 0.5,
            mfc="None",
            mec=cmap(colval),
            mew=0.5,
            ls="",
            zorder=4,
            label="sim",
        )
        ax.legend(loc="best", fontsize=le_size)
        ax.set_ylabel("O-C ({})".format(u[1]), fontsize=ll_size)
        # ax.set_xticks([])
        xlims += ax.get_xlim()
        axs.append(ax)

        irow += 2
        # plot res
        ax = plt.subplot2grid((nrows, ncols), (irow, icol), rowspan=1)
        ax.ticklabel_format(useOffset=False)
        # ax.tick_params(axis="both", labelsize=tl_size, labelbottom=False)

        ax.axhline(0.0, color="black", ls="-", lw=0.8, zorder=2)
        ax.errorbar(
            tlin,
            res * u[0],
            yerr=etra_obs * u[0],
            color=cmap(colval),
            marker="o",
            ms=ms,
            mfc=cmap(colval),
            mec="white",
            mew=0.5,
            ls="",
            ecolor=cmap(colval),
            elinewidth=0.7,
            capsize=0,
            zorder=3,
        )
        ax.set_ylabel("res ({})".format(u[1]), fontsize=ll_size)
        ax.set_xlabel("$\mathrm{BJD_{TDB}}$", fontsize=ll_size)
        axs.append(ax)

        irow += 1
        colval += dmap

    xmin = np.min(xlims)
    xmax = np.max(xlims)
    for ax in axs:
        ax.set_xlim([xmin, xmax])

    for ax in axs[:-1]:
        ax.tick_params(axis="both", labelsize=tl_size, labelbottom=False)
    axs[-1].tick_params(axis="both", labelsize=tl_size, labelbottom=True)

    plt.subplots_adjust(hspace=0.05)
    fig.align_ylabels(axs)
    if output_folder is not None:
        output_file = os.path.join(output_folder, "OC_diagram.png")
        fig.savefig(output_file, dpi=300, bbox_inches="tight")
    if show_plot:
        plt.show()
    plt.close(fig)

    return


def set_photometry_portion(
    time, flux, flux_err, n_oversample=1, t_exp_d=60.0 * cst.sec2day, ancillary=None
):

    ndata = len(time)
    portion = {}
    portion["time"] = time
    portion["flux"] = flux
    portion["flux_err"] = flux_err
    portion["ancillary"] = ancillary
    if ancillary is not None:
        portion["ancillary_interp"] = {
            k: interp1d(time, a, bounds_error=False, fill_value=(a[0],a[-1])) for k, a in ancillary.items()
        }
    else:
        portion["ancillary_interp"] = ancillary
    portion["ndata"] = ndata
    portion["t_exp_d"] = t_exp_d
    portion["n_oversample"] = n_oversample
    time_min = np.min(time)
    time_max = np.max(time)
    portion["time_min"] = time_min
    portion["time_max"] = time_max
    portion["time_med"] = np.median(time)
    t_exp_d_full = min(t_exp_d, 60.0*cst.sec2day)
    portion["t_exp_d_full"] = t_exp_d_full
    portion["time_full"] = np.arange(time_min, time_max+(0.5*t_exp_d_full), t_exp_d_full)
    portion["ndata_full"] = len(portion["time_full"])
    portion["n_oversample_full"] = 1

    return portion


def get_photometry_residuals(obs_photometry, sim_photometry):
    res, err2 = [], []
    for ko, o in obs_photometry.items():
        s = sim_photometry[ko]
        r = {}
        wr = {}
        for kon, op in o.items():
            sp = s[kon]
            rp = op["flux"] - sp
            err = op["flux_err"]
            res = np.concatenate([res, rp])
            err2 = np.concatenate([err2, err * err])

    return res, err2


def get_transit_times_residuals(obs_transits, sim_transits):

    res, err2 = [], []
    for bd, otra in obs_transits.items():
        oT0 = otra["T0s"]
        eT0 = otra["err_T0s"]
        if bd in sim_transits:
            sT0 = sim_transits[bd]
            rT0 = oT0 - sT0
        else:
            rT0 = oT0
        res = np.concatenate([res, rT0])
        err2 = np.concatenate([err2, eT0 * eT0])

    return res, err2


def get_radial_velocities_residuals(
    t_epoch,
    t_rv_obs,
    rv_obs,
    rv_obs_err,
    rv_set_id,
    gammas,
    jitters,
    sim_rv,
    ctrends=[0.0],
):

    trv = t_rv_obs
    orv = rv_obs
    erv = rv_obs_err
    nrv = len(trv)

    tt = trv - t_epoch
    ntrend = len(ctrends)
    trend = 0.0
    if ntrend > 0:
        for i, c in enumerate(ctrends):
            trend += c * (tt ** (i + 1))

    srv = sim_rv["rv"] + trend

    res = np.zeros((nrv))
    jit = np.zeros((nrv))
    err2 = np.zeros((nrv))

    ids = np.unique(rv_set_id)
    # for i, n in enumerate(ids):
    #     sel = rv_set_id == n
    #     gg = gammas[i]
    #     srv[sel] += gg
    #     jit[sel] = jitters[i]
    for n in ids:
        # print(n)
        sel = rv_set_id == n
        gg = gammas[n]
        srv[sel] += gg
        jit[sel] = jitters[n]

    res = orv - srv
    err2 = erv * erv + jit * jit

    return res, err2


class PhotoTRADES:
    def __init__(self, n_body, t_epoch, t_start, t_int, duration_check=1):

        self.n_body = n_body
        self.t_epoch = t_epoch
        self.t_start = t_start
        self.t_int = t_int
        self.duration_check = duration_check

        args_init(
            n_body,
            duration_check,
            t_epoch=t_epoch,
            t_start=t_start,
            t_int=t_int,
        )

        self.Tref = {}
        self.Pref = {}

        self.n_data = 0

        self.transits = {}
        self.n_transits = {}
        self.n_all_transits = 0

        self.photometry = {}
        self.photometry_flux_trend = {}
        self.ancillary_coeff = {}
        self.n_photometry = 0
        self.n_photometry_points = 0
        self.tm = pytransit.QuadraticModel()

        self.radial_velocity = {}
        self.id_radial_velocity = {}
        self.nrv = 0
        self.n_rv_set = 0
        self.t_rv_obs = np.array([], dtype=float)
        self.rv_obs = np.array([], dtype=float)
        self.rv_obs_err = np.array([], dtype=float)
        self.rv_set_id = np.array([], dtype=int)
        self.rv_set_name = []
        self.idx_rv_obs = np.array([], dtype=int)

        return

    def add_default_physical_parameters(
        self, mass, radius, period, ecc, argp, meana, inc, longn
    ):

        self.mass = mass
        self.radius = radius
        self.period = period
        self.ecc = ecc
        self.argp = argp
        self.meana = meana
        self.inc = inc
        self.longn = longn

        return

    def add_linear_ephem(self, body_id, Tref, Pref):

        if isinstance(Tref, (float, int)):
            self.Tref[body_id] = list([Tref])
        else:
            self.Tref[body_id] = list(Tref)
        if isinstance(Pref, (float, int)):
            self.Pref[body_id] = list([Pref])
        else:
            self.Pref[body_id] = list(Pref)

        return

    def add_transits(self, body_id, T0s, err_T0s, epo=None, sources=None):

        n = len(T0s)
        if sources is None:
            sources = np.array(["-"] * n)
        self.n_transits[body_id] = n
        self.transits[body_id] = {
            "nT0:": n,
            "epoch": epo,
            "T0s": T0s,
            "err_T0s": err_T0s,
            "sources": sources,
        }

        self.n_all_transits = np.sum([i for i in self.n_transits.values()])

        return

    def update_linear_ephem(self):

        for bd, tra_bd in self.transits.items():
            nT0 = self.n_transits[bd]
            if len(self.Tref) == 0:
                Tx = tra_bd[nT0 // 2]
            else:
                Tx = self.Tref[bd][0]

            if len(self.Pref) == 0:
                Px = np.median(np.diff(tra_bd))
            else:
                Px = self.Pref[bd][0]

            sort_tra = np.argsort(tra_bd["T0s"])
            T0s = tra_bd["T0s"][sort_tra]
            eT0s = tra_bd["err_T0s"][sort_tra]
            sources = tra_bd["sources"][sort_tra]
            epox = np.rint((T0s - Tx) / Px)

            Tr, Pr, _ = linear_fit(epox, T0s, ey=eT0s)
            epo = np.rint((T0s - Tr[0]) / Pr[0])
            self.add_linear_ephem(bd, Tr, Pr)
            self.add_transits(bd, T0s, eT0s, epo=epo, sources=sources)

        return

    def get_simulated_transits(self, transits, body_flag):

        transits_out = {}
        body_ids = np.unique(body_flag)
        for body_id in body_ids:
            sel = body_id == body_flag
            tra_s = transits[sel]
            epo_s = np.rint(
                (tra_s - self.Tref[body_id][0]) / self.Pref[body_id][0]
            ).astype(int)
            obs = self.transits[body_id]
            tra_obs = obs["T0s"]
            epo_obs = obs["epoch"]
            s_in_o = np.isin(epo_s, epo_obs)
            tra_c = tra_s[s_in_o]
            epo_c = epo_s[s_in_o]
            tra_out = np.zeros_like(tra_obs)
            for i, epo in enumerate(epo_c):
                check = epo == epo_obs
                tra_out[check] = tra_c[i]
                # print(i, epo, check, tra_c[i])
            transits_out[body_id] = tra_out

        return transits_out

    def add_photometry(self, name, data, flux_time_trend, anc_coeff=None):
        self.photometry[name] = data
        self.photometry_flux_trend[name] = flux_time_trend
        if anc_coeff is None:
            self.ancillary_coeff[name] = {k: None for k in data.keys()}
        else:
            self.ancillary_coeff[name] = anc_coeff
        self.n_photometry = len(self.photometry)

        self.n_photometry_points = 0
        for ph in self.photometry.values():
            for d in ph.values():
                self.n_photometry_points += len(d["time"])

        return

    def add_radial_velocity(self, name, data, id_num):
        self.radial_velocity[name] = data
        self.id_radial_velocity[name] = id_num
        self.n_rv_set = len(self.radial_velocity)
        t = data["time"]
        n = len(t)
        self.t_rv_obs = np.concatenate((self.t_rv_obs, t))
        self.rv_obs = np.concatenate((self.rv_obs, data["rv"]))
        self.rv_obs_err = np.concatenate((self.rv_obs_err, data["rv_err"]))
        self.rv_set_id = np.concatenate(
            (self.rv_set_id, np.zeros(n).astype(int) + id_num)
        )
        self.rv_set_name += [name] * n
        self.idx_rv_obs = np.concatenate(
            (self.idx_rv_obs, np.arange(0, n, 1) + self.nrv)
        )
        self.nrv += n
        return

    def set_radial_velocity_sorting(self):

        sort_t_rv = np.argsort(self.t_rv_obs)
        self.t_rv_obs = self.t_rv_obs[sort_t_rv]
        self.rv_obs = self.rv_obs[sort_t_rv]
        self.rv_obs_err = self.rv_obs_err[sort_t_rv]
        self.rv_set_id = self.rv_set_id[sort_t_rv]
        self.rv_set_name = np.array(self.rv_set_name)[sort_t_rv]
        self.idx_rv_obs = self.idx_rv_obs[sort_t_rv]

        return

    def update_n_data(self):

        self.n_data = self.n_photometry_points + self.n_all_transits + self.nrv

        return

    def kelements_to_orbits_full(
        self,
        mass,
        radius,
        period,
        ecc,
        w,
        ma,
        inc,
        long,
        step_size=None,
        n_steps_smaller_orbits=10.0,
    ):

        time_steps, orbits = kelements_to_orbits_full(
            self.t_start,
            self.t_epoch,
            self.t_int,
            mass,
            radius,
            period,
            ecc,
            w,
            ma,
            inc,
            long,
            specific_times=self.t_rv_obs,  # RV TIMES
            step_size=step_size,
            n_steps_smaller_orbits=n_steps_smaller_orbits,
        )

        return time_steps, orbits

    def orbital_parameters_to_transits(
        self, mass, radius, period, ecc, w, ma, inc, long
    ):

        (
            time_steps,
            orbits,
            transits,
            durations,
            kep_elem,
            body_flag,
            rv_sim,
        ) = orbital_parameters_to_transits(
            self.t_start,
            self.t_epoch,
            self.t_int,
            mass,
            radius,
            period,
            ecc,
            w,
            ma,
            inc,
            long,
            self.t_rv_obs,
        )
        return time_steps, orbits, transits, durations, kep_elem, body_flag, rv_sim

    def get_simulate_flux(
        self, radius, ld_quads, transits, durations, body_flag, kep_elem, time_key="time"
    ):

        rp_rs, per, aRs, inc, ecc, w = set_transit_parameters(
            radius, transits, body_flag, kep_elem
        )
        sim_photometry = {}
        for phot_name, phot in self.photometry.items():
            ld_quad = ld_quads[phot_name]
            sim_phot = get_simulate_flux(
                self.tm,
                phot,
                transits,
                durations,
                rp_rs,
                ld_quad,
                per,
                aRs,
                inc,
                ecc,
                w,
                time_key=time_key
            )
            sim_photometry[phot_name] = sim_phot
        return sim_photometry

    def full_photodyn(self, mass, radius, period, ecc, w, ma, inc, long, ld_quads):

        (
            time_steps,
            orbits,
            transits,
            durations,
            kep_elem,
            body_flag,
            rv_sim,
        ) = self.orbital_parameters_to_transits(
            mass, radius, period, ecc, w, ma, inc, long
        )
        sim_photometry = self.get_simulate_flux(
            radius, ld_quads, transits, durations, body_flag, kep_elem, time_key="time"
        )
        if len(self.transits) > 0:
            sim_transits = self.get_simulated_transits(transits, body_flag)
        else:
            sim_transits = {}

        return sim_photometry, rv_sim, sim_transits

    def plot_photometry(
        self,
        sim_photometry,
        mod_photometry=None,
        trend_photometry=None,
        figsize=(3, 3),
        show_plot=True,
        output_folder=None,
        return_rms=False,
    ):

        rms = plot_photometry(
            self.photometry,
            sim_photometry,
            mod_photometry=mod_photometry,
            trend_photometry=trend_photometry,
            figsize=figsize,
            show_plot=show_plot,
            output_folder=output_folder,
            return_rms=return_rms,
        )

        return rms

    def plot_oc(
        self, transits_sim, pl_names, figsize=(8, 4), show_plot=True, output_folder=None
    ):

        phototrades_plot_oc(
            self.transits,
            self.Tref,
            self.Pref,
            transits_sim,
            pl_names,
            figsize=figsize,
            show_plot=show_plot,
            output_folder=output_folder,
        )

        return

    def plot_rv(
        self,
        rv_sim,
        gamma_rv,
        jitter_rv,
        markers,
        rv_trend_coeff=None,
        print_rv=False,
        figsize=(4, 2),
        show_plot=True,
        output_folder=None,
        remove_dataset=None,
    ):

        t_epoch = self.t_epoch
        ms = 3
        fig = plt.figure(figsize=figsize)

        colors = plt.get_cmap("nipy_spectral")(
            np.linspace(0.1, 0.9, endpoint=True, num=len(gamma_rv))
        )

        if remove_dataset is not None and len(remove_dataset) > 0:
            additional_panel = True
        else:
            additional_panel = False
        nrows, ncols = 3, 1
        if additional_panel:
            nrows += 1

        axs = []

        ax1 = plt.subplot2grid((nrows, ncols), (0, 0), rowspan=2)
        ax1.ticklabel_format(useOffset=False)
        ax1.get_xaxis().set_ticks([])
        ax1.tick_params(axis="both", labelsize=plt.rcParams["xtick.labelsize"] - 4)

        ax2 = plt.subplot2grid((nrows, ncols), (2, 0))
        ax2.ticklabel_format(useOffset=False)
        ax2.tick_params(axis="both", labelsize=plt.rcParams["xtick.labelsize"] - 4)

        if additional_panel:
            ax3 = plt.subplot2grid((nrows, ncols), (3, 0))
            ax3.ticklabel_format(useOffset=False)
            ax3.tick_params(axis="both", labelsize=plt.rcParams["xtick.labelsize"] - 4)

        rvs_full = np.zeros_like(self.t_rv_obs)

        for i_k, k in enumerate(self.radial_velocity.keys()):
            kid = self.id_radial_velocity[k]
            gg = gamma_rv[kid]
            sel = k == self.rv_set_name
            t_o = self.t_rv_obs[sel]
            rv_o = self.rv_obs[sel]
            erv_o = self.rv_obs_err[sel]
            t_s = rv_sim["time"][sel]
            rv_s = rv_sim["rv"][sel]

            trend = 0.0
            if rv_trend_coeff is not None:
                w = t_o - t_epoch
                if isinstance(rv_trend_coeff, (float, int)):
                    trend += rv_trend_coeff * w
                else:
                    if len(rv_trend_coeff) > 0:
                        for idx_c, coeff in enumerate(rv_trend_coeff):
                            o = idx_c + 1
                            trend += coeff * (w**o)

            rvs_full[sel] = rv_s + gg + trend

            res = rv_o - gg - trend - rv_s
            eres = np.sqrt(erv_o * erv_o + jitter_rv[kid] * jitter_rv[kid])

            ax1.axhline(0.0, color="black", ls="-", lw=0.7)
            ax1.errorbar(
                t_o,
                rv_o - gg,
                yerr=erv_o,
                color=colors[i_k],
                marker=markers[k],
                ms=ms,
                mec="white",
                mew=0.4,
                ls="",
                elinewidth=0.7,
                capsize=0,
                label=k.upper(),
            )
            ax1.plot(
                t_s,
                rv_s + trend,
                color=colors[i_k],
                marker=markers[k],
                ms=ms + 1.0,
                mfc="None",
                mew=0.4,
                ls="",
            )

            ax2.axhline(0.0, color="black", ls="-", lw=0.7)
            ax2.errorbar(
                t_o,
                res,
                yerr=eres,
                color=colors[i_k],
                marker=markers[k],
                ms=ms,
                mec="white",
                mew=0.4,
                ls="",
                elinewidth=0.7,
                capsize=0,
            )

            if additional_panel:
                ax3.axhline(0.0, color="black", ls="-", lw=0.7)
                if k.lower() not in remove_dataset:
                    ax3.errorbar(
                        t_o,
                        res,
                        yerr=eres,
                        color=colors[i_k],
                        marker=markers[k],
                        ms=ms,
                        mec="white",
                        mew=0.4,
                        ls="",
                        elinewidth=0.7,
                        capsize=0,
                    )

        ax1.legend(
            bbox_to_anchor=(1.05, 0.5), loc="upper left", fontsize=6, frameon=False
        )
        ax1.set_ylabel("RV (m/s) $- \gamma$", fontsize=6)
        xlims = ax1.get_xlim()

        ax2.set_ylabel("res (m/s)", fontsize=6)
        ax2.set_xlim(xlims)

        if additional_panel:
            ax2.get_xaxis().set_ticks([])
            ax3.set_ylabel("res (m/s)", fontsize=6)
            ax3.set_xlim(xlims)
            ax3.set_xlabel("time", fontsize=6)
            axs = [ax1, ax2, ax3]
        else:
            ax2.set_xlabel("time", fontsize=6)
            axs = [ax1, ax2]

        plt.subplots_adjust(hspace=0.05)
        fig.align_ylabels()

        if output_folder is not None:
            output_file = os.path.join(output_folder, "RV_plot.png")
            fig.savefig(output_file, dpi=300, bbox_inches="tight")
        if show_plot:
            plt.show()
        plt.close(fig)

        if print_rv:
            for trvo, rvo, idx, iset, nn, trvs, rvs in zip(
                self.t_rv_obs,
                self.rv_obs,
                self.idx_rv_obs,
                self.rv_set_id,
                self.rv_set_name,
                rv_sim["time"],
                rvs_full,
            ):
                print(
                    "{:11.5f} {:9.2f} {:3d} {:2d} {:10s} {:11.5f} {:9.2f}".format(
                        trvo, rvo, idx, iset, nn, trvs, rvs
                    )
                )

        return

    def get_photometry_residuals(self, sim_photometry):

        res, err2 = get_photometry_residuals(self.photometry, sim_photometry)

        return res, err2

    def get_transit_times_residuals(self, sim_transits):

        res, err2 = get_transit_times_residuals(self.transits, sim_transits)

        return res, err2

    def get_radial_velocities_residuals(self, gammas, jitters, rv_sim, ctrends=[0.0]):

        res, err2 = get_radial_velocities_residuals(
            self.t_epoch,
            self.t_rv_obs,
            self.rv_obs,
            self.rv_obs_err,
            self.rv_set_id,
            gammas,
            jitters,
            rv_sim,
            ctrends=ctrends,
        )

        return res, err2


# =============================================================================
# TRADES: SET-UP FROM FOLDER/FILES
# USED IN trades_emcee.py/trades_emcee_analysis.property
# IT USES THE LOGLIKELIHOOD IN FORTRAN90 AND DEFAULT PARAMETERIZATION
# =============================================================================

class TRADESfolder:
    def __init__(self, input_path, sub_folder="", nthreads=1, seed=42, m_type="e"):

        self.work_path = os.path.join(input_path, "")
        self.sub_folder = sub_folder
        self.full_path = os.path.join(self.work_path, self.sub_folder, "")
        self.nthreads = nthreads
        self.seed = seed
        self.m_type = m_type

        return

    def init_trades_from_path(self, work_path, sub_path):

        print("Initialising TRADES by f90 modules ...")
        sys.stdout.flush()

        self.reset()

        w_path = os.path.join(work_path, "")

        # INITIALISE TRADES WITH SUBROUTINE WITHIN TRADES_LIB -> PARAMETER NAMES, MINMAX, INTEGRATION ARGS, READ DATA ...
        f90trades.initialize_trades(w_path, sub_path, self.nthreads)

        print("... done")
        sys.stdout.flush()

        self.full_path = f90trades.get_path().decode("utf-8").strip()
        self.full_path = os.path.join(os.path.abspath(self.full_path), "")

        self.n_bodies = (
            f90trades.n_bodies
        ).item()  # NUMBER OF TOTAL BODIES OF THE SYSTEM
        self.n_planets = self.n_bodies - 1  # NUMBER OF PLANETS IN THE SYSTEM
        self.ndata = (f90trades.ndata).item()  # TOTAL NUMBER OF DATA AVAILABLE
        self.nkel = (f90trades.nkel).item()  # NUMBER OF FITTING KEPLERIAN ELEMENTS
        self.nfit = (f90trades.nfit).item()  # NUMBER OF PARAMETERS TO FIT
        self.nfree = (f90trades.nfree).item()  # NUMBER OF FREE PARAMETERS (ie nrvset)
        self.dof = (f90trades.dof).item()  # NUMBER OF DEGREES OF FREEDOM = NDATA - NFIT
        self.inv_dof = (f90trades.inv_dof).item()

        self.tstart = (f90trades.tstart).item()
        self.tepoch = (f90trades.tepoch).item()
        self.tint   = (f90trades.tint).item()

        # fitting parameters
        str_len = f90trades.str_len
        temp_names = f90trades.get_parameter_names(self.nfit, str_len)
        # names
        self.fitting_names = anc.convert_fortran_charray2python_strararray(temp_names)
        # fitting values
        self.fitting_parameters = f90trades.fitting_parameters.copy()
        # boundaries
        self.fitting_minmax = f90trades.parameters_minmax.copy()

        # all system_parameters
        self.system_parameters = f90trades.system_parameters.copy()

        # priors
        self.n_priors = (f90trades.n_priors).item()

        # rv trend
        self.rv_trend_order = (f90trades.rv_trend_order).item()

        # Stellar Mass and Radius
        self.MR_star = f90trades.mr_star.copy()
        # Mass conversion factor for output
        # m_factor_0, self.mass_unit = anc.mass_type_factor(
        #     Ms=1.0, mtype=self.m_type, mscale=False
        # )
        (
            m_factor_0,
            self.mass_unit,
            r_factor_0,
            self.radius_unit,
        ) = anc.mass_radius_type_factor(mtype=self.m_type)
        self.m_factor = m_factor_0 * self.MR_star[0, 0]
        self.m_factor_post = self.m_factor  # temporary
        self.r_factor = r_factor_0 * self.MR_star[1, 0]
        self.r_factor_post = self.r_factor  # temporary

        # units of the fitted parameters
        self.fitting_units = anc.get_units(self.fitting_names, self.mass_unit)

        # get information on how to decompose fitted parameters in physical/derived ones
        # (
        #     _,
        #     _,
        #     self.bodies_file,
        #     self.id_fit,
        #     self.id_all,
        #     self.nfit_list,
        #     self.cols_list,
        #     self.case,
        # ) = anc.get_fitted(self.work_path)

        # set kepler elements, divided by names ...
        self.mass = np.zeros((self.n_bodies))
        (
            self.radius,
            self.period,
            self.sma,
            self.ecc,
            self.argp,
            self.meanA,
            self.inc,
            self.longN,
        ) = (
            np.zeros_like(self.mass),
            np.zeros_like(self.mass),
            np.zeros_like(self.mass),
            np.zeros_like(self.mass),
            np.zeros_like(self.mass),
            np.zeros_like(self.mass),
            np.zeros_like(self.mass),
            np.zeros_like(self.mass),
        )
        (
            self.mass,
            self.radius,
            self.period,
            self.sma,
            self.ecc,
            self.argp,
            self.meanA,
            self.inc,
            self.longN,
            self.checkpar,
        ) = f90trades.convert_trades_par_to_kepelem(
            self.system_parameters, self.fitting_parameters, self.n_bodies
        )
        # and in array, without the start!
        # M_Msun R_Rsun P_d a_AU e w_deg mA_deg inc_deg lN_deg
        self.kep_elem = np.zeros((self.n_planets, 9))
        self.kep_elem[:, 0] = self.mass[1:]
        self.kep_elem[:, 1] = self.radius[1:]
        self.kep_elem[:, 2] = self.period[1:]
        self.kep_elem[:, 3] = self.sma[1:]
        self.kep_elem[:, 4] = self.ecc[1:]
        self.kep_elem[:, 5] = self.argp[1:]
        self.kep_elem[:, 6] = self.meanA[1:]
        self.kep_elem[:, 7] = self.inc[1:]
        self.kep_elem[:, 8] = self.longN[1:]

        # computes the derived/physical parameters and their names
        (
            self.physical_names,
            self.physical_parameters,
            self.physical_posterior,
            self.phys_type,
            self.physical_chains,
            self.physical_chains_posterior,
        ) = anc.compute_physical_parameters(
            self.n_bodies,
            self.fitting_parameters,
            self.fitting_names,
            self.system_parameters,
            mass_conv_factor=self.m_factor,
            radius_conv_factor=self.r_factor,
        )
        # units of the derived/physical parameters
        self.physical_units = anc.get_units(self.physical_names, self.mass_unit)
        self.nphy = len(self.physical_parameters)

        # RADIAL VELOCITIES SET
        self.n_rv = f90trades.nrv
        self.n_set_rv = f90trades.nrvset  # number of jitter parameters

        # TRANSITS SET
        self.n_t0 = f90trades.nt0
        self.n_t0_sum = f90trades.ntts
        self.n_set_t0 = 0
        for i in range(0, self.n_planets):
            if self.n_t0[i] > 0:
                self.n_set_t0 += 1

        self.pephem = f90trades.pephem
        self.tephem = f90trades.tephem
        self.wrttime = f90trades.wrttime

        return

    def init_trades(self):

        self.init_trades_from_path(self.work_path, self.sub_folder)

        return

    def set_one_fit_par_boundaries(self, ifit, min_val, max_val):

        self.fitting_minmax[ifit, 0] = min_val
        self.fitting_minmax[ifit, 1] = max_val
        set_one_fit_par_boundaries(ifit+1, min_val, max_val)

        return
    
    def reset_all_fit_boundaries(self):

        reset_all_fit_boundaries()
        self.fitting_minmax = f90trades.parameters_minmax.copy()

        return

    def path_change(self, new_path):

        f90trades.path_change(os.path.join(new_path, ""))

        return

    def write_time_change(self, new_wrttime):

        f90trades.wrttime = new_wrttime
        self.wrttime = f90trades.wrttime

        return

    def update_system_fitting_parameters_from_keplerian_input(
        self,
        mass,
        radius,
        period,
        ecc,
        argp,
        meanA,
        inc,
        longN,
    ):

        n_all_par = len(self.system_parameters)
        all_par = np.zeros((n_all_par))
        fit_par = np.zeros((self.nfit))
        (all_par, fit_par,) = f90trades.update_parameters_from_keplerian(
            mass,
            radius,
            period,
            ecc,
            argp,
            meanA,
            inc,
            longN,
        )
        self.system_parameters = all_par.copy()
        self.fitting_parameters = fit_par.copy()

        return

    def update_system_fitting_parameters(self):

        self.update_system_fitting_parameters_from_keplerian_input(
            self.mass,
            self.radius,
            self.period,
            self.ecc,
            self.argp,
            self.meanA,
            self.inc,
            self.longN,
        )

        return

    def run_and_get_stats_from_parameters(self, fit_pars):

        (
            chi_square,
            reduced_chi_square,
            lgllhd,
            lnprior,
            ln_const,
            bic,
            check,
        ) = f90trades.fortran_fitness_function(np.asarray(fit_pars, dtype=float))

        return chi_square, reduced_chi_square, lgllhd, lnprior, ln_const, bic, check

    def run_and_get_stats(self):

        (
            chi_square,
            reduced_chi_square,
            lgllhd,
            lnprior,
            ln_const,
            bic,
            check,
        ) = self.run_and_get_stats_from_parameters(self.fitting_parameters)

        return chi_square, reduced_chi_square, lgllhd, lnprior, ln_const, bic, check

    def get_all_observables_from_parameters(self, fit_par):

        # run TRADES to get all TTs and RVs during the integration time
        nt0_full, nrv_nmax = f90trades.get_ntt_nrv(fit_par)
        # anc.print_both("full nRV = {:d} and nT0s = {}".format(nrv_nmax, nt0_full))
        (
            ttra_full,
            dur_full,
            id_ttra_full,
            stats_ttra,
            time_rv_nmax,
            rv_nmax,
            stats_rv,
        ) = f90trades.fit_par_to_ttra_rv(fit_par, nt0_full, nrv_nmax)

        return (
            ttra_full,
            dur_full,
            id_ttra_full,
            stats_ttra,
            time_rv_nmax,
            rv_nmax,
            stats_rv,
        )

    def save_models_from_parameters(self, fit_par, smp_h5, smp_name):

        (
            ttra_full,
            dur_full,
            id_ttra_full,
            stats_ttra,
            time_rv_nmax,
            rv_nmax,
            stats_rv,
        ) = self.get_all_observables_from_parameters(fit_par)

        gr = smp_h5.create_group(smp_name)
        gr.create_dataset(
            "fitting_parameters", data=fit_par, dtype=np.float64, compression="gzip"
        )
        gr["fitting_parameters"].attrs["fitting_names"] = self.fitting_names
        gr.attrs["tepoch"] = self.tepoch

        stats_rv = np.array(stats_rv).astype(bool)
        time_rvx, rvx = time_rv_nmax[stats_rv], rv_nmax[stats_rv]
        id_rv = np.argsort(time_rvx)
        time_rv, rv = time_rvx[id_rv], rvx[id_rv]

        # anc.print_both("time_rv min = {:.5f} max = {:.5f}".format(np.min(time_rv), np.max(time_rv)))
        gr.create_dataset(
            "time_rv_mod", data=time_rv, dtype=np.float64, compression="gzip"
        )
        gr.create_dataset("rv_mod", data=rv, dtype=np.float64, compression="gzip")

        stats_ttra = np.array(stats_ttra).astype(bool)
        for ipl in range(2, self.n_bodies + 1):
            sel_tra = np.logical_and(
                np.logical_and(id_ttra_full == ipl, stats_ttra), ttra_full > -9.0e10
            )
            ttra = np.array(ttra_full)[sel_tra]
            dur = np.array(dur_full)[sel_tra]
            idx_tra = np.argsort(ttra)
            ttra = ttra[idx_tra]
            dur = dur[idx_tra]
            gr.create_dataset(
                "TTs_{:d}".format(ipl), data=ttra, dtype=np.float64, compression="gzip"
            )
            gr.create_dataset(
                "T41s_{:d}".format(ipl), data=dur, dtype=np.float64, compression="gzip"
            )
            # if np.sum(sel_tra) > 0:
            #     anc.print_both("T0_{:02d} min = {:.5f} max = {:.5f}".format(ipl, np.min(ttra), np.max(ttra)))

        return

    def run_and_write_summary_files(
        self,
        sim_name,
        fit_pars,
        phy_pars,
        id_sim=1,
        sigma_fit=None,
        sigma_phy=None,
        par_description="",
        return_stats=False,
    ):

        print("RUNNING sim {}".format(sim_name))
        # print("with parameters:")
        # print("{}".format(fit_pars))

        out_folder = os.path.join(self.full_path, sim_name, "")
        os.makedirs(out_folder, exist_ok=True)
        f90trades.path_change(out_folder)

        for iname, name in enumerate(self.fitting_names):
            if (name[0] == "w"
                or name[0:2] == "mA"
                or name[0:2] == "lN"
                or "lambda" in name
            ):
                fit_pars[iname] %= 360.0
                
        (
            chi_square,
            reduced_chi_square,
            lgllhd,
            lnprior,
            ln_const,
            bic,
            check,
        ) = f90trades.write_summary_files(id_sim, fit_pars)

        # smp_name = "{:04d}_{:s}".format(id_sim, sim_name)
        smp_name = sim_name
        smp_file = os.path.join(out_folder, "{}_models.hdf5".format(smp_name))
        with h5py.File(smp_file, "w") as smp_h5:

            self.save_models_from_parameters(fit_pars, smp_h5, smp_name)

            gr = smp_h5[smp_name]
            gr.attrs["chi_square"] = chi_square
            gr.attrs["reduced_chi_square"] = reduced_chi_square
            gr.attrs["lgllhd"] = lgllhd
            gr.attrs["lnprior"] = lnprior
            gr.attrs["ln_const"] = ln_const
            gr.attrs["bic"] = bic

        mass = np.zeros_like(self.mass)
        radius, period, sma, ecc, argp, meanA, inc, longN = (
            np.zeros_like(mass),
            np.zeros_like(mass),
            np.zeros_like(mass),
            np.zeros_like(mass),
            np.zeros_like(mass),
            np.zeros_like(mass),
            np.zeros_like(mass),
            np.zeros_like(mass),
        )

        (
            mass,
            radius,
            period,
            sma,
            ecc,
            argp,
            meanA,
            inc,
            longN,
            checkpar,
        ) = f90trades.convert_trades_par_to_kepelem(
            self.system_parameters, fit_pars, self.n_bodies
        )
        # M_Msun R_Rsun P_d a_AU e w_deg mA_deg inc_deg lN_deg
        kep_elem = np.zeros((self.n_planets, 9))
        kep_elem[:, 0] = mass[1:]
        kep_elem[:, 1] = radius[1:]
        kep_elem[:, 2] = period[1:]
        kep_elem[:, 3] = sma[1:]
        kep_elem[:, 4] = ecc[1:]
        kep_elem[:, 5] = argp[1:]
        kep_elem[:, 6] = meanA[1:]
        kep_elem[:, 7] = inc[1:]
        kep_elem[:, 8] = longN[1:]

        top_header, header = anc.get_header(anc.percentile_val)
        out_file = os.path.join(out_folder, "parameters_summary.txt")
        with open(out_file, "w") as out:

            if sigma_fit is None:
                sigma_print = np.zeros((8, self.nfit))
            else:
                sigma_print = sigma_fit
            anc.print_both("", out)
            anc.print_both("# " + "=" * 40, out)
            anc.print_both(
                "\n# SUMMARY: {:s} ==> {:s}".format(sim_name, par_description), out
            )
            anc.print_both("# FITTED PARAMETERS", out)
            anc.print_parameters(
                top_header,
                header,
                self.fitting_names,
                self.fitting_units,
                fit_pars,
                sigma_parameters=sigma_print,
                output=out,
            )

            if sigma_phy is None:
                sigma_print = np.zeros((8, self.nphy))
            else:
                sigma_print = sigma_phy
            anc.print_both("# DERIVED PARAMETERS", out)
            anc.print_parameters(
                top_header,
                header,
                self.physical_names,
                self.physical_units,
                phy_pars,
                sigma_parameters=sigma_print,
                output=out,
            )
            if not bool(check):
                anc.print_both(
                    "WRITING WARNING FILE: {:s}".format(
                        os.path.join(out_folder, "WARNING.txt")
                    ),
                    out,
                )
                warn_o = open(os.path.join(out_folder, "WARNING.txt"), "w")
                warn_o.write(
                    "*******\nWARNING: FITTED PARAMETERS COULD NOT BE PHYSICAL!\nWARNING: BE VERY CAREFUL WITH THIS PARAMETER SET!\n*******"
                )
                warn_o.close()

            anc.print_both("", out)
            anc.print_both("# " + "=" * 40, out)
            anc.print_both("", out)

        f90trades.path_change(self.full_path)

        if return_stats:
            return chi_square, reduced_chi_square, lgllhd, lnprior, ln_const, bic, check

        return

    def reset(self):

        f90trades.deallocate_variables()

        return
