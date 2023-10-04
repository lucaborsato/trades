#!/usr/bin/env python
# -*- coding: utf-8 -*-

## TO DO LIST
# 0) simulate the grid, starting from the t_epoch to end of the CHEOPS mission
# A)
# 1) read grid_summary file with fitted/grid sim_id,Mc,Pc,ec,ic & fitness: 'summary_grid_sims_0.dat'
# 2) read all the TTs from sim_id_0_NB2_simT0.dat
# 3) provide linear ephemeris of wasp-8b ... or what ==> T0lin
# 4) compute O-C for each sim_id and amplitude rms_oc = sqrt( sum(x^2)/nx ) = sqrt( <x^2> )
# 5) plot correlation plot (triangle plot) of Mc vs Pc vs ec vs ic with rms_oc as colorbar
#    should I smooth (kernel or what) stuff or what to plot them??
# 6) overplot (may in a new plot) the rms_oc above the threshold == err_T0s (decide how: gray scale below, viridis scale above?)
# 7) save file with sims id, parameters, fitness and rms_oc
#
# B)
# 1) Define N_TT during cheops mission --> N_TT = [3, 5, 6, 10, 15]
# 2) for each sim_id with rms_oc above err_T0s:
#    a) for each N_TT value:
#      i) select N_TT transits for N_mc times
#      ii) compute rms_oc for each N_mc
#      iii) compute average rms_oc of mc: rms_oc_mc
# 3) new triangle plot Mc vs Pc vs ec vs ic with rms_oc_mc as colorbar (showing above err_T0s)
# 4) save file with sims id, parameters, fitness, rms_oc, rms_oc_mc

# ==============================================================================
# IMPORT

# no more "zero" integer division bugs!:P
import argparse
import numpy as np  # array
import h5py
import astropy.time as atime
import os
import sys

# import emcee.interruptible_pool as emcee_pool
import gc

import gls

# import pp

import matplotlib as mpl

# from matplotlib import use as mpluse
mpl.use("Agg", warn=False)
# mpluse("Qt4Agg")
# mpluse("TkAgg")
import matplotlib.pyplot as plt

# matplotlib rc params
plt.rcParams["text.usetex"] = True
plt.rcParams["font.family"] = "serif"
plt.rcParams["figure.figsize"] = [6, 6]
plt.rcParams["figure.facecolor"] = "white"
plt.rcParams["savefig.facecolor"] = "white"
plt.rcParams["figure.dpi"] = 200
plt.rcParams["savefig.dpi"] = 300
plt.rcParams["font.size"] = 10

# import matplotlib.mlab as mlab
import matplotlib.cm as cm

# import matplotlib.colors as mc
# from mpl_toolkits.axes_grid1 import make_axes_locatable
from matplotlib.ticker import ScalarFormatter

script_path = os.path.realpath(__file__)
module_path = os.path.abspath(os.path.join(os.path.dirname(script_path), "../pytrades"))
sys.path.append(module_path)

import ancillary as anc
import constants as cst
from pytrades_lib import f90trades
import pytrades

# import emcee.interruptible_pool as emcee_pool

import multiprocessing as mpro

# from multiprocessing.pool import ThreadPool

# import threading
import _thread

import time

# ==============================================================================

# read cli arguments
def get_args():

    parser = argparse.ArgumentParser(description="CHEOPS TTV LEVEL ANALYSIS")
    parser.add_argument(
        "-p",
        "--path",
        action="store",
        dest="full_path",
        required=True,
        help="The path (absolute or relative) with simulation files for TRADES.",
    )

    # parser.add_argument('-lm', '--lm-flag', action='store', dest='lm_flag', default='0', help='LM flag. If you want to read files w/o LM set to 0, otherwise 1. Default is 0')

    # SUBFOLDER TO SAVE ALL SIMULATION DATA
    parser.add_argument(
        "-sf",
        "--sub-folder",
        "--sub-f",
        action="store",
        dest="sub_folder",
        default="grid2ttv",
        help="Sub-folder name, without full path. Default = grid2ttv",
    )

    parser.add_argument(
        "-e",
        "--e-T0",
        "--err-T0",
        action="store",
        dest="err_T0",
        default="0.0",
        help="Mean uncertainty on the Transit Times (T0s)."
        + 'Provide value with unit "d" "m" "s" for days, minutes, and seconds in this way:'
        + '"0.001528d" or "2.2m" or "132.s". Default unit is days and value set to 0.',
    )

    parser.add_argument(
        "-k",
        "--K",
        "--K-rv",
        action="store",
        dest="K_rv",
        default=None,
        help="Radial Velocity semi amplitude (K) threshold in m/s (no check is performed). "
        + "Only the simulations with rms_RV <= K will be taken into account for plots. "
        + "The rms_Rv is computed as the rms of the RV calculated during the orbital integration",
    )

    parser.add_argument(
        "-mr",
        "--mr-type",
        "--mass-radius-type",
        action="store",
        dest="mr_type",
        default="e",
        help="Mass and radius type to convert the grid. "
        + '"e ea ear eart earth" or "j ju jup jupi jupit jupite jupiter" or '
        + '"n ne nep nept neptu neptun neptune".',
    )

    parser.add_argument(
        "-s",
        "--seed",
        action="store",
        dest="seed",
        default="None",
        help="Seed for random number generator. Defauto is None.",
    )

    parser.add_argument(
        "-c",
        "--cpu",
        "--threads",
        action="store",
        dest="nthreads",
        default=1,
        type=int,
        help="Number of threads/cpu to use. Default is 1.",
    )

    parser.add_argument(
        "-id",
        "--id-perturber",
        action="store",
        dest="perturber_id",
        default=3,
        type=int,
        help="Id of the perturber body: star is 1, so it has to be from 2 to Number of bodies. Default is 3.",
    )

    parser.add_argument(
        "-P-tra",
        "--P-tra",
        "-P-transiting",
        "--P-transiting",
        action="store",
        dest="period_transiting",
        default=None,
        help="Period of the transiting body. If provided the internal and external MMR of "
        + "2:1, 3:1, 4:1, 2:3, 3:4, 4:5 will be plotted. Default is None.",
    )

    parser.add_argument(
        "-ts",
        "--ts",
        "-time-start",
        "--time-start",
        action="store",
        dest="t_start",
        default=0.0,
        help="Enter starting time of selection for the plot in same unit of the simulation. Default is 0.0",
    )

    parser.add_argument(
        "-te",
        "--te",
        "-time-end",
        "--time-end",
        action="store",
        dest="t_end",
        default=365.25,
        help="Enter ending time of selection for the plot in same unit of the simulation. Default is 365.25",
    )

    cli = parser.parse_args()
    cli.full_path = os.path.join(os.path.abspath(cli.full_path), "")
    cli.sub_folder = os.path.join(os.path.relpath(cli.sub_folder), "")
    # cli.lm_flag = str(cli.lm_flag)
    try:
        cli.seed = int(cli.seed)
    except:
        cli.seed = None

    try:
        cli.perturber_id = int(cli.perturber_id)
    except:
        cli.perturber_id = 3

    try:
        cli.K_rv = np.float64(cli.K_rv)
    except:
        cli.K_rv = None

    try:
        cli.period_transiting = np.float64(cli.period_transiting)
    except:
        # cli.K_rv = None
        cli.period_transiting = None

    try:
        cli.t_start = np.float64(cli.t_start)
    except:
        cli.t_start = 0.0

    try:
        cli.t_end = np.float64(cli.t_end)
    except:
        cli.t_end = 365.25

    return cli


# -------------------------------------


def set_err_T0(err_T0_in):

    if "d" in err_T0_in:
        err_T0 = np.float64(err_T0_in.split("d")[0])
    elif "m" in err_T0_in:
        err_T0 = np.float64(err_T0_in.split("m")[0]) / 1440.0
    elif "s" in err_T0_in:
        err_T0 = np.float64(err_T0_in.split("s")[0]) / 86400.0
    else:
        anc.print_both("WARNING: wrong mean err_T0: set to 30s")
        err_T0 = 30.0 / 86400.0

    return err_T0


# -------------------------------------


def set_analysis_folder_log(full_path, sub_folder):

    out_folder = os.path.join(full_path, sub_folder)
    if not os.path.isdir(out_folder):
        os.makedirs(out_folder)

    anc.copy_simulation_files(full_path, out_folder)

    log_file = os.path.join(out_folder, "analysis_log.txt")
    of_log = open(log_file, "w")

    return out_folder, log_file, of_log


# -------------------------------------


def mass_factor(mr_type):

    if mr_type.lower() in "j ju jup jupi jupit jupite jupiter".split():
        # mr_convert = Msun2jup
        mr_convert = [cst.Msjup, cst.Rsjup]
        mr_unit = ["M_jup", "R_jup"]
    elif mr_type.lower() in "n ne nep nept neptu neptun neptune".split():
        # mr_convert = Msun2nep
        mr_convert = [cst.Msnep, cst.Rsnep]
        mr_unit = ["M_nep", "R_nep"]
    else:
        # mr_convert = Msun2ear
        mr_convert = [cst.Msear, cst.Rsear]
        mr_unit = ["M_ear", "R_ear"]

    return mr_convert, mr_unit


# -------------------------------------


def set_names_parameters(id_body, mr_unit, of_log):

    body_num = "{0:d}".format(id_body)
    base_names_lst = "M R P a e w mA tau i lN".split()
    units = "{0:s} {1:s} day au - deg deg day deg deg".format(mr_unit[0], mr_unit[1])
    units = units.split()
    n_names = len(base_names_lst)
    names_lst = [
        "{0:s}{1:s}".format(base_names_lst[i], body_num) for i in range(0, n_names)
    ]
    temp_lst = [
        "{0:s} [{1:s}]".format(names_lst[i], units[i].strip())
        for i in range(0, n_names)
    ]
    names_str = " ".join(["{0:23s}".format(temp_lst[i]) for i in range(0, n_names)])

    return names_lst, names_str, units


# -------------------------------------


def select_and_sort_ttra(
    ttra_full,
    id_ttra_full,
    stats_ttra,
    id_transiting=2,
    t_min_sel=None,
    t_max_sel=None,
    of_log=None,
):

    sel_tra = np.logical_and(id_ttra_full == id_transiting, stats_ttra)
    ttra_all = np.sort(ttra_full[sel_tra])
    nttra_all = np.shape(ttra_all)[0]
    # anc.print_both(' | number of transits nT0 = %d' %(nttra_all), of_log)

    if nttra_all == 0:
        return None, 0, None, 0, None, 0

    if t_min_sel is None:
        sel_min = 0.0
    else:
        sel_min = t_min_sel
    if t_max_sel is None:
        sel_max = ttra_all[-1]
    else:
        sel_max = t_max_sel

    sel_pre = (
        ttra_all < sel_min
    )  # T0s before CHEOPS mission, they are use to define Tref, Pref

    nttra_pre = int(np.sum(np.asarray(sel_pre, dtype=np.int)))
    # anc.print_both(' | selected pre-CHEOPS transits: nT0 (pre) = %d' %(nttra_pre), of_log)
    if nttra_pre == 0:
        return None, 0, None, 0, None, 0

    sel_cheops = np.logical_and(ttra_all >= sel_min, ttra_all <= sel_max)
    nttra_cheops = int(np.sum(np.asarray(sel_cheops, dtype=np.int)))
    # anc.print_both(' | selected CHEOPS transits: nT0 (cheops) = %d' %(nttra_cheops), of_log)
    if nttra_cheops == 0:
        return None, 0, None, 0, None, 0

    return ttra_all, nttra_all, sel_pre, nttra_pre, sel_cheops, nttra_cheops


# -------------------------------------


def select_and_sort_rv(time_rv_nmax, rv_nmax, stats_rv):

    time_rv_tosort = time_rv_nmax[stats_rv].copy()
    rv_tosort = rv_nmax[stats_rv].copy()
    id_s = np.argsort(time_rv_tosort)
    time_rv = time_rv_tosort[id_s]
    rv = rv_tosort[id_s]

    return time_rv, rv


# -------------------------------------


def dirty_transits(T0, err_T0):

    sig_err = 0.1 * err_T0
    min_err = 0.5 * err_T0
    nT0 = np.shape(T0)[0]

    dirty = []
    ierr = 0
    while True:
        # eee = np.random.normal(loc=0., scale=sig_err)
        eee = np.random.normal(loc=err_T0, scale=sig_err)
        if np.abs(eee) > min_err:
            dirty.append(eee)
            ierr += 1
            if ierr == nT0:
                break

    # eT0 = np.array(dirty, dtype=np.float64) + err_T0
    eT0 = np.abs(np.array(dirty, dtype=np.float64))
    eT0[eT0 == 0.0] = err_T0
    T0dirty = T0 + np.array(dirty, dtype=np.float64)

    return T0dirty, eT0


# -------------------------------------


def rms(vec):

    rms_vec = np.sqrt(np.mean(vec * vec))

    return rms_vec


# -------------------------------------


def compute_oc(Tref, Pref, T0):

    epo = anc.calculate_epoch(T0, Tref, Pref)
    T0lin = Tref + epo * Pref
    oc = T0 - T0lin
    # rms_oc = np.sqrt(np.mean(oc*oc))
    rms_oc = rms(oc)

    return epo, oc, rms_oc


# -------------------------------------


def compute_rms_mc(oc, nT0_sel, nMC):

    if np.shape(oc)[0] < nT0_sel:
        return 0.0, 0.0
    rms_all = [
        rms(np.random.choice(oc, nT0_sel, replace=False)) for imc in range(0, nMC)
    ]
    rms_mean = np.mean(rms_all)
    rms_med = np.percentile(rms_all, 50.0, interpolation="midpoint")

    return rms_mean, rms_med


# -------------------------------------


def do_chi2r(res_in, err_in, dof):

    if np.shape(res_in)[0] == 0 or np.shape(err_in)[0] == 0:

        chi2r = -1.0

    else:

        if np.any(err_in == 0.0):
            err0 = err_in != 0.0
            err = err_in[err0]
            res = res_in[err0]
            werr = res / err
        else:
            werr = res_in / err_in

        chi2r = np.sum(np.power(werr, 2)) / np.float64(dof)

    return chi2r


# -------------------------------------


def compute_chi2r_mc(oc, eoc, nT0_sel, nMC):

    noc = np.shape(oc)[0]
    dof = nT0_sel - 2
    if noc < nT0_sel:
        return -1.0, -1.0
    if dof <= 0:
        print("WARNING dof <= 0! Set to dof = 1")
        dof = 1
    for _ in range(0, nMC):
        id_sel = np.sort(np.random.choice(noc, nT0_sel, replace=False))
        chi2r_all = do_chi2r(oc[id_sel], eoc[id_sel], dof)
    chi2r_mean = np.mean(chi2r_all)
    chi2r_med = np.percentile(chi2r_all, 50.0, interpolation="midpoint")

    return chi2r_mean, chi2r_med


# -------------------------------------

# ==============================================================================

# -------------------------------------


def create_hdf5_file(
    full_path,
    seed,
    id_body,
    perturber_grid,
    mr_type,
    t_sel_min,
    t_sel_max,
    err_T0,
    of_log=None,
):

    mr_convert, mr_unit = mass_factor(mr_type)
    names_lst, names_str, units = set_names_parameters(id_body, mr_unit, of_log)
    ngrid, _ = np.shape(perturber_grid)
    write_grid = perturber_grid.copy()
    write_grid[:, 0] = write_grid[:, 0] * mr_convert[0]  # mass
    write_grid[:, 1] = write_grid[:, 1] * mr_convert[1]  # radius

    file_name = os.path.join(full_path, "rms_analysis_summary.hdf5")
    f_h5 = h5py.File(file_name, "w")

    f_h5.create_dataset(
        "grid/grid_names",
        data=anc.encode_list(names_lst),
        dtype="S20",
        compression="gzip",
    )
    f_h5.create_dataset(
        "grid/grid_units", data=anc.encode_list(units), dtype="S20", compression="gzip"
    )
    f_h5.create_dataset(
        "grid/perturber_grid", data=write_grid, dtype=np.float64, compression="gzip"
    )

    f_h5["grid"].attrs["ngrid"] = ngrid
    f_h5["grid"].attrs["seed"] = seed
    f_h5["grid"].attrs["t_sel_min"] = t_sel_min
    f_h5["grid"].attrs["t_sel_max"] = t_sel_max
    f_h5["grid"].attrs["err_T0_day"] = err_T0
    f_h5["grid"].attrs["err_T0_min"] = err_T0 * 1440.0

    return file_name, f_h5, mr_convert, mr_unit, names_lst, names_str, units


# -------------------------------------


def save_results_sim(f_h5, results_sim, mr_type):

    mr_convert, _ = mass_factor(mr_type)
    temp_par = results_sim["perturber_par"].copy()
    temp_par[0] = temp_par[0] * mr_convert[0]
    temp_par[1] = temp_par[1] * mr_convert[1]

    sim_num = results_sim["sim_num"]
    f_h5.create_dataset(
        "sim_{0:06d}/perturber_par_names".format(sim_num),
        data=anc.encode_list(np.array(results_sim["perturber_par_names"], dtype="S10")),
        dtype="S20",
        compression="gzip",
    )
    f_h5.create_dataset(
        "sim_{0:06d}/perturber_par_units".format(sim_num),
        data=anc.encode_list(np.array(results_sim["perturber_par_units"], dtype="S10")),
        dtype="S20",
        compression="gzip",
    )
    f_h5.create_dataset(
        "sim_{0:06d}/perturber_par".format(sim_num),
        data=temp_par,
        dtype=np.float64,
        compression="gzip",
    )
    f_h5.create_dataset(
        "sim_{0:06d}/MR_sun2unit".format(sim_num),
        data=np.array(mr_convert),
        dtype=np.float64,
        compression="gzip",
    )

    f_h5.create_dataset(
        "sim_{0:06d}/keplerian_elements_header".format(sim_num),
        data=anc.encode_list(
            np.array(results_sim["keplerian_elements_header"], dtype="S10")
        ),
        dtype="S20",
        compression="gzip",
    )
    f_h5.create_dataset(
        "sim_{0:06d}/keplerian_elements".format(sim_num),
        data=results_sim["keplerian_elements"],
        dtype=np.float64,
        compression="gzip",
    )

    # f_h5.create_dataset('sim_%06d/ttra_all'%(sim_num), data=results_sim['ttra_all'], dtype=np.float64, compression="gzip")
    # f_h5.create_dataset('sim_%06d/sel_pre_cheops'%(sim_num), data=results_sim['sel_pre_cheops'], dtype=bool, compression="gzip")
    # f_h5.create_dataset('sim_%06d/sel_cheops'%(sim_num), data=results_sim['sel_cheops'], dtype=bool, compression="gzip")
    # f_h5.create_dataset('sim_%06d/T0_eT0'%(sim_num), data=results_sim['T0_eT0'], dtype=np.float64, compression="gzip")
    sel_cheops = results_sim["sel_cheops"]
    if np.sum(sel_cheops.astype(int)) > 0:
        f_h5.create_dataset(
            "sim_{0:06d}/T0_eT0".format(sim_num),
            data=results_sim["T0_eT0"][sel_cheops, :],
            dtype=np.float64,
            compression="gzip",
        )
    else:
        f_h5.create_dataset(
            "sim_{0:06d}/T0_eT0".format(sim_num),
            data=np.zeros((2)),
            dtype=np.float64,
            compression="gzip",
        )

    # f_h5.create_dataset('sim_%06d/time_rv'%(sim_num), data=results_sim['time_rv'], dtype=np.float64, compression="gzip")
    f_h5.create_dataset(
        "sim_{0:06d}/rms_rv".format(sim_num),
        data=np.array([results_sim["rms_rv"]]),
        dtype=np.float64,
        compression="gzip",
    )

    f_h5.create_dataset(
        "sim_{0:06d}/linear_ephem".format(sim_num),
        data=results_sim["linear_ephem"],
        dtype=np.float64,
        compression="gzip",
    )

    f_h5.create_dataset(
        "sim_{0:06d}/rms_ttv".format(sim_num),
        data=np.array([results_sim["rms_ttv"]]),
        dtype=np.float64,
        compression="gzip",
    )
    f_h5.create_dataset(
        "sim_{0:06d}/nT0_sel".format(sim_num),
        data=results_sim["nT0_sel"],
        dtype=np.int,
        compression="gzip",
    )
    f_h5.create_dataset(
        "sim_{0:06d}/rms_mc".format(sim_num),
        data=results_sim["rms_mc"],
        dtype=np.float64,
        compression="gzip",
    )

    f_h5.create_dataset(
        "sim_{0:06d}/chi2r_cheops".format(sim_num),
        data=np.array([results_sim["chi2r_cheops"]]),
        dtype=np.float64,
        compression="gzip",
    )
    f_h5.create_dataset(
        "sim_{0:06d}/chi2r_mc".format(sim_num),
        data=results_sim["chi2r_mc"],
        dtype=np.float64,
        compression="gzip",
    )

    f_h5["sim_{0:06d}".format(sim_num)].attrs["sim_num"] = results_sim["sim_num"]
    f_h5["sim_{0:06d}".format(sim_num)].attrs["nMC"] = results_sim["nMC"]
    f_h5["sim_{0:06d}".format(sim_num)].attrs["A_TTV"] = results_sim["A_TTV"]
    f_h5["sim_{0:06d}".format(sim_num)].attrs["P_TTV"] = results_sim["P_TTV"]

    return f_h5


# -------------------------------------


def save_rms_analysis(f_h5, rms_grid, above_threshold, nT0_sel, rms_rv):

    rms_header_str = "rms_ttv_d {0:s}".format(
        " ".join(
            [
                "rms_mean_N{0:d}_d rms_median_N{1:d}_d".format(
                    nT0_sel[i_n], nT0_sel[i_n]
                )
                for i_n in range(len(nT0_sel))
            ]
        )
    )
    rms_header_lst = rms_header_str.split()

    f_h5.create_dataset(
        "grid/rms_grid", data=rms_grid, dtype=np.float64, compression="gzip"
    )

    f_h5.create_dataset(
        "grid/rms_header",
        data=anc.encode_list(np.array(rms_header_lst, dtype="S20")),
        dtype="S20",
        compression="gzip",
    )

    f_h5.create_dataset(
        "grid/above_threshold", data=above_threshold, dtype=bool, compression="gzip"
    )

    f_h5.create_dataset(
        "grid/rms_rv", data=rms_rv, dtype=np.float64, compression="gzip"
    )

    return f_h5


# -------------------------------------

# ==============================================================================

# -------------------------------------


def simulation_and_ttv_analysis(
    i_grid,
    transit_id,
    perturber_id,
    n_bodies,
    perturber_par,
    t_start,
    t_end,
    err_T0,
    nT0_sel,
    names_lst,
    units,
    mr_type,
    nMC=1000,
    of_log=None,
):

    pid = os.getpid()
    tid = _thread.get_ident()
    results_sim = {}
    sim_num = i_grid + 1

    kel_header = "Mass_Msun Radius_Rsun Period_d sma_au ecc argp_deg meanA_deg inc_deg longN_deg".split()
    # anc.print_both('\n - sim %06d: START\n' %(sim_num), of_log)

    dof = 1

    try:
        # mass, radius, period, sma, ecc, argp, mean_an, inc, long_n, nt0_full, nrv_nmax = f90trades.wrapper_set_grid_orbelem(sim_num, perturber_id, perturber_grid, n_bodies)
        (
            mass,
            radius,
            period,
            sma,
            ecc,
            argp,
            mean_an,
            inc,
            long_n,
            nt0_full,
            nrv_nmax,
        ) = f90trades.wrapper_set_grid_orbelem(
            sim_num, perturber_id, perturber_par, n_bodies
        )
        # print ' | Set grid keplerian elements'
        kep_elem = np.column_stack(
            (mass, radius, period, sma, ecc, argp, mean_an, inc, long_n)
        )
    except:
        print(
            " | PID {0:d} || TID {1:d} sim {2:d}: error setting grip par".format(
                pid, tid, sim_num
            )
        )
        kep_elem = np.zeros((n_bodies, 9))

    results_sim["sim_num"] = sim_num
    results_sim["perturber_par_names"] = names_lst
    results_sim["perturber_par_units"] = units
    results_sim["perturber_par"] = perturber_par

    results_sim["keplerian_elements_header"] = kel_header
    results_sim["keplerian_elements"] = kep_elem

    # 2018-03-08 ADDING GLS AMPLITUDE AND PERIOD OF TTV/OC

    try:
        # anc.print_both(' - Run integration ...', of_log)
        (
            ttra_full,
            id_ttra_full,
            stats_ttra,
            time_rv_nmax,
            rv_nmax,
            stats_rv,
        ) = f90trades.wrapper_run_grid_combination(
            mass,
            radius,
            period,
            sma,
            ecc,
            argp,
            mean_an,
            inc,
            long_n,
            nt0_full,
            nrv_nmax,
        )
        # print ' | Completed orbit intergration'

        # anc.print_both(' - Select transit times:', of_log)
        (
            ttra_all,
            nttra_all,
            sel_pre,
            nttra_pre,
            sel_cheops,
            nttra_cheops,
        ) = select_and_sort_ttra(
            ttra_full,
            id_ttra_full,
            stats_ttra.astype(bool),
            id_transiting=transit_id,
            t_min_sel=t_start,
            t_max_sel=t_end,
            of_log=of_log,
        )
        # print ' | Selected and sorted TTs and RVs'
    except:
        ttra_all = None
        print(
            " | PID {0:d} || TID {1:d} sim {2:d}: error integrating".format(
                pid, tid, sim_num
            )
        )

    sel_kel = perturber_id - 1
    rms_rv = anc.compute_Kms(
        mass[0], mass[sel_kel] * cst.Msjup, inc[sel_kel], period[sel_kel], ecc[sel_kel]
    )

    if ttra_all is not None:
        T0, eT0 = dirty_transits(ttra_all, err_T0)
        T0_eT0 = np.column_stack((T0, eT0))
        print(
            " | PID {0:d} || TID {1:d} sim {2:d} shape(T0_eT0): {3:}".format(
                pid, tid, sim_num, np.shape(T0_eT0)
            )
        )
        # epo_pre, Tref_pre, Pref_pre, _    = anc.compute_lin_ephem(T0[sel_pre], eT0=eT0[sel_pre])
        Tin = T0[sel_pre][0]
        Pin = np.median(np.diff(T0[sel_pre]))
        epo_pre = np.rint( (T0[sel_pre] - Tin) / Pin)
        Tref_pre, Pref_pre, _ = pytrades.linear_fit(epo_pre, T0[sel_pref], ey=eT0[sel_pre])
        epo_pre = np.rint( (T0[sel_pre] - Tref[0]) / Pref[0])
        epo_cheops, oc_cheops, rms_cheops = compute_oc(
            Tref_pre[0], Pref_pre[0], T0[sel_cheops]
        )
        T0cheops = T0[sel_cheops]
        eT0cheops = eT0[sel_cheops]

        Ptra = period[transit_id - 1]
        print("periods: {}".format(period))
        print("transit id -1 = {} ==> P = {}".format(transit_id - 1, Ptra))
        if Pref_pre < Ptra:
            Pgls0 = 3.0 * Ptra
        else:
            Pgls0 = 3.0 * Pref_pre
        Pgls1 = 1.0 * (T0cheops[-1] - T0cheops[0])
        ogls = gls.Gls(
            (T0cheops, oc_cheops, eT0cheops),
            #  Pbeg=3.*Pref_pre, Pend=2.*(T0cheops[-1]-T0cheops[0]),
            Pbeg=Pgls0,
            Pend=Pgls1,
            ofac=5,
            verbose=False,
        )
        attv = ogls.hpstat["amp"] * 1440.0
        pttv = 1.0 / ogls.hpstat["fbest"]

        dof = np.shape(oc_cheops)[0] - 2
        if dof <= 0:
            dof = 1

        print(" | PID {0:d} || TID {1:d} sim {2:d} OK TT".format(pid, tid, sim_num))
        chi2r_cheops = do_chi2r(oc_cheops, eT0[sel_cheops], dof)

        rms_mc = [
            list(compute_rms_mc(oc_cheops, nT0_sel[i_n], nMC))
            for i_n in range(len(nT0_sel))
        ]
        chi2r_mc = [
            list(compute_chi2r_mc(oc_cheops, eT0[sel_cheops], nT0_sel[i_n], nMC))
            for i_n in range(len(nT0_sel))
        ]

        time_rv_all, rv_all = select_and_sort_rv(
            time_rv_nmax, rv_nmax, stats_rv.astype(bool)
        )

        if time_rv_all is None:
            time_rv = np.zeros((1, 2))
            # rms_rv = -1.

        else:
            time_rv = np.column_stack((time_rv_all, rv_all))
            # sel_kel = perturber_id-1
            # rms_rv = anc.compute_Kms(mass[0],mass[sel_kel]*cst.Msjup,inc[sel_kel],period[sel_kel],ecc[sel_kel])

    else:  # if the transit file is empty
        ttra_all = np.zeros((1)).reshape((-1, 1))
        sel_pre = np.zeros((1)).astype(bool)
        sel_cheops = np.zeros((1)).astype(bool)
        T0_eT0 = np.zeros((1, 2))
        Tref_pre = 0.0
        Pref_pre = 0.0
        attv = 0.0
        pttv = 0.0
        rms_cheops = -1.0
        rms_mc = -np.ones((len(nT0_sel), 2))
        chi2r_cheops = -1.0
        chi2r_mc = -np.ones((len(nT0_sel), 2))
        time_rv = np.zeros((1, 2))
        # rms_rv = anc.compute_Kms(mass[0],mDass[sel_kel]*cst.Msjup,inc[sel_kel],period[sel_kel],ecc[sel_kel])
        print(" | PID {0:d} || TID {1:d} sim {2:d} NO TT".format(pid, tid, sim_num))

    results_sim["ttra_all"] = ttra_all
    results_sim["sel_pre_cheops"] = sel_pre
    results_sim["sel_cheops"] = sel_cheops
    results_sim["T0_eT0"] = T0_eT0

    results_sim["A_TTV"] = attv
    results_sim["P_TTV"] = pttv

    results_sim["time_rv"] = time_rv
    results_sim["rms_rv"] = rms_rv

    results_sim["linear_ephem"] = [Tref_pre[0], Pref_pre[0]]
    results_sim["rms_ttv"] = rms_cheops
    results_sim["nT0_sel"] = nT0_sel
    results_sim["nMC"] = nMC
    results_sim[
        "rms_mc"
    ] = rms_mc  # mean and median OC rms for different nT0_sel during CHEOPS mission

    results_sim["chi2r_cheops"] = chi2r_cheops
    results_sim["chi2r_mc"] = chi2r_mc

    # anc.print_both('\n | sim %06d: COMPLETED\n' %(sim_num), of_log)

    del sim_num, names_lst, units, perturber_par, dof
    del mass, radius, period, sma, ecc, argp, mean_an, inc, long_n
    del ttra_all, sel_pre, sel_cheops, T0_eT0, time_rv, rms_rv, Tref_pre, Pref_pre
    del rms_cheops, nT0_sel, nMC, rms_mc, chi2r_cheops, chi2r_mc
    if "ogls" in locals():
        del ogls
    gc.collect()

    return results_sim


# -------------------------------------


def run_the_simulation(sim_args):

    i_grid = sim_args["i_grid"]
    transit_id = sim_args["transit_id"]
    perturber_id = sim_args["perturber_id"]
    perturber_par = sim_args["perturber_par"]
    n_bodies = sim_args["n_bodies"]
    t_start = sim_args["t_start"]
    t_end = sim_args["t_end"]
    err_T0 = sim_args["err_T0"]
    nT0_sel = sim_args["nT0_sel"]
    names_lst = sim_args["names_lst"]
    units = sim_args["units"]
    mr_type = sim_args["mr_type"]
    nMC = sim_args["nMC"]

    pid = os.getpid()
    tid = _thread.get_ident()
    print(" PID {0:d} || TID {1:d} RUNNING sim {2:d} ...".format(pid, tid, i_grid + 1))

    results_sim = simulation_and_ttv_analysis(
        i_grid,
        transit_id,
        perturber_id,
        n_bodies,
        perturber_par,
        t_start,
        t_end,
        err_T0,
        nT0_sel,
        names_lst,
        units,
        mr_type,
        nMC,
        of_log=None,
    )

    print(" PID {0:d} || TID {1:d} COMPLETED sim {2:d}".format(pid, tid, i_grid + 1))

    del i_grid, transit_id, perturber_id, perturber_par, n_bodies, t_start, t_end
    del err_T0, nT0_sel, names_lst, units, mr_type, nMC
    gc.collect()

    return results_sim


# -------------------------------------


def get_rms_from_results(
    results_sim, rms_grid, above_threshold, nT0_sel, err_T0, rms_rv
):

    i_grid = results_sim["sim_num"] - 1
    rms_ttv = results_sim["rms_ttv"]
    rms_mc = results_sim["rms_mc"]
    rms_grid[i_grid, 0] = rms_ttv
    rms_grid[i_grid, 1:] = np.reshape(rms_mc, newshape=(len(nT0_sel) * 2))
    above_threshold[i_grid, :] = rms_grid[i_grid, :] > err_T0
    rms_rv[i_grid] = results_sim["rms_rv"]

    return


# -------------------------------------

# ==============================================================================

# -------------------------------------


def compute_limits(xx, exx=None, delta=0.0):

    if exx is not None:
        xxmin = np.min(xx - exx)
        xxmax = np.max(xx + exx)
    else:
        xxmin = np.min(xx)
        xxmax = np.max(xx)
    dxx = np.abs(xxmax - xxmin)
    inflim = xxmin - dxx * delta
    suplim = xxmax + dxx * delta

    return inflim, suplim


# -------------------------------------


def insert_rms_lines(ax, rms_ref_m, rms_max, colors, lstyle):

    nref = np.shape(rms_ref_m)[0]

    for i_r in range(0, nref):
        if rms_ref_m[i_r] / 1440.0 <= rms_max:
            ax.axvline(
                rms_ref_m[i_r] / 1440.0,
                color=colors[i_r],
                ls=lstyle[i_r],
                lw=0.7,
                alpha=0.99,
                zorder=3,
                label="rms = {0:.0f} min".format(rms_ref_m[i_r]),
            )

    return ax


# -------------------------------------


def insert_grey_area(ax, min_x, max_x, z_order=8):

    ax.axvspan(min_x, max_x, color="gray", alpha=0.45, zorder=z_order)

    return ax


# -------------------------------------


def do_own_histogram(ax, vals, vrange, hist_color, z_order=5):

    kbins = anc.get_auto_bins(vals)
    cnts, bin_edges = np.histogram(vals, bins=kbins, range=vrange)
    # ncnts = cnts / np.shape(vals)[0]
    ncnts = cnts / np.max(cnts)
    bin_width = 0.9 * np.abs(bin_edges[1] - bin_edges[0])
    bin_pos = 0.5 * (bin_edges[:-1] + bin_edges[1:])
    ax.bar(bin_pos, ncnts, width=bin_width, color=hist_color, zorder=z_order)

    return ax, kbins


# -------------------------------------


def plot_rms(
    out_folder,
    body_id,
    grid_par,
    p_names,
    pl_lab,
    err_T0_d,
    rms,
    above,
    rms_head,
    rms_max_in=None,
    K_rv=None,
    rms_rv=None,
    of_log=None,
):

    ngrid, npar = np.shape(grid_par)
    sims = np.arange(0, ngrid).astype(int) + 1

    if rms_max_in is None:
        rmslim1, rmslim2 = compute_limits(rms, delta=0.05)
        # rmsmin = np.min(rms)
        rmsmin, rmslim1 = 0.0, 0.0
        rmsmax = np.max(rms)
    else:
        rmsmin = 0.0
        rmsmax = rms_max_in
        rmslim1, rmslim2 = compute_limits(np.array([rmsmin, rmsmax]), delta=0.05)
        rmslim1 = 0.0

    # rms_ref_m = np.array([1., 5., 10., 15., 20., 25., 30., 45., 60.])
    rms_ref_m = np.array([1.0, 5.0, 10.0, 15.0])
    nref = np.shape(rms_ref_m)[0]

    sel_rms = np.logical_and(rms > rmsmin, rms <= rmsmax)
    sel_rms_b = np.logical_and(rms > rmsmin, rms <= err_T0_d)
    if K_rv is not None and rms_rv is not None:
        if K_rv <= 0.0:
            anc.print_both("**WARNING: input K_rv <= 0", of_log)
        sel_rv = np.isnan(rms_rv) == False
        sel_rv = np.logical_and(rms_rv[sel_rv] > 0.0, rms_rv[sel_rv] <= K_rv)
        if np.sum(np.asarray(sel_rv, dtype=np.int)) == 0:
            anc.print_both(
                " There are no simulations with rms(RV) within 0.0 and {0:.2f} m/s.\n No plot to do.".format(
                    K_rv
                ),
                of_log,
            )
            return None
        sel_rms = np.logical_and(sel_rms, sel_rv)
        sel_rms_b = np.logical_and(sel_rms_b, sel_rv)
    sel_rms_a = np.logical_and(sel_rms, above)

    color_map = cm.viridis
    # color_map = cm.winter

    # ncmap = np.shape(color_map.colors)[0]
    # csel = np.random.choice(ncmap, size=nref, replace=False)
    # xvcolors = [color_map.colors[csel[ic]] for ic in range(0, nref)]

    xvcolors = ["black"] * nref
    lstyle = "-- -. : -".split()

    fig = plt.figure()
    fig.subplots_adjust(hspace=0.05, wspace=0.05)

    h_yes = 1  # set to 1 if the first plot is the histogram

    ax = plt.subplot2grid((npar + h_yes, 1), (0, 0))

    err_T0_min = err_T0_d * 1440.0
    rmsmax_min = rmsmax * 1440.0
    ncomb = np.shape(rms[sel_rms_a])[0]
    fig_title = r"$N_\textrm{{combinations}} = {0:d}$ for ${1:.1f} \leq \textrm{{rms(TTV)}} \leq {2:.1f}$ min".format(
        ncomb, err_T0_min, rmsmax_min
    )
    if K_rv is not None and rms_rv is not None:
        fig_title = r"{0:s} and $\textrm{{rms(TTV)}} \leq K_\textrm{{RV}} = {1:.2f}$ m/s".format(
            fig_title, K_rv
        )

    # fig.text(0.5, 0.85, fig_title, fontsize=11)
    ax.set_title(fig_title)

    ax.ticklabel_format(useOffset=False)

    # ax.axvline(0., color='darkgray', ls='--', lw=0.5, alpha=0.33, zorder=1)
    ax.axvline(err_T0_d, color="C0", ls="-", lw=1.7, alpha=0.99, zorder=9)
    ax = insert_rms_lines(ax, rms_ref_m, rmsmax, xvcolors, lstyle)

    # kbins = get_bins(rms[sel_rms], rule='sq')
    ax, kbins = do_own_histogram(ax, rms[sel_rms], [rmsmin, rmsmax], "lightgray")
    anc.print_both(
        " | plotting {0:s} histogram witn nbins (all) = {1:d}".format(rms_head, kbins),
        of_log,
    )
    # ax.hist(rms[sel_rms], bins=kbins, range=[rmsmin, rmsmax], histtype='stepfilled', color='gray', edgecolor='lightgray', align='mid', orientation='vertical') # plot all in gray

    # kbins = get_bins(rms[sel_rms_a], rule='sq')
    # anc.print_both(' | plotting %s histogram witn nbins (above) = %d' %(rms_head, kbins), of_log)
    # ax.hist(rms[sel_rms_a], bins=kbins, range=[err_T0_d, rmsmax], histtype='stepfilled', color='C0', edgecolor='lightgray', align='mid', orientation='vertical') # plot above in azure
    ax, kbins = do_own_histogram(ax, rms[sel_rms_a], [err_T0_d, rmsmax], "C0")
    anc.print_both(
        " | plotting {0:s} histogram witn nbins (above) = {1:d}".format(
            rms_head, kbins
        ),
        of_log,
    )

    ax = insert_grey_area(ax, rmslim1, err_T0_d)

    ax.set_xlim([rmslim1, rmslim2])
    ax.get_xaxis().set_visible(False)

    ax.get_yaxis().set_visible(True)
    # ylabel = 'Normalised Counts'
    ylabel = r"$N / N_\textrm{{max}}$"
    ax.set_ylabel(ylabel)
    ax.get_yaxis().set_label_coords(-0.075, 0.5)

    plt.draw()

    for iy in range(0, npar):
        yy = grid_par[:, iy]
        # yya = yy[above]
        # yylim1, yylim2 = compute_limits(yy, delta=0.5)

        anc.print_both(
            " | plotting %s vs %s%d" % (rms_head, p_names[iy], body_id), of_log
        )

        ax = plt.subplot2grid((npar + h_yes, 1), (iy + h_yes, 0))
        ax.ticklabel_format(useOffset=False)

        ax.axvline(
            err_T0_d,
            color="C0",
            ls="-",
            lw=1.7,
            alpha=0.99,
            zorder=9,
            label=r"rms = $\sigma_{{\textrm{{T}}_0}}$ = {0:.3f} min".format(err_T0_min),
        )
        ax = insert_rms_lines(ax, rms_ref_m, rms_max_in, xvcolors, lstyle)

        ax.scatter(
            rms[sel_rms_a],
            yy[sel_rms_a],
            marker="o",
            s=25,
            c="None",
            edgecolors="C0",
            lw=1.1,
            alpha=1.0,
            zorder=8,
            label="observable",
        )
        ax.scatter(
            rms[sel_rms_a],
            yy[sel_rms_a],
            c=sims[sel_rms_a],
            cmap=color_map,
            marker="o",
            s=24,
            lw=0,
            alpha=1.0,
            zorder=8,
        )

        # ax.scatter(rms, yy, color='black', marker='o', s=20, lw=0, alpha=0.05, zorder=5)
        # sax = ax.scatter(rms[sel_rms], yy[sel_rms], c=sims, cmap=color_map, marker='o', s=25, lw=0, alpha=0.25, zorder=7)
        ax.scatter(
            rms[sel_rms_b],
            yy[sel_rms_b],
            c=sims[sel_rms_b],
            cmap=color_map,
            marker="o",
            s=25,
            lw=0,
            alpha=0.25,
            zorder=8,
        )
        # ax.scatter(rms[above], yya, c=sims[above], cmap=color_map, marker='o', s=20, lw=0, alpha=0.1, zorder=8)

        ax = insert_grey_area(ax, rmslim1, err_T0_d)

        ax.set_xlim([rmslim1, rmslim2])
        ax.get_xaxis().set_visible(False)
        if iy == npar - 1:
            ax.get_xaxis().set_visible(True)
            xlabel = "rms(O-C) [days]"
            ax.set_xlabel(xlabel)

        ax.get_yaxis().set_visible(True)
        ylabel = pl_lab[iy]
        ax.set_ylabel(ylabel)
        # ax.yaxis.labelpad = 20
        ax.get_yaxis().set_label_coords(-0.075, 0.5)

        if iy == 0:
            ax.legend(
                loc="upper left", fontsize=9, bbox_to_anchor=(1, 1), ncol=1
            )  # one to legend them all

        plt.draw()

    # create ad-hoc colorbar
    ngrid_full = np.shape(rms)[0]
    ticks_lst = list(np.linspace(1, ngrid_full, 5, endpoint=True).astype(int))
    fig.subplots_adjust(hspace=0.05)
    cbaxes = fig.add_axes([0.15, 0.04, 0.7, 0.02])  # [left, bottom, width, height]
    nbar = mpl.colors.Normalize(vmin=1, vmax=ngrid_full)
    cbar = mpl.colorbar.ColorbarBase(
        cbaxes, cmap=color_map, norm=nbar, orientation="horizontal", ticks=ticks_lst
    )
    cbar.set_label("simulation number")

    plt.draw()

    fig_file = os.path.join(out_folder, "plot_%s_max%.1fmin" % (rms_head, rmsmax_min))
    if K_rv is not None and rms_rv is not None:
        fig_file = "%s_K%.2fmps" % (fig_file, K_rv)
    fig.savefig("%s.pdf" % (fig_file), bbox_inches="tight", dpi=300)
    plt.close(fig)
    anc.print_both(" ++ saved plot %s.pdf" % (fig_file), of_log)

    return fig_file


# -------------------------------------


def plot_summary(
    out_folder,
    hdf5_file,
    body_id=3,
    mr_type="e",
    rms_max_in=None,
    K_rv=None,
    of_log=None,
):

    letters = "a b c d e f g h i j k l m n o p q r s t u v x y z".split()
    p_names = "M R P a e w mA tau i lN".split()
    npar = np.shape(p_names)[0]

    l_names = r"M R P a e \omega \nu \tau i \Omega".split()
    units = r"[$M_\oplus$] [$R_\oplus$] [day] [au] [] [$^\circ$] [$^\circ$] [BJD] [$^\circ$] [$^\circ$]".split()
    units[2] = ""
    # if(mr_type in 'j ju jup jupi jupit jupit jupite jupiter'):
    jupl = "j ju jup jupi jupit jupit jupite jupiter".split()
    nepl = "n ne nep nept neptu neptun neptune".split()
    if np.any([mr_type == jj for jj in jupl]):
        units[0] = r"[$M_\textrm{{Jup}}$]"
        units[1] = r"[$R_\textrm{{Jup}}$]"
    # elif(mr_type in 'n ne nep nept neptu neptun neptune'):
    elif np.any([mr_type == nn for nn in nepl]):
        units[0] = r"[$M_\textrm{{Nep}}$]"
        units[1] = r"[$R_\textrm{{Nep}}$]"
    else:
        units[0] = r"[$M_\oplus$]"
        units[1] = r"[$R_\oplus$]"

    pl_lab = [
        r"${0:s}_\textrm{{1:s}}$ {2:s}".format(
            l_names[i_n], letters[body_id - 1], units[i_n].replace("[]", " ").strip()
        )
        for i_n in range(0, npar)
    ]

    f_h5 = h5py.File(hdf5_file, "r")
    grid_par = f_h5["grid/perturber_grid"][...]
    err_T0_d = f_h5["grid"].attrs["err_T0_day"]

    rms_grid = f_h5["grid/rms_grid"][...]
    rms_header = f_h5["grid/rms_header"][...]
    above_threshold = f_h5["grid/above_threshold"][...]
    try:
        rms_rv = f_h5["grid/rms_rv"][...]
    except:
        rms_rv = None
    f_h5.close()

    plot_folder = os.path.join(out_folder, "plots")
    if not os.path.isdir(plot_folder):
        os.makedirs(plot_folder)

    # ngrid, nrms = np.shape(rms_grid)
    _, nrms = np.shape(rms_grid)

    id_noRtau = np.array([0, 2, 3, 4, 5, 6, 8, 9]).astype(int)

    fig_files = []
    for irms in range(0, nrms):
        fig_file = plot_rms(
            plot_folder,
            body_id,
            grid_par[:, id_noRtau],
            np.array(p_names)[id_noRtau],
            np.array(pl_lab)[id_noRtau],
            err_T0_d,
            rms_grid[:, irms],
            above_threshold[:, irms],
            rms_header[irms],
            rms_max_in,
            K_rv,
            rms_rv,
            of_log,
        )  # here call function that plot single rms vs parameters
        fig_files.append(fig_file)

    return fig_files


# ==============================================================================


def compute_grid_to_freq_obs(
    grid_par,
    rms_all,
    lower_rms=0.0,
    upper_rms=15.0 / 1440.0,
    K_rv=None,
    rms_rv=None,
    of_log=None,
):

    # ngrid, npar = np.shape(grid_par)
    ngrid, _ = np.shape(grid_par)
    #            0 1 2 3 4 5 6  7   8 9
    # p_names = 'M R P a e w mA tau i lN'.split() # columns of grid_par (npar)

    m_grid = np.unique(grid_par[:, 0])
    n_m = np.shape(m_grid)[0]

    p_grid = np.unique(grid_par[:, 2])
    n_p = np.shape(p_grid)[0]

    # n_mp = n_m * n_p

    n_e = np.shape(np.unique(grid_par[:, 4]))[0]
    n_mA = np.shape(np.unique(grid_par[:, 6]))[0]
    n_i = np.shape(np.unique(grid_par[:, 8]))[0]

    n_emAi = n_e * n_mA * n_i
    # n_emAi = ngrid / n_mp

    nrms = np.shape(rms_all)[1]

    mass, period = [], []
    freq_obs, freq_rv = [], []
    for iim in range(0, n_m):
        for iip in range(0, n_p):

            mm = m_grid[iim]
            pp = p_grid[iip]
            mass.append(mm)
            period.append(pp)

            freq_otemp = np.zeros((nrms)).astype(int)
            if K_rv is not None and rms_rv is not None:
                freq_rvtemp = np.zeros((1)).astype(int)
            for iig in range(0, ngrid):
                if mm == grid_par[iig, 0] and pp == grid_par[iig, 2]:
                    rms_sel = np.logical_and(
                        rms_all[iig, :] >= lower_rms, rms_all[iig, :] <= upper_rms
                    )
                    freq_otemp += rms_sel.astype(int)
                    if K_rv is not None and rms_rv is not None:
                        rv_sel = np.logical_and(rms_rv[iig] >= 0.0, rms_rv[iig] <= K_rv)
                        freq_rvtemp += rv_sel.astype(int)

            freq_obs.append(np.float64(freq_otemp) / np.float64(n_emAi))
            if K_rv is not None and rms_rv is not None:
                freq_rv.append(np.float64(freq_rvtemp) / np.float64(n_emAi))

    mass = np.array(mass)
    period = np.array(period)
    freq_obs = np.array(freq_obs)

    if K_rv is None or rms_rv is None:
        freq_rv = None
    else:
        freq_rv = np.array(freq_rv)

    # if (K_rv is not None and rms_rv is not None):
    # sel_rv = np.logical_and(rms_rv >= 0., rms_rv <= K_rv)
    # else:
    # sel_rv = np.ones((ngrid)).astype(bool)

    return mass, period, freq_obs, freq_rv


# ==============================================================================


def plot_m_vs_p_vs_freq_obs(
    plot_file, mass, period, freq_obs, freq_rv, n_label="all 3 5 6 10 15".split()
):

    # nplt = 6 # plot rms mean or median
    nplt = len(n_label)

    fig = plt.figure()
    fig.subplots_adjust(hspace=0.25, wspace=0.065)

    # xmin, xmax = anc.compute_limits(period, delta=0.025)
    log_p = np.log10(period)
    # xmin, xmax = anc.compute_limits(log_p, delta=0.025)

    # ymin, ymax = anc.compute_limits(mass, delta=0.025)
    log_m = np.log10(mass)
    # ymin, ymax = anc.compute_limits(log_m, delta=0.025)

    pp = np.unique(log_p)
    n_p = np.shape(pp)[0]
    mm = np.unique(log_m)
    n_m = np.shape(mm)[0]

    # ntks = 5
    # lp_ticks = np.linspace(log_p[0], log_p[-1], num=ntks, endpoint=True)

    # if(freq_rv is not None):
    #   gray_rv = freq_rv <= 0.0

    clev = [0.0 + ilev * 0.1 for ilev in range(0, 11)]

    ix = 0
    iy = 0
    for iplt in range(0, nplt):

        iy = iplt % 2
        ax = plt.subplot2grid((nplt, 2), (ix, iy))
        print(" Plot n. {:d}: ".format(iplt))

        # ax.scatter(period, mass, c=freq_obs[:,iplt], cmap=cm.viridis, marker='s', s=20, zorder=5, label='T0 CHEOPS N = %s' %(n_label[iplt]))
        # ax.legend(loc='upper right', fontsize=9, numpoints=0)

        # ssax = ax.scatter(log_p, log_m, c=freq_obs[:,iplt], cmap=cm.magma, marker='s', s=30, alpha=0.5, zorder=10)

        # ax.text(0.90, 0.90, 'T0 CHEOPS N = %s' %(n_label[iplt]), fontsize=8, horizontalalignment='right', verticalalignment='center',transform=ax.transAxes)

        # print(' | scatter done')

        mesh_p, mesh_m = np.meshgrid(pp, mm, indexing="ij")
        mesh_f = np.zeros((n_p, n_m))
        mesh_f_im = np.zeros((n_p, n_m))
        ipm = 0
        for i_p in range(0, n_p):
            for i_m in range(0, n_m):
                mesh_f[i_p, i_m] = freq_obs[ipm, iplt]
                mesh_f_im[i_m, i_p] = freq_obs[ipm, iplt]
                ipm += 1

        # scax = ax.imshow(np.flipud(mesh_f_im.T), cmap=cm.viridis, extent=[np.min(log_p), np.max(log_p), np.min(log_m), np.max(log_m)], vmin=0.0, vmax=1.0, interpolation='bilinear', aspect='auto', zorder=4)
        # print(' | mesh/imshow done')

        scax = ax.contourf(
            mesh_p, mesh_m, mesh_f_im, levels=clev, cmap=cm.viridis, zorder=4
        )
        print(" | mesh/contourf done")

        # colorbar
        cbar = plt.colorbar(scax, ax=ax, orientation="vertical", pad=0.01)
        cbar.ax.set_ylabel(r"$f(\textrm{{TTV}}) >$ threshold")
        cbar.ax.tick_params(labelsize=10)

        # ccax = ax.contour(mesh_p, mesh_m, mesh_f_im, levels=clev, colors='k', linewidths=0.5, zorder=5)
        # plt.clabel(ccax, clev[1::2], fontsize=6, inline=True, zorder=10) # not writing labels with plot of imshow!
        # print(' | mesh/contour done')

        if freq_rv is not None:
            # scax_rv = ax.scatter(log_p[gray_rv], log_m[gray_rv], c='darkgray', marker='s', s=20, alpha=0.66, zorder=5)
            mesh_rv = np.zeros((n_p, n_m))
            ipm = 0
            for i_p in range(0, n_p):
                for i_m in range(0, n_m):
                    if freq_rv[ipm] > 0.0:
                        mesh_rv[i_p, i_m] = np.ma.masked
                    else:
                        mesh_rv[i_p, i_m] = 1.0
                    ipm += 1
            ax.contourf(
                mesh_p,
                mesh_m,
                mesh_rv.T,
                colors="gray",
                # cmap=cm.gray,
                ls="",
                alpha=0.90,
                zorder=6,
            )
            # scax_rv = ax.scatter(log_p[gray_rv], log_m[gray_rv], c='darkgray', marker='s', s=20, alpha=0.88, zorder=5)
            print(" | excluded by RV done")

        ax.axhline(
            np.log10(1.0), color="lightgray", ls="-", lw=0.8, alpha=0.66, zorder=6
        )
        ax.axhline(
            np.log10(5.0), color="lightgray", ls="-", lw=0.8, alpha=0.66, zorder=6
        )
        ax.axhline(
            np.log10(10.0), color="lightgray", ls="-", lw=0.8, alpha=0.66, zorder=6
        )

        # change ticklabels to base
        # ax.xaxis.set_ticks([10.**lp_ticks[i] for i in range(0, ntks)])
        # ax.set_xticks([10.**lp_ticks[i] for i in range(0, ntks)])

        # ax.xaxis.set_major_formatter(ScalarFormatter())
        # ax.get_xaxis().set_major_formatter(ScalarFormatter())

        ax.xaxis.set_tick_params(labelsize=10)
        ax.set_xlabel(
            r"$\log P [\textrm{{days}}] , N_{{T0}} = {:s}$".format(n_label[iplt])
        )
        # ax.set_xlim([np.log10(xmin), np.log10(xmax)])

        # ax.yaxis.set_ticks([mass[i] for i in range(0, np.shape(period)[0])])
        # ax.yaxis.set_major_formatter(ScalarFormatter())

        ax.yaxis.set_tick_params(labelsize=10)
        ax.set_ylabel(r"$\log M [M_\oplus]$")
        # ax.set_ylim([np.log10(ymin), np.log10(ymax)])

        ax.tick_params(direction="inout")
        ix += iy
        plt.draw()

    plt.tight_layout()
    fig.savefig(plot_file)  # , bbox_inches='tight')
    plt.close(fig)

    return


# ==============================================================================


def read_and_plot_p_vs_m_vs_freq(
    out_folder, hdf5_file, rms_min=None, rms_max=15.0 / 1440.0, K_rv=None, of_log=None
):

    f_h5 = h5py.File(hdf5_file, "r")
    grid_par = f_h5["grid/perturber_grid"][...]
    err_T0_d = f_h5["grid"].attrs["err_T0_day"]
    rms_grid = f_h5["grid/rms_grid"][...]
    # rms_header = f_h5['grid/rms_header'][...]
    # above_threshold = f_h5['grid/above_threshold'][...]
    try:
        rms_rv = f_h5["grid/rms_rv"][...]
    except:
        rms_rv = None
    f_h5.close()

    if rms_min is not None:
        rms_lower = rms_min
    else:
        rms_lower = err_T0_d

    mass, period, freq_obs, freq_rv = compute_grid_to_freq_obs(
        grid_par,
        rms_grid,
        lower_rms=rms_lower,
        upper_rms=rms_max,
        K_rv=K_rv,
        rms_rv=rms_rv,
        of_log=of_log,
    )

    _, nrms = np.shape(rms_grid)

    mean_sel = np.ones((nrms)).astype(bool)
    mean_sel[1:] = np.arange(1, nrms, 1) % 2 == 0

    median_sel = np.ones((nrms)).astype(bool)
    median_sel[1:] = np.arange(1, nrms, 1) % 2 == 1

    plt_file_mean = os.path.join(
        out_folder, "p_vs_m_rms_{0:s}_max{1:.2f}min".format("mean", rms_max * 1440.0)
    )
    plt_file_med = os.path.join(
        out_folder, "p_vs_m_rms_{0:s}_max{1:.2f}min".format("median", rms_max * 1440.0)
    )
    if K_rv is not None and rms_rv is not None:
        plt_file_mean = "{0:s}_K{1:.2f}mps".format(plt_file_mean, K_rv)
        plt_file_med = "{0:s}_K{1:.2f}mps".format(plt_file_med, K_rv)
    plt_file_mean = "{0:s}.pdf".format(plt_file_mean)
    plt_file_med = "{0:s}.pdf".format(plt_file_med)

    # call function that plots P vs M with color scale from freq_obs[irms] and sel by K_rv, where irms is the sing rms type (TTV, Nx, Ny, ...Nz)
    plot_m_vs_p_vs_freq_obs(plt_file_mean, mass, period, freq_obs[:, mean_sel], freq_rv)
    print(plt_file_mean)
    plot_m_vs_p_vs_freq_obs(
        plt_file_med, mass, period, freq_obs[:, median_sel], freq_rv
    )
    print(plt_file_med)

    return plt_file_mean, plt_file_med


# ==============================================================================
# ==============================================================================


def get_chi2r_rms_from_results(
    results_sim, chi2r_oc_grid, rms_oc_grid, above_threshold, nT0_sel, err_T0, rms_rv
):

    i_grid = results_sim["sim_num"] - 1
    chi2r_cheops = results_sim["chi2r_cheops"]
    chi2r_mc = results_sim["chi2r_mc"]
    chi2r_oc_grid[i_grid, 0] = chi2r_cheops
    chi2r_oc_grid[i_grid, 1:] = np.reshape(chi2r_mc, newshape=(len(nT0_sel) * 2))

    rms_ttv = results_sim["rms_ttv"]
    rms_mc = results_sim["rms_mc"]
    rms_oc_grid[i_grid, 0] = rms_ttv
    rms_oc_grid[i_grid, 1:] = np.reshape(rms_mc, newshape=(len(nT0_sel) * 2))
    above_threshold[i_grid, :] = rms_oc_grid[i_grid, :] > err_T0

    rms_rv[i_grid] = results_sim["rms_rv"]

    return


# -------------------------------------


def save_chi2r_rms_analysis(
    f_h5,
    chi2r_oc_grid,
    rms_oc_grid,
    above_threshold,
    nT0_sel,
    rms_rv,
    attv_grid=None,
    pttv_grid=None,
):

    # save chi2r_oc_grid
    chi2r_header_str = "chi2r_oc {0:s}".format(
        " ".join(
            [
                "chi2r_mean_N{0:d}_d chi2r_median_N{1:d}_d".format(
                    nT0_sel[i_n], nT0_sel[i_n]
                )
                for i_n in range(len(nT0_sel))
            ]
        )
    )
    chi2r_header_lst = chi2r_header_str.split()

    f_h5.create_dataset(
        "grid/chi2r_oc_grid", data=chi2r_oc_grid, dtype=np.float64, compression="gzip"
    )
    f_h5.create_dataset(
        "grid/chi2r_header",
        data=np.array(chi2r_header_lst, dtype="S20"),
        dtype="S20",
        compression="gzip",
    )

    # save rms_oc_grid
    rms_header_str = "rms_ttv_d {0:s}".format(
        " ".join(
            [
                "rms_mean_N{0:d}_d rms_median_N{1:d}_d".format(
                    nT0_sel[i_n], nT0_sel[i_n]
                )
                for i_n in range(len(nT0_sel))
            ]
        )
    )
    rms_header_lst = rms_header_str.split()

    f_h5.create_dataset(
        "grid/rms_oc_grid", data=rms_oc_grid, dtype=np.float64, compression="gzip"
    )
    f_h5.create_dataset(
        "grid/rms_header",
        data=np.array(rms_header_lst, dtype="S20"),
        dtype="S20",
        compression="gzip",
    )
    f_h5.create_dataset(
        "grid/above_threshold", data=above_threshold, dtype=bool, compression="gzip"
    )

    # save rms_RV
    f_h5.create_dataset(
        "grid/rms_rv", data=rms_rv, dtype=np.float64, compression="gzip"
    )

    if not attv_grid is None:
        f_h5.create_dataset(
            "grid/A_TTV", data=attv_grid, dtype=np.float64, compression="gzip"
        )
        f_h5.create_dataset(
            "grid/P_TTV", data=pttv_grid, dtype=np.float64, compression="gzip"
        )

    return f_h5


# -------------------------------------


def compute_grid_to_chi2r_rms(
    grid_par, chi2r_all, rms_oc_all, rms_rv=None, of_log=None
):

    ngrid, _ = np.shape(grid_par)
    #            0 1 2 3 4 5 6  7   8 9
    # p_names = 'M R P a e w mA tau i lN'.split() # columns of grid_par (npar)

    m_grid = np.unique(grid_par[:, 0])
    n_m = np.shape(m_grid)[0]

    p_grid = np.unique(grid_par[:, 2])
    n_p = np.shape(p_grid)[0]

    # n_mp = n_m * n_p

    # n_e  = np.shape(np.unique(grid_par[:,4]))[0]
    # n_mA = np.shape(np.unique(grid_par[:,6]))[0]
    # n_i  = np.shape(np.unique(grid_par[:,8]))[0]

    # n_emAi = n_e * n_mA * n_i

    # nrms = np.shape(rms_oc_all)[1] # columns mean-median

    mass, period = [], []
    chi2r_MP, rms_oc_MP, rv_MP = [], [], []
    eps64bit = np.finfo(1.0).eps

    for iim in range(0, n_m):
        for iip in range(0, n_p):

            mm = m_grid[iim]
            pp = p_grid[iip]
            mass.append(mm)
            period.append(pp)

            chi2r_temp, rms_oc_temp, rv_temp = [], [], []

            for iig in range(0, ngrid):
                adm = np.abs(grid_par[iig, 0] - mm)
                adp = np.abs(grid_par[iig, 2] - pp)
                if adm <= eps64bit and adp <= eps64bit:

                    chi2r_temp.append(chi2r_all[iig, :])
                    rms_oc_temp.append(rms_oc_all[iig, :])
                    if rms_rv is not None:
                        # rv_temp.append(rms_rv[iig, :])
                        rv_temp.append(rms_rv[iig])

            # chi2r_MP.append(np.mean(chi2r_temp, axis=0))
            chi2r_MP.append(
                np.percentile(
                    np.array(chi2r_temp), 50.0, interpolation="midpoint", axis=0
                )
            )

            # rms_oc_MP.append(np.mean(rms_oc_temp, axis=0))
            rms_oc_MP.append(
                np.percentile(
                    np.array(rms_oc_temp), 50.0, interpolation="midpoint", axis=0
                )
            )

            if rms_rv is not None:
                # rv_MP.append(np.mean(rv_temp, axis=0))
                rv_MP.append(
                    np.percentile(
                        np.array(rv_temp), 50.0, interpolation="midpoint", axis=0
                    )
                )

    mass = np.array(mass)
    period = np.array(period)

    chi2r_MP = np.array(chi2r_MP)
    rms_oc_MP = np.array(rms_oc_MP)

    if rms_rv is None:
        rv_MP = None
    else:
        rv_MP = np.array(rv_MP)

    return mass, period, chi2r_MP, rms_oc_MP, rv_MP


# -------------------------------------


def compute_grid_to_chi2r_rms_MP_e0copmA90(
    grid_par, chi2r_all, rms_oc_all, rms_rv=None, of_log=None
):

    ngrid, _ = np.shape(grid_par)
    #            0 1 2 3 4 5 6  7   8 9
    # p_names = 'M R P a e w mA tau i lN'.split() # columns of grid_par (npar)

    m_grid = np.unique(grid_par[:, 0])
    n_m = np.shape(m_grid)[0]

    p_grid = np.unique(grid_par[:, 2])
    n_p = np.shape(p_grid)[0]

    # n_mp = n_m * n_p

    # n_e  = np.shape(np.unique(grid_par[:,4]))[0]
    # n_mA = np.shape(np.unique(grid_par[:,6]))[0]
    # n_i  = np.shape(np.unique(grid_par[:,8]))[0]

    # n_emAi = n_e * n_mA * n_i

    # nrms = np.shape(rms_oc_all)[1] # columns mean-median

    mass, period = [], []
    chi2r_MP, rms_oc_MP, rv_MP = [], [], []
    eps64bit = np.finfo(1.0).eps

    for iim in range(0, n_m):
        for iip in range(0, n_p):

            mm = m_grid[iim]
            pp = p_grid[iip]
            mass.append(mm)
            period.append(pp)

            chi2r_temp, rms_oc_temp, rv_temp = [], [], []

            for iig in range(0, ngrid):
                adm = np.abs(grid_par[iig, 0] - mm)
                adp = np.abs(grid_par[iig, 2] - pp)
                admA = np.abs(grid_par[iig, 6] - 90.0)

                if adm <= eps64bit and adp <= eps64bit and admA <= eps64bit:
                    chi2r_temp.append(chi2r_all[iig, :])
                    rms_oc_temp.append(rms_oc_all[iig, :])
                    if rms_rv is not None:
                        # rv_temp.append(rms_rv[iig, :])
                        rv_temp.append(rms_rv[iig])

            # chi2r_MP.append(np.mean(chi2r_temp, axis=0))
            # chi2r_MP.append(np.percentile(np.array(chi2r_temp), 50., interpolation='midpoint', axis=0))
            chi2r_MP.append(chi2r_temp)

            # rms_oc_MP.append(np.mean(rms_oc_temp, axis=0))
            # rms_oc_MP.append(np.percentile(np.array(rms_oc_temp), 50., interpolation='midpoint', axis=0))
            rms_oc_MP.append(rms_oc_temp)

            if rms_rv is not None:
                # rv_MP.append(np.mean(rv_temp, axis=0))
                # rv_MP.append(np.percentile(np.array(rv_temp), 50., interpolation='midpoint', axis=0))
                rv_MP.append(rv_temp)

    mass = np.array(mass)
    period = np.array(period)

    # chi2r_MP = np.array(chi2r_MP)
    # testing reshape
    n1, n2, n3 = np.shape(chi2r_MP)
    chi2r_MP = np.array(chi2r_MP).reshape((n1 * n2, n3))
    # print ' chi2r n1 x n2 x n3 = ', n1, n2 , n3, ' ==> new shape: ',np.shape(chi2r_MP)
    # rms_oc_MP = np.array(rms_oc_MP)
    # testing reshape
    n1, n2, n3 = np.shape(rms_oc_MP)
    rms_oc_MP = np.array(rms_oc_MP).reshape((n1 * n2, n3))
    # print ' rmsOC n1 x n2 x n3 = ', n1, n2 , n3, ' ==> new shape: ',np.shape(rms_oc_MP)

    if rms_rv is None:
        rv_MP = None
    else:
        rv_MP = np.array(rv_MP)

    return mass, period, chi2r_MP, rms_oc_MP, rv_MP


# ==============================================================================


def ticklabels_logtolin(tickslabels_obj):

    ticks_log_str = [
        tickslabels_obj[ii].get_text() for ii in range(len(tickslabels_obj))
    ]
    nticks = len(ticks_log_str)
    ticks_lin_str = [""] * nticks
    for ii in range(nticks):
        if len(ticks_log_str[ii]) > 0:
            tick_temp = np.power(10.0, np.float64(ticks_log_str[ii].strip("$")))
            ticks_lin_str[ii] = r"${.2f}$".format(tick_temp)

    return ticks_lin_str


# ==============================================================================


def plot_MMR(ax, P_tra):

    MMR_1st = [2.0, 3.0 / 2.0, 4.0 / 3.0, 5.0 / 4.0]
    MMR_2nd = [3.0, 5.0 / 3.0]
    MMR_3rd = [4.0, 5.0 / 2.0]

    # MMR j:j-1 // j:j+1
    for ii in range(0, len(MMR_1st)):
        Pin = P_tra / MMR_1st[ii]
        ax.axvline(Pin, color="red", ls="-", lw=1.2, alpha=0.66, zorder=7)
        # lgPin = np.log10(Pin)
        # ax.axvline(lgPin, color='red', ls='-', lw=1.2, alpha=0.66, zorder=7)
        Pout = P_tra * MMR_1st[ii]
        ax.axvline(Pout, color="red", ls="-", lw=1.2, alpha=0.66, zorder=7)
        # lgPout = np.log10(Pout)
        # ax.axvline(lgPout, color='red', ls='-', lw=1.2, alpha=0.66, zorder=7)

    # MMR j:j-2 // j:j+2
    for ii in range(0, len(MMR_2nd)):
        Pin = P_tra / MMR_2nd[ii]
        ax.axvline(Pin, color="orange", ls="-", lw=0.95, alpha=0.66, zorder=7)
        # lgPin = np.log10(Pin)
        # ax.axvline(lgPin, color='orange', ls='-', lw=0.95, alpha=0.66, zorder=7)

        Pout = P_tra * MMR_2nd[ii]
        ax.axvline(Pout, color="orange", ls="-", lw=0.95, alpha=0.66, zorder=7)
        # lgPout = np.log10(P_tra*MMR_2nd[ii])
        # ax.axvline(lgPout, color='orange', ls='-', lw=0.95, alpha=0.66, zorder=7)

    # MMR j:j-3 // j:j+3
    for ii in range(0, len(MMR_3rd)):
        Pin = P_tra / MMR_3rd[ii]
        ax.axvline(Pin, color="cyan", ls="-", lw=0.8, alpha=0.66, zorder=7)
        # lgPin = np.log10(P_tra/MMR_3rd[ii])
        # ax.axvline(lgPin, color='cyan', ls='-', lw=0.8, alpha=0.66, zorder=7)
        Pout = P_tra * MMR_3rd[ii]
        ax.axvline(Pout, color="cyan", ls="-", lw=0.8, alpha=0.66, zorder=7)
        # lgPout = np.log10(P_tra*MMR_3rd[ii])
        # ax.axvline(lgPout, color='cyan', ls='-', lw=0.8, alpha=0.66, zorder=7)

    return ax


# ==============================================================================

# THE GOOD PLOT ONE!!
# THE GOOD PLOT ONE!!
# THE GOOD PLOT ONE!!
# THE GOOD PLOT ONE!!


def plot_m_vs_p_vs_chi2r(
    plot_file,
    mass,
    period,
    chi2r_MP,
    rms_oc_MP,
    rms_oc_max=5.0 / 1440.0,
    K_rv=5.0,
    rv_MP=None,
    P_tra=None,
    n_label="all 3 5 6 10 15".split(),
):

    print("*** plot_m_vs_p_vs_chi2r START ***")

    n_p = np.shape(np.unique(period))[0]
    n_m = np.shape(np.unique(mass))[0]

    mesh_p = period.reshape((n_m, n_p))
    if P_tra is not None:
        mesh_p = mesh_p / P_tra
    mesh_m = mass.reshape((n_m, n_p))

    # contour levels
    cmap_get = cm.get_cmap("viridis")

    clev05 = [0.0, 1.0, 2.0, 3.0, 4.0, 5.0]
    maxchi2r = np.max(np.max(chi2r_MP))
    if maxchi2r <= 6.0:
        clev = clev05 + [6.0]
    else:
        lval6max = np.linspace(np.log10(6.0), np.log10(maxchi2r), num=4, endpoint=True)
        clev = clev05 + [np.power(10.0, lval6max[ic]) for ic in range(4)]
    nlev = np.shape(clev)[0]
    levcols = cmap_get(np.linspace(0.0, 1.0, num=nlev, endpoint=True))

    lev_minutes = [1.0, 5.0, 10.0]
    rms_lev = [ilev / 1400.0 for ilev in lev_minutes]

    origin_mesh = None

    nplt = len(n_label)
    print("*** nplt = {}".format(nplt))
    nrow = int(nplt / 2)
    ncol = 2
    # fig = plt.figure()
    # fig.subplots_adjust(hspace=1.0, wspace=1.0)
    fig, axs = plt.subplots(nrow, ncol)
    fig.subplots_adjust(hspace=0.45, wspace=0.55)

    ix = 0
    iy = 0
    for iplt in range(0, nplt):

        iy = iplt % 2
        # ax = plt.subplot2grid((nplt, 2), (ix,iy))
        ax = axs[ix, iy]
        print((" Plot n. {:d}: ".format(iplt)))

        # testing withou cycle:
        mesh_chi2r_im = chi2r_MP[:, iplt].reshape((n_m, n_p))
        mesh_rms_oc = rms_oc_MP[:, iplt].reshape((n_m, n_p))
        mesh_rv = rv_MP.reshape((n_m, n_p))

        scax = ax.contourf(
            mesh_p,
            mesh_m,
            mesh_chi2r_im,
            # linestyles='solid',
            # levels=clev, cmap=cm.viridis, zorder=4)
            levels=clev,
            colors=levcols,
            origin=origin_mesh,
            zorder=4,
        )
        print(" | mesh/contourf done")

        # colorbar
        cbar = plt.colorbar(
            scax, ax=ax, orientation="vertical", pad=0.01, format="%4.2e"
        )
        cbar.ax.set_ylabel(r"$\chi^{{2}}_\textrm{{r}}(O-C)$")
        cbar.ax.tick_params(labelsize=5)

        # cnrms = ax.contour(mesh_p, mesh_m, mesh_rms_oc,
        ax.contour(
            mesh_p,
            mesh_m,
            mesh_rms_oc,
            levels=rms_lev,
            # colors='black',
            colors="white",
            linewidths=0.5,
            # linestyles='solid',
            linestyles=[(0.0, (1.0, 1.0)), (0.0, (3.0, 2.0)), "-"],
            origin=origin_mesh,
            zorder=6,
        )
        print(" | mesh/contour rms done")

        # to be removed or use at my will!!
        # mesh_rv = anc.compute_Kms(1.030, mesh_m*cst.Mears*cst.Msjup, 88.55, mesh_p, 0.0)

        # inverted Greys cmap with alpha from 0 (transparent) to 1 (opaque)
        greys_r_t = cm.Greys_r(np.arange(cm.Greys_r.N))
        greys_r_t[:, -1] = np.linspace(0.0, 1.0, cm.Greys_r.N)
        greys_r_t = mpl.colors.ListedColormap(greys_r_t)

        v_rvmin = np.min(np.min(mesh_rv))
        if v_rvmin < 0.0:
            v_rvmin = 0.0
        v_rvmax = np.max(np.max(mesh_rv))
        # mesh_rv[mesh_rv < 3.] = np.ma.masked
        print((" | RV min = {0:.6f} m/s max = {1:.6f} m/s".format(v_rvmin, v_rvmax)))
        # cnrv = ax.contourf(mesh_p, mesh_m, mesh_rv,
        ax.contourf(
            mesh_p,
            mesh_m,
            mesh_rv,
            # levels=[3., 5., np.max(np.max(mesh_rv))],
            # colors=['lightgray', 'gray'], alpha=0.8,
            levels=[1.0, 2.0, 3.0, 5.0, v_rvmax],
            # cmap=cm.Greys_r, alpha=0.444,
            cmap=greys_r_t,
            # cmap=cm.Greys, alpha=0.7,
            # linestyles=[(0.,(1.,1.)), (0.,(3.,2.)),'-'],
            # origin=origin_mesh,
            vmin=v_rvmin,
            vmax=v_rvmax,
            zorder=7,
        )
        print(" | mesh/contourf rv done")

        ax.xaxis.set_tick_params(labelsize=10)
        ax.yaxis.set_tick_params(labelsize=10)
        # ax.set_ylabel('$\log M [M_\oplus]$', fontsize=6)
        ax.set_ylabel(r"$M (M_\oplus)$")
        ax.tick_params(direction="inout")
        # plt.draw()

        if P_tra is not None:
            # ax = plot_MMR(ax, P_tra)
            ax = plot_MMR(ax, 1.0)  # P_tra / P_tra
            ax.set_xlim([np.min(period) / P_tra, np.max(period) / P_tra])
            xlab = r"$P_\textrm{{per}}/P_\textrm{{tra}} , N_{{T0}} = \textrm{{{0:s}}}$".format(
                n_label[iplt]
            )

            # ll = [item.get_text() for item in ax.get_xticklabels()]
            # print '0', ll
            # tlabels = [ill.replace(u'\u2212', '-').replace('$','') for ill in ll]
            # print '1', tlabels
            # tlfloats = [u'%.2f' %(np.float64(ill)/P_tra) if ill != u'' else ill \
            # for ill in tlabels
            # ]
            # print '2', tlfloats
            # ax.set_xticklabels(tlfloats)

        else:
            ax.set_xlim([np.min(period), np.max(period)])
            xlab = r"$P (\textrm{{days}}) , N_{{T0}} = \textrm{{{0:s}}}$".format(
                n_label[iplt]
            )

        ax.set_ylim([np.min(mass), np.max(mass)])

        # ax.set_xlabel('$\log P [\\textrm{days}] , N_{T_0} = \\textrm{%s}$' %(n_label[iplt]), fontsize=6)

        ax.set_xlabel(xlab)

        # use only if log-scale
        plt.draw()
        # xtickslabels_obj = ax.get_xmajorticklabels()
        # xticks_lin_str = ticklabels_logtolin(xtickslabels_obj)
        # ax.set_xticklabels(xticks_lin_str)

        # ytickslabels_obj = ax.get_ymajorticklabels()
        # yticks_lin_str = ticklabels_logtolin(ytickslabels_obj)
        # ax.set_yticklabels(yticks_lin_str)
        # plt.draw()

        ix += iy

    # plt.tight_layout()
    fig.savefig(plot_file)
    plt.close(fig)
    print("*** plot_m_vs_p_vs_chi2r END ***")

    return


# ==============================================================================
def read_and_plot_p_vs_m_vs_chi2r(
    out_folder, hdf5_file, rms_oc_max=5.0 / 1440.0, K_rv=None, P_tra=None, of_log=None
):

    print(" - Read hdf5 file ... ", end=" ")
    f_h5 = h5py.File(hdf5_file, "r")
    grid_par = f_h5["grid/perturber_grid"][...]
    err_T0_d = f_h5["grid"].attrs["err_T0_day"]

    chi2r_oc_grid = f_h5["grid/chi2r_oc_grid"][...]
    rms_oc_grid = f_h5["grid/rms_oc_grid"][...]
    # rms_header = f_h5['grid/rms_header'][...]

    # above_threshold = f_h5['grid/above_threshold'][...]

    try:
        rms_rv = f_h5["grid/rms_rv"][...]
    except:
        rms_rv = None
    f_h5.close()
    print("done")

    # if(rms_min is not None):
    # rms_lower = rms_min
    # else:
    # rms_lower = err_T0_d

    # mass, period, freq_obs, freq_rv = compute_grid_to_freq_obs(grid_par, rms_oc_grid, lower_rms= rms_lower, upper_rms=rms_oc_max, K_rv=K_rv, rms_rv=rms_rv, of_log=of_log)

    print(" - Create mass, period grid vs chi2r/rmsOC/rmsRV ...", end=" ")
    mass, period, chi2r_MP, rms_oc_MP, rv_MP = compute_grid_to_chi2r_rms(
        grid_par, chi2r_oc_grid, rms_oc_grid, rms_rv=rms_rv, of_log=of_log
    )
    # mass, period, chi2r_MP, rms_oc_MP, rv_MP = compute_grid_to_chi2r_rms_MP_e0copmA90(grid_par,
    # chi2r_oc_grid, rms_oc_grid, rms_rv=rms_rv,
    # of_log=of_log)
    print("done")

    _, nrms = np.shape(rms_oc_grid)

    mean_sel = np.ones((nrms)).astype(bool)
    mean_sel[1:] = np.arange(1, nrms, 1) % 2 == 0

    median_sel = np.ones((nrms)).astype(bool)
    median_sel[1:] = np.arange(1, nrms, 1) % 2 == 1

    if not os.path.isdir(out_folder):
        os.makedirs(out_folder)

    plt_file_mean = os.path.join(
        out_folder,
        "p_vs_m_chi2r_{0:s}_max{1:.2f}min".format("mean", rms_oc_max * 1440.0),
    )
    plt_file_med = os.path.join(
        out_folder,
        "p_vs_m_chi2r_{0:s}_max{1:.2f}min".format("median", rms_oc_max * 1440.0),
    )

    if K_rv is not None and rms_rv is not None:
        plt_file_mean = "{0:s}_K{1:.2f}mps".format(plt_file_mean, K_rv)
        plt_file_med = "{0:s}_K{1:.2f}mps".format(plt_file_med, K_rv)

    plt_file_mean = "{0:s}.pdf".format(plt_file_mean)
    plt_file_med = "{0:s}.pdf".format(plt_file_med)

    print(" - plot mean MC ...")
    # plot_m_vs_p_vs_freq_obs(plt_file_mean, mass, period, freq_obs[:,mean_sel], freq_rv)
    plot_m_vs_p_vs_chi2r(
        plt_file_mean,
        mass,
        period,
        chi2r_MP[:, mean_sel],
        rms_oc_MP[:, mean_sel],
        rms_oc_max=rms_oc_max,
        K_rv=K_rv,
        rv_MP=rv_MP,
        P_tra=P_tra,
        n_label="all 3 5 6 10 15".split(),
    )
    print(plt_file_mean)

    print(" - plot median MC ...")
    # plot_m_vs_p_vs_freq_obs(plt_file_med, mass, period, freq_obs[:,median_sel], freq_rv)
    plot_m_vs_p_vs_chi2r(
        plt_file_med,
        mass,
        period,
        chi2r_MP[:, median_sel],
        rms_oc_MP[:, median_sel],
        rms_oc_max=rms_oc_max,
        K_rv=K_rv,
        rv_MP=rv_MP,
        P_tra=P_tra,
        n_label="all 3 5 6 10 15".split(),
    )
    print(plt_file_med)

    del (
        grid_par,
        err_T0_d,
        chi2r_oc_grid,
        rms_oc_grid,
        mass,
        period,
        chi2r_MP,
        rms_oc_MP,
        rv_MP,
    )
    if K_rv is not None and rms_rv is not None:
        del rms_rv
    gc.collect()

    print(" - COMPLETED PLOTS")

    return plt_file_mean, plt_file_med


# ==============================================================================
# ==============================================================================

# needed by pp
class object_sim:
    def __init__(self, sim_args):
        self.sim_args = sim_args


class object_res:
    def __init(self, sim_res):
        self.results_sim = sim_res


def run_the_simulation_pp(osim):

    sim_res = run_the_simulation(osim.sim_args)

    out_res = object_res(sim_res)

    return out_res


# ==============================================================================

# ==============================================================================
# MAIN


def main():

    cli = get_args()

    np.random.seed(cli.seed)

    out_folder, log_file, of_log = set_analysis_folder_log(
        cli.full_path, cli.sub_folder
    )

    anc.print_both("# THIS IS: {:s}".format(log_file), of_log)
    anc.print_both("Main folder: {:s}".format(cli.full_path), of_log)

    err_T0 = set_err_T0(cli.err_T0)
    anc.print_both("Mean T0 uncertainty = {0:.6f} days".format(err_T0), of_log)

    # t_launch = atime.Time('2018-06-30T06:00:00.0', format='isot', scale='tdb')
    # t_launch = atime.Time('2019-10-15T07:00:00.0', format='isot', scale='tdb')
    # t_end = t_launch.jd + 365.25*3.5

    anc.print_both(
        "CHEOPS DATE: LAUNCH = {0:.4f} [BJD] END = {1:.4f} [BJD]".format(
            cli.t_start, cli.t_end
        ),
        of_log,
    )
    anc.print_both(" ", of_log)

    # init TRADES
    anc.print_both(" - INIT TRADES", of_log)

    f90trades.initialize_trades(cli.full_path, cli.sub_folder, cli.nthreads)
    n_bodies = f90trades.n_bodies  # NUMBER OF TOTAL BODIES OF THE SYSTEM
    # n_planets = n_bodies - 1 # NUMBER OF PLANETS IN THE SYSTEM

    # read grid parameters ==> output ngrid and ncol
    ngrid, ncol = f90trades.wrapper_read_grid(1)
    # create grid
    perturber_grid = f90trades.wrapper_grid_init(1, ngrid, ncol)
    anc.print_both(" - SET PERTURBER GRID", of_log)
    # perturber_grid = f90trades.wrapper_grid_init(1)
    anc.print_both(" | To run {:d} simulation/s.\n".format(ngrid), of_log)

    # mr_convert, mr_unit = mass_factor(cli.mr_type)
    hdf5_name, f_h5, _, _, names_lst, _, units = create_hdf5_file(
        out_folder,
        cli.seed,
        cli.perturber_id,
        perturber_grid,
        cli.mr_type,
        cli.t_start,
        cli.t_end,
        err_T0,
        of_log,
    )
    anc.print_both(" - created hdf5 file: {:s}".format(hdf5_name), of_log)

    # define the number of transits to select:
    nT0_sel = [3, 5, 6, 10, 15]
    # define the number of Monte-Carlo repeatition for each nT0_sel
    # nMC = 1000
    nMC = 100
    n_rms = 1 + len(nT0_sel) * 2  # rms_ttv + nT0_sel[mean, median]

    transit_id = 2

    niter = ngrid
    # niter = 500

    chi2r_oc_grid = np.zeros((niter, n_rms))
    rms_oc_grid = np.zeros((niter, n_rms))
    above_threshold = np.zeros((niter, n_rms)).astype(bool)
    rms_rv_grid = -np.ones((niter))

    # 2018-03-08 ADDING GLS AMPLITUDE AND PERIOD OF TTV/OC
    attv_grid = np.zeros((niter))
    pttv_grid = np.zeros((niter))

    ## WORKING IN SERIAL
    # for igrid in range(0, niter):

    # results_sim = simulation_and_ttv_analysis(igrid, transit_id, cli.perturber_id, n_bodies, perturber_grid[igrid,:], t_launch.jd, t_end, err_T0, nT0_sel, names_lst, units, cli.mr_type, nMC=nMC, of_log=of_log)
    # f_h5 = save_results_sim(f_h5, results_sim, cli.mr_type)
    ##get_rms_from_results(results_sim, rms_oc_grid, above_threshold, nT0_sel, err_T0, rms_rv_grid)
    # get_chi2r_rms_from_results(results_sim, chi2r_oc_grid, rms_oc_grid, above_threshold, nT0_sel, err_T0, rms_rv_grid)
    # anc.print_both(' - Saved sim %06d into hdf5 file.' %(results_sim['sim_num']), of_log)

    print(" - Creating sim_args list ...", end=" ")

    list_sim_args = []
    # mp Queue
    # qq = mpro.Queue()
    # oo = mpro.Queue()
    list_osim_args = []

    for igrid in range(0, niter):
        sim_args = {}
        sim_args["i_grid"] = igrid
        sim_args["transit_id"] = transit_id
        sim_args["perturber_id"] = cli.perturber_id
        sim_args["perturber_par"] = perturber_grid[igrid, :]
        sim_args["n_bodies"] = n_bodies
        sim_args["t_start"] = cli.t_start
        sim_args["t_end"] = cli.t_end
        sim_args["err_T0"] = err_T0
        sim_args["nT0_sel"] = nT0_sel
        sim_args["names_lst"] = names_lst
        sim_args["units"] = units
        sim_args["mr_type"] = cli.mr_type
        sim_args["nMC"] = nMC
        list_sim_args.append(sim_args)
        list_osim_args.append(object_sim(sim_args))
        # qq.put(sim_args)
    print("done")

    print(" - *** RUN THEM ALL *** ...\n")

    if cli.nthreads <= 1:
        # serial
        list_results_sim = [
            run_the_simulation(list_sim_args[igrid]) for igrid in range(niter)
        ]
    else:

        ## **POOL** ##
        ## try InterruptiblePool from emcee to reduce process spawn...
        # threads_pool = emcee_pool.InterruptiblePool(cli.nthreads)
        ## TRY ThreadPool
        # threads_pool = ThreadPool(processes=cli.nthreads)

        # threads_pool = mpro.Pool(cli.nthreads)
        # list_results_sim = threads_pool.map(run_the_simulation, list_sim_args)
        # # close the pool of threads
        # threads_pool.close()
        # # threads_pool.terminate()
        # threads_pool.join()
        with mpro.Pool(cli.nthreads) as threads_pool:
            list_results_sim = threads_pool.map(run_the_simulation, list_sim_args)

    print("done\n")

    # save results of each sim in hdf5 file
    for igrid in range(0, niter):
        results_sim = list_results_sim[igrid]
        attv_grid[igrid] = results_sim["A_TTV"]
        pttv_grid[igrid] = results_sim["P_TTV"]

        print(
            " - Saving results for sim #{0:06d} ... into hdf5 file ...".format(
                results_sim["sim_num"]
            ),
            end=" ",
        )
        f_h5 = save_results_sim(f_h5, results_sim, cli.mr_type)

        print(" | and populating the chi2r and rms grid ...", end=" ")
        get_chi2r_rms_from_results(
            results_sim,
            chi2r_oc_grid,
            rms_oc_grid,
            above_threshold,
            nT0_sel,
            err_T0,
            rms_rv_grid,
        )
        print(" | done")
    print()

    # f_h5 = save_chi2r_rms_analysis(f_h5, chi2r_oc_grid, rms_oc_grid, above_threshold, nT0_sel, rms_rv_grid)
    f_h5 = save_chi2r_rms_analysis(
        f_h5,
        chi2r_oc_grid,
        rms_oc_grid,
        above_threshold,
        nT0_sel,
        rms_rv_grid,
        attv_grid,
        pttv_grid,
    )

    # f_h5 = save_rms_analysis(f_h5, rms_oc_grid, above_threshold, nT0_sel, rms_rv_grid)

    f_h5.close()

    anc.print_both(" ", of_log)
    anc.print_both(" ===============================", of_log)
    anc.print_both(" - Saved chi2r/rms_oc_grid into hdf5 file", of_log)
    anc.print_both(" ", of_log)

    # --------------------

    for imin in [15.0]:
        rms_threshold = np.float64(imin) / 1440.0
        _, _ = read_and_plot_p_vs_m_vs_chi2r(
            out_folder,
            hdf5_name,
            rms_oc_max=rms_threshold,
            K_rv=cli.K_rv,
            P_tra=cli.period_transiting,
            of_log=of_log,
        )

    of_log.close()

    gc.collect()

    return


# ==============================================================================
# RUN MAIN

if __name__ == "__main__":
    main()

# ==============================================================================
