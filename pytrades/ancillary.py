#!/usr/bin/env python
# -*- coding: utf-8 -*-

# no more "zero" integer division bugs!:P
import sys
import argparse
import os
import numpy as np  # array
import h5py
import sys
import yaml

# import trades_lib
# import random

import constants as cst  # local constants module

# from emcee import autocorr as acor
import emcee

# import acor
import matplotlib.pyplot as plt

import glob
import shutil

# import scipy.optimize as sciopt
# import scipy.odr as sodr
# try:
#     import sklearn.linear_model as sklm
# except:
#     print("sklearn not installed ... ")
# # import scipy.stats as scits
# try:
#     import statsmodels.api as sm
# except:
#     print("statsmodels not installed ...")

# common variables needed to create labels and parameter names
kel_fmt = ["%.3f", "%.4f", "%.3f", "%.1f", "%.1f", "%.1f", "%.3f", "%.3f"]

letters = ["a", "b", "c", "d", "e", "f", "g", "h", "i", "l", "m", "n", "o", "p"]

# deg2rad = np.pi / 180.0
# rad2deg = 180.0 / np.pi
deg2rad = cst.deg2rad
rad2deg = cst.rad2deg

eps64bit = np.finfo(np.float64(1.0)).eps
eps32bit = np.finfo(np.float32(1.0)).eps

# sigma =          1      0    -1       1     -2       2     -3       3
# sigma pos        0      1     2       3      4       5      6       7
percentile_val = [68.27, 50.0, 15.865, 84.135, 2.265, 97.735, 0.135, 99.865]

parameters_folders = [
    "0000_sim_initial",
    "0004_sim_pso",
    "0006_sim_de",
    "0777_sim_from_file",
    "1050_sim_median",
    "1051_sim_median_to_physical",
    "2050_sim_map",
    "2051_sim_map_to_physical",
    "3050_sim_mode",
    "0668_sim_map_hdi",
    "0669_sim_map_hdi_to_physical",
]

# ==============================================================================
def set_rcParams():

    my_dpi = 300
    plt.rcParams["text.usetex"] = False  # True
    plt.rcParams["font.family"] = "serif"  #'sans-serif'
    plt.rcParams["font.serif"] = ["Computer Modern Roman", "Palatino", "DejaVu Serif"]
    plt.rcParams["mathtext.fontset"] = "cm"
    plt.rcParams["figure.figsize"] = [5, 5]
    plt.rcParams["figure.facecolor"] = "white"
    plt.rcParams["savefig.facecolor"] = "white"
    plt.rcParams["figure.dpi"] = my_dpi
    plt.rcParams["savefig.dpi"] = 300
    plt.rcParams["font.size"] = 10
    plt.rcParams["xtick.labelsize"] = 9
    plt.rcParams["ytick.labelsize"] = 9
    plt.rcParams["animation.html"] = "jshtml"
    plt.rcParams["axes.formatter.useoffset"] = False

    return


def print_both(line, output=None):

    print(line)
    if output is not None:
        output.write(line + "\n")

    return


def decode_list(alist):

    if not "st" in str(type(alist[0])):
        blist = [a.decode("utf-8") for a in alist]
    else:
        blist = alist

    return blist


def encode_list(alist):

    if not "bytes" in str(type(alist[0])):
        blist = [a.encode("utf-8") for a in alist]
    else:
        blist = alist

    return blist


# ==============================================================================
# INITIALISE FOLDER AND LOG FILE
# ==============================================================================
def copy_simulation_files(dir_src, dir_dest):
    # copy all the simulation files from the source directory to destination directory
    ## arg.in
    arg_file = os.path.join(dir_src, "arg.in")
    shutil.copy(arg_file, os.path.join(dir_dest, ""))
    ## bodies.lst
    bodies_file = os.path.join(dir_src, "bodies.lst")
    shutil.copy(bodies_file, os.path.join(dir_dest, ""))
    ## bodies file (inside bodies.lst)
    obd = open(bodies_file, "r")
    for line in obd.readlines():
        shutil.copy(
            os.path.join(dir_src, line.strip().split()[0]), os.path.join(dir_dest, "")
        )
    obd.close()
    ## TT files
    t0files = glob.glob(os.path.join(dir_src, "NB*_observations.dat"))
    for t0f in t0files:
        shutil.copy(t0f, os.path.join(dir_dest, ""))
    ## RV file
    if os.path.exists(os.path.join(dir_src, "obsRV.dat")):
        shutil.copy(os.path.join(dir_src, "obsRV.dat"), os.path.join(dir_dest, ""))

    ## lm.opt pso.opt pik.opt
    opt_files = "lm.opt pso.opt pik.opt priors.in de.yml".split()
    for opt in opt_files:
        if os.path.exists(os.path.join(dir_src, opt)):
            shutil.copy(os.path.join(dir_src, opt), os.path.join(dir_dest, ""))

    return


# ==============================================================================


def init_folder(working_path, sub_folder):
    working_folder = os.path.join(working_path, sub_folder)
    if not os.path.isdir(working_folder):
        os.makedirs(working_folder)

    copy_simulation_files(working_path, working_folder)

    run_log = os.path.join(working_folder, "trades_run.log")
    of_run = open(run_log, "w")
    print_both("# pyTRADES LOG FILE", of_run)
    print_both("# working_path = %s" % (working_path), of_run)
    print_both("# working_folder = %s" % (working_folder), of_run)
    print_both("# run_log = %s" % (run_log), of_run)

    return working_folder, run_log, of_run


# ==============================================================================


def set_bool_argument(arg_in):
    if str(arg_in).lower() in ["t", "t", "tru", "true", "y", "ye", "yes", "1"]:
        arg_out = True
    else:
        arg_out = False
    return arg_out


# ==============================================================================


def set_int_argument(arg_in, default=0):
    try:
        arg_out = int(arg_in)
    except:
        arg_out = default
    return arg_out


# ==============================================================================


# def set_overplot(arg_in):
#     if str(arg_in).lower() == "none":
#         arg_out = None
#     else:
#         arg_out = os.path.abspath(arg_in)
#         if not os.path.isfile(arg_out):
#             try:
#                 arg_out = int(arg_in)
#             except:
#                 arg_out = None
#             if arg_out not in [0, 666, 667, 668, 777, 1050, 1051, 2050, 3050, 3051]:
#                 arg_out = None
#     return arg_out


def set_overplot(arg_in):

    arg_str = str(arg_in).lower().strip()
    arg_out = None
    if arg_str != "none":
        for p_folder in parameters_folders:
            q_folder = p_folder.split("sim_")[1]
            if arg_str == q_folder:
                arg_out = p_folder

    return arg_out


# ==============================================================================


def set_adhoc_file(arg_in):
    if str(arg_in).lower() == "none":
        arg_out = None
    else:
        arg_out = os.path.abspath(arg_in)
        if not os.path.isfile(arg_out):
            arg_out = None

    return arg_out


# ==============================================================================


def set_int_or_none(arg_in):
    try:
        arg_out = int(arg_in)
    except:
        arg_out = None
    return arg_out


# ==============================================================================

# read command line (cli) arguments
def get_args():
    parser = argparse.ArgumentParser(description="TRADES+EMCEE PLOT")

    parser.add_argument(
        "-p",
        "--path",
        action="store",
        dest="full_path",
        required=True,
        help="The path (absolute or relative) with simulation files for TRADES.",
    )

    parser.add_argument(
        "-nb",
        "--nb",
        "--nburn",
        "-nburn",
        "--nburnin",
        "-nburnin",
        action="store",
        dest="nburnin",
        required=True,
        help="The number of posterior/burn in steps to discard at the beginning of each chain. It has to be > 0",
    )

    parser.add_argument(
        "-m",
        "--mtype",
        "--mass-type",
        action="store",
        dest="m_type",
        default="e",
        help="Mass type: j = Jupiter, e = Earth, s = Sun. Default is Earth = e.",
    )

    parser.add_argument(
        "-t",
        "--temp",
        action="store",
        dest="temp_status",
        default=False,
        help="If you want to read temporary emcee analysis from file, because simulation is not finished yet. Default is False",
    )

    parser.add_argument(
        "-s",
        "--steps",
        "--gelmanrubin-steps",
        action="store",
        dest="gr_steps",
        default=10,
        help="Set to positive integer. It will allow to compute GelmanRubin statistics of only gr_steps, and not for each step. Just to debug and overall view of the chains. Default is 10.",
    )
    parser.add_argument(
        "-gk-steps",
        "--gk-steps",
        "-gelweke-steps",
        "--gelmanrubin-steps",
        action="store",
        dest="gk_steps",
        default=10,
        help="Set to positive integer. It will allow to compute Geweke statistics of only gk_steps, and not for each step. Just to debug and overall view of the chains. Default is 10.",
    )

    parser.add_argument(
        "-u",
        "--thin",
        "--use-thin",
        action="store",
        dest="use_thin",
        default=False,
        help="Set if you want to use or not the thinning parameter computed from the autocorrelation time. Default False",
    )

    parser.add_argument(
        "--sample",
        action="store",
        dest="sample_str",
        default="None",
        help='Set a list of parameter names to select a parameter sample within the Credible Intervals, first will be the "pivoting" parameters, then the parameters to not take into account. Default = None, it means first parameters fitted as "pivot", all other have to be within the Credible Intervals.',
    )

    parser.add_argument(
        "--seed",
        action="store",
        dest="seed",
        default="None",
        help="Set the seed for random number generator. Default is None.",
    )

    parser.add_argument(
        "--overplot",
        action="store",
        dest="overplot",
        default="None",
        help="""
            Define which parameter set to be overplotted on the correlation plot.\n
            Input the proper string:\n
            initial,\n
            adhoc,\n
            pso,\n
            de,\n
            median,\n
            map,\n
            mode,\n
            map_hdi.
            """,
    )

    parser.add_argument(
        "--ad-hoc",
        action="store",
        dest="adhoc",
        default="None",
        help="Define adhoc file with fitted parameters to be overplotted as adhoc",
    )

    parser.add_argument(
        "-n-samples",
        "--n-samples",
        action="store",
        dest="n_samples",
        default=0,
        help="Number of sample parameter within CI to generate T0s and RVs to be overplotted in the O-C plots. Defautl is 0",
    )

    parser.add_argument(
        "-corner-type",
        "--corner-type",
        action="store",
        dest="corner_type",
        default="both",
        help="Corner plot type: custom, pygtc, both. Default both.",
    )

    # which analysis and plot to do arguments
    parser.add_argument(
        "-all",
        "--all",
        nargs="?",
        type=str,
        const="True",
        dest="all",
        help="Do the full analysis and plots",
    )

    parser.add_argument(
        "-save-posterior",
        "--save-posterior",
        "-save-post",
        "--save-post",
        nargs="?",
        type=str,
        const="True",
        dest="save_posterior",
        help="Save posterior (fitted and physical) to file.",
    )

    parser.add_argument(
        "-save-parameters",
        "--save-parameters",
        "-save-par",
        "--save-par",
        nargs="?",
        type=str,
        const="True",
        dest="save_parameters",
        help="Save best-fit parameters to file. They will be read by overplot keyword.",
    )

    parser.add_argument(
        "-chain",
        "--chain",
        "-chains",
        "--chains",
        "-trace",
        "--trace",
        nargs="?",
        type=str,
        const="False",
        dest="chain",
        help="Chains/trace plots",
    )

    parser.add_argument(
        "-gelman-rubin",
        "--gelman-rubin",
        "-gr",
        "--gr",
        nargs="?",
        type=str,
        const="False",
        dest="gelman_rubin",
        help="Gelman-Rubin plots",
    )

    parser.add_argument(
        "-geweke",
        "--geweke",
        "-ge",
        "--ge",
        nargs="?",
        type=str,
        const="False",
        dest="geweke",
        help="Geweke plots",
    )

    parser.add_argument(
        "-correlation-fitted",
        "--correlation-fitted",
        "-triangle-fitted",
        "--triangle-fitted",
        nargs="?",
        type=str,
        const="False",
        dest="correlation_fitted",
        help="Correlation plot of fitted parameters",
    )

    parser.add_argument(
        "-correlation-physical",
        "--correlation-physical",
        "-triangle-physical",
        "--triangle-physical",
        nargs="?",
        type=str,
        const="False",
        dest="correlation_physical",
        help="Correlation plot of physical parameters",
    )

    cli = parser.parse_args()

    cli.full_path = os.path.abspath(cli.full_path)
    cli.nburnin = set_int_or_none(cli.nburnin)
    cli.m_type = cli.m_type.lower()
    cli.temp_status = set_bool_argument(cli.temp_status)
    # cli.boot_id = set_int_argument(cli.boot_id)
    cli.use_thin = set_int_argument(cli.use_thin, default=False)
    # cli.use_thin = set_bool_argument(cli.use_thin)
    cli.seed = set_int_or_none(cli.seed)
    # cli.overplot = set_overplot(cli.overplot)
    cli.adhoc = set_adhoc_file(cli.adhoc)
    # cli.pyscript = os.path.abspath(cli.pyscript)
    cli.n_samples = set_int_argument(cli.n_samples, default=0)

    cli.all = set_bool_argument(cli.all)
    if cli.all:
        cli.save_posterior = True
        cli.save_parameters = True
        cli.chain = True
        cli.gelman_rubin = True
        cli.geweke = True
        cli.correlation_fitted = True
        cli.correlation_physical = True
    else:
        cli.save_posterior = set_bool_argument(cli.save_posterior)
        cli.save_parameters = set_bool_argument(cli.save_parameters)
        cli.chain = set_bool_argument(cli.chain)
        cli.gelman_rubin = set_bool_argument(cli.gelman_rubin)
        cli.geweke = set_bool_argument(cli.geweke)
        cli.correlation_fitted = set_bool_argument(cli.correlation_fitted)
        cli.correlation_physical = set_bool_argument(cli.correlation_physical)

    return cli


# ==============================================================================


def get_input_file():

    parser = argparse.ArgumentParser(description="TRADES")

    parser.add_argument(
        "-input",
        "--input",
        "-input-file",
        "--input-file",
        action="store",
        dest="input",
        required=True,
        help="The path (absolute or relative) with yaml configuration file (for run or analysis) for TRADES.",
    )

    cli = parser.parse_args()

    input_file = os.path.abspath(cli.input)

    return input_file


# ==============================================================================
class CLI_OC:
    def __init__(
        self,
        full_path=os.path.abspath("."),
        idsim=1,
        lmflag=0,
        tscale=2440000.5,
        ocunit="d",
        samples_file=None,
        limits="obs",
        kep_ele=False,
    ):
        self.full_path = os.path.abspath(full_path)
        self.idsim = int(idsim)
        self.sim_type = os.path.basename(full_path)
        self.lmflag = int(lmflag)
        self.tscale = tscale
        self.ocunit = ocunit
        self.samples_file = samples_file
        self.limits = limits
        self.kep_ele = kep_ele

        return


class CLI_RV:
    def __init__(
        self,
        full_path=os.path.abspath("."),
        idsim=1,
        lmflag=0,
        tscale=2440000.5,
        labels="None",
        samples_file=None,
        limits="obs",
        color_map="viridis",
    ):
        self.full_path = os.path.abspath(full_path)
        self.sim_type = os.path.basename(full_path)
        self.idsim = int(idsim)
        self.lmflag = int(lmflag)
        self.tscale = tscale
        if str(labels).lower() != "none":
            if isinstance(labels, str):
                self.labels = labels.split()
            else:
                self.labels = labels
        else:
            self.labels = None
        self.samples_file = samples_file
        self.limits = limits
        self.color_map = color_map

        return


# ==============================================================================


class ConfigurationRun:
    def __init__(self, yaml_file):

        with open(yaml_file) as f:
            conf_all = yaml.load(f, Loader=yaml.FullLoader)

        conf = {
            "full_path": os.path.abspath("."),
            "sub_folder": "",
            "seed": 42,
            "mass_type": "e",
            "delta_sigma": 1.0e-4,
            "nthreads": 1,
            "trades_previous": None,
        }
        conf_keys = conf.keys()

        cemcee = {
            "nwalkers": 42,
            "nruns": 1,
            "thin_by": 1,
            "emcee_restart": False,
            "emcee_progress": True,
            "move": {"type": ["ai"], "fraction": [1.0]}
        }
        emcee_keys = cemcee.keys()

        cpyde = {
            "type": False,  # False/no, run, resume, to_emcee
            "npop": 42,
            "ngen": 4200,
            "save": 100,
            "f": 0.5,
            "c": 0.5,
            "maximize": True,
        }
        pyde_keys = cpyde.keys()

        cultranest = {
            "ultranest_live_points": 400,
            "ultranest_resume_flag": "resume",  # 'resume', 'resume-similar', 'overwrite' or 'subfolder'
            "ultranest_dlogz": 0.5,  # desired accuracy on logz
        }

        cdyne = {
            "dynesty_live_points": 500,
            "dynesty_bound": "multi",
            "dynesty_sample": "auto",
            "dynesty_restore": False,
            "dynesty_dlogz": 0.01,  # desired accuracy on logz
            "dynesty_pfrac": 0.8,
            "dynesty_just_plots": False,
        }
        dyne_keys = cdyne.keys()

        if "run" in conf_all:
            conf_input = conf_all["run"]
            for k, v in conf_input.items():
                if k in conf_keys:
                    conf[k] = v
                elif k == "pyde":
                    for kp, vp in conf_input["pyde"].items():
                        cpyde[kp] = vp
                elif k == "emcee":
                    for ke, ve in conf_input["emcee"].items():
                        cemcee[ke] = ve
                elif k == "ultranest":
                    for ku, vu in conf_input["ultranest"].items():
                        cultranest["ultranest_{}".format(ku)] = vu
                elif k == "dynesty":
                    for kd, vd in conf_input["dynesty"].items():
                        cdyne["dynesty_{}".format(kd)] = vd

        self.full_path = os.path.join(os.path.abspath(conf["full_path"]), "")
        self.sub_folder = os.path.join(self.full_path, conf["sub_folder"], "")
        self.seed = set_int_or_none(conf["seed"])
        self.m_type = conf["mass_type"]
        self.delta_sigma = conf["delta_sigma"]
        self.nthreads = set_int_argument(conf["nthreads"], default=1)
        self.trades_previous = set_adhoc_file(conf["trades_previous"])

        # pyde
        self.de_type = str(cpyde["type"]).lower()
        self.npop_de = set_int_argument(cpyde["npop"], default=42)
        self.ngen_de = set_int_argument(cpyde["ngen"], default=4200)
        self.nsave_de = set_int_argument(cpyde["save"], default=100)
        self.de_f = cpyde["f"]
        self.de_c = cpyde["c"]
        self.de_maximize = set_bool_argument(cpyde["maximize"])

        # emcee
        self.nwalkers = set_int_argument(cemcee["nwalkers"], default=42)
        self.nruns = set_int_argument(cemcee["nruns"], default=1)
        self.thin_by = set_int_argument(cemcee["thin_by"], default=1)
        self.emcee_restart = set_bool_argument(cemcee["emcee_restart"])
        self.emcee_progress = set_bool_argument(cemcee["emcee_progress"])

        # print()
        # print("====================")
        # print("cemcee['move']", cemcee["move"])
        if cemcee["move"] is None:
            self.emcee_move = [(emcee.moves.StretchMove(), 1.0)]
        else:
            emoves = []
            nmove = len(cemcee["move"]["type"])
            for i_em in range(nmove):
                em = cemcee["move"]["type"][i_em]
                # print("em = ", em)
                if em.lower() == "de":
                    m = emcee.moves.DEMove()
                elif em.lower() == "desnooker":
                    m = emcee.moves.DESnookerMove()
                else:
                    m = emcee.moves.StretchMove()
                fm = cemcee["move"]["fraction"][i_em]
                # print("==> ", (m, fm))
                emoves.append((m, fm))
            self.emcee_move = emoves
        # print("====================")
        # print()


        # ultranest
        self.ultranest_live_points = set_int_argument(
            cultranest["ultranest_live_points"], default=500
        )
        self.resume_flag = cultranest["ultranest_resume_flag"]
        self.dlogz = cultranest["ultranest_dlogz"]

        # dynesty
        self.dynesty_live_points = set_int_argument(
            cdyne["dynesty_live_points"], default=500
        )
        self.dynesty_bound = cdyne["dynesty_bound"]
        self.dynesty_sample = cdyne["dynesty_sample"]
        self.dynesty_restore = cdyne["dynesty_restore"]
        if str(cdyne["dynesty_dlogz"]).lower() == "none":
            self.dynesty_dlogz = 0.01
        else:
            self.dynesty_dlogz = float(cdyne["dynesty_dlogz"])
        self.dynesty_pfrac = float(cdyne["dynesty_pfrac"])
        self.dynesty_just_plots = set_bool_argument(cdyne["dynesty_just_plots"])

        return


# ==============================================================================


class ConfigurationAnalysis:
    def __init__(self, yaml_file):
        self.yaml_file = yaml_file
        with open(yaml_file) as f:
            conf_all = yaml.load(f, Loader=yaml.FullLoader)

        conf_analysis = {
            "full_path": os.path.abspath("."),
            "m_type": "e",
            "nburnin": 0,
            "emcee_old_save": False,
            "temp_status": True,
            "use_thin": 1,
            "seed": None,
            "from_file": None,
            "n_samples": 0,
            "corner_type": "pygtc",
            "overplot": "map_hdi",
            "all_analysis": False,
            "save_posterior": False,
            "save_parameters": False,
            "chain": False,
            "gelman_rubin": False,
            "gr_steps": 10,
            "geweke": False,
            "gk_steps": 10,
            "correlation_fitted": False,
            "correlation_physical": False,
        }
        conf_keys = conf_analysis.keys()

        # conf = conf_all["analysis"]
        if "analysis" in conf_all.keys():
            conf_input = conf_all["analysis"]
            for k, v in conf_input.items():
                if k in conf_keys:
                    conf_analysis[k] = v

        self.full_path = os.path.abspath(conf_analysis["full_path"])

        self.m_type = str(conf_analysis["m_type"]).lower()
        self.nburnin = set_int_or_none(conf_analysis["nburnin"])
        self.emcee_old_save = set_bool_argument(conf_analysis["emcee_old_save"])
        self.temp_status = set_bool_argument(conf_analysis["temp_status"])
        self.use_thin = set_int_argument(conf_analysis["use_thin"], default=1)
        self.seed = set_int_or_none(conf_analysis["seed"])
        self.from_file = set_adhoc_file(conf_analysis["from_file"])
        self.n_samples = set_int_argument(conf_analysis["n_samples"], default=0)
        self.corner_type = str(conf_analysis["corner_type"]).lower().strip()

        self.overplot = str(conf_analysis["overplot"]).lower().strip()
        if self.overplot == "mle":
            self.overplot = "map_hdi"

        self.gr_steps = set_int_argument(conf_analysis["gr_steps"], default=10)
        self.gk_steps = set_int_argument(conf_analysis["gk_steps"], default=10)

        self.all = set_bool_argument(conf_analysis["all_analysis"])
        if self.all:
            self.save_posterior = True
            self.save_parameters = True
            self.chain = True
            self.gelman_rubin = True
            self.geweke = True
            self.correlation_fitted = True
            self.correlation_physical = True
        else:
            self.save_posterior = set_bool_argument(conf_analysis["save_posterior"])
            self.save_parameters = set_bool_argument(conf_analysis["save_parameters"])
            self.chain = set_bool_argument(conf_analysis["chain"])
            self.gelman_rubin = set_bool_argument(conf_analysis["gelman_rubin"])
            self.geweke = set_bool_argument(conf_analysis["geweke"])
            self.correlation_fitted = set_bool_argument(
                conf_analysis["correlation_fitted"]
            )
            self.correlation_physical = set_bool_argument(
                conf_analysis["correlation_physical"]
            )

        # === OC === #
        conf_oc = {
            "plot_oc": False,
            "full_path": os.path.abspath("."),
            "sim_name": ["map_hdi", "map", "median", "initial"],
            "idplanet_name": None,
            "lmflag": 0,
            "tscale": None,
            "unit": "auto",
            "samples_file": None,
            "limits": "obs",
            "kep_ele": False,
        }
        # print("conf_oc[sim_name]    = {}".format(conf_oc["sim_name"]))

        conf_keys = conf_oc.keys()
        if "OC" in conf_all.keys():
            conf_input = conf_all["OC"]
            for k, v in conf_input.items():
                if k in conf_keys:
                    conf_oc[k] = v

        # print("conf_input[sim_name] = {}".format(conf_input["sim_name"]))
        # print("conf_oc[sim_name]    = {}".format(conf_oc["sim_name"]))
        # print("conf_oc[samples_file]    = {} is {}".format(conf_oc["samples_file"], type(conf_oc["samples_file"])))

        self.plot_oc = conf_oc["plot_oc"]
        self.idplanet_name = conf_oc["idplanet_name"]
        fpath = os.path.abspath(conf_oc["full_path"])
        self.ocs = []
        # check if mle and change it to map_hdi
        for i_name, sname in enumerate(conf_oc["sim_name"]):
            if sname == "mle":
                conf_oc["sim_name"][i_name] = "map_hdi"
        for keyn in ["map", "map_hdi", "median"]:
            keyp = "{}_to_physical".format(keyn)
            if (keyn in conf_oc["sim_name"]) and not (keyp in conf_oc["sim_name"]):
                conf_oc["sim_name"].append(keyp)

        if str(conf_oc["samples_file"]).lower() != "none":
            samples_file = os.path.join(fpath, conf_oc["samples_file"])
        else:
            samples_file = None
        for sname in conf_oc["sim_name"]:
            xname = "sim_{}".format(sname)
            for pfolder in parameters_folders:
                if xname.split("sim_")[1] == pfolder.split("sim_")[1]:
                    fsim = pfolder
                    isim = int(pfolder.split("_sim")[0])
                    self.ocs.append(
                        CLI_OC(
                            full_path=os.path.join(fpath, fsim),
                            idsim=isim,
                            lmflag=conf_oc["lmflag"],
                            tscale=conf_oc["tscale"],
                            ocunit=conf_oc["unit"],
                            samples_file = samples_file,
                            limits=conf_oc["limits"],
                            kep_ele=conf_oc["kep_ele"],
                        )
                    )
        # print("OC OBJ")
        # for oco in self.ocs:
        #     print(oco.idsim, oco.full_path)

        # === RV === #
        conf_rv = {
            "plot_rv": False,
            "full_path": os.path.abspath("."),
            "sim_name": ["map_hdi", "map", "median", "initial"],
            "lmflag": 0,
            "tscale": None,
            "samples_file": None,
            "limits": "obs",
            "labels": None,
            "color_map": "nipy_spectral",
        }
        conf_keys = conf_rv.keys()
        if "RV" in conf_all.keys():
            conf_input = conf_all["RV"]
            for k, v in conf_input.items():
                if k in conf_keys:
                    conf_rv[k] = v

        self.plot_rv = conf_rv["plot_rv"]
        if self.plot_rv:
            fpath = os.path.abspath(conf_rv["full_path"])
            self.rvs = []
            # check if mle and change it to map_hdi
            for i_name, sname in enumerate(conf_rv["sim_name"]):
                if sname == "mle":
                    conf_rv["sim_name"][i_name] = "map_hdi"
            for keyn in ["map", "map_hdi", "median"]:
                keyp = "{}_to_physical".format(keyn)
                if (keyn in conf_rv["sim_name"]) and not (keyp in conf_rv["sim_name"]):
                    conf_rv["sim_name"].append(keyp)

            if str(conf_rv["samples_file"]).lower() != "none":
                samples_file=os.path.join(fpath, conf_rv["samples_file"])
            else:
                samples_file = None
            for sname in conf_rv["sim_name"]:
                xname = "sim_{}".format(sname)
                for pfolder in parameters_folders:
                    if xname.split("sim_")[1] == pfolder.split("sim_")[1]:
                        fsim = pfolder
                        isim = int(pfolder.split("_sim")[0])

                        self.rvs.append(
                            CLI_RV(
                                full_path=os.path.join(fpath, fsim),
                                idsim=isim,
                                lmflag=conf_rv["lmflag"],
                                tscale=conf_rv["tscale"],
                                samples_file=samples_file,
                                limits=conf_rv["limits"],
                                labels=conf_rv["labels"],
                                color_map=conf_rv["color_map"],
                            )
                        )
            
            # print("RV OBJ")
            # for rvo in self.rvs:
            #     print(rvo.idsim, rvo.full_path)

        return


# ==============================================================================
def period_to_semimajoraxis(Ms_Msun, Mp_Msun, P_days):

    mu = cst.Giau * (Ms_Msun + Mp_Msun)
    P2 = P_days * P_days
    twopi2 = cst.dpi * cst.dpi
    sma = ((mu * P2) / twopi2) ** cst.onethird

    return sma


# ==============================================================================

# given the mass flag letter it computes the proper mass conversion factor
def mass_type_factor(Ms=1.0, mtype="earth", mscale=True):
    if mtype.lower() in ["s", "su", "sun"]:
        conv_factor = 1.0
        scale_factor = Ms
        mass_unit = "M_Sun"
    elif mtype.lower() in ["e", "ea", "ea", "eart", "earth"]:
        conv_factor = cst.Msear
        scale_factor = Ms * cst.Msear
        mass_unit = "M_Earth"
    elif mtype.lower() in ["n", "ne", "nep", "nept", "neptu", "neptun", "neptune"]:
        conv_factor = cst.Msnep
        scale_factor = Ms * cst.Msnep
        mass_unit = "M_Nep"
    else:
        conv_factor = cst.Msjup
        scale_factor = Ms * cst.Msjup
        mass_unit = "M_Jup"
    if mscale:
        return scale_factor, mass_unit
    else:
        return conv_factor, mass_unit


# given the mass/radius flag letter it computes the proper mass/radius conversion factor
def mass_radius_type_factor(mtype="earth"):
    if mtype.lower() in ["s", "su", "sun"]:
        mass_conv_factor = 1.0
        mass_unit = "M_Sun"
        radius_conv_factor = 1.0
        radius_unit = "R_Sun"
    elif mtype.lower() in ["e", "ea", "ea", "eart", "earth"]:
        mass_conv_factor = cst.Msear
        mass_unit = "M_Earth"
        radius_conv_factor = 1.0 / cst.Rears
        radius_unit = "R_Earth"
    elif mtype.lower() in ["n", "ne", "nep", "nept", "neptu", "neptun", "neptune"]:
        mass_conv_factor = cst.Msnep
        mass_unit = "M_Nep"
        radius_conv_factor = 1.0 / cst.Rneps
        radius_unit = "R_Nep"
    else:
        mass_conv_factor = cst.Msjup
        mass_unit = "M_Jup"
        radius_conv_factor = 1.0 / cst.Rjups
        radius_unit = "R_Jup"

    return mass_conv_factor, mass_unit, radius_conv_factor, radius_unit


# ==============================================================================

# prepare the labels of the Keplerian orbital elements for the legend
def keplerian_legend(parameter_names, m_type):
    # _, m_unit = mass_type_factor(1.0, m_type, mscale=False)
    _, m_unit, _, _ = mass_radius_type_factor(mtype=m_type)

    if m_unit == "M_Jup":
        unit_mass = "M_\mathrm{Jup}"
        unit_radius = "R_\mathrm{Jup}"
    elif m_unit == "M_Earth":
        unit_mass = "M_\oplus"
        unit_radius = "R_\oplus"
    elif m_unit == "M_Sun":
        unit_mass = "M_\odot"
        unit_radius = "R_\odot"
    elif m_unit == "M_Nep":
        unit_mass = "M_\mathrm{Nep}"
        unit_radius = "R_\mathrm{Nep}"
    else:
        unit_mass = " "
        unit_radius = " "

    nfit = np.shape(parameter_names)[0]
    # kel_legends = np.zeros((nfit), dtype="|S256")
    kel_legends = [""] * nfit

    for i in range(0, nfit):
        parameter_names[i] = parameter_names[i].strip()

        if "Ms" in parameter_names[i]:
            planet_id = int(parameter_names[i].split("m")[1].split("Ms")[0]) - 1
            # kel_legends[i] = r'%s$_\mathrm{%s}$ %s' %(kel_id_2[0], letters[planet_id], kel_units_2[0])
            kel_legends[i] = (
                # "$\displaystyle\frac{M_\mathrm{%s}}{M_\star} \left[\frac{%s}{M_\odot}\right]$"
                # "$\frac{M_\mathrm{%s}}{M_\star} \left(\frac{%s}{M_\odot}\right)$"
                "$M_\mathrm{%s}/M_\star ({%s}/M_\odot)$"
                % (letters[planet_id], unit_mass)
            )

        elif parameter_names[i][0] == "m" and parameter_names[i][-1] != "s":
            planet_id = int(parameter_names[i][1:]) - 1
            kel_legends[i] = "$M_\mathrm{%s} [%s]$" % (letters[planet_id], unit_mass)

        elif parameter_names[i][0] == "r":
            planet_id = int(parameter_names[i][1:]) - 1
            kel_legends[i] = "$R_\mathrm{%s}$ [%s]$" % (
                letters[planet_id],
                unit_radius,
            )

        elif parameter_names[i][0] == "P":
            planet_id = int(parameter_names[i][1:]) - 1
            kel_legends[i] = "$P_\mathrm{%s}$ [days]" % (letters[planet_id])

        elif "e" in parameter_names[i]:
            if parameter_names[i][0:4] == "ecos":
                planet_id = int(parameter_names[i].split("w")[1].strip()) - 1
                kel_legends[i] = "$e\cos\omega_\mathrm{%s}$" % (letters[planet_id])
            elif parameter_names[i][0:6] == "sqrtec":
                planet_id = int(parameter_names[i].split("w")[1].strip()) - 1
                kel_legends[i] = "$\sqrt{e}\cos\omega_\mathrm{%s}$" % (
                    letters[planet_id]
                )
            elif parameter_names[i][0:4] == "seco":
                planet_id = int(parameter_names[i].split("w")[1].strip()) - 1
                kel_legends[i] = "$\sqrt{e}\cos\omega_\mathrm{%s}$" % (
                    letters[planet_id]
                )
            elif parameter_names[i][0:4] == "esin":
                planet_id = int(parameter_names[i].split("w")[1].strip()) - 1
                kel_legends[i] = "$e\sin\omega_\mathrm{%s}$" % (letters[planet_id])
            elif parameter_names[i][0:6] == "sqrtes":
                planet_id = int(parameter_names[i].split("w")[1].strip()) - 1
                kel_legends[i] = "$\sqrt{e}\sin\omega_\mathrm{%s}$" % (
                    letters[planet_id]
                )
            elif parameter_names[i][0:4] == "sesi":
                planet_id = int(parameter_names[i].split("w")[1].strip()) - 1
                kel_legends[i] = "$\sqrt{e}\sin\omega_\mathrm{%s}$" % (
                    letters[planet_id]
                )
            else:
                planet_id = int(parameter_names[i].split("e")[1].strip()) - 1
                kel_legends[i] = "$e_\mathrm{%s}$" % (letters[planet_id])

        elif parameter_names[i][0] == "w":
            planet_id = int(parameter_names[i].split("w")[1].strip()) - 1
            kel_legends[i] = "$\omega_\mathrm{%s} [^\circ]$" % (letters[planet_id])

        elif parameter_names[i][0:2] == "mA":
            planet_id = int(parameter_names[i][2:]) - 1
            kel_legends[i] = "$\mathcal{M}_\mathrm{%s} [^\circ]$" % (letters[planet_id])

        elif "lambda" in parameter_names[i]:
            planet_id = int(parameter_names[i].split("lambda")[1].strip()) - 1
            kel_legends[i] = "$\lambda_\mathrm{%s} [^\circ]$" % (letters[planet_id])

        # elif parameter_names[i][0] == "i":
        #     if "icos" in parameter_names[i]:
        #         planet_id = int(parameter_names[i].split("lN")[1].strip()) - 1
        #         kel_legends[i] = "$i\cos\Omega_\mathrm{%s}$" % (letters[planet_id])
        #     elif "isin" in parameter_names[i]:
        #         planet_id = int(parameter_names[i].split("lN")[1].strip()) - 1
        #         kel_legends[i] = "$\i\sin\Omega_\mathrm{%s}$" % (letters[planet_id])
        #     else:
        #         planet_id = int(parameter_names[i][1:]) - 1
        #         kel_legends[i] = "$i_\mathrm{%s} [^\circ]$" % (letters[planet_id])
        elif parameter_names[i][0:2] == "ci" or "cosi" in parameter_names[i]:
            if "icos" in parameter_names[i]:
                planet_id = int(parameter_names[i].split("lN")[1].strip()) - 1
                kel_legends[i] = "$\cosi\cos\Omega_\mathrm{%s}$" % (letters[planet_id])
            elif "isin" in parameter_names[i]:
                planet_id = int(parameter_names[i].split("lN")[1].strip()) - 1
                kel_legends[i] = "$\cosi\sin\Omega_\mathrm{%s}$" % (letters[planet_id])
        elif parameter_names[i][0] == "i":
            planet_id = int(parameter_names[i][1:]) - 1
            kel_legends[i] = "$i_\mathrm{%s} [^\circ]$" % (letters[planet_id])

        elif parameter_names[i][0:2] == "lN":
            planet_id = int(parameter_names[i][2:]) - 1
            kel_legends[i] = "$\Omega_\mathrm{%s} [^\circ]$" % (letters[planet_id])

        else:
            kel_legends[i] = str(parameter_names[i])

    kel_legends = [
        "{0}".format(
            # kl.replace("[", "(").replace("]", ")").strip()
            kl.replace("[", "(").replace("]", ")")
        )
        for kl in decode_list(kel_legends)
    ]
    kel_legends = np.array(kel_legends)

    return kel_legends


# ==============================================================================


def derived_labels(derived_names, m_type):
    # _, m_unit = mass_type_factor(1.0, m_type, mscale=False)
    _, m_unit, _, _ = mass_radius_type_factor(mtype=m_type)

    if m_unit == "M_Jup":
        unit_mass = "M_\mathrm{Jup}"
        unit_radius = "R_\mathrm{Jup}"
    elif m_unit == "M_Earth":
        unit_mass = "M_\oplus"
        unit_radius = "R_\oplus"
    elif m_unit == "M_Sun":
        unit_mass = "M_\odot"
        unit_radius = "R_\odot"
    elif m_unit == "M_Nep":
        unit_mass = "M_\mathrm{Nep}"
        unit_radius = "R_\mathrm{Nep}"
    else:
        unit_mass = " "
        unit_radius = " "

    labels_list = []
    n_der = np.shape(derived_names)[0]
    for ider in range(0, n_der):

        # mass
        if derived_names[ider][0] == "m" and derived_names[ider][1] != "A":
            planet_id = int(derived_names[ider].split("m")[1]) - 1
            labels_list.append("$M_\mathrm{%s} [%s]$" % (letters[planet_id], unit_mass))

        # radius
        elif derived_names[ider][0] == "r":
            planet_id = int(derived_names[ider].split("r")[1]) - 1
            labels_list.append(
                "$R_\mathrm{%s} [%s]$" % (letters[planet_id], unit_radius)
            )

        # eccentricity
        elif derived_names[ider][0] == "e":
            planet_id = int(derived_names[ider].split("e")[1]) - 1
            labels_list.append("$e_\mathrm{%s}$" % (letters[planet_id]))

        # argument of pericentre
        elif derived_names[ider][0] == "w":
            planet_id = int(derived_names[ider].split("w")[1]) - 1
            # labels_list.append(r'$\omega_\mathrm{%s} [\mathrm{deg}]$' %(letters[planet_id]))
            labels_list.append("$\omega_\mathrm{%s} [^\circ]$" % (letters[planet_id]))

        # mean anomaly
        elif derived_names[ider][0:2] == "mA":
            planet_id = int(derived_names[ider].split("mA")[1]) - 1
            # labels_list.append(r'$\mathcal{M}_\mathrm{%s} [\mathrm{deg}]$' %(letters[planet_id]))
            labels_list.append(
                "$\mathcal{M}_\mathrm{%s} [^\circ]$" % (letters[planet_id])
            )

        # inclination
        elif derived_names[ider][0] == "i":
            planet_id = int(derived_names[ider].split("i")[1]) - 1
            # labels_list.append(r'$i_\mathrm{%s} [\mathrm{deg}]$' %(letters[planet_id]))
            labels_list.append("$i_\mathrm{%s} [^\circ]$" % (letters[planet_id]))

        # longitude of node
        elif derived_names[ider][0:2] == "lN":
            planet_id = int(derived_names[ider].split("lN")[1]) - 1
            # labels_list.append(r'$\Omega_\mathrm{%s} [\mathrm{deg}]$' %(letters[planet_id]))
            labels_list.append("$\Omega_\mathrm{%s} [^\circ]$" % (letters[planet_id]))

        else:
            labels_list.append(derived_names[ider])

    labels_list = [
        "{0}".format(l.replace("[", "(").replace("]", ")"))
        for l in decode_list(labels_list)
    ]

    return labels_list


# ==============================================================================

# only needed by check_good_parameters to convert Mjup to flagged mass
def check_good_parameters(good_id, good_parameters_in, m_factor, nfit):

    good_parameters_out = np.zeros(nfit) + good_parameters_in
    for i in range(0, nfit):
        if good_id[i].strip()[0] == "m" and good_id[i].strip()[1] != "A":
            good_parameters_out[i] = good_parameters_in[i] * cst.Mjups * m_factor

    return good_parameters_out


# ==============================================================================


def read_check_good_parameters(full_path, good_status, m_factor, nfit):

    # read good parameters file
    good_parameters_0, good_parameters_1 = False, False
    if good_status:
        good_file = os.path.join(full_path, "good_parameters.dat")
        if os.path.exists(good_file):
            good_parameters_0 = np.genfromtxt(good_file, usecols=(1), dtype=np.float64)
            good_id = np.genfromtxt(good_file, usecols=(0), dtype="|S10")
            good_parameters_1 = check_good_parameters(
                good_id, good_parameters_0, m_factor, nfit
            )

    return good_parameters_0, good_parameters_1


# ==============================================================================

# sometimes the angle distribution are bimodal because they are centered in the wrong position
# e.g.: -1° = 359° etc...
def rescale_angle(angle):

    # with arctan2
    cosa = np.cos(angle * deg2rad)
    sina = np.sin(angle * deg2rad)
    new_angle = np.arctan2(sina, cosa) * 180.0 / np.pi

    return new_angle


# ==============================================================================

# renormalize the angular parameter mean Anomaly mA
def renormalize_parameters(parameters, parameter_names):

    new_parameters = parameters.copy()
    # nfit = np.array(new_parameters).shape[0]
    nfit, _ = np.shape(new_parameters)
    for i in range(0, nfit):
        if parameter_names[i][:2] == "mA":
            new_parameters[i] = (parameters[i] + 360.0) % 360.0

    return new_parameters


# ==============================================================================

# get indices of max of an array
def get_max_indices(array_values):

    lmax = np.max(array_values)
    xx = np.where(array_values == lmax)
    idx = [xx[0][0], xx[1][0]]

    return idx[0], idx[1]


# ==============================================================================

# prepare file and best folder emcee
def get_emcee_file_and_best(emcee_folder, temp_status):

    if temp_status:
        folder_best = "best_temp"
        emcee_file = os.path.join(emcee_folder, "emcee_temp.hdf5")
        if not os.path.exists(emcee_file):
            emcee_file = os.path.join(emcee_folder, "emcee_summary.hdf5")
        emcee_best = os.path.join(emcee_folder, folder_best)
    else:
        folder_best = "best"
        emcee_file = os.path.join(emcee_folder, "emcee_summary.hdf5")
        emcee_best = os.path.join(emcee_folder, folder_best)

    return emcee_file, emcee_best, folder_best


# ==============================================================================


def get_devol_file(emcee_folder):

    devol_file = os.path.join(emcee_folder, "best_devol.hdf5")

    return devol_file


# ==============================================================================


def get_percentile_angle(angle_posterior):

    cosa_posterior = np.cos(angle_posterior * np.pi / 180.0)
    sina_posterior = np.sin(angle_posterior * np.pi / 180.0)
    temp_cos = np.percentile(cosa_posterior, 50.0)
    temp_sin = np.percentile(sina_posterior, 50.0)
    median_angle = (np.arctan2(temp_sin, temp_cos) * 180.0 / np.pi) % 360.0
    lower_angle = np.percentile(angle_posterior, 16.0)
    upper_angle = np.percentile(angle_posterior, 84.0)

    return median_angle, lower_angle, upper_angle


# ==============================================================================


def get_data(emcee_file, temp_status):

    ## read data from hdf5 file
    completed_steps = 0
    # print ' reading file', emcee_file
    # f_read = h5py.File(emcee_file, "r")
    with h5py.File(emcee_file, mode="r", swmr=True) as f_read:
        data_names = [str(i) for i in list(f_read.keys())]
        parameter_names_emcee, parameter_boundaries = [], []
        chains = []
        acceptance_fraction, autocor_time, lnprobability = [], [], []

        if "parameter_names" in data_names:
            # parameter_names_emcee = np.array(f_read['parameter_names'], dtype='S15').astype(str)
            names_temp = f_read["parameter_names"][...]
            parameter_names_emcee = decode_list(names_temp)
            # parameter_names_emcee = names_temp

        if "boundaries" in data_names:
            parameter_boundaries = np.array(f_read["boundaries"][...], dtype=np.float64)
        # if ('final_parameters' in data_names):  final_parameters = np.array(f_read['final_parameters'], dtype=np.float64)
        if "chains" in data_names:
            chains = np.array(
                f_read["chains"], dtype=np.float64
            )  # shape (nwalkers, nruns, nfit)
            chains = np.swapaxes(chains, 1, 0)  # nruns x nwalkers x nfit
        if "acceptance_fraction" in data_names:
            acceptance_fraction = np.array(
                f_read["acceptance_fraction"], dtype=np.float64
            )
        if "autocor_time" in data_names:
            autocor_time = np.array(f_read["autocor_time"], dtype=np.float64)
        if "lnprobability" in data_names:
            lnprobability = np.array(
                f_read["lnprobability"], dtype=np.float64
            )  # nwalkers x nruns
            lnprobability = np.swapaxes(lnprobability, 1, 0)

        try:
            ln_err_const = f_read["lnprobability"].attrs["ln_err_const"]
        except:
            ln_err_const = 0.0
        if temp_status:
            completed_steps = f_read["chains"].attrs["completed_steps"]
        else:
            completed_steps = np.shape(chains)[1]
    # close hdf5 file
    # f_read.close()

    return (
        parameter_names_emcee,
        parameter_boundaries,
        chains,
        acceptance_fraction,
        autocor_time,
        lnprobability,
        ln_err_const,
        completed_steps,
    )


# ==============================================================================


def get_last_emcee_iteration(emcee_file, nwalkers):

    with h5py.File(emcee_file, "r", swmr=True) as f_read:
        chains = f_read["chains"][...]
        nw_c, nr_c, _ = np.shape(chains)
        done = True

        try:
            last_step = f_read["chains"].attrs["completed_steps"]
        except:
            last_step = nr_c
    # f_read.close()

    if nwalkers == nw_c:
        last_p0 = [chains[iw, last_step - 1, :].copy() for iw in range(0, nw_c)]
        done = True
    else:
        last_p0 = None
        done = False

    del chains

    return last_p0, nw_c, done


# ==============================================================================


def compute_autocor_time(chains, walkers_transposed=True):

    if walkers_transposed:

        _, nw, nfit = np.shape(chains)
        autocor_time = np.zeros((nfit), dtype=np.float64)
        for ifit in range(0, nfit):
            x = chains[:, :, ifit]
            # tau, mean, sigma = acor.acor(x)
            # temp_acor = np.mean(np.array([acor.acor(x[:, iw]) for iw in range(0, nw)]), axis=0)
            temp_acor = np.mean(
                np.array(
                    [
                        emcee.autocorr.integrated_time(x[:, iw], c=1)
                        for iw in range(0, nw)
                    ]
                ),
                axis=0,
            )
            autocor_time[ifit] = temp_acor[0]

    else:

        nw, _, nfit = chains.shape
        autocor_time = np.zeros((nfit), dtype=np.float64)
        for ifit in range(0, nfit):
            x = chains[:, :, ifit]
            # temp_acor = np.mean(np.array([acor.acor(x[iw, :]) for iw in range(0, nw)]), axis=1)
            temp_acor = np.mean(
                np.array(
                    [
                        emcee.autocorr.integrated_time(x[iw, :], c=1)
                        for iw in range(0, nw)
                    ]
                ),
                axis=1,
            )
            autocor_time[ifit] = temp_acor[0]

    return autocor_time


# ==============================================================================
def compute_acor_time(sampler, steps_done=None):

    if steps_done is None:
        actime_const = 100
    else:
        if steps_done < 100:
            actime_const = steps_done
        else:
            actime_const = 100

    _, _, nfit = np.shape(sampler.chain)
    try:
        acor_time = sampler.acor
    except:
        acor_time = np.zeros((nfit), dtype=np.float64) + actime_const

    return acor_time


# ==============================================================================


def get_emcee_parameters(chains, temp_status, nburnin_in, completed_steps):

    # determine nwalkers, nruns, nburnin, and nfit parameters
    # nwalkers = chains.shape[0]
    # nruns = chains.shape[1]
    # nwalkers, nruns, nfit = np.shape(chains)
    nruns, nwalkers, nfit = np.shape(chains)
    if temp_status:
        nruns = int(completed_steps)

    # select posterior chains, without burn in steps
    if nburnin_in < 0:
        nburnin = 0
    elif nburnin_in >= nruns or nburnin_in is None:
        nburnin = np.rint(0.5 * nruns).astype(int)
    else:
        nburnin = int(nburnin_in)
    nruns_sel = nruns - nburnin

    return nfit, nwalkers, nruns, nburnin, nruns_sel


# ==============================================================================


def print_memory_usage(array_values, output=None):

    print_both(
        " MEMORY USAGE: array_values = {:d} bytes = {:.2f} MBytes = {:.4f} GBytes".format(
            array_values.nbytes,
            array_values.nbytes / (1024.0**2),
            array_values.nbytes / (1024.0**3),
        ),
        output=output,
    )

    return


# ==============================================================================


def select_transpose_convert_chains(
    nfit,
    nwalkers,
    npost,
    nruns,
    nruns_sel,
    m_factor,
    parameter_names_emcee,
    parameter_boundaries_in,
    chains,
):

    # chain is transposed: needed to plot quicker chains for each walker: nruns vs value of parameter
    parameter_boundaries = parameter_boundaries_in.copy()
    chains_T = np.zeros((nruns, nwalkers, nfit))

    for ii in range(0, nfit):
        chains_T[:, :, ii] = chains[:, :nruns, ii].T  # transpose
        if parameter_names_emcee[ii][0] == "m" and parameter_names_emcee[ii][1] != "A":

            if "Ms" in parameter_names_emcee[ii]:
                m_conv = m_factor
            else:
                m_conv = np.float64(1.0)

            chains_T[:, :, ii] = chains_T[:, :, ii] * m_conv
            parameter_boundaries[ii, :] = parameter_boundaries[ii, :] * m_conv

    return chains_T, parameter_boundaries


# ==============================================================================


def posterior_back_to_msun(m_factor, parameter_names_emcee, flatchain_posterior_in):

    nfit = flatchain_posterior_in.shape[1]
    flatchain_posterior_out = flatchain_posterior_in.copy()
    for ifit in range(0, nfit):
        if (
            parameter_names_emcee[ifit][0] == "m"
            and parameter_names_emcee[ifit][1] != "A"
        ):
            flatchain_posterior_out[:, ifit] = (
                flatchain_posterior_out[:, ifit] / m_factor
            )

    return flatchain_posterior_out


# ==============================================================================


def prepare_plot_folder(full_path):

    plot_folder = prepare_emcee_plot_folder(full_path)

    return plot_folder


# ==============================================================================


def prepare_emcee_plot_folder(full_path):

    emcee_plots = os.path.join(full_path, "plots")
    os.makedirs(emcee_plots, exist_ok=True)

    return emcee_plots


# ==============================================================================


def computation_time(elapsed):

    elapsed_d = elapsed / 86400.0
    elapsed_h = (elapsed_d - int(elapsed_d)) * 24.0
    elapsed_m = (elapsed_h - int(elapsed_h)) * 60.0
    elapsed_s = (elapsed_m - int(elapsed_m)) * 60.0

    return int(elapsed_d), int(elapsed_h), int(elapsed_m), elapsed_s


# ==============================================================================


def get_pso_data(pso_file):

    with h5py.File(pso_file, "r", swmr=True) as of_pso:
        population = np.array(of_pso["population"], dtype=np.float64)
        population_fitness = np.array(of_pso["population_fitness"], dtype=np.float64)
        pso_parameters = np.array(of_pso["pso_parameters"], dtype=np.float64)
        pso_fitness = np.array(of_pso["pso_fitness"], dtype=np.float64)
        if "pso_best_evolution" in list(of_pso.keys()):
            pso_best_evolution = np.array(
                of_pso["pso_best_evolution"], dtype=np.float64
            )
        else:
            pso_best_evolution = False
        if "parameters_minmax" in list(of_pso.keys()):
            parameters_minmax = np.array(of_pso["parameters_minmax"], dtype=np.float64)
        else:
            parameters_minmax = False
        if "parameter_names" in list(of_pso.keys()):
            temp_names = np.array(of_pso["parameter_names"], dtype="S10")
            parameter_names = [pn.decode("utf-8") for pn in temp_names]
        else:
            parameter_names = False
    # of_pso.close()
    pop_shape = population.shape

    return (
        population,
        population_fitness,
        pso_parameters,
        pso_fitness,
        pso_best_evolution,
        parameters_minmax,
        parameter_names,
        pop_shape,
    )


def pso_to_emcee(nfit, nwalkers, pso_population, names_par, sqrt_par=True):

    # _, npop, niter = np.shape(pso_population)
    pop_flat = np.reshape(np.swapaxes(pso_population, 2, 1), newshape=(nfit, -1))
    p0trades = pop_flat[:, -nwalkers:]
    p0 = []
    if sqrt_par:
        for i_w in range(0, nwalkers):
            p0.append(e_to_sqrte_fitting(p0trades[:, i_w], names_par))
    else:
        p0 = p0trades

    return p0


# ==============================================================================


def compute_limits(vec_a, delta=0.05):

    a_min = np.min(np.array(vec_a))
    a_max = np.max(np.array(vec_a))
    da = np.abs(a_max - a_min)
    lim_min = a_min - da * delta
    lim_max = a_max + da * delta

    return lim_min, lim_max


# ==============================================================================


def thin_the_chains(
    use_thin,
    nburnin,
    nruns,
    chains,
    lnprobability,
    burnin_done=False,
    full_chains_thinned=False,
):

    _, _, nfit = np.shape(chains)
    print(" nburnin   = ", nburnin)
    print(" nruns     = ", nruns)
    print(" thin      = ", use_thin)
    nr, nc = np.shape(lnprobability)

    if not burnin_done:
        chains_posterior = chains[nburnin:nruns, :, :].copy()
        lnprobability_posterior = lnprobability[nburnin:nruns, :].copy()
    else:
        chains_posterior = chains[:nruns, :, :].copy()
        lnprobability_posterior = lnprobability[:nruns, :].copy()

    nr_post, nw_post, _ = np.shape(chains_posterior)

    if use_thin <= 1:  # use_thin <= 1

        fitting_posterior = (
            chains_posterior[:, :, :].copy().reshape((nr_post * nw_post, nfit))
        )
        print(" fitting_posterior not thinned shape = ", np.shape(fitting_posterior))
        lnprob_posterior = (
            lnprobability_posterior[:, :].copy().reshape((nr_post * nw_post))
        )
        print(" lnprob_posterior not thinned shape  = ", np.shape(lnprob_posterior))
        thin_steps = 0
        chains_full_thinned = chains[:nruns, :, :]
        lnprobability_full_thinned = lnprobability[:nruns, :]

    else:  # use_thin > 0
        thin_steps = use_thin

        sel_thin_steps = np.arange(0, nr_post + thin_steps, thin_steps).astype(int)
        if sel_thin_steps[-1] >= nr_post:
            sel_thin_steps[-1] = nr_post - 1
        n_thin = np.shape(sel_thin_steps)[0]

        print(" n_thin = ", n_thin)

        chains_posterior = chains_posterior[sel_thin_steps, :, :]
        # create a flat array of the posterior: from (nruns_sel, nwalkers, nfit) -> (n_thin * nwalkers, nfit)
        nr_thin, nw_thin, _ = np.shape(chains_posterior)
        fitting_posterior = (
            chains_posterior[:, :, :].copy().reshape((nr_thin * nw_thin, nfit))
        )
        print(" fitting_posterior thinned shape = ", np.shape(fitting_posterior))
        lnprobability_posterior = lnprobability_posterior[sel_thin_steps, :]
        lnprob_posterior = lnprobability_posterior.copy().reshape((nr_thin, nw_thin))
        print(" lnprob_posterior thinned shape  = ", np.shape(lnprob_posterior))

        if full_chains_thinned:
            tempsel = np.arange(0, nruns + thin_steps, thin_steps).astype(int)
            if tempsel[-1] >= nruns:
                tempsel[-1] = nruns - 1
            chains_full_thinned = chains[tempsel, :, :]
            lnprobability_full_thinned = lnprobability[tempsel, :]
        else:
            chains_full_thinned = chains[:nruns, :, :]
            lnprobability_full_thinned = lnprobability[:nruns, :]

    return (
        chains_posterior,
        fitting_posterior,
        lnprobability_posterior,
        lnprob_posterior,
        chains_full_thinned,
        lnprobability_full_thinned,
    )


# ==============================================================================


def get_sigmas(best_parameters, flatchain_posterior):

    sigmas_percentiles = [15.87, 2.28, 0.13, 84.13, 97.72, 99.87]
    sigma_perc68 = np.percentile(
        np.abs(flatchain_posterior - best_parameters),
        68.27,
        axis=0,
        # ,
    )
    # retrieve confidence intervals of the residual distribution
    sigma_confint = np.percentile(
        flatchain_posterior - best_parameters,
        sigmas_percentiles,
        axis=0,
        # ,
    )

    return sigma_perc68, sigma_confint


# ==============================================================================


def get_maxlnprob_parameters(lnprob_burnin, chains_T, flatchain_posterior):

    print(" -- shape(chains_T)      = ", np.shape(chains_T))
    print(" -- shape(lnprob_burnin) = ", np.shape(lnprob_burnin))
    maxlnprob_row, maxlnprob_col = get_max_indices(lnprob_burnin)
    maxlnprob = lnprob_burnin[maxlnprob_row, maxlnprob_col]
    # retrieve best parameters <-> best lnprobability
    maxlnprob_parameters = chains_T[maxlnprob_row, maxlnprob_col, :]
    # retrieve 1sigma as 68.27th percentile of the absolute residual distribution
    # retrieve confidence intervals of the residual distribution
    maxlnprob_perc68, maxlnprob_confint = get_sigmas(
        maxlnprob_parameters, flatchain_posterior
    )

    return maxlnprob, maxlnprob_parameters, maxlnprob_perc68, maxlnprob_confint


# ==============================================================================


def get_median_parameters(flatchain_posterior):

    median_parameters = np.percentile(
        flatchain_posterior, 50.0, axis=0, methodo="midpoint"
    )
    median_perc68, median_confint = get_sigmas(median_parameters, flatchain_posterior)

    return median_parameters, median_perc68, median_confint


# ==============================================================================


def get_parameters_median_fitness(
    nwalkers, npost, nruns, lnprobability, flatchain_posterior, ln_err_const
):

    nruns_sel = nruns - npost
    lnprob_burnin = lnprobability[:, npost:nruns]
    flat_lnprob = lnprob_burnin.T.reshape((nwalkers * nruns_sel))
    flat_fitness = -2.0 * (flat_lnprob - ln_err_const)
    n_med = int(nwalkers * nruns_sel * 0.5)
    id_med = np.argsort(flat_fitness)[n_med]
    # retrieve median fitness
    median_fitness = flat_fitness[id_med]
    # retrieve parameters at id_med
    medfit_parameters = flatchain_posterior[id_med, :]
    medfit_perc68, medfit_confint = get_sigmas(medfit_parameters, flatchain_posterior)

    return median_fitness, medfit_parameters, medfit_perc68, medfit_confint


# ==============================================================================


def compute_max_mean(data_vec, k):

    hist_counts, bin_edges = np.histogram(data_vec, bins=k)
    max_bin = np.argmax(np.array(hist_counts))
    if max_bin == 0 or max_bin == k:
        ext_bin = 0
    elif max_bin == 1 or max_bin == k - 1:
        ext_bin = 1
    else:
        ext_bin = 2

    data_max = data_vec[
        np.logical_and(
            data_vec >= bin_edges[max_bin - ext_bin],
            data_vec < bin_edges[max_bin + ext_bin],
        )
    ]
    if np.shape(data_max)[0] < 1:
        print("WARNING: n(data_max) < 1!")
        return np.nan, 0
        sys.stdout.flush()
    max_mean = np.mean(data_max)

    return max_mean, max_bin


# ==============================================================================


def get_mode_parameters_full(flatchain_posterior, k):

    _, nfit = np.shape(flatchain_posterior)
    mode_parameters = np.zeros((nfit))
    mode_bin = np.zeros((nfit)).astype(int)

    for i_fit in range(0, nfit):
        data_vec = flatchain_posterior[:, i_fit]
        mode_parameters[i_fit], mode_bin[i_fit] = compute_max_mean(data_vec, k)

    mode_perc68, mode_confint = get_sigmas(mode_parameters, flatchain_posterior)

    return mode_bin, mode_parameters, mode_perc68, mode_confint


# ==============================================================================


def get_mode_parameters(flatchain_posterior, k):

    _, nfit = np.shape(flatchain_posterior)
    mode_parameters = np.zeros((nfit))
    mode_bin = np.zeros((nfit)).astype(int)

    print("get_mode_parameters")
    sys.stdout.flush()

    for i_fit in range(0, nfit):
        print("i_fit = {0:d}".format(i_fit), end=" ")
        data_vec = flatchain_posterior[:, i_fit]
        print("computing mode and bin ... ", end=" ")
        mode_parameters[i_fit], mode_bin[i_fit] = compute_max_mean(data_vec, k)
        print("({:.5f} , {:d}) done".format(mode_parameters[i_fit], mode_bin[i_fit]))
        sys.stdout.flush()

    return mode_bin, mode_parameters


# ==============================================================================


def get_sample_list(sample_str, parameter_names):

    nfit = np.shape(parameter_names)[0]
    print(" input sample to select = ", sample_str)
    str_list = sample_str.strip().split()

    name_par = parameter_names[0]
    name_excluded = None
    if len(str_list) == 1:
        if str_list[0].lower() == "none":
            name_par = parameter_names[0]
            name_excluded = None
        else:
            for ifit in range(0, nfit):
                if str_list[0].lower() == parameter_names[ifit].lower():
                    name_par = parameter_names[ifit]
                    name_excluded = None
    else:
        name_excluded = []
        for ilist in range(0, len(str_list)):
            for ifit in range(0, nfit):
                if str_list[ilist].lower() == parameter_names[ifit].lower():
                    if ilist == 0:
                        name_par = parameter_names[ifit]
                    else:
                        name_excluded.append(parameter_names[ifit])

    return name_par, name_excluded


# ==============================================================================


def check_sample_parameters(
    parameter_names, sample_parameters, post_ci, name_excluded=None
):

    nfit = np.shape(sample_parameters)[0]

    check_sample = True
    for ifit in range(0, nfit):
        if name_excluded is None:
            if sample_parameters[ifit] < post_ci[0, ifit]:
                check_sample = False
                break
            elif sample_parameters[ifit] > post_ci[1, ifit]:
                check_sample = False
                break
        else:
            if parameter_names[ifit] not in name_excluded:
                if sample_parameters[ifit] < post_ci[0, ifit]:
                    check_sample = False
                    break
                elif sample_parameters[ifit] > post_ci[1, ifit]:
                    check_sample = False
                    break

    return check_sample


# ==============================================================================


def pick_sample_parameters(
    posterior, parameter_names, name_par=None, name_excluded=None, post_ci=None
):

    if name_par is not None:
        npost, nfit = np.shape(posterior)

        if post_ci is None:
            post_ci = compute_hdi_full(posterior).T[0:2, :]

        sel_par = 0
        for ipar in range(0, nfit):
            if name_par.lower() == parameter_names[ipar].lower():
                sel_par = ipar
                break
        # use median
        # get idx sorted of the selected parameter-posterior
        idx_posterior = np.argsort(posterior[:, sel_par])
        # define the number of testing sample given the values within credible intervals
        n_test = int(
            0.5
            * np.sum(
                np.logical_and(
                    posterior[:, sel_par] >= post_ci[0, sel_par],
                    posterior[:, sel_par] <= post_ci[1, sel_par],
                ).astype(int)
            )
        )
        # create the list with the sorted idx, starting from the 50th and moving of +1,-1 every time
        testing = []
        testing.append(0)
        for itest in range(0, n_test):
            testing.append(itest + 1)
            testing.append(-(itest + 1))
        n_testing = len(testing)
        sel_idx = [int(0.5 * npost) + testing[ii] for ii in range(0, n_testing)]

        for ii in range(0, n_testing):
            # w/ median
            idx = idx_posterior[sel_idx[ii]]

            sample_parameters = posterior[idx, :]
            check_sample = True
            check_sample = check_sample_parameters(
                parameter_names, sample_parameters, post_ci, name_excluded
            )
            if check_sample:
                print(" found good sample parameters at index: ", idx)
                print(sample_parameters)
                return sample_parameters, idx

        return None, None


# =============================================================================

# 1) select parameters and lgllhd within ci
def select_within_all_ci(posterior, post_ci, lnprobability):

    npost, nfit = np.shape(posterior)
    use_ci = post_ci

    # use_ci should have: nfit x nci,
    # where nci:
    # -1sigma(0) +1sigma(1) -2sigma(2) +2sigma(3) -3sigma(4) +3sigma(5)
    ok_sel = np.ones((npost)).astype(bool)
    for ifit in range(0, nfit):
        ok_temp = np.logical_and(
            posterior[:, ifit] >= use_ci[ifit, 0], posterior[:, ifit] <= use_ci[ifit, 1]
        )
        ok_sel = np.logical_and(ok_sel, ok_temp)

    post_sel = posterior[ok_sel, :]
    lnprob_sel = lnprobability[ok_sel]

    return post_sel, lnprob_sel


# ==============================================================================


def select_maxlglhd_with_hdi(posterior, post_ci, lnprobability, return_idmax=False):

    npost, _ = np.shape(posterior)
    post_sel, lnprob_sel = select_within_all_ci(
        posterior, post_ci, lnprobability  # .reshape((npost))
    )
    idmax = np.argmax(lnprob_sel)

    sample_par = post_sel[idmax, :]
    sample_lgllhd = lnprob_sel[idmax]

    if return_idmax:
        return sample_par, sample_lgllhd, idmax

    return sample_par, sample_lgllhd


# ==============================================================================

# 2) determine the median lgllhd lgllhd_med of the selected pars
# 3) sort by lgllhd-lgllhd_med and take the first set of pars
def get_sample_by_sorted_lgllhd(posterior, lnprobability, post_ci=None):

    if post_ci is None:
        post_ci = compute_hdi_full(posterior)

    npost, _ = np.shape(posterior)

    use_lnprob = lnprobability.reshape((npost))

    post_sel, lnprob_sel = select_within_all_ci(posterior, post_ci, use_lnprob)

    lgllhd_med = np.percentile(
        use_lnprob,
        50.0,
    )
    abs_dlg = np.abs(lnprob_sel - lgllhd_med)

    idx_sample = np.argmin(abs_dlg)
    sample_parameters = post_sel[idx_sample, :]
    sample_lgllhd = lnprob_sel[idx_sample]

    return sample_parameters, sample_lgllhd


# ==============================================================================


def get_sample_by_par_and_lgllhd(
    posterior, lnprobability, parameter_names, post_ci=None, name_par=None
):

    if post_ci is None:
        post_ci = compute_hdi_full(posterior)

    npost, nfit = np.shape(posterior)

    if name_par is None:
        sel_par = 0
    else:
        for ifit in range(0, nfit):
            if name_par == parameter_names[ifit].strip():
                sel_par = ifit
                break

    use_lnprob = lnprobability.reshape((npost))

    post_sel, lnprob_sel = select_within_all_ci(posterior, post_ci, use_lnprob)

    lgllhd_med = np.percentile(
        use_lnprob,
        50.0,
    )
    abs_dlg = np.abs(lnprob_sel - lgllhd_med)
    lgllhd_mad = np.percentile(
        abs_dlg,
        50.0,
    )
    lgllhd_min = lgllhd_med - lgllhd_mad
    lgllhd_max = lgllhd_med + lgllhd_mad

    par_posterior = posterior[:, sel_par]
    par_med = np.percentile(
        par_posterior,
        50.0,
    )

    abs_dpar = np.abs(post_sel[:, sel_par] - par_med)
    ids_par = np.argsort(abs_dpar)
    lgllhd_check = np.logical_and(lnprob_sel >= lgllhd_min, lnprob_sel <= lgllhd_max)
    sample_parameters, sample_lgllhd = None, None
    for ii in range(0, np.shape(ids_par)[0]):
        idx = ids_par[ii]
        if lgllhd_check[idx]:
            sample_parameters = post_sel[idx, :]
            sample_lgllhd = lnprob_sel[idx]
            break

    return sample_parameters, sample_lgllhd


# ==============================================================================


def take_n_samples(posterior, lnprob=None, post_ci=None, n_samples=100):

    npost, nfit = np.shape(posterior)

    print("take_n_samples")
    print("shape posterior ", np.shape(posterior))
    print("shape lbprob    ", np.shape(lnprob))
    print("shape post_ci   ", np.shape(post_ci))

    if lnprob is not None:
        idx_post = np.argsort(lnprob)[::-1]  # sort: descending
    else:
        idx_post = np.arange(0, npost).astype(int)

    post_sort = posterior[
        idx_post, :
    ]  # sort properly so then selection matches indexes

    sel_within_ci = np.ones((npost)).astype(bool)
    if post_ci is not None:
        # set selection withing hdi (ci) to True
        # sel_within_ci = np.ones((npost)).astype(bool)
        # loop on parameters and check if within hdi (ci)
        for ifit in range(0, nfit):
            sel_par = np.logical_and(
                post_sort[:, ifit] >= post_ci[0, ifit],
                post_sort[:, ifit] <= post_ci[1, ifit],
            )
            sel_within_ci = np.logical_and(sel_within_ci, sel_par)

    #     idx_sample = np.random.choice(idx_post[sel_within_ci], n_samples, replace=False)
    # else:
    #     idx_sample = np.random.choice(idx_post, n_samples, replace=False)
    nsel = np.sum(sel_within_ci)
    if n_samples > nsel:
        nmissing = n_samples - nsel
        not_sel = sel_within_ci == False
        not_sel_idx = np.random.choice(idx_post[not_sel], nmissing, replace=False)
        sel_within_ci[not_sel_idx] = True
    # else:
    #     sel = sel_within_ci
    idx_sample = np.random.choice(idx_post[sel_within_ci], n_samples, replace=False)

    sample_parameters = posterior[idx_sample, :]

    return sample_parameters


# =============================================================================


def print_parameters_nolog(parameter_names, parameters, perc68, confint, par_type):

    sigmas_percentiles = [15.87, 2.28, 0.13, 84.13, 97.72, 99.87]
    nsigmas = ["-1", "-2", "-3", "+1", "+2", "+3"]

    f_w = 23

    print()
    line_68 = "+/-1sigma(68.27|res|)".rjust(f_w)
    line_s = "%s" % (
        " ".join(
            [
                "%2ssigma(%5.2fres)".rjust(f_w) % (nsigmas[j], sigmas_percentiles[j])
                for j in range(len(sigmas_percentiles))
            ]
        )
    )
    header = "%s %23s %s %s" % ("#".ljust(10), par_type, line_68, line_s)
    print(header)

    for i_fit in range(0, parameters.shape[0]):
        sigma_line = "".join(
            [
                "%23.16f" % (confint[j_perc, i_fit])
                for j_perc in range(len(sigmas_percentiles))
            ]
        )
        line = "%10s %23.16f %23.16f %s" % (
            parameter_names[i_fit],
            parameters[i_fit],
            perc68[i_fit],
            sigma_line,
        )
        print(line)
    print()

    return


# ==============================================================================


def print_parameters_logtxt(
    of_run, parameter_names, parameters, perc68, confint, par_type
):

    sigmas_percentiles = [15.87, 2.28, 0.13, 84.13, 97.72, 99.87]
    nsigmas = ["-1", "-2", "-3", "+1", "+2", "+3"]

    f_w = 23

    print()
    of_run.write("\n")
    line_68 = "+/-1sigma(68.27|res|)".rjust(f_w)
    line_s = "%s" % (
        " ".join(
            [
                "%2ssigma(%5.2fres)".rjust(f_w) % (nsigmas[j], sigmas_percentiles[j])
                for j in range(len(sigmas_percentiles))
            ]
        )
    )
    header = "%s %23s %s %s" % ("#".ljust(10), par_type, line_68, line_s)
    print(header)
    of_run.write(header + "\n")

    for i_fit in range(0, parameters.shape[0]):
        sigma_line = "".join(
            [
                "%23.16f" % (confint[j_perc, i_fit])
                for j_perc in range(len(sigmas_percentiles))
            ]
        )
        line = "%10s %23.16f %23.16f %s" % (
            parameter_names[i_fit],
            parameters[i_fit],
            perc68[i_fit],
            sigma_line,
        )
        print(line)
        of_run.write(line + "\n")
    print()
    of_run.write("\n")

    return


# ==============================================================================


def print_parameters_logger(
    logger, parameter_names, parameters, perc68, confint, par_type
):

    sigmas_percentiles = [15.87, 2.28, 0.13, 84.13, 97.72, 99.87]
    nsigmas = ["-1", "-2", "-3", "+1", "+2", "+3"]

    f_w = 23

    logger.info("")
    line_68 = "+/-1sigma(68.27|res|)".rjust(f_w)
    line_s = "%s" % (
        " ".join(
            [
                "%2ssigma(%5.2fres)".rjust(f_w) % (nsigmas[j], sigmas_percentiles[j])
                for j in range(len(sigmas_percentiles))
            ]
        )
    )
    header = "%s %23s %s %s" % ("#".ljust(10), par_type, line_68, line_s)
    logger.info(header)

    for i_fit in range(0, parameters.shape[0]):
        sigma_line = "".join(
            [
                "%23.16f" % (confint[j_perc, i_fit])
                for j_perc in range(len(sigmas_percentiles))
            ]
        )
        line = "%10s %23.16f %23.16f %s" % (
            parameter_names[i_fit],
            parameters[i_fit],
            perc68[i_fit],
            sigma_line,
        )
        logger.info(line)
    logger.info("")

    return


# =============================================================================


def GelmanRubin(chains_T):

    n, M = np.shape(chains_T)

    theta_m = np.mean(chains_T, axis=0)

    B_n = np.var(theta_m, ddof=1)

    arg_W = np.var(chains_T, ddof=1, axis=0)
    W = np.mean(arg_W)

    n_frac = (n - 1) / n
    var_plus = n_frac * W + B_n
    Var = var_plus + (B_n / M)

    Rc = np.sqrt(Var / W)

    return Rc


# ==============================================================================

# Geweke 1992 test 1
def geweke_test(chains_T, start_frac=0.05, n_sel_steps=5):

    n_steps, n_chains = chains_T.shape
    half = int(0.5 * n_steps)

    int_b = chains_T[half:, :]
    mean_b = np.mean(int_b, axis=0)
    var_b = np.var(int_b, axis=0, ddof=1)

    start_step = int(start_frac * n_steps)
    sel_a = np.linspace(
        start=start_step, stop=half, num=n_sel_steps, endpoint=False, dtype=np.int
    )

    zscore = np.zeros((n_sel_steps, n_chains))
    for i_a in range(0, n_sel_steps):
        int_a = chains_T[sel_a[i_a] : half, :]
        mean_a = np.mean(int_a, axis=0)
        var_a = np.var(int_a, axis=0, ddof=1)
        z_temp = (mean_a - mean_b) / np.sqrt(var_a + var_b)
        zscore[i_a, :] = z_temp

    return sel_a, zscore


# =============================================================================

# ==============================================================================


def get_units(names, mass_unit):

    n_names = len(names)
    units_par = []

    for i in range(0, n_names):

        if str(names[i])[0] == "m":
            if "Ms" in names[i]:
                units_par.append("(M_sun/M_star)")
            elif "mA" in names[i]:
                units_par.append("(deg)")
            elif "m1" == names[i]:
                units_par.append("(M_sun)")
            else:
                units_par.append("(%s)" % (mass_unit))

        elif "R" in str(names[i]):
            units_par.append("(R_sun)")

        elif "P" in str(names[i]):
            units_par.append("(days)")

        elif str(names[i])[0] == "e":  # the same for e, ecosw, esinw
            units_par.append(" ")

        elif str(names[i])[0] == "s":  # sqrtecosw, sqrtesinw
            units_par.append(" ")

        elif str(names[i])[0] == "w":
            units_par.append("(deg)")

        elif "lambda" in str(names[i]):
            units_par.append("(deg)")

        elif "cosi" in names[i]:
            units_par.append(" ")
        elif str(names[i])[0] == "i":
            units_par.append("(deg)")

        elif str(names[i])[0:2] == "lN":
            units_par.append("(deg)")

        elif "jitter" in str(names[i]):
            units_par.append("(m/s)")
        else:
            units_par.append(" ")

    return units_par


# ==============================================================================


def read_initial_elements(fpath, idsim, lmf=0):

    kel_file = os.path.join(fpath, "{:d}_{:d}_initialElements.dat".format(idsim, lmf))
    try:
        kep_elem = np.genfromtxt(
            kel_file
        )  # (NB-1) x (M_Msun R_Rsun P_d a_AU ecc w_deg mA_deg inc_deg lN_deg)
    except:
        print(" KEPLERIAN ELEMENTS FILE NOT FOUND {:s}".format(kel_file))
        kel_file, kep_elem = None, None

    return kel_file, kep_elem


# ==============================================================================


def get_case(id_body, fit_body):

    fit_n = np.arange(1, 11, 1)
    id_fit = [fit_n[j] for j in range(len(fit_body)) if (fit_body[j] == "1")]
    id_all = [8 * (id_body - 1) + 2 + id_fit[j] for j in range(len(id_fit))]

    if 6 not in id_fit:  # not fitting lambda
        case = [0]

    else:  # (6 in id_fit) # fitting lambda
        if all(x in id_fit for x in [4, 5, 6, 7, 8]):  # fit secosw, sesinw, lambda, lN
            case = [4]
        elif all(x in id_fit for x in [4, 5, 6]):  # fit secosw, sesinw, lambda
            case = [2]
        elif all(x in id_fit for x in [6, 7, 8]):  # fit lambda, cicoslN, cisinlN
            case = [3]
        else:
            case = [1]  # fit lambda & w || lambda & lN || lamda & w & lN

    return id_fit, id_all, case


# ==============================================================================


def get_fitted(full_path):

    of = open(os.path.join(full_path, "bodies.lst"), "r")
    lines = of.readlines()
    of.close()
    NB = len(lines)
    bodies_file = []
    fit_string = ""
    fit_list = []
    for line in lines:
        fname = line.strip().split()[0]
        bodies_file.append(fname)
        temp_string = line.strip().split(fname)[1].split("#")[0].strip()
        if "F" in temp_string:
            temp_string = temp_string.split("F")[0].strip()
        elif "f" in temp_string:
            temp_string = temp_string.split("f")[0].strip()
        elif "T" in temp_string:
            temp_string = temp_string.split("T")[0].strip()
        elif "t" in temp_string:
            temp_string = temp_string.split("t")[0].strip()

        fit_string = "%s %s" % (fit_string, temp_string)
        fit_list.append(temp_string.split())

    nfit = np.sum(np.array(fit_string.split(), dtype=np.int))

    case = []
    id_fit = []
    id_all = []
    nfit_list = []
    cols_list = []
    nfit_cnt = 0
    for j in range(1, NB + 1):
        id_fit_j, id_all_j, case_j = get_case(j, fit_list[j - 1])
        id_fit.append(id_fit_j)
        id_all.append(id_all_j)
        case.append(case_j)
        nfit_j = len(id_fit_j)
        nfit_list.append(nfit_j)
        cols_list.append([cc for cc in range(nfit_cnt, nfit_cnt + nfit_j)])
        nfit_cnt += nfit_j

    return nfit, NB, bodies_file, id_fit, id_all, nfit_list, cols_list, case


# ==============================================================================


def compute_intervals(flatchain, parameters, percentiles):

    sigma_par = np.percentile(
        np.subtract(flatchain, parameters),
        percentiles,
        axis=0,
        # ,
    )  # (n_percentile x nfit)
    sigma_par[0] = np.percentile(
        np.abs(np.subtract(flatchain, parameters)),
        percentiles[0],
        axis=0,
        # ,
    )  # 68.27th
    sigma_par[1] = np.percentile(
        np.abs(np.subtract(flatchain, parameters)),
        percentiles[1],
        axis=0,
        # ,
    )  # MAD

    return sigma_par


# ==============================================================================


def compute_hdi_full(flatchains):

    # alpha=[0.3173, 0.0456, 0.0026]
    _, npar = np.shape(flatchains)
    hdi_full = [calculate_hpd(flatchains[:, ipar]) for ipar in range(npar)]

    hdi_l = [
        [
            hdi_full[ipar][0][0],
            hdi_full[ipar][0][1],  # -1sigma +1sigma
            hdi_full[ipar][1][0],
            hdi_full[ipar][1][1],  # -2sigma +2sigma
            hdi_full[ipar][2][0],
            hdi_full[ipar][2][1],  # -3sigma +3sigma
        ]
        for ipar in range(npar)
    ]

    hdi_a = np.array(hdi_l, dtype=np.float64)

    hdi = np.reshape(hdi_a, newshape=((npar, -1)))

    return hdi


# ==============================================================================
def hdi_to_sigma(par, hdi):

    # print("par shape = {}".format(np.shape(par)))
    # print("hdi shape = {}".format(np.shape(hdi)))
    # sigma = np.zeros_like(hdi)
    # sigma[:, 0] = hdi[:, 0] - par
    # sigma[:, 1] = hdi[:, 1] - par
    sigma = hdi - par

    return sigma


def posterior_to_rms_mad(posterior, parameters):

    npar = len(parameters)
    ares = np.abs(posterior - parameters)
    mad = np.percentile(
        ares,
        50.00,
        # ,
        axis=0,
    ).reshape(npar, 1)
    rms = np.percentile(
        ares,
        68.27,
        # ,
        axis=0,
    ).reshape(npar, 1)

    return mad, rms


# ==============================================================================


def compute_sigma_hdi(flatchains, parameters):

    # alpha=[0.3173, 0.0456, 0.0026]

    _, npar = np.shape(flatchains)
    hdi_full = [calculate_hpd(flatchains[:, ipar]) for ipar in range(npar)]

    sigma_par = np.array(
        [
            [
                hdi_full[ipar][0][0] - parameters[ipar],
                hdi_full[ipar][0][1] - parameters[ipar],
                hdi_full[ipar][1][0] - parameters[ipar],
                hdi_full[ipar][1][1] - parameters[ipar],
                hdi_full[ipar][2][0] - parameters[ipar],
                hdi_full[ipar][2][1] - parameters[ipar],
            ]
            for ipar in range(npar)
        ]
    ).reshape((npar, -1))

    delta_flat_T = np.abs(
        [flatchains[:, ipar] - parameters[ipar] for ipar in range(npar)]
    ).T

    sigma_r68p = np.percentile(
        delta_flat_T,
        68.27,
        axis=0,
    )  # 68.27th

    sigma_mad = np.percentile(
        delta_flat_T,
        50.0,
        axis=0,
    )  # MAD

    sigma_par = np.column_stack((sigma_r68p, sigma_mad, sigma_par))

    return sigma_par


# ==============================================================================


def get_good_distribution(posterior_scale, posterior_mod, type_out=False, debug=False):

    par_scale = np.median(posterior_scale)
    p68_scale = np.percentile(
        np.abs(posterior_scale - par_scale),
        68.27,
    )
    std_scale = np.std(posterior_scale, ddof=1)
    s_scale = max(p68_scale, std_scale)

    par_mod = np.median(posterior_mod)
    p68_mod = np.percentile(
        np.abs(posterior_mod - par_mod),
        68.27,
    )
    std_mod = np.std(posterior_mod, ddof=1)
    s_mod = max(p68_mod, std_mod)

    if debug:
        print(
            "scaled [{:7.2f}, {:7.2f}] std: {} ==> {}".format(
                np.min(posterior_scale),
                np.max(posterior_scale),
                (p68_scale, std_scale),
                s_scale,
            )
        )
        print(
            "mod    [{:7.2f}, {:7.2f}] std: {} ==> {}".format(
                np.min(posterior_mod), np.max(posterior_mod), (p68_mod, std_mod), s_mod
            )
        )

    if s_scale < s_mod:
        if debug:
            print("Recenter as SCALE")
        par_out = par_scale
        posterior_out = posterior_scale
        par_type = "scale"
    else:
        if debug:
            print("Recenter as MOD")
        par_out = par_mod
        posterior_out = posterior_mod
        par_type = "mod"

    if type_out:
        return par_out, posterior_out, par_type
    else:
        return par_out, posterior_out


# ==============================================================================


def get_proper_posterior_correlated(posterior, col=None, type_out=False, debug=False):

    if col is not None:
        posterior_scale = np.arctan2(posterior[:, col + 1], posterior[:, col]) * rad2deg
        posterior_mod = posterior_scale % 360.0
    else:
        posterior_scale = (
            np.arctan2(np.sin(posterior * deg2rad), np.cos(posterior * deg2rad))
            * rad2deg
        )
        posterior_mod = posterior_scale % 360.0

    # pass type_out=True, but return only if requested
    par_out, post_out, par_type = get_good_distribution(
        posterior_scale, posterior_mod, type_out=True, debug=debug
    )
    if type_out:
        return par_out, post_out, par_type
    else:
        return par_out, post_out


# ==============================================================================
# ==============================================================================
def get_arctan_angle(alpha_deg):

    alpha_rad = alpha_deg * cst.deg2rad
    cosa = np.cos(alpha_rad)
    sina = np.sin(alpha_rad)
    beta_deg = np.arctan2(sina, cosa) * cst.rad2deg

    return beta_deg

# ==============================================================================
# ==============================================================================
# fix lambda only in flatchain posterior
def fix_lambda(flatchain_post, names_par):

    _, nfit = np.shape(flatchain_post)
    flatchain_post_out = flatchain_post.copy()
    for ifit in range(0, nfit):
        if "lambda" in names_par[ifit]:
            _, xout = get_proper_posterior_correlated(flatchain_post[:, ifit])
            flatchain_post_out[:, ifit] = xout.copy()

    return flatchain_post_out


# ==============================================================================
# ==============================================================================
# recenter single angle distribution
def recenter_angle_distribution(alpha_deg, debug=False, type_out=False):

    beta_deg = get_arctan_angle(alpha_deg)

    _, recentered, rec_type = get_good_distribution(beta_deg, alpha_deg, type_out=True, debug=debug)
    if type_out:
        return recentered, rec_type
    else:
        return recentered


# ==============================================================================
# ==============================================================================


# def get_trigonometric_posterior(id_NB, id_l, col, idpar, posterior):
#     # ecosw,esinw // sqrt(e)cosw,sqrt(e)sinw -> (e,w) || No lambda

#     cosboot = posterior[:, col]
#     sinboot = posterior[:, col + 1]

#     sum_square = (cosboot * cosboot) + (sinboot * sinboot)

#     if (idpar[id_l][0:6] == "sqrtec") or (idpar[id_l][0:2] == "se"):
#         # ecc = (sqrt(e)*cosw)^2 + (sqrt(e)*sinw)^2
#         aboot = sum_square
#     else:
#         # ecc = sqrt(ecosw^2+esinw^2)
#         aboot = np.sqrt(sum_square)

#     _, bboot = get_proper_posterior_correlated(posterior, col)

#     if "e" in idpar[id_l]:
#         sel_ezero = aboot <= eps64bit
#         bboot[sel_ezero] = 90.0

#     if "ecosw" in idpar[id_l]:
#         id_aa = "e"
#         id_bb = "w"
#     # else:
#     #     id_aa = "i"
#     #     id_bb = "lN"

#     id_a0 = "{:s}{}".format(id_aa, id_NB + 1)
#     id_b0 = "{:s}{}".format(id_bb, id_NB + 1)

#     return id_a0, aboot, id_b0, bboot


# ==============================================================================


# def get_trigonometric_parameters(
#     id_NB, id_l, col, idpar, parameters
# ):  # ecosw,esinw // sqrt(e)cosw,sqrt(e)sinw -> (e,w) && icoslN,isinlN -> (i,lN) || No lambda

#     cosp = parameters[col]
#     sinp = parameters[col + 1]
#     sum_square = (cosp * cosp) + (sinp * sinp)

#     if (idpar[id_l][0:6] == "sqrtec") or (idpar[id_l][0:2] == "se"):
#         # ecc = (sqrt(e)*cosw)^2 + (sqrt(e)*sinw)^2
#         a_par = sum_square
#     elif ("cosi" in idpar[id_l]) or ("ci" in idpar[id_l]):
#         # i = arccos(sqrt(cosicoslN^2+cosisinlN^2))
#         a_par = np.arccos(np.sqrt(sum_square)) * cst.rad2deg
#     else:
#         # ecc = sqrt(ecosw^2+esinw^2) || inc = sqrt(icoslN^2+isinlN^2)
#         a_par = np.sqrt(sum_square)

#     b_par = np.arctan2(sinp, cosp) * rad2deg

#     if "e" in idpar[id_l]:
#         if a_par <= eps64bit:
#             b_par = 90.0

#     if "ecosw" in idpar[id_l]:
#         id_aa = "e"
#         id_bb = "w"
#     else:
#         id_aa = "i"
#         id_bb = "lN"

#     id_a0 = "{:s}{}".format(id_aa, id_NB + 1)
#     id_b0 = "{:s}{}".format(id_bb, id_NB + 1)

#     return id_a0, a_par, id_b0, b_par


# ==============================================================================


# def derived_posterior_case0(
#     idpar_NB, id_fit_NB, i_NB, cols, posterior
# ):  # not fitting lambda

#     name_der = []
#     der_posterior = []

#     nlist = len(id_fit_NB)
#     for i_l in range(nlist):
#         # this condition work for sqrt(e)cos/sinw, ecos/sinw, cosicos/sinlN, icos/sinlN
#         if "ecosw" in idpar_NB[i_l] or "icoslN" in idpar_NB[i_l]:
#             id_a0, aboot, id_b0, bboot = get_trigonometric_posterior(
#                 i_NB, i_l, cols[i_l], idpar_NB, posterior
#             )
#             name_der.append(id_a0)
#             name_der.append(id_b0)
#             der_posterior.append(aboot)
#             der_posterior.append(bboot)

#     return name_der, der_posterior


# ==============================================================================


# def derived_parameters_case0(
#     idpar_NB, id_fit_NB, i_NB, cols, parameters
# ):  # not fitting lambda

#     name_der = []
#     der_par = []

#     nlist = len(id_fit_NB)
#     for i_l in range(nlist):
#         # this condition work for sqrt(e)cos/sinw, ecos/sinw, cosicos/sinlN, icos/sinlN
#         if "ecosw" in idpar_NB[i_l] or "icoslN" in idpar_NB[i_l]:
#             id_a0, a_par, id_b0, b_par = get_trigonometric_parameters(
#                 i_NB, i_l, cols[i_l], idpar_NB, parameters
#             )
#             name_der.append(id_a0)
#             name_der.append(id_b0)
#             der_par.append(a_par)
#             der_par.append(b_par)

#     return name_der, der_par


# ==============================================================================


# def derived_posterior_case1(
#     idpar_NB, id_fit_NB, i_NB, cols, kep_elem, posterior
# ):  # fitting lambda || lambda & w || lambda & lN || lamda & w & lN

#     name_der = []
#     der_posterior = []

#     nlist = len(id_fit_NB)

#     for i_l in range(0, nlist):
#         if id_fit_NB[i_l] == 6:  # lambda
#             lambda_post = posterior[:, cols[i_l]]
#             if i_l > 0:

#                 if i_l - 1 == 5:  # id 5 == argp
#                     argp = posterior[:, cols[i_l - 1]]  # argp as posterior
#                 else:
#                     argp = kep_elem[5]  # argp is fixed

#                 if (
#                     i_l < nlist - 1 and i_l + 1 == nlist - 1
#                 ):  # id 6 is not the last --> 8 is the next (longn)
#                     longn = posterior[:, cols[i_l + 1]]  # long node as posterior
#                 else:
#                     # if(i_l == nlist-1):
#                     longn = kep_elem[8]  # long of node is fixed

#             else:
#                 argp = kep_elem[5]  # argp is fixed
#                 longn = kep_elem[8]  # longn is fixed

#             aboot = (lambda_post - argp - longn) % 360.0
#             _, aboot = get_proper_posterior_correlated(aboot)
#             name_der.append("{:s}{}".format("mA", i_NB + 1))
#             der_posterior.append(aboot)

#     return name_der, der_posterior


# ==============================================================================


# def derived_parameters_case1(
#     idpar_NB, id_fit_NB, i_NB, cols, kep_elem, parameters
# ):  # fitting lambda || lambda & w || lambda & lN || lamda & w & lN
#     name_der = []
#     der_par = []

#     nlist = len(id_fit_NB)

#     for i_l in range(0, nlist):
#         if id_fit_NB[i_l] == 6:  # lambda
#             lambda_par = parameters[cols[i_l]]
#             if i_l > 0:

#                 if i_l - 1 == 5:  # id 5 == argp
#                     argp = parameters[cols[i_l - 1]]  # argp as parameter
#                 else:
#                     argp = kep_elem[5]  # argp is fixed

#                 if (
#                     i_l < nlist - 1 and i_l + 1 == nlist - 1
#                 ):  # id 6 is not the last --> 8 is the next (longn)
#                     longn = parameters[cols[i_l + 1]]  # long node as posterior
#                 else:
#                     # if(i_l == nlist-1):
#                     longn = kep_elem[8]  # long of node is fixed

#             else:
#                 argp = kep_elem[5]  # argp is fixed
#                 longn = kep_elem[8]  # longn is fixed

#             a_par = (lambda_par - argp - longn) % 360.0
#             name_der.append("{:s}{}".format("mA", i_NB + 1))
#             der_par.append(a_par)

#     return name_der, der_par


# ==============================================================================


# def derived_posterior_case2(
#     idpar_NB, id_fit_NB, i_NB, cols, kep_elem, posterior
# ):  # fit ecosw//sqrt(e)cosw, esinw//sqrt(e)sinw, lambda

#     name_der = []
#     der_posterior = []

#     nlist = len(id_fit_NB)
#     for i_l in range(0, nlist):

#         if id_fit_NB[i_l] == 4:
#             id_e0, eboot, id_w0, wboot = get_trigonometric_posterior(
#                 i_NB, i_l, cols[i_l], idpar_NB, posterior
#             )
#             name_der.append(id_e0)
#             name_der.append(id_w0)
#             der_posterior.append(eboot)
#             der_posterior.append(wboot)

#             if idpar_NB[i_l + 2][0:2] != "mA":
#                 mAboot = (
#                     posterior[:, cols[i_l + 2]] - wboot - kep_elem[8]
#                 ) % 360.0  # mA = lambda - w - lN
#                 _, mAboot = get_proper_posterior_correlated(mAboot)
#                 name_der.append("{:s}{}".format("mA", i_NB + 1))
#                 der_posterior.append(mAboot)

#     return name_der, der_posterior


# ==============================================================================


# def derived_parameters_case2(
#     idpar_NB, id_fit_NB, i_NB, cols, kep_elem, parameters
# ):  # fit ecosw//sqrt(e)cosw, esinw//sqrt(e)sinw, lambda

#     name_der = []
#     der_par = []

#     nlist = len(id_fit_NB)
#     for i_l in range(0, nlist):

#         if id_fit_NB[i_l] == 4:
#             id_e0, e_par, id_w0, w_par = get_trigonometric_parameters(
#                 i_NB, i_l, cols[i_l], idpar_NB, parameters
#             )
#             name_der.append(id_e0)
#             name_der.append(id_w0)
#             der_par.append(e_par)
#             der_par.append(w_par)

#             if idpar_NB[i_l + 2][0:2] != "mA":
#                 mA_par = (
#                     parameters[cols[i_l + 2]] - w_par - kep_elem[8]
#                 ) % 360.0  # mA = lambda - w - lN
#                 name_der.append("{:s}{}".format("mA", i_NB + 1))
#                 der_par.append(mA_par)
#     return name_der, der_par


# ==============================================================================


# def derived_posterior_case3(
#     idpar_NB, id_fit_NB, i_NB, cols, kep_elem, posterior
# ):  # fit lambda, icoslN, isinlN

#     name_der = []
#     der_posterior = []

#     nlist = len(id_fit_NB)
#     for i_l in range(0, nlist):

#         if id_fit_NB[i_l] == 7:
#             id_i0, iboot, id_lN0, lNboot = get_trigonometric_posterior(
#                 i_NB, i_l, cols[i_l], idpar_NB, posterior
#             )
#             name_der.append(id_i0)
#             name_der.append(id_lN0)
#             der_posterior.append(iboot)
#             der_posterior.append(lNboot)

#             if idpar_NB[i_l + 2][0:2] != "mA":
#                 mAboot = (
#                     posterior[:, cols[i_l - 1]] - kep_elem[5] - lNboot
#                 ) % 360.0  # mA = lambda - w - lN
#                 _, mAboot = get_proper_posterior_correlated(mAboot)
#                 name_der.append("{:s}{}".format("mA", i_NB + 1))
#                 der_posterior.append(mAboot)

#     return name_der, der_posterior


# ==============================================================================


# def derived_parameters_case3(
#     idpar_NB, id_fit_NB, i_NB, cols, kep_elem, parameters
# ):  # fit lambda, cosicoslN, cosisinlN

#     name_der = []
#     der_par = []

#     nlist = len(id_fit_NB)
#     for i_l in range(0, nlist):

#         if id_fit_NB[i_l] == 7:
#             id_i0, i_par, id_lN0, lN_par = get_trigonometric_parameters(
#                 i_NB, i_l, cols[i_l], idpar_NB, parameters
#             )
#             name_der.append(id_i0)
#             name_der.append(id_lN0)
#             der_par.append(i_par)
#             der_par.append(lN_par)

#             if idpar_NB[i_l + 2][0:2] != "mA":
#                 mA_par = (
#                     parameters[cols[i_l - 1]] - kep_elem[5] - lN_par
#                 ) % 360.0  # mA = lambda - w - lN
#                 name_der.append("{:s}{}".format("mA", i_NB + 1))
#                 der_par.append(mA_par)

#     return name_der, der_par


# ==============================================================================


# def derived_posterior_case4(
#     idpar_NB, id_fit_NB, i_NB, cols, kep_elem, posterior
# ):  # fit ecosw//sqrt(e)cosw, esinw//sqrt(e)sinw, lambda, cosicoslN, cosisinlN

#     name_der = []
#     der_posterior = []

#     nlist = len(id_fit_NB)
#     for i_l in range(0, nlist):

#         if id_fit_NB[i_l] == 4:
#             id_e0, eboot, id_w0, wboot = get_trigonometric_posterior(
#                 i_NB, i_l, cols[i_l], idpar_NB, posterior
#             )
#             name_der.append(id_e0)
#             name_der.append(id_w0)
#             der_posterior.append(eboot)
#             der_posterior.append(wboot)

#             id_i0, iboot, id_lN0, lNboot = get_trigonometric_posterior(
#                 i_NB, i_l + 3, cols[i_l + 3], idpar_NB, posterior
#             )
#             name_der.append(id_i0)
#             name_der.append(id_lN0)
#             der_posterior.append(iboot)
#             der_posterior.append(lNboot)

#             if idpar_NB[i_l + 2][0:2] != "mA":
#                 mAboot = (
#                     posterior[:, cols[i_l + 2]] - wboot - lNboot
#                 ) % 360.0  # mA = lambda - w - lN
#                 _, mAboot = get_proper_posterior_correlated(mAboot)
#                 name_der.append("{:s}{}".format("mA", i_NB + 1))
#                 der_posterior.append(mAboot)

#     return name_der, der_posterior


# ==============================================================================


# def derived_parameters_case4(
#     idpar_NB, id_fit_NB, i_NB, cols, kep_elem, parameters
# ):  # fit ecosw//sqrt(e)cosw, esinw//sqrt(e)sinw, lambda, cosicoslN, cosisinlN

#     name_der = []
#     der_par = []

#     nlist = len(id_fit_NB)
#     for i_l in range(0, nlist):

#         if id_fit_NB[i_l] == 4:
#             id_e0, e_par, id_w0, w_par = get_trigonometric_parameters(
#                 i_NB, i_l, cols[i_l], idpar_NB, parameters
#             )
#             name_der.append(id_e0)
#             name_der.append(id_w0)
#             der_par.append(e_par)
#             der_par.append(w_par)

#             id_i0, i_par, id_lN0, lN_par = get_trigonometric_parameters(
#                 i_NB, i_l + 3, cols[i_l + 3], idpar_NB, parameters
#             )
#             name_der.append(id_i0)
#             name_der.append(id_lN0)
#             der_par.append(i_par)
#             der_par.append(lN_par)

#             if idpar_NB[i_l + 2][0:2] != "mA":
#                 mA_par = (
#                     parameters[cols[i_l + 2]] - w_par - lN_par
#                 ) % 360.0  # mA = lambda - w - lN
#                 name_der.append("{:s}{}".format("mA", i_NB + 1))
#                 der_par.append(mA_par)

#     return name_der, der_par


# ==============================================================================


def derived_posterior_check(derived_names_in, derived_posterior_in):
    # check distribution of angles...
    derived_posterior = derived_posterior_in.copy()
    derived_names = decode_list(derived_names_in)
    n_der = np.shape(derived_names)[0]
    for ider in range(0, n_der):
        if (
            "w" in derived_names[ider]
            or "mA" in derived_names[ider]
            or "lN" in derived_names[ider]
        ):
            cosp = np.cos(derived_posterior_in[:, ider] * deg2rad)
            sinp = np.sin(derived_posterior_in[:, ider] * deg2rad)
            p_scale = np.arctan2(sinp, cosp) * rad2deg
            p_mod = (p_scale + 360.0) % 360.0
            _, temp_post = get_good_distribution(p_scale, p_mod)
            derived_posterior[:, ider] = temp_post

    return derived_posterior


# ==============================================================================


def derived_parameters_check(derived_names, derived_parameters_in, derived_posterior):
    derived_parameters = derived_parameters_in.copy()
    n_der = np.shape(derived_names)[0]
    for ider in range(0, n_der):
        # print derived_names[ider]
        if (
            "w" in derived_names[ider]
            or "mA" in derived_names[ider]
            or "lN" in derived_names[ider]
        ):
            if np.min(derived_posterior[:, ider]) < 0.0:
                if np.max(derived_posterior[:, ider]) < 0.0:
                    derived_parameters[ider] = derived_parameters[ider] % -360.0
            else:
                derived_parameters[ider] = derived_parameters[ider] % 360.0

    return derived_parameters


# ==============================================================================
# IT HAS TO BE MODIFIED!!!


# def compute_derived_posterior(
#     idpar, kep_elem_in, id_fit, case_list, cols_list, posterior, conv_factor=1.0
# ):
#     NB = len(case_list)

#     if NB == 2:
#         kep_elem = np.array(kep_elem_in).reshape((1, -1))
#     else:
#         kep_elem = kep_elem_in

#     id_derived = []
#     der_posterior = []

#     cc_fit = 0

#     for i_NB in range(0, NB):
#         nfit_NB = len(id_fit[i_NB])

#         if nfit_NB > 0:

#             idpar_NB = idpar[
#                 cols_list[i_NB][0] : cols_list[i_NB][-1] + 1
#             ]  # names of the parameter for the body i_NB
#             id_fit_NB = id_fit[
#                 i_NB
#             ]  # integers that identify the proper type of the fitted parameters

#             for i_fit in range(0, nfit_NB):

#                 if (
#                     id_fit[i_NB][i_fit] == 1
#                 ):  # convert Mp and Mp/Ms into mass unit specified by the user <-> mtype

#                     if "Ms" in idpar[cc_fit]:
#                         mboot = posterior[:, cc_fit] * conv_factor
#                         xid = "m{:d}".format(i_NB + 1)
#                         id_derived.append(xid)
#                         der_posterior.append(mboot)

#                 cc_fit += 1

#             id_temp, der_temp = [], []

#             if case_list[i_NB][0] == 0:
#                 id_temp, der_temp = derived_posterior_case0(
#                     idpar_NB, id_fit_NB, i_NB, cols_list[i_NB], posterior
#                 )

#             elif case_list[i_NB][0] == 1:
#                 id_temp, der_temp = derived_posterior_case1(
#                     idpar_NB,
#                     id_fit_NB,
#                     i_NB,
#                     cols_list[i_NB],
#                     kep_elem[i_NB - 1, :],
#                     posterior,
#                 )

#             elif case_list[i_NB][0] == 2:
#                 id_temp, der_temp = derived_posterior_case2(
#                     idpar_NB,
#                     id_fit_NB,
#                     i_NB,
#                     cols_list[i_NB],
#                     kep_elem[i_NB - 1, :],
#                     posterior,
#                 )

#             elif case_list[i_NB][0] == 3:
#                 id_temp, der_temp = derived_posterior_case3(
#                     idpar_NB,
#                     id_fit_NB,
#                     i_NB,
#                     cols_list[i_NB],
#                     kep_elem[i_NB - 1, :],
#                     posterior,
#                 )

#             elif case_list[i_NB][0] == 4:
#                 id_temp, der_temp = derived_posterior_case4(
#                     idpar_NB,
#                     id_fit_NB,
#                     i_NB,
#                     cols_list[i_NB],
#                     kep_elem[i_NB - 1, :],
#                     posterior,
#                 )

#             id_derived.append(id_temp)
#             der_posterior.append(der_temp)

#     n_der = 0
#     n_xder = len(der_posterior)
#     if n_xder > 0:
#         n_der_single = []
#         derived_post = []
#         names_derived = []

#         for ii in range(0, n_xder):
#             temp_names = id_derived[ii]
#             ntemp = np.size(temp_names)
#             nls = len(np.shape(temp_names))
#             n_der_single.append(ntemp)
#             n_der += ntemp
#             temp_der = der_posterior[ii]

#             if ntemp > 0:

#                 if ntemp == 1 and nls == 0:
#                     derived_post.append(temp_der)
#                     names_derived.append(temp_names)
#                 else:
#                     for jj in range(0, ntemp):
#                         derived_post.append(temp_der[jj])
#                         names_derived.append(temp_names[jj])

#     # convert l2j_x to jitter_x
#     for name_par in idpar:
#         if "l2j_" in name_par:
#             names_derived.append(name_par.replace("l2j_", "jitter_"))
#             derived_post.append(2.0 ** posterior[:, idpar.index(name_par)])

#     return (
#         np.array(names_derived, dtype=str),
#         np.array(derived_post, dtype=np.float64).T,
#     )


# ==============================================================================

# def compute_derived_parameters(
#     idpar, kep_elem_in, id_fit, case_list, cols_list, parameters, conv_factor=1.0
# ):
#     NB = len(case_list)

#     if NB == 2:
#         kep_elem = np.array(kep_elem_in).reshape((1, -1))
#     else:
#         kep_elem = kep_elem_in

#     id_derived = []
#     der_par = []

#     cc_fit = 0

#     for i_NB in range(0, NB):
#         nfit_NB = len(id_fit[i_NB])

#         if nfit_NB > 0:

#             idpar_NB = idpar[
#                 cols_list[i_NB][0] : cols_list[i_NB][-1] + 1
#             ]  # names of the parameter for the body i_NB
#             id_fit_NB = id_fit[
#                 i_NB
#             ]  # integers that identify the proper type of the fitted parameters

#             for i_fit in range(0, nfit_NB):

#                 if (
#                     id_fit[i_NB][i_fit] == 1
#                 ):  # convert Mp and Mp/Ms into mass unit specified by the user <-> mtype

#                     if "Ms" in idpar[cc_fit]:
#                         m_par = parameters[cc_fit] * conv_factor
#                         xid = "m{:d}".format(i_NB + 1)
#                         id_derived.append(xid)
#                         der_par.append(m_par)

#                 cc_fit += 1

#             id_temp, der_temp = [], []

#             if case_list[i_NB][0] == 0:
#                 id_temp, der_temp = derived_parameters_case0(
#                     idpar_NB, id_fit_NB, i_NB, cols_list[i_NB], parameters
#                 )

#             elif case_list[i_NB][0] == 1:
#                 id_temp, der_temp = derived_parameters_case1(
#                     idpar_NB,
#                     id_fit_NB,
#                     i_NB,
#                     cols_list[i_NB],
#                     kep_elem[i_NB - 1, :],
#                     parameters,
#                 )

#             elif case_list[i_NB][0] == 2:
#                 id_temp, der_temp = derived_parameters_case2(
#                     idpar_NB,
#                     id_fit_NB,
#                     i_NB,
#                     cols_list[i_NB],
#                     kep_elem[i_NB - 1, :],
#                     parameters,
#                 )

#             elif case_list[i_NB][0] == 3:
#                 id_temp, der_temp = derived_parameters_case3(
#                     idpar_NB,
#                     id_fit_NB,
#                     i_NB,
#                     cols_list[i_NB],
#                     kep_elem[i_NB - 1, :],
#                     parameters,
#                 )

#             elif case_list[i_NB][0] == 4:
#                 id_temp, der_temp = derived_parameters_case4(
#                     idpar_NB,
#                     id_fit_NB,
#                     i_NB,
#                     cols_list[i_NB],
#                     kep_elem[i_NB - 1, :],
#                     parameters,
#                 )

#             id_derived.append(id_temp)
#             der_par.append(der_temp)

#     n_der = 0
#     n_xder = len(der_par)
#     if n_xder > 0:
#         n_der_single = []
#         derived_par = []
#         names_derived = []

#         for ii in range(0, n_xder):
#             temp_names = id_derived[ii]
#             ntemp = np.size(temp_names)
#             nls = len(np.shape(temp_names))
#             n_der_single.append(ntemp)
#             n_der += ntemp
#             temp_der = der_par[ii]

#             if ntemp > 0:

#                 if ntemp == 1 and nls == 0:
#                     derived_par.append(temp_der)
#                     names_derived.append(temp_names)
#                 else:
#                     for jj in range(0, ntemp):
#                         derived_par.append(temp_der[jj])
#                         names_derived.append(temp_names[jj])

#     # convert l2j_x to jitter_x
#     for name_par in idpar:
#         if "l2j_" in name_par:
#             names_derived.append(name_par.replace("l2j_", "jitter_"))
#             derived_par.append(2.0 ** parameters[idpar.index(name_par)])

#     return np.array(names_derived, dtype=str), np.array(derived_par, dtype=np.float64).T


def update_parameterisation_parameters(fitting_names, fitting_parameters):

    new_fitting_names = fitting_names.copy()
    new_fitting_parameters = fitting_parameters.copy()

    for i_p, name in enumerate(fitting_names):
        if "icoslN" in name:
            id_body = name.split("lN")[1]
            new_fitting_names[i_p] = "i{}".format(id_body)
            new_fitting_names[i_p + 1] = "lN{}".format(id_body)

            # take current parameters and posterior
            cosln = fitting_parameters[i_p]
            sinln = fitting_parameters[i_p + 1]
            lon = (np.arctan2(sinln, cosln) * cst.rad2deg) % 360.0

            # if the parameter was icoslN, isinlN
            inc = np.sqrt(cosln * cosln + sinln * sinln)

            # or if it was cosicoslN, cosisinlN
            if "cosi" == name[0:4]:
                inc = np.arccos(inc) * cst.rad2deg

            new_fitting_parameters[i_p] = inc
            new_fitting_parameters[i_p + 1] = lon

    return new_fitting_names, new_fitting_parameters


def update_parameterisation_posterior(fitting_names, fitting_posterior):

    new_fitting_names = fitting_names.copy()
    new_fitting_posterior = fitting_posterior.copy()

    for i_p, name in enumerate(fitting_names):
        if "icoslN" in name:
            id_body = name.split("lN")[1]
            new_fitting_names[i_p] = "i{}".format(id_body)
            new_fitting_names[i_p + 1] = "lN{}".format(id_body)

            # take current parameters and posterior
            cosln = fitting_posterior[:, i_p]
            sinln = fitting_posterior[:, i_p + 1]
            lon = (np.arctan2(sinln, cosln) * cst.rad2deg) % 360.0

            # if the parameter was icoslN, isinlN
            inc = np.sqrt(cosln * cosln + sinln * sinln)

            # or if it was cosicoslN, cosisinlN
            if "cosi" == name[0:4]:
                inc = np.arccos(inc) * cst.rad2deg

            new_fitting_posterior[:, i_p] = inc
            new_fitting_posterior[:, i_p + 1] = lon

    return new_fitting_names, new_fitting_posterior


def update_parameterisation_chains(fitting_names, fitting_chains):

    new_fitting_names = fitting_names.copy()
    new_fitting_chains = fitting_chains.copy()

    for i_p, name in enumerate(fitting_names):
        if "icoslN" in name:
            id_body = name.split("lN")[1]
            new_fitting_names[i_p] = "i{}".format(id_body)
            new_fitting_names[i_p + 1] = "lN{}".format(id_body)

            # take current parameters and posterior
            cosln = fitting_chains[:, :, i_p]
            sinln = fitting_chains[:, :, i_p + 1]
            lon = (np.arctan2(sinln, cosln) * cst.rad2deg) % 360.0

            # if the parameter was icoslN, isinlN
            inc = np.sqrt(cosln * cosln + sinln * sinln)

            # or if it was cosicoslN, cosisinlN
            if "cosi" == name[0:4]:
                inc = np.arccos(inc) * cst.rad2deg

            new_fitting_chains[:, :, i_p] = inc
            new_fitting_chains[:, :, i_p + 1] = lon

    return new_fitting_names, new_fitting_chains

def scale_angle(val_deg):

    val_rad = val_deg * cst.deg2rad
    cv = np.cos(val_rad)
    sv = np.sin(val_rad)
    scale_rad = np.arctan2(sv, cv)
    scale_deg = scale_rad * cst.rad2deg

    return scale_deg

def search_scale_parameter(par, par_type):

    scale_par = par.copy()
    for i_p, p in enumerate(par):
        if "mod" in par_type[i_p]:
            scale_par[i_p] = p%360.0
        elif "scale" in par_type[i_p]:
            scale_par[i_p] = scale_angle(p)

    return scale_par

def compute_physical_parameters(
    n_bodies,
    par_fit,
    names_fit,
    all_system_parameters,
    mass_conv_factor=1.0,
    radius_conv_factor=1.0,
    chains_full_fit=None,
    chains_posterior_fit=None,
    posterior_fit=None,
    mass_post_conv_factor=[1.0, 1.0, 1.0],
    radius_post_conv_factor=[1.0, 1.0, 1.0],
):

    names_phys = []
    par_phys = []
    chains_full_phys = []
    chains_posterior_phys = []
    posterior_phys = []
    phys_type = []

    nfit = len(par_fit)

    for i_body in range(1, n_bodies):

        idx_body = i_body + 1

        mass = None
        radius = None
        ecc = None
        argp = None
        meana = None
        meanl = None
        inc = None
        longn = None

        mass_post = None
        radius_post = None
        ecc_post = None
        argp_post = None
        meana_post = None
        meanl_post = None
        inc_post = None
        longn_post = None

        idx_p = 0
        # mass: mXMs to Mp
        keyp = "m{:d}Ms".format(idx_body)
        if keyp in names_fit:
            idx_p = names_fit.index(keyp)
            mass = par_fit[idx_p] * mass_conv_factor
            names_phys.append("m{:d}".format(idx_body))
            par_phys.append(mass)
            phys_type.append("-")
            if posterior_fit is not None:
                mass_post = posterior_fit[:, idx_p] * mass_post_conv_factor[2]
                posterior_phys.append(mass_post)
            if chains_full_fit is not None:
                mass_f_post = chains_full_fit[:,:, idx_p] * mass_post_conv_factor[0]
                chains_full_phys.append(mass_f_post)
            if chains_posterior_fit is not None:
                mass_p_post = chains_posterior_fit[:,:, idx_p] * mass_post_conv_factor[1]
                chains_posterior_phys.append(mass_p_post)

        idx_p = 0
        # radius: rXRs to Rp (but you shouldn't fit Radius of planet)
        keyp = "r{:d}Rs".format(idx_body)
        if keyp in names_fit:
            idx_p = names_fit.index(keyp)
            radius = par_fit[idx_p] * radius_conv_factor
            names_phys.append("r{:d}".format(idx_body))
            par_phys.append(radius)
            phys_type.append("-")
            if posterior_fit is not None:
                radius_post = posterior_fit[:, idx_p] * radius_post_conv_factor[2]
                posterior_phys.append(radius_post)
            if chains_full_fit is not None:
                radius_f_post = chains_full_fit[:,:, idx_p] * radius_post_conv_factor[0]
                chains_full_phys.append(radius_f_post)
            if chains_posterior_fit is not None:
                radius_p_post = chains_posterior_fit[:,:, idx_p] * radius_post_conv_factor[1]
                chains_posterior_phys.append(radius_p_post)

        idx_p, idx_q = 0, 0
        # ecc, argp: secoswX, sesinwX (also check sqrte) to eX, wX
        keyp = "secosw{:d}".format(idx_body)
        if keyp in names_fit:
            idx_p = names_fit.index(keyp)
        elif keyp.replace("se", "sqrte") in names_fit:
            idx_p = names_fit.index(keyp)
        if idx_p > 0:
            idx_q = idx_p + 1
            ecc = par_fit[idx_p] * par_fit[idx_p] + par_fit[idx_q] * par_fit[idx_q]
            if ecc <= 1.0e-9:
                ecc = 0.0
                argp_scale = 90.0
                argp = argp_scale
            else:
                argp_scale = np.arctan2(par_fit[idx_q], par_fit[idx_p]) * cst.rad2deg
                argp = argp_scale % 360.0

            names_phys.append("e{:d}".format(idx_body))
            par_phys.append(ecc)
            phys_type.append("-")
            names_phys.append("w{:d}".format(idx_body))
            par_phys.append(argp)
            phys_type.append("mod")

            if posterior_fit is not None:
                secwpost = posterior_fit[:, idx_p]
                seswpost = posterior_fit[:, idx_q]
                ecc_post = secwpost * secwpost + seswpost * seswpost
                argp_post = (np.arctan2(seswpost, secwpost) * cst.rad2deg) % 360.0
                _, argp_post, argp_type = get_proper_posterior_correlated(
                    argp_post, type_out=True
                )
                sel_ecc = ecc_post <= 1.0e-9
                ecc_post[sel_ecc] = 0.0
                posterior_phys.append(ecc_post)
                argp_post[sel_ecc] = 90.0
                posterior_phys.append(argp_post)
                # substitute proper value
                if "mod" in argp_type:
                    argp %= 360.0
                else:
                    argp = argp_scale
                par_phys[-1] = argp
                phys_type[-1] = argp_type
                
            if chains_full_fit is not None:
                cc = chains_full_fit[:,:, idx_p]
                ss = chains_full_fit[:,:, idx_q]
                f_post = cc*cc + ss*ss
                sel_ecc = f_post <= 1.0e-9
                f_post[sel_ecc] = 0.0
                chains_full_phys.append(f_post)
                f_post = np.arctan2(ss, cc) * cst.rad2deg
                if "mod" in phys_type[-1]:
                    f_post %= 360.0
                f_post[sel_ecc] = 90.0
                chains_full_phys.append(f_post)
                argp_f_post = f_post

            if chains_posterior_fit is not None:
                cc = chains_posterior_fit[:,:, idx_p]
                ss = chains_posterior_fit[:,:, idx_q]
                p_post = cc*cc + ss*ss
                sel_ecc = p_post <= 1.0e-9
                p_post[sel_ecc] = 0.0
                chains_posterior_phys.append(p_post)
                p_post = np.arctan2(ss, cc) * cst.rad2deg
                if "mod" in phys_type[-1]:
                    p_post %= 360.0
                p_post[sel_ecc] = 90.0
                chains_posterior_phys.append(p_post)
                argp_p_post = p_post

        # ecosw/esinw for backward compatibility
        idx_p, idx_q = 0, 0
        keyp = "ecosw{:d}".format(idx_body)
        if keyp in names_fit:
            idx_p = names_fit.index(keyp)
            idx_q = idx_p + 1
            ecc = np.sqrt(
                par_fit[idx_p] * par_fit[idx_p] + par_fit[idx_q] * par_fit[idx_q]
            )
            if ecc <= 1.0e9:
                ecc = 0.0
                argp_scale = 90.0
                argp = argp_scale
            else:
                argp_scale = np.arctan2(par_fit[idx_q], par_fit[idx_p]) * cst.rad2deg
                argp = argp_scale % 360.0

            names_phys.append("e{:d}".format(idx_body))
            par_phys.append(ecc)
            phys_type.append("-")
            names_phys.append("w{:d}".format(idx_body))
            par_phys.append(argp)
            phys_type.append("mod")

            if posterior_fit is not None:
                cc = posterior_fit[:, idx_p]
                ss = posterior_fit[:, idx_q]
                ecc_post = np.sqrt(cc * cc + ss * ss)
                argp_post = (np.arctan2(ss, cc) * cst.rad2deg) % 360.0
                _, argp_post, argp_type = get_proper_posterior_correlated(
                    argp_post, type_out=True
                )
                sel_ecc = ecc_post <= 1.0e9
                ecc_post[sel_ecc] = 0.0
                argp_post[sel_ecc] = 90.0
                posterior_phys.append(ecc_post)
                posterior_phys.append(argp_post)
                # substitute proper value
                if "mod" in argp_type:
                    argp %= 360.0
                else:
                    argp = argp_scale
                par_phys[-1] = argp
                phys_type[-1] = argp_type

            if chains_full_fit is not None:
                cc = chains_full_fit[:,:, idx_p]
                ss = chains_full_fit[:,:, idx_q]
                f_post = np.sqrt(cc*cc + ss*ss)
                sel_ecc = f_post <= 1.0e-9
                f_post[sel_ecc] = 0.0
                chains_full_phys.append(f_post)
                f_post = np.arctan2(ss, cc) * cst.rad2deg
                if "mod" in phys_type[-1]:
                    f_post %= 360.0
                f_post[sel_ecc] = 90.0
                chains_full_phys.append(f_post)
                argp_f_post = f_post

            if chains_posterior_fit is not None:
                cc = chains_posterior_fit[:,:, idx_p]
                cc = chains_posterior_fit[:,:, idx_q]
                p_post = np.sqrt(cc*cc + ss*ss)
                sel_ecc = p_post <= 1.0e-9
                p_post[sel_ecc] = 0.0
                chains_posterior_phys.append(p_post)
                p_post = np.arctan2(ss, cc) * cst.rad2deg
                if "mod" in phys_type[-1]:
                    p_post %= 360.0
                p_post[sel_ecc] = 90.0
                chains_posterior_phys.append(p_post)
                argp_p_post = p_post

        # icoslN/isinlN for backward compatibility
        idx_p, idx_q = 0, 0
        keyp = "icoslN{:d}".format(idx_body)
        if keyp in names_fit:
            idx_p = names_fit.index(keyp)
            idx_q = idx_p + 1
            inc = np.sqrt(
                par_fit[idx_p] * par_fit[idx_p] + par_fit[idx_q] * par_fit[idx_q]
            )
            longn_scale = np.arctan2(par_fit[idx_q], par_fit[idx_p]) * cst.rad2deg
            longn = longn_scale % 360.0

            names_phys.append("i{:d}".format(idx_body))
            par_phys.append(inc)
            phys_type.append("-")
            names_phys.append("lN{:d}".format(idx_body))
            par_phys.append(longn)
            phys_type.append("mod")

            if posterior_fit is not None:
                iclNpost = posterior_fit[:, idx_p]
                islNpost = posterior_fit[:, idx_q]
                inc_post = np.sqrt(iclNpost * iclNpost + islNpost * islNpost)
                longn_post = (np.arctan2(islNpost, iclNpost) * cst.rad2deg) % 360.0
                _, longn_post, longn_type = get_proper_posterior_correlated(
                    longn_post, type_out=True
                )
                posterior_phys.append(inc_post)
                posterior_phys.append(longn_post)
                # substitute proper value
                if "mod" in argp_type:
                    longn %= 360.0
                else:
                    longn =longn_scale
                par_phys[-1] = longn
                phys_type[-1] = longn_type
            
            if chains_full_fit is not None:
                cc = chains_full_fit[:,:, idx_p]
                ss = chains_full_fit[:,:, idx_q]
                f_post = np.sqrt(cc*cc + ss*ss)
                chains_full_phys.append(f_post)
                f_post = np.arctan2(ss, cc) * cst.rad2deg
                if "mod" in phys_type[-1]:
                    f_post %= 360.0
                chains_full_phys.append(f_post)
                longn_f_post = f_post

            if chains_posterior_fit is not None:
                cc = chains_posterior_fit[:,:, idx_p]
                ss = chains_posterior_fit[:,:, idx_q]
                p_post = np.sqrt(cc*cc + ss*ss)
                chains_posterior_phys.append(p_post)
                sel_ecc = p_post <= 1.0e-9
                p_post[sel_ecc] = 0.0
                chains_posterior_phys.append(p_post)
                p_post = np.arctan2(ss, cc) * cst.rad2deg
                if "mod" in phys_type[-1]:
                    p_post %= 360.0
                p_post[sel_ecc] = 90.0
                chains_posterior_phys.append(p_post)
                longn_p_post = p_post

        idx_p = 0
        # meanA: lambdaX to mAX
        keyp = "lambda{:d}".format(idx_body)
        if keyp in names_fit:
            idx_p = names_fit.index(keyp)
            meanl = par_fit[idx_p]
            # print("DEBUG: body {} meanl = {}".format(idx_body, meanl))
            if posterior_fit is not None:
                meanl_post = posterior_fit[:, idx_p]
            if chains_full_fit is not None:
                meanl_f_post = chains_full_fit[:, :, idx_p]
            if chains_posterior_fit is not None:
                meanl_p_post = chains_posterior_fit[:, :, idx_p]

            idx_q = 0
            # get to proper wX
            if argp is None:
                keyq = "w{:d}".format(idx_body)
                if keyq in names_fit:
                    idx_q = names_fit.index(keyq)
                    argp = par_fit[idx_q]
                    if posterior_fit is not None:
                        argp_post = posterior_fit[:, idx_q]
                    if chains_full_fit is not None:
                        argp_f_post = chains_full_fit[:, :, idx_q]
                    if chains_posterior_fit is not None:
                        argp_p_post = chains_posterior_fit[:, :, idx_q]
                else:
                    idx_q = (7 + (idx_body - 2) * 8) - 1  # fortran index to python: -1
                    argp = all_system_parameters[idx_q]
                    if posterior_fit is not None:
                        argp_post = argp
                    if chains_full_fit is not None:
                        argp_f_post = argp
                    if chains_posterior_fit is not None:
                        argp_p_post = argp
            # else:
            #     if posterior_fit is not None:
            #         argp_post = posterior_fit[:, idx_q]
            #     if chains_full_fit is not None:
            #         argp_f_post = chains_full_fit[:, :, idx_q]
            #     if chains_posterior_fit is not None:
            #         argp_p_post = chains_posterior_fit[:, :, idx_q]

            # print("DEBUG: body {} argp = {}".format(idx_body, argp))

            idx_q = 0
            # get proper lNX
            if longn is None:
                keyq = "lN{:d}".format(idx_body)
                if keyq in names_fit:
                    idx_q = names_fit.index(keyq)
                    longn = par_fit[idx_q]
                    if posterior_fit is not None:
                        longn_post = posterior_fit[:, idx_q]
                    if chains_full_fit is not None:
                        lf_post = chains_full_fit[:, :, idx_q]
                    if chains_posterior_fit is not None:
                        lp_post = chains_posterior_fit[:, :, idx_q]
                else:
                    idx_q = (10 + (idx_body - 2) * 8) - 1  # fortran index to python: -1
                    longn = all_system_parameters[idx_q]
                    if posterior_fit is not None:
                        longn_post = longn
                    if chains_full_fit is not None:
                        longn_f_post = longn
                    if chains_posterior_fit is not None:
                        longn_p_post = longn
            # else:
            #     if posterior_fit is not None:
            #         longn_post = posterior_fit[:, idx_q]
            #     if chains_full_fit is not None:
            #         lf_post = chains_full_fit[:, :, idx_q]
            #     if chains_posterior_fit is not None:
            #         lp_post = chains_posterior_fit[:, :, idx_q]

            # print("DEBUG: body {} longn = {}".format(idx_body, longn))

            meana = (meanl - longn - argp) % 360.0
            # meana_scale = np.arctan2(np.sin(meana*cst.deg2rad), np.cos(meana*cst.deg2rad)) * cst.rad2deg
            meana_scale = get_arctan_angle(meana)
            names_phys.append("mA{:d}".format(idx_body))
            par_phys.append(meana)
            phys_type.append("mod")

            # print("DEBUG: body {} meana = {}".format(idx_body, meana))

            if posterior_fit is not None:
                meana_post = (meanl_post - longn_post - argp_post) % 360.0
                _, meana_post, meana_type = get_proper_posterior_correlated(
                    meana_post, type_out=True
                )
                posterior_phys.append(meana_post)
                # substitute proper value
                if "mod" in meana_type:
                    meana %= 360.0
                else:
                    meana = meana_scale
                par_phys[-1] = meana
                phys_type[-1] = meana_type
            if chains_full_fit is not None:
                meana_f_post = (meanl_f_post - longn_f_post - argp_f_post) % 360.0
                # _, meana_post, _ = get_proper_posterior_correlated(
                #     meana_post, type_out=True
                # )
                if "scale" in meana_type:
                    meana_f_post = get_arctan_angle(meana_f_post)
                chains_full_phys.append(meana_f_post)
            if chains_posterior_fit is not None:
                meana_p_post = (meanl_p_post - longn_p_post - argp_p_post) % 360.0
                # _, meana_post, _ = get_proper_posterior_correlated(
                #     meana_post, type_out=True
                # )
                if "scale" in meana_type:
                    meana_p_post = get_arctan_angle(meana_p_post)
                chains_posterior_phys.append(meana_p_post)

    # check jitter
    # convert l2j_x to jitter_x
    idx_p = 0
    # l2j_sel = np.array(["l2j_" in n for n in names_fit])
    for idx_p, n_fit in enumerate(names_fit):
        if "l2j_" in n_fit:
            # idx_p = names_fit.index(n_fit)
            jitter = 2.0 ** par_fit[idx_p]
            names_phys.append(n_fit.replace("l2j_", "jitter_"))
            par_phys.append(jitter)
            phys_type.append("-")
            if posterior_fit is not None:
                jitter_post = 2.0 ** posterior_fit[:, idx_p]
                posterior_phys.append(jitter_post)
            if chains_full_fit is not None:
                f_post = 2.0 ** chains_full_fit[:, :, idx_p]
                chains_full_phys.append(f_post)
            if chains_posterior_fit is not None:
                p_post = 2.0 ** chains_posterior_fit[:, :, idx_p]
                chains_posterior_phys.append(p_post)

    names_phys = np.array(names_phys)
    par_phys = np.array(par_phys)
    phys_type = np.array(phys_type)


    if posterior_fit is not None:
        posterior_phys = np.column_stack(posterior_phys)
    
    if chains_full_fit is not None:
        chains_full_phys = np.array(chains_full_phys).swapaxes(0,1).swapaxes(1,2)

    if chains_posterior_fit is not None:
        chains_posterior_phys = np.array(chains_posterior_phys).swapaxes(0,1).swapaxes(1,2)

    return names_phys, par_phys, posterior_phys, phys_type, chains_full_phys, chains_posterior_phys


# ==============================================================================


def compare_par_post_angle(der_par, der_post):
    cospar = np.cos(der_par * deg2rad)
    sinpar = np.sin(der_par * deg2rad)

    par_scale = np.arctan2(sinpar, cospar) * rad2deg
    par_mod = (par_scale + 360.0) % 360.0

    cospos = np.cos(der_post * deg2rad)
    sinpos = np.sin(der_post * deg2rad)

    pos_scale = np.arctan2(sinpos, cospos) * rad2deg
    pos_mod = (pos_scale + 360.0) % 360.0

    delta_scale = np.abs(np.mean(pos_scale) - par_scale)
    delta_mod = np.abs(np.mean(pos_mod) - par_mod)

    if np.abs(delta_scale - delta_mod) <= eps64bit:
        return par_mod, pos_mod
    else:
        if delta_scale > delta_mod:
            return par_mod, pos_mod
        else:
            return par_scale, pos_scale
    return par_mod, pos_mod


# ==============================================================================


def adjust_derived_parameters(derived_names, derived_par, derived_post):
    n_der = np.shape(derived_par)[0]
    derived_par_adj = derived_par.copy()
    derived_post_adj = derived_post.copy()

    for i_d in range(0, n_der):
        if derived_names[i_d][0] in ["w", "i"] or derived_names[i_d][0:2] in [
            "lN",
            "mA",
        ]:
            derived_par_adj[i_d], derived_post_adj[:, i_d] = compare_par_post_angle(
                derived_par[i_d], derived_post[:, i_d]
            )

    return derived_par_adj, derived_post_adj


# ==============================================================================


def get_header(perc_val):

    top_header = "# {:15s} {:20s} {:<23s} {:<23s} {:<23s} {:<23s} {:<23s} {:<23s} {:<23s} {:<23s} {:<23s}".format(
        "-",
        "-",
        "-",
        "+/-1sigma",
        "MAD",
        "-1sigma",
        "+1sigma",
        "-2sigma",
        "+2sigma",
        "-3sigma",
        "+3sigma",
    )
    header = "# {:15s} {:20s} {:23s}".format("name", "unit", "parameter")
    perc_str = " ".join(
        [
            "{:23s}".format("{:4.2f}-th".format(perc_val[i_p]))
            for i_p in range(len(perc_val))
        ]
    )
    header = "{:s} {:s}".format(header, perc_str)

    return top_header, header


# ==============================================================================


def print_parameters(
    top_header,
    header,
    name_parameters,
    unit_parameters,
    parameters,
    sigma_parameters=None,
    output=None,
):

    print_both(top_header, output)
    print_both(header, output)

    if sigma_parameters is not None:
        n_sig, n_par = np.shape(sigma_parameters)
        for i_p in range(0, n_par):
            unit = unit_parameters[i_p].strip()
            if len(unit) == 0:
                unit = "-"
            sigma_line = " ".join(
                ["{:+23.15e}".format(sigma_parameters[ii, i_p]) for ii in range(n_sig)]
            )
            line = "  {:15s} {:20s} {:+23.15e} {:s}".format(
                name_parameters[i_p],
                unit,
                parameters[i_p],
                sigma_line,
            )
            print_both(line, output)
    else:
        print_both("PARAMETERS NOT AVAILABLE", output)

    return


# ==============================================================================


def print_confidence_intervals(
    percentiles,
    conf_interv=None,
    name_parameters=None,
    unit_parameters=None,
    output=None,
):

    if conf_interv is not None:
        header = "# {:15s} {:20s}".format("name", "unit")
        perc_str = " ".join(
            [
                "{:23s}".format("{:4.2f}-th".format(percentiles[i_p]))
                for i_p in range(len(percentiles))
            ]
        )
        header = "{:s} {:s}".format(header, perc_str)
        print_both(header, output)

        n_par = len(name_parameters)
        for i_p in range(0, n_par):
            unit = unit_parameters[i_p].strip()
            if len(unit) == 0:
                unit = "-"
            ci_line = " ".join(
                [
                    "{:23.16e}".format(conf_interv[ii, i_p])
                    for ii in range(0, len(percentiles))
                ]
            )
            line = " {:15s} {:23s} {:s}".format(name_parameters[i_p], unit, ci_line)
            print_both(line, output)

    else:
        print_both("EMPTY", output)

    return


# ==============================================================================


def print_hdi(
    conf_interv=None, name_parameters=None, unit_parameters=None, output=None
):

    if conf_interv is not None:
        header = "# {:15s} {:20s} {:23s} {:23s} {:23s} {:23s} {:23s} {:23s}".format(
            "name",
            "unit",
            "HDI(68.27%)lower",
            "HDI(68.27%)upper",
            "HDI(95.44%)lower",
            "HDI(95.44%)upper",
            "HDI(99.74%)lower",
            "HDI(99.74%)upper",
        )
        print_both(header, output)

        # n_par, n_hdi = np.shape(conf_interv)
        n_hdi, n_par = np.shape(conf_interv)
        for i_p in range(0, n_par):
            unit = unit_parameters[i_p].strip()
            if len(unit) == 0:
                unit = "-"
            ci_line = " ".join(
                ["{:23.16e}".format(conf_interv[ii, i_p]) for ii in range(n_hdi)]
            )
            line = "  {:15s} {:23s} {:s}".format(name_parameters[i_p], unit, ci_line)
            print_both(line, output)

    else:
        print_both("EMPTY", output)

    return


# =============================================================================


def read_fitted_file(fitted_file):
    of = open(fitted_file, "r")
    lines = of.readlines()
    of.close()

    names, fitted_par = [], []
    for ll in lines:
        line = ll.strip()
        if line[0] != "#":
            lsplit = line.split()
            names.append(lsplit[0].strip())
            fitted_par.append(np.float64(lsplit[1].strip()))

    return names, np.array(fitted_par, dtype=np.float64)


# ==============================================================================
# STURGES RULE TO COMPUTE THE NUMBER OF BINS FOR HISTOGRAM
# ==============================================================================


def sturges_nbins(npost):
    """
    Following Sturges' rule:
    nbins = log2(n) + 1
    """
    nbins = np.ceil(np.log2(np.float64(npost))).astype(int) + 1

    return nbins


# ==============================================================================
# FREEDMAN-DIACONIS RULE TO COMPUTE THE NUMBER OF BINS FOR HISTOGRAM
# ==============================================================================


def freedman_diaconis_nbins(x):
    """
    Following Freedman-Diaconis rule:
    width = 2 * IRQ / n^1/3
    nbins = [ (max-min) / width ]
    """
    q75 = np.percentile(
        x,
        75.0,
        #
    )
    q25 = np.percentile(
        x,
        25.0,
        #
    )
    irq = np.abs(q75 - q25)
    nx = np.float64(np.shape(x)[0])
    width = 2.0 * irq / np.power(nx, 1.0 / 3.0)
    nbins = np.ceil(np.abs(np.max(x) - np.min(x)) / width).astype(int)

    return nbins


# ==============================================================================
# DOANE'S FORMULA TO COMPUTE THE NUMBER OF BINS FOR HISTOGRAM
# ==============================================================================


def doane_nbins(x):
    """
    Doane's formula:
      nbins = 1 + lon2(n) + log2(1 + |g1|/s_g1)
      n: number of points/sample
      g1 = 3rd-moment-skewness
      s_g1 = sqrt( (6 x (n-2)) / ((n+1) x (n+3)) )
    """
    nxf = np.float64(np.shape(x)[0])
    mux = np.mean(x)
    stdx = np.std(x, ddof=1)
    g1 = np.mean(((x - mux) / stdx) ** 3)  # skew
    s_g1 = np.sqrt(6.0 * (nxf - 2.0) / ((nxf + 1.0) * (nxf + 3.0)))
    # print nxf
    # print mux,stdx
    # print g1, s_g1
    nbins = int(1.0 + np.log2(nxf) + np.log2(1.0 + np.abs(g1) / s_g1))

    return nbins


# ==============================================================================
# GIVEN THE POSTERIOR IT SELECT THE PROPER NBINS FOR ALL THE PARAMETER DISTRIBUTIONS
# ==============================================================================


def get_auto_bins(posterior):

    npost, nfit = np.shape(posterior)
    # call Freedman-Diaconis rule
    nbins_fd = [freedman_diaconis_nbins(posterior[:, ifit]) for ifit in range(nfit)]
    # doane's formula
    nbins_doa = [doane_nbins(posterior[:, ifit]) for ifit in range(nfit)]
    # and Sturges' rule
    nbins = int(np.mean([min(nbins_fd), min(nbins_doa), sturges_nbins(npost)]))

    return nbins


# ==============================================================================
# GIVEN THE POSTERIOR AND THE RULE IT COMPUTES THE NBINS
# ==============================================================================


def get_old_bins(x, rule="sq"):

    nx = np.shape(x)[0]
    nxf = np.float64(nx)

    if "sq" in rule.lower():

        nbins = int(np.sqrt(np.float64(nx)))

    elif "stu" in rule.lower():

        nbins = sturges_nbins(nxf)

    elif "doa" in rule.lower():

        nbins = doane_nbins(x)

    elif "fd" in rule.lower():

        nbins = freedman_diaconis_nbins(x)

    else:  # rice

        nbins = int(np.ceil(2.0 * np.power(nxf, 1.0 / 3.0)))

    return nbins


# ==============================================================================
# HIGH DENSITY INTERVALS - DOING BAYESIAN DATA ANALYSIS by Kruschke
# ==============================================================================


def calculate_hdi(x, nbins, alpha=[0.05], mode_output=False):

    # 68.27% (15.87th-84.13th) ==> alpha = 1. - 0.6827 = 0.3173
    # 95.44% ( 2.28th-97.72th) ==> alpha = 1. - 0.9544 = 0.0456
    # 99.74% ( 0.13th-99.87th) ==> alpha = 1. - 0.9974 = 0.0026

    counts, bin_edges = np.histogram(x, bins=nbins)

    bin_width = bin_edges[1] - bin_edges[0]
    hwidth = 0.5 * bin_width
    bin_mid = [bin_edges[ibin] + hwidth for ibin in range(nbins)]

    # probability mass distribution = (counts / ndata)
    # probability density = probability mass distribution / bin_width
    pmd = np.asarray(counts) / np.float64(np.shape(x)[0])
    # pdd = pmd / bin_width
    spmd = np.sort(pmd)[::-1]
    hdi_ci = []
    for ialpha in alpha:
        lowestHeightIdx = np.min(np.where(np.cumsum(spmd) > (1.0 - ialpha)))
        lowestHeight = spmd[lowestHeightIdx]
        bin_sel = np.asarray(bin_mid)[pmd >= lowestHeight]
        hdi_l = np.min(bin_sel) - hwidth
        hdi_u = np.max(bin_sel) + hwidth
        hdi_ci.append([hdi_l, hdi_u])

    if mode_output:
        max_bin = np.argmax(np.array(counts))
        sel_bin = np.logical_and(x >= bin_edges[max_bin], x < bin_edges[max_bin + 1])

        mode = np.mean(x[sel_bin])
        hdi_ci.append(mode)

    return hdi_ci


# ==============================================================================
# COMPUTES THE HDI/HPD FROM PYASTRONOMY
# ==============================================================================
def hpd(trace, cred=0.6827):
    """
    Estimate the highest probability density interval.

    This function determines the shortest, continuous interval
    containing the specified fraction (cred) of steps of
    the Markov chain. Note that multi-modal distribution
    may require further scrutiny.

    Parameters
    ----------
    trace : array
        The steps of the Markov chain.
    cred : float
        The probability mass to be included in the
        interval (between 0 and 1).

    Returns
    -------
    start, end : float
        The start and end points of the interval.
    """
    cred_def = 0.6827
    if (cred > 1.0) or (cred < 0.0):
        print("CRED HAS TO BE: 0 < cred < 1 ==> setting to cred = {}".format(cred_def))
        cred = cred_def

    # Sort the trace steps in ascending order
    st = np.sort(trace)

    # Number of steps in the chain
    n = len(st)
    # Number of steps to be included in the interval
    nin = int(n * cred)

    # All potential intervals must be 1) continuous and 2) cover
    # the given number of trace steps. Potential start and end
    # points of the HPD are given by
    starts = st[0:-nin]
    ends = st[nin:]
    # All possible widths are
    widths = ends - starts
    # The density is highest in the shortest one
    imin = np.argmin(widths)
    return starts[imin], ends[imin]


# ==============================================================================
def calculate_hpd(trace):

    # cred = 0.6827 <-> -/+ 1 sigma
    # cred = 0.9544 <-> -/+ 2 sigma
    # cred = 0.9974 <-> -/+ 3 sigma

    credval = [0.6827, 0.9544, 0.9974]
    hdi = [hpd(trace, cred=c) for c in credval]

    return hdi


# ==============================================================================
# ==============================================================================

# use in emcee script sqrt(e)cos(w),sqrt(e)sin(w), while in trades use (e)cos(w),(e)sin(w)


def convert_fortran2python_strarray(array_str_in, nfit, str_len=10):

    temp = [
        array_str_in[r].decode("utf-8")[c] for r in range(nfit) for c in range(str_len)
    ]
    temp = np.asarray(temp).reshape((nfit, str_len), order="F")

    array_str_out = ["".join(temp[i, :]).strip() for i in range(0, nfit)]

    return array_str_out


# ==============================================================================


def convert_fortran_charray2python_strararray(array_str_in):

    ndim, str_len = np.shape(array_str_in)
    temp = np.asarray(
        [
            array_str_in[r, c].decode("utf-8")
            for r in range(ndim)
            for c in range(str_len)
        ]
    ).reshape((ndim, str_len))
    array_str_out = ["".join(temp[i, :]).strip() for i in range(0, ndim)]

    return array_str_out


# ==============================================================================


def trades_names_to_emcee(trades_names):
    nfit = np.shape(trades_names)[0]
    emcee_names = list(trades_names)
    for ifit in range(0, nfit):
        if trades_names[ifit][:2] == "ec":
            emcee_names[ifit] = "sqrt{:s}".format(trades_names[ifit])
            emcee_names[ifit + 1] = "sqrt{:s}".format(trades_names[ifit + 1])

    return emcee_names


# ==============================================================================


def emcee_names_to_trades(emcee_names):
    nfit = np.shape(emcee_names)[0]
    trades_names = list(emcee_names)
    for ifit in range(0, nfit):
        if emcee_names[ifit][0:6] == "sqrtec":
            trades_names[ifit] = "{:s}".format(emcee_names[ifit].split("sqrt")[1])
            trades_names[ifit + 1] = "{:s}".format(
                emcee_names[ifit + 1].split("sqrt")[1]
            )

    return trades_names


# ==============================================================================


def e_to_sqrte_boundaries(boundaries_in, names_par):
    nfit = np.shape(boundaries_in)[0]
    boundaries_out = np.array(boundaries_in, dtype=np.float64).copy()
    for ifit in range(0, nfit):
        if names_par[ifit][0:2] == "ec":
            ec = boundaries_in[ifit, :]
            es = boundaries_in[ifit + 1, :]
            boundaries_out[ifit, :] = np.sign(ec) * np.sqrt(np.abs(ec))
            boundaries_out[ifit + 1, :] = np.sign(es) * np.sqrt(np.abs(es))

    return boundaries_out


# ==============================================================================


def e_to_sqrte_fitting(fitting_in, names_par):
    nfit = np.shape(fitting_in)[0]
    fitting_out = np.array(fitting_in, dtype=np.float64).copy()
    for ifit in range(0, nfit):
        if names_par[ifit][:2] == "ec":
            ec = fitting_in[ifit]
            es = fitting_in[ifit + 1]
            ss = ec * ec + es * es
            sqrte = np.power(ss, 0.25)
            ww_rad = np.arctan2(es, ec)
            fitting_out[ifit] = sqrte * np.cos(ww_rad)
            fitting_out[ifit + 1] = sqrte * np.sin(ww_rad)

    return fitting_out


# ==============================================================================


def e_to_sqrte_flatchain(fitting_in, names_par):
    _, nfit = np.shape(fitting_in)
    fitting_out = np.array(fitting_in, dtype=np.float64).copy()
    for ifit in range(0, nfit):
        if names_par[ifit][:2] == "ec":
            ec = fitting_in[:, ifit]
            es = fitting_in[:, ifit + 1]
            ss = ec * ec + es * es
            sqrte = np.power(ss, 0.25)
            ww_rad = np.arctan2(es, ec)
            fitting_out[:, ifit] = sqrte * np.cos(ww_rad)
            fitting_out[:, ifit + 1] = sqrte * np.sin(ww_rad)

    return fitting_out


# ==============================================================================


def e_to_sqrte_chain(fitting_in, names_par):
    _, _, nfit = np.shape(fitting_in)
    fitting_out = np.array(fitting_in, dtype=np.float64).copy()
    for ifit in range(0, nfit):
        if names_par[ifit][:2] == "ec":
            ec = fitting_in[:, :, ifit]
            es = fitting_in[:, :, ifit + 1]
            ss = ec * ec + es * es
            sqrte = np.power(ss, np.float64(0.25))
            ww_rad = np.arctan2(es, ec)
            fitting_out[:, :, ifit] = sqrte * np.cos(ww_rad)
            fitting_out[:, :, ifit + 1] = sqrte * np.sin(ww_rad)

    return fitting_out


# ==============================================================================


def e_to_sqrte_parameters(fitting_in, names_par):
    if len(np.shape(fitting_in)) == 1:  # fitting parameters
        fitting_out = e_to_sqrte_fitting(fitting_in, names_par)

    elif len(np.shape(fitting_in)) == 2:  # flatchain posterior
        fitting_out = e_to_sqrte_flatchain(fitting_in, names_par)

    elif (
        len(np.shape(fitting_in)) == 3
    ):  # full chain (original, nw x nr x nfit, transposed nr x nw x nfit)
        fitting_out = e_to_sqrte_chain(fitting_in, names_par)

    return fitting_out


# ==============================================================================


def sqrte_to_e_fitting(fitting_in, names_par):
    nfit = np.shape(fitting_in)[0]
    fitting_out = np.array(fitting_in, dtype=np.float64).copy()
    for ifit in range(0, nfit):
        if "sqrtec" in names_par[ifit] or "se" in names_par[ifit]:
            sec = fitting_in[ifit]
            ses = fitting_in[ifit + 1]
            ee = sec * sec + ses * ses
            ww_rad = np.arctan2(ses, sec)
            fitting_out[ifit] = ee * np.cos(ww_rad)
            fitting_out[ifit + 1] = ee * np.sin(ww_rad)

    return fitting_out


# ==============================================================================


def sqrte_to_e_flatchain(fitting_in, names_par):
    _, nfit = np.shape(fitting_in)
    fitting_out = np.array(fitting_in, dtype=np.float64).copy()
    for ifit in range(0, nfit):
        if "sqrtec" in names_par[ifit] or "se" in names_par[ifit]:
            sec = fitting_in[:, ifit]
            ses = fitting_in[:, ifit + 1]
            ee = sec * sec + ses * ses
            ww_rad = np.arctan2(ses, sec)
            fitting_out[:, ifit] = ee * np.cos(ww_rad)
            fitting_out[:, ifit + 1] = ee * np.sin(ww_rad)

    return fitting_out


# ==============================================================================


def sqrte_to_e_chain(fitting_in, names_par):
    _, _, nfit = np.shape(fitting_in)
    fitting_out = np.array(fitting_in, dtype=np.float64).copy()
    for ifit in range(0, nfit):
        if "sqrtec" in names_par[ifit] or "se" in names_par[ifit]:
            sec = fitting_in[:, :, ifit]
            ses = fitting_in[:, :, ifit + 1]
            ee = sec * sec + ses * ses
            ww_rad = np.arctan2(ses, sec)
            fitting_out[:, :, ifit] = ee * np.cos(ww_rad)
            fitting_out[:, :, ifit + 1] = ee * np.sin(ww_rad)

    return fitting_out


# ==============================================================================


def sqrte_to_e_parameters(fitting_in, names_par):
    if len(np.shape(fitting_in)) == 1:  # fitting parameters
        fitting_out = sqrte_to_e_fitting(fitting_in, names_par)

    elif len(np.shape(fitting_in)) == 2:  # flatchain posterior
        fitting_out = sqrte_to_e_flatchain(fitting_in, names_par)

    elif (
        len(np.shape(fitting_in)) == 3
    ):  # full chain (original, nw x nr x nfit, transposed nr x nw x nfit)
        fitting_out = sqrte_to_e_chain(fitting_in, names_par)

    return fitting_out


# ==============================================================================
def check_wrapped_parameters(names_par):

    nfit = len(names_par)
    wrapped = [False] * nfit

    for ifit, name in enumerate(names_par):
        if (
            # ("cos" in name)
            # or ("sin" in name)
            # or
            ("lambda" in name)
            or (name[0] == "w")
            or ("mA" in name)
        ):
            wrapped[ifit] = True

    return wrapped


# ==============================================================================

# ==============================================================================

# ==============================================================================

# RV semi-amplitude K in m/s
def compute_Kms(Ms_sun, Mp_jup, inc_deg, P_day, ecc):

    Ms_jup = Ms_sun * cst.Msjup
    sini = np.sin(inc_deg * cst.deg2rad)
    G_m_mj_s = cst.Gsi * cst.Mjup
    P_sec = P_day * cst.day2sec

    P_factor = np.power((cst.dpi * G_m_mj_s / P_sec), 1.0 / 3.0)
    M_factor = (Mp_jup * sini) / np.power((Ms_jup + Mp_jup), 2.0 / 3.0)
    e_factor = 1.0 / np.sqrt(1.0 - (ecc * ecc))

    Kms = P_factor * M_factor * e_factor

    return Kms


# ==============================================================================

# ==============================================================================

# epoch or transit number for each T0 given a T0ref and a Pref
def calculate_epoch(T0, Tref, Pref):

    epo = np.rint((T0 - Tref) / Pref).astype(int)

    return epo


# ==============================================================================

# 2) linear fit function with errors on y-axis.
# It returns also the errors.
# y = b_ls + m_ls * x
def lstsq_fit(x, y, yerr):

    A = np.vstack((np.ones_like(x), x)).T
    C = np.diag(yerr * yerr)
    cov = np.linalg.inv(np.dot(A.T, np.linalg.solve(C, A)))
    b_ls, m_ls = np.dot(cov, np.dot(A.T, np.linalg.solve(C, y)))
    coeff_err = np.sqrt(np.diag(cov))

    return b_ls, m_ls, coeff_err


# ==============================================================================


def compute_lin_ephem_lstsq(T0, eT0):

    nT0 = np.shape(T0)[0]
    # Tref0 = T0[0]
    Tref0 = T0[nT0 // 2]
    # dT = [np.abs(T0[i + 1] - T0[i]) for i in range(nT0 - 1)]
    dT = np.diff(T0)
    Pref0 = np.percentile(
        np.array(dT),
        50.0,
        #
    )

    epo = calculate_epoch(T0, Tref0, Pref0)
    Tref, Pref, _ = lstsq_fit(epo, T0, eT0)
    epo = calculate_epoch(T0, Tref, Pref)

    return epo, Tref, Pref


# ==============================================================================


def linear_model(par, x):

    linmod = par[0] + x * par[1]

    return linmod


# ==============================================================================


def linear_model_curve_fit(x, q, m):

    linmod = q + x * m

    return linmod


# ==============================================================================


def res_linear_model(par, x, y, ey=None):

    linmod = linear_model(par, x)
    if ey is not None:
        wres = (y - linmod) / ey
    else:
        wres = y - linmod

    return wres


# ==============================================================================


def chi2r_linear_model(par, x, y, ey=None):

    wres = res_linear_model(par, x, y, ey)
    nfit = np.shape(par)[0]
    ndata = np.shape(x)[0]
    dof = ndata - nfit
    chi2r = np.sum(np.power(wres, 2)) / np.float64(dof)

    return chi2r


# ==============================================================================


# def compute_lin_ephem(T0, eT0=None, epoin=None, modefit="wls"):

#     nT0 = np.shape(T0)[0]

#     if eT0 is None:
#         errTT = np.ones((nT0)) / 86400.0
#     else:
#         errTT = eT0

#     if epoin is None:
#         Tref0 = T0[nT0 // 2]  # T0[int(0.5*nT0)]
#         dT = np.diff(T0)  # [np.abs(T0[i+1]-T0[i]) for i in range(nT0-1)]
#         Pref0 = np.min(np.abs(dT))
#         epo = calculate_epoch(T0, Tref0, Pref0)
#     else:
#         epo = epoin
#         Tref0, Pref0, _ = lstsq_fit(epo, T0, errTT)

#     if modefit in ["optimize", "minimize"]:
#         # SCIPY.OPTIMIZE.MINIMIZE
#         optres = sciopt.minimize(
#             chi2r_linear_model,
#             [Tref0, Pref0],
#             method="nelder-mead",
#             args=(epo, T0, errTT),
#         )
#         Tref, Pref = optres.x[0], optres.x[1]
#         TP_err = [0.0, 0.0]
#         epo = calculate_epoch(T0, Tref, Pref)

#     elif modefit == "curve_fit":
#         optres = sciopt.curve_fit(linear_model_curve_fit, xdata=epo, ydata=T0)[0]
#         Tref, Pref = optres[0], optres[1]
#         TP_err = [0.0, 0.0]
#         epo = calculate_epoch(T0, Tref, Pref)

#     elif modefit in ["sklearn", "linear_model"]:
#         # SKLEARN.LINEAR_MODEL
#         sk_mod = sklm.LinearRegression()
#         sk_mod.fit(epo.reshape((-1, 1)), T0)
#         Tref, Pref = np.asscalar(np.array(sk_mod.intercept_)), np.asscalar(
#             np.array(sk_mod.coef_)
#         )
#         epo = calculate_epoch(T0, Tref, Pref)

#     elif modefit == "odr":
#         # SCIPY.ODR
#         odr_lm = sodr.Model(linear_model)
#         if eT0 is not None:
#             odr_data = sodr.RealData(epo, T0, sy=eT0)
#         else:
#             odr_data = sodr.RealData(epo, T0)
#         init_coeff = [Tref0, Pref0]
#         odr_mod = sodr.ODR(odr_data, odr_lm, beta0=init_coeff)
#         odr_out = odr_mod.run()
#         Tref, Pref = odr_out.beta[0], odr_out.beta[1]
#         TP_err = odr_out.sd_beta

#     else:  # wls
#         X = sm.add_constant(epo)
#         if eT0 is not None:
#             wls = sm.WLS(T0, X, weights=1.0 / (eT0 * eT0)).fit()
#         else:
#             wls = sm.WLS(T0, X).fit()
#         Tref, Pref = wls.params[0], wls.params[1]
#         epo = calculate_epoch(T0, Tref, Pref)
#         TP_err = wls.bse

#     return epo, Tref, Pref, TP_err


# ==============================================================================
# ==============================================================================


def read_priors(full_path):

    prior_file = os.path.join(full_path, "priors.in")
    priors = {}
    if os.path.exists(prior_file):
        of = open(prior_file, "r")
        lines = of.readlines()
        for line_raw in lines:
            line = line_raw.strip()
            if len(line) > 0:
                if line[0] != "#":
                    l = line.split("#")[0].split()
                    key_p = l[0]
                    val_p = [np.float64(x) for x in l[1:]]
                    priors[key_p] = val_p
        of.close()

    return priors


# =============================================================================
# Computes the log-penalty due to a Gaussian (asymmetric) prior
# =============================================================================
def set_gaussian_lnpar(par, pv):  # par = fitted, pv = prior value [value, nerr, perr]

    lp = par - pv[0]
    if lp < 0.0:
        lnpar = -0.5 * (lp * lp / (pv[1] * pv[1]))
    else:
        lnpar = -0.5 * (lp * lp / (pv[2] * pv[2]))

    return lnpar


# ==============================================================================


def all_parameters_to_kep_elem(all_par, NB):

    # based on how I wrote with trades the XXX_initialElements.dat file
    # M_Msun R_Rsun P_d a_AU e w_deg mA_deg inc_deg lN_deg
    # planet id 2
    # planet id ...
    # no star line!

    kep_elem = np.zeros((NB - 1, 9))
    # # M and R star
    # kep_elem[0,0] = all_par[0]
    # kep_elem[0,1] = all_par[1]

    for i in range(1, NB):
        kep_elem[i - 1, 0] = all_par[2 + (i - 1) * 8]  # Mplanet
        kep_elem[i - 1, 1] = all_par[3 + (i - 1) * 8]  # Rplanet
        kep_elem[i - 1, 2] = all_par[4 + (i - 1) * 8]  # Pplanet
        kep_elem[i - 1, 3] = period_to_semimajoraxis(
            all_par[0], kep_elem[i - 1, 0], all_par[4 + (i - 1) * 8]
        )  # aplanet
        kep_elem[i - 1, 4] = all_par[5 + (i - 1) * 8]  # eplanet
        kep_elem[i - 1, 5] = all_par[6 + (i - 1) * 8]  # wplanet
        kep_elem[i - 1, 6] = all_par[7 + (i - 1) * 8]  # mAplanet
        kep_elem[i - 1, 7] = all_par[8 + (i - 1) * 8]  # iplanet
        kep_elem[i - 1, 8] = all_par[9 + (i - 1) * 8]  # lNplanet

    return kep_elem


# ==============================================================================


def compute_proper_sigma(nfit, delta_sigma, parameter_names):

    delta_sigma_out = np.ones((nfit)) * delta_sigma
    for ifit in range(0, nfit):
        if delta_sigma > 1.0e-6:
            if "Ms" in parameter_names[ifit]:
                delta_sigma_out[ifit] = delta_sigma * 1.0e-3
            if delta_sigma > 1.0e-2:
                if "P" in parameter_names[ifit]:
                    delta_sigma_out[ifit] = delta_sigma * 1.0e-2
                elif "mA" in parameter_names[ifit] or "lambda" in parameter_names[ifit]:
                    delta_sigma_out[ifit] = 1.0e-4
        if "gamma" in parameter_names[ifit]:
            delta_sigma_out[ifit] = 10.0

    return delta_sigma_out


# ==============================================================================


def compute_initial_walkers(
    lnprob,
    nfit,
    nwalkers,
    fitting_parameters,
    parameters_minmax,
    delta_sigma,
    parameter_names,
    args=(),
    of_run=None,
):
    print_both(
        " Inititializing walkers with delta_sigma = %s" % (str(delta_sigma).strip()),
        of_run,
    )
    p0 = []
    i_p0 = 0

    print_both(" good p0:", of_run)

    # 2017-02-03 LUCA --0--
    try:
        d_sigma = np.float64(delta_sigma)
    except:
        d_sigma = np.float64(1.0e-4)

    delta_sigma_out = compute_proper_sigma(nfit, d_sigma, parameter_names)
    print(" ", end=" ")
    # init all initial walkers
    while True:
        test_p0 = np.array(
            [
                fitting_parameters[ifit]
                + np.random.normal(loc=0.0, scale=delta_sigma_out[ifit])
                for ifit in range(0, nfit)
            ],
            dtype=np.float64,
        )
        for p_n in parameter_names:
            if ("Ms" in p_n) or ("P" in p_n):
                test_p0[parameter_names.index(p_n)] = np.abs(
                    test_p0[parameter_names.index(p_n)]
                )
        test_lg = lnprob(test_p0, *args)
        if not np.isinf(test_lg):
            i_p0 += 1
            p0.append(test_p0)
            print(i_p0, end=" ")
            if i_p0 == nwalkers:
                break
    # I want the original fitting paramameters in the initial walkers
    p0[-1] = fitting_parameters
    print()
    # if 'random' opt ==> create other Gaussian starting points (<->nwalkers)
    if "ran" in str(delta_sigma).strip().lower():
        delta_parameters = np.abs(
            parameters_minmax[:, 1] - parameters_minmax[:, 0]
        )  # DELTA BETWEEN MAX AND MIN OF BOUNDARIES
        # a new Gaussian starting point each nw_min walkers,
        # keeping at least nw_min walkers Gaussian to the original fitting parameters
        # nw_min = 30
        nw_min = nwalkers // 2
        # n_gpts = int(
        #     (nwalkers - nw_min) / nw_min
        # )
        n_gpts = nwalkers - nw_min
        print(" new gaussian starting points: ", n_gpts)
        if n_gpts > 0:
            print(" doing random-gaussian points ... ")
            for i_gpt in range(0, n_gpts):
                # create new starting point, but check if lnL != -inf
                new_start = fitting_parameters.copy()
                sel_fit = int(
                    np.random.random() * (nfit - 1)
                )  # change only parameter...
                print("gpt ", i_gpt + 1)
                print("selected sel_fit = ", sel_fit, " ==> ", parameter_names[sel_fit])
                print(
                    "val = ",
                    new_start[sel_fit],
                    " with min = ",
                    parameters_minmax[sel_fit, 0],
                    " and delta = ",
                    delta_parameters[sel_fit],
                )
                while True:
                    new_start[sel_fit] = (
                        parameters_minmax[sel_fit, 0]
                        + delta_parameters[sel_fit] * np.random.random()
                    )
                    test_lg = lnprob(new_start, *args)
                    if not np.isinf(test_lg):
                        break
                i_pos = nw_min * i_gpt
                print("i_pos = ", end=" ")
                while True:
                    test_p0 = np.array(
                        [
                            new_start[ifit]
                            + np.random.normal(loc=0.0, scale=delta_sigma_out[ifit])
                            for ifit in range(0, nfit)
                        ],
                        dtype=np.float64,
                    )
                    test_lg = lnprob(test_p0, *args)
                    if not np.isinf(test_lg):
                        p0[i_pos] = test_p0
                        print(i_pos, end=" ")
                        i_pos += 1
                        if i_pos % nw_min == 0:
                            break
            print()
        print()

    print_both(" done initial walkers.", of_run)

    return p0


# =============================================================================
def set_automatic_unit_time(val_d):

    val_h = val_d * cst.day2hour
    val_m = val_d * cst.day2min
    val_s = val_d * cst.day2sec

    if val_s < 60.0:
        scale = [86400.0, "s"]
    elif 1.0 <= val_m < 60.0:
        scale = [1440.0, "min"]
    elif 1.0 <= val_h < 24.0:
        scale = [24.0, "hours"]
    else:
        scale = [1.0, "d"]

    return scale


# =============================================================================


def de_hdf5_save_one_dataset(de_file, data_name, data, data_type, hdf5_mode="a"):

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


def de_hdf5_update_attr(de_file, data_name, attrs_name, attrs_value):

    de_hdf5 = h5py.File(
        de_file,
        "r+",
        # libver='latest'
    )
    # de_hdf5.swmr_mode = True
    de_hdf5[data_name].attrs[attrs_name] = attrs_value
    de_hdf5.close()

    return


def de_save_evolution(
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
    de_maximize=True,
):

    # de_file = os.path.join(de_path, 'de_run.hdf5')
    de_hdf5_save_one_dataset(de_file, "population", de_pop, "f8")
    de_hdf5_update_attr(de_file, "population", "npop", npop_de)
    de_hdf5_update_attr(de_file, "population", "ngen", ngen_de)
    de_hdf5_update_attr(de_file, "population", "iter_de", iter_de)
    de_hdf5_update_attr(de_file, "population", "iter_global", iter_global)
    de_hdf5_update_attr(de_file, "population", "nfit", nfit)
    de_hdf5_update_attr(de_file, "population", "ndata", ndata)

    de_hdf5_save_one_dataset(de_file, "population_fitness", de_fit, "f8")

    de_hdf5_save_one_dataset(de_file, "best_population", de_pop_best, "f8")

    de_hdf5_save_one_dataset(de_file, "best_fitness", de_fit_best, "f8")

    de_hdf5_save_one_dataset(de_file, "parameters_minmax", de_bounds, "f8")

    # best_loc = np.argmax(de_fit_best[:iter_de])
    # de_hdf5_save_one_dataset(de_file, 'de_parameters', de_pop_best[best_loc,:].copy(), 'f8')
    de_parameters = de_get_best_parameters(
        de_fit_best, de_pop_best, iter_de, de_maximize=de_maximize
    )
    de_hdf5_save_one_dataset(de_file, "de_parameters", de_parameters, "f8")

    de_hdf5_save_one_dataset(de_file, "parameter_names", parameter_names, "S10")

    return


# =============================================================================


def de_get_best_parameters(de_fit_best, de_pop_best, iter_de, de_maximize=True):

    loc_best = (
        np.argmax(de_fit_best[: iter_de + 1])
        if de_maximize
        else np.argmin(de_fit_best[: iter_de + 1])
    )

    de_parameters = de_pop_best[loc_best, :].copy()

    return de_parameters


# =============================================================================
def de_load_parameters(de_file):

    with h5py.File(
        de_file,
        mode="r",
        # libver='latest',
        swmr=True,
    ) as de_hdf5:
        de_parameters = de_hdf5["de_parameters"][...]
    # de_hdf5.close()

    return de_parameters


# =============================================================================
def compute_epoch(Tref, Pref, TTs):

    epo = np.rint((TTs - Tref) / Pref)

    return epo

# =============================================================================
def u1u2_to_q1q2(u1, u2):

    u1u2 = u1+u2
    q1 = u1u2*u1u2
    q2 = u1/(2.0*u1u2)

    return q1, q2

def q1q2_to_u1u2(q1, q2):

    sq1 = np.sqrt(q1)
    u1 = 2.0*sq1*q2
    u2 = sq1*(1.0-(2.0*q2))

    return u1, u2

# =============================================================================
def normalization_standard(x):

    xs = (x - np.mean(x))/np.std(x, ddof=1)

    return xs

def normalization_range(x):

    xs = (2.0*x - (np.amin(x)+np.amax(x))) / np.ptp(x)

    return xs

def normalization_max(x):

    xs = (x - np.amin(x))/np.ptp(x)

    return xs

def normalization_constant(x, xc=1.0):

    xs = x / xc

    return xs

# =============================================================================
def angle_to_harmonics(phi, n_harmonics=3):

    phir = phi*cst.deg2rad

    harmonics = {}
    for nh in range(1,n_harmonics+1):
        harmonics["cos{:02d}phi".format(nh)] = np.cos(nh*phir)
        harmonics["sin{:02d}phi".format(nh)] = np.sin(nh*phir)

    return harmonics