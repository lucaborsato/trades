"""Generates TRADES configuration files from a system dict."""
import ast
import json
import os
import re


PLANET_LABELS = "bcdefghijklmnopqrstuvwxyz"

# These are the eight fit/prior parameters represented by the planet columns in
# bodies.lst. The derived semi-major-axis and pericentre-time rows stay in the
# planet .dat file but do not have bodies.lst fit flags.
FIT_PARAMETERS = (
    ("mass", "m", "Mass"),
    ("radius", "R", "Radius"),
    ("period", "P", "Period"),
    ("ecc", "e", "eccentricity"),
    ("omega", "w", "argument of pericenter"),
    ("M0", "mA", "mean anomaly"),
    ("inc", "i", "inclination"),
    ("lnode", "lN", "longitude of ascending node"),
)


def get_planet_label(idx):
    """Get letter label for planet index (0-based)."""
    if idx < len(PLANET_LABELS):
        return PLANET_LABELS[idx]
    return str(idx + 2)


def get_body_id(label, n_planets):
    """Get body ID from planet label. Body 1=star, 2=first planet (b), etc."""
    if label == "star":
        return 1
    planet_idx = PLANET_LABELS.index(label) if label in PLANET_LABELS else -1
    if planet_idx >= 0 and planet_idx < n_planets:
        return planet_idx + 2
    return 2


def format_number(value):
    """Shortest round-trip representation, keeping full user precision."""
    return repr(float(value))


def _format_min_decimals(value, min_decimals):
    """User precision when it has enough digits, otherwise padded to min_decimals."""
    number = float(value)
    text = repr(number)
    if "e" in text or "E" in text:
        return f"{number:.{min_decimals}f}"
    decimals = text.partition(".")[2]
    if len(decimals) < min_decimals:
        return f"{number:.{min_decimals}f}"
    return text


def format_time(value):
    """Times keep user precision, with at least 7 decimal digits."""
    return _format_min_decimals(value, 7)


def format_rv(value):
    """RV values keep user precision, with at least 3 decimal digits."""
    return _format_min_decimals(value, 3)


def generate_bodies_lst(n_planets, star_fit_mass=False, star_fit_radius=False, planets=None):
    """Generate bodies.lst content."""
    lines = []
    lines.append(f"star.dat {int(bool(star_fit_mass))} {int(bool(star_fit_radius))}          #filename Mass Radius [1=to fit, 0=fixed]")
    for i in range(n_planets):
        label = get_planet_label(i)
        planet = planets[i] if planets and i < len(planets) else {}
        flags = " ".join(str(int(bool(planet.get(f"fit_{key}", True)))) for key, _, _ in FIT_PARAMETERS)
        transit_flag = "T" if planet.get("transiting", True) else "F"
        lines.append(f"{planet.get('label', label)}.dat {flags} {transit_flag}          #filename Mass Radius Period eccentricity Arg.ofPericenter MeanAnomaly inclination Long.ofNodes transit_flag(F or f if it should not transit)")
    return "\n".join(lines) + "\n"


def generate_star_dat(star_mass, star_mass_err, star_radius, star_radius_err, transits):
    """Generate the two-line star.dat file."""
    lines = []
    lines.append(f"{format_number(star_mass)} {format_number(star_mass_err)}          # Mstar [Msun]")
    lines.append(f"{format_number(star_radius)} {format_number(star_radius_err)}          # Rstar [Rsun]")
    return "\n".join(lines) + "\n"


def generate_planet_dat(p):
    """Generate N.dat content for a planet."""
    lines = []
    rows = [
        ("Mass", p["mass"], p["mass_min"], p["mass_max"], "[Mjup]"),
        ("Radius", p["radius"], p["radius_min"], p["radius_max"], "[Rjup]"),
        ("Period", p["period"], p["period_min"], p["period_max"], "[day]"),
        ("semi major axis", p["sma"], p["sma_min"], p["sma_max"], "[AU]"),
        ("eccentricity", p["ecc"], p["ecc_min"], p["ecc_max"], ""),
        ("argument of the pericenter", p["omega"], p["omega_min"], p["omega_max"], "[deg]"),
        ("mean anomaly", p["M0"], p["M0_min"], p["M0_max"], "[deg]"),
        ("time of pericenter passage", p["tp"], p["tp_min"], p["tp_max"], "[JD]"),
        ("orbit inclination", p["inc"], p["inc_min"], p["inc_max"], "[deg]"),
        ("longitude of the ascending node", p["lnode"], p["lnode_min"], p["lnode_max"], "[deg]"),
    ]
    for name, val, vmin, vmax, unit in rows:
        u = f" {unit}" if unit else ""
        lines.append(f"{format_number(val)} {format_number(vmin)} {format_number(vmax)}      # {name}{u}")
    return "\n".join(lines) + "\n"


def generate_transit_observations(obs_list, include_duration=True):
    """Generate an NB*_observations.dat file with one header row."""
    lines = []
    if include_duration:
        lines.append("# Epoch Transit_time(days) err_Transit_time(days) Duration(min) err_Duration(min) LC_id")
        for o in obs_list:
            lines.append(f"{int(o['epoch'])} {format_time(o['tt'])} {format_time(o['ett'])} {format_time(o.get('dur', 0.0))} {format_time(o.get('edur', 0.0))} {o.get('lc_id', 1)}")
    else:
        lines.append("# Epoch Transit_time(days) err_Transit_time(days) LC_id")
        for o in obs_list:
            lines.append(f"{int(o['epoch'])} {format_time(o['tt'])} {format_time(o['ett'])} {o.get('lc_id', 1)}")
    return "\n".join(lines) + "\n"


def generate_rv_observations(obs_list):
    """Generate obsRV.dat content with one header row."""
    lines = ["# Time_RV(days) RV(m/s) err_RV(m/s) RV_set"]
    for o in obs_list:
        lines.append(f"{format_time(o['jd'])}  {format_rv(o['rv'])}  {format_rv(o['erv'])} {o['rvset']}")
    return "\n".join(lines) + "\n"


# configuration.yml parameter definitions: (key, type, default, comment, choices).
# Types: text, int, float, bool, choice. Comments come from the base_2p example
# and are reported both in the GUI and in the generated file.
RUN_PARAMS = [
    ("full_path", "text", ".", None, None),
    ("sub_folder", "text", "test_emcee", None, None),
    ("seed", "int", 42, None, None),
    ("mass_type", "choice", "e", "suggested: j, e, n", ["j", "e", "n"]),
    ("delta_sigma", "float", 1.0e-4, None, None),
    ("trades_previous", "text", None, "/full/path/to/run/trades/finalXpar.dat", None),
    ("nthreads", "int", 34, None, None),
]

PYDE_PARAMS = [
    ("type", "choice", "run", "False/no, run, resume, to_emcee", ["no", "run", "resume", "to_emcee"]),
    ("npop", "int", 68, "population size == number of different configurations for each iteration of the DE; if it matches nwalkers is better", None),
    ("ngen", "int", 42000, "generation size == number of iterations", None),
    ("save", "int", 4200, "saving steps == number of steps for the save points ==> it creates a de_run.hdf5 file", None),
    ("f", "float", 0.5, "Difference amplification factor. Values between 0.5-0.8 are good in most cases.", None),
    ("c", "float", 0.5, "Cross-over probability. Use 0.9 to test for fast convergence, and smaller values (~0.1) for a more elaborate search.", None),
    ("maximize", "bool", True, "maximize == True: maximize the function, maximize == False: minimize the function", None),
]

EMCEE_PARAMS = [
    ("nwalkers", "int", 68, None, None),
    ("nruns", "int", 2000, None, None),
    ("thin_by", "int", 100, "!! the true nruns will be nruns x thin_by, this is also the keyword used by emcee and it will keep only thin_by steps", None),
    ("emcee_restart", "bool", False, None, None),
    ("emcee_progress", "bool", True, None, None),
    ("pre_optimise", "bool", True, "run a local optimiser before passing initial configuration to emcee", None),
    ("move_type", "text", "de, desnooker", "move types: ai == affine invariant, de == differential evolution, desnooker == de with different implementation", None),
    ("move_fraction", "text", "0.8, 0.2", "fractions associated to each move type; they have to sum to 1", None),
]

ANALYSIS_PARAMS = [
    ("full_path", "text", "", "empty = auto from run.sub_folder", None),
    ("m_type", "choice", "e", "suggested: j, e, n", ["j", "e", "n"]),
    ("nburnin", "int", 100, "number of steps to discard as burn-in", None),
    ("emcee_old_save", "bool", False, None, None),
    ("temp_status", "bool", False, "True if you are analysing a temporary emcee run", None),
    ("use_thin", "int", 1, "remember if you used thin_by > 1 or not", None),
    ("seed", "int", 42, None, None),
    ("from_file", "text", None, "trades file of a simulation output to simulate and plot. Default = None.", None),
    ("n_samples", "int", 42, "set to 0 if you don't want to simulate and save samples --> samples_ttra_rv.hdf5", None),
    ("overplot", "choice", "map_hdi", "initial, pso, de, median, map_hdi, map, mode, adhoc", ["initial", "pso", "de", "median", "map_hdi", "map", "mode", "adhoc"]),
    ("corner_type", "choice", "pygtc", "pygtc, custom, both", ["pygtc", "custom", "both"]),
    ("lnProb_selection", "text", "[None, None]", "further selection on ln-probability, default [None, None] = [-np.inf, np.inf]", None),
    ("all_analysis", "bool", True, "if set to False you should set by hand the following arguments, otherwise it run all of them", None),
    ("save_posterior", "bool", True, None, None),
    ("save_parameters", "bool", True, None, None),
    ("chain", "bool", True, None, None),
    ("gelman_rubin", "bool", True, None, None),
    ("gr_steps", "int", 10, None, None),
    ("geweke", "bool", True, None, None),
    ("gk_steps", "int", 10, None, None),
    ("correlation_fitted", "bool", True, None, None),
    ("correlation_physical", "bool", True, None, None),
]

OC_PARAMS = [
    ("plot_oc", "bool", True, None, None),
    ("full_path", "text", "", "empty = auto from run.sub_folder", None),
    ("sim_name", "text", "de, initial, median, map_hdi, map", "list of [initial, pso, de, median, map_hdi, map, mode, from_file]", None),
    ("idplanet_name", "text", "{}", "e.g {2: 'b', 3: 'c'}", None),
    ("lmflag", "int", 0, "if LM has been run (1) or not (0). Default 0.0", None),
    ("tscale", "text", None, "time scale to remove from x-axis plot; suggested [value, '.Xf']; Default None means removing nothing, like [0.0, '.0f']", None),
    ("plot_title", "bool", True, None, None),
    ("unit", "choice", "auto", "unit of the O-C: auto-detect best unit for each body, or d, h, m, s", ["auto", "d", "h", "m", "s"]),
    ("samples_file", "text", None, "hdf5 file name with samples created by analysis/n_samples > 0: i.e. samples_ttra_rv.hdf5", None),
    ("plot_obs_sim", "bool", False, "plot simulations point corrisponding to observations? True/False. Default=True.", None),
    ("plot_samples", "choice", "ci", "how to plot samples: 'ci' == 'sigma', or 'all' (each sample plotted)", ["ci", "all"]),
    ("limits", "choice", "obs", "x-y limits based on observations (obs) or samples (sam)", ["obs", "sam"]),
    ("kep_ele", "bool", False, "plot of Keplerian Elements of each transiting body during transits", None),
    ("legend", "choice", "in", "in or out, everything else means no legend", ["in", "out"]),
    ("linear_ephemeris", "text", None, "default None: compute the linear ephemeris from data for each body; otherwise provide per body as {2: {Tref: [val, err], Pref: [val, err]} }", None),
    ("idsource_name", "text", None, "e.g {1: 'TESS', 2: 'CHEOPS'} to map sub-data set with telescope -> different color and marker for each data-set", None),
    ("color_map", "text", "nipy_spectral", "or as {1: 'C0', 2: 'C1'} where 1, 2, etc match the ids in idsource_name", None),
]

RV_PARAMS = [
    ("plot_rv", "bool", True, None, None),
    ("full_path", "text", "", "empty = auto from run.sub_folder", None),
    ("sim_name", "text", "de, initial, median, map_hdi, map", "list of [initial, pso, de, median, map_hdi, map, mode, from_file]", None),
    ("lmflag", "int", 0, "if LM has been run (1) or not (0). Default 0.0", None),
    ("tscale", "text", None, "time scale to remove from x-axis plot; suggested [value, '.Xf']; Default None means removing nothing, like [0.0, '.0f']", None),
    ("plot_title", "bool", True, None, None),
    ("samples_file", "text", None, "hdf5 file name with samples created by analysis/n_samples > 0: i.e. samples_ttra_rv.hdf5", None),
    ("limits", "choice", "obs", "x-y limits based on observations (obs) or samples (sam)", ["obs", "sam"]),
    ("legend", "choice", "in", "in or out, everything else means no legend", ["in", "out"]),
    ("labels", "text", "RVdataset#1, RVdataset#2", "list of labels of different RV dataset (order has to match the numbering of obsRV.dat file) or dictionary: {1: 'RVdataset#1', 2: 'RVdataset#2'}", None),
    ("color_map", "text", "nipy_spectral", "color map name, color name or {1: 'C0', 2: 'C1'} matching the RV set ids", None),
]

DEFAULT_CONFIG = {}
for _prefix, _params in [
    ("run", RUN_PARAMS), ("pyde", PYDE_PARAMS), ("emcee", EMCEE_PARAMS),
    ("analysis", ANALYSIS_PARAMS), ("oc", OC_PARAMS), ("rv", RV_PARAMS),
]:
    for _key, _ptype, _default, _comment, _choices in _params:
        DEFAULT_CONFIG[f"{_prefix}.{_key}"] = _default


def _parse_scalar_text(value):
    """Parse a free-text YAML-ish value: empty -> None, brackets -> list/dict."""
    if value is None:
        return None
    text = str(value).strip()
    if not text or text.lower() in ("none", "null"):
        return None
    if (text.startswith("[") and text.endswith("]")) or (text.startswith("{") and text.endswith("}")):
        try:
            return ast.literal_eval(text)
        except (ValueError, SyntaxError):
            return text
    return text


def _split_list(value, cast=str):
    """Parse a comma-separated text field into a list."""
    if value is None:
        return None
    if isinstance(value, (list, tuple)):
        return [cast(item) for item in value]
    text = str(value).strip()
    if not text or text.lower() in ("none", "null"):
        return None
    items = [item.strip() for item in text.split(",") if item.strip()]
    return [cast(item) for item in items] if items else None


def _yaml_scalar(value):
    """Format a python value as a YAML scalar."""
    if value is None:
        return "None"
    if isinstance(value, bool):
        return "True" if value else "False"
    if isinstance(value, float):
        return format_number(value)
    if isinstance(value, int):
        return str(value)
    if isinstance(value, (list, tuple)):
        return "[" + ", ".join(_yaml_scalar(item) for item in value) + "]"
    if isinstance(value, dict):
        return "{" + ", ".join(f"{_yaml_scalar(k)}: {_yaml_scalar(v)}" for k, v in value.items()) + "}"
    text = str(value).strip()
    if re.match(r"^[A-Za-z_][A-Za-z0-9_.\-]*$", text) and text.lower() not in (
        "true", "false", "null", "none", "nan", "yes", "no", "on", "off", "y", "n",
    ):
        return text
    return json.dumps(text)


def _yaml_line(indent, key, value, comment=None):
    line = f"{' ' * indent}{key}: {value}"
    if comment:
        line += f" # {comment}"
    return line


def generate_configuration_yml(system):
    """Generate configuration.yml content including the original comments."""
    cfg = dict(DEFAULT_CONFIG)
    cfg.update(system.get("config", {}))

    def value_for(prefix, key, ptype):
        value = cfg.get(f"{prefix}.{key}", DEFAULT_CONFIG.get(f"{prefix}.{key}"))
        if ptype != "text":
            return value
        if key == "full_path" and prefix in ("analysis", "oc", "rv") and (value is None or not str(value).strip()):
            sub = str(cfg.get("run.sub_folder", "test_emcee") or "").strip()
            return f"./{sub}" if sub and sub.lower() != "none" else None
        if key in ("sim_name", "labels"):
            return _split_list(value, str)
        return _parse_scalar_text(value)

    lines = []
    lines.append("# RUN SECTION e.g. for EMCEE, DE+EMCEE, etc.")
    lines.append("run:")
    for key, ptype, _default, comment, _choices in RUN_PARAMS:
        lines.append(_yaml_line(4, key, _yaml_scalar(value_for("run", key, ptype)), comment))
    lines.append("    pyde:")
    for key, ptype, _default, comment, _choices in PYDE_PARAMS:
        lines.append(_yaml_line(8, key, _yaml_scalar(value_for("pyde", key, ptype)), comment))
    lines.append("    emcee:")
    for key, ptype, _default, comment, _choices in EMCEE_PARAMS:
        if key.startswith("move_"):
            continue
        lines.append(_yaml_line(8, key, _yaml_scalar(value_for("emcee", key, ptype)), comment))
    move_comment = "where ai == affine invariant, de == differential evolution, desnooker == de with different implementation"
    move_types = _split_list(cfg.get("emcee.move_type"), str) or ["de", "desnooker"]
    move_fractions = _split_list(cfg.get("emcee.move_fraction"), float) or [0.8, 0.2]
    lines.append("        move:")
    lines.append(f"            type: {_yaml_scalar(move_types)} # {move_comment}")
    lines.append(f"            fraction: {_yaml_scalar(move_fractions)}")
    lines.append("# ANALYSIS SECTION, TO RUN AFTER EMCEE SAVED AT LEAST ONE SAVE POINT, YOU CAN RUN IT ALSO FOR DE/PSO OUTPUT")
    lines.append("analysis:")
    for key, ptype, _default, comment, _choices in ANALYSIS_PARAMS:
        lines.append(_yaml_line(4, key, _yaml_scalar(value_for("analysis", key, ptype)), comment))
    lines.append("# O-C PLOTS, if missing the SECTION NO PLOT AT ALL")
    lines.append("OC:")
    for key, ptype, _default, comment, _choices in OC_PARAMS:
        lines.append(_yaml_line(4, key, _yaml_scalar(value_for("oc", key, ptype)), comment))
    lines.append("RV:")
    for key, ptype, _default, comment, _choices in RV_PARAMS:
        lines.append(_yaml_line(4, key, _yaml_scalar(value_for("rv", key, ptype)), comment))
    return "\n".join(lines) + "\n"


# arg.in parameter definitions: (key, type, default, comment, choices).
# The comment is the description line of the original arg.in file and is
# written above each parameter in the generated file.
ARG_PARAMETERS = [
    ("progtype", "choice", 2, "1=grid search, 2=integration/Levenberg-Marquardt, 3=PIKAIA (GA), 4=PSO, 5=PolyChord.", [1, 2, 3, 4, 5]),
    ("nboot", "int", 0, "bootstrap: <=0 no bootstrap, >0 yes bootstrap (Nboot set to 100 if <100)", None),
    ("bootstrap_scaling", "bool", True, "bootstrap_scaling = .true. or T, .false. or F, default .false.", None),
    ("tepoch", "time", 0.0, "epoch of the elements [JD].", None),
    ("tstart", "time", 0.0, "time start of the integration [JD]. Set this to a value > 9.e7 if you want to use the tepoch (next parameter)", None),
    ("tint", "time", 3652.5, "time duration of the integration in days", None),
    ("step", "time", 1.0e-3, "initial time step size in days", None),
    ("wrttime", "time", 0.04167, "time interval in days of write data in files (if < stepsize the program will write every step)", None),
    ("idtra", "int", 1, "number of the body to check if it transits (from 2 to N, 1 for everyone, 0 for no check).", None),
    ("durcheck", "int", 0, "0/1=no/yes duration fit", None),
    ("tol_int", "float", 1.0e-13, "tolerance in the integration (stepsize selection etc)", None),
    ("wrtorb", "int", 0, "write orbit condition: 1[write], 0[do not write]", None),
    ("wrtconst", "int", 0, "write constants condition: 1[write], 0[do not write]", None),
    ("wrtel", "int", 0, "write orbital elements condition: 1[write], 0[do not write]", None),
    ("rvcheck", "int", 0, "check of Radial Velocities condition: 1[check, read from obsRV.dat], 0[do not check]", None),
    ("rv_res_gls", "bool", False, "check of Radial Velocities residuals with GLS (look for added signals close to planetary periods): Set as F or T. Default is False: F.", None),
    ("rv_trend_order", "int", 0, "if you want to fit a trend to RV define the order (default = 0)", None),
    ("idpert", "int", 3, "grid option: id of the perturber body [integer >1, <= tot bodies; else no perturber]", None),
    ("lmon", "int", 0, "lmon: Levenberg-Marquardt off = 0[no LM], on=1[yes LM]. Default lmon = 0", None),
    ("secondary_parameters", "int", 0, "secondary_parameters: define if the program has to check only boundaries for derived parameters (1) or if it has also to fix values due to derived parameters (2) or do nothing (0) with derived parameters. Default secondary_parameters = 0.", None),
    ("close_encounter_check", "bool", True, "close encounters", None),
    ("do_hill_check", "bool", False, "Hill check", None),
    ("amd_hill_check", "bool", False, "AMD Hill check", None),
    ("ncpu", "int", 1, "number of cpu to use with opemMP. Default is 1.", None),
]

NB_COMMENT = "number of bodies to use (from 2 to N, where the max N is the number of files in bodies.lst)"


def _arg_scalar(value, ptype):
    """Format one arg.in value according to its declared type."""
    if ptype == "bool":
        return "T" if value in (True, 1, "T", "t", "true", "True") else "F"
    if ptype == "time":
        return format_time(value)
    if ptype == "float":
        return format_number(value)
    return str(int(value))


def generate_arg_in(args):
    """Generate arg.in content: comment line above each parameter."""
    lines = []
    for key, ptype, default, comment, _choices in ARG_PARAMETERS:
        lines.append(f"# {comment}")
        lines.append(f"{key} = {_arg_scalar(args.get(key, default), ptype)}")
    lines.append(f"# {NB_COMMENT}")
    lines.append(f"NB = {args.get('NB', args.get('n_planets', 3))}")
    return "\n".join(lines) + "\n"


def generate_priors(planets, star=None):
    """Generate priors.in content."""
    lines = [
        "# priors of physical parameters",
        "# label val -1sigma +1sigma",
        "# possible label: mX, PX, eX, wX, mAX, iX, lNX, jitter_N, gamma_N, c_N",
        "# X the id number of the body, where the first planet has X == 2, and last X = NB: X = 2, 3, 4, ..., NB",
        "# m1 : stellar mass in Msun; R1 : stellar radius in Rsun",
        "# mX : mass of planet X in Mearth!",
        "# PX : period of planet X in days",
        "# angles wX, mAX, iX, lNX in degree",
        "# jitter_N in m/s, gamma_N in m/s, c_N between -1 and 1.",
        "# example:",
        "# m2 1.0 -0.3 +0.3",
    ]
    star = star or {}
    if star.get("mass_prior"):
        sigma = float(star["mass_err"])
        lines.append(f"m1 {format_number(star['mass'])} {format_number(-sigma)} {format_number(sigma)}")
    if star.get("radius_prior"):
        sigma = float(star["radius_err"])
        lines.append(f"R1 {format_number(star['radius'])} {format_number(-sigma)} {format_number(sigma)}")

    for i, p in enumerate(planets):
        body_id = i + 2
        for key, prior_label, _ in FIT_PARAMETERS:
            if p.get(f"prior_{key}"):
                value = float(p.get(f"prior_{key}_value", p[key]))
                sigma = float(p[f"prior_{key}_sigma"])
                lines.append(f"{prior_label}{body_id} {format_number(value)} {format_number(-sigma)} {format_number(sigma)}")
    return "\n".join(lines) + "\n"


def generate_lm_opt(lm):
    """Generate lm.opt content."""
    lines = [
        f"{lm.get('max_eval', -1)}    # max function evaluation (if <=0 it will be set to 200*(n+1), where n=number of parameters to be fitted)",
        f"{lm.get('ftol', 1e-8)}    # ftol (must be > 0., if <=0 it will be setted with an internal TOLERANCE) DETERMINED IN THE PROGRAM",
        f"{lm.get('xtol', 1e-8)}    # xtol (must be > 0., if <=0 it will be setted with an internal TOLERANCE) DETERMINED IN THE PROGRAM",
        f"{lm.get('gtol', 1e-8)}    # gtol (must be >= 0., if <0 it will be setted to zero)",
        f"{lm.get('step_factor', -1.0)}    # value that multiply the parameters to detect a parameter step in lmdif (usually to be setted to sqrt( machine precision ), if <0 it will be set to this value) DETERMINED IN THE PROGRAM",
        f"{lm.get('nprint', 0)}    # nprint: it prints some results every nprint function evaluations (if <0 it will be set to 0)",
    ]
    return "\n".join(lines) + "\n"


def generate_pikaia_opt(pk):
    """Generate pikaia.opt content."""
    lines = [
        f"{pk.get('nindiv', 10.)}     # ctrl(1) = number of individuals",
        f"{pk.get('ngen', 10)}     # ctrl(2) = number of generations",
        f"{pk.get('ndigits', 7.)}      # ctrl(3) = number of significant digits",
        f"{pk.get('cross', 0.85)}    # ctrl(4) = crossover probability; must be  <= 1.0 (default is 0.85)",
        f"{pk.get('mut_mode', 2)}      # ctrl(5) = mutation mode; 1/2=steady/variable (default is 2)",
        f"{pk.get('mut_init', 0.02)}    # ctrl(6) = initial mutation rate; should be small (default is 0.005) (Note: the mutation rate is the probability that any one gene locus will mutate in any one generation.)",
        f"{pk.get('mut_min', 0.0005)}  # ctrl(7) = minimum mutation rate; must be >= 0.0 (default is 0.0005)",
        f"{pk.get('mut_max', 0.3)}     # ctrl(8) = maximum mutation rate; must be <= 1.0 (default is 0.25)",
        f"{pk.get('rel_diff', 1.)}      # ctrl(9) = relative fitness differential; range from 0 (none) to 1 (maximum).  (default is 1.)",
        f"{pk.get('reprod', 1)}        # ctrl(10) = reproduction plan; 1/2/3=Full generational replacement/Steady-state-replace-random/Steady-state-replace-worst (default is 3)",
        f"{pk.get('elitism', 0)}       # ctrl(11) = elitism flag; 0/1=off/on (default is 0) (Applies only to reproduction plans 1 and 2)",
        f"{pk.get('output', 0)}        # ctrl(12) = printed output 0/1/2=None/Minimal/Verbose (default is 0)",
        f"{pk.get('seed', 123456)}  # seed",
        f"{pk.get('wrtAll', 0)}       # wrtAll = 0 [not writing all individuals for each iteration] 1 [writing all individuals for each iteration]",
        f"{pk.get('nGlobal', 1)}        # nGlobal = number of global search, number of times that GA will be run, i.e., Nindividual x Ngeneration x nGlobal",
    ]
    return "\n".join(lines) + "\n"


def generate_pso_opt(ps):
    """Generate pso.opt content."""
    lines = [
        f"{ps.get('nparticles', 10)}     # number of particle",
        f"{ps.get('niter', 10)}     # number of iterations to do",
        f"{ps.get('screen', 0)}      # counter to write to screen PSO run...if 0 = no write",
        f"{ps.get('wrtAll', 0)}      # wrtAll = 0 [not writing all individuals for each iteration] 1 [writing all individuals for each iteration]",
        f"{ps.get('nGlobal', 1)}      # nGlobal = number of global search, number of times that PSO will be run, i.e., Npart x Niter x nGlobal",
        f"{ps.get('seed', 123456)}  # seed",
        f"{ps.get('inertia', 0.9)}     # inertia = inertia parameter. Reccomended 0.9",
        f"{ps.get('self_param', 2.0)}     # self = self intention parameter. Reccomended 2.0",
        f"{ps.get('swarm_param', 2.0)}     # swarm = swarm intention parameter. Reccomended 2.0",
        f"{ps.get('randsearch', 1e-5)}  # randsearch = random search parameter. Reccomended very small 1.0e-5",
        f"{ps.get('vmax', 0.5)}     # vmax = limit of vector length. Reccomended between 0.5 and 1.0",
        f"{ps.get('vrand', 0.07)}    # vrand = velocity perturbation parameter. Reccomended between 0.0 and 0.1",
    ]
    return "\n".join(lines) + "\n"


def write_system(base_dir, system):
    """Write all TRADES files to base_dir. Returns dict of created files."""
    os.makedirs(base_dir, exist_ok=True)
    n_planets = system.get("n_planets", 2)
    planets = system.get("planets", [])
    star = system.get("star", {})
    nb = n_planets + 1

    created = {}

    # bodies.lst
    path = os.path.join(base_dir, "bodies.lst")
    with open(path, "w") as f:
        f.write(generate_bodies_lst(
            n_planets,
            star.get("fit_mass", False),
            star.get("fit_radius", False),
            planets,
        ))
    created["bodies.lst"] = path

    # star.dat
    transits = system.get("transits", [])
    path = os.path.join(base_dir, "star.dat")
    with open(path, "w") as f:
        f.write(generate_star_dat(
            star.get("mass", 1.0), star.get("mass_err", 0.01),
            star.get("radius", 1.0), star.get("radius_err", 0.01),
            transits
        ))
    created["star.dat"] = path

    # N.dat for each planet
    for i, p in enumerate(planets):
        label = p.get("label", get_planet_label(i))
        path = os.path.join(base_dir, f"{label}.dat")
        with open(path, "w") as f:
            f.write(generate_planet_dat(p))
        created[f"{label}.dat"] = path

    # Transit observations for each transiting planet with at least one observation
    for i, p in enumerate(planets):
        if p.get("transiting") and p.get("transit_obs"):
            body_id = i + 2
            fname = f"NB{body_id}_observations.dat"
            path = os.path.join(base_dir, fname)
            with open(path, "w") as f:
                f.write(generate_transit_observations(
                    p.get("transit_obs", []), p.get("include_duration", True)
                ))
            created[fname] = path

    # RV observations
    rv_obs = system.get("rv_obs", [])
    if rv_obs:
        path = os.path.join(base_dir, "obsRV.dat")
        with open(path, "w") as f:
            f.write(generate_rv_observations(rv_obs))
        created["obsRV.dat"] = path

    # configuration.yml
    path = os.path.join(base_dir, "configuration.yml")
    with open(path, "w") as f:
        f.write(generate_configuration_yml(system))
    created["configuration.yml"] = path

    # arg.in
    args = system.get("arg_in", {})
    args["NB"] = nb
    args["n_planets"] = n_planets
    path = os.path.join(base_dir, "arg.in")
    with open(path, "w") as f:
        f.write(generate_arg_in(args))
    created["arg.in"] = path

    # priors.in
    path = os.path.join(base_dir, "priors.in")
    with open(path, "w") as f:
        f.write(generate_priors(planets, star))
    created["priors.in"] = path

    # lm.opt
    path = os.path.join(base_dir, "lm.opt")
    with open(path, "w") as f:
        f.write(generate_lm_opt(system.get("lm_opt", {})))
    created["lm.opt"] = path

    # pikaia.opt
    path = os.path.join(base_dir, "pikaia.opt")
    with open(path, "w") as f:
        f.write(generate_pikaia_opt(system.get("pikaia_opt", {})))
    created["pikaia.opt"] = path

    # pso.opt
    path = os.path.join(base_dir, "pso.opt")
    with open(path, "w") as f:
        f.write(generate_pso_opt(system.get("pso_opt", {})))
    created["pso.opt"] = path

    return created


def load_example(system_dir):
    """Load an existing TRADES system directory into a system dict."""
    import shutil
    system = {
        "n_planets": 0,
        "star": {"mass": 1.0, "mass_err": 0.01, "radius": 1.0, "radius_err": 0.01, "fixed_mass": True, "fixed_radius": True},
        "planets": [],
        "transits": [],
        "transit_obs": [],
        "rv_obs": [],
        "run_full_path": ".",
        "run_sub_folder": "test_emcee",
        "run_seed": 42,
        "mass_type": "e",
        "nthreads": 34,
        "pyde_npop": 68,
        "pyde_ngen": 42000,
        "pyde_save": 4200,
        "pyde_f": 0.5,
        "pyde_c": 0.5,
        "emcee_nwalkers": 68,
        "emcee_nruns": 2000,
        "emcee_thin_by": 100,
        "arg_in": {"progtype": 2, "NB": 3, "n_planets": 3},
    }

    # Read bodies.lst
    bodies_file = os.path.join(system_dir, "bodies.lst")
    if os.path.exists(bodies_file):
        with open(bodies_file) as f:
            for line in f:
                line = line.split("#")[0].strip()
                if not line:
                    continue
                parts = line.split()
                fname = parts[0]
                if fname == "star.dat":
                    system["n_planets"] = 0
                elif fname.endswith(".dat"):
                    system["n_planets"] += 1

    # Read star.dat
    star_file = os.path.join(system_dir, "star.dat")
    if os.path.exists(star_file):
        with open(star_file) as f:
            lines = [l.split("#")[0].strip() for l in f if l.strip()]
            if len(lines) >= 2:
                parts = lines[0].split()
                if len(parts) >= 2:
                    system["star"]["mass"] = float(parts[0])
                    system["star"]["mass_err"] = float(parts[1])
                parts = lines[1].split()
                if len(parts) >= 2:
                    system["star"]["radius"] = float(parts[0])
                    system["star"]["radius_err"] = float(parts[1])

    # Read planet .dat files
    for i in range(system["n_planets"]):
        label = get_planet_label(i)
        pfile = os.path.join(system_dir, f"{label}.dat")
        p = {
            "name": label,
            "transiting": True,
            "mass": 0.01, "mass_min": 0.0, "mass_max": 20.0,
            "radius": 0.1, "radius_min": 0.0, "radius_max": 20.0,
            "period": 10.0, "period_min": 1.0, "period_max": 1000.0,
            "sma": 999.0, "sma_min": 0.0, "sma_max": 0.0,
            "ecc": 0.0, "ecc_min": 0.0, "ecc_max": 1.0,
            "omega": 90.0, "omega_min": 0.0, "omega_max": 360.0,
            "M0": 0.0, "M0_min": 0.0, "M0_max": 360.0,
            "tp": 9e8, "tp_min": 0.0, "tp_max": 0.0,
            "inc": 90.0, "inc_min": 0.0, "inc_max": 180.0,
            "lnode": 180.0, "lnode_min": 0.0, "lnode_max": 360.0,
            "transit_obs": [],
        }
        if os.path.exists(pfile):
            with open(pfile) as f:
                dat_lines = [l.split("#")[0].strip() for l in f if l.strip() and not l.startswith("# epoch")]
                for j, dl in enumerate(dat_lines[:10]):
                    parts = dl.split()
                    if len(parts) >= 3:
                        vals = [float(x) for x in parts[:3]]
                        key_map = ["mass", "radius", "period", "sma", "ecc", "omega", "M0", "tp", "inc", "lnode"]
                        if key_map[j] == "mass":
                            p["mass"] = vals[0]; p["mass_min"] = vals[1]; p["mass_max"] = vals[2]
                        elif key_map[j] == "radius":
                            p["radius"] = vals[0]; p["radius_min"] = vals[1]; p["radius_max"] = vals[2]
                        elif key_map[j] == "period":
                            p["period"] = vals[0]; p["period_min"] = vals[1]; p["period_max"] = vals[2]
                        elif key_map[j] == "sma":
                            p["sma"] = vals[0]; p["sma_min"] = vals[1]; p["sma_max"] = vals[2]
                        elif key_map[j] == "ecc":
                            p["ecc"] = vals[0]; p["ecc_min"] = vals[1]; p["ecc_max"] = vals[2]
                        elif key_map[j] == "omega":
                            p["omega"] = vals[0]; p["omega_min"] = vals[1]; p["omega_max"] = vals[2]
                        elif key_map[j] == "M0":
                            p["M0"] = vals[0]; p["M0_min"] = vals[1]; p["M0_max"] = vals[2]
                        elif key_map[j] == "tp":
                            p["tp"] = vals[0]; p["tp_min"] = vals[1]; p["tp_max"] = vals[2]
                        elif key_map[j] == "inc":
                            p["inc"] = vals[0]; p["inc_min"] = vals[1]; p["inc_max"] = vals[2]
                        elif key_map[j] == "lnode":
                            p["lnode"] = vals[0]; p["lnode_min"] = vals[1]; p["lnode_max"] = vals[2]
        system["planets"].append(p)

    # Read arg.in for progtype and NB
    arg_file = os.path.join(system_dir, "arg.in")
    if os.path.exists(arg_file):
        with open(arg_file) as f:
            for line in f:
                if "progtype" in line and "=" in line:
                    system["arg_in"]["progtype"] = int(line.split("=")[1].strip())
                if "NB" in line and "=" in line:
                    try:
                        system["arg_in"]["NB"] = int(line.split("=")[1].strip())
                    except ValueError:
                        pass

    return system
