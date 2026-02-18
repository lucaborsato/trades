# pytrades module, it imports the fortran-python library and creates a trades object
import numpy as np
import os
import h5py
import sys

from . import constants as cst
from . import ancillary as anc
from .pytrades_lib import f90trades

from scipy.interpolate import interp1d

import pytransit
import numba

# os.environ["OMP_NUM_THREADS"] = "1"

# numba.set_num_threads(1)
# numba.config.THREADING_LAYER = "tbb"

# numba.config.DISABLE_JIT = 1


set_unit_base = anc.set_unit_base

import matplotlib.pyplot as plt

# =============================================================================
anc.set_rcParams()

# =============================================================================
get_data_info = f90trades.get_data_info
set_fitness = f90trades.set_fitness
convert_trades_par_to_kepelem = f90trades.convert_trades_par_to_kepelem
compute_ln_priors = f90trades.compute_ln_priors
check_boundaries = f90trades.check_boundaries
set_one_fit_par_boundaries = f90trades.set_one_fit_par_boundaries
reset_all_fit_boundaries = f90trades.reset_all_fit_boundaries
orbits_to_elements = f90trades.orbits_to_elements
set_close_encounter_check = f90trades.set_close_encounter_check
set_hill_check = f90trades.set_hill_check
set_amd_hill_check = f90trades.set_amd_hill_check
set_rv_res_gls = f90trades.set_rv_res_gls
period_to_sma = f90trades.f90_period_to_sma
period_to_sma_vec = f90trades.f90_period_to_sma_vec
sma_to_period = f90trades.f90_sma_to_period
astrocentric_to_barycentric_orbits = f90trades.astrocentric_to_barycentric_orbits
convert_trades_par_to_kepelem = f90trades.convert_trades_par_to_kepelem
print_obsdata = f90trades.print_obsdata
angular_momentum_deficit = f90trades.angular_momentum_deficit
angular_momentum_deficit_fit_parameters = (
    f90trades.angular_momentum_deficit_fit_parameters
)
check_rv_res_periodogram = f90trades.check_rv_res_periodogram
# angular_momentum_deficit_posterior = f90trades.angular_momentum_deficit_posterior
trueanomaly_to_eccanomaly = f90trades.trueanomaly_to_eccanomaly
meananomaly_to_eccanomaly = f90trades.meananomaly_to_eccanomaly
meananomaly_to_trueanomaly = f90trades.meananomaly_to_trueanomaly
# =============================================================================


def args_init(
    n_body,
    duration_check,
    t_epoch=None,
    t_start=None,
    t_int=None,
    encounter_check=True,
    do_hill_check=False,
    amd_hill_check=False,
    rv_res_gls=False,
):
    """
    Initialize the arguments for the function.

    Parameters:
        n_body (int): Number of bodies.
        duration_check (int): Duration check value.
        t_epoch (int, optional): Optional epoch time value.
        t_start (int, optional): Optional start time value.
        t_int (int, optional): Optional interval time value.
        encounter_check (bool, optional): Optional flag for close encounter check. Default True.
        do_hill_check (bool, optional): Optional flag for Hill check. Default False.
        amd_hill_check (bool, optional): Optional flag for AMD Hill check. Default False.
        rv_res_gls (bool, optional): Optional flag for RV GLS check of inserted periods during fit.
    Returns:
        None
    """

    f90trades.args_init(
        n_body,
        duration_check,
    )

    if t_epoch is not None:
        f90trades.set_epoch_time(t_epoch)

    if t_start is not None:
        f90trades.set_starting_time(t_start)

    if t_int is not None:
        f90trades.set_integration_time(t_int)

    if encounter_check is not None:
        set_close_encounter_check(encounter_check)

    if do_hill_check is not None:
        set_hill_check(do_hill_check)
    if amd_hill_check is not None:
        set_amd_hill_check(amd_hill_check)
    if rv_res_gls is not None:
        set_rv_res_gls(rv_res_gls)

    return


# =============================================================================
# get_priors = f90trades.get_priors
def get_priors():
    """
    Retrieves the priors from f90trades, decodes the names to utf-8, and returns the names and values.

    Parameters:
    None

    Returns:
    names (list): A list of prior names decoded to utf-8.
    values (numpy array): An array of prior values.
    """

    n_priors = (f90trades.n_priors).item()
    names, values = f90trades.get_priors(n_priors)
    names = [n.decode("utf-8") for n in names]

    return names, values


# =============================================================================
def deallocate_rv_dataset():
    """
    deallocate_rv_dataset function deallocates the RV dataset.
    """

    f90trades.deallocate_rv_dataset()

    return


# =============================================================================
def set_rv_dataset(t_rv, rv_obs, erv_obs, rv_setid=None, n_rvset=1):
    """
    Set the RV dataset for a given set of RV observations.

    Parameters:
    t_rv (array): Array of time values for the RV observations.
    rv_obs (array): Array of RV observations.
    erv_obs (array): Array of errors on the RV observations.
    rv_setid (array, optional): Array defining the RV set IDs. If not provided, default set IDs are generated.
    n_rvset (int, optional): Number of RV sets.

    Returns:
    None
    """

    if rv_setid is not None:
        n = len(np.unique(rv_setid))
        if n != n_rvset:
            n_rvset = n
    else:
        rv_setid = np.ones((len(t_rv)))
    f90trades.set_rv_dataset(t_rv, rv_obs, erv_obs, rv_setid.astype(int), n_rvset)

    return


# =============================================================================
def set_t0_dataset(body_id, epo, t0, et0, t14=None, et14=None, sources_id=None):
    """
    A function to set the t0 dataset with the given parameters.

    Parameters:
    body_id (int): The ID of the body for which the t0 dataset needs to be set.
    epo (array): Array of epoch values for the t0 observations.
    t0 (array): Array of t0 values for the observations.
    et0 (array): Array of errors on the t0 observations.
    t14 (array): Array of t14 (durations) values for the observations.
    et14 (array): Array of errors on the t14 (durations) observations.
    sources_id (array, optional): Array defining the sources IDs. If not provided, default sources IDs are generated.

    Returns:
    None
    """
      

    # If sources_id is not provided, generate default sources IDs.
    if sources_id is None:
        sources_id = np.ones((len(epo)), dtype=int)

    # Set the t0 dataset using the F90 wrapper.
    f90trades.set_t0_dataset(body_id, epo, t0, et0, sources_id)
    # Set durations if provided
    if (t14 is not None) and (et14 is not None):
        f90trades.set_t14_dataset(body_id, t14, et14)

    return


# =============================================================================
def deallocate_t0_dataset(body_id):
    """
    Deallocates a T0 dataset for a given body ID.

    Parameters:
    body_id (int): The ID of the body for which the T0 dataset needs to be deallocated.

    Returns:
    None
    """

    f90trades.deallocate_t0_dataset(body_id)

    return


# =============================================================================
def kelements_to_observed_rv_and_t0s(
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
    transit_flag,
):
    """
    Computes observed radial velocities (RVs) and T0s based on given Keplerian elements.
    
    Parameters:
        t_epoch (float): epoch time
        t_start (float): start time
        t_int (float): integration time step
        M_msun (float): mass of the sun
        R_rsun (float): radius of the sun
        P_day (float): orbital period in days
        ecc (float): eccentricity
        argp_deg (float): argument of periapsis in degrees
        mA_deg (float): mean anomaly in degrees
        inc_deg (float): inclination in degrees
        lN_deg (float): longitude of the ascending node in degrees
        transit_flag (list): flag for transiting bodies: [False, True, True, ...]
    
    Returns:
        tuple: Contains rv_sim (RVs), body_T0_sim (T0s), epo_sim (epochs), t0_sim, t14_sim, lambda_rm_sim, kel_sim, stable
    """
   
    n_kep = 8  # number of keplerian elements in output for each T0s
    n_rv = f90trades.nrv
    n_T0s = f90trades.ntts
    # print("Computing RVs and T0s ...", flush=True)
    # print("n_rv    = {}".format(n_rv), flush=True)
    # print("n_T0s   = {}".format(n_T0s), flush=True)
    # print("t_epoch = {}".format(t_epoch), flush=True)
    # print("t_start = {}".format(t_start), flush=True)
    # print("t_int   = {}".format(t_int), flush=True)
    # print("", flush=True)

    (rv_sim, body_T0_sim, epo_sim, t0_sim, t14_sim, lambda_rm_sim, kel_sim, stable) = (
        f90trades.kelements_to_rv_and_t0s(
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
            transit_flag,
            n_rv,
            n_T0s,
            n_kep,
        )
    )

    return rv_sim, body_T0_sim, epo_sim, t0_sim, t14_sim, lambda_rm_sim, kel_sim, stable


# =============================================================================
def kelements_to_observed_rv(
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
    # n_rv,
):
    """
    Computes observed radial velocities (RVs) based on given Keplerian elements.
    
    Parameters:
        t_epoch (float): epoch time
        t_start (float): start time
        t_int (float): integration time step
        M_msun (float): mass of the sun
        R_rsun (float): radius of the sun
        P_day (float): orbital period in days
        ecc (float): eccentricity
        argp_deg (float): argument of periapsis in degrees
        mA_deg (float): mean anomaly in degrees
        inc_deg (float): inclination in degrees
        lN_deg (float): longitude of the ascending node in degrees
    
    Returns:
        tuple: Contains rv_sim (RVs) and stable flag
    """

    n_rv = f90trades.nrv
    rv_sim, stable = f90trades.kelements_to_rv(
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
        n_rv,
    )

    return rv_sim, stable


# =============================================================================
def kelements_to_observed_t0s(
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
    transit_flag,
):
    """
    Computes orbits from a given set of Keplerian elements and generates RVs and T0s.

    Parameters:
        t_epoch (float): epoch time
        t_start (float): start time
        t_int (float): integration time step
        M_msun (float): mass of the sun
        R_rsun (float): radius of the sun
        P_day (float): period in days
        ecc (float): eccentricity
        argp_deg (float): argument of periapsis in degrees
        mA_deg (float): mean anomaly in degrees
        inc_deg (float): inclination in degrees
        lN_deg (float): longitude of the ascending node in degrees
        transit_flag (int): flag for transit

    Returns:
        tuple: Contains rv_sim (RVs), body_T0_sim (T0s), epo_sim (epochs), t0_sim, t14_sim, lambda_rm_sim, kel_sim, stable
    """

    n_kep = 8  # number of keplerian elements in output for each T0s
    n_T0s = f90trades.ntts
    body_T0_sim, epo_sim, t0_sim, t14_sim, lambda_rm_sim, kel_sim, stable = (
        f90trades.kelements_to_t0s(
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
            transit_flag,
            n_T0s,
            n_kep,
        )
    )

    return body_T0_sim, epo_sim, t0_sim, t14_sim, lambda_rm_sim, kel_sim, stable


# =============================================================================
def kelements_to_orbits(
    steps, M_msun, R_rsun, P_day, ecc, argp_deg, mA_deg, inc_deg, lN_deg
):
    """
    Calculate orbits from Keplerian elements.

    Parameters:
        steps (list): time steps.
        M_msun (list): List of masses in solar masses.
        R_rsun (list): List of radii in solar radii.
        P_day (list): List of orbital periods in days.
        ecc (list): List of eccentricities.
        argp_deg (list): List of argument of pericenter in degrees.
        mA_deg (list): List of mean anomalies in degrees.
        inc_deg (list): List of inclinations in degrees.
        lN_deg (list): List of longitude of ascending node in degrees.

    Returns:
        orbits (list): List of calculated orbits.
        check (int): Check value (stable or not).
    """

    # Calculate the number of bodies
    nb_dim = len(M_msun) * 6

    # Call the fortran routine to calculate the orbits
    orbits, check = f90trades.kelements_to_orbits(
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

    # Return the calculated orbits and the check value
    return orbits, check


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
    """
    Generate orbits for a given set of keplerian elements over a specified time interval.

    Parameters:
    - t_epoch: float, epoch time
    - t_start: float, start time of the interval
    - t_int: float, length of the time interval
    - M_msun: array-like, masses of the bodies in solar masses
    - R_rsun: array-like, radii of the bodies in solar radii
    - P_day: array-like, orbital periods in days
    - ecc: array-like, eccentricities of the orbits
    - argp_deg: array-like, argument of periastron in degrees
    - mA_deg: array-like, mean anomaly at epoch in degrees
    - inc_deg: array-like, inclination in degrees
    - lN_deg: array-like, longitude of the ascending node in degrees
    - specific_times: array-like or None, specific times to include in the orbit computation
    - step_size: float or None, step size for the orbit computation
    - n_steps_smaller_orbits: float, number of steps for the smallest orbit

    Returns:
    - time_steps: array, time steps for the orbits
    - orbits: array, calculated orbits
    - check: bool, check flag for the orbits
    """

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
        orbits_backward, check_backward = kelements_to_orbits(
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
        check_backward = True

    if np.abs(dt_end_epoch) > 0.0:
        # forward integration
        steps_forward = np.arange(t_epoch, t_end, step_size * np.sign(dt_end_epoch))
        if steps_forward[-1] < t_end:
            steps_forward = np.concatenate((steps_forward, [t_end]))
        # add specific times to compute the orbits, such as observed RV times
        t_in = t_add[t_add >= t_epoch]
        steps_forward = np.sort(np.concatenate((steps_forward, t_in)))

        # n_steps_forward = len(steps_forward)
        orbits_forward, check_forward = kelements_to_orbits(
            steps_forward, M_msun, R_rsun, P_day, ecc, argp_deg, mA_deg, inc_deg, lN_deg
        )
    else:
        steps_forward = np.array([])
        orbits_forward = np.empty((0, nb_dim))
        check_forward = True

    time_steps = np.concatenate((steps_backward, steps_forward))
    orbits = np.concatenate((orbits_backward, orbits_forward))
    check = check_backward & check_forward

    return time_steps, orbits, check


# =============================================================================


def orbits_to_rvs(M_msun, orbits):
    """
    Calculate radial velocity from stellar mass and orbital state vectors.
    Parameters:
    - M_msun: float, shape (n_body), mass of the star in solar masses
    - orbits: array-like, shape (n_steps, n_body * 6),orbital state vectors
    Returns:
    - rvs: array-like, radial velocity
    """

    rvs = f90trades.orbits_to_rvs(M_msun, orbits)

    return rvs


# =============================================================================


def orbits_to_transits(
    n_all_transits, time_steps, M_msun, R_rsun, orbits, transiting_body=1
):
    """
    Convert orbital state vectors to transit times and durations.

    Parameters:
        n_all_transits (int): Total number of transits
        time_steps (ndarray): Array of time steps
        M_msun (float): Masses in solar masses
        R_rsun (float): Radii in solar radii
        orbits (ndarray): Array of orbital state vectors
        transiting_body (str): Name of the transiting body, 1 means all bodies

    Returns:
        tuple: A tuple containing the transit times, durations, spin-orbit misalignment, Keplerian elements, and body flags
    """

    # Pass the input arguments to the fortran function
    transits, durations, lambda_rm, kep_elem, body_flag = f90trades.orbits_to_transits(
        n_all_transits, time_steps, M_msun, R_rsun, orbits, transiting_body
    )

    # Remove elements with body_flag < 2
    # This is because body_flag is set to 1 for bodies that do not transit
    # and set to 2 for bodies that do transit
    sel = body_flag > 1
    transits, durations, lambda_rm, kep_elem, body_flag = (
        transits[sel],
        durations[sel],
        lambda_rm[sel],
        kep_elem[sel, :],
        body_flag[sel],
    )

    # Return the results
    return transits, durations, lambda_rm, kep_elem, body_flag


# =============================================================================


def linear_fit(x, y, ey=None):
    """
    Perform a linear fit on the given data points.

    Parameters:
    x (array-like): The x-values of the data points.
    y (array-like): The y-values of the data points.
    ey (array-like, optional): The errors associated with the y-values. Defaults to None.

    Returns:
    tuple: Tuple containing two tuples. The first tuple contains the intercept and its error (q, err_q).
    The second tuple contains the slope and its error (m, err_m).
    float: The chi-squared value of the fit. Be aware that if ey is None, the chi-squared value is not normalized and has no meaning.
    """

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
    t_epoch, t_start, t_int, mass, radius, period, ecc, w, ma, inc, long, t_rv_obs
):
    """
    From orbital parameters compute orbits of the planets and return
    the transits (times, durations, kep. elements. and body flags) and radial velocity at specific times.
    Similar to kelements_to_orbits_full + orbits_to_transits + orbits_to_rvs

    Parameters:
    - t_epoch: float, epoch time
    - t_start: float, start time
    - t_int: float, interval time
    - mass: array-like, mass values
    - radius: array-like, radius values
    - period: array-like, period values
    - ecc: array-like, eccentricity values
    - w: array-like, w values
    - ma: array-like, ma values
    - inc: array-like, inclination values
    - long: array-like, longitude values
    - t_rv_obs: array-like, observed RV times

    Returns:
    - time_steps: array-like, time steps
    - orbits: array-like, orbital parameters
    - transits: array-like, transit information
    - lambda_rm: array-like, spin-orbit misalignment
    - durations: array-like, duration values
    - kep_elem: array-like, Kepler elements
    - body_flag: array-like, body flags
    - rv_sim: dictionary, simulated RV values
    - stable: int, stability flag: 1=stable/0=unstable
    """

    # Compute the orbits from the orbital parameters
    time_steps, orbits, stable = kelements_to_orbits_full(
        t_epoch,
        t_start,
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

    # Select the time steps that correspond to the observed RV times
    sel_t_rv = np.isin(time_steps, t_rv_obs)
    # Compute the radial velocity at the selected times
    rv = orbits_to_rvs(mass, orbits[sel_t_rv, :])
    # Create a dictionary with the simulated RV values
    rv_sim = {"time": time_steps[sel_t_rv], "rv": rv}

    # Set the transiting body to 1 (all planets)
    transiting_body = 1
    n_body = len(mass)
    # Determine the number of transits ... it has to be done in advance
    # n_transits = (t_int / period[1:]).astype(int)
    # n_transits = len(time_steps)
    # n_all_transits = np.sum(n_transits) + (
    #     n_body - 1
    # )  # star has no transits by definition
    n_all_transits = len(time_steps)*(n_body-1)
    # Compute the transits, durations, lambda_rm, Kepler elements, and body flags
    transits, durations, lambda_rm, kep_elem, body_flag = orbits_to_transits(
        n_all_transits, time_steps, mass, radius, orbits, transiting_body
    )
    # kep_elem == period 0, sma 1, ecc 2, inc 3, meana 4, argp 5, truea 6, longn 7
    return (
        time_steps,
        orbits,
        transits,
        durations,
        lambda_rm,
        kep_elem,
        body_flag,
        rv_sim,
        stable,
    )


def set_transit_parameters(radius, transits, body_flag, kep_elem):
    """
    Set the transit parameters based on the given inputs.

    Parameters:
    - radius: array-like, the radius values
    - transits: array-like, the transit values
    - body_flag: int, flag indicating the body
    - kep_elem: array-like, the Keplerian elements

    Returns:
    - rp_rs: array-like, ratio of radius to the reference radius
    - per: array-like, period values
    - aRs: array-like, semi-major axis values
    - inc: array-like, inclination values
    - ecc: array-like, eccentricity values
    - w: array-like, argument of periastron values
    """

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
    rp_rs,
    ld_quad,
    per,
    aRs,
    inc,
    ecc,
    w,
    #  body_flag,
    time_key="time",
):
    """
    Generate simulated flux for given parameters and photometry data.

    Parameters:
    - tm: The PyTransit model to use for simulation
    - photometry: Dictionary containing photometry data
    - transits: Array of transit times
    - durations: Array of transit durations
    - rp_rs: Array of planet-to-star radius ratios
    - ld_quad: Limb darkening coefficients
    - per: Array of planet orbital periods
    - aRs: Array of semi-major axis to stellar radius ratios
    - inc: Array of orbital inclinations
    - ecc: Array of eccentricities
    - w: Array of arguments of periastron
    # - body_flag: Body number of the planet/bodies, first planet is 2 (start is 1)
    - time_key: Key to access time values in the photometry dictionary

    Returns:
    - sim_photometry: Dictionary containing simulated flux for each dataset
    """
    if "full" in time_key:
        n_over_key = "n_oversample_full"
        t_exp_key = "t_exp_d_full"
    else:
        n_over_key = "n_oversample"
        t_exp_key = "t_exp_d"
    sim_photometry = {}
    for k, vis in photometry.items():
        t = vis[time_key]

        # select transiting planets in the time range
        tra_in_t = np.logical_and(transits >= t.min(), transits <= t.max())
        n_tra = np.sum(tra_in_t)
        # select partial transits in the time range
        tra_dur_in_t = np.logical_and(
            transits - 0.5 * durations * cst.min2day >= t.min(),
            transits + 0.5 * durations * cst.min2day <= t.max(),
        )
        n_dur = np.sum(tra_dur_in_t)
        # number of events based on the max between n_tra and n_dur
        # n = max(n_tra, n_dur)
        if n_tra >= n_dur:
            n = n_tra
            sel_tra = tra_in_t
        else:
            n = n_dur
            sel_tra = tra_dur_in_t
        # sel_body = body_flag[sel_tra]

        if n > 0:
            # u_nb = np.unique(sel_body)
            # flux = []
            # for bd in u_nb:
            #     sel = np.logical_and(sel_tra, sel_body == bd)
            #     nsel = np.sum(sel)
            #     tra = np.atleast_1d(transits[sel])
            #     Px  = np.atleast_1d(per[sel])[0]
            #     epo = np.rint( (tra - tra[0]) / Px).astype(int)
            #     lid = np.rint( (t   - tra[0]) / Px).astype(int)
            #     tm.set_data(
            #         t,
            #         lcids=lid,
            #         epids=epo,
            #         nsamples=vis[n_over_key],
            #         exptimes=vis[t_exp_key]
            #     )
            #     ff = tm.evaluate(
            #         k=rp_rs[sel],
            #         ldc=np.array([ld_quad] * nsel),
            #         t0=tra,
            #         p=per[sel],
            #         a=aRs[sel],
            #         i=inc[sel],
            #         e=ecc[sel],
            #         w=w[sel],
            #     )
            #     flux.append(ff)

            tra_sel = np.atleast_1d(transits[sel_tra])
            rp_rs_sel = np.atleast_1d(rp_rs[sel_tra])
            per_sel = np.atleast_1d(per[sel_tra])
            aRs_sel = np.atleast_1d(aRs[sel_tra])
            inc_sel = np.atleast_1d(inc[sel_tra])
            ecc_sel = np.atleast_1d(ecc[sel_tra])
            w_sel = np.atleast_1d(w[sel_tra])

            flux = []
            for itra, tra in enumerate(tra_sel):
                # sel_t = np.logical_and(
                sel_t = (t >= tra - 0.5 * per[itra],)
                tm.set_data(t[sel_t], nsamples=vis[n_over_key], exptimes=vis[t_exp_key])
                ff = tm.evaluate(
                    k=rp_rs_sel[itra],
                    ldc=ld_quad,
                    t0=tra,
                    p=per_sel[itra],
                    a=aRs_sel[itra],
                    i=inc_sel[itra],
                    e=ecc_sel[itra],
                    w=w_sel[itra],
                )
                flux.append(ff)
            f2d = np.atleast_2d(flux)
            flux_ = np.sum(f2d - 1.0, axis=0) + 1.0
        else:
            flux_ = np.ones((len(t)))
        sim_photometry[k] = flux_

    return sim_photometry


# TODO:
# GET PHOTOMETRY TRENDS
# GET PHOTOMETRY ANCILLARY
# def get_photometry_trend(obs_photometry, trend_coeff, time_key="time"):

#     out_trend = {}
#     for name_phot, phot in obs_photometry.items():
#         out_trend[name_phot] = {}
#         for k, vis in phot.items():
#             t = vis[time_key]
#             nt = len(t)
#             ctrend = trend_coeff[name_phot][k]
#             tscale = t - t.min()
#             ftrend = np.ones((nt))

#     return

def add_photometry_trend(
    obs_photometry,
    sim_photometry,
    photometry_coeff_trend,
    ancillary_coeff,
    time_key="time",
):
    """
    Generate photometry trend for observed and simulated photometry data.

    Parameters:
    - obs_photometry (dict): A dictionary containing observed photometry data.
    - sim_photometry (dict): A dictionary containing simulated photometry data.
    - photometry_coeff_trend (dict): A dictionary containing photometry coefficient trends.
    - ancillary_coeff (dict): A dictionary containing ancillary coefficients.
    - time_key (str): The key for the time values in the photometry data. Defaults to "time".

    Returns:
    - out_photometry (dict): A dictionary containing the photometry data with trends applied.
    - out_trend (dict): A dictionary containing the trend values for each photometry data.
    """

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
            ctrend = photometry_coeff_trend[name_phot][k]
            tscale = t - vis["time_min"]
            
            if ctrend is not None:
                # ftrend = np.zeros((nt))
                ftrend = 0.0
                for o, c in enumerate(ctrend):
                    ftrend += c * (tscale**o)
            else:
                ftrend = np.ones((nt))
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
get_photometry_full_trend = add_photometry_trend


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
    """
    A function to plot photometry data along with simulated, model, and trend photometry.

    Parameters:
    - photometry: dictionary containing photometry data
    - sim_photometry: dictionary containing simulated photometry data
    - mod_photometry: dictionary containing model photometry data (optional)
    - trend_photometry: dictionary containing trend photometry data (optional)
    - figsize: tuple specifying the figure size (default is (3, 3))
    - show_plot: boolean indicating whether to display the plot (default is True)
    - output_folder: string specifying the output folder for saving plots
    - return_rms: boolean indicating whether to return the root mean square values

    Returns:
    - rms_photometry: dictionary containing root mean square values if return_rms is True, otherwise None
    """

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
    """
    Plot the observed and simulated O-C (Observed minus Calculated) values for transits of multiple planets.

    Parameters:
    - transits_obs: Dictionary containing observed transit data for each planet.
    - Tref: Dictionary containing reference transit times for each planet.
    - Pref: Dictionary containing reference periods for each planet.
    - transits_sim: Dictionary containing simulated transit times for each planet.
    - pl_names: List of planet names.
    - figsize: Tuple representing the figure size (default is (8, 4)).
    - show_plot: Boolean indicating whether to display the plot (default is True).
    - output_folder: Path to the output folder for saving the plot (default is None).

    Returns:
    - None
    """

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
    """
    Creates a portion dictionary containing time, flux, flux error, ancillary data, and additional derived information.

    Parameters:
        time (array-like): Array of time values.
        flux (array-like): Array of flux values.
        flux_err (array-like): Array of flux error values.
        n_oversample (int, optional): Number of oversampling points. Default is 1.
        t_exp_d (float, optional): Exposure time in days. Default is 60.0 sec in days.
        ancillary (dict, optional): Dictionary with ancillary data. Default is None.

    Returns:
        dict: A dictionary containing time, flux, flux error, ancillary data, and additional information including time statistics and interpolated ancillary data.
    """

    ndata = len(time)
    portion = {}
    portion["time"] = time
    portion["flux"] = flux
    portion["flux_err"] = flux_err
    portion["ancillary"] = ancillary
    if ancillary is not None:
        portion["ancillary_interp"] = {
            k: interp1d(time, a, bounds_error=False, fill_value=(a[0], a[-1]))
            for k, a in ancillary.items()
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
    t_exp_d_full = min(t_exp_d, 60.0 * cst.sec2day)
    portion["t_exp_d_full"] = t_exp_d_full
    portion["time_full"] = np.arange(
        time_min, time_max + (0.5 * t_exp_d_full), t_exp_d_full
    )
    portion["ndata_full"] = len(portion["time_full"])
    portion["n_oversample_full"] = 1

    return portion


def get_photometry_residuals(obs_photometry, sim_photometry):
    """
    Calculate photometry residuals between observed and simulated photometry data.

    Parameters:
        obs_photometry (dict): A dictionary containing observed photometry data.
        sim_photometry (dict): A dictionary containing simulated photometry data.

    Returns:
        res (array): Array of photometry residuals.
        err2 (array): Array of squared photometry errors.
    """
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
    """
    Calculate photometry residuals between observed and simulated transit times data.

    Parameters:
        obs_photometry (dict): A dictionary containing observed transit times data.
        sim_photometry (dict): A dictionary containing simulated transit times data.

    Returns:
        res (array): Array of photometry residuals.
        err2 (array): Array of squared photometry errors.
    """
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
    """
    Calculate the residuals and errors for radial velocities data.

    Parameters:
        t_epoch (array): Array of epoch times.
        t_rv_obs (array): Array of observed radial velocities times.
        rv_obs (array): Array of observed radial velocities.
        rv_obs_err (array): Array of errors in observed radial velocities.
        rv_set_id (array): Array of IDs for radial velocity sets.
        gammas (array): Array of gamma values.
        jitters (array): Array of jitter values.
        sim_rv (dict): Dictionary containing simulated radial velocities.
        ctrends (list, optional): List of trend coefficients. Defaults to [0.0].

    Returns:
        res (array): Array of residuals for radial velocities.
        err2 (array): Array of squared errors for radial velocities.
    """
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
    def __init__(
        self,
        n_body,
        t_epoch,
        t_start,
        t_int,
        encounter_check=True,
        duration_check=1,
        do_hill_check=False,
        amd_hill_check=False,
        rv_res_gls=False,
    ):
        """
        Initializes the PhotoTRADES object with the specified parameters.

        Parameters:
            n_body (int): Number of bodies.
            t_epoch (int): Epoch time.
            t_start (int): Start time.
            t_int (int): Interval time.
            duration_check (int, optional): Duration check value (default is 1).
            do_hill_check (bool, optional): Flag for Hill check.
            amd_hill_check (bool, optional): Flag for AMD Hill check.
            rv_res_gls (bool, optional): Flag for RV GLS.

        Returns:
            None
        """
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
            encounter_check=encounter_check,
            do_hill_check=do_hill_check,
            amd_hill_check=amd_hill_check,
            rv_res_gls=rv_res_gls,
        )

        self.Tref = {}
        self.Pref = {}

        self.n_data = 0

        self.transits = {}
        self.n_transits = {}
        self.n_all_transits = 0

        self.photometry = {}
        self.photometry_flux_trend_coeff = {}
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
        """
        Set the default physical parameters for the object.

        Parameters:
            mass (float): The mass of the object.
            radius (float): The radius of the object.
            period (float): The period of the object.
            ecc (float): The eccentricity of the object.
            argp (float): The argument of the pericenter of the object.
            meana (float): The mean anomaly of the object.
            inc (float): The inclination of the object.
            longn (float): The longitude of the ascending node of the object.

        Returns:
            None
        """
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
        """
        Adds linear ephemeris data for a specific body.
        
        Parameters:
            self: The object instance.
            body_id (int): The ID of the body for which the ephemeris data is being added.
            Tref (float or list): The reference time for the ephemeris data.
            Pref (float or list): The reference position for the ephemeris data.
        
        Returns:
            None
        """
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
        """
        Adds transit data for a specific body.
        
        Parameters:
            self: The object instance.
            body_id (int): The ID of the body for which the transit data is being added.
            T0s (array): Array of transit times.
            err_T0s (array): Array of errors on the transit times.
            epo (optional): The epoch value.
            sources (array, optional): Array of sources for the transit data.

        Returns:
            None
        """
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
        """
        Updates linear ephemeris data based on transit times and reference values for each body.
        This function iterates over transits, calculates Tx and Px values, sorts the transit times, and computes the ephemeris.
        
        Parameters:
            self: The object instance.
            
        Returns:
            None
        """
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
        """
        Generates simulated transits based on input transits and body flags.
        
        Parameters:
            self: The object instance.
            transits: The input transits data.
            body_flag: The flag indicating the body.
            
        Returns:
            transits_out: A dictionary containing the simulated transits.
        """
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

    def add_photometry(self, name, data, flux_time_trend_coeff, anc_coeff=None):
        """
        Adds photometry data to the object instance.

        Parameters:
            self: The object instance.
            name (str): The name of the photometry data.
            data (dict): The photometry data to be added.
            flux_time_trend_coeff (dict): The flux time trend coefficients to be added.
            anc_coeff (dict, optional): The ancillary coefficients data. Defaults to None.

        Returns:
            None
        """
        self.photometry[name] = data
        self.photometry_flux_trend_coeff[name] = flux_time_trend_coeff
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
        """
        Adds radial velocity data to the object instance.

        Parameters:
            self: The object instance.
            name (str): The name of the radial velocity data.
            data (dict): The radial velocity data to be added.
            id_num (int): The ID number associated with the radial velocity data.

        Returns:
            None
        """
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
        """
        Sorts the radial velocity observations based on the time of observation.
        No parameters.
        Returns None.
        """
        sort_t_rv = np.argsort(self.t_rv_obs)
        self.t_rv_obs = self.t_rv_obs[sort_t_rv]
        self.rv_obs = self.rv_obs[sort_t_rv]
        self.rv_obs_err = self.rv_obs_err[sort_t_rv]
        self.rv_set_id = self.rv_set_id[sort_t_rv]
        self.rv_set_name = np.array(self.rv_set_name)[sort_t_rv]
        self.idx_rv_obs = self.idx_rv_obs[sort_t_rv]

        return

    def update_n_data(self):
        """
        Updates the total number of data points by summing the number of photometry points, all transits, and radial velocity points.
        No parameters.
        Returns None.
        """
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
        """
        Calculates the time steps, orbits, and a check value based on the provided mass, radius, period, eccentricity, argument of periapsis, mean anomaly, inclination, and longitude of the ascending node. 
        Optionally, a step size and number of steps for smaller orbits can be specified. Returns the calculated time steps, orbits, and stable flag.
        """

        time_steps, orbits, stable = kelements_to_orbits_full(
            self.t_epoch,
            self.t_start,
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

        return time_steps, orbits, stable

    def orbital_parameters_to_transits(
        self, mass, radius, period, ecc, w, ma, inc, long
    ):
        """
        Calculate transits based on the orbital parameters.

        Parameters:
            - mass (float): Mass of the celestial body.
            - radius (float): Radius of the celestial body.
            - period (float): Orbital period of the celestial body.
            - ecc (float): Eccentricity of the orbit.
            - w (float): Argument of periastron.
            - ma (float): Mean anomaly.
            - inc (float): Inclination of the orbit.
            - long (float): Longitude of the ascending node.

        Returns:
            - time_steps (list): List of time steps.
            - orbits (list): List of orbits.
            - transits (list): List of transits.
            - durations (list): List of durations.
            - lambda_rm (list): List of lambda_rm values.
            - kep_elem (list): List of keplerian elements.
            - body_flag (list): List of body flags.
            - rv_sim (list): List of simulated radial velocities.
            - stable (int): 1=stable/0=unstable.
        """

        (
            time_steps,
            orbits,
            transits,
            durations,
            lambda_rm,
            kep_elem,
            body_flag,
            rv_sim,
            stable,
        ) = orbital_parameters_to_transits(
            self.t_epoch,
            self.t_start,
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
        return (
            time_steps,
            orbits,
            transits,
            durations,
            lambda_rm,
            kep_elem,
            body_flag,
            rv_sim,
            stable,
        )

    def get_simulate_flux(
        self,
        radius,
        ld_quads,
        transits,
        durations,
        body_flag,
        kep_elem,
        time_key="time",
    ):
        """
        A function to simulate flux based on input parameters and return simulated photometry.
        
        Parameters:
            radius (float): The radius parameter.
            ld_quads (dict): Dictionary of LD (limb darkening) quadratics.
            transits (list): List of transits.
            durations (list): List of durations.
            body_flag (int): Body flag parameter.
            kep_elem (dict): Dictionary of Keplerian elements.
            time_key (str, optional): Key for time parameter. Defaults to "time".
        
        Returns:
            dict: A dictionary containing simulated photometry for each photo name.
        """
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
                # body_flag,
                time_key=time_key,
            )
            sim_photometry[phot_name] = sim_phot
        return sim_photometry

    def full_photodyn(self, mass, radius, period, ecc, w, ma, inc, long, ld_quads):
        """
        A function to perform full photodynamic calculations based on the given parameters.

        Parameters:
            mass (float): The mass parameter.
            radius (float): The radius parameter.
            period (float): The period parameter.
            ecc (float): The eccentricity parameter.
            w (float): The argument of pericenter parameter.
            ma (float): The mean anomaly parameter.
            inc (float): The inclination parameter.
            long (float): The longitude of the node parameter.
            ld_quads (dict): Dictionary of limb darkening quadratics.

        Returns:
            tuple: A tuple containing simulated photometry, radial velocity simulation, and simulated transits.
        """
        (
            time_steps,
            orbits,
            transits,
            durations,
            lambda_rm,
            kep_elem,
            body_flag,
            rv_sim,
            stable,
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
        """
        Plot the photometry data along with optional model and trend photometry.
        
        Parameters:
            sim_photometry (dict): A dictionary containing simulated photometry data.
            mod_photometry (dict, optional): A dictionary containing model photometry data. Defaults to None.
            trend_photometry (dict, optional): A dictionary containing trend photometry data. Defaults to None.
            figsize (tuple, optional): The figure size for the plot. Defaults to (3, 3).
            show_plot (bool, optional): A flag to show the plot. Defaults to True.
            output_folder (str, optional): The output folder path. Defaults to None.
            return_rms (bool, optional): A flag to return RMS value. Defaults to False.
        
        Returns:
            rms: The root mean square value calculated from the photometry data.
        """
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
        """
        Plot the O-C diagram for the transits simulation data.

        Parameters:
            transits_sim (dict): A dictionary containing the simulated transits data.
            pl_names (list): A list of planet names.
            figsize (tuple, optional): A tuple specifying the figure size. Default is (8, 4).
            show_plot (bool, optional): A boolean indicating whether to display the plot. Default is True.
            output_folder (str, optional): A string specifying the output folder path. Default is None.
        """
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
        """
        Plot the RV observations along with the simulated RV data and residuals.

        Parameters:
            rv_sim (dict): A dictionary containing the simulated RV data.
            gamma_rv (array): Array of gamma RV values.
            jitter_rv (dict): A dictionary containing jitter RV values.
            markers (dict): A dictionary of markers for each dataset.
            rv_trend_coeff (float or list, optional): Coefficients for RV trend calculation.
            print_rv (bool, optional): Flag to print RV data. Default is False.
            figsize (tuple, optional): A tuple specifying the figure size. Default is (4, 2).
            show_plot (bool, optional): Flag to display the plot. Default is True.
            output_folder (str, optional): Path to the output folder. Default is None.
            remove_dataset (list, optional): List of datasets to remove.

        Returns:
            None
        """
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
        """
        Calculate photometry residuals between observed and simulated photometry data.

        Parameters:
            obs_photometry (dict): A dictionary containing observed photometry data.
            sim_photometry (dict): A dictionary containing simulated photometry data.

        Returns:
            res (array): Array of photometry residuals.
            err2 (array): Array of squared photometry errors.
        """
        res, err2 = get_photometry_residuals(self.photometry, sim_photometry)

        return res, err2

    def get_transit_times_residuals(self, sim_transits):
        """
        Calculate photometry residuals between observed and simulated transit times data.

        Parameters:
            obs_photometry (dict): A dictionary containing observed transit times data.
            sim_photometry (dict): A dictionary containing simulated transit times data.

        Returns:
            res (array): Array of photometry residuals.
            err2 (array): Array of squared photometry errors.
        """
        res, err2 = get_transit_times_residuals(self.transits, sim_transits)

        return res, err2

    def get_radial_velocities_residuals(self, gammas, jitters, rv_sim, ctrends=[0.0]):
        """
        Calculate radial velocities residuals between observed and simulated radial velocity data.

        Parameters:
            gammas (array): Array of gamma parameters.
            jitters (array): Array of jitter parameters.
            rv_sim (array): Array of simulated radial velocities.
            ctrends (list, optional): List of radial velocity correction trends. Defaults to [0.0].

        Returns:
            res (array): Array of radial velocity residuals.
            err2 (array): Array of squared radial velocity errors.
        """
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
        """
        Initializes the TRADESfolder class with input parameters.

        Parameters:
            input_path (str): The input path for the folder.
            sub_folder (str, optional): The sub-folder name. Defaults to "".
            nthreads (int, optional): The number of threads to use. Defaults to 1.
            seed (int, optional): The seed value. Defaults to 42.
            m_type (str, optional): The type identifier. Defaults to "e".

        Returns:
            None
        """
        self.work_path = os.path.join(input_path, "")
        self.sub_folder = sub_folder
        self.full_path = os.path.join(self.work_path, self.sub_folder, "")
        self.nthreads = nthreads
        self.seed = seed
        self.m_type = m_type

        return

    def init_trades_from_path(self, work_path, sub_path):
        """
        Initializes the TRADES by calling a Fortran90 module, resetting the state, setting the work path, 
        and initializing various parameters related to the system, fitting, priors, RV trend, stellar mass, 
        and radius. It also computes physical parameters, sets up radial velocities and transit information,
        and prepares arrays for Kepler elements. No parameters are passed in. No explicit return value.
        """
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
        self.npar = (f90trades.npar).item()  # NUMBER OF SYSTEM PARAMETERS
        self.nfit = (f90trades.nfit).item()  # NUMBER OF PARAMETERS TO FIT
        self.nfree = (f90trades.nfree).item()  # NUMBER OF FREE PARAMETERS (ie nrvset)
        self.dof = (f90trades.dof).item()  # NUMBER OF DEGREES OF FREEDOM = NDATA - NFIT
        self.inv_dof = (f90trades.inv_dof).item()

        self.transit_flag = f90trades.get_do_transit_flag(self.n_bodies)

        self.tstart = (f90trades.tstart).item()
        self.tepoch = (f90trades.tepoch).item()
        self.tint = (f90trades.tint).item()

        # fitting parameters
        str_len = f90trades.str_len
        temp_names = f90trades.get_parameter_names(self.nfit, str_len)
        # names
        self.fitting_names = anc.convert_fortran_charray2python_strararray(temp_names)
        self.fitting_parameters, self.fitting_minmax = (
            f90trades.get_default_fitting_parameters_and_minmax(self.nfit)
        )
        (
            self.system_parameters,
            self.system_parameters_min,
            self.system_parameters_max,
        ) = f90trades.get_system_parameters_and_minmax(self.npar)

        self.all_names = anc.convert_fortran_charray2python_strararray(
            f90trades.get_all_keplerian_names(self.npar, str_len)
        )
        self.to_fit = f90trades.get_tofit(self.npar)
        sel_fixed = self.to_fit == 0
        self.fixed_names = anc.decode_list(np.array(self.all_names)[sel_fixed])
        self.fixed_parameters = self.system_parameters[sel_fixed]

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
        self.n_t0, self.pephem, self.tephem = f90trades.get_observed_transits_info(
            self.n_planets
        )

        self.n_t0_sum = f90trades.ntts
        self.n_set_t0 = 0
        for i in range(0, self.n_planets):
            if self.n_t0[i] > 0:
                self.n_set_t0 += 1

        self.wrttime = f90trades.wrttime

        return

    def init_trades(self):
        """
        Initializes the TRADES class instance by calling init_trades_from_path with work_path and sub_folder.
        No parameters.
        Returns None.
        """
        self.init_trades_from_path(self.work_path, self.sub_folder)

        return

    def set_one_fit_par_boundaries(self, ifit, min_val, max_val):
        """
        Sets the boundaries for a single fitting parameter.

        Parameters:
            self: The instance of the class.
            ifit: The index of the fitting parameter.
            min_val: The minimum value for the parameter.
            max_val: The maximum value for the parameter.

        Returns:
            None
        """
        # self.fitting_minmax[ifit, 0] = min_val
        # self.fitting_minmax[ifit, 1] = max_val
        self.fitting_minmax = f90trades.set_one_fit_par_boundaries(
            ifit + 1, min_val, max_val, self.fitting_minmax
        )

        return

    def reset_all_fit_boundaries(self):
        """
        Resets all fit boundaries using f90trades.reset_all_fit_boundaries() and assigns the result to self.fitting_minmax.
        No parameters.
        Returns None.
        """
        # reset_all_fit_boundaries()
        # self.fitting_minmax = f90trades.parameters_minmax.copy()
        self.fitting_minmax = f90trades.reset_all_fit_boundaries()

        return

    def path_change(self, new_path):
        """
        Changes the path using f90trades.path_change.

        Parameters:
            self: The instance of the class.
            new_path: The new path to change to.

        Returns:
            None
        """
        f90trades.path_change(os.path.join(new_path, ""))

        return

    def write_time_change(self, new_wrttime):
        """
        Updates the write time for TRADES and the local object.

        Parameters:
            new_wrttime: The new write time value.

        Returns:
            None
        """

        f90trades.wrttime = new_wrttime
        self.wrttime = f90trades.wrttime

        return

    def print_summary_boundaries(self, par_print="all", return_boundaries=False):
        """
        Print summary boundaries based on the specified parameter print options.

        Parameters:
            par_print (str): Specifies which boundaries to print. Default is "all".

        Returns:
            None
        """

        if return_boundaries:
            out = {}
        else:
            out = None

        sys.stdout.flush()
        if ("fix" in par_print) or ("all" in par_print):
            print("# fixed value", flush=True)
            for n, p in zip(self.fixed_names, self.fixed_parameters):
                print("{:12s} = {:16.8f}".format(n, p), flush=True)
            if return_boundaries:
                out["fixed"] = {n: p for n, p in zip(self.fixed_names, self.fixed_parameters)}
        if ("fit" in par_print) or ("all" in par_print):
            print("# fitted min -- max", flush=True)
            for n, p in zip(self.fitting_names, self.fitting_minmax):
                print("{:12s} = {:16.8f} -- {:16.8f}".format(n, p[0], p[1]), flush=True)
            if return_boundaries:
                out["fitted"] = {n: p for n, p in zip(self.fitting_names, self.fitting_minmax)}

        if ("sys" in par_print) or ("all" in par_print):
            print("# system min -- max", flush=True)
            for n, pl, pu in zip(
                self.all_names, self.system_parameters_min, self.system_parameters_max
            ):
                print("{:12s} = {:16.8f} -- {:16.8f}".format(n, pl, pu), flush=True)
            if return_boundaries:
                out["physical"] = {n: [p1, p2] for n, p1, p2 in zip(self.all_names, self.system_parameters_min, self.system_parameters_max)}

        print("", flush=True)

        return out

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
        """
        Update system and fitting parameters based on keplerian input parameters.
        Parameters:
            mass: Mass of the system
            radius: Radius of the system
            period: Period of the system
            ecc: Eccentricity of the system
            argp: Argument of pericenter of the system
            meanA: Mean anomaly of the system
            inc: Inclination of the system
            longN: Longitude of the ascending node of the system
        Returns:
            None
        """
        n_all_par = len(self.system_parameters)
        all_par = np.zeros((n_all_par))
        fit_par = np.zeros((self.nfit))
        (
            all_par,
            fit_par,
        ) = f90trades.update_parameters_from_keplerian(
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
        """
        Update system and fitting parameters based on default keplerian parameters.
        Parameters:
            mass: Mass of the system
            radius: Radius of the system
            period: Period of the system
            ecc: Eccentricity of the system
            argp: Argument of pericenter of the system
            meanA: Mean anomaly of the system
            inc: Inclination of the system
            longN: Longitude of the ascending node of the system
        Returns:
            None
        """
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
        """
        Run and get statistics from the input parameters.

        Parameters:
            fit_pars: List of fitting parameters for the function

        Returns:
            chi_square: Chi-square value
            reduced_chi_square: Reduced chi-square value
            lgllhd: Log-likelihood value
            lnprior: Log-prior value
            ln_const: Log-constant value
            bic: Bayesian Information Criterion value
            check: Check value
        """
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
        """
        Run and get statistics from the input parameters.
        
        Returns:
            chi_square: Chi-square value
            reduced_chi_square: Reduced chi-square value
            lgllhd: Log-likelihood value
            lnprior: Log-prior value
            ln_const: Log-constant value
            bic: Bayesian Information Criterion value
            check: Check value
        """
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
        """
        Run TRADES to get all TTs and RVs during the integration time.
        
        Parameters:
            fit_par: The input parameters for the function.
        
        Returns:
            Tuple containing ttra_full, dur_full, lambda_rm_full, id_ttra_full, stats_ttra, time_rv_nmax, rv_nmax, stats_rv.
        """
        # run TRADES to get all TTs and RVs during the integration time
        nt0_full, nrv_nmax = f90trades.get_ntt_nrv(fit_par)
        # anc.print_both("full nRV = {:d} and nT0s = {}".format(nrv_nmax, nt0_full))
        (
            ttra_full,
            dur_full,
            lambda_rm_full,
            id_ttra_full,
            stats_ttra,
            time_rv_nmax,
            rv_nmax,
            stats_rv,
        ) = f90trades.fit_par_to_ttra_rv(fit_par, nt0_full, nrv_nmax)

        return (
            ttra_full,
            dur_full,
            lambda_rm_full,
            id_ttra_full,
            stats_ttra,
            time_rv_nmax,
            rv_nmax,
            stats_rv,
        )

    def save_models_from_parameters(self, fit_par, smp_h5, smp_name):
        """
        Run TRADES to save models based on the input parameters.
        
        Parameters:
            fit_par: The fitting parameters.
            smp_h5: The HDF5 file to save the models.
            smp_name: The name of the sample.
        
        Returns:
            None
        """
        (
            ttra_full,
            dur_full,
            lambda_rm_full,
            id_ttra_full,
            stats_ttra,
            time_rv_nmax,
            rv_nmax,
            stats_rv,
        ) = self.get_all_observables_from_parameters(fit_par)

        gr = smp_h5.create_group(smp_name)
        gr.create_dataset(
            "fitting_parameters", data=fit_par, dtype=float, compression="gzip"
        )
        gr["fitting_parameters"].attrs["fitting_names"] = self.fitting_names
        gr.attrs["tepoch"] = self.tepoch

        stats_rv = np.array(stats_rv).astype(bool)
        time_rvx, rvx = time_rv_nmax[stats_rv], rv_nmax[stats_rv]
        id_rv = np.argsort(time_rvx)
        time_rv, rv = time_rvx[id_rv], rvx[id_rv]

        # anc.print_both("time_rv min = {:.5f} max = {:.5f}".format(np.min(time_rv), np.max(time_rv)))
        gr.create_dataset(
            "time_rv_mod", data=time_rv, dtype=float, compression="gzip"
        )
        gr.create_dataset("rv_mod", data=rv, dtype=float, compression="gzip")

        stats_ttra = np.array(stats_ttra).astype(bool)
        for ipl in range(2, self.n_bodies + 1):
            sel_tra = np.logical_and(
                np.logical_and(id_ttra_full == ipl, stats_ttra), ttra_full > -9.0e10
            )
            ttra = np.array(ttra_full)[sel_tra]
            dur = np.array(dur_full)[sel_tra]
            lambda_rm = lambda_rm_full[sel_tra]
            idx_tra = np.argsort(ttra)
            ttra = ttra[idx_tra]
            dur = dur[idx_tra]
            lambda_rm = lambda_rm[idx_tra]
            gr.create_dataset(
                "TTs_{:d}".format(ipl), data=ttra, dtype=float, compression="gzip"
            )
            gr.create_dataset(
                "T41s_{:d}".format(ipl), data=dur, dtype=float, compression="gzip"
            )
            gr.create_dataset(
                "lambda_rm_{:d}".format(ipl),
                data=lambda_rm,
                dtype=float,
                compression="gzip",
            )
            # if np.sum(sel_tra) > 0:
            #     anc.print_both("T0_{:02d} min = {:.5f} max = {:.5f}".format(ipl, np.min(ttra), np.max(ttra)))

        return

    def fit_pars_to_keplerian_elements(self, fit_pars):
        """
        Converts fitted parameters to Keplerian elements; returns mass, radius, period, semi-major axis, eccentricity, argument of pericenter, mean anomaly, inclination, and longitude of nodes.
        """
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

        return (mass, radius, period, sma, ecc, argp, meanA, inc, longN)

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
        """
        Runs the simulation and writes summary files based on the given parameters.

        Parameters:
            sim_name (str): The name of the simulation.
            fit_pars (array): Fitted parameters for the simulation.
            phy_pars (array): Physical parameters for the simulation.
            id_sim (int): The simulation ID (default is 1).
            sigma_fit (array): Sigma for fitted parameters (default is None).
            sigma_phy (array): Sigma for physical parameters (default is None).
            par_description (str): Description of the parameters (default is an empty string).
            return_stats (bool): Flag to return statistics (default is False).

        Returns:
            Tuple: If return_stats is True, returns chi_square, reduced_chi_square, lgllhd, lnprior, ln_const, bic, and check. Otherwise, returns None.
        """
        print("RUNNING sim {}".format(sim_name))
        # print("with parameters:")
        # print("{}".format(fit_pars))

        out_folder = os.path.join(self.full_path, sim_name, "")
        os.makedirs(out_folder, exist_ok=True)
        f90trades.path_change(out_folder)

        for iname, name in enumerate(self.fitting_names):
            if (
                name[0] == "w"
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
                    "# WRITING WARNING FILE: {:s}".format(
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

    def computes_observables_from_keplerian_elements(
        self,
        t_epoch,
        t_start,
        t_int,
        mass,
        radius,
        period,
        ecc,
        argp,
        meana,
        inc,
        longn,
        transit_flag=None,
    ):
        """
        A function that computes observables from Keplerian elements.

        Parameters:
            self: the object instance
            t_epoch: the epoch time
            t_start: the start time
            t_int: the integration time
            mass: the mass
            radius: the radius
            period: the period
            ecc: the eccentricity
            argp: the argument of pericenter
            meana: the mean anomaly
            inc: the inclination
            longn: the longitude of nodes
            transit_flag: the transit flag (default is None)

        Returns:
            Tuple containing computed observables:
                - rv_sim: the simulated radial velocity
                - body_t0_sim: the simulated body T0
                - epo_sim: the simulated epoch
                - t0_sim: the simulated T0
                - t14_sim: the simulated T14
                - lambda_rm_sim: the simulated lambda_rm
                - kel_sim: the simulated kel
                - stable: the stability flag
        """
        if transit_flag is None:
            transit_flag = self.transit_flag

        n_kep = 8
        (
            rv_sim,
            body_t0_sim,
            epo_sim,
            t0_sim,
            t14_sim,
            lambda_rm_sim,
            kel_sim,
            stable,
        ) = f90trades.kelements_to_rv_and_t0s(
            t_epoch,
            t_start,
            t_int,
            mass,
            radius,
            period,
            ecc,
            argp,
            meana,
            inc,
            longn,
            transit_flag,
            self.n_rv,
            self.n_t0_sum,
            n_kep,
        )

        return (
            rv_sim,
            body_t0_sim,
            epo_sim,
            t0_sim,
            t14_sim,
            lambda_rm_sim,
            kel_sim,
            stable,
        )

    def computes_observables_from_default_keplerian_elements(
        self, t_epoch, t_start, t_int, transit_flag=None
    ):
        """
        A function that computes observables from default Keplerian elements.

        Parameters:
            self: the object instance
            t_epoch: the epoch time
            t_start: the start time
            t_int: the integration time
            transit_flag: the transit flag (default is None)

        Returns:
            Tuple containing computed observables:
                - rv_sim: the simulated radial velocity
                - body_t0_sim: the simulated body T0
                - epo_sim: the simulated epoch
                - t0_sim: the simulated T0
                - t14_sim: the simulated T14
                - lambda_rm_sim: the simulated lambda_rm
                - kel_sim: the simulated kel
                - stable: the stability flag
        """
        if transit_flag is None:
            transit_flag = self.transit_flag

        (
            rv_sim,
            body_t0_sim,
            epo_sim,
            t0_sim,
            t14_sim,
            lambda_rm_sim,
            kel_sim,
            stable,
        ) = self.computes_observables_from_keplerian_elements(
            t_epoch,
            t_start,
            t_int,
            self.mass,
            self.radius,
            self.period,
            self.ecc,
            self.argp,
            self.meana,
            self.inc,
            self.longn,
            in_transit_flag=transit_flag,
        )

        return (
            rv_sim,
            body_t0_sim,
            epo_sim,
            t0_sim,
            t14_sim,
            lambda_rm_sim,
            kel_sim,
            stable,
        )

    def reset(self):

        f90trades.deallocate_variables()

        return


# =============================================================================
def base_plot_orbits(
    time_steps,
    orbits,
    radius,
    n_body,
    body_names,
    figsize=(4, 4),
    sky_scale="star",
    side_scale="star",
    title=None,
    show_plot=False,
):
    """
    Generate a plot displaying orbits in the sky, orbit, and side planes.

    Parameters:
    - time_steps: array of time steps
    - orbits: array of orbital positions
    - radius: array of radii
    - n_body: number of bodies
    - body_names: array of body names
    - figsize: tuple specifying figure size (default is (4, 4))
    - sky_scale: string specifying sky scale (default is "star")
    - side_scale: string specifying side scale (default is "star")
    - title: string specifying title of the plot (default is None)
    - show_plot: boolean to determine if the plot should be displayed (default is False)

    Returns:
    - fig: the generated plot figure
    """
    n_steps = len(time_steps)
    # base_colors = plt.get_cmap("nipy_spectral")(
    #     np.linspace(0.1, 0.9, endpoint=True, num=n_body - 1)
    # )
    base_colors = anc.set_colors(
        n_body - 1, vmin=0.05, vmax=0.95, colormap="nipy_spectral"
    )

    # colors with alpha based on increasing time (alpha(t_start) = 0.2, alpha(t_end) = 1.0)
    alphas = np.linspace(0.05, 0.95, endpoint=True, num=n_steps)
    colors = []
    for c in base_colors:
        colors.append([[c[0], c[1], c[2], a] for a in alphas])

    idx_X = [i * 6 for i in range(1, n_body)]
    idx_Y = [i + 1 for i in idx_X]
    idx_Z = [i + 1 for i in idx_Y]

    X = orbits[:, idx_X]
    Y = orbits[:, idx_Y]
    Z = orbits[:, idx_Z]

    rstar = radius[0] * cst.RsunAU
    scale = 1.50
    alpha_star = 1.0
    # sky_scale = "star" # "X"
    # side_scale = "star"
    leg_size = 4
    lab_size = 12
    tic_size = 7

    zo_s = 2
    zo_0 = 3
    zo_z = zo_s + n_body

    ms_scatter = 5
    ms_z = np.sqrt(ms_scatter) + 0.1
    # ---------------------------
    # X vs Y == sky plane
    fig = plt.figure(figsize=figsize)
    axs = []
    nrows, ncols = 2, 2
    if title is not None:
        irow, icol = 0, 1
        ax = plt.subplot2grid((nrows, ncols), (irow, icol))
        ax.set_axis_off()
        ax.text(0.5, 0.5, title, ha="center", va="center", fontsize=lab_size)

    irow, icol = 0, 0
    ax = plt.subplot2grid((nrows, ncols), (irow, icol))
    ax.tick_params(axis="both", labelsize=tic_size)
    ax.tick_params(axis="x", labelrotation=45)
    ax.axhline(0.0, color="black", ls="-", lw=0.6, zorder=zo_0)
    ax.axvline(0.0, color="black", ls="-", lw=0.6, zorder=zo_0)
    star = plt.Circle(
        (0, 0),
        radius=rstar,
        facecolor="C1",
        zorder=zo_s,
        edgecolor="None",
        alpha=alpha_star,
    )
    ax.add_artist(star)

    leg_elem = []
    for i_pl in range(0, n_body - 1):

        # ax.scatter(
        #     X[:, i_pl],
        #     Y[:, i_pl],
        #     c=colors[i_pl],
        #     marker="o",
        #     s=ms_scatter,
        #     edgecolors="None",
        #     zorder=zo_z + i_pl + 1,
        # )

        Zpos = Z[:, i_pl] > 0.0
        ax.scatter(
            X[Zpos, i_pl],
            Y[Zpos, i_pl],
            c=np.array(colors[i_pl])[Zpos],
            marker="o",
            s=ms_scatter,
            edgecolors="None",
            zorder=zo_z + i_pl + 1,
        )
        ax.plot(
            X[Zpos, i_pl],
            Y[Zpos, i_pl],
            marker="o",
            ms=ms_z,
            mfc="None",
            mec="black",
            mew=0.4,
            ls="",
            zorder=zo_z + i_pl + 1,
        )
        lp = plt.Line2D(
            [0],
            [0],
            color=base_colors[i_pl],
            ls="",
            marker="o",
            ms=ms_z,
            mec="black",
            mew=0.3,
            label="planet {} Z > 0".format(body_names[i_pl + 1]),
        )
        leg_elem.append(lp)

        Zneg = Zpos == False
        ax.scatter(
            X[Zneg, i_pl],
            Y[Zneg, i_pl],
            c=np.array(colors[i_pl])[Zneg],
            marker="o",
            s=ms_scatter,
            edgecolors="None",
            zorder=zo_s - 1,
        )
        ax.plot(
            X[Zneg, i_pl],
            Y[Zneg, i_pl],
            marker="o",
            ms=ms_z,
            mfc="None",
            mec="white",
            mew=0.3,
            ls="",
            # zorder=zo_z + i_pl + 0,
            zorder=zo_s - 1,
        )
        ln = plt.Line2D(
            [0],
            [0],
            color=base_colors[i_pl],
            ls="",
            marker="o",
            ms=ms_z,
            mec="white",
            mew=0.3,
            label="planet {} Z < 0".format(body_names[i_pl + 1]),
        )
        leg_elem.append(ln)

    ax.legend(
        handles=leg_elem,
        # loc='best',
        bbox_to_anchor=(1.05, 1),
        loc="upper left",
        borderaxespad=0.0,
        fontsize=leg_size,
        facecolor=[0.8] * 3 + [1.0],  #'lightgray',
        edgecolor="None",
    )
    ax.set_title("sky plane", fontsize=tic_size)
    ax.set_xlabel("X (au)")
    ax.set_ylabel("Y (au)")
    if sky_scale.lower() == "star":
        lims = [-scale * rstar, +scale * rstar]
        xlims = lims
        ylims = lims
    else:
        xlims = ax.get_xlim()
        ylims = xlims

    ax.set_xlim(xlims)
    ax.set_ylim(ylims)

    axs.append(ax)
    # ---------------------------
    # X vs Z == orbit plane
    irow, icol = 1, 0

    ax = plt.subplot2grid((nrows, ncols), (irow, icol))
    ax.tick_params(axis="both", labelsize=tic_size)
    ax.tick_params(axis="x", labelrotation=45)
    ax.axhline(0.0, color="black", ls="-", lw=0.7, zorder=1)
    ax.axvline(0.0, color="black", ls="-", lw=0.7, zorder=1)
    zo_s = 3
    star = plt.Circle(
        (0, 0),
        radius=rstar,
        facecolor="C1",
        zorder=zo_s,
        edgecolor="None",
        alpha=alpha_star,
    )
    ax.add_artist(star)

    for i_pl in range(0, n_body - 1):
        ax.scatter(
            X[:, i_pl],
            Z[:, i_pl],
            c=colors[i_pl],
            marker="o",
            s=ms_scatter,
            edgecolors="None",
            zorder=zo_z + i_pl + 1,
        )
    ax.set_title("orbits plane - observer at top", fontsize=tic_size)
    ax.set_xlabel("X (au)")
    ax.set_ylabel("Z (au)")

    axs.append(ax)

    # ---------------------------
    # Z vs Y == side plane
    irow, icol = 1, 1

    ax = plt.subplot2grid((nrows, ncols), (irow, icol))
    ax.tick_params(axis="both", labelsize=tic_size)
    ax.tick_params(axis="x", labelrotation=45)
    ax.axhline(0.0, color="black", ls="-", lw=0.7, zorder=1)
    ax.axvline(0.0, color="black", ls="-", lw=0.7, zorder=1)
    zo_s = 3
    star = plt.Circle(
        (0, 0),
        radius=rstar,
        facecolor="C1",
        zorder=zo_s,
        edgecolor="None",
        alpha=alpha_star,
    )
    ax.add_artist(star)

    for i_pl in range(0, n_body - 1):
        ax.scatter(
            Z[:, i_pl],
            Y[:, i_pl],
            c=colors[i_pl],
            marker="o",
            s=ms_scatter,
            edgecolors="None",
            # label="body {}".format(i_pl+2),
            zorder=zo_z + i_pl + 1,
        )
        Zpos = Z[:, i_pl] > 0.0
        # ax.plot(
        #     Z[Zpos, i_pl],
        #     Y[Zpos, i_pl],
        #     marker="o",
        #     ms=ms_z,
        #     mfc="None",
        #     mec="black",
        #     mew=0.4,
        #     ls="",
        #     zorder=zo_z + i_pl + 1,
        # )
        Zneg = Zpos == False
        # ax.plot(
        #     Z[Zneg, i_pl],
        #     Y[Zneg, i_pl],
        #     marker="o",
        #     ms=ms_z,
        #     mfc="None",
        #     mec="white",
        #     mew=0.3,
        #     ls="",
        #     # zorder=zo_z + i_pl + 0,
        #     zorder=zo_s - 1,
        # )
    # ax.legend(loc='best', fontsize=leg_size)
    ax.set_title("side plane - observer at right", fontsize=tic_size)
    ax.set_xlabel("Z (au)")
    ax.set_ylabel("Y (au)")
    if side_scale.lower() == "star":
        lims = [-scale * rstar, +scale * rstar]
        xlims = lims
        ylims = lims
    elif "pos" in side_scale.lower():
        dZp = np.ptp(Z[Zpos, :])
        xlims = [-scale * rstar, np.max(Z)+ dZp*0.03]
        ylims = [-scale * rstar, +scale * rstar]
    else:
        xlims = ax.get_xlim()
        # ylims = ax.get_ylim()
        ylims = xlims

    ax.set_xlim(xlims)
    ax.set_ylim(ylims)

    axs.append(ax)
    fig.align_ylabels(axs)
    plt.tight_layout()

    if show_plot:
        plt.show()
    # plt.close(fig)
    return fig


# =============================================================================
# =============================================================================
def angular_momentum_deficit_posterior(n_bodies, post_fit, all_pars):
    """
    Calculates the posterior for the angular momentum deficit given the number of bodies, post-fit data, and all parameters.

    Parameters:
    - n_bodies (int): The total number of bodies.
    - post_fit (numpy.ndarray): The post-fit data.
    - all_pars (list): A list of all parameters.

    Returns:
    - amd_names (list): A list of names corresponding to different parameters.
    - amd_full (numpy.ndarray): A 2D array containing calculated values related to the angular momentum deficit.
    - amd_stable: The stability of the angular momentum deficit.
    """
    n_post, n_fit = np.shape(post_fit)
    n_all = len(all_pars)
    n_pairs = n_bodies - 2

    lambdas_bodies, amd_bodies, amd, amd_r_pairs, amd_h_pairs, amd_stable = (
        f90trades.angular_momentum_deficit_posterior(
            n_bodies,
            # n_post, n_fit,
            post_fit,
            # n_all,
            all_pars,
            n_pairs,
        )
    )
    namd = amd / np.sum(lambdas_bodies[:, 1:], axis=1)

    amd_names = (
        ["Lambda{}".format(i + 1) for i in range(1, n_bodies)]
        + ["AMD{}".format(i + 1) for i in range(1, n_bodies)]
        + ["AMD", "NAMD"]
        + ["rAMD_pair{}".format(i + 1) for i in range(0, n_pairs)]
        + ["hAMD_pair{}".format(i + 1) for i in range(0, n_pairs)]
    )

    amd_full = np.column_stack(
        (lambdas_bodies[:, 1:], amd_bodies[:, 1:], amd, namd, amd_r_pairs, amd_h_pairs)
    )

    return amd_names, amd_full, amd_stable


# =============================================================================
# =============================================================================
# rename some functions for backward compatibility
kelements_to_rv_and_t0s = kelements_to_observed_rv_and_t0s
kelements_to_rv = kelements_to_observed_rv
kelements_to_t0s = kelements_to_observed_t0s
# =============================================================================
