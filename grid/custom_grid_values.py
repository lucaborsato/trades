import numpy as np
import os
import sys
import scipy.stats as sts

# ==============================================================================
# custom modules
script_path = os.path.realpath(__file__)
module_path = os.path.abspath(os.path.join(os.path.dirname(script_path), "../pytrades"))
sys.path.append(module_path)

import constants as cst


def set_grid_values_01(n_sims=1):

    # DICT WITH MULTIPLE-PLANET SYSTEMS
    sort_ecc = [False, False, True]
    grid_values = {}

    # STAR
    body_id = 0  # 0 == star, 1 ==  planet b, 2 == perturber x
    mass = np.zeros((n_sims)) + 1.066
    radius = np.zeros((n_sims)) + 0.962
    period = np.zeros((n_sims))
    ecc = np.zeros((n_sims))
    if sort_ecc[body_id]:
        ecc = np.sort(ecc)
    argp = np.zeros((n_sims))
    meanA = np.zeros((n_sims))
    inc = np.zeros((n_sims))
    longN = np.zeros((n_sims))
    params = np.column_stack((mass, radius, period, ecc, argp, meanA, inc, longN))
    grid_values[body_id] = params

    # PLANET b
    body_id = 1  # 0 == star, 1 ==  planet b, 2 == perturber x
    mass = np.zeros((n_sims)) + 2.244 * cst.Mjups
    radius = np.zeros((n_sims)) + 1.11237 * cst.Rjups
    period = np.zeros((n_sims)) + 8.15874
    ecc = np.zeros((n_sims)) + 0.31
    if sort_ecc[body_id]:
        ecc = np.sort(ecc)
    argp = np.zeros((n_sims)) + 94.27
    meanA = np.zeros((n_sims)) + (360.0 * 2109.52371 / 8.15874) % 360.0
    inc = np.zeros((n_sims)) + 88.40
    longN = np.zeros((n_sims)) + 180.0
    params = np.column_stack((mass, radius, period, ecc, argp, meanA, inc, longN))
    grid_values[body_id] = params

    # PLANET x by DGa
    body_id = 2  # 0 == star, 1 ==  planet b, 2 == perturber x
    mass = np.random.normal(loc=0.1 * cst.Mjups, scale=1.0 * cst.Mears, size=n_sims)
    radius = np.zeros((n_sims)) + 0.2
    period = np.random.normal(loc=32.8, scale=0.01, size=n_sims)
    ecc = np.random.random(n_sims)
    if sort_ecc[body_id]:
        ecc = np.sort(ecc)
    argp = np.random.random(n_sims) * 360.0
    meanA = np.random.random(n_sims) * 360.0
    inc = 30.0 + np.random.random(n_sims) * 120.0
    longN = np.zeros((n_sims)) + 180.0
    params = np.column_stack((mass, radius, period, ecc, argp, meanA, inc, longN))
    grid_values[body_id] = params

    return grid_values


# =============================================================================


def set_grid_values_02(n_sims=1):

    # DICT WITH MULTIPLE-PLANET SYSTEMS
    sort_ecc = [False, False, True]
    grid_values = {}

    # STAR
    body_id = 0  # 0 == star, 1 ==  planet b, 2 == perturber x
    mass = np.zeros((n_sims)) + 1.066
    radius = np.zeros((n_sims)) + 0.962
    period = np.zeros((n_sims))
    ecc = np.zeros((n_sims))
    if sort_ecc[body_id]:
        ecc = np.sort(ecc)
    argp = np.zeros((n_sims))
    meanA = np.zeros((n_sims))
    inc = np.zeros((n_sims))
    longN = np.zeros((n_sims))
    params = np.column_stack((mass, radius, period, ecc, argp, meanA, inc, longN))
    grid_values[body_id] = params

    # PLANET b
    body_id = 1  # 0 == star, 1 ==  planet b, 2 == perturber x
    mass = np.abs(
        np.random.normal(loc=2.13 * cst.Mjups, scale=0.11 * cst.Mjups, size=n_sims)
    )
    radius = np.zeros((n_sims)) + 1.11237 * cst.Rjups
    period = np.zeros((n_sims)) + 8.15874
    ecc = np.zeros((n_sims)) + 0.312
    if sort_ecc[body_id]:
        ecc = np.sort(ecc)
    argp = np.zeros((n_sims)) + (274.5 - 180.0)
    meanA = (
        np.zeros((n_sims))
        + (360.0 * (0.0 - (2459105.62952 - cst.btjd)) / period[0]) % 360.0
    )
    inc = np.zeros((n_sims)) + 88.40
    # mass from RV analysis, that is mass * sin(inc)
    mass = mass / np.sin(inc * cst.deg2rad)
    longN = np.zeros((n_sims)) + 180.0
    params = np.column_stack((mass, radius, period, ecc, argp, meanA, inc, longN))
    grid_values[body_id] = params

    # PLANET x by DGa
    body_id = 2  # 0 == star, 1 ==  planet b, 2 == perturber x
    mass = np.abs(np.random.normal(loc=0.101 * cst.Mjups, scale=0.014, size=n_sims))
    radius = np.zeros((n_sims)) + 0.2
    period = np.random.normal(loc=32.827, scale=0.04, size=n_sims)
    ecc = np.random.random(n_sims)
    if sort_ecc[body_id]:
        ecc = np.sort(ecc)
    argp = np.random.random(n_sims) * 360.0
    tperi = np.random.normal(loc=4504.0, scale=1.7, size=n_sims)  # days, input
    meanA = (360.0 * (0.0 - tperi) / period) % 360.0
    inc = 30.0 + np.random.random(n_sims) * 120.0
    # mass from RV analysis, that is mass * sin(inc)
    mass = mass / np.sin(inc * cst.deg2rad)
    longN = np.zeros((n_sims)) + 180.0
    params = np.column_stack((mass, radius, period, ecc, argp, meanA, inc, longN))
    grid_values[body_id] = params

    return grid_values


# =============================================================================
def get_truncnorm(l=0.0, s=1.0, vmin=-100.0, vmax=100.0, size=1):

    a, b = (vmin - l) / s, (vmax - l) / s
    r = sts.truncnorm.rvs(a, b, loc=l, scale=s, size=size)

    return r


# =============================================================================


def set_grid_values_03(n_sims=1):

    # DICT WITH MULTIPLE-PLANET SYSTEMS
    sort_ecc = [False, False, False]
    grid_values = {}

    # STAR
    body_id = 0  # 0 == star, 1 ==  planet b, 2 == perturber x
    mass = get_truncnorm(
        l=1.066, s=0.044, vmin=1.0e-2, vmax=2.0, size=n_sims
    )  # bounded to 0.01-2 Msun
    radius = get_truncnorm(
        l=0.962, s=0.030, vmin=1.0e-2, vmax=2.0, size=n_sims
    )  # bounded to 0.01-2 Rsun
    period = np.zeros((n_sims))
    ecc = np.zeros((n_sims))
    if sort_ecc[body_id]:
        ecc = np.sort(ecc)
    argp = np.zeros((n_sims))
    meanA = np.zeros((n_sims))
    inc = np.zeros((n_sims))
    longN = np.zeros((n_sims))
    params = np.column_stack((mass, radius, period, ecc, argp, meanA, inc, longN))
    grid_values[body_id] = params

    # PLANET b
    body_id = 1  # 0 == star, 1 ==  planet b, 2 == perturber x
    mass = get_truncnorm(
        l=2.13 * cst.Mjups,
        s=0.11 * cst.Mjups,
        vmin=1.0 * cst.Mears,
        vmax=10.0 * cst.Mjups,
        size=n_sims,
    )  # bounded to 1 Mearth and 10 Mjups
    radius = get_truncnorm(
        l=1.112 * cst.Rjups,
        s=0.035 * cst.Rjups,
        vmin=0.01 * cst.Rears,
        vmax=10.0 * cst.Rjups,
        size=n_sims,
    )  # np.zeros((n_sims)) + 1.11237 * cst.Rjups # Rp =  1.112 +/- 0.035 Rjup
    period = get_truncnorm(l=8.1587364, s=0.000003, vmin=0.0, vmax=1.0e5, size=n_sims)
    ecc = get_truncnorm(
        l=0.312, s=0.003, vmin=0.0, vmax=1.0, size=n_sims
    )  # 0.312 +/- 0.003
    if sort_ecc[body_id]:
        ecc = np.sort(ecc)
    argp = np.zeros((n_sims)) + (274.5 - 180.0)
    meanA = (
        np.zeros((n_sims))
        + (360.0 * (0.0 - (2459105.62952 - cst.btjd)) / 8.1587364) % 360.0
    )
    inc = get_truncnorm(
        l=88.40, s=0.02, vmin=60.0, vmax=120.0, size=n_sims
    )  # bad for some reason ...
    # inc = 88.40 + sts.norm.rvs(size=n_sims)*0.02
    # mass from RV analysis, that is mass * sin(inc)
    mass = mass / np.sin(inc * cst.deg2rad)
    longN = np.zeros((n_sims)) + 180.0
    params = np.column_stack((mass, radius, period, ecc, argp, meanA, inc, longN))
    grid_values[body_id] = params

    # PLANET x by DGa
    body_id = 2  # 0 == star, 1 ==  planet b, 2 == perturber x
    mass = get_truncnorm(
        l=0.101 * cst.Mjups,
        s=0.014 * cst.Mjups,
        vmin=1.0 * cst.Mears,
        vmax=10.0 * cst.Mjups,
        size=n_sims,
    )
    radius = np.zeros((n_sims)) + 0.2
    period = np.random.normal(loc=32.827, scale=0.04, size=n_sims)
    # ecc    = get_truncnorm(l=0.09, s=0.02, vmin=0.0, vmax=1.0, size=n_sims) # Van Eylen 2019: half-Gaussian single: 0.32 +/- 0.06; multi: 0.083 +/- 0.020
    ecc = get_truncnorm(
        l=0.32, s=0.06, vmin=0.0, vmax=1.0, size=n_sims
    )  # Van Eylen 2019: half-Gaussian single: 0.32 +/- 0.06; multi: 0.083 +/- 0.020
    if sort_ecc[body_id]:
        ecc = np.sort(ecc)
    argp = np.random.random(n_sims) * 360.0
    tperi = np.random.normal(loc=4504.0, scale=1.7, size=n_sims)  # days, input
    meanA = (360.0 * (0.0 - tperi) / period) % 360.0
    # inc    = get_truncnorm(l=88.40, s=0.02, vmin=30.0, vmax=150.0, size=n_sims) # 30.0 + np.random.random(n_sims) * 120.0
    inc = get_truncnorm(l=90.0, s=2.0, vmin=30.0, vmax=150.0, size=n_sims)
    # mass from RV analysis, that is mass * sin(inc)
    mass = mass / np.sin(inc * cst.deg2rad)
    longN = np.zeros((n_sims)) + 180.0
    params = np.column_stack((mass, radius, period, ecc, argp, meanA, inc, longN))
    grid_values[body_id] = params

    return grid_values


# =============================================================================


def set_grid_values_v1298Tau_01(n_sims=1):

    # DICT WITH MULTIPLE-PLANET SYSTEMS
    sort_ecc = [False, False, False, False, False]
    grid_values = {}

    # STAR
    body_id = 0  # 0 == star, 1 ==  planet b, 2 == perturber x
    mass = get_truncnorm(
        l=1.17, s=0.06, vmin=1.0e-2, vmax=2.0, size=n_sims
    )  # bounded to 0.01-2 Msun
    radius = get_truncnorm(
        l=1.28, s=0.07, vmin=1.0e-2, vmax=2.0, size=n_sims
    )  # bounded to 0.01-2 Rsun
    period = np.zeros((n_sims))
    ecc = np.zeros((n_sims))
    if sort_ecc[body_id]:
        ecc = np.sort(ecc)
    argp = np.zeros((n_sims))
    meanA = np.zeros((n_sims))
    inc = np.zeros((n_sims))
    longN = np.zeros((n_sims))
    params = np.column_stack((mass, radius, period, ecc, argp, meanA, inc, longN))
    grid_values[body_id] = params

    # PLANET b
    body_id = 1  # 0 == star, 1 ==  planet b, 2 == perturber c, ...
    mass = get_truncnorm(
        l=0.64 * cst.Mjups,
        s=0.19 * cst.Mjups,
        vmin=0.01 * cst.Mears,
        vmax=10.0 * cst.Mjups,
        size=n_sims,
    )  # bounded to 0.01 Mearth and 10 Mjups
    radius = get_truncnorm(
        l=0.868 * cst.Rjups,
        s=0.056 * cst.Rjups,
        vmin=0.01 * cst.Rears,
        vmax=10.0 * cst.Rjups,
        size=n_sims,
    )
    period = get_truncnorm(l=24.1399, s=0.0015, vmin=0.0, vmax=1.0e5, size=n_sims)
    ecc = get_truncnorm(
        l=0.134, s=0.075, vmin=0.0, vmax=1.0, size=n_sims
    )
    if sort_ecc[body_id]:
        ecc = np.sort(ecc)
    argp = np.random.random(n_sims) * 360.0
    meanA = (360.0 * (0.0 - np.random.normal(7067.0486, 0.0015, size=n_sims)) / period) % 360.0
    inc = 88.7 + np.random.random(n_sims) * ((90.0-88.7)*2)
    longN = np.zeros((n_sims)) + 180.0
    params = np.column_stack((mass, radius, period, ecc, argp, meanA, inc, longN))
    grid_values[body_id] = params

    # PLANET c
    body_id = 2 
    mass = 0.01*cst.Mears + np.random.random(n_sims) * (0.24*cst.Mjups-0.01*cst.Mears)
    radius = get_truncnorm(
        l=0.460 * cst.Rjups,
        s=0.034 * cst.Rjups,
        vmin=0.01 * cst.Rears,
        vmax=10.0 * cst.Rjups,
        size=n_sims,
    )
    period = get_truncnorm(l=8.24892, s=0.00083, vmin=0.0, vmax=1.0e5, size=n_sims)
    ecc = np.random.random(n_sims) * 0.3
    if sort_ecc[body_id]:
        ecc = np.sort(ecc)
    argp = np.random.random(n_sims) * 360.0
    meanA = (360.0 * (0.0 - np.random.normal(7064.2801, 0.0041, size=n_sims)) / period) % 360.0
    inc = 87.5 + np.random.random(n_sims) * ((90.0-87.5)*2)
    longN = np.zeros((n_sims)) + 180.0
    params = np.column_stack((mass, radius, period, ecc, argp, meanA, inc, longN))
    grid_values[body_id] = params

    # PLANET d
    body_id = 3
    mass = 0.01*cst.Mears + np.random.random(n_sims) * (0.31*cst.Mjups-0.01*cst.Mears)
    radius = get_truncnorm(
        l=0.574 * cst.Rjups,
        s=0.041 * cst.Rjups,
        vmin=0.01 * cst.Rears,
        vmax=10.0 * cst.Rjups,
        size=n_sims,
    )
    period = get_truncnorm(l=12.4058, s=0.0018, vmin=0.0, vmax=1.0e5, size=n_sims)
    ecc = np.random.random(n_sims) * 0.2
    if sort_ecc[body_id]:
        ecc = np.sort(ecc)
    argp = np.random.random(n_sims) * 360.0
    meanA = (360.0 * (0.0 - np.random.normal(7072.3907, 0.0063, size=n_sims)) / period) % 360.0
    inc = 88.3 + np.random.random(n_sims) * ((90.0-88.3)*2)
    longN = np.zeros((n_sims)) + 180.0
    params = np.column_stack((mass, radius, period, ecc, argp, meanA, inc, longN))
    grid_values[body_id] = params

    # PLANET e
    body_id = 4
    mass = get_truncnorm(
        l=1.16 * cst.Mjups,
        s=0.30 * cst.Mjups,
        vmin=0.01 * cst.Mears,
        vmax=10.0 * cst.Mjups,
        size=n_sims,
    )  # bounded to 0.01 Mearth and 10 Mjups
    radius = get_truncnorm(
        l=0.735 * cst.Rjups,
        s=0.072 * cst.Rjups,
        vmin=0.01 * cst.Rears,
        vmax=10.0 * cst.Rjups,
        size=n_sims,
    )
    period = get_truncnorm(l=40.2, s=1.0, vmin=0.0, vmax=1.0e5, size=n_sims)
    ecc = get_truncnorm(
        l=0.10, s=0.09, vmin=0.0, vmax=1.0, size=n_sims
    )
    if sort_ecc[body_id]:
        ecc = np.sort(ecc)
    argp = np.random.random(n_sims) * 360.0
    meanA = (360.0 * (0.0 - np.random.normal(7096.6226, 0.0031, size=n_sims)) / period) % 360.0
    inc = 89.0 + np.random.random(n_sims) * ((90.0-89.0)*2)
    longN = np.zeros((n_sims)) + 180.0
    params = np.column_stack((mass, radius, period, ecc, argp, meanA, inc, longN))
    grid_values[body_id] = params

    return grid_values

def set_grid_values_v1298Tau_02(n_sims=1):

    # DICT WITH MULTIPLE-PLANET SYSTEMS
    sort_ecc = [False, False, False, False, False]
    grid_values = {}

    # STAR
    body_id = 0  # 0 == star, 1 ==  planet b, 2 == perturber x
    mass = get_truncnorm(
        l=1.17, s=0.06, vmin=1.0e-2, vmax=2.0, size=n_sims
    )  # bounded to 0.01-2 Msun
    radius = get_truncnorm(
        l=1.28, s=0.07, vmin=1.0e-2, vmax=2.0, size=n_sims
    )  # bounded to 0.01-2 Rsun
    period = np.zeros((n_sims))
    ecc = np.zeros((n_sims))
    if sort_ecc[body_id]:
        ecc = np.sort(ecc)
    argp = np.zeros((n_sims))
    meanA = np.zeros((n_sims))
    inc = np.zeros((n_sims))
    longN = np.zeros((n_sims))
    params = np.column_stack((mass, radius, period, ecc, argp, meanA, inc, longN))
    grid_values[body_id] = params

    # PLANET b
    body_id = 1  # 0 == star, 1 ==  planet b, 2 == perturber c, ...
    mass = get_truncnorm(
        l=0.64 * cst.Mjups,
        s=0.19 * cst.Mjups,
        vmin=0.01 * cst.Mears,
        vmax=10.0 * cst.Mjups,
        size=n_sims,
    )  # bounded to 0.01 Mearth and 10 Mjups
    radius = get_truncnorm(
        l=0.868 * cst.Rjups,
        s=0.056 * cst.Rjups,
        vmin=0.01 * cst.Rears,
        vmax=10.0 * cst.Rjups,
        size=n_sims,
    )
    period = get_truncnorm(l=24.1399, s=0.0015, vmin=0.0, vmax=1.0e5, size=n_sims)
    ecc = get_truncnorm(
        l=0.08, s=0.02, vmin=0.0, vmax=1.0, size=n_sims
    ) # van eylen 2019
    if sort_ecc[body_id]:
        ecc = np.sort(ecc)
    argp = np.random.random(n_sims) * 360.0
    meanA = (360.0 * (0.0 - np.random.normal(7067.0486, 0.0015, size=n_sims)) / period) % 360.0
    inc = 88.7 + np.random.random(n_sims) * ((90.0-88.7)*2)
    longN = np.zeros((n_sims)) + 180.0
    params = np.column_stack((mass, radius, period, ecc, argp, meanA, inc, longN))
    grid_values[body_id] = params

    # PLANET c
    body_id = 2 
    mass = 0.01*cst.Mears + np.random.random(n_sims) * (0.24*cst.Mjups-0.01*cst.Mears)
    radius = get_truncnorm(
        l=0.460 * cst.Rjups,
        s=0.034 * cst.Rjups,
        vmin=0.01 * cst.Rears,
        vmax=10.0 * cst.Rjups,
        size=n_sims,
    )
    period = get_truncnorm(l=8.24892, s=0.00083, vmin=0.0, vmax=1.0e5, size=n_sims)
    ecc = get_truncnorm(
        l=0.08, s=0.02, vmin=0.0, vmax=0.3, size=n_sims
    ) # van eylen 2019
    if sort_ecc[body_id]:
        ecc = np.sort(ecc)
    argp = np.random.random(n_sims) * 360.0
    meanA = (360.0 * (0.0 - np.random.normal(7064.2801, 0.0041, size=n_sims)) / period) % 360.0
    inc = 87.5 + np.random.random(n_sims) * ((90.0-87.5)*2)
    longN = np.zeros((n_sims)) + 180.0
    params = np.column_stack((mass, radius, period, ecc, argp, meanA, inc, longN))
    grid_values[body_id] = params

    # PLANET d
    body_id = 3
    mass = 0.01*cst.Mears + np.random.random(n_sims) * (0.31*cst.Mjups-0.01*cst.Mears)
    radius = get_truncnorm(
        l=0.574 * cst.Rjups,
        s=0.041 * cst.Rjups,
        vmin=0.01 * cst.Rears,
        vmax=10.0 * cst.Rjups,
        size=n_sims,
    )
    period = get_truncnorm(l=12.4058, s=0.0018, vmin=0.0, vmax=1.0e5, size=n_sims)
    ecc = get_truncnorm(
        l=0.08, s=0.02, vmin=0.0, vmax=0.2, size=n_sims
    ) # van eylen 2019
    if sort_ecc[body_id]:
        ecc = np.sort(ecc)
    argp = np.random.random(n_sims) * 360.0
    meanA = (360.0 * (0.0 - np.random.normal(7072.3907, 0.0063, size=n_sims)) / period) % 360.0
    inc = 88.3 + np.random.random(n_sims) * ((90.0-88.3)*2)
    longN = np.zeros((n_sims)) + 180.0
    params = np.column_stack((mass, radius, period, ecc, argp, meanA, inc, longN))
    grid_values[body_id] = params

    # PLANET e
    body_id = 4
    mass = get_truncnorm(
        l=1.16 * cst.Mjups,
        s=0.30 * cst.Mjups,
        vmin=0.01 * cst.Mears,
        vmax=10.0 * cst.Mjups,
        size=n_sims,
    )  # bounded to 0.01 Mearth and 10 Mjups
    radius = get_truncnorm(
        l=0.735 * cst.Rjups,
        s=0.072 * cst.Rjups,
        vmin=0.01 * cst.Rears,
        vmax=10.0 * cst.Rjups,
        size=n_sims,
    )
    period = get_truncnorm(l=40.2, s=1.0, vmin=0.0, vmax=1.0e5, size=n_sims)
    ecc = get_truncnorm(
        l=0.10, s=0.09, vmin=0.0, vmax=1.0, size=n_sims
    )
    # ecc = get_truncnorm(
    #     l=0.08, s=0.02, vmin=0.0, vmax=1.0, size=n_sims
    # ) # van eylen 2019
    if sort_ecc[body_id]:
        ecc = np.sort(ecc)
    argp = np.random.random(n_sims) * 360.0
    meanA = (360.0 * (0.0 - np.random.normal(7096.6226, 0.0031, size=n_sims)) / period) % 360.0
    inc = 89.0 + np.random.random(n_sims) * ((90.0-89.0)*2)
    longN = np.zeros((n_sims)) + 180.0
    params = np.column_stack((mass, radius, period, ecc, argp, meanA, inc, longN))
    grid_values[body_id] = params

    return grid_values
