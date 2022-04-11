#!/usr/bin/env python
# -*- coding: utf-8 -*-

# no more "zero" integer division bugs!:P
#import argparse
import os
from hypothesis import example
import numpy as np  # array
from pytrades_lib import pytrades
# from constants import Mjups
import constants as cst
import ancillary as anc


def main():
    # integration parameters
    t_start = np.float64(2454965.0)
    t_epoch = np.float64(2455088.212)
    t_int = np.float64(500.0)
    step_in = np.float64(1.0e-3)

    print('set integration parameters')

    # number of bodies in the system
    n_body = 3

    print('set n_body')

    # example_folder = os.path.join(
    #     os.path.abspath('..'), 'trades_example', 'Kepler-9_example'
    # )

    # read rv data
    # rv_path = os.path.join(
    #     example_folder, 'obsRV.dat'
    # )
    # rv_data = np.genfromtxt(rv_path, usecols=(0, 1, 2)
    #                         )  # t_rv, rv_obs, erv_obs
    # n_rv = rv_data.shape[0]

    # t_rv = rv_data[:, 0]
    # rv_obs = rv_data[:, 1]
    # erv_obs = rv_data[:, 2]
    t_rv = np.array([
        2455342.947836,
        2455344.061419,
        2455351.037415,
        2455352.011289,
        2455367.083543,
        2455376.930532,
    ])
    rv_obs = np.array([
         17.15,
          2.23,
         -7.89,
         -4.61,
        -25.45,
         17.39,
    ])
    erv_obs = np.array([
        2.51,
        2.74,
        2.36,
        2.21,
        2.49,
        2.61,
    ])
    rv_data = np.column_stack((t_rv, rv_obs, erv_obs))
    n_rv = len(t_rv)

    print('set rv data')

    # read t0 data
    # t0b_path = os.path.join(
    #     example_folder, "NB2_observations.dat"
    # )
    # t0b_data = np.genfromtxt(t0b_path)
    # n_t0b = t0b_data.shape[0]
    epo_b = np.array([
        -5,  
        -4,  
        -2,  
        -1,  
        0,  
        2,  
        3,  
        4,  
        5, 
    ]) 
    t0_b = np.array([
        2454977.24875,
        2454996.48240,
        2455034.95437,
        2455054.19058,
        2455073.43412,
        2455111.92589,
        2455131.17167,
        2455150.42951,
        2455169.68103,
    ])
    et0_b = np.array([
        0.00059,
        0.00062,
        0.00052,
        0.00052,
        0.00072,
        0.00050,
        0.00048,
        0.00048,
        0.00046,
    ])
    t0b_data = np.column_stack((epo_b, t0_b, et0_b))
    n_t0b = len(t0_b)


    # t0c_path = os.path.join(
    #     example_folder, "NB3_observations.dat"
    # )
    # t0c_data = np.genfromtxt(t0c_path)
    # n_t0c = t0c_data.shape[0]
    epo_c = np.array([
        -3,
        -2,
        -1,
        0,
        1,
        2,
    ]) 
    t0_c = np.array([
        2454969.30577,
        2455008.33086,
        2455047.33560,
        2455086.31251,
        2455125.26284,
        2455164.18168,
    ])
    et0_c = np.array([
        0.00070,
        0.00061,
        0.00058,
        0.00059,
        0.00053,
        0.00055,
    ])
    t0c_data = np.column_stack((epo_c, t0_c, et0_c))
    n_t0c = len(t0_c)

    n_t0 = np.zeros((n_body)).astype(int)
    n_t0[1] = n_t0b
    n_t0[2] = n_t0c

    n_max_t0 = np.max(n_t0)

    t0_num = np.zeros((n_max_t0, n_body)).astype(int)
    t0_obs = np.zeros((n_max_t0, n_body))
    et0_obs = np.zeros((n_max_t0, n_body))

    t0_num[0:n_t0[1], 1] = t0b_data[:, 0].astype(int)
    t0_obs[0:n_t0[1], 1] = t0b_data[:, 1]
    et0_obs[0:n_t0[1], 1] = t0b_data[:, 2]

    t0_num[0:n_t0[2], 2] = t0c_data[:, 0].astype(int)
    t0_obs[0:n_t0[2], 2] = t0c_data[:, 1]
    et0_obs[0:n_t0[2], 2] = t0c_data[:, 2]

    print('set T0 data')
    print('n_t0 = ', n_t0)
    print('n_max_t0 = ', n_max_t0)

    pytrades.args_init(t_start, t_epoch, t_int, n_body,
                       n_t0, t0_num, t0_obs, et0_obs
                       )

    print('initialised args for trades')

    transit_flag = np.array([False, True, True])
    dur_check = 1
    # dur_check:
    # 1 yes computes the Total Transit Duration as T4 - T1
    # 0 no, but it returns a vector with same size of t0_sim, but with zero values

    M_msun   = np.array([1.0, 0.136985*cst.Mjups, 0.094348*cst.Mjups])
    R_rsun   = np.array([1.1, 0.842*cst.Rjups   , 0.823*cst.Rjups])
    P_day    = np.array([0.0, 19.238756         , 38.986097])
    ecc      = np.array([0.0, 0.058             , 0.068])
    argp_deg = np.array([0.0, 356.1             , 167.6])
    mA_deg   = np.array([0.0, 3.8               , 307.4])
    inc_deg  = np.array([0.0, 88.55             , 89.12])
    lN_deg   = np.array([0.0,  0.0              ,  0.0])

    print('set orbital elements')
    print('running trades:')

    print('- once ...')
    # old
    # rv_sim, t0_sim = pytrades.kelements_to_data(t_start,t_epoch,step_in,t_int,
    # M_msun,R_rsun,P_day,ecc,argp_deg,mA_deg,inc_deg,lN_deg,
    # t_rv,transit_flag,n_t0,t0_num)
    rv_sim, t0_sim, t14_sim, kel_sim = pytrades.kelements_to_data(
        t_start, t_epoch, step_in, t_int,
        M_msun, R_rsun,
        P_day, ecc, argp_deg, mA_deg, inc_deg, lN_deg,
        t_rv,
        transit_flag, dur_check, n_max_t0
    )

    print('- twice ...')
    # run it twice as check
    rv_sim, t0_sim, t14_sim, kel_sim = pytrades.kelements_to_data(
        t_start, t_epoch, step_in, t_int,
        M_msun, R_rsun,
        P_day, ecc, argp_deg, mA_deg, inc_deg, lN_deg,
        t_rv,
        transit_flag, dur_check, n_max_t0
    )

    print('end trades runs.')
    print()
    print('trades output summary:')
    print()
    print('## RV')
    print('# time_rv rv_obs rv_sim')
    for irv in range(0, n_rv):
        print('%18.6f %12.6f %12.6f' % (t_rv[irv], rv_obs[irv], rv_sim[irv]))
    print()

    print('## T0')
    for inb in range(1, n_body):
        print('# BODY %d (planet %d)' % (inb + 1, inb))
        print('# {:>6s} {:>20s} {:>20s} {:>16s} {:>16s} {:>16s} {:>16s} {:>16s} {:>16s} {:>16s} {:>16s} {:>16s}'.format(
            "t0_num", "t0_obs", "t0_sim", "t14_sim", "period_d", "sma_au", "ecc", "inc_deg", "meana_deg", "argp_deg", "truea_deg", "longn_deg"
        ))
        for itt in range(0, n_t0[inb]):
            kel_str = ' '.join(["{:16.6f}".format(kel) for kel in kel_sim[itt, inb, :]])
            print('{:8d} {:20.6f} {:20.6f} {:16.6f} {:s}'.format(
                t0_num[itt, inb], t0_obs[itt, inb],
                t0_sim[itt, inb], t14_sim[itt, inb], kel_str
            ))
        print()
    return


if __name__ == "__main__":
    main()
