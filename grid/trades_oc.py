#!/usr/bin/env python
# -*- coding: utf-8 -*-

import argparse
import numpy as np  # array
import h5py

# import astropy.time as atime
import os
import sys

import matplotlib as mpl

mpl.use("Agg")
import matplotlib.pyplot as plt

# matplotlib rc params
plt.rcParams['text.usetex']       = False
# plt.rcParams['font.family']       = 'sans-serif'
plt.rcParams['font.family']       = 'serif'
plt.rcParams['font.serif']        = ['Computer Modern Roman', 'Palatino', 'DejaVu Serif']
plt.rcParams['mathtext.fontset']  = 'cm'
plt.rcParams['figure.figsize']    = [5, 5]
plt.rcParams["figure.facecolor"]  = 'white'
plt.rcParams["savefig.facecolor"] = 'white'
plt.rcParams["figure.dpi"]        = 200
plt.rcParams["savefig.dpi"]       = 300
plt.rcParams["font.size"]         = 12
plt.rcParams["xtick.labelsize"]   = 10
plt.rcParams["ytick.labelsize"]   = plt.rcParams["xtick.labelsize"]


# ==============================================================================
# custom modules
script_path = os.path.realpath(__file__)
module_path = os.path.abspath(os.path.join(os.path.dirname(script_path), "../pytrades"))
sys.path.append(module_path)

import ancillary as anc
import constants as cst
import gls
from pytrades_lib import f90trades

from custom_grid_values import *

# =============================================================================
class Configuration:
    def __init__(self, sim_folder, grid_values, sub_folder="output", int_time=365.25, n_sims=1, n_cpu=1, seed=42):
        self.sim_folder  = sim_folder
        self.sub_folder  = sub_folder
        self.int_time    = int_time
        self.n_sims      = n_sims
        self.grid_values = grid_values
        self.n_cpu       = n_cpu
        self.seed        = seed

        return

    def update_with_system_info(self, system_id, system_name):
        self.system_id = system_id
        self.system_name = system_name

        return

# ==============================================================================

class exoSystem:
    def __init__(self,
        system_id,
        system_name,
        n_bodies,
        output_folder = "./output",
        seed=42,
        id_sim=1
        ):
        self.system_id = system_id  # numeric id of the system, starting from 0
        self.system_name = system_name  # name of the system
        self.id_sim = id_sim  # numeric id of the simulation
        self.seed = seed
        self.output_folder = os.path.abspath(output_folder)
        self.n_bodies = n_bodies  # number of the bodies in the system. Not planets!
        self.n_planets = n_bodies - 1  # number of planets
        self.mass = np.ones((n_bodies))
        self.radius = np.ones((n_bodies))
        self.period = np.zeros((n_bodies))
        self.period[1:] = 365.25
        self.sma = np.zeros((n_bodies))
        self.sma[1:] = 1.0
        self.ecc = np.zeros((n_bodies))
        self.argp = np.zeros((n_bodies))
        self.argp[1:] = 90.0
        self.meanA = np.zeros((n_bodies))
        self.inc = np.zeros((n_bodies))
        self.inc[1:] = 90.0
        self.longN = np.zeros((n_bodies))
        self.longN[1:] = 180.0
        self.stable = True

    # data_in structure: example 2-planet system
    # 0         1         2           3             4               5               6             7
    # planet_id mass_star radius_star mass_planet_b radius_planet_b period_planet_b tref_planet_b inclination_planet_b
    # planet_id mass_star radius_star mass_planet_c radius_planet_c period_planet_c tref_planet_c inclination_planet_c
    def update_system(self, body_id, mass, radius, period, ecc, argp, meanA, inc, longN):
        self.mass[body_id]   = mass * cst.Mears
        self.radius[body_id] = radius * cst.Rears
        self.period[body_id] = period
        self.sma[body_id]    = anc.period_to_semimajoraxis(
            self.mass[0],
            self.mass[body_id],
            self.period[body_id]
        )
        self.ecc[body_id]   = ecc
        self.argp[body_id]  = argp
        self.meanA[body_id] = meanA
        self.inc[body_id]   = inc
        self.longN[body_id] = longN

    def set_transits(self, transits):
        self.transits = transits

    def run_simulation(self, int_time):

        # n_planets = n_bodies - 1
        # system = exoSystem(id_system, id_sim, n_bodies)
        # system.update_system(parameters_in)
        # system.meanA[2:] = np.random.random(size = n_planets-1)*360.

        n_bodies = self.n_bodies
        nTT_max, nRV_max = f90trades.get_max_nt0_nrv(self.period, n_bodies)

        # TTs_full, id_TT_full, stats_TTs, time_rv_nmax, rv_nmax, stats_rv = f90trades.wrapper_run_grid_combination(system.mass,
        TTs_full, id_TT_full, stats_TTs, _, _, _ = f90trades.wrapper_run_grid_combination(
            self.mass,
            self.radius,
            self.period,
            self.sma,
            self.ecc,
            self.argp,
            self.meanA,
            self.inc,
            self.longN,
            nTT_max,
            nRV_max,
        )

        transits = [
            set_complete_transits(int_time, ibd, TTs_full, id_TT_full, stats_TTs)
            for ibd in range(2, n_bodies + 1)
        ]
        if np.any(stats_TTs):
            idx_max_TT = np.argmax(TTs_full)
            max_TT = TTs_full[idx_max_TT]
            if max_TT + self.period[id_TT_full[idx_max_TT] - 1] < int_time:
                self.stable = False
        else:
            self.stable = False

        self.set_transits(transits)

    def save_plot_oc(self, int_time):

        # creates output folders
        if(not os.path.isdir(self.output_folder)):
            os.makedirs(self.output_folder, exist_ok=True)

        # save self data into hdf5 file
        output_file = os.path.join(
            self.output_folder,
            '{:04d}_{:s}_{:04d}.hdf5'.format(
                self.system_id,
                self.system_name,
                self.id_sim
            )
        )
        if(os.path.exists(output_file) and os.path.isfile(output_file)):
            os.remove(output_file)

        h5f = h5py.File(output_file, 'w')
        gg = h5f.create_group('{:04d}_{:s}_{:04d}'.format(
            self.system_id,
            self.system_name,
            self.id_sim
            )
        )
        gg.attrs["system_name"] = str(self.system_name)
        gg.attrs["system_id"] = self.system_id
        gg.attrs["id_sim"] = self.id_sim
        gg.attrs['n_bodies'] = self.n_bodies
        gg.attrs['n_planets'] = self.n_planets
        gg.attrs['stable'] = self.stable
        gg.attrs['seed'] = self.seed

        gg.create_dataset('mass', data=self.mass, dtype=np.float64)
        gg['mass'].attrs['unit'] = 'Msun'
        gg.create_dataset('radius', data=self.radius, dtype=np.float64)
        gg['radius'].attrs['unit'] = 'Rsun'
        gg.create_dataset('period', data=self.period, dtype=np.float64)
        gg['period'].attrs['unit'] = 'days'
        gg.create_dataset('sma', data=self.sma, dtype=np.float64)
        gg['sma'].attrs['unit'] = 'au'
        gg.create_dataset('ecc', data=self.ecc, dtype=np.float64)
        gg.create_dataset('argp', data=self.argp, dtype=np.float64)
        gg['argp'].attrs['unit'] = 'deg'
        gg.create_dataset('meanA', data=self.meanA, dtype=np.float64)
        gg['meanA'].attrs['unit'] = 'deg'
        gg.create_dataset('inc', data=self.inc, dtype=np.float64)
        gg['inc'].attrs['unit'] = 'deg'
        gg.create_dataset('longN', data=self.longN, dtype=np.float64)
        gg['longN'].attrs['unit'] = 'deg'

        # save planet TTs (and stuff) and create plot
        output_fig = os.path.splitext(output_file)[0]
        if(self.stable):
            sys_stat = 'stable'
        else:
            sys_stat = 'unstable'


        deltax = 0.03*int_time
        xminl = 0 - deltax
        xmaxl = int_time + deltax

        # fig = plt.figure(figsize=(12,12))
        fig = plt.figure(figsize=(10,10))
        fig.subplots_adjust(hspace=0.3)
        fig.suptitle(r'{:s} ({:04d}): sim {:04d} ({})'.format(
            self.system_name,
            self.system_id,
            self.id_sim,
            sys_stat,
            str(self.stable)
            )
        )

        normal_font = plt.rcParams["font.size"]
        small_font = plt.rcParams["xtick.labelsize"]
        axs = []

        id_name = anc.letters
        for ipl in range(0,self.n_planets):
            planet_str = id_name[ipl+1]
            oo = self.transits[ipl]
            og = gg.create_group('planet_{:s}'.format(planet_str))
            og.attrs['planet_id'] = oo.planet_id
            og.create_dataset('TTs', data=oo.TTs, dtype=np.float64)
            og['TTs'].attrs['nTTs'] = oo.nTTs
            og.create_dataset('epo', data=oo.epo, dtype=np.float64)
            og.create_dataset('linTTs', data=oo.linTTs, dtype=np.float64)
            og['linTTs'].attrs['TTref'] = oo.TTref
            og['linTTs'].attrs['Pref'] = oo.Pref
            og.create_dataset('OCs', data=oo.OCs, dtype=np.float64)
            og['OCs'].attrs['unit'] = 'days'
            og['OCs'].attrs['P_TTV_d'] = oo.P_TTV_d
            og['OCs'].attrs['amp_TTV_d'] = oo.amp_TTV_d
            og['OCs'].attrs['amp_TTV_m'] = oo.amp_TTV_m

            if oo.amp_TTV_d < cst.min2day:
                u_s, u_v = "sec", cst.day2sec
            elif cst.min2day <= oo.amp_TTV_d < cst.hour2day:
                u_s, u_v = "min", cst.day2min
            elif cst.hour2day <= oo.amp_TTV_d < 1.0:
                u_s, u_v = "hour", cst.day2hour
            else:
                u_s, u_v = "day", 1.0


            ax = plt.subplot2grid((self.n_planets, 1), (ipl,0))
            tlin_str = r'planet {:s} ($P={:.6f}$): $TT_\mathrm{{lin}} = {:.6f} + N \times {:.6f}$'.format(
                planet_str,
                self.period[ipl+1],
                oo.TTref,
                oo.Pref
            )
            ttv_str = r'$P_\mathrm{{TTV}}={:.6f}\ \mathrm{{day}}\ \mathrm{{Amp}}_\mathrm{{TTV}}={:.6f}\ \mathrm{{{:s}}}$'.format(
                oo.P_TTV_d,
                oo.amp_TTV_d*u_v, u_s
            )
            # ax.set_title(r'%s , %s' %(tlin_str, ttv_str), fontsize=9)
            ax.text(0.5, 1.11, tlin_str,
                fontsize=normal_font, #small_font,
                horizontalalignment='center',
                verticalalignment='bottom',
                transform=ax.transAxes
                )
            ax.text(0.5, 1.01, ttv_str,
                fontsize=normal_font, #small_font,
                horizontalalignment='center',
                verticalalignment='bottom',
                transform=ax.transAxes
                )

            ax.axhline(0., color='black')
            if(oo.nTTs > 0):
                xTT = oo.TTs - oo.TTs[0]
                ax.errorbar(
                    xTT,
                    oo.OCs*u_v,
                    color='gray',
                    marker='o', ms=4, mfc='C0', mec='None',
                    ls='-',
                    ecolor='lightgray', capsize=0
                )
            else:
                ax.text(0.5, 0.5, "Stable: {}".format(str(self.stable)),
                    fontsize=normal_font, #small_font,
                    horizontalalignment='center',
                    verticalalignment='center',
                    transform=ax.transAxes
                    )

            ax.set_xlim([xminl, xmaxl])
            ax.set_xlabel('Time (days)')
            ax.set_ylabel('O-C ({:s})'.format(u_s))
            axs.append(ax)

        h5f.close()

        fig.align_ylabels(axs)
        # rect = [0, 0.03, 1, 0.95] # ok
        rect = [0, 0.01, 1, 0.97]
        plt.tight_layout(rect=rect) # rect = [left, bottom, right, top]

        fig.savefig('{:s}.png'.format(output_fig), bbox_inches='tight')
        plt.close(fig)

        print('WRITTEN HDF5 FILE: {:s}'.format(output_file))
        print('WRITTEN IMAGE: {:s}.png'.format(output_fig))

        return

# ==============================================================================
class transit_times:
    def __init__(self, planet_id, TTs_full, id_TT_full, stats_TT):
        self.planet_id = planet_id
        sel_TTs = np.logical_and(
            np.logical_and(id_TT_full == planet_id, stats_TT), TTs_full > -9.0e10
        )
        self.TTs = np.sort(TTs_full[sel_TTs])
        # self.nTTs = np.shape(self.TTs)[0]
        self.nTTs = np.sum(sel_TTs)
        self.epo = 0
        self.TTref = 0.0
        self.Pref = 1.0
        self.TP_err = []
        self.linTTs = []
        self.OCs = np.zeros((self.nTTs))
        self.P_TTV_d = 0.0
        self.amp_TTV_d = 0.0
        self.amp_TTV_m = 0.0

    def set_linearephem(self):
        self.epo, self.TTref, self.Pref, self.TP_err = anc.compute_lin_ephem(self.TTs)

        # self.epo, self.TTref, self.Pref, self.TP_err = anc.compute_lin_ephem(self.TTs,
        #     modefit='minimize'
        # )
        
        # self.epo, self.TTref, self.Pref, self.TP_err = anc.compute_lin_ephem(self.TTs,
        #     modefit='curve_fit'
        # )
        
        # self.epo, self.TTref, self.Pref, self.TP_err = anc.compute_lin_ephem(self.TTs,
        #     modefit='sklearn'
        # ) # BAD

        # self.TTref = self.TTs[self.nTTs//2]
        # dT = np.diff(self.TTs-self.TTref)
        # self.Pref = dT[(self.nTTs-1)//2]

        # self.TTref = self.TTs[0]
        # self.Pref = np.median(np.diff(self.TTs))
        # self.epo = np.rint((self.TTs-self.TTref)/self.Pref)
        # self.TP_err = [0.0, 0.0]

        # self.TTref = self.TTs[0]
        # dT = np.diff(self.TTs-self.TTref)
        # self.Pref = dT[0]
        # self.epo = np.arange(0,self.nTTs)
        # self.epo, self.TTref, self.Pref, self.TP_err = anc.compute_lin_ephem(self.TTs,
        #     epoin=self.epo,
        #     modefit='wls'
        # )


    def set_oc(self, int_time):
        self.linTTs = self.TTref + self.epo * self.Pref
        self.OCs = self.TTs - self.linTTs

        Pstart = 1.5 * self.Pref
        Pfinish = 1.5 * int_time # 0.75 * int_time
        # print("Pstart = ", Pstart, " Pfinish = ", Pfinish, "( Pref = ", self.Pref, ")")
        fstart = 1.0 / Pfinish
        fend = 1.0 / Pstart
        # print("fstar = ", fstart, " ffinish = ", fend)
        freq = np.linspace(start=fstart, stop=fend, num=1000, endpoint=True)
        gls_o = gls.Gls(
            (self.TTs, self.OCs, np.ones(self.nTTs) / 86400.0),
            fbeg=fstart,
            fend=fend,
            freq=freq,
            ofac=1,
            verbose=False,
        )

        self.P_TTV_d = 1.0 / gls_o.hpstat["fbest"]
        self.amp_TTV_d = gls_o.hpstat["amp"]
        self.amp_TTV_m = self.amp_TTV_d * 1440.0
        del gls_o


# ==============================================================================


def set_complete_transits(int_time, planet_id, TTs_full, id_TT_full, stats_TT):

    TT_obj = transit_times(planet_id, TTs_full, id_TT_full, stats_TT)
    print("{:s} nTTs = {:d}".format(anc.letters[planet_id - 1], TT_obj.nTTs))
    if TT_obj.nTTs > 2:
        TT_obj.set_linearephem()
        if TT_obj.Pref > 0.0:
            TT_obj.set_oc(int_time)

    return TT_obj

# =============================================================================

# def run_multiple_sims(grid_values, sim_folder, sub_folder="output", int_time=365.25, n_sims=1):
def run_multiple_sims(conf):
    print()
    grid_values = conf.grid_values
    sim_folder  = conf.sim_folder
    sub_folder  = conf.sub_folder
    int_time    = conf.int_time
    n_cpu       = conf.n_cpu
    n_sims      = conf.n_sims


    output_folder = os.path.join(sim_folder, sub_folder)
    os.makedirs(output_folder, exist_ok=True)
    print("OUTPUT FOLDER: {}".format(output_folder))

    # saving input grid_values:
    for body_id, params in grid_values.items():
        file_out = os.path.join(output_folder, "body_{:02d}_grid_values.dat".format(body_id))
        head="mass_Msun radius_Rsun period_d ecc argp_deg meanA_deg inc_deg longN_deg"
        fmt="%23.16e"
        np.savetxt(file_out, params, fmt=fmt, header=head)

    # init trades with proper base folder
    print('SET TRADES')
    f90trades.initialize_trades(sim_folder, '', n_cpu)

    n_bodies = f90trades.n_bodies
    # int_time = f90trades.tint

    f90trades.tint = int_time
    n_bodies = f90trades.n_bodies

    # save simulation info into file:
    summary_file = os.path.join(
        output_folder,
        "summary.dat"
    )
    with open(summary_file, 'w') as of:
        of.write("# index_sim body_id Tref Pref_d P_TTV_d A_TTV_d Stable_flag\n")

        for i_sim in range(0, n_sims):
            exosys = exoSystem(
                conf.system_id,
                conf.system_name,
                n_bodies,
                output_folder=output_folder,
                id_sim=i_sim
            )
            print("\n{:s} (id {:04d}) simulation n. {:d}".format(
                    exosys.system_name,
                    exosys.system_id,
                    exosys.id_sim
                )
            )
            for body_id, params in grid_values.items():
                print("body = {}".format(body_id))
                mass, radius, period, ecc, argp, meanA, inc, longN = params[i_sim, :]
                print(mass, radius, period, ecc, argp, meanA, inc, longN)
                exosys.update_system(body_id, mass, radius, period, ecc, argp, meanA, inc, longN)
            exosys.run_simulation(int_time)
            print("The system {} - sim {} is stable? {}".format(
                exosys.system_name,
                i_sim,
                str(exosys.stable)
                )
            )
            exosys.save_plot_oc(int_time)
            # save each line 
            for i_body in range(1, n_bodies):
                tra = exosys.transits[i_body-1]
                line = "{:05d} {:03d} {:20.6f} {:13.6f} {:13.6f} {:13.6f} {:d}\n".format(
                    i_sim, 
                    i_body,
                    tra.TTref,
                    tra.Pref,
                    tra.P_TTV_d,
                    tra.amp_TTV_d,
                    exosys.stable
                )
                of.write(line)

    return


# =============================================================================



# =============================================================================
# =============================================================================
def main():

    seed = 42
    np.random.seed(seed)
    n_sims = 100
    
    # grid_values = set_grid_values_01(n_sims=n_sims)
    # grid_values = set_grid_values_02(n_sims=n_sims)
    # grid_values = set_grid_values_03(n_sims=n_sims)

    grid_values = set_grid_values_v1298Tau_01(n_sims=n_sims)
    # grid_values = set_grid_values_v1298Tau_02(n_sims=n_sims)

    # SIM FOLDER FOR TRADES INITIALISATION BASED ON THE SYSTEM AND OUTPUT
    # WASP-8
    # sim_folder = os.path.join(
    #     os.path.abspath("/home/borsato/Dropbox/Research/exoplanets/objects/WASP/WASP-8/sims/2021-10-22_WASP-8_grid_OC"),
    #     ""
    #     )
    # sub_folder = "2021-10-25_output_004"
    # v1298Tau
    sim_folder = os.path.join(
        os.path.abspath("/home/borsato/Dropbox/Research/exoplanets/objects/V1298TAu/sims/2021-11-03_V1298Tau_oc_001/"),
        ""
        )
    sub_folder = "2021-11-05_output_001"
    # sub_folder = "2021-11-05_output_002"

    #  INTEGRATION TIME OF X YEARS
    int_time = 365.25*3.5

    conf = Configuration(
        sim_folder,
        grid_values,
        sub_folder=sub_folder,
        int_time=int_time,
        n_sims=n_sims,
        n_cpu=1,
        seed=seed
    )
    # conf.update_with_system_info(0, "WASP-8_bd")
    conf.update_with_system_info(0, "V1298Tau")

    run_multiple_sims(conf)
    return


# =============================================================================
# =============================================================================

if __name__ == '__main__':
    main()
