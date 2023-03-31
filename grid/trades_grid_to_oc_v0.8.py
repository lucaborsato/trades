#!/usr/bin/env python
# -*- coding: utf-8 -*-

## TO DO LIST
# 0) simulate the grid, starting from the t_epoch to end of the CHEOPS mission
# A)
# 1) read grid_summary file with fitted/grid sim_id,Mc,Pc,ec,ic & fitness: 'summary_grid_sims_0.dat'
# 2) run simulation for each grid and computes TTs and RVs
# 3) compute O-C for time range (linear ephemeris computed from all TTs)
# 4) compute O-C parameters (sigma_TTV, P_TTV) and plot
# 5) RV plot - optional
# ==============================================================================
# IMPORT

 # no more "zero" integer division bugs!:P
import argparse
import numpy as np # array
import h5py
import astropy.time as atime
import os
import sys
#import emcee.interruptible_pool as emcee_pool
import gc

import gls

# import pp

import matplotlib as mpl
#from matplotlib import use as mpluse
mpl.use("Agg", warn=False)
#mpluse("Qt4Agg")
#mpluse("TkAgg")
import matplotlib.pyplot as plt
# matplotlib rc params
plt.rcParams['text.usetex'] = True
plt.rcParams['font.family'] = 'serif'
plt.rcParams['figure.figsize'] = [6, 6]
plt.rcParams["figure.facecolor"] = 'white'
plt.rcParams["savefig.facecolor"] = 'white'
plt.rcParams["figure.dpi"]  = 200
plt.rcParams["savefig.dpi"] = 300
plt.rcParams["font.size"]   = 12

import matplotlib.cm as cm
from matplotlib.ticker import ScalarFormatter

script_path = os.path.realpath(__file__)
module_path = os.path.abspath(os.path.join(os.path.dirname(script_path),
                                           '../pytrades'))
sys.path.append(module_path)

import ancillary as anc
import constants as cst
from pytrades_lib import f90trades

import multiprocessing as mpro

import _thread

import time

# ==============================================================================

# read cli arguments
def get_args():
  
  parser = argparse.ArgumentParser(description='CHEOPS TTV LEVEL ANALYSIS')
  parser.add_argument('-p', '--path',
    action='store', dest='full_path', required=True,
    help='The path (absolute or relative) with simulation files for TRADES.'
    )
  
  #parser.add_argument('-lm', '--lm-flag', action='store', dest='lm_flag', default='0', help='LM flag. If you want to read files w/o LM set to 0, otherwise 1. Default is 0')
  
  # SUBFOLDER TO SAVE ALL SIMULATION DATA
  parser.add_argument('-sf', '--sub-folder', '--sub-f',
    action='store', dest='sub_folder', default='grid2ttv',
    help='Sub-folder name, without full path. Default = grid2ttv'
    )
  
  parser.add_argument('-e', '--e-T0', '--err-T0',
    action='store', dest='err_T0', default='0.0',
    help='Mean uncertainty on the Transit Times (T0s).'+\
          'Provide value with unit "d" "m" "s" for days, minutes, and seconds in this way:'+\
          '"0.001528d" or "2.2m" or "132.s". Default unit is days and value set to 0.'
    )
  
  parser.add_argument('-mr', '--mr-type', '--mass-radius-type',
                      action='store', dest='mr_type', default='e',
                      help='Mass and radius type to convert the grid. '+\
                      '"e ea ear eart earth" or "j ju jup jupi jupit jupite jupiter" or '+\
                      '"n ne nep nept neptu neptun neptune".'
                      )

  parser.add_argument('-ts', '--ts', '-time-start','--time-start',  
                      action='store', dest='t_start', default=0.0, 
                      help='Enter starting time of selection for the plot in same unit of the simulation. Default is 0.0')

  parser.add_argument('-te', '--te', '-time-end','--time-end',  
                      action='store', dest='t_end', default=365.25, 
                      help='Enter ending time of selection for the plot in same unit of the simulation. Default is 365.25')

  parser.add_argument('-s', '--seed', 
    action='store', dest='seed', default='None',
    help='Seed for random number generator. Defauto is None.'
    )
  
  parser.add_argument('-c', '--cpu', '--threads', 
    action='store', dest='nthreads', default=1, type=int, 
    help='Number of threads/cpu to use. Default is 1.'
    )


  cli = parser.parse_args()
  cli.full_path = os.path.join(os.path.abspath(cli.full_path), '')
  cli.sub_folder = os.path.join(os.path.relpath(cli.sub_folder), '')

    
  try:
    cli.t_start = np.float64(cli.t_start)
  except:
    cli.t_start = 0.0

  try:
    cli.t_end = np.float64(cli.t_end)
  except:
    cli.t_end = 365.25

  try:
    cli.seed = int(cli.seed)
  except:
    cli.seed = None


  return cli

# -------------------------------------

def set_err_T0(err_T0_in):
  
  if('d' in err_T0_in):
    err_T0 = np.float64(err_T0_in.split('d')[0])
  elif('m' in err_T0_in):
    err_T0 = np.float64(err_T0_in.split('m')[0])/1440.
  elif('s' in err_T0_in):
    err_T0 = np.float64(err_T0_in.split('s')[0])/86400.
  else:
    anc.print_both('WARNING: wrong mean err_T0: set to 30s')
    err_T0 = 30.0/86400.0
    
  return err_T0

# -------------------------------------

def set_analysis_folder_log(full_path, sub_folder):
  
  out_folder = os.path.join(full_path, sub_folder)
  if(not os.path.isdir(out_folder)): os.makedirs(out_folder)
  
  anc.copy_simulation_files(full_path, out_folder)
  
  log_file = os.path.join(out_folder, 'analysis_log.txt')
  of_log = open(log_file, 'w')
  
  return out_folder, log_file, of_log

# -------------------------------------

def mass_factor(mr_type):
  
  if(mr_type.lower() in 'j ju jup jupi jupit jupite jupiter'.split()):
    #mr_convert = Msun2jup
    mr_convert = [cst.Msjup, cst.Rsjup]
    mr_unit = ['M_jup', 'R_jup']
  elif(mr_type.lower() in 'n ne nep nept neptu neptun neptune'.split()):
    #mr_convert = Msun2nep
    mr_convert = [cst.Msnep, cst.Rsnep]
    mr_unit = ['M_nep', 'R_nep']
  else:
    #mr_convert = Msun2ear
    mr_convert = [cst.Msear, cst.Rsear]
    mr_unit = ['M_ear', 'R_ear']
  
  return mr_convert, mr_unit
    
# -------------------------------------

def set_names_parameters(id_body, mr_unit):
  
  body_num = '{0:d}'.format(id_body)
  base_names_lst = 'M R P a e w mA tau i lN'.split()
  units = '{0:s} {1:s} day au - deg deg day deg deg'.format(mr_unit[0], mr_unit[1])
  units = units.split()
  n_names = len(base_names_lst)
  names_lst = [r'{0:s}{1:s}'.format(base_names_lst[i], body_num) for i in range(0, n_names)]
  temp_lst = [r'{0:s} [{1:s}]'.format(names_lst[i], units[i].strip()) for i in range(0, n_names)]
  names_str = ' '.join([r'{0:23s}'.format(temp_lst[i]) for i in range(0, n_names)])
  
  return names_lst, names_str, units

# -------------------------------------

def select_and_sort_rv(time_rv_nmax, rv_nmax, stats_rv):
  
  time_rv_tosort = time_rv_nmax[stats_rv].copy()
  rv_tosort      = rv_nmax[stats_rv].copy()
  id_s           = np.argsort(time_rv_tosort)
  time_rv        = time_rv_tosort[id_s]
  rv             = rv_tosort[id_s]
  
  return time_rv, rv


# ==============================================================================

# # -------------------------------------

# def create_hdf5_file(full_path, seed, id_body, perturber_grid, mr_type, t_sel_min, t_sel_max, err_T0, of_log=None):
  
#   mr_convert, mr_unit = mass_factor(mr_type)
#   names_lst, names_str, units = set_names_parameters(id_body, mr_unit, of_log)
#   ngrid, _ = np.shape(perturber_grid)
#   write_grid = perturber_grid.copy()
#   write_grid[:,0] = write_grid[:,0] * mr_convert[0] # mass
#   write_grid[:,1] = write_grid[:,1] * mr_convert[1] # radius
  
#   file_name = os.path.join(full_path, 'rms_analysis_summary.hdf5')
#   f_h5 = h5py.File(file_name, 'w')
  
#   f_h5.create_dataset('grid/grid_names',
#                       data=anc.encode_list(names_lst),
#                       dtype='S20', compression="gzip")
#   f_h5.create_dataset('grid/grid_units',
#                       data=anc.encode_list(units),
#                       dtype='S20', compression="gzip")
#   f_h5.create_dataset('grid/perturber_grid',
#                       data=write_grid,
#                       dtype=np.float64, compression="gzip")
  
#   f_h5['grid'].attrs['ngrid'] = ngrid
#   f_h5['grid'].attrs['seed'] = seed
#   f_h5['grid'].attrs['t_sel_min'] = t_sel_min
#   f_h5['grid'].attrs['t_sel_max'] = t_sel_max
#   f_h5['grid'].attrs['err_T0_day'] = err_T0
#   f_h5['grid'].attrs['err_T0_min'] = err_T0*1440.
  
#   return file_name, f_h5, mr_convert, mr_unit, names_lst, names_str, units


# # -------------------------------------

# def save_results_sim(f_h5, results_sim, mr_type):
  
#   mr_convert, _ = mass_factor(mr_type)
#   temp_par = results_sim['perturber_par'].copy()
#   temp_par[0] = temp_par[0]*mr_convert[0]
#   temp_par[1] = temp_par[1]*mr_convert[1]
  
#   sim_num = results_sim['sim_num']
#   f_h5.create_dataset('sim_{0:06d}/perturber_par_names'.format(sim_num),
#                       data=anc.encode_list(np.array(results_sim['perturber_par_names'],
#                                                     dtype='S10')),
#                       dtype='S20', compression="gzip")
#   f_h5.create_dataset('sim_{0:06d}/perturber_par_units'.format(sim_num),
#                       data=anc.encode_list(np.array(results_sim['perturber_par_units'],
#                                                     dtype='S10')),
#                       dtype='S20', compression="gzip")
#   f_h5.create_dataset('sim_{0:06d}/perturber_par'.format(sim_num),
#                       data=temp_par,
#                       dtype=np.float64, compression="gzip")
#   f_h5.create_dataset('sim_{0:06d}/MR_sun2unit'.format(sim_num),
#                       data=np.array(mr_convert),
#                       dtype=np.float64, compression="gzip")
  

#   f_h5.create_dataset('sim_{0:06d}/keplerian_elements_header'.format(sim_num),
#                       data=anc.encode_list(np.array(results_sim['keplerian_elements_header'], 
#                                                     dtype='S10')),
#                       dtype='S20', compression="gzip")
#   f_h5.create_dataset('sim_{0:06d}/keplerian_elements'.format(sim_num),
#                       data=results_sim['keplerian_elements'],
#                       dtype=np.float64, compression="gzip")
  
#   #f_h5.create_dataset('sim_%06d/ttra_all'%(sim_num), data=results_sim['ttra_all'], dtype=np.float64, compression="gzip")
#   #f_h5.create_dataset('sim_%06d/sel_pre_cheops'%(sim_num), data=results_sim['sel_pre_cheops'], dtype=bool, compression="gzip")
#   #f_h5.create_dataset('sim_%06d/sel_cheops'%(sim_num), data=results_sim['sel_cheops'], dtype=bool, compression="gzip")
#   #f_h5.create_dataset('sim_%06d/T0_eT0'%(sim_num), data=results_sim['T0_eT0'], dtype=np.float64, compression="gzip")
#   sel_cheops = results_sim['sel_cheops']
#   if(np.sum(sel_cheops.astype(int)) > 0):
#     f_h5.create_dataset('sim_{0:06d}/T0_eT0'.format(sim_num),
#                         data=results_sim['T0_eT0'][sel_cheops,:],
#                         dtype=np.float64, compression="gzip")
#   else:
#     f_h5.create_dataset('sim_{0:06d}/T0_eT0'.format(sim_num),
#                         data=np.zeros((2)),
#                         dtype=np.float64, compression="gzip")
  
#   #f_h5.create_dataset('sim_%06d/time_rv'%(sim_num), data=results_sim['time_rv'], dtype=np.float64, compression="gzip")
#   f_h5.create_dataset('sim_{0:06d}/rms_rv'.format(sim_num),
#                       data=np.array([results_sim['rms_rv']]),
#                       dtype=np.float64, compression="gzip")
  
#   f_h5.create_dataset('sim_{0:06d}/linear_ephem'.format(sim_num),
#                       data=results_sim['linear_ephem'],
#                       dtype=np.float64, compression="gzip")
  
#   f_h5.create_dataset('sim_{0:06d}/rms_ttv'.format(sim_num),
#                       data=np.array([results_sim['rms_ttv']]),
#                       dtype=np.float64, compression="gzip")
#   f_h5.create_dataset('sim_{0:06d}/nT0_sel'.format(sim_num),
#                       data=results_sim['nT0_sel'],
#                       dtype=np.int, compression="gzip")
#   f_h5.create_dataset('sim_{0:06d}/rms_mc'.format(sim_num),
#                       data=results_sim['rms_mc'], 
#                       dtype=np.float64, compression="gzip")

#   f_h5.create_dataset('sim_{0:06d}/chi2r_cheops'.format(sim_num), 
#                       data=np.array([results_sim['chi2r_cheops']]),
#                       dtype=np.float64, compression="gzip")
#   f_h5.create_dataset('sim_{0:06d}/chi2r_mc'.format(sim_num),
#                       data=results_sim['chi2r_mc'],
#                       dtype=np.float64, compression="gzip")

#   f_h5['sim_{0:06d}'.format(sim_num)].attrs['sim_num'] = results_sim['sim_num']
#   f_h5['sim_{0:06d}'.format(sim_num)].attrs['nMC'] = results_sim['nMC']
#   f_h5['sim_{0:06d}'.format(sim_num)].attrs['A_TTV'] = results_sim['A_TTV']
#   f_h5['sim_{0:06d}'.format(sim_num)].attrs['P_TTV'] = results_sim['P_TTV']
  
#   return f_h5

# # -------------------------------------

# def save_rms_analysis(f_h5, rms_grid, above_threshold, nT0_sel, rms_rv):
  
#   rms_header_str = 'rms_ttv_d {0:s}'.format(' '.join(
#                    ['rms_mean_N{0:d}_d rms_median_N{1:d}_d'.format(nT0_sel[i_n], nT0_sel[i_n])
#                     for i_n in range(len(nT0_sel))]
#                     ))
#   rms_header_lst = rms_header_str.split()
  
#   f_h5.create_dataset('grid/rms_grid',
#                       data=rms_grid,
#                       dtype=np.float64, compression="gzip")
  
#   f_h5.create_dataset('grid/rms_header', 
#                       data= anc.encode_list(np.array(rms_header_lst, dtype='S20')),
#                       dtype='S20', compression="gzip")
  
#   f_h5.create_dataset('grid/above_threshold', 
#                       data=above_threshold, 
#                       dtype=bool, compression="gzip")
  
#   f_h5.create_dataset('grid/rms_rv', 
#                       data=rms_rv, 
#                       dtype=np.float64, compression="gzip")
  
#   return f_h5

# # -------------------------------------

# ==============================================================================

# -------------------------------------

def simulation_and_ttv_analysis(i_grid, n_bodies, perturber_par, t_start, t_end, \
  err_T0, names_lst, units, mr_type, plot_folder, of_log=None):
  
  pid = os.getpid()
  tid = _thread.get_ident()
  results_sim = {}
  sim_num = i_grid + 1
  
  tsize=9
  tsize2=tsize-3
  lsize=10
  msize=4
  
  kel_header = 'Mass_Msun Radius_Rsun Period_d sma_au ecc argp_deg meanA_deg inc_deg longN_deg'.split()
  #anc.print_both('\n - sim %06d: START\n' %(sim_num), of_log)

  perturber_id = f90trades.idpert
  dof = 1

  try:
    #mass, radius, period, sma, ecc, argp, mean_an, inc, long_n, nt0_full, nrv_nmax = f90trades.wrapper_set_grid_orbelem(sim_num, perturber_id, perturber_grid, n_bodies)
    mass, radius, period, sma, ecc, argp, mean_an, inc, long_n, nt0_full, nrv_nmax = \
      f90trades.wrapper_set_grid_orbelem(sim_num, perturber_id, perturber_par, n_bodies)
    #print ' | Set grid keplerian elements'
    kep_elem = np.column_stack((mass, radius, period, sma, ecc, argp, mean_an, inc, long_n))
    print(' | PID {0:d} || TID {1:d} sim {2:d}: set grid par'
          .format(pid, tid, sim_num))
  except:
    print(' | PID {0:d} || TID {1:d} sim {2:d}: error setting grid par'
          .format(pid, tid, sim_num))
    kep_elem = np.zeros((n_bodies, 9))
    
  results_sim['sim_num']             = sim_num
  results_sim['perturber_par_names'] = names_lst
  results_sim['perturber_par_units'] = units
  results_sim['perturber_par']       = perturber_par
  
  results_sim['keplerian_elements_header'] = kel_header
  results_sim['keplerian_elements']        = kep_elem
  
  # RV semi-amplitude from Keplerian ...
  K_rv = {}
  for i in range(1,n_bodies):
    K_rv[i+1] = anc.compute_Kms(mass[0], mass[i]*cst.Msjup,
                                inc[i], period[i], ecc[i]
                                )
  results_sim['K_rv'] = K_rv

  # run trades simulations
  try:
    ttra_full, id_ttra_full, stats_ttra, time_rv_nmax, rv_nmax, stats_rv = \
      f90trades.wrapper_run_grid_combination(mass, radius, period, sma, ecc, argp, mean_an, inc, long_n,
                                            nt0_full, nrv_nmax
                                            )
    print(' | PID {0:d} || TID {1:d} sim {2:d}: completed integration'
          .format(pid, tid, sim_num))
    results_sim['ttra_full']    = ttra_full
    results_sim['id_ttra_full'] = id_ttra_full
    results_sim['stats_ttra']   = stats_ttra

    results_sim['time_rv_nmax'] = time_rv_nmax
    results_sim['rv_nmax']      = rv_nmax
    results_sim['stats_rv']     = stats_rv

  except:
    ttra_full = None
    print(' | PID {0:d} || TID {1:d} sim {2:d}: error integrating'
          .format(pid, tid, sim_num))
    results_sim['ttra_full']    = -np.ones(1)
    results_sim['id_ttra_full'] = -np.ones(1)
    results_sim['stats_ttra']   = -np.ones(1)

    results_sim['time_rv_nmax'] = -np.ones(1)
    results_sim['rv_nmax']      = -np.ones(1)
    results_sim['stats_rv']     = -np.ones(1)  


  A_TTV        = {}
  P_TTV        = {}
  linear_ephem = {}
  statistics   = {}
  oc_unit      = [1.0, "day"]
  
  nfreq = 1000
  if(ttra_full is not None):
    for i in range(1,n_bodies): # star == 0
      i_body  = i + 1 # first planet has id == 2
      # select transits of the planet
      sel_tra = np.logical_and(np.logical_and(id_ttra_full == i_body, stats_ttra),
                               ttra_full > -9.0e10
                               )
      TTs     = np.sort(ttra_full[sel_tra])
      # TTs     = ttra_full[sel_tra]
      nTTs    = np.shape(TTs)[0]

      plot_file = os.path.join(plot_folder, 'oc_{0:05d}_NB{1:02d}'.format(sim_num, i_body))
      if(nTTs > 3):
        print(' | PID {0:d} || TID {1:d} sim {2:d}: found {4:d} transits for body {3:d}.'
              .format(pid, tid, sim_num, i_body, nTTs))
        # computes O-C
        eTTs  = np.zeros(nTTs) + err_T0
        epo, Tref, Pref, _ = anc.compute_lin_ephem(TTs, eTTs)
        print(' | PID {0:d} || TID {1:d} sim {2:d}: Tref = {3:.5f} Pref = {4:.5f}.'
              .format(pid, tid, sim_num, Tref, Pref))
        TTlin = Tref + epo * Pref
        oc    = TTs - TTlin
        # computes GLS
        Pgls0 = 3.0*Pref
        Pgls1 = TTs[-1] - TTs[0]
        print(' | PID {0:d} || TID {1:d} sim {2:d}: Pgls0 = {3:.5f} Pgls1 = {4:.5f}.'
              .format(pid, tid, sim_num, Pgls0, Pgls1))
        print(' | PID {0:d} || TID {1:d} sim {2:d}: TTs[-1] = {3:.5f} TTs[0] = {4:.5f}.'
              .format(pid, tid, sim_num, TTs[-1], TTs[0]))
        fmin  = 1.0/Pgls1
        fmax  = 1.0/Pgls0
        frequencies = np.linspace(fmin, fmax, num=nfreq, endpoint=True)
        ogls  = gls.Gls((TTs, oc, eTTs), 
                        # Pbeg=Pgls0, Pend=Pgls1,
                        fbeg=fmin, fend=fmax,
                        freq=frequencies,
                        ofac=1,
                        verbose=False
                       )
        attv = ogls.hpstat['amp']
        
        if(attv <= 60*cst.sec2day):
          oc_unit = [cst.day2sec,  "sec"]
        elif(attv > cst.min2day  and attv <= 60*cst.min2day):
          oc_unit = [cst.day2min,  "min"]
        elif(attv > cst.hour2day and attv <= 24*cst.hour2day):
          oc_unit = [cst.day2hour, "hour"]
        else:
          oc_unit = [1.0,          "day"]

        pttv = 1.0/ogls.hpstat['fbest']

        chisquare     = np.sum( np.power(oc/eTTs, 2))
        dof           = nTTs - 2 # 2 parameters of the linear ephemeris
        red_chisquare = chisquare/dof

        # plot O-C
        
        fig = plt.figure()

        plt.axhline(0., color='black', ls='-', lw=0.33, alpha=0.77)
        plt.plot(TTs, oc*oc_unit[0],
                 color='lightgray', marker='o', ms=msize, mfc='C0',
                 ls='-', lw=0.45
                )

        plt.title('{0:05d} (O-C) body {1:02d}'.format(sim_num, i_body), loc='left',
                  fontsize=tsize
                  )
       
        try:
          title1='Tlin = {:.6f} + N x {:.6f}'.format(Tref, Pref)
        except:
          title1='ERROR lin eph.'
        # title2=r'$\sigma_\textrm{{TTV}}={0:.3f}$ {1:s} $P_\textrm{{TTV}}={2:.4f}$ d'.format(
        #        attv*oc_unit[0], oc_unit[1], pttv
        #        )
        try:
          title2='A(TTV) = {0:.3f} {1:s} P(TTV) = {2:.4f} d'.format(
                attv*oc_unit[0], oc_unit[1], pttv
                )
        except:
          title2='ERROR TTV'
        plt.title('{:s}\n{:s}'.format(title1, title2),
                  loc='right', fontsize=tsize2
                  )
        
        plt.xlabel('days',
                   fontsize=lsize
                  )
        plt.ylabel('O-C ({:s})'.format(oc_unit[1]), 
                   fontsize=lsize
                  )
        
        fig.savefig('{:s}.png'.format(plot_file), bbox_inches='tight')
        plt.close(fig)

      else:
        print(' | PID {0:d} || TID {1:d} sim {2:d}: found less than 4 transits for body {3:d}.'
              .format(pid, tid, sim_num, i_body))
        
        attv = -1
        pttv = -1
        Tref = -1
        Pref = -1

        chisquare     = -1
        dof           = -1
        red_chisquare = -1

        # plot EMPTY O-C
        fig = plt.figure()
        plt.title('{0:05d} (O-C) body {1:02d}'.format(sim_num, i_body), loc='left',
                    fontsize=tsize
                    )
        plt.text(0.5, 0.5, 'NO O-C PLOT: not enough transits',
                 horizontalalignment='center', verticalalignment='center',
                 fontsize=lsize
                 )
        fig.savefig('{:s}.png'.format(plot_file), bbox_inches='tight')
        plt.close(fig)

      A_TTV[i_body] = attv
      P_TTV[i_body] = pttv
      linear_ephem[i_body] = [Tref, Pref]
      statistics[i_body] = [chisquare, dof, red_chisquare]

      print(' | PID {0:d} || TID {1:d} sim {2:d}: body {3:d} A_TTV = {4:.2f} {5:s} , P_TTV = {6:.4f} days , Red Chi Square = {7:.1f}/{8:d} = {9:.3f}.'
              .format(pid, tid, sim_num, i_body, attv*oc_unit[0], oc_unit[1], pttv, chisquare, dof, red_chisquare))

  results_sim['A_TTV']        = A_TTV
  results_sim['P_TTV']        = P_TTV
  results_sim['linear_ephem'] = linear_ephem
  results_sim['statistics']   = statistics
    
  #anc.print_both('\n | sim %06d: COMPLETED\n' %(sim_num), of_log)
   
  return results_sim

# -------------------------------------

def run_the_simulation(sim_args):
  
  i_grid        = sim_args['i_grid']
  perturber_par = sim_args['perturber_par']
  n_bodies      = sim_args['n_bodies']
  t_start       = sim_args['t_start']
  t_end         = sim_args['t_end']
  err_T0        = sim_args['err_T0']
  names_lst     = sim_args['names_lst']
  units         = sim_args['units']
  mr_type       = sim_args['mr_type']
  plot_folder   = sim_args['plot_folder']
  
  pid = os.getpid()
  tid = _thread.get_ident()
  print(' PID {0:d} || TID {1:d} RUNNING sim {2:d} ...'.format(pid, tid, i_grid+1))
  
  results_sim = simulation_and_ttv_analysis(i_grid, n_bodies,
                perturber_par, t_start, t_end, err_T0, names_lst, units, mr_type,
                plot_folder,
                of_log=None
                )
  
  print(' PID {0:d} || TID {1:d} COMPLETED sim {2:d}'.format(pid, tid, i_grid+1))

  # del i_grid, transit_id, perturber_id, perturber_par, n_bodies, t_start, t_end
  # del err_T0, nT0_sel, names_lst, units, mr_type, 
  # gc.collect()

  return results_sim

# -------------------------------------

# ==============================================================================

# -------------------------------------

def compute_limits(xx, exx=None, delta=0.):
  
  if(exx is not None):
    xxmin = np.min(xx-exx)
    xxmax = np.max(xx+exx)
  else:
    xxmin = np.min(xx)
    xxmax = np.max(xx)
  dxx = np.abs(xxmax-xxmin)
  inflim = xxmin - dxx*delta
  suplim = xxmax + dxx*delta
  
  return inflim, suplim

# -------------------------------------

def insert_rms_lines(ax, rms_ref_m, rms_max, colors, lstyle):
  
  nref = np.shape(rms_ref_m)[0]
  
  for i_r in range(0, nref):
    if(rms_ref_m[i_r]/1440. <= rms_max):
      ax.axvline(rms_ref_m[i_r]/1440.0,
        color=colors[i_r], ls=lstyle[i_r], lw=0.7,
        alpha=0.99, zorder=3,
        label='rms = {0:.0f} min'.format(rms_ref_m[i_r])
        )
  
  return ax

# -------------------------------------

def insert_grey_area(ax, min_x, max_x, z_order=8):
  
  ax.axvspan(min_x, max_x,
    color='gray', alpha=0.45,
    zorder=z_order
    )

  return ax

# -------------------------------------


# ==============================================================================
# ==============================================================================


def ticklabels_logtolin(tickslabels_obj):
  
  ticks_log_str = [tickslabels_obj[ii].get_text() for ii in range(len(tickslabels_obj))]
  nticks = len(ticks_log_str)
  ticks_lin_str = ['']*nticks
  for ii in range(nticks):
    if(len(ticks_log_str[ii]) > 0):
      tick_temp = np.power(10., np.float64(ticks_log_str[ii].strip('$')))
      ticks_lin_str[ii] = r'${.2f}$'.format(tick_temp)
  
  return ticks_lin_str

# ==============================================================================

def plot_MMR(ax, P_tra):
  
  MMR_1st = [2.0, 3.0/2.0, 4.0/3.0, 5.0/4.0]
  MMR_2nd = [3.0, 5.0/3.0]
  MMR_3rd = [4.0, 5.0/2.0]
  
  # MMR j:j-1 // j:j+1
  for ii in range(0, len(MMR_1st)):
    Pin = P_tra/MMR_1st[ii]
    ax.axvline(Pin, color='red', ls='-', lw=1.2, alpha=0.66, zorder=7)
    #lgPin = np.log10(Pin)
    #ax.axvline(lgPin, color='red', ls='-', lw=1.2, alpha=0.66, zorder=7)
    Pout = P_tra*MMR_1st[ii]
    ax.axvline(Pout, color='red', ls='-', lw=1.2, alpha=0.66, zorder=7)
    #lgPout = np.log10(Pout)
    #ax.axvline(lgPout, color='red', ls='-', lw=1.2, alpha=0.66, zorder=7)
    
  # MMR j:j-2 // j:j+2
  for ii in range(0, len(MMR_2nd)):
    Pin = P_tra/MMR_2nd[ii]
    ax.axvline(Pin, color='orange', ls='-', lw=0.95, alpha=0.66, zorder=7)
    #lgPin = np.log10(Pin)
    #ax.axvline(lgPin, color='orange', ls='-', lw=0.95, alpha=0.66, zorder=7)
    
    Pout = P_tra*MMR_2nd[ii]
    ax.axvline(Pout, color='orange', ls='-', lw=0.95, alpha=0.66, zorder=7)
    #lgPout = np.log10(P_tra*MMR_2nd[ii])
    #ax.axvline(lgPout, color='orange', ls='-', lw=0.95, alpha=0.66, zorder=7)
    
  # MMR j:j-3 // j:j+3
  for ii in range(0, len(MMR_3rd)):
    Pin = P_tra/MMR_3rd[ii]
    ax.axvline(Pin, color='cyan', ls='-', lw=0.8, alpha=0.66, zorder=7)
    #lgPin = np.log10(P_tra/MMR_3rd[ii])
    #ax.axvline(lgPin, color='cyan', ls='-', lw=0.8, alpha=0.66, zorder=7)
    Pout = P_tra*MMR_3rd[ii]
    ax.axvline(Pout, color='cyan', ls='-', lw=0.8, alpha=0.66, zorder=7)
    #lgPout = np.log10(P_tra*MMR_3rd[ii])
    #ax.axvline(lgPout, color='cyan', ls='-', lw=0.8, alpha=0.66, zorder=7)
  
  return ax

# ==============================================================================


# needed by pp
class object_sim:
  def __init__(self, sim_args):
    self.sim_args=sim_args

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
  
  out_folder, log_file, of_log = set_analysis_folder_log(cli.full_path, cli.sub_folder)
  
  anc.print_both('# THIS IS: {:s}'.format(log_file), of_log)
  anc.print_both('Main folder: {:s}'.format(cli.full_path), of_log)
  
  err_T0 = set_err_T0(cli.err_T0)
  anc.print_both('Mean T0 uncertainty = {0:.6f} days'.format(err_T0), of_log)
  
  anc.print_both('TIME SELECTION: START = {0:.4f} (days) END = {1:.4f} (days)'.format(cli.t_start, cli.t_end), of_log)
  anc.print_both(' ', of_log)
  
  # init TRADES
  anc.print_both(' - INIT TRADES', of_log)
  
  f90trades.initialize_trades(cli.full_path, cli.sub_folder, cli.nthreads)
  plot_folder = os.path.join(out_folder, "plots")
  if(not os.path.isdir(plot_folder)):
    os.makedirs(plot_folder)

  n_bodies = f90trades.n_bodies # NUMBER OF TOTAL BODIES OF THE SYSTEM
  
  # read grid parameters ==> output ngrid and ncol
  ngrid, ncol = f90trades.wrapper_read_grid(1)
  # create grid
  perturber_grid = f90trades.wrapper_grid_init(1, ngrid, ncol)
  anc.print_both(' - SET PERTURBER GRID', of_log)
  #perturber_grid = f90trades.wrapper_grid_init(1)
  anc.print_both(' | To run {:d} simulation/s.\n'.format(ngrid), of_log)
  
  # hdf5_name, f_h5, _, _, names_lst, _, units = create_hdf5_file(out_folder, cli.seed, cli.perturber_id, perturber_grid, cli.mr_type, cli.t_start, cli.t_end, err_T0, of_log)
  # anc.print_both(' - created hdf5 file: {:s}'.format(hdf5_name), of_log)
  names_lst = {}
  names_str = {}
  units     = {}
  mr_convert, mr_unit = mass_factor(cli.mr_type)
  for i in range(1,n_bodies):
    i_body = i + 1
    n_lst, n_str, n_u = set_names_parameters(i_body, mr_unit)
    names_lst[i_body] = n_lst
    names_str[i_body] = n_str
    units[i_body]     = n_u

  niter = ngrid
  
  # # 2018-03-08 ADDING GLS AMPLITUDE AND PERIOD OF TTV/OC
  # attv_grid = np.zeros((niter))
  # pttv_grid = np.zeros((niter))
  
  
  print(' - Creating sim_args list ...', end=' ')
  
  list_sim_args = []
  list_osim_args = []
  
  for igrid in range(0, niter):
    sim_args                  = {}
    sim_args['i_grid']        = igrid
    sim_args['perturber_par'] = perturber_grid[igrid,:]
    sim_args['n_bodies']      = n_bodies
    sim_args['t_start']       = cli.t_start
    sim_args['t_end']         = cli.t_end
    sim_args['err_T0']        = err_T0
    sim_args['names_lst']     = names_lst
    sim_args['names_str']     = names_str
    sim_args['units']         = units
    sim_args['mr_type']       = cli.mr_type
    sim_args['plot_folder']   = plot_folder
    list_sim_args.append(sim_args)
    list_osim_args.append(object_sim(sim_args))
  print('done')
  
  print(' - *** RUN THEM ALL *** ...\n')
  
  if(cli.nthreads <= 1):
  # serial
    list_results_sim = [run_the_simulation(list_sim_args[igrid])
                        for igrid in range(niter)
                        ]
  else:
    
    with mpro.Pool(cli.nthreads) as threads_pool:
      list_results_sim = threads_pool.map(run_the_simulation, list_sim_args)

  print('done\n')

  print('Saving results ...')  
  
  # prepare files for summary
  n_planets         = n_bodies - 1 # NUMBER OF PLANETS IN THE SYSTEM
  summary_grid      = np.zeros((n_planets, ngrid, 12))
  rv_perturber_grid = np.zeros((ngrid, 4))

  for i in range(1,n_bodies):
    i_body = i + 1
    s_file = os.path.join(out_folder, "summary_NB{:02d}_grid.dat".format(i_body))
    header = "sim_number mass_NB{1:02d}_{0:s} period_NB{1:02d}_d mass_perturber_{0:s} period_perturber_d Tref Pref A_TTV_d P_TTV_d chi_square dof red_chi_square \n".format(mr_unit[0], i_body)

    for i_grid in range(0,niter):
      r = list_results_sim[i_grid]
          
      summary_grid[i-1, i_grid, 0]  = r['sim_num']
      summary_grid[i-1, i_grid, 1]  = r["keplerian_elements"][i,0]*mr_convert[0] # mass
      summary_grid[i-1, i_grid, 2]  = r["keplerian_elements"][i,2] # period
      summary_grid[i-1, i_grid, 3]  = r["perturber_par"][0]*mr_convert[0] # mass
      summary_grid[i-1, i_grid, 4]  = r["perturber_par"][2] # period
      summary_grid[i-1, i_grid, 5]  = r["linear_ephem"][i_body][0] # Tref
      summary_grid[i-1, i_grid, 6]  = r["linear_ephem"][i_body][1] # Pref
      summary_grid[i-1, i_grid, 7]  = r["A_TTV"][i_body] # A_TTV days
      summary_grid[i-1, i_grid, 8]  = r["P_TTV"][i_body] # P_TTV days
      summary_grid[i-1, i_grid, 9]  = r["statistics"][i_body][0] # chi square
      summary_grid[i-1, i_grid, 10] = r["statistics"][i_body][1] # dof
      summary_grid[i-1, i_grid, 11] = r["statistics"][i_body][2] # red chi square
      
      if(i_body == f90trades.idpert):
        rv_perturber_grid[i_grid, 0] = r['sim_num']
        rv_perturber_grid[i_grid, 1] = r["perturber_par"][0]*mr_convert[0] # mass
        rv_perturber_grid[i_grid, 2] = r["perturber_par"][2] # period
        rv_perturber_grid[i_grid, 3] = r["K_rv"][i_body]

    np.savetxt(s_file, summary_grid[i-1, :, :],
              fmt='%05.0f %23.16e %23.16e %23.16e %23.16e %23.16e %23.16e %23.16e %23.16e %23.16e %4d %23.16e ',
              header=header
              )
    anc.print_both('Written file {}'.format(s_file))
  
  rv_file   = os.path.join(out_folder, "summary_perturber_Krv.dat")
  header_rv = "sim_num mass_perturber_{0:s} period_perturber_d K_rv_mps".format(mr_unit[0])
  np.savetxt(rv_file, rv_perturber_grid, 
             fmt="%05.0f %23.16e %23.16e %23.16e ",
             header=header_rv
             )
  anc.print_both('Written file {}'.format(rv_file))

  # #f_h5 = save_chi2r_rms_analysis(f_h5, chi2r_oc_grid, rms_oc_grid, above_threshold, nT0_sel, rms_rv_grid)
  # f_h5 = save_chi2r_rms_analysis(f_h5, chi2r_oc_grid, rms_oc_grid, above_threshold,
  #        nT0_sel, rms_rv_grid, attv_grid, pttv_grid
  #        )
  
  # f_h5.close()
    
  anc.print_both(' ', of_log)
  anc.print_both(' ===============================', of_log)
  anc.print_both(' ', of_log)
  
  # --------------------

  of_log.close()
  
  return

# ==============================================================================
# RUN MAIN

if __name__ == "__main__":
  main()

# ==============================================================================

