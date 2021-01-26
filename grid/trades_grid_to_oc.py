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
mpl.use("Agg") #, warn=False)
#mpluse("Qt4Agg")
#mpluse("TkAgg")
import matplotlib.pyplot as plt
# matplotlib rc params
plt.rcParams['text.usetex'] = False
plt.rcParams['font.family'] = 'sans-serif'
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
from pytrades_lib import pytrades

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

  parser.add_argument('-n', '-n-T0', '--n-T0',
    action='store', dest='n_T0', default='None',
    help='Number of Transit Times (T0) to select, between time start and end, to compute TTV amplitude and period.'+\
         'Default is None, meaning it will take all the T0s in the time range.'
    )
  parser.add_argument('-mc', '-n-MC', '--n-MC',
    action='store', dest='n_mc', default=1,
    help='Number of Monte-Carlo repetition of selection of T0s. Default 1.'
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

  try:
    cli.n_T0 = int(cli.n_T0)
  except:
    cli.n_T0 = None
  
  try:
    cli.n_mc = int(cli.n_mc)
  except:
    cli.n_mc = 1


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

def simulation_and_ttv_analysis(i_grid, n_bodies, perturber_par, t_start, t_end, \
  err_T0, n_T0, n_mc, names_lst, units, mr_type, plot_folder, of_log=None):
  
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

  perturber_id = pytrades.idpert
  dof = 1

  try:
    #mass, radius, period, sma, ecc, argp, mean_an, inc, long_n, nt0_full, nrv_nmax = pytrades.wrapper_set_grid_orbelem(sim_num, perturber_id, perturber_grid, n_bodies)
    mass, radius, period, sma, ecc, argp, mean_an, inc, long_n, nt0_full, nrv_nmax = \
      pytrades.wrapper_set_grid_orbelem(sim_num, perturber_id, perturber_par, n_bodies)
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
      pytrades.wrapper_run_grid_combination(mass, radius, period, sma, ecc, argp, mean_an, inc, long_n,
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
  
  A_TTV_MC     = {} # median
  P_TTV_MC     = {}
  A_MC         = {} # MC
  P_MC         = {}

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

      if(n_T0 is not None):
        nTTsel = n_T0
      else:
        nTTsel = nTTs

      A_MC[i_body] = np.zeros((n_mc)) - 1.0
      P_MC[i_body] = np.zeros((n_mc)) - 1.0

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

        print(' | PID {0:d} || TID {1:d} sim {2:d}: running {3:d} MC for body {4:d}'
              .format(pid, tid, sim_num, n_mc, i_body))
        # MC analysis
        for i_mc in range(n_mc):
          if(nTTsel >= nTTs):
            selTTs = np.arange(0,nTTs).astype(int)
          else:
            selTTs = np.random.choice(
              np.arange(0,nTTs).astype(int),
              size=nTTsel,
              replace=False
            )
          
          # too few data-points to compute properly the GLS amplitude and period
          # ogls  = gls.Gls((TTs[selTTs], oc[selTTs], eTTs[selTTs]), 
          #                 # Pbeg=Pgls0, Pend=Pgls1,
          #                 fbeg=fmin, fend=fmax,
          #                 freq=frequencies,
          #                 ofac=1,
          #                 verbose=False
          #               )
          # amp = ogls.hpstat['amp']
          # per = 1.0/ogls.hpstat['fbest']

          # computing amp of TTV as semi-amplitude of refitted OC
          epox, Trefx, Prefx, _ = anc.compute_lin_ephem(
            TTs[selTTs],
            eT0=eTTs[selTTs],
            epoin=epo[selTTs],
            modefit='wls'
          )
          TTlinx = Trefx + epox * Prefx
          ocx    = TTs[selTTs] - TTlinx
          

          A_MC[i_body][i_mc] = 0.5*(np.max(ocx)-np.min(ocx))
          P_MC[i_body][i_mc] = pttv # setting equal to full gls analysis

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
      A_TTV_MC[i_body] = np.median(A_MC[i_body])
      P_TTV_MC[i_body] = np.median(P_MC[i_body])

      print(' | PID {0:d} || TID {1:d} sim {2:d}: body {3:d} A_TTV = {4:.2f} {5:s} , P_TTV = {6:.4f} days , Red Chi Square = {7:.1f}/{8:d} = {9:.3f}.'
              .format(pid, tid, sim_num, i_body, attv*oc_unit[0], oc_unit[1], pttv, chisquare, dof, red_chisquare))

  results_sim['A_TTV']        = A_TTV
  results_sim['P_TTV']        = P_TTV
  results_sim['linear_ephem'] = linear_ephem
  results_sim['statistics']   = statistics
  results_sim['A_TTV_MC']        = A_TTV_MC
  results_sim['P_TTV_MC']        = P_TTV_MC
    
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
  n_T0          = sim_args['n_T0']
  n_mc        = sim_args['n_mc']
  names_lst     = sim_args['names_lst']
  units         = sim_args['units']
  mr_type       = sim_args['mr_type']
  plot_folder   = sim_args['plot_folder']
  
  pid = os.getpid()
  tid = _thread.get_ident()
  print(' PID {0:d} || TID {1:d} RUNNING sim {2:d} ...'.format(pid, tid, i_grid+1))
  
  results_sim = simulation_and_ttv_analysis(i_grid, n_bodies,
                perturber_par, t_start, t_end, err_T0, n_T0, n_mc,
                names_lst, units, mr_type,
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
  anc.print_both('Number of T0 to select: {}'.format(cli.n_T0))
  anc.print_both('Number of MC analysis: {}'.format(cli.n_mc))
  anc.print_both(' ', of_log)
  
  # init TRADES
  anc.print_both(' - INIT TRADES', of_log)
  
  pytrades.initialize_trades(cli.full_path, cli.sub_folder, cli.nthreads)
  plot_folder = os.path.join(out_folder, "plots")
  if(not os.path.isdir(plot_folder)):
    os.makedirs(plot_folder)

  n_bodies = pytrades.n_bodies # NUMBER OF TOTAL BODIES OF THE SYSTEM
  
  # read grid parameters ==> output ngrid and ncol
  ngrid, ncol = pytrades.wrapper_read_grid(1)
  # create grid
  perturber_grid = pytrades.wrapper_grid_init(1, ngrid, ncol)
  anc.print_both(' - SET PERTURBER GRID', of_log)
  #perturber_grid = pytrades.wrapper_grid_init(1)
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
    sim_args['n_T0']          = cli.n_T0
    sim_args['n_mc']          = cli.n_mc
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
  summary_grid      = np.zeros((n_planets, ngrid, 14))
  rv_perturber_grid = np.zeros((ngrid, 4))

  for i in range(1,n_bodies):
    i_body = i + 1
    s_file = os.path.join(out_folder, "summary_NB{:02d}_grid.dat".format(i_body))
    header = "sim_number mass_NB{1:02d}_{0:s} period_NB{1:02d}_d mass_perturber_{0:s} period_perturber_d Tref Pref A_TTV_d P_TTV_d chi_square dof red_chi_square A_TTV_MC_d P_TTV_MC_d\n".format(mr_unit[0], i_body)

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
      summary_grid[i-1, i_grid, 12]  = r["A_TTV_MC"][i_body] # A_TTV_MC days
      summary_grid[i-1, i_grid, 13]  = r["P_TTV_MC"][i_body] # P_TTV_MC days
      
      if(i_body == pytrades.idpert):
        rv_perturber_grid[i_grid, 0] = r['sim_num']
        rv_perturber_grid[i_grid, 1] = r["perturber_par"][0]*mr_convert[0] # mass
        rv_perturber_grid[i_grid, 2] = r["perturber_par"][2] # period
        rv_perturber_grid[i_grid, 3] = r["K_rv"][i_body]

    np.savetxt(s_file, summary_grid[i-1, :, :],
              fmt='%05.0f %23.16e %23.16e %23.16e %23.16e %23.16e %23.16e %23.16e %23.16e %23.16e %4d %23.16e %23.16e %23.16e ',
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

# # ==============================================================================
# # ==============================================================================
# # READ AND PLOT FUNCTIONS

# # TTV data
# class TTVdata:
#   def __init__(self, TTV_file):
#     self.sim_num, _, _, \
#       self.Mpert, self.Ppert, _, _, \
#         self.ATTV, self.PTTV, _, self.dof, self.rchisqr, \
#           self.ATTV_MC, self.PTTV_MC = np.genfromtxt(TTV_file, unpack=True)

# # RV data
# class RVdata:
#   def __init__(self, RV_file):
#     self.sim_num, self.Mpert, self.Ppert, self.Krv_mps = np.genfromtxt(RV_file, unpack=True)

# # Plot/s
# def get_limits(v, scale):
#   vmin = np.min(v)
#   vmax = np.max(v)
#   dv   = vmax - vmin
#   lim_min = vmin - dv * scale
#   lim_max = vmax + dv * scale
#   return [lim_min, lim_max]

# def scatter_plot_TTV(TTV):
  
#   fig = plt.figure()
  
#   vir_map = copy.copy(plt.cm.get_cmap("viridis_r"))
  
#   plt.scatter(TTV.Ppert, TTV.Mpert, c=TTV.ATTV*cst.day2min,
#               vmin=0.0,
#               vmax=10.0,
#               cmap=vir_map,
#              )
#   plt.colorbar(label=r'A$_\mathrm{TTV}$ (min)')
#   plt.xlabel(r"$P_\mathrm{perturber}$ (days)")
#   plt.ylabel(r"$M_\mathrm{perturber}$ ($M_\oplus$)")
#   plt.show()
#   return fig

# def grid_plot(TTV, RV, grid_on=True, rchisq=False):
    
#   n_P = np.shape(np.unique(TTV.Ppert))[0]
#   n_M = np.shape(np.unique(TTV.Mpert))[0]
  
#   mesh_P   = TTV.Ppert.reshape((n_M, n_P))
#   mesh_M   = TTV.Mpert.reshape((n_M, n_P))
#   mesh_TTV = TTV.ATTV.reshape((n_M, n_P))*cst.day2min
#   mesh_RV  = RV.Krv_mps.reshape((n_M, n_P))
#   mesh_c2r = TTV.rchisqr.reshape((n_M, n_P))
  
#   fig, ax = plt.subplots()

#   # TTV amplitude in min
#   TTVlevels_min = [np.float(i) for i in range(0,11)] + [np.max(TTV.ATTV)*cst.day2min]
#   vir_map = copy.copy(plt.cm.get_cmap("viridis_r"))
#   vir_map.set_under("white")
#   vir_map.set_over("black")
  
#   nTTVlev = np.shape(TTVlevels_min)[0]
#   TTV_colors = vir_map(np.linspace(0.2,1.0,num=nTTVlev,endpoint=True))
  
#   TTV_ax = ax.contourf(mesh_P, mesh_M, mesh_TTV,
#                        levels=TTVlevels_min,
#                        colors=TTV_colors,
#                        origin=None,
#                        extend='max',
#                        zorder=5
#                       )
  
#   TTV_bar = fig.colorbar(TTV_ax)
#   TTV_bar.set_label(r'A$_\mathrm{TTV}$ (min)')
  
#   # RV amplitude in m/s
# #   RVlevels_mps = [1.0, 2.0, 3.0, 4.0, 5.0, 6.0] # m/s
#   RVlevels_mps = [float(i) for i in range(1, 11)] + [int(0.5*(10.0+np.max(RV.Krv_mps)))]  # m/s
# #   print(RVlevels_mps)
#   nRVlev       = len(RVlevels_mps)
#   gval         = 0.5
#   RV_colors    = [(gval, gval, gval, i) for i in np.linspace(0.0, 1.0, num=nRVlev, endpoint=True)]
  
# #   RV_ax = ax.contourf(mesh_P, mesh_M, mesh_RV,
# #                       levels=RVlevels_mps,
# #                       colors=RV_colors,
# #                       origin=None,
# #                       extend='max',
# #                       zorder=6
# #                      )
  
# #   RV_bar = fig.colorbar(RV_ax)
# #   RV_bar.set_label(r'$K_\mathrm{RV}$ (m/s)')
  
#   RV_lines = ax.contour(mesh_P, mesh_M, mesh_RV,
#     levels=RVlevels_mps,
#     #                         colors=RV_colors,
#     colors='white',
#     origin=None,
#     extend='max',
#     linewidths=0.5,
#     antialiased=True,
#     zorder=7
#   )
#   clabels = ax.clabel(RV_lines, inline=True, fontsize=8, colors='white', fmt='%.0f m/s')
#   for clabel in clabels:
#     clabel.set_zorder(8)
  
#   if(rchisq):
    
# #     # \chi^_r as lines contour
#     c2rlevels = [0.0, 1.0, 2.0, 3.0, 5.0]
#     nc2rlev   = len(c2rlevels)  

#     gval         = 0.3
#     c2r_colors    = [(gval, gval, gval, i) for i in np.linspace(0.0, 0.7, num=nc2rlev, endpoint=True)]
#     c2r = ax.contourf(mesh_P, mesh_M, mesh_c2r,
#                        levels=c2rlevels,
#                        colors=c2r_colors,
#                        origin=None,
#                        extend='max',
#                        zorder=6
#                       )
# #     cbaxes  = fig.add_axes([0.0, 1.0, 0.8, 0.03]) # [left, bottom, width, height]
#     c2r_bar = fig.colorbar(c2r, orientation="horizontal",
# #                            cax=cbaxes
#                            anchor=(0.0, 0.0)
#     )
#     c2r_bar.set_label(r'$\chi^2_r$')
    
  
#   if(grid_on): ax.plot(TTV.Ppert, TTV.Mpert, color='gray', marker='o', ms=1, mec='None', ls='', alpha=0.5, zorder=8)
  
  
#   scale = 0.01
#   ax.set_xlim(get_limits(TTV.Ppert, scale))
#   ax.set_ylim(get_limits(TTV.Mpert, scale))
  
#   ax.set_xlabel(r"$P_\mathrm{perturber}$ (days)")
#   ax.set_ylabel(r"$M_\mathrm{perturber}$ ($M_\oplus$)")
#   plt.tight_layout()
#   plt.show()
  
#   return fig

# # ==============================================================================
# # ==============================================================================

# ==============================================================================
# RUN MAIN

if __name__ == "__main__":
  main()

# ==============================================================================

