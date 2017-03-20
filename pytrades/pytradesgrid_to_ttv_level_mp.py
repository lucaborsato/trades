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

from __future__ import division # no more "zero" integer division bugs!:P
import argparse
import numpy as np # array
import h5py
import astropy.time as atime
import os
import sys

from matplotlib import use as mpluse
mpluse("Agg")
#mpluse("Qt4Agg")
#mpluse("TkAgg")
import matplotlib.pyplot as plt
plt.rc('font',**{'family':'serif','serif':['Computer Modern Roman']})
plt.rc('text', usetex=True)

#import matplotlib.mlab as mlab
import matplotlib.cm as cm
#import matplotlib.colors as mc
#from mpl_toolkits.axes_grid1 import make_axes_locatable

import ancillary as anc
import constants as cst
from pytrades_lib import pytrades
#import pytrades_lib.pytrades as pytrades

import time
import multiprocessing as mp

# ==============================================================================

# read cli arguments
def get_args():
  
  parser = argparse.ArgumentParser(description='CHEOPS TTV LEVEL ANALYSIS')
  parser.add_argument('-p', '--path', action='store', dest='full_path', required=True, help='The path (absolute or relative) with simulation files for TRADES.')
  
  #parser.add_argument('-lm', '--lm-flag', action='store', dest='lm_flag', default='0', help='LM flag. If you want to read files w/o LM set to 0, otherwise 1. Default is 0')
  
  # SUBFOLDER TO SAVE ALL SIMULATION DATA
  parser.add_argument('-sf', '--sub-folder', '--sub-f', action='store', dest='sub_folder', default='grid2ttv', help='Sub-folder name, without full path. Default = grid2ttv')
  
  parser.add_argument('-e', '--e-T0', '--err-T0', action='store', dest='err_T0', default='0.', help='Mean uncertainty on the Transit Times (T0s). Provide value with unit "d" "m" "s" for days, minutes, and seconds in this way: "0.001528d" or "2.2m" or "132.s". Default unit is days and value set to 0.')
  
  parser.add_argument('-mr', '--mr-type', '--mass-radius-type', action='store', dest='mr_type', default='e', help='Mass and radius type to convert the grid. "e ea ear eart earth" or "j ju jup jupi jupit jupite jupiter" or "n ne nep nept neptu neptun neptune".')
  
  parser.add_argument('-s', '--seed', action='store', dest='seed', default='None', help='Seed for random number generator. Defauto is None.')
  
  parser.add_argument('-c', '--cpu', '--threads', action='store', dest='nthreads', default=1, type=int, help='Number of threads/cpu to use. Default is 1.')
  
  parser.add_argument('-id', '--id-perturber', action='store', dest='perturber_id', default=3, type=int, help='Id of the perturber body: star is 1, so it has to be from 2 to Number of bodies. Default is 3.')
  
  cli = parser.parse_args()
  cli.full_path = os.path.join(os.path.abspath(cli.full_path), '')
  cli.sub_folder = os.path.join(os.path.relpath(cli.sub_folder), '')
  #cli.lm_flag = str(cli.lm_flag)
  try:
    cli.seed = int(cli.seed)
  except:
    cli.seed = None
  try:
    cli.perturber_id = int(cli.perturber_id)
  except:
    cli.perturber_id = 3
  
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
    err_T0 = 0.
    
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

def set_names_parameters(id_body, mr_unit, of_log):
  
  body_num = '%d' %(id_body)
  base_names_lst = 'M R P a e w mA tau i lN'.split()
  units = '%s %s day au - deg deg day deg deg' %(mr_unit[0], mr_unit[1])
  units = units.split()
  n_names = len(base_names_lst)
  names_lst = ['%s%s' %(base_names_lst[i], body_num) for i in range(0, n_names)]
  temp_lst = ['%s [%s]' %(names_lst[i], units[i].strip()) for i in range(0, n_names)]
  names_str = ' '.join(['%23s' %(temp_lst[i]) for i in range(0, n_names)])
  
  return names_lst, names_str, units

# -------------------------------------

def select_and_sort_ttra(ttra_full, id_ttra_full, stats_ttra, id_transiting=2, t_min_sel=None, t_max_sel=None, of_log=None):
  
  sel_tra = np.logical_and(id_ttra_full == id_transiting, stats_ttra)
  ttra_all = np.sort(ttra_full[sel_tra])
  nttra_all = np.shape(ttra_all)[0]
  #anc.print_both(' | number of transits nT0 = %d' %(nttra_all), of_log)

  if(nttra_all == 0):
    return None, 0, None, 0, None, 0

  if(t_min_sel is None):
    sel_min = 0.
  else:
    sel_min = t_min_sel
  if(t_max_sel is None):
    sel_max = ttra_all[-1]
  else:
    sel_max = t_max_sel
    
  sel_pre = ttra_all < sel_min # T0s before CHEOPS mission, they are use to define Tref, Pref
  
  nttra_pre = int(np.sum(np.asarray(sel_pre, dtype=np.int)))
  #anc.print_both(' | selected pre-CHEOPS transits: nT0 (pre) = %d' %(nttra_pre), of_log)
  if(nttra_pre == 0):
    return None, 0, None, 0, None, 0
  
  sel_cheops = np.logical_and(ttra_all >= sel_min, ttra_all <= sel_max)
  nttra_cheops = int(np.sum(np.asarray(sel_cheops, dtype=np.int)))
  #anc.print_both(' | selected CHEOPS transits: nT0 (cheops) = %d' %(nttra_cheops), of_log)
  if(nttra_cheops == 0):
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
    eee = np.random.normal(loc=0., scale=sig_err)
    if(eee < min_err):
      dirty.append(eee)
      ierr += 1
      if(ierr == nT0):
        break
      
  eT0 = np.array(dirty, dtype=np.float64) + err_T0
  T0dirty = T0 + np.array(dirty, dtype=np.float64)
  
  return T0dirty, eT0

# -------------------------------------

# epoch or transit number for each T0 given a T0ref and a Pref
def calculate_epoch(T0, Tref, Pref):
  
  epo = np.rint((T0-Tref)/Pref).astype(int)
  
  return epo

# -------------------------------------

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

# -------------------------------------

def compute_lin_ephem(T0, eT0):
  
  nT0 = np.shape(T0)[0]
  Tref0 = T0[0]
  dT = [np.abs(T0[i+1]-T0[i]) for i in range(nT0-1)]
  Pref0 = np.percentile(np.array(dT), 50., interpolation='midpoint')
  
  epo = calculate_epoch(T0, Tref0, Pref0)
  Tref, Pref, coeff_err = lstsq_fit(epo, T0, eT0)
  epo = calculate_epoch(T0, Tref, Pref)

  return epo, Tref, Pref

# -------------------------------------

def rms(vec):
  
  rms_vec = np.sqrt(np.mean(vec*vec))
  
  return rms_vec
# -------------------------------------

def compute_oc(Tref, Pref, T0):
  
  epo = calculate_epoch(T0, Tref, Pref)
  T0lin = Tref + epo * Pref
  oc = T0 - T0lin
  #rms_oc = np.sqrt(np.mean(oc*oc))
  rms_oc = rms(oc)
  
  return epo, oc, rms_oc

# -------------------------------------

def compute_rms_mc(oc, nT0_sel, nMC):
  
  if(np.shape(oc)[0] < nT0_sel):
    return 0.0, 0.0
  rms_all = [rms(np.random.choice(oc, nT0_sel, replace=False)) for imc in range(0, nMC)]
  rms_mean = np.mean(rms_all)
  rms_med = np.percentile(rms_all, 50., interpolation='midpoint')
    
  return rms_mean, rms_med

# -------------------------------------

# ==============================================================================

# -------------------------------------

def create_hdf5_file(full_path, seed, id_body, perturber_grid, mr_type, t_sel_min, t_sel_max, err_T0, of_log=None):
  
  mr_convert, mr_unit = mass_factor(mr_type)
  names_lst, names_str, units = set_names_parameters(id_body, mr_unit, of_log)
  ngrid, ncol = np.shape(perturber_grid)
  write_grid = perturber_grid.copy()
  write_grid[:,0] = write_grid[:,0] * mr_convert[0] # mass
  write_grid[:,1] = write_grid[:,1] * mr_convert[1] # radius
  
  file_name = os.path.join(full_path, 'rms_analysis_summary.hdf5')
  f_h5 = h5py.File(file_name, 'w')
  
  f_h5.create_dataset('grid/grid_names', data=names_lst, dtype='S20', compression="gzip")
  f_h5.create_dataset('grid/grid_units', data=units, dtype='S20', compression="gzip")
  f_h5.create_dataset('grid/perturber_grid', data=write_grid, dtype=np.float64, compression="gzip")
  
  f_h5['grid'].attrs['ngrid'] = ngrid
  f_h5['grid'].attrs['seed'] = seed
  f_h5['grid'].attrs['t_sel_min'] = t_sel_min
  f_h5['grid'].attrs['t_sel_max'] = t_sel_max
  f_h5['grid'].attrs['err_T0_day'] = err_T0
  f_h5['grid'].attrs['err_T0_min'] = err_T0*1440.
  
  return file_name, f_h5, mr_convert, mr_unit, names_lst, names_str, units


# -------------------------------------

def save_results_sim(f_h5, results_sim, mr_type):
  
  mr_convert, mr_unit = mass_factor(mr_type)
  temp_par = results_sim['perturber_par'].copy()
  temp_par[0] = temp_par[0]*mr_convert[0]
  temp_par[1] = temp_par[1]*mr_convert[1]
  
  sim_num = results_sim['sim_num']
  f_h5.create_dataset('sim_%06d/perturber_par_names' %(sim_num), data=np.array(results_sim['perturber_par_names'], dtype='S10'), dtype='S20', compression="gzip")
  f_h5.create_dataset('sim_%06d/perturber_par_units' %(sim_num), data=np.array(results_sim['perturber_par_units'], dtype='S10'), dtype='S20', compression="gzip")
  f_h5.create_dataset('sim_%06d/perturber_par'%(sim_num), data=temp_par, dtype=np.float64, compression="gzip")
  f_h5.create_dataset('sim_%06d/MR_sun2unit'%(sim_num), data=np.array(mr_convert), dtype=np.float64, compression="gzip")
  

  f_h5.create_dataset('sim_%06d/keplerian_elements_header' %(sim_num), data=np.array(results_sim['keplerian_elements_header'], dtype='S10'), dtype='S20', compression="gzip")
  f_h5.create_dataset('sim_%06d/keplerian_elements'%(sim_num), data=results_sim['keplerian_elements'], dtype=np.float64, compression="gzip")
  
  #f_h5.create_dataset('sim_%06d/ttra_all'%(sim_num), data=results_sim['ttra_all'], dtype=np.float64, compression="gzip")
  #f_h5.create_dataset('sim_%06d/sel_pre_cheops'%(sim_num), data=results_sim['sel_pre_cheops'], dtype=bool, compression="gzip")
  #f_h5.create_dataset('sim_%06d/sel_cheops'%(sim_num), data=results_sim['sel_cheops'], dtype=bool, compression="gzip")
  #f_h5.create_dataset('sim_%06d/T0_eT0'%(sim_num), data=results_sim['T0_eT0'], dtype=np.float64, compression="gzip")
  
  #f_h5.create_dataset('sim_%06d/time_rv'%(sim_num), data=results_sim['time_rv'], dtype=np.float64, compression="gzip")
  
  f_h5.create_dataset('sim_%06d/linear_ephem' %(sim_num), data=results_sim['linear_ephem'] , dtype=np.float64, compression="gzip")
  f_h5.create_dataset('sim_%06d/rms_ttv' %(sim_num), data=np.array([results_sim['rms_ttv']]) , dtype=np.float64, compression="gzip")
  f_h5.create_dataset('sim_%06d/nT0_sel' %(sim_num), data=results_sim['nT0_sel'] , dtype=np.int, compression="gzip")
  f_h5.create_dataset('sim_%06d/rms_mc' %(sim_num), data=results_sim['rms_mc'] , dtype=np.float64, compression="gzip")

  f_h5['sim_%06d' %(sim_num)].attrs['sim_num'] = results_sim['sim_num']
  f_h5['sim_%06d' %(sim_num)].attrs['nMC'] = results_sim['nMC']
  
  return f_h5

# -------------------------------------

def save_rms_analysis(f_h5, rms_grid, above_threshold, nT0_sel):
  
  rms_header_str = 'rms_ttv_d %s' %(' '.join(['rms_mean_N%d_d rms_median_N%d_d' %(nT0_sel[i_n], nT0_sel[i_n]) for i_n in range(len(nT0_sel))]))
  rms_header_lst = rms_header_str.split()
  
  f_h5.create_dataset('grid/rms_grid', data=rms_grid, dtype=np.float64, compression="gzip")
  
  f_h5.create_dataset('grid/rms_header', data=np.array(rms_header_lst, dtype='S20'), dtype='S20', compression="gzip")
  
  f_h5.create_dataset('grid/above_threshold', data=above_threshold, dtype=bool, compression="gzip")
  
  return f_h5

# -------------------------------------

# ==============================================================================

# -------------------------------------

def simulation_and_ttv_analysis(i_grid, transit_id, perturber_id, n_bodies, perturber_par, t_launch, t_end, err_T0, nT0_sel, names_lst, units, mr_type, nMC=1000, of_log=None):
  
  results_sim = {}
  sim_num = i_grid + 1
  anc.print_both('\n - sim %06d: START\n' %(sim_num), of_log)
  
  #mass, radius, period, sma, ecc, argp, mean_an, inc, long_n, nt0_full, nrv_nmax = pytrades.wrapper_set_grid_orbelem(sim_num, perturber_id, perturber_grid, n_bodies)
  mass, radius, period, sma, ecc, argp, mean_an, inc, long_n, nt0_full, nrv_nmax = pytrades.wrapper_set_grid_orbelem(sim_num, perturber_id, perturber_par, n_bodies)
  #anc.print_both(' | Get orbital elements', of_log)
  
  #anc.print_both(' - Run integration ...', of_log)
  ttra_full, id_ttra_full, stats_ttra, time_rv_nmax, rv_nmax, stats_rv = pytrades.wrapper_run_grid_combination(mass, radius, period, sma, ecc, argp, mean_an, inc, long_n, nt0_full, nrv_nmax)
  #anc.print_both(' | ...done. Computed all transit times and rvs.', of_log)
  
  #anc.print_both(' - Select rv', of_log)
  time_rv_all, rv_all = select_and_sort_rv(time_rv_nmax, rv_nmax, stats_rv.astype(bool))
  if(time_rv_all is None):
    time_rv = np.zeros((2))
  else:
    time_rv = np.column_stack((time_rv_all, rv_all))

  #anc.print_both(' - Select transit times:', of_log)
  ttra_all, nttra_all, sel_pre, nttra_pre, sel_cheops, nttra_cheops = select_and_sort_ttra(ttra_full, id_ttra_full, stats_ttra.astype(bool), id_transiting=transit_id,  t_min_sel=t_launch, t_max_sel=t_end, of_log=of_log)
  
  if(ttra_all is not None):
    T0, eT0 = dirty_transits(ttra_all, err_T0)
    #anc.print_both(' | dirty transits done', of_log)
    #anc.print_both(' | new mean T0 error = %.6f day = %.6f min' %(np.mean(eT0), np.mean(eT0)*1440.), of_log)
    T0_eT0 = np.column_stack((T0, eT0))
    
    epo_pre, Tref_pre, Pref_pre = compute_lin_ephem(T0[sel_pre], eT0[sel_pre])
    #anc.print_both(' - computed the linear ephemeris on pre-CHEOPS transits', of_log)
    #anc.print_both(' | T0lin = Tref + N x Pref = %.8f + N x %.8f' %(Tref_pre, Pref_pre), of_log)
    
    epo_cheops, oc_cheops, rms_cheops = compute_oc(Tref_pre, Pref_pre, T0[sel_cheops])
    #anc.print_both(' - computed O-C = T0cheops - T0lin', of_log)
    #anc.print_both(' | rms_TTV (cheops) = %.10f day = %.8f hour = %.6f min %.3f sec' %(rms_cheops, rms_cheops*24., rms_cheops*1440., rms_cheops*86400.), of_log)
    
    #anc.print_both(' - computing Monte-Carlo rms of O-C for different number of T0s', of_log)
    rms_mc = [list(compute_rms_mc(oc_cheops, nT0_sel[i_n], nMC)) for i_n in range(len(nT0_sel))]
    
  else: # if the transit file is empty
    #anc.print_both(' | NO transits at all ...')
    ttra_all = np.zeros((1))
    sel_pre = np.zeros((1)).astype(bool)
    sel_cheops = np.zeros((1)).astype(bool)
    T0_eT0 = np.column_stack((T0, eT0))
    Tref_pre = 0.
    Pref_pre = 0.
    rms_cheops = -1.
    rms_mc = -np.ones((len(nT0_sel), 2))
  
  results_sim['sim_num']             = sim_num
  results_sim['perturber_par_names'] = names_lst
  results_sim['perturber_par_units'] = units
  results_sim['perturber_par']       = perturber_par
  
  results_sim['keplerian_elements_header'] = 'Mass_Msun Radius_Rsun Period_d sma_au ecc argp_deg meanA_deg inc_deg longN_deg'.split()
  results_sim['keplerian_elements']        = np.column_stack((mass, radius, period, sma, ecc, argp, mean_an, inc, long_n))
  
  results_sim['ttra_all']            = ttra_all
  results_sim['sel_pre_cheops']      = sel_pre
  results_sim['sel_cheops']          = sel_cheops
  results_sim['T0_eT0']              = T0_eT0
  
  results_sim['time_rv'] = time_rv
  
  results_sim['linear_ephem'] = [Tref_pre, Pref_pre]
  results_sim['rms_ttv']      = rms_cheops
  results_sim['nT0_sel']      = nT0_sel
  results_sim['nMC']          = nMC
  results_sim['rms_mc']       = rms_mc  # mean and median OC rms for different nT0_sel during CHEOPS mission
  
  #f_h5 = save_results_sim(f_h5, results_sim, mr_type)
  
  #rms_grid[i_grid, 0]  = rms_cheops
  #rms_grid[i_grid, 1:] = np.reshape(rms_mc, newshape=(len(nT0_sel)*2))
  #above_threshold[i_grid, :] = rms_grid[i_grid, :] > err_T0
  
  anc.print_both('\n | sim %06d: COMPLETED\n' %(sim_num), of_log)
  
  return results_sim

# -------------------------------------

def get_rms_from_results(results_sim, rms_grid, above_threshold, nT0_sel, err_T0):
  
  i_grid = results_sim['sim_num'] - 1
  rms_ttv = results_sim['rms_ttv']
  rms_mc = results_sim['rms_mc']
  rms_grid[i_grid, 0]  = rms_ttv
  rms_grid[i_grid, 1:] = np.reshape(rms_mc, newshape=(len(nT0_sel)*2))
  above_threshold[i_grid, :] = rms_grid[i_grid, :] > err_T0
  
  return

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

def get_bins(x, rule='rice'):
  
  nx = np.shape(x)[0]
  nxf = np.float64(nx)
  
  if ('sq' in rule.lower()):
    
    k_bins = int(np.sqrt(np.float64(nx)))
    
  elif ('stu' in rule.lower()):
    
    k_bins = int(1. + np.log2(np.float64(nx)))
    
  elif ('doa' in rule.lower()):
    
    mux = np.mean(x, axis=0)
    stdx = np.std(x, ddof=1, axis=0)
    g1 = np.max(np.mean( ((x-mux)/stdx)**3 , axis=0)) # skew
    #g1 = np.min(np.mean( ((x-mux)/stdx)**3 , axis=0)) # skew
    s_g1 = np.sqrt( 6.*(nxf-2) / ((nxf+1.)*(nxf+3.)) )
    k_bins = int(1. + np.log2(nxf) + np.log2( 1. + np.abs(g1)/s_g1 ))
  
  else: # rice
    
    k_bins = int(np.ceil(2. * nxf**(1./3.)))
    
  return k_bins

# -------------------------------------

def insert_rms_lines(ax, rms_ref_m, rms_max, colors):
  
  nref = np.shape(rms_ref_m)[0]
  
  for i_r in range(0, nref):
    if(rms_ref_m[i_r]/1440. <= rms_max):
      ax.axvline(rms_ref_m[i_r]/1440., color=colors[i_r], ls='-', lw=0.7, alpha=0.99, zorder=3, label='rms = %.0f min' %(rms_ref_m[i_r]))
  
  return ax

# -------------------------------------

def plot_summary(out_folder, hdf5_file, body_id=2, mass_type='e', rms_id_in=None, rms_max_in = None, of_log=None):
  
  letters = 'a b c d e f g h i j k l m n o p q r s t u v x y z'.split()
  p_names = 'm P e w mA i lN'.split()
  npar = np.shape(p_names)[0]
  if(rms_id_in is None):
    rms_id = npar
  else:
    rms_id = rms_id_in
        
  l_names = 'M P e \\omega \\nu i \\Omega'.split()
  units = '[$M_\oplus$] [day] [] [$^\circ$] [$^\circ$] [$^\circ$] [$^\circ$]'.split()
  units[2] = ''
  if(mass_type in 'j ju jup jupi jupit jupit jupite jupiter'):
    units[0] = '[$M_\\textrm{Jup}$]'
  elif(mass_type in 'n ne nep nept neptu neptun neptune'):
    units[0] = '[$M_\\textrm{Nep}$]'
  
  pl_lab = ['$%s_\\textrm{%s}$ %s' %(l_names[i_n], letters[body_id], units[i_n]) for i_n in range(0, npar)]
  
  f_h5 = h5py.File(hdf5_file, 'r')
  final_summary = f_h5['summary/final_summary'][...]
  final_summary_above = f_h5['summary/final_summary_above'][...]
  above = f_h5['summary/above_threshold'][...]
  err_T0_d = f_h5['summary'].attrs['err_T0']
  f_h5.close()

  xx = final_summary[:,rms_id] # rms_ttv is the npar-th col of the final_summary grid
  xxa = final_summary_above[:,rms_id]
  
  if(rms_max_in is None):
    xxlim1, xxlim2 = compute_limits(xx, delta=0.05)
    xxmin = np.min(xx)
    xxmax = np.max(xx)
  else:
    xxmin = 0.
    xxmax = rms_max_in
    xxlim1, xxlim2 = compute_limits(np.array([xxmin, xxmax]), delta=0.05)
  
  rms_ref_m = np.array([1., 5., 10., 15., 20., 25., 30., 45., 60.])
  nref = np.shape(rms_ref_m)[0]
  
  sel_rms = np.logical_and(xx > xxmin, xx <= xxmax)
  ngrid = np.shape(xx[sel_rms])[0]
  sims = np.arange(0,ngrid) + 1
  
  color_map = cm.viridis
  #color_map = cm.winter
  ncmap = np.shape(color_map.colors)[0]
  csel = np.random.choice(ncmap, size=nref, replace=False)
  xvcolors = [color_map.colors[csel[ic]] for ic in range(0, nref)]
  
  fig = plt.figure(figsize=(12,12))
  fig.subplots_adjust(hspace=0.05, wspace=0.05)
  
  h_yes = 1 # set to 1 if the first plot is the histogram
  ax = plt.subplot2grid((npar+h_yes , 1), (0,0))
  ax.ticklabel_format(useOffset=False)
  
  ax.axvline(err_T0_d, color='black', ls='-', lw=0.77, alpha=0.99, zorder=3, label='rms threshold = %.3f min' %(err_T0_d*1440.))
  ax = insert_rms_lines(ax, rms_ref_m, rms_max_in, xvcolors)
  
  #kbins = get_bins(xx, rule='stu')
  kbins = get_bins(xx[sel_rms], rule='sq')
  anc.print_both(' | plotting rms_ttv histogram witn nbins = %d' %(kbins), of_log)

  #ax.hist(xx, bins=kbins, histtype='stepfilled', color='#7db0d3', edgecolor='lightgray', align='mid', orientation='vertical', normed=True, stacked=True)
  ax.hist(xx[sel_rms], bins=kbins, range=[xxmin, xxmax], histtype='stepfilled', color='#7db0d3', edgecolor='lightgray', align='mid', orientation='vertical')
  
  ax.set_xlim([xxlim1, xxlim2])
  ax.get_xaxis().set_visible(False)

  ax.get_yaxis().set_visible(True)
  #ylabel = 'Normalised Counts'
  ylabel = 'Counts'
  ax.set_ylabel(ylabel)
  ax.get_yaxis().set_label_coords(-0.075, 0.5)
  
  ax.legend(loc='upper left', fontsize=9, bbox_to_anchor=(1, -2.), ncol=1) # one to legend them all
  plt.draw()

  for iy in range(0,npar):
    yy = final_summary[:,iy]
    yya = final_summary_above[:,iy]
    yylim1, yylim2 = compute_limits(yy, delta=0.5)
    
    anc.print_both(' | plotting rms_ttv vs %s%d' %(p_names[iy], body_id+1), of_log)
  
    ax = plt.subplot2grid((npar+h_yes , 1), (iy+h_yes,0))
    ax.ticklabel_format(useOffset=False)
    
    ax.axvline(err_T0_d, color='black', ls='-', lw=0.77, alpha=0.99, zorder=3)
    ax = insert_rms_lines(ax, rms_ref_m, rms_max_in, xvcolors)
    
    #ax.scatter(xx, yy, color='black', marker='o', s=20, lw=0, alpha=0.05, zorder=5)
    ax.scatter(xx[sel_rms], yy[sel_rms], c=sims, cmap=color_map, marker='o', s=50, lw=0, alpha=0.1, zorder=8)
    #ax.scatter(xxa, yya, c=sims[above], cmap=color_map, marker='o', s=20, lw=0, alpha=0.1, zorder=8)
    
    ax.set_xlim([xxlim1, xxlim2])
    ax.get_xaxis().set_visible(False)
    if(iy == npar-1):
      ax.get_xaxis().set_visible(True)
      xlabel = 'rms(O-C) [days]'
      ax.set_xlabel(xlabel)

    ax.get_yaxis().set_visible(True)
    ylabel = pl_lab[iy]
    ax.set_ylabel(ylabel)
    #ax.yaxis.labelpad = 20
    ax.get_yaxis().set_label_coords(-0.075, 0.5)
    
    
    plt.draw()
  
  plt.draw()
  
  fig_file = os.path.join(out_folder, 'plot_final_summary')
  fig.savefig('%s.pdf' %(fig_file), bbox_inches='tight', dpi=300)
  plt.close(fig)

  return fig_file

# ==============================================================================

# -------------------------------------

def worker(igrid_queue, results_queue):
  
  while igrid_queue.qsize() != 0:
    init_sim = igrid_queue.get()
    igrid = init_sim['igrid']
    transit_id = init_sim['transit_id']
    perturber_par = init_sim['perturber_par']
    perturber_id = init_sim['perturber_id']
    n_bodies = init_sim['n_bodies']
    sel_time = init_sim['sel_time']
    err_T0    = init_sim['err_T0']
    nT0_sel   = init_sim['nT0_sel']
    names_lst = init_sim['names_lst']
    units     = init_sim['units']
    mr_type   = init_sim['mr_type']
    nMC       = init_sim['nMC']
    
    print mp.current_process().name, ' -> igrid = ',igrid,' doing ... (id pert = ',perturber_id,')'
    results_sim = simulation_and_ttv_analysis(igrid, transit_id, perturber_id, n_bodies, perturber_par, sel_time[0], sel_time[1], err_T0, nT0_sel, names_lst, units, mr_type, nMC=nMC)
    #results_sim = {'sim_num': igrid+1}
    print mp.current_process().name, ' -> igrid = ',igrid, ' created results_sim'
    
    results_queue.put(results_sim)
    print mp.current_process().name, ' -> igrid = ',igrid, ' put in results_queue'
    
  return results_queue
  #return

# -------------------------------------

# Define the function for checking the process status
def status(proc):
  
  #print proc.is_alive(), type(proc.is_alive())
  if (proc.is_alive()==True):
    return 'alive'
  elif (proc.is_alive()==False):
    return 'dead'
  else:
    return proc.is_alive()

# -------------------------------------

def get_results_sim_and_save(results_queue, rms_grid, above_threshold, nT0_sel, err_T0, f_h5, mr_type, of_log=None):

  nq_size = results_queue.qsize()
  for i_q in range(0, nq_size):
    results_sim = results_queue.get()
    f_h5 = save_results_sim(f_h5, results_sim, mr_type)
    get_rms_from_results(results_sim, rms_grid, above_threshold, nT0_sel, err_T0)
    anc.print_both(' - Saved sim %06d into hdf5 file.' %(results_sim['sim_num']), of_log)
  
  #f_h5 = save_rms_analysis(f_h5, rms_grid, above_threshold, nT0_sel)
  
  return f_h5

# -------------------------------------


# ==============================================================================
# ==============================================================================
# MAIN

def main():
  
  cli = get_args()
  
  np.random.seed(cli.seed)
  
  
  out_folder, log_file, of_log = set_analysis_folder_log(cli.full_path, cli.sub_folder)
  
  anc.print_both('# THIS IS: %s' %(log_file), of_log)
  anc.print_both('Main folder: %s' %(cli.full_path), of_log)
  
  err_T0 = set_err_T0(cli.err_T0)
  anc.print_both('Mean T0 uncertainty = %.6f days' %(err_T0), of_log)
  
  t_launch = atime.Time('2018-06-30T06:00:00.0', format='isot', scale='tdb')
  t_end = t_launch.jd + 365.25*3.5
  
  anc.print_both('CHEOPS DATE: LAUNCH = %.4f [BJD] END = %.4f [BJD]' %(t_launch.jd, t_end), of_log)
  anc.print_both(' ', of_log)
  
  anc.print_both('NUMBER OF CPUS = %d' %(cli.nthreads), of_log)
  anc.print_both(' ', of_log)
  
  # init TRADES
  anc.print_both(' - INIT TRADES', of_log)
  pytrades.initialize_trades(cli.full_path, cli.sub_folder, 1)
  n_bodies = pytrades.n_bodies # NUMBER OF TOTAL BODIES OF THE SYSTEM
  n_planets = n_bodies - 1 # NUMBER OF PLANETS IN THE SYSTEM
  
  # read grid parameters ==> output ngrid and ncol
  ngrid, ncol = pytrades.wrapper_read_grid(1)
  # create grid
  perturber_grid = pytrades.wrapper_grid_init(1, ngrid, ncol)
  anc.print_both(' - SET PERTURBER GRID', of_log)
  #perturber_grid = pytrades.wrapper_grid_init(1)
  anc.print_both(' | To run %d simulation/s.\n' %(ngrid), of_log)
  
  #mr_convert, mr_unit = mass_factor(cli.mr_type)
  hdf5_name, f_h5, mr_convert, mr_unit, names_lst, names_str, units = create_hdf5_file(out_folder, cli.seed, cli.perturber_id, perturber_grid, cli.mr_type, t_launch.jd, t_end, err_T0, of_log)
  anc.print_both(' - created hdf5 file: %s' %(hdf5_name), of_log)
  
  # define the number of transits to select:
  nT0_sel = [3, 5, 6, 10, 15]
  # define the number of Monte-Carlo repeatition for each nT0_sel
  nMC = 1000
  n_rms = 1 + len(nT0_sel)*2 # rms_ttv + nT0_sel[mean, median]
  
  transit_id = 2
  
  #niter = ngrid
  niter = 4
  
  igrid_queue = mp.Queue()
  results_queue = mp.Queue()
  for igrid in range(0, niter):
    init_sim = {'igrid':igrid, 
                 'perturber_par':perturber_grid[igrid, :],
                 'transit_id': transit_id,
                 'perturber_id': cli.perturber_id,
                 'n_bodies': n_bodies,
                 'sel_time':[t_launch.jd, t_end],
                 'err_T0': err_T0, 
                 'nT0_sel': nT0_sel, 
                 'names_lst': names_lst, 
                 'units': units, 
                 'mr_type': cli.mr_type, 
                 'nMC': nMC
                 }
    igrid_queue.put(init_sim)
  
  procs = []
  for i_p in range(cli.nthreads):
    pp = mp.Process(target=worker, args=(igrid_queue, results_queue))
    pp.start()
    procs.append(pp)
    
  for pp in procs:
    pp.join()
    
  for pp in procs:
    print "Process ",pp, " @ ", pp.pid, " is ", status(pp)
  while len(procs) != 0:
    for pp in procs:
      if status(pp) == "dead":
        procs.pop(procs.index(pp))
    # loose 2 seconds before checking again
    time.sleep(2)
  
  anc.print_both(' - Completed all the simulations. Now save simulation analysis.', of_log)
  
  rms_grid = np.zeros((niter, n_rms))
  above_threshold = np.zeros((niter, n_rms)).astype(bool)
  f_h5 = get_results_sim_and_save(results_queue, rms_grid, above_threshold, nT0_sel, err_T0, f_h5, cli.mr_type, of_log)
  f_h5 = save_rms_analysis(f_h5, rms_grid, above_threshold, nT0_sel)  
  f_h5.close()
  
  anc.print_both(' ', of_log)
  anc.print_both(' ===============================', of_log)
  anc.print_both(' - Saved rms_grid into hdf5 file', of_log)
  anc.print_both(' ', of_log)
  
  ## create the plot:
  #fig_file = plot_summary(out_folder, hdf5_name, body_id=2, rms_id_in=npar, of_log=of_log)

  of_log.close()
  
  return

# ==============================================================================
# RUN MAIN

if __name__ == "__main__":
  main()

# ==============================================================================

