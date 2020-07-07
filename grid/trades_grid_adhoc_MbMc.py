#!/usr/bin/env python
# -*- coding: utf-8 -*-

# ==============================================================================
# IMPORT

 # no more "zero" integer division bugs!:P
import argparse
import numpy as np # array
#import h5py
#import astropy.time as atime
import os
import sys
import h5py

import emcee.interruptible_pool as emcee_pool
#import multiprocessing as mp
#import gc
# ==============================================================================
import matplotlib as mpl
#from matplotlib import use as mpluse
mpl.use("Agg", warn=False)
#mpluse("Qt4Agg")
#mpluse("TkAgg")
import matplotlib.pyplot as plt
plt.rcParams['text.usetex'] = True
plt.rcParams['font.family'] = 'serif'
#import matplotlib.mlab as mlab
import matplotlib.cm as cm
# ==============================================================================
import gls
# ==============================================================================
# custom modules
script_path = os.path.realpath(__file__)
module_path = os.path.abspath(os.path.join(os.path.dirname(script_path),
                                           '../pytrades'))
sys.path.append(module_path)

import ancillary as anc
import constants as cst
from pytrades_lib import pytrades
# ==============================================================================

def set_analysis_folder_log(full_path, sub_folder):
  
  out_folder = os.path.join(full_path, sub_folder)
  if(not os.path.isdir(out_folder)): os.makedirs(out_folder)
  
  anc.copy_simulation_files(full_path, out_folder)
  
  log_file = os.path.join(out_folder, 'analysis_log.txt')
  of_log = open(log_file, 'w')
  
  return out_folder, log_file, of_log

# ==============================================================================

def set_err_T0(err_T0_in):
  
  if('d' in err_T0_in):
    err_T0 = np.float64(err_T0_in.split('d')[0])
  elif('m' in err_T0_in):
    err_T0 = np.float64(err_T0_in.split('m')[0])/1440.
  elif('s' in err_T0_in):
    err_T0 = np.float64(err_T0_in.split('s')[0])/86400.
  else:
    err_T0 = 10./86400.
    
  return err_T0

# ==============================================================================

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

# ==============================================================================

def run_simulations(sim_args):
  
  sim_id = sim_args['sim_id']
  mb_s, mc_s = sim_args['m_grid']
  mass_0, radius, period, sma, ecc, argp, meanA, inc, longN = sim_args['kep_elem_0']
  mass = mass_0.copy()
  mass[1] = mb_s
  mass[2] = mc_s
  n_bodies = sim_args['n_bodies']
  nt0_full = sim_args['nt0_full']
  nrv_nmax = sim_args['nrv_nmax']
  t_start, t_end = sim_args['time_sel']
  err_T0 = sim_args['err_T0']
  
  # set default output_sim:
  output_sim = {'sim_id': sim_id,
                'Tb_ephem': 0., 'Pb_ephem': 0.,
                'T0b_cheops': None, 'ocb_cheops': None,
                'sigma_ocb': 0., 'rms_ocb': 0., 'amp_b': 0.,
                'PTTVb': 0., 'FAPb': 0.,
                'Tc_ephem': 0., 'Pc_ephem': 0.,
                'T0c_cheops': None, 'occ_cheops': None,
                'sigma_occ': 0., 'rms_occ': 0., 'amp_c': 0.,
                'PTTVc': 0., 'FAPc': 0.
                }
  
  pid = os.getpid()
  print(' Process %d RUNNING sim %d ...' %(pid, sim_id))
  
  #ttra_full, id_ttra_full, stats_ttra, time_rv_nmax, rv_nmax, stats_rv = fit_par_to_ttra_rv(fit_parameters)
  ttra_full, id_ttra_full, stats_ttra, time_rv_nmax, rv_nmax, stats_rv = pytrades.wrapper_run_grid_combination(mass, radius, period, sma, ecc, argp, meanA, inc, longN, nt0_full, nrv_nmax)
  
  # check if transits of planet b exist
  T0b_all = np.sort(ttra_full[np.logical_and(id_ttra_full==2, stats_ttra)])
  nT0b_all = np.shape(T0b_all)[0]
  # check if transits of planet c exist
  T0c_all = np.sort(ttra_full[np.logical_and(id_ttra_full==3, stats_ttra)])
  nT0c_all = np.shape(T0c_all)[0]
  
  if(nT0b_all == 0 or nT0c_all == 0):
    return output_sim

  # select transits within CHEOPS mission
  
  T0b = T0b_all[np.logical_and(T0b_all >= t_start, T0b_all <= t_end)]
  nT0b = np.shape(T0b)[0]
  
  T0c = T0c_all[np.logical_and(T0c_all >= t_start, T0c_all <= t_end)]
  nT0c = np.shape(T0c)[0]
  
  if(nT0b == 0 or nT0c_all == 0):
    return output_sim
  
  # fit linear ephemeris, computes O-C, std(O-C), rms(O-C)
  eT0b = np.zeros((nT0b))+err_T0
  epob, Tb_ephem, Pb_ephem = compute_lin_ephem(T0b, eT0b)
  ocb_d = T0b - (Tb_ephem + epob * Pb_ephem)
  sigma_ocb = np.std(ocb_d, ddof=1)
  rms_ocb = np.sqrt(np.mean(ocb_d*ocb_d))
  
  eT0c = np.zeros((nT0c))+err_T0
  epoc, Tc_ephem, Pc_ephem = compute_lin_ephem(T0c, eT0c)
  occ_d = T0c - (Tc_ephem + epoc * Pc_ephem)
  sigma_occ = np.std(occ_d, ddof=1)
  rms_occ = np.sqrt(np.mean(occ_d*occ_d))
  
  # computes GLS on the O-C
  time_span_b = T0b[-1] - T0b[0]
  #glsb_obj = gls.Gls((T0b, ocb_d, eT0b), Pbeg=3.*Pb_ephem, Pend=2.*time_span_b, ofac=10, verbose=False)
  glsb_obj = gls.Gls((T0b, ocb_d, eT0b), Pbeg=3.*Pb_ephem, Pend=4.*time_span_b, ofac=10, verbose=False)
  max_pwb, fbestb = glsb_obj.pmax, glsb_obj.hpstat['fbest']
  PTTVb, FAPb, amp_b = 1./fbestb, glsb_obj.FAP(max_pwb), glsb_obj.hpstat['amp']
  # comment!!
  #glsb_obj.plot(block=False)
  
  time_span_c = T0c[-1] - T0c[0]
  #glsc_obj = gls.Gls((T0c, occ_d, eT0c), Pbeg=3.*Pc_ephem, Pend=2.*time_span_c, ofac=10, verbose=False)
  glsc_obj = gls.Gls((T0c, occ_d, eT0c), Pbeg=3.*Pc_ephem, Pend=4.*time_span_c, ofac=10, verbose=False)
  max_pwc, fbestc = glsc_obj.pmax, glsc_obj.hpstat['fbest']
  PTTVc, FAPc, amp_c = 1./fbestc, glsc_obj.FAP(max_pwc), glsc_obj.hpstat['amp']
  # comment!!
  #glsc_obj.plot(block=False)

  output_sim['Tb_ephem'] = Tb_ephem
  output_sim['Pb_ephem'] = Pb_ephem
  output_sim['T0b_cheops'] = T0b
  output_sim['ocb_cheops'] = ocb_d
  output_sim['sigma_ocb'] = sigma_ocb
  output_sim['rms_ocb'] = rms_ocb
  output_sim['amp_b'] = amp_b
  output_sim['PTTVb'] = PTTVb
  output_sim['FAPb'] = FAPb
  output_sim['Tc_ephem'] = Tc_ephem
  output_sim['Pc_ephem'] = Pc_ephem
  output_sim['T0c_cheops'] = T0c
  output_sim['occ_cheops'] = occ_d
  output_sim['sigma_occ'] = sigma_occ
  output_sim['rms_occ'] = rms_occ
  output_sim['amp_c'] = amp_c
  output_sim['PTTVc'] = PTTVc
  output_sim['FAPc'] = FAPc
  
  print(' Process %d COMPLETED sim %d' %(pid, sim_id))
  
  return output_sim

# ==============================================================================
def plot_meshgrid(Mb_mesh, Mc_mesh, amp_b, amp_c, PTTV_b, PTTV_c, out_folder):
  
  #amp_min = min(np.min(amp_b), np.min(amp_c))
  amp_min = 0.
  amp_max = max(np.max(amp_b), np.max(amp_c))
  lev_amp = np.arange(0., amp_max+1., 1.)
  
  plot_folder = os.path.join(out_folder, 'plots')
  if(not os.path.isdir(plot_folder)): os.makedirs(plot_folder)
  plot_file = os.path.join(plot_folder, 'MbMc_grid')
  
  fig = plt.figure(figsize=(12,12))
  fig.subplots_adjust(hspace=0.13, wspace=0.11)
  
  # Mb vs Mc vs amp_b
  ax = plt.subplot2grid((2, 2), (0,0))
  #scax = ax.pcolormesh(Mb_mesh, Mc_mesh, amp_b,
  scax = ax.contourf(Mb_mesh, Mc_mesh, amp_b,
                    levels=lev_amp,
                    cmap=cm.viridis,
                    vmin=amp_min, vmax=amp_max,
                    zorder=4)

  # colorbar
  cbar = plt.colorbar(scax, ax=ax, orientation='vertical', pad=0.01, format='%4.1f')
  cbar.ax.set_ylabel('$\sigma_\\textrm{TTV} (\\textrm{min})$', size=7)
  cbar.ax.tick_params(labelsize=6)
  
  ax.tick_params(direction='in')
  #ax.set_xticklabels([])
  ax.set_xlabel('$M_\\textrm{b} (M_\oplus)$')
  ax.set_ylabel('$M_\\textrm{c} (M_\oplus)$')
  ax.set_title('planet b')
  
  # Mb vs Mc vs amp_c
  ax = plt.subplot2grid((2, 2), (0,1))
  #scax = ax.pcolormesh(Mb_mesh, Mc_mesh, amp_c,
  scax = ax.contourf(Mb_mesh, Mc_mesh, amp_c,
                    levels=lev_amp,
                    cmap=cm.viridis,
                    vmin=amp_min, vmax=amp_max,
                    zorder=4)
  
  ax.tick_params(direction='in')
  #ax.set_xticklabels([])
  #ax.set_yticklabels([])
  ax.set_xlabel('$M_\\textrm{b} (M_\oplus)$')
  ax.set_ylabel('$M_\\textrm{c} (M_\oplus)$')
  ax.set_title('planet c')
  
  # colorbar
  cbar = plt.colorbar(scax, ax=ax, orientation='vertical', pad=0.01, format='%4.1f')
  cbar.ax.set_ylabel('$\sigma_\\textrm{TTV} (\\textrm{min})$', size=7)
  cbar.ax.tick_params(labelsize=6)

  # Mb vs Mc vs PTTV_b
  lev_pb = np.linspace(np.min(PTTV_b), np.max(PTTV_b), num=10, endpoint=True)
  ax = plt.subplot2grid((2, 2), (1,0))
  #scax = ax.pcolormesh(Mb_mesh, Mc_mesh, PTTV_b,
  scax = ax.contourf(Mb_mesh, Mc_mesh, PTTV_b,
                    levels=lev_pb,
                    cmap=cm.winter,
                    #vmin=PTTV_min, vmax=PTTV_max,
                    zorder=4)

  ax.tick_params(direction='in')
  ax.set_xlabel('$M_\\textrm{b} (M_\oplus)$')
  ax.set_ylabel('$M_\\textrm{c} (M_\oplus)$')
  ax.set_title('planet b')

  # colorbar
  cbar = plt.colorbar(scax, ax=ax, orientation='vertical', pad=0.01, format='%4.1f')
  cbar.ax.set_ylabel('$P_\\textrm{TTV} (\\textrm{days})$', size=7)
  cbar.ax.tick_params(labelsize=6)
  
  # Mb vs Mc vs PTTV_c
  lev_pc = np.linspace(np.min(PTTV_c), np.max(PTTV_c), num=10, endpoint=True)
  ax = plt.subplot2grid((2, 2), (1,1))
  #scax = ax.pcolormesh(Mb_mesh, Mc_mesh, PTTV_c,
  scax = ax.contourf(Mb_mesh, Mc_mesh, PTTV_c,
                    levels=lev_pc,
                    cmap=cm.autumn,
                    #vmin=PTTV_min, vmax=PTTV_max,
                    zorder=4)
  
  ax.tick_params(direction='in')
  #ax.set_yticklabels([])
  ax.set_title('planet c')
  
  # colorbar
  cbar = plt.colorbar(scax, ax=ax, orientation='vertical', pad=0.01, format='%4.1f')
  cbar.ax.set_ylabel('$P_\\textrm{TTV} (\\textrm{days})$', size=7)
  cbar.ax.tick_params(labelsize=6)

  ax.set_xlabel('$M_\\textrm{b} (M_\oplus)$')
  ax.set_ylabel('$M_\\textrm{c} (M_\oplus)$')
  
  fig.savefig('%s.png' %(plot_file), bbox_inches='tight', dpi=300)
  plt.close(fig)

  return

# ==============================================================================

def read_grid_hdf5(h5_file):
  
  h5f = h5py.File(h5_file, 'r')
  Mb_e = h5f['Mb_e'][...]
  Mb_efull = h5f['Mb_efull'][...]
  Mb_mesh = h5f['Mb_mesh'][...]
  Mc_e = h5f['Mc_e'][...]
  Mc_efull = h5f['Mc_efull'][...]
  Mc_mesh = h5f['Mc_mesh'][...]
  
  amp_b = h5f['amp_b'][...]
  amp_b_full = h5f['amp_b_full'][...]
  amp_c = h5f['amp_c'][...]
  amp_c_full = h5f['amp_c_full'][...]
  
  PTTV_b = h5f['PTTV_b'][...]
  PTTV_b_full = h5f['PTTV_b_full'][...]
  PTTV_c = h5f['PTTV_c'][...]
  PTTV_c_full = h5f['PTTV_c_full'][...]
  h5f.close()

  return Mb_e, Mb_efull, Mb_mesh, Mc_e, Mc_efull, Mc_mesh, amp_b, amp_b_full, amp_c, amp_c_full, PTTV_b, PTTV_b_full, PTTV_c, PTTV_c_full
  
# ==============================================================================
# read cli arguments
def get_args():
  
  parser = argparse.ArgumentParser(description='CHEOPS TTV LEVEL ANALYSIS')
  parser.add_argument('-p', '--path', action='store', dest='full_path', required=True, help='The path (absolute or relative) with simulation files for TRADES.')
  
  #parser.add_argument('-lm', '--lm-flag', action='store', dest='lm_flag', default='0', help='LM flag. If you want to read files w/o LM set to 0, otherwise 1. Default is 0')
  
  # SUBFOLDER TO SAVE ALL SIMULATION DATA
  parser.add_argument('-sf', '--sub-folder', '--sub-f', action='store', dest='sub_folder', default='grid2ttv', help='Sub-folder name, without full path. Default = grid2ttv')
  
  parser.add_argument('-e', '--e-T0', '--err-T0', action='store', dest='err_T0', default='0.', help='Mean uncertainty on the Transit Times (T0s). Provide value with unit "d" "m" "s" for days, minutes, and seconds in this way: "0.001528d" or "2.2m" or "132.s". Default unit is days and value set to 0.')
  
  parser.add_argument('-c', '--cpu', '--threads', action='store', dest='nthreads', default=1, type=int, help='Number of threads/cpu to use. Default is 1.')
    
  cli = parser.parse_args()
  cli.full_path = os.path.join(os.path.abspath(cli.full_path), '')
  cli.sub_folder = os.path.join(os.path.relpath(cli.sub_folder), '')
    
  return cli

# ==============================================================================
def main():
  
  cli = get_args()
  
  out_folder, log_file, of_log = set_analysis_folder_log(cli.full_path, cli.sub_folder)
  err_T0 = set_err_T0(cli.err_T0)
  
  ncpu = int(cli.nthreads)
  # init TRADES
  pytrades.initialize_trades(cli.full_path, cli.sub_folder, ncpu)
  all_par_0 = pytrades.system_parameters
  fit_par_0 = pytrades.fitting_parameters
  n_bodies = pytrades.n_bodies # NUMBER OF TOTAL BODIES OF THE SYSTEM
  n_planets = n_bodies - 1 # NUMBER OF PLANETS IN THE SYSTEM
  MRs = pytrades.mr_star
  
  # hip41378
  #Mb_e = np.linspace(0.1, 21.4, num=100, endpoint=True) # Mearth
  #Mc_e = np.linspace(0.1, 27.15, num=100, endpoint=True) # Mearth
  
  # K2-24
  Mb_e = np.linspace(7.30, 39.1, num=100, endpoint=True) # Mearth
  Mc_e = np.linspace(11.8, 50.2, num=100, endpoint=True) # Mearth
  
  nMb = np.shape(Mb_e)[0]
  nMc = np.shape(Mc_e)[0]
  nsims = nMb * nMc
  
  Mb_efull = np.repeat(Mb_e, nMc)
  Mc_efull = np.tile(Mc_e, nMb)

  Mb_s = Mb_efull * cst.Mears
  Mc_s = Mc_efull * cst.Mears
  #MbMs = Mb_efull * cst.Mears / MRs[0,0]
  #McMs = Mc_efull * cst.Mears / MRs[0,0]
  #fit_par_grid = [[MbMs[isim], McMs[isim]] for isim in range(nsims)]
  
  mass, radius, period, sma, ecc, argp, meanA, inc, longN = pytrades.convert_trades_par_to_kepelem(all_par_0, fit_par_0, n_bodies)
  nt0_full, nrv_nmax = pytrades.get_max_nt0_nrv(period)
  kep_elem_0 = [mass, radius, period, sma, ecc, argp, meanA, inc, longN]
  
  #t_launch = atime.Time('2018-06-30T06:00:00.0', format='isot', scale='tdb')
  #t_end = t_launch.jd + 365.25*3.5
  t_launch = 2.*period[2]
  t_end = t_launch + 365.25*3.5
    
  list_sim_args = []
  for isim in range(0, nsims):
    sim_args = {'sim_id': isim+1,
                'm_grid': [Mb_s[isim], Mc_s[isim]], 
                'kep_elem_0': kep_elem_0,
                'n_bodies': n_bodies,
                'nt0_full': nt0_full,
                'nrv_nmax': nrv_nmax,
                'time_sel': [t_launch, t_end],
                'err_T0': err_T0
                }
    list_sim_args.append(sim_args)
    
    #output_sim = run_simulations(sim_args)
    #print
    #print 'sim: iter+1 = %d' %(isim+1)
    #print 'Mb = %12.4f Mearth || Mc = %12.4f Mearth' %(Mb_efull[isim], Mc_efull[isim])
    #print 'Planet b sigma(O-C) = %12.6f d rms(O-C) = %12.6f d amp(O-C) = %12.6f d PTTV = %12.4f d FAP = %23.16e' %(output_sim['sigma_ocb'], output_sim['rms_ocb'], output_sim['amp_b'], output_sim['PTTVb'], output_sim['FAPb'])
    #print 'Planet c sigma(O-C) = %12.6f d rms(O-C) = %12.6f d amp(O-C) = %12.6f d PTTV = %12.6f d FAP = %23.16e' %(output_sim['sigma_occ'], output_sim['rms_occ'], output_sim['amp_c'], output_sim['PTTVc'], output_sim['FAPc'])

  #print
  #raw_input('HIT ENTER TO CLOSE EVERYTHING!')

  # create pool
  threads_pool = emcee_pool.InterruptiblePool(ncpu)
  #threads_pool = mp.Pool(ncpu)
  # run simulations with the map
  list_output_sims = threads_pool.map(run_simulations, list_sim_args)
  # close the pool of threads
  threads_pool.close()
  threads_pool.terminate()
  threads_pool.join()
  
  #print 'len(list_output_sims) = %d' %(len(list_output_sims))
  #for isim in range(0, nsims):
    #print 'sim: iter+1 = %d' %(isim+1)
    #print 'Mb = %12.4f Mearth || Mc = %12.4f Mearth' %(Mb_efull[isim], Mc_efull[isim])
    #print 'Planet b sigma(O-C) = %12.6f d rms(O-C) = %12.6f d amp(O-C) = %12.6f d PTTV = %12.4f d FAP = %23.16e' %(list_output_sims[isim]['sigma_ocb'], list_output_sims[isim]['rms_ocb'], list_output_sims[isim]['amp_b'], list_output_sims[isim]['PTTVb'], list_output_sims[isim]['FAPb'])
    #print 'Planet c sigma(O-C) = %12.6f d rms(O-C) = %12.6f d amp(O-C) = %12.6f d PTTV = %12.6f d FAP = %23.16e' %(list_output_sims[isim]['sigma_occ'], list_output_sims[isim]['rms_occ'], list_output_sims[isim]['amp_c'], list_output_sims[isim]['PTTVc'], list_output_sims[isim]['FAPc'])
  
  Mb_mesh = np.reshape(Mb_efull, (nMc, nMb))
  Mc_mesh = np.reshape(Mc_efull, (nMc, nMb))
  
  amp_b_full = np.array([list_output_sims[isim]['amp_b'] for isim in range(0,nsims)]) * 1440.
  amp_b = np.reshape([list_output_sims[isim]['amp_b'] for isim in range(0,nsims)], (nMc, nMb)) * 1440.
  amp_c_full = np.array([list_output_sims[isim]['amp_c'] for isim in range(0,nsims)]) * 1440.
  amp_c = np.reshape([list_output_sims[isim]['amp_c'] for isim in range(0,nsims)], (nMc, nMb)) * 1440.
  
  PTTV_b = np.reshape([list_output_sims[isim]['PTTVb'] for isim in range(0,nsims)], (nMc, nMb))
  PTTV_b_full = np.array([list_output_sims[isim]['PTTVb'] for isim in range(0,nsims)])
  PTTV_c = np.reshape([list_output_sims[isim]['PTTVc'] for isim in range(0,nsims)], (nMc, nMb))
  PTTV_c_full = np.array([list_output_sims[isim]['PTTVc'] for isim in range(0,nsims)])
  
  h5_file = os.path.join(out_folder, 'MbMc_grid.hdf5')
  h5f = h5py.File(h5_file, 'w')
  h5f.create_dataset('Mb_e', data=Mb_e, dtype=np.float64, compression='gzip')
  h5f.create_dataset('Mb_efull', data=Mb_e, dtype=np.float64, compression='gzip')
  h5f.create_dataset('Mb_mesh', data=Mb_mesh, dtype=np.float64, compression='gzip')
  h5f.create_dataset('Mc_e', data=Mc_e, dtype=np.float64, compression='gzip')
  h5f.create_dataset('Mc_efull', data=Mb_e, dtype=np.float64, compression='gzip')
  h5f.create_dataset('Mc_mesh', data=Mc_mesh, dtype=np.float64, compression='gzip')
  
  h5f.create_dataset('amp_b', data=amp_b, dtype=np.float64, compression='gzip')
  h5f.create_dataset('amp_b_full', data=amp_b_full, dtype=np.float64, compression='gzip')
  h5f.create_dataset('amp_c', data=amp_c, dtype=np.float64, compression='gzip')
  h5f.create_dataset('amp_c_full', data=amp_c_full, dtype=np.float64, compression='gzip')
  
  h5f.create_dataset('PTTV_b', data=PTTV_b, dtype=np.float64, compression='gzip')
  h5f.create_dataset('PTTV_b_full', data=PTTV_b_full, dtype=np.float64, compression='gzip')
  h5f.create_dataset('PTTV_c', data=PTTV_c, dtype=np.float64, compression='gzip')
  h5f.create_dataset('PTTV_c_full', data=PTTV_c_full, dtype=np.float64, compression='gzip')
  h5f.close()
  
  plot_meshgrid(Mb_mesh, Mc_mesh, amp_b, amp_c, PTTV_b, PTTV_c, out_folder)
  
  pytrades.deallocate_variables()
  
  return
# ==============================================================================
# ==============================================================================
if __name__ == "__main__":
  main()



