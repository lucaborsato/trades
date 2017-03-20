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
import sys
import os
import argparse
import numpy as np # array
import h5py
import astropy.time as atime

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

# ==============================================================================

# masses conversions
Msun2mer = 6.0236e6 # Msun to Mmer
Mmer2sun = 1./Msun2mer # Mmer to Msun
Msun2ven = 4.08523719e5 # Msun to Mven
Mven2sun = 1./Msun2ven   #  Mven to Msun
Msun2ear = 332946.0487 # Msun to Mear
Mear2sun = 1./Msun2ear  #  Mear to Msun
Msun2mar = 3.09870359e6 # Msun to Mmar
Mmar2sun = 1./Msun2mar   #  Mmar to Msun
Msun2jup = 1.047348644e3 # Msun to Mjup
Mjup2sun = 1./Msun2jup    #  Mjup to Msun
Msun2sat = 3.4979018e3 # Msun to Msat
Msat2sun = 1./Msun2sat  #  Msat to Msun
Msun2ura = 2.290298e4 # Msun to Mura
Mura2sun = 1./Msun2ura #  Mura to Msun
Msun2nep = 1.941226e4 # Msun to Mnep
Mnep2sun = 1./Msun2nep #  Mnep to Msun

# ==============================================================================
# read cli arguments
def get_args():
  
  parser = argparse.ArgumentParser(description='CHEOPS TTV LEVEL ANALYSIS')
  parser.add_argument('-p', '--path', action='store', dest='full_path', required=True, help='The path (absolute or relative) with simulation files for TRADES.')
  parser.add_argument('-lm', '--lm-flag', action='store', dest='lm_flag', default='0', help='LM flag. If you want to read files w/o LM set to 0, otherwise 1. Default is 0')
  parser.add_argument('-e', '--e-T0', '--err-T0', action='store', dest='err_T0', default='0.', help='Mean uncertainty on the Transit Times (T0s). Provide value with unit "d" "m" "s" for days, minutes, and seconds in this way: "0.001528d" or "2.2m" or "132.s". Default unit is days and value set to 0.')
  parser.add_argument('-m', '--m-type', '--mass-type', action='store', dest='mass_type', default='e', help='Mass type to convert the grid. "e ea ear eart earth" or "j ju jup jupi jupit jupite jupiter" or "n ne nep nept neptu neptun neptune".')
  parser.add_argument('-s', '--seed', action='store', dest='seed', default='None', help='Seed for random number generator. Defauto is None.')
  
  cli = parser.parse_args()
  cli.full_path = os.path.abspath(cli.full_path)
  cli.lm_flag = str(cli.lm_flag)
  try:
    cli.seed = int(cli.seed)
  except:
    cli.seed = None
  
  return cli

# -------------------------------------

def print_both(line, of_log=None):
  
  print line
  if(of_log is not None):
    of_log.write(line + '\n')
  
  return

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

def set_analysis_folder_log(full_path):
  
  out_folder = os.path.join(full_path, 'level_analysis')
  if(not os.path.isdir(out_folder)): os.makedirs(out_folder)
  
  log_file = os.path.join(out_folder, 'analysis_log.txt')
  of_log = open(log_file, 'w')
  
  return out_folder, log_file, of_log

# -------------------------------------

def mass_factor(mass_type):
  
  if(mass_type.lower() in 'j ju jup jupi jupit jupite jupiter'.split()):
    mass_convert = Msun2jup
    mass_unit = 'M_jup'
  elif(mass_type.lower() in 'n ne nep nept neptu neptun neptune'.split()):
    mass_convert = Msun2nep
    mass_unit = 'M_nep'
  else:
    mass_convert = Msun2ear
    mass_unit = 'M_ear'
  
  return mass_convert, mass_unit
    
# -------------------------------------

def read_grid_summary(full_path, lm_flag, of_log=None):
  
  grid_file = os.path.join(full_path, 'summary_grid_sims_%s.dat' %(lm_flag))
  try:
    grid = np.genfromtxt(grid_file, names=True)
    grid_names = grid.dtype.names
    ngrid = np.shape(grid)[0]
  except:
    print_both('%s NOT FOUND' %(grid_file))
    grid, grid_names, ngrid = None, None, 0
  
  return grid, grid_names, ngrid

# -------------------------------------

def print_perturber_parameters(grid, sim_id, id_body, mass_convert, mass_unit, of_log):
  
  body_num = '%d' %(id_body + 1)
  base_names_lst = 'm P e w mA i lN'.split()
  units = '%s day - deg deg deg deg' %(mass_unit)
  units = units.split()
  n_names = len(base_names_lst)
  names_lst = ['%s%s' %(base_names_lst[i], body_num) for i in range(0, n_names)]
  temp_lst = ['%s [%s]' %(names_lst[i], units[i].strip()) for i in range(0, n_names)]
  names_str = ' '.join(['%23s' %(temp_lst[i]) for i in range(0, n_names)])
  
  mass   = grid[names_lst[0]][sim_id]*mass_convert
  period = grid[names_lst[1]][sim_id]
  ecc    = grid[names_lst[2]][sim_id]
  argp   = grid[names_lst[3]][sim_id]
  meanA  = grid[names_lst[4]][sim_id]
  inc    = grid[names_lst[5]][sim_id]
  longN  = grid[names_lst[6]][sim_id]
  
  perturber_parameters = [mass, period, ecc, argp, meanA, inc, longN]
  param_str = ' '.join(['%23.16e' %(perturber_parameters[i]) for i in range(0, n_names)])
  
  print_both(' - Perturber parameters', of_log)
  print_both(' | %s' %(names_str), of_log)
  print_both(' | %s' %(param_str), of_log)

  return names_lst, names_str, units, perturber_parameters

# -------------------------------------

def read_transit_file(full_path, lm_flag, id_sim, id_body, t_min_sel=None, t_max_sel=None, of_log=None):
  
  transit_file = os.path.join(full_path, '%d_%s_NB%d_tra.dat' %(id_sim, lm_flag, id_body))
  try:
    data = np.genfromtxt(transit_file, usecols=(0,1))
    print_both(' | %s' %(transit_file), of_log)
    T0_all = data[:,0] + data[:,1]
    nT0_all = np.shape(T0_all)[0]
    print_both(' | number of transits nT0 = %d' %(nT0_all), of_log)
    
    if(t_min_sel is None):
      sel_min = 0.
    else:
      sel_min = t_min_sel
    if(t_max_sel is None):
      sel_max = T0_all[-1]
    else:
      sel_max = t_max_sel
      
    sel_pre = T0_all < sel_min # T0s before CHEOPS mission, they are use to define Tref, Pref
    nT0_pre = int(np.sum(np.asarray(sel_pre, dtype=np.int)))
    print_both(' | selected pre-CHEOPS transits: nT0 (pre) = %d' %(nT0_pre) , of_log)
    
    sel_cheops = np.logical_and(T0_all >= sel_min, T0_all <= sel_max)
    nT0_cheops = int(np.sum(np.asarray(sel_cheops, dtype=np.int)))
    print_both(' | selected CHEOPS transits: nT0 (cheops) = %d' %(nT0_cheops) , of_log)
  except:
    return None, 0, None, 0, None, 0
  
  
  return T0_all, nT0_all, sel_pre, nT0_pre, sel_cheops, nT0_cheops

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
  
  # with numpy array
  #rms_all = np.zeros(nMC)
  #for imc in nrange(0, nMC):
    #sel_oc = np.random.choice(oc, nT0_sel, replace=False) # select nT0_sel from the oc without replacement
    #rms_all[imc] = rms(sel_oc)
  rms_all = [rms(np.random.choice(oc, nT0_sel, replace=False)) for imc in range(0, nMC)]
  rms_mean = np.mean(rms_all)
  rms_med = np.percentile(rms_all, 50., interpolation='midpoint')
    
  return rms_mean, rms_med

# -------------------------------------

# ==============================================================================

# -------------------------------------

def create_hdf5_file(full_path, seed, grid, grid_names, ngrid, t_sel_min, t_sel_max, err_T0):
  
  new_grid = np.array([grid[grid_names[i_n]][:] for i_n in range(len(grid_names))], dtype=np.float64).T
  
  file_name = os.path.join(full_path, 'rms_analysis_summary.hdf5')
  f_h5 = h5py.File(file_name, 'w')
  
  f_h5.create_dataset('grid_header', data=grid_names, dtype='S20')
  f_h5.create_dataset('grid', data=new_grid, dtype=np.float64)
  f_h5['grid'].attrs['ngrid'] = ngrid
  f_h5['grid'].attrs['seed'] = seed
  f_h5['grid'].attrs['t_sel_min'] = t_sel_min
  f_h5['grid'].attrs['t_sel_max'] = t_sel_max
  f_h5['grid'].attrs['err_T0_day'] = err_T0
  f_h5['grid'].attrs['err_T0_min'] = err_T0*1440.
  
  return file_name, f_h5

# -------------------------------------

def save_results_sim(f_h5, results_sim):
  
  sim_num = results_sim['sim_num']
  f_h5.create_dataset('sim_%06d/parameter_names'     %(sim_num), data=np.array(results_sim['parameter_names'], dtype='S10')     , dtype='S20')
  f_h5.create_dataset('sim_%06d/parameter_units'     %(sim_num), data=np.array(results_sim['parameter_units'], dtype='S10')     , dtype='S20')
  f_h5.create_dataset('sim_%06d/perturber_parameters'%(sim_num), data=results_sim['perturber_parameters'], dtype=np.float64)
  f_h5.create_dataset('sim_%06d/linear_ephem'        %(sim_num), data=results_sim['linear_ephem']        , dtype=np.float64)
  f_h5.create_dataset('sim_%06d/rms_ttv'             %(sim_num), data=results_sim['rms_ttv']             , dtype=np.float64)
  f_h5.create_dataset('sim_%06d/nT0_sel'             %(sim_num), data=results_sim['nT0_sel']             , dtype=np.int)
  f_h5.create_dataset('sim_%06d/rms_mc'              %(sim_num), data=results_sim['rms_mc']              , dtype=np.float64)

  f_h5['sim_%06d' %(sim_num)].attrs['sim_num'] = results_sim['sim_num']
  f_h5['sim_%06d' %(sim_num)].attrs['nMC'] = results_sim['nMC']
  
  return f_h5

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
  print_both(' | plotting rms_ttv histogram witn nbins = %d' %(kbins), of_log)

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
    
    print_both(' | plotting rms_ttv vs %s%d' %(p_names[iy], body_id+1), of_log)
  
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
# MAIN

def main():
  
  cli = get_args()
  
  np.random.seed(cli.seed)
  
  out_folder, log_file, of_log = set_analysis_folder_log(cli.full_path)
  
  print_both('# THIS IS: %s' %(log_file), of_log)
  print_both('Main folder: %s' %(cli.full_path), of_log)
  print_both('LM flag: %s' %(cli.lm_flag), of_log)
  
  err_T0 = set_err_T0(cli.err_T0)
  print_both('Mean T0 uncertainty = %.6f days' %(err_T0), of_log)
  
  grid, grid_names, ngrid = read_grid_summary(cli.full_path, cli.lm_flag, of_log)
  print_both('FOUND %d simulation/s.' %(ngrid), of_log)
  
  mass_convert, mass_unit = mass_factor(cli.mass_type)
  
  t_launch = atime.Time('2018-06-30T06:00:00.0', format='isot', scale='tdb')
  t_end = t_launch.jd + 365.25*3.5
  
  print_both('CHEOPS DATE: LAUNCH = %.4f [BJD] END = %.4f [BJD]' %(t_launch.jd, t_end), of_log)
  print_both(' ', of_log)
  
  hdf5_name, f_h5 = create_hdf5_file(out_folder, cli.seed, grid, grid_names, ngrid, t_launch.jd, t_end, err_T0)
  print_both(' - created hdf5 file: %s' %(hdf5_name), of_log)
  
  # define the number of transits to select:
  nT0_sel = [3, 5, 6, 10, 15]
  # define the number of Monte-Carlo repeatition for each nT0_sel
  nMC = 1000
  
  results_full = []
  
  niter = ngrid
  #niter = 20
  for igrid in range(0, niter):
    results_sim = {}
    sim_num = igrid + 1
    print_both('ANALYSING SIMULATION NUMBER %06d' %(sim_num), of_log)
    names_lst, names_str, units, perturber_parameters = print_perturber_parameters(grid, igrid, 2, mass_convert, mass_unit, of_log)
    # read file
    print_both(' - read trasit file' , of_log)
    T0_all, nT0_all, sel_pre, nT0_pre, sel_cheops, nT0_cheops = read_transit_file(cli.full_path, cli.lm_flag, sim_num, 2, t_min_sel=t_launch.jd, t_max_sel=t_end, of_log=of_log)
    
    if(T0_all is not None):
      T0, eT0 = dirty_transits(T0_all, err_T0)
      print_both(' - dirty transits done', of_log)
      print_both(' | new mean T0 error = %.6f day = %.6f min' %(np.mean(eT0), np.mean(eT0)*1440.), of_log)
    
      epo_pre, Tref_pre, Pref_pre = compute_lin_ephem(T0[sel_pre], eT0[sel_pre])
      print_both(' - computed the linear ephemeris on pre-CHEOPS transits', of_log)
      print_both(' | T0lin = Tref + N x Pref = %.8f + N x %.8f' %(Tref_pre, Pref_pre), of_log)
      
      epo_cheops, oc_cheops, rms_cheops = compute_oc(Tref_pre, Pref_pre, T0[sel_cheops])
      print_both(' - computed O-C = T0cheops - T0lin', of_log)
      print_both(' | rms_TTV (cheops) = %.10f day = %.8f hour = %.6f min %.3f sec' %(rms_cheops, rms_cheops*24., rms_cheops*1440., rms_cheops*86400.), of_log)
      
      print_both(' - computing Monte-Carlo rms of O-C for different number of T0s', of_log)
      rms_mc = [list(compute_rms_mc(oc_cheops, nT0_sel[i_n], nMC)) for i_n in range(len(nT0_sel))]
      
    else: # if the transit file is empty
      print_both(' | transit file is empty ...')
      Tref_pre = 0.
      Pref_pre = 0.
      rms_cheops = 0.
      rms_mc = np.zeros((len(nT0_sel), 2))
    
    
    results_sim['sim_num']              = sim_num
    results_sim['parameter_names']      = names_lst
    results_sim['parameter_units']      = units
    results_sim['perturber_parameters'] = perturber_parameters
    results_sim['linear_ephem']         = [Tref_pre, Pref_pre]
    results_sim['rms_ttv']              = rms_cheops
    results_sim['nT0_sel']              = nT0_sel
    results_sim['nMC']                  = nMC
    results_sim['rms_mc']               = rms_mc  # mean and median OC rms for different nT0_sel during CHEOPS mission
    
    f_h5 = save_results_sim(f_h5, results_sim)
    
    print_both(' | %7s %16s %16s' %('nT0_sel', 'rms_mean[day]', 'rms_median[day]'), of_log)
    for i_sel in range(0,len(nT0_sel)):
      print_both(' | %7d %16.6f %16.6f' %(nT0_sel[i_sel], rms_mc[i_sel][0], rms_mc[i_sel][1]), of_log)
    
    results_full.append(results_sim)
    
    print_both(' ', of_log)
    
  
  # use last sim to determine number of parameters of the perturber
  npar = np.shape(perturber_parameters)[0]
  
  print_both(' - Save final summary', of_log)
  
  # put perturber_parameters / rms_ttv and / rms_mc in a numpy array
  all_perturber_par = np.array([results_full[i_s]['perturber_parameters'] for i_s in range(0,niter)])
  all_rms_ttv = np.array([results_full[i_s]['rms_ttv'] for i_s in range(0,niter)])
  all_rms_mc = np.array([np.reshape(results_full[i_s]['rms_mc'], newshape = (len(nT0_sel)*2)) for i_s in range(0,niter)])
  above_threshold = all_rms_ttv > err_T0
  
  final_summary = np.column_stack((all_perturber_par, all_rms_ttv, all_rms_mc))
  final_summary_above = np.column_stack((all_perturber_par[above_threshold,:], all_rms_ttv[above_threshold], all_rms_mc[above_threshold,:]))
  # very long header... ...
  header = '%s rms_ttv_d %s' %(' '.join('%s_%s'.strip('-').strip('_') %(names_lst[i_x], units[i_x]) for i_x in range(len(units))), ' '.join(['rms_mean_N%d_d rms_median_N%d_d' %(nT0_sel[i_n], nT0_sel[i_n]) for i_n in range(len(nT0_sel))]))
  header_lst = header.split()

  # save into txt files
  np.savetxt(os.path.join(out_folder, 'final_summary.dat'), final_summary, fmt='%23.16e', header=header)
  print_both(' | saved final_summary.dat', of_log)
  
  np.savetxt(os.path.join(out_folder, 'final_summary_above.dat'), final_summary_above, fmt='%23.16e', header=header)
  print_both(' | saved final_summary_above.dat', of_log)

  # save into hdf5 file both final array
  gsum = f_h5.create_group('summary')
  gsum.attrs['npar'] = npar
  gsum.attrs['err_T0'] = err_T0
  gsum.attrs['header'] = header
  print_both(' | summary attrs added', of_log)
  gsum.create_dataset('final_summary', data=final_summary, dtype=np.float64)
  gsum.create_dataset('final_summary_above', data=final_summary_above, dtype=np.float64)
  gsum.create_dataset('above_threshold', data=above_threshold, dtype=bool)
  gsum.create_dataset('header_lst', data=np.array(header_lst, dtype='S20'), dtype='S20')
  
  f_h5.close()
  print_both(' | saved final summary into hdf5 file', of_log)

  # create the plot:
  fig_file = plot_summary(out_folder, hdf5_name, body_id=2, rms_id_in=npar, of_log=of_log)

  of_log.close()
  
  return

# ==============================================================================
# RUN MAIN

if __name__ == "__main__":
  main()

# ==============================================================================

