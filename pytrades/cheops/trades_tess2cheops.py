#!/usr/bin/env python
# -*- coding: utf-8 -*-

# ==============================================================================
# ==============================================================================
# IT READS SYNTHETIC MULTIPLE PLANET SYSTEMS FROM TESS
# IT COMPUTES IF THE SYSTEM IS STABLE AND IF IT SHOWS TTV
# ==============================================================================
# ==============================================================================

# ==============================================================================
# IMPORT

from __future__ import division # no more "zero" integer division bugs!:P
import argparse
import numpy as np # array
import h5py
#import astropy.time as atime
import os
import sys
import gc
import time
#import emcee.interruptible_pool as emcee_pool
#import multiprocessing as mp

import matplotlib as mpl
mpl.use("Agg", warn=False)
import matplotlib.pyplot as plt
plt.rc('font',**{'family':'serif','serif':['Computer Modern Roman']})
plt.rc('text', usetex=True)
#import matplotlib.cm as cm
#from  matplotlib import colors as mc

# custom modules
import gls

script_path = os.path.realpath(__file__)
module_path = os.path.abspath(os.path.join(os.path.dirname(script_path), '../'))
sys.path.append(module_path)

import ancillary as anc
import constants as cst
from pytrades_lib import pytrades

# ==============================================================================
# ==============================================================================

def get_args():
  
  parser = argparse.ArgumentParser(description='TESS-CHEOPS')
  
  parser.add_argument('-f', '--f', '-i', '--i', '-input-file', '--input-file',
                      action='store', dest='input_file', required=True,
                       help='The input csv file with TESS synthetic systems.')
  
  parser.add_argument('-n', '--n', '-np', '--np', '-n-planets', '--n-planets',
                      action='store', dest='n_planets', required=True,
                       help='Number of planets in the system (The planets, not the bodies!).')
  
  parser.add_argument('-b', '--b', '-base', '--base', '-base-folder', '--base-folder',
                      action='store', dest='main_base_folder', required=True,
                       help='Full folder path containing the proper base_#p folder.')
  
  parser.add_argument('-o', '--o', '-output-folder', '--output-folder',
                      action='store', dest='output_folder', required=True,
                       help='Full folder path to store outputs.')
  
  parser.add_argument('-c', '--c', '-cpu', '--cpu',
                      action='store', dest='n_cpu', default='1',
                       help='Number of cpus/threads to use. Default = 1.')
  
  parser.add_argument('-s', '--s', '-seed', '--seed',
                      action='store', dest='seed', default='None',
                       help='Seed for random number generator. Default is None')
  
  parser.add_argument('-r', '--r', '-repeat', '--repeat',
                      action='store', dest='repeat', default='1',
                       help='Number of simulations per system. This allows to compute different mean anomalies for external planets. Default is 1')
  
  cli = parser.parse_args()
  cli.input_file = os.path.abspath(cli.input_file)
  cli.main_base_folder = os.path.abspath(cli.main_base_folder)
  cli.output_folder = os.path.abspath(cli.output_folder)
  
  try:
    cli.n_cpu = int(cli.n_cpu)
    if(cli.n_cpu < 1 ): cli.n_cpu = 1
  except:
    cli.n_cpu = 1
  
  try:
    cli.seed = int(cli.seed)
  except:
    cli.seed = None
  
  try:
    cli.repeat = int(cli.repeat)
  except:
    cli.repeat = 1
  
  return cli

# ==============================================================================

class set_inout:
  
  def __init__(self, cli):
    self.input_file = cli.input_file
    self.main_base_folder = cli.main_base_folder
    self.n_planets_str = cli.n_planets
    self.n_planets = int(cli.n_planets)
    self.output_folder = cli.output_folder
    self.seed = cli.seed
    self.base_folder = os.path.join(os.path.join(
      cli.main_base_folder, 'base_%sp' %(self.n_planets_str)), '')
    self.sub_output_folder = os.path.join(os.path.join(
      cli.output_folder, '%sp_systems_seed%d' %(cli.n_planets, cli.seed)), '')

  def print_io(self):
    print
    print 'INPUT FILE'
    print self.input_file
    print 'NUMBER OF PLANETS = ', self.n_planets_str
    print 'MAIN BASE FOLDER'
    print self.main_base_folder
    print 'BASE FOLDER'
    print self.base_folder
    print 'OUTPUT FOLDER'
    print self.output_folder
    print 'SUB OUTPUT FOLDER'
    print self.sub_output_folder
    print
    
# ==============================================================================

class input_data:
  
  def __init__(self, csv_read, n_planets):
    self.csv_read = csv_read
    self.col_names = list(self.csv_read.dtype.names)
    self.n_cols = len(self.col_names)
    self.n_rows = np.shape(self.csv_read[self.col_names[0]])[0]
    self.n_systems = int(self.n_rows / n_planets)
    self.data = np.array([self.csv_read[self.col_names[icol]] for icol in range(self.n_cols)], dtype=np.float64).T

# ==============================================================================

def read_inputfile(input_file, n_planets):
  
  # checks if file exists
  if(os.path.isfile(input_file)):
    csv_read = np.genfromtxt(input_file, delimiter=',', names=True)
    data_obj = input_data(csv_read, n_planets)
  else:
    sys.exit('THE INPUT DOES NOT EXIST')
  
  return data_obj

# ==============================================================================

class transit_times:
  
  def __init__(self, planet_id, TTs_full, id_TT_full, stats_TT):
    self.planet_id = planet_id
    sel_TTs = np.logical_and(id_TT_full == planet_id, stats_TT)
    self.TTs = np.sort(TTs_full[sel_TTs])
    #self.nTTs = np.shape(self.TTs)[0]
    self.nTTs = np.sum(sel_TTs)
    self.epo = 0
    self.TTref = 0.
    self.Pref = 1.
    self.linTTs = 0.
    self.OCs = np.zeros((self.nTTs))
    self.P_TTV_d = 0.
    self.amp_TTV_d = 0.
    self.amp_TTV_m = 0.
    
  def set_linearephem(self):
    #self.epo, self.TTref, self.Pref = anc.compute_lin_ephem(self.TTs, np.ones(self.nTTs)/86400.)
    self.epo, self.TTref, self.Pref = anc.compute_lin_ephem(self.TTs)
    #print self.planet_id, self.TTref, self.Pref
    
  def set_oc(self, int_time):
    self.linTTs = self.TTref + self.epo * self.Pref
    self.OCs = self.TTs - self.linTTs
    Pstart = 1.5*self.Pref
    Pfinish = 2. * int_time
    #print Pstart, Pfinish, '(', self.Pref, ')'
    gls_o = gls.Gls((self.TTs, self.OCs, np.ones(self.nTTs)/86400.),
                    Pbeg=Pstart,
                    Pend=Pfinish, 
                    ofac=10,
                    verbose = False)
    self.P_TTV_d = 1./gls_o.hpstat['fbest']
    self.amp_TTV_d = gls_o.hpstat['amp']
    self.amp_TTV_m = self.amp_TTV_d * 1440.
    del gls_o
    
# -------------------------------------

def set_complete_transits(int_time, planet_id, TTs_full, id_TT_full, stats_TT):
  
  TT_obj = transit_times(planet_id, TTs_full, id_TT_full, stats_TT)
  print '%s nTTs = %d' %(anc.letters[planet_id-1], TT_obj.nTTs)
  if(TT_obj.nTTs > 2):
    TT_obj.set_linearephem()
    if(TT_obj.Pref != 0.):
      TT_obj.set_oc(int_time)
  
  return TT_obj
    
# ==============================================================================

class exoSystem:
  
  def __init__(self, id_system, id_sim, n_bodies):
    self.id_system = id_system
    self.id_sim = id_sim
    self.n_bodies = n_bodies
    self.n_planets = n_bodies - 1
    self.mass = np.ones((n_bodies))
    self.radius = np.ones((n_bodies))
    self.period = np.zeros((n_bodies))
    self.period[1:] = 365.25
    self.sma = np.zeros((n_bodies))
    self.sma[1:] = 1.
    self.ecc = np.zeros((n_bodies))
    self.argp = np.zeros((n_bodies))
    self.argp[1:] = 90.
    self.meanA = np.zeros((n_bodies))
    self.inc = np.zeros((n_bodies))
    self.inc[1:] = 90.
    self.longN = np.zeros((n_bodies))
    self.longN[1:] = 180.
    self.stable = True
  
  def update_system(self, data_in):
    self.mass[0] = data_in[0,0]
    self.mass[1:] = data_in[0:, 2]*cst.Mears
    self.radius[0] = data_in[0,1]
    self.radius[1:] = data_in[0:, 3]*cst.Rears
    self.period[1:] = data_in[0:, 4]
    self.sma[1:] = data_in[0:, 5]
    self.inc[1:] = np.arccos(data_in[0:, 6])*180./np.pi

  def set_transits(self, transits):
    self.transits = transits
    
# ==============================================================================

def run_simulation(int_time, id_system, id_sim, n_bodies, parameters_in):
  
  print 'SYSTEM NUMBER %d: SIM %d' %(id_system, id_sim)
  n_planets = n_bodies - 1
  system = exoSystem(id_system, id_sim, n_bodies)
  system.update_system(parameters_in)
  system.meanA[2:] = np.random.random(size = n_planets-1)*360.
  nTT_max, nRV_max = pytrades.get_max_nt0_nrv(system.period, n_bodies)
  
  TTs_full, id_TT_full, stats_TTs, time_rv_nmax, rv_nmax, stats_rv = pytrades.wrapper_run_grid_combination(system.mass,
                                        system.radius,
                                        system.period, 
                                        system.sma, 
                                        system.ecc, 
                                        system.argp, 
                                        system.meanA, 
                                        system.inc, 
                                        system.longN, 
                                        nTT_max, nRV_max)
  
  transits = [set_complete_transits(int_time, ibd, TTs_full, id_TT_full, stats_TTs)
              for ibd in range(2, n_bodies+1)]
  if(np.any(stats_TTs)):
    idx_max_TT = np.argmax(TTs_full)
    max_TT = TTs_full[idx_max_TT]
    if(max_TT+system.period[id_TT_full[idx_max_TT]-1] < int_time): system.stable = False
  else:
    system.stable = False
  
  system.set_transits(transits)
  
  return system

# ==============================================================================

def save_plot_oc(inout, int_time, system):
  
  # creates output folders
  if(not os.path.isdir(inout.sub_output_folder)):
    os.makedirs(inout.sub_output_folder)
  
  # save system data into hdf5 file
  output_file = os.path.join(inout.sub_output_folder,
                             '%04d_%04d_system.hdf5' %(system.id_system, system.id_sim))
  h5f = h5py.File(output_file, 'w')
  gg = h5f.create_group('%04d_%04d' %(system.id_system, system.id_sim))
  gg.attrs['n_bodies'] = system.n_bodies
  gg.attrs['n_planets'] = system.n_planets
  gg.attrs['stable'] = system.stable
  gg.attrs['seed'] = inout.seed
  
  gg.create_dataset('mass', data=system.mass, dtype=np.float64)
  gg['mass'].attrs['unit'] = 'Msun'
  gg.create_dataset('radius', data=system.radius, dtype=np.float64)
  gg['radius'].attrs['unit'] = 'Rsun'
  gg.create_dataset('period', data=system.period, dtype=np.float64)
  gg['period'].attrs['unit'] = 'days'
  gg.create_dataset('sma', data=system.sma, dtype=np.float64)
  gg['sma'].attrs['unit'] = 'au'
  gg.create_dataset('ecc', data=system.ecc, dtype=np.float64)
  gg.create_dataset('argp', data=system.argp, dtype=np.float64)
  gg['argp'].attrs['unit'] = 'deg'
  gg.create_dataset('meanA', data=system.meanA, dtype=np.float64)
  gg['meanA'].attrs['unit'] = 'deg'
  gg.create_dataset('inc', data=system.inc, dtype=np.float64)
  gg['inc'].attrs['unit'] = 'deg'
  gg.create_dataset('longN', data=system.longN, dtype=np.float64)
  gg['longN'].attrs['unit'] = 'deg'
  
  # save planet TTs (and stuff) and create plot
  output_fig = os.path.splitext(output_file)[0]
  if(system.stable):
    sys_stat = 'stable'
  else:
    sys_stat = 'unstable'
    
    
  deltax = 0.03*int_time
  xminl = 0 - deltax
  xmaxl = int_time + deltax
  
  fig = plt.figure(figsize=(12,12))
  fig.subplots_adjust(hspace=0.2)
  fig.suptitle('system %04d\_%04d: %s' %(system.id_system, system.id_sim, sys_stat))
  
  id_name = anc.letters
  for ipl in range(0,system.n_planets):
    planet_str = id_name[ipl+1]
    oo = system.transits[ipl]
    og = gg.create_group('planet_%s' %(planet_str))
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
    
    ax = plt.subplot2grid((system.n_planets, 1), (ipl,0))
    tlin_str = r'planet %s ($P=%.6f$): $TT_\mathrm{lin} = %.6f + N \times %.6f$' %(planet_str, system.period[ipl+1], oo.TTref, oo.Pref)
    #tlin_str = 'planet %s: TT(lin) = %.6f + N x %.6f' %(planet_str, oo.TTref, oo.Pref)
    ttv_str = r'$P_\mathrm{TTV}=%.6f\ \mathrm{day}\ \mathrm{Amp}_\mathrm{TTV}=%.6f\ \mathrm{min}$' %(oo.P_TTV_d, oo.amp_TTV_m)
    #ttv_str = 'P(TTV)=%.6f day Amp(TTV)=%.6f min' %(oo.P_TTV_d, oo.amp_TTV_m)
    ax.set_title(r'%s , %s' %(tlin_str, ttv_str), fontsize=9)
    
    ax.axhline(0., color='black')
    if(oo.nTTs > 0):
      xTT = oo.TTs - oo.TTs[0]
      ax.plot(xTT, oo.OCs*1440., 
              color='lightgray', marker='o', ms=4, mfc='#1f77b4', mec='None',ls='-')
    
    ax.set_xlim([xminl, xmaxl])
    ax.set_xlabel('Time (days)')
    ax.set_ylabel('O-C (min)')
    
  h5f.close()
    
  #fig.savefig('%s.png' %(output_fig), bbox_inches='tight', dpi=300)
  fig.tight_layout(rect=[0, 0.03, 1, 0.95])
  fig.savefig('%s.png' %(output_fig), dpi=300)
  plt.close(fig)
  
  print 'WRITTEN HDF5 FILE: %s' %(output_file)
  print 'WRITTEN IMAGE: %s.png' %(output_fig)

  return
  

# ==============================================================================

def run_sim_and_plot(int_time, id_system, id_sim, inout, n_bodies, data_in):
  
  system = run_simulation(int_time, id_system, id_sim, n_bodies, data_in)
  save_plot_oc(inout, int_time, system)
  
  return system

# ==============================================================================

def main():
  
  # command line input arguments
  cli = get_args()
  
  np.random.seed(cli.seed)
  
  # set paths and folders
  inout = set_inout(cli)
  inout.print_io()
  
  # init trades with proper base folder
  print 'SET TRADES'
  pytrades.initialize_trades(inout.base_folder, '', cli.n_cpu)
  n_bodies = pytrades.n_bodies
  int_time = pytrades.tint
  print 'NUMBER OF BODIES = ', n_bodies
  print 'INTEGRATION TIME = ', int_time, ' days'
  print
  
  # reads input data
  data_obj = read_inputfile(inout.input_file, inout.n_planets)
  
  ## testing
  #test = data_obj.data[0:inout.n_planets, :]
  #run_sim_and_plot(int_time, 1, inout, n_bodies, test)
  
  systems = [run_sim_and_plot(int_time, (ii/inout.n_planets)+1, jj+1, inout, n_bodies, data_obj.data[ii:ii+inout.n_planets, :]) for ii in range(0, data_obj.n_rows, inout.n_planets) for jj in range(0,cli.repeat)]
 
 
  return

# ==============================================================================
# ==============================================================================

if __name__ == "__main__":
  main()

# ==============================================================================
