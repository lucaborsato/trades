#!/usr/bin/env python
# -*- coding: utf-8 -*-

from __future__ import division # no more "zero" integer division bugs!:P
import argparse
import os
import numpy as np # array
import pytrades_lib
import emcee
import h5py
import sys
import time
import shutil
from constants import Mjups, Mears, Msear
import ancillary as anc

#from matplotlib import use as mpluse
##mpluse("Agg")
#mpluse("Qt4Agg")
#import matplotlib.pyplot as plt
#plt.rc('font',**{'family':'serif','serif':['Computer Modern Roman']})
#plt.rc('text', usetex=True)

# ---
# path argument
def get_args():
  parser = argparse.ArgumentParser(description='TRADES+PSO+EMCEE')
  
  # PATH FOLDER: full_path
  parser.add_argument('-p', '--path', action='store', dest='full_path', required=True, help='The path (absolute or relative) with simulation files for TRADES.')

  # SUBFOLDER TO SAVE ALL SIMULATION DATA
  parser.add_argument('-s', '--sub-folder', '--sb', action='store', dest='sub_folder', default='emcee_run', help='Sub-folder name, without full path. Default = emcee_run')

  cli = parser.parse_args()
  
  cli.full_path = os.path.join(os.path.abspath(cli.full_path), '')
  return cli


#
# INITIALISE FOLDER AND LOG FILE
#
def init_folder(working_path, sub_folder):
  working_folder = os.path.join(working_path, sub_folder)
  if (not os.path.isdir(working_folder)):
      os.makedirs(working_folder)
  run_log = os.path.join(working_folder, "trades_run.log")
  of_run = open(run_log, 'w')
  of_run.write("# pyTRADES LOG FILE\n")
  of_run.write("# working_path = %s\n" %(working_path))
  of_run.write("# working_folder = %s\n" %(working_folder))
  of_run.write("# run_log = %s\n" %(run_log))
  
  return working_folder, run_log, of_run

# set names of all parameters
def set_all_parameters_names(n_planets):
  id_star = ['M_a_sun', 'R_a_sun']
  #mass_id = [0]
  mass_id = []
  letter_planets = ['b', 'c', 'd', 'e', 'f', 'g', 'h', 'i']
  elements = ['M', 'R', 'P', 'e', 'w', 'mA', 'i', 'lN']
  units_elements = ['earth', 'earth', 'days', '', 'deg', 'deg', 'deg', 'deg']
  all_parameters_names = np.zeros((2+n_planets*len(elements))).astype(str)
  all_parameters_names[:2] = np.asarray(id_star)
  #print all_parameters_names
  for i in range(n_planets):
    temp = ['%s_%s_%s' %(elements[j], letter_planets[i], units_elements[j]) for j in range(len(elements))]
    #print temp
    aaa = 2+8*i
    bbb = 2+8*(i+1)
    #print aaa, bbb
    all_parameters_names[aaa:bbb] = np.asarray(temp)
    mass_id.append(aaa)
  #print all_parameters_names
  return all_parameters_names, mass_id

# ---

#
# --- main --- #
#

# STARTING TIME
start = time.time()

# ----
# -- 1
# get command arguments

cli = get_args()

# ----
# -- 2
# read MCMC-RV parameters
RV_sims_file = '/home/borsato/Research/Kepler/Kepler-19/RV/2016-02-08_output_byRanger.dat'
RV_names = np.asarray(['P2', 'phi2', 'm2', 'P3', 'phi3', 'm3', 'P4', 'phi4', 'm4'], dtype='S15')
RV_names_units = np.asarray(['P_b_days', 'phi_b_rad', 'm_b_jup', 'P_c_days', 'phi_c_rad', 'm_c_jup', 'P_d_days', 'phi_d_rad', 'm_d_jup'], dtype='S15')
m_id = [2, 5, 8]
RV_sims = np.genfromtxt(RV_sims_file)
#RV_sims[:,m_id] = RV_sims[:,m_id]*Mears
RV_sims[:,m_id] = RV_sims[:,m_id]*Mjups
n_RV_sims = RV_sims.shape[0]

# ----
# -- 3
# initialise TRADES:
#   * init TRADES arguments -> read integration arguments, T0 data, no RV?

working_path = cli.full_path

# INITIALISE TRADES WITH SUBROUTINE WITHIN TRADES_LIB -> PARAMETER NAMES, MINMAX, INTEGRATION ARGS, READ DATA ...
pytrades_lib.pytrades.initialize_trades(working_path, cli.sub_folder)
# RETRIEVE DATA AND VARIABLES FROM TRADES_LIB MODULE
parameters_minmax = pytrades_lib.pytrades.parameters_minmax # PARAMETER BOUNDARIES
n_bodies = pytrades_lib.pytrades.n_bodies # NUMBER OF TOTAL BODIES OF THE SYSTEM
n_planets = n_bodies - 1 # NUMBER OF PLANETS IN THE SYSTEM
ndata = pytrades_lib.pytrades.ndata # TOTAL NUMBER OF DATA AVAILABLE
npar  = pytrades_lib.pytrades.npar # NUMBER OF TOTAL PARAMATERS ~n_planets X 6
nfit  = pytrades_lib.pytrades.nfit # NUMBER OF PARAMETERS TO FIT
dof   = pytrades_lib.pytrades.dof # NUMBER OF DEGREES OF FREEDOM = NDATA - NFIT

# RADIAL VELOCITIES SET
n_rv = pytrades_lib.pytrades.nrv
n_set_rv = pytrades_lib.pytrades.nrvset

# TRANSITS SET
n_t0 = pytrades_lib.pytrades.nt0
n_t0_sum = np.sum(n_t0)
n_set_t0 = 0
for i in range(0, n_bodies):
  if (np.sum(n_t0[i]) > 0): n_set_t0 += 1

# READ THE NAMES OF THE PARAMETERS FROM THE TRADES_LIB AND CONVERT IT TO PYTHON STRINGS
reshaped_names = pytrades_lib.pytrades.parameter_names.reshape((10,nfit), order='F').T
parameter_names = [''.join(reshaped_names[i,:]).strip() for i in range(0,nfit)]
  
# INITIALISE SCRIPT FOLDER/LOG FILE
working_folder, run_log, of_run = init_folder(working_path, cli.sub_folder)

print 
print ' =================== '
print ' pyTRADES - MCMC2PSO '
print ' =================== '
print
print ' WORKING PATH = %s' %(working_path)
#print ' NUMBER OF THREADS = %d' %(nthreads)
print ' dof = ndata(%d) - nfit(%d) = %d' %(ndata, nfit, dof)
print ' npar = %d' %(npar)
print ' Total N_RV = %d for %d set(s)' %(n_rv, n_set_rv)
print ' Total N_T0 = %d for %d out of %d planet(s)' %(n_t0_sum, n_set_t0, n_planets)
print
#print ' %s = %.7f' %('log constant error = ', ln_err_const)
#of_run.write(' NUMBER OF THREADS = %d\n' %(nthreads))
of_run.write(' dof = ndata(%d) - nfit(%d) = %d\n' %(ndata, nfit, dof))
of_run.write(' npar = %d\n' %(npar))
of_run.write(' Total N_RV = %d for %d set(s)\n' %(n_rv, n_set_rv))
of_run.write(' Total N_T0 = %d for %d out of %d planet(s)\n' %(n_t0_sum, n_set_t0, n_planets))
#of_run.write(' %s = %.7f' %('log constant error = ', ln_err_const))
of_run.write('\n')

all_parameters_names, mass_id = set_all_parameters_names(n_planets)
#sys.exit()

# ----
# -- 4
# PSO & FITNESS
# save final parameter set for each simulation and fitness

all_parameters_summary = np.zeros((n_RV_sims, npar+1)) + 1.e10
parameters_summary = np.zeros((n_RV_sims, nfit+1)) + 1.e10
# INITIALISE PSO ARGUMENTS FROM pso.opt FILE
pytrades_lib.pytrades.init_pso(1,working_path) # read PSO options
# PSO VARIABLES
np_pso = pytrades_lib.pytrades.np_pso
nit_pso = pytrades_lib.pytrades.nit_pso
n_global = pytrades_lib.pytrades.n_global

# prepare to save all the parameters into file
parameters_hdf5 = h5py.File(os.path.join(working_folder, 'parameters_summary.hdf5'), 'w')
# save original RV simulation parameters
parameters_hdf5.create_dataset('RV_names', data=RV_names_units, dtype='S15')
parameters_hdf5.create_dataset('RV_sims', data=RV_sims, dtype=np.float64)
# prepare datasets for the PSO fitting parameters
parameters_hdf5.create_dataset('all_parameters_names', data=all_parameters_names, dtype='S15')
parameters_hdf5.create_dataset('all_parameters_summary', (n_RV_sims, npar+1), dtype=np.float64)
parameters_hdf5.create_dataset('parameter_names', data=parameter_names, dtype='S15')
parameters_hdf5.create_dataset('parameters_summary', (n_RV_sims, nfit+1), dtype=np.float64)
parameters_hdf5.close()

# run PSO for each combination of the MCMC-RV parameters
for i in range(0,n_RV_sims):
#for i in range(0,2):
  pso_start = time.time()
  i_sim = i + 1
  print ' SIM # %07d: npop = %d niter = %d' %(i_sim, np_pso, nit_pso)
  #print RV_sims[i,:]
  pytrades_lib.pytrades.init_fix_parameters(RV_names, RV_sims[i,:])
  pso_path = os.path.join(os.path.join(working_folder, '%07d_sim' %(i_sim)), '')
  pytrades_lib.pytrades.path_change(pso_path)
  if(not os.path.exists(pso_path)): os.makedirs(pso_path)
  pso_parameters, pso_fitness = pytrades_lib.pytrades.run_pso(i_sim, nfit)
  pso_best_evolution = np.asarray(pytrades_lib.pytrades.pso_best_evolution[:,:], dtype=np.float64)
  #best_fitness, best_lgllhd, best_check  = pytrades_lib.pytrades.write_summary_files(i_sim,pso_parameters)
  
  # save all the parameters into array
  temp_all_parameters_in = np.asarray(pytrades_lib.pytrades.system_parameters, dtype=np.float64)
  temp_all_parameters_out = pytrades_lib.pytrades.update_system_parameters(temp_all_parameters_in, pso_parameters)
  temp_all_parameters_earth = temp_all_parameters_out
  temp_all_parameters_earth[mass_id] = temp_all_parameters_out[mass_id]*Msear
  all_parameters_summary[i,0:npar] = temp_all_parameters_earth
  all_parameters_summary[i,-1] = pso_fitness
  parameters_summary[i,0:nfit] = pso_parameters
  parameters_summary[i,-1] = pso_fitness
  
  # store best parameters and fitness into file
  parameters_hdf5 = h5py.File(os.path.join(working_folder, 'parameters_summary.hdf5'), 'a')
  temporary_all_hdf5 = parameters_hdf5['all_parameters_summary']
  temporary_all_hdf5[i,0:npar]= temp_all_parameters_out
  temporary_all_hdf5[i,-1] = pso_fitness
  temporary_par_hdf5 = parameters_hdf5['parameters_summary']
  temporary_par_hdf5[i,0:nfit]= pso_parameters
  temporary_par_hdf5[i,-1] = pso_fitness
  parameters_hdf5.close()
  
  # SAVE PSO SIMULATION IN pso_run.hdf5 FILE
  pso_hdf5 = h5py.File(os.path.join(pso_path, 'pso_run.hdf5'), 'w')
  pso_hdf5.create_dataset('population', data=pytrades_lib.pytrades.population, dtype=np.float64)
  pso_hdf5.create_dataset('population_fitness', data=pytrades_lib.pytrades.population_fitness, dtype=np.float64)
  pso_hdf5.create_dataset('pso_parameters', data=pso_parameters, dtype=np.float64)
  pso_hdf5.create_dataset('pso_fitness', data=np.array(pso_fitness), dtype=np.float64)
  pso_hdf5.create_dataset('pso_best_evolution', data=pso_best_evolution, dtype=np.float64)
  pso_hdf5.create_dataset('parameters_minmax', data=parameters_minmax, dtype=np.float64)
  pso_hdf5.create_dataset('parameter_names', data=parameter_names, dtype='S10')
  pso_hdf5.close()

  pso_end = time.time()
  pso_time_d, pso_time_h, pso_time_m, pso_time_s = anc.computation_time(pso_end-pso_start)
  print ' PSO FINISHED in %2d day %02d hour %02d min %.2f sec - bye bye' %(int(pso_time_d), int(pso_time_h), int(pso_time_m), pso_time_s)
  print

## save all the parameters into file
#parameters_hdf5 = h5py.File(os.path.join(working_folder, 'parameters_summary.hdf5'), 'w')
#parameters_hdf5.create_dataset('parameters_summary', data=parameters_summary, dtype=np.float64)
#parameters_hdf5.close()

pytrades_lib.pytrades.path_change(os.path.join(working_folder,''))
sort_idx = np.argsort(parameters_summary[:,-1])
best_id = sort_idx[0] + 1

best_all_parameters = all_parameters_summary[sort_idx[0],0:npar]
best_parameters = parameters_summary[sort_idx[0],0:nfit]
best_pso_fitness = parameters_summary[sort_idx[0],-1]

all_parameters_good = all_parameters_summary[all_parameters_summary[:,-1]<1.e10, 0:npar]
#print all_parameters_good.shape
residuals = all_parameters_good - best_all_parameters
sigmas_parameters = np.percentile(residuals, [0.13, 2.28, 15.87, 50., 84.13, 97.72, 99.87], axis=0)
#print sigmas_parameters.shape

pytrades_lib.pytrades.init_fix_parameters(RV_names, RV_sims[sort_idx[0],:])
best_fitness, best_lgllhd, best_check  = pytrades_lib.pytrades.write_summary_files(best_id, best_parameters)

print
line = ' sort id = %d => best id  = %d' %(sort_idx[0], best_id)
print line
of_run.write('%s\n' %(line))
line = ' best fitness = %.7f (%.7f)' %(best_pso_fitness, best_fitness)
print line
of_run.write('%s\n' %(line))
line = ' %15s = %19s - %19s + %19s' %('id', 'best_value', '-1sigma', '+1sigma')
print line
of_run.write('%s\n' %(line))
for i in range(0, npar):
  line = '%15s = %19.12f - %19.12f + %19.12f' %(all_parameters_names[i], best_all_parameters[i], sigmas_parameters[2,i], sigmas_parameters[4,i])
  print line
  of_run.write('%s\n' %(line))
print

# ----
# -- 5
# closing stuff: deallocate and print time elapsed
elapsed = time.time() - start
elapsed_d, elapsed_h, elapsed_m, elapsed_s = anc.computation_time(elapsed)
print
print ' mcmc2fitness FINISHED in %2d day %02d hour %02d min %.2f sec - bye bye' %(int(elapsed_d), int(elapsed_h), int(elapsed_m), elapsed_s)
of_run.write(' mcmc2fitness FINISHED in %2d day %02d hour %02d min %.2f sec - bye bye\n' %(int(elapsed_d), int(elapsed_h), int(elapsed_m), elapsed_s))
of_run.close()
print 
pytrades_lib.pytrades.deallocate_variables()

# ----



