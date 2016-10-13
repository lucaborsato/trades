#!/usr/bin/env python
# -*- coding: utf-8 -*-

from __future__ import division # no more "zero" integer division bugs!:P
import argparse
import os
import numpy as np # array
import pytrades_lib
#import emcee
import h5py
import sys
import time
import shutil
#from emcee.utils import MPIPool
from constants import Mjups
from ancillary import *
#from pyde.de import DiffEvol as devol
import pymultinest
from pymultinest.solve import solve

#from matplotlib import use as mpluse
##mpluse("Agg")
#mpluse("Qt4Agg")
#import matplotlib.pyplot as plt
#plt.rc('font',**{'family':'serif','serif':['Computer Modern Roman']})
#plt.rc('text', usetex=True)

# 
# 
#
def get_args():
  parser = argparse.ArgumentParser(description='TRADES+pyMultiNest')
  
  # PATH FOLDER: full_path
  parser.add_argument('-p', '--path', action='store', dest='full_path', required=True, help='The path (absolute or relative) with simulation files for TRADES.')

  # SUBFOLDER TO SAVE ALL SIMULATION DATA
  parser.add_argument('-s', '--sub-folder', '--sb', action='store', dest='sub_folder', default='emcee_run', help='Sub-folder name, without full path. Default = emcee_run')
  
  ## HOW TO RUN PSO OR SKIP
  #parser.add_argument('-r', '--pso', '--pso-type', action='store', dest='pso_type', default='run', help='Define PSO run type: "skip" = do not run PSO, only emcee with random walkers; "run" = run PSO normally; "exists" = do not run PSO, but read previous PSO (pso_run.hdf5) file and start emcee from its population')
  
  ## NUMBER OF CPU TO USE WITH EMCE !! WARNING !! TO USE PARALLEL IN PSO BEFORE THIS SCRIPT TYPE:
  ## export OMP_NUM_THREADS=CPUNUMBER
  #parser.add_argument('-c', '--cpu', '--nthreads', action='store', dest='nthreads', default=1, help='Number of threads to use. default nthreads = 1.')
  
  ## NUMBER OF WALKERS TO USE WITH EMCEE
  #parser.add_argument('-nw', '--nwalkers', '-np', '--npop', action='store', dest='nwalkers', default=1, help='Number of walkers (or number of chains) to use. default nwalkers = nfit*2')
  
  ## NUMBER OF STEPS/RUNS TO DO FOR EACH WALKER OF EMCEE
  #parser.add_argument('-nr', '--nruns', '-ns', '--nsteps', action='store', dest='nruns', default=10000, help='Number of runs/steps to use for each chain. default nruns = 10000.')
  
  ## NUMBER OF STEPS/RUNS TO SAVE TEMPORARY EMCEE SIMULATION
  #parser.add_argument('--isave', '--iter-save', '--iterations-save', action='store', dest='nsave', default='False', help='Number of iterations to do for save temporary chain. default each 0.1 of nruns.')
  
  ## NUMBER OF BURNIN TO SKIP FOR PRELIMINARY ANALYSIS
  #parser.add_argument('-nb', '--nburn', '--npost', action='store', dest='npost', default=1000, help='Number of burn in, or number of posterior to discard at the beginning of the each chain. default npost = 1000.')
  
  # COMPUTE OR NOT CONSTANT TO ADD TO THE LGLIKELIHOOD: - (1/2dof) * SUM( ln 2pi * sigma_obs^2 )
  #parser.add_argument('-l', '--ln-err', '--ln-err-const', action='store', dest='ln_flag', default=True, help='Computes or not constant to add to the lglikelihood: - (1/2dof) * SUM( ln 2pi * sigma_obs^2 ). Default = True. Set to False if not use this value.')
  parser.add_argument('-l', '--ln-err', '--ln-err-const', action='store', dest='ln_flag', default=True, help='Computes or not constant to add to the lglikelihood: SUM( ln 2pi * sigma_obs^2 ). Default = True. Set to False if not use this value.')

  
  cli = parser.parse_args()
  
  cli.full_path = os.path.join(os.path.abspath(cli.full_path), '')
  cli.sub_folder = os.path.join(os.path.relpath(cli.sub_folder), '')
  
  #if (cli.pso_type not in ['skip', 'run', 'exists']):
    #cli.pso_type = 'run'
  
  #cli.nthreads = int(cli.nthreads)
  
  #cli.nwalkers = int(cli.nwalkers)
  #cli.nruns = int(cli.nruns)
  #cli.npost = int(cli.npost)
  
  cli.ln_flag = set_bool_argument(cli.ln_flag)
    
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

def compute_ln_err_const(ndata, dof, e_RVo, e_T0o, ln_flag):
  if (ln_flag):
    eRV = e_RVo[e_RVo > 0.]
    eT0 = e_T0o[e_T0o > 0.]
    
    ln_e_RVo = np.abs(np.sum(np.log(eRV*eRV)))
    ln_e_T0o = np.abs(np.sum(np.log(eT0*eT0)))

    #ln_err_const = -(0.5/dof) * (ndata*np.log(2.*np.pi) + ln_e_RVo + ln_e_T0o)
    ln_err_const = np.log(2.*np.pi) + ln_e_RVo + ln_e_T0o
  else:
    ln_err_const = 0.
  return ln_err_const


#
# PSO SIMULATION TO EMCEE p0
#


# MAIN -- TRADES + pyMultiNest

# READ COMMAND LINE ARGUMENTS
cli = get_args()

# STARTING TIME
start = time.time()

# RENAME 
working_path = cli.full_path
#nthreads=cli.nthreads

# INITIALISE TRADES WITH SUBROUTINE WITHIN TRADES_LIB -> PARAMETER NAMES, MINMAX, INTEGRATION ARGS, READ DATA ...
pytrades_lib.pytrades.initialize_trades(working_path, cli.sub_folder)

# RETRIEVE DATA AND VARIABLES FROM TRADES_LIB MODULE
fitting_parameters = pytrades_lib.pytrades.fitting_parameters # INITIAL PARAMETER SET (NEEDED ONLY TO HAVE THE PROPER ARRAY/VECTOR)
parameters_minmax = pytrades_lib.pytrades.fitting_parametersameters_minmax # PARAMETER BOUNDARIES
delta_parameters = np.abs(parameters_minmax[:,1] - parameters_minmax[:,0]) # DELTA BETWEEN MAX AND MIN OF BOUNDARIES

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


# compute global constant for the loglhd
global ln_err_const

try:
  e_RVo = np.asarray(pytrades_lib.pytrades.ervobs[:], dtype=np.float64) # fortran variable RV in python will be rv!!!
except:
  e_RVo = np.asarray([0.], dtype=np.float64)
try:
  e_T0o = np.asarray(pytrades_lib.pytrades.et0obs[:,:], dtype=np.float64).reshape((-1))
except:
  e_T0o = np.asarray([0.], dtype=np.float64)
ln_err_const = compute_ln_err_const(ndata, dof, e_RVo, e_T0o, cli.ln_flag)

# READ THE NAMES OF THE PARAMETERS FROM THE TRADES_LIB AND CONVERT IT TO PYTHON STRINGS
reshaped_names = pytrades_lib.pytrades.fitting_parametersameter_names.reshape((10,nfit), order='F').T
parameter_names = [''.join(reshaped_names[i,:]).strip() for i in range(0,nfit)]

# INITIALISE SCRIPT FOLDER/LOG FILE
working_folder, run_log, of_run = init_folder(working_path, cli.sub_folder)

print 
print ' ==================== '
print ' pyTRADES-pyMultiNest'
print ' ==================== '
print
print ' WORKING PATH = %s' %(working_path)
#print ' NUMBER OF THREADS = %d' %(nthreads)
print ' dof = ndata(%d) - nfit(%d) = %d' %(ndata, nfit, dof)
print ' Total N_RV = %d for %d set(s)' %(n_rv, n_set_rv)
print ' Total N_T0 = %d for %d out of %d planet(s)' %(n_t0_sum, n_set_t0, n_planets)
print ' %s = %.7f' %('log constant error = ', ln_err_const)
#of_run.write(' NUMBER OF THREADS = %d\n' %(nthreads))
of_run.write(' dof = ndata(%d) - nfit(%d) = %d\n' %(ndata, nfit, dof))
of_run.write(' Total N_RV = %d for %d set(s)\n' %(n_rv, n_set_rv))
of_run.write(' Total N_T0 = %d for %d out of %d planet(s)\n' %(n_t0_sum, n_set_t0, n_planets))
of_run.write(' %s = %.7f\n' %('log constant error = ', ln_err_const))

# save parameter_names and boundaries to be read by a script
trades_hdf5 = h5py.File(os.path.join(working_folder, 'system_summary.hdf5'), 'w')
trades_hdf5.create_dataset('parameter_names', data=parameter_names, dtype='S10')
trades_hdf5.create_dataset('parameters_minmax', data=parameters_minmax, dtype=np.float64)
trades_hdf5.create_dataset('ln_err_const', data=np.asarray([ln_err_const], dtype=np.float64), dtype=np.float64)
trades_hdf5.close()

# MULTINEST HERE
#output_mnest = os.path.join(working_folder, 'trades_mnest_')
#output_mnest = os.path.join(cli.sub_folder, 'trades_mnest_')
output_mnest = 'trades_mnest_'
os.chdir(working_folder)
log_zero_value = -0.5*1.e8
seed_value = 292929
n_pop = nfit+int(nfit*0.5)
n_update = 2
resume_flag = False
multi_node_flag = True
mpi_onoff = False
print ' pyMultiNest parameters:'
print ' folder = %s\n seed_value = %d , n_pop = %d , n_update = %d\n resume_flag = %r , multi_node_flag = %r\n' %(output_mnest, seed_value, n_pop, n_update, resume_flag, multi_node_flag)
of_run.write(' pyMultiNest parameters:\n')
of_run.write(' folder = %s\n seed_value = %d , n_pop = %d , n_update = %d\n resume_flag = %r , multi_node_flag = %r\n' %(output_mnest, seed_value, n_pop, n_update, resume_flag, multi_node_flag))

#
# RESCALE PARAMETERS FUNCTION NEEDED BY LNLIKE
#
def trades_rescale(fitting_parameters, ndim, nparams):
  for i in range(0,ndim):
    fitting_parameters[i] = parameters_minmax[i,0] + fitting_parameters[i]*delta_parameters[i]
  return fitting_parameters

# LNPRIOR TO BE ADDED TO LOGLHD
def lnprior(fitting_parameters, ndim):
  lnprior_value = 0.
  for i in range(0, ndim):
    if('ecosw' in parameter_names[i]):
      ecc = np.sqrt(fitting_parameters[i]**2 + fitting_parameters[i+1]**2)
      if (ecc <= 0.):
        ecc = np.finfo(float).eps
      elif (ecc > 1.):
        ecc = 1.
      # ln(0.,1.) <= 0. --> -ln(0.,1.) > 0.
      lnprior_value = lnprior_value + np.abs(-np.log(ecc))
  return lnprior_value

# LNLIKEHOOD FUNCTION NEEDED BY MULTINEST
def lnlike(fitting_parameters, ndim, nparams):
  loglhd = 0.
  check = 1
  #print 'loglhd, check', loglhd, check
  trades_parameters = np.asarray([fitting_parameters[i] for i in range(0,ndim)], dtype=np.float64)
  loglhd, check = pytrades_lib.pytrades.fortran_loglikelihood(trades_parameters)
  #print 'loglhd, check', loglhd, check
  lnprior_value = lnprior(fitting_parameters, ndim)
  loglhd = loglhd + ln_err_const + lnprior_value # ln_err_const (global variable) & lnprior_value
  #print 'ln_err_const, lnprior_value, loglhd', ln_err_const, lnprior_value, loglhd
  if ( check == 0 ):
    loglhd = -0.5e10
  #sys.exit()
  return loglhd

# run MultiNest
pymultinest.run(LogLikelihood = lnlike, Prior = trades_rescale, n_dims = nfit, n_params = nfit, outputfiles_basename=output_mnest, multimodal = multi_node_flag, log_zero=log_zero_value, seed = seed_value, n_live_points = n_pop, n_iter_before_update = n_update, resume = resume_flag, verbose = True, init_MPI = mpi_onoff)

elapsed = time.time() - start
elapsed_d, elapsed_h, elapsed_m, elapsed_s = computation_time(elapsed)

print
print ' pyTRADES: pyMultiNest FINISHED in %2d day %02d hour %02d min %.2f sec - bye bye' %(int(elapsed_d), int(elapsed_h), int(elapsed_m), elapsed_s)
of_run.write(' pyTRADES: pyMultiNest FINISHED in %2d day %02d hour %02d min %.2f sec - bye bye\n' %(int(elapsed_d), int(elapsed_h), int(elapsed_m), elapsed_s))
of_run.close()
print 
pytrades_lib.pytrades.deallocate_variables()




