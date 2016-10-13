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
#from emcee.utils import MPIPool
from constants import Mjups
from ancillary import *
import logging
#from pyde.de import DiffEvol as devol

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
  parser = argparse.ArgumentParser(description='TRADES+PSO+EMCEE')
  
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
  
  ## COMPUTE OR NOT CONSTANT TO ADD TO THE LGLIKELIHOOD: - (1/2dof) * SUM( ln 2pi * sigma_obs^2 )
  #parser.add_argument('-l', '--ln-err', '--ln-err-const', action='store', dest='ln_flag', default=True, help='Computes or not constant to add to the lglikelihood: - (1/2dof) * SUM( ln 2pi * sigma_obs^2 ). Default = True. Set to False if not use this value.')
  
  cli = parser.parse_args()
  
  cli.full_path = os.path.join(os.path.abspath(cli.full_path), '')
  cli.sub_folder = os.path.join(os.path.relpath(cli.sub_folder), '')
  
  #if (cli.pso_type not in ['skip', 'run', 'exists']):
    #cli.pso_type = 'run'
  
  #cli.nthreads = int(cli.nthreads)
  
  #cli.nwalkers = int(cli.nwalkers)
  #cli.nruns = int(cli.nruns)
  #cli.npost = int(cli.npost)
  
  #cli.ln_flag = set_bool_argument(cli.ln_flag)
    
  return cli

# 
# LOGPROBABILITY FUNCTION NEEDED BY EMCEE
#
def lnprob(fitting_parameters):
  loglhd = 0.
  check = 1
  loglhd, check = pytrades_lib.pytrades.fortran_loglikelihood(fitting_parameters)
  #loglhd = loglhd + ln_err_const # ln_err_const: global variable
  if ( check == 0 ):
    loglhd = -np.inf
  return loglhd

#
# INITIALISE FOLDER AND LOG FILE
#
def init_folder(working_path, sub_folder):
  working_folder = os.path.join(working_path, sub_folder)
  if (not os.path.isdir(working_folder)):
      os.makedirs(working_folder)
  
  return working_folder

#def compute_ln_err_const(ndata, dof, e_RVo, e_T0o, ln_flag):
  #if (ln_flag):
    #eRV = e_RVo[e_RVo > 0.]
    #eT0 = e_T0o[e_T0o > 0.]
    
    #ln_e_RVo = np.sum(np.log(eRV*eRV))
    #ln_e_T0o = np.sum(np.log(eT0*eT0))

    #ln_err_const = -(0.5/dof) * (ndata*np.log(2.*np.pi) + ln_e_RVo + ln_e_T0o)
  #else:
    #ln_err_const = 0.
  #return ln_err_const


# MAIN -- TRADES + EMCEE

# ---
# initialize logger
logger = logging.getLogger("Main_log")
logger.setLevel(logging.DEBUG)
formatter = logging.Formatter("%(asctime)s - %(message)s")


# READ COMMAND LINE ARGUMENTS
cli = get_args()

# STARTING TIME
start = time.time()

# RENAME 
working_path = cli.full_path

log_file = os.path.join(working_path, '%s_log.txt' %(os.path.dirname(cli.sub_folder)))

flog = logging.FileHandler(log_file, 'w')
flog.setLevel(logging.DEBUG)
flog.setFormatter(formatter)
logger.addHandler(flog)
# log screen
slog = logging.StreamHandler()
slog.setLevel(logging.DEBUG)
slog.setFormatter(formatter)
logger.addHandler(slog)


# INITIALISE TRADES WITH SUBROUTINE WITHIN TRADES_LIB -> PARAMETER NAMES, MINMAX, INTEGRATION ARGS, READ DATA ...
pytrades_lib.pytrades.initialize_trades(working_path, cli.sub_folder)

# RETRIEVE DATA AND VARIABLES FROM TRADES_LIB MODULE
fitting_parameters = pytrades_lib.pytrades.fitting_parameters[...] # INITIAL PARAMETER SET (NEEDED ONLY TO HAVE THE PROPER ARRAY/VECTOR)
parameters_minmax = pytrades_lib.pytrades.parameters_minmax # PARAMETER BOUNDARIES
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


## compute global constant for the loglhd
#global ln_err_const

#try:
  #e_RVo = np.asarray(pytrades_lib.pytrades.ervobs[:], dtype=np.float64) # fortran variable RV in python will be rv!!!
#except:
  #e_RVo = np.asarray([0.], dtype=np.float64)
#try:
  #e_T0o = np.asarray(pytrades_lib.pytrades.et0obs[:,:], dtype=np.float64).reshape((-1))
#except:
  #e_T0o = np.asarray([0.], dtype=np.float64)
#ln_err_const = compute_ln_err_const(ndata, dof, e_RVo, e_T0o, cli.ln_flag)


# READ THE NAMES OF THE PARAMETERS FROM THE TRADES_LIB AND CONVERT IT TO PYTHON STRINGS
reshaped_names = pytrades_lib.pytrades.parameter_names.reshape((10,nfit), order='F').T
parameter_names = [''.join(reshaped_names[i,:]).strip() for i in range(0,nfit)]

# INITIALISE SCRIPT FOLDER/LOG FILE
working_folder = init_folder(working_path, cli.sub_folder)

logger.info('         ')
logger.info('======== ')
logger.info('pyTRADES' )
logger.info('======== ')
logger.info('         ')
logger.info('WORKING PATH = %s' %(working_path))
#logger.info('NUMBER OF THREADS = %d' %(nthreads))
logger.info('dof = ndata(%d) - nfit(%d) = %d' %(ndata, nfit, dof))
logger.info('Total N_RV = %d for %d set(s)' %(n_rv, n_set_rv))
logger.info('Total N_T0 = %d for %d out of %d planet(s)' %(n_t0_sum, n_set_t0, n_planets))


# P2 ecosw2 esinw2 mA2
fitting_parameters[0] = 9.286764153758
fitting_parameters[1] = -0.040366268120
fitting_parameters[2] = 0.059594807318
fitting_parameters[3] = 76.957476236345
# M3 P3 ecosw3 esinw3 mA3
fitting_parameters[4] = 0.046876837442 * Mjups
fitting_parameters[5] = 28.722915000165
fitting_parameters[6] = 0.020439443979
fitting_parameters[7] = -0.116064013999
fitting_parameters[8] = 246.552195809933
# ecosw4 esinw4 mA4
fitting_parameters[9] = -0.000579175704
fitting_parameters[10] = -0.058476768919
fitting_parameters[11] = 264.880970723681

# GOOD FIT
good_fit = fitting_parameters.copy()
#good_fit=np.zeros((nfit))

logger.info('')
logger.info('%15s %23s %23s %23s %23s' %('name', 'min', 'max', 'parameter', 'fit'))
for ii in range(0,nfit):
  logger.info('%15s %23.16E %23.16E %23.16E %23.16E' %(parameter_names[ii], parameters_minmax[ii,0], parameters_minmax[ii,1],good_fit[ii],fitting_parameters[ii]))

logger.info('')

lglhd, check_lg = pytrades_lib.pytrades.fortran_loglikelihood(fitting_parameters)
logger.info('loglikelihood = %23.16E with Check = %r' %(lglhd, check_lg))

logger.info('')

lg_prob = lnprob(fitting_parameters)
logger.info('lnprob = %23.16E' %(lg_prob))

# TEST PSO

# INITIALISE PSO ARGUMENTS FROM pso.opt FILE
logger.info('')
logger.info('INIT PSO')
pytrades_lib.pytrades.init_pso(1,working_path) # read PSO options

pso_path = os.path.join(os.path.join(working_folder, 'test_pso'), '')
pytrades_lib.pytrades.path_change(pso_path)
if(not os.path.exists(pso_path)): os.makedirs(pso_path)
logger.info('CHANGED PSO DIRECTORY TO %s' %(pso_path))
logger.info('')

all_parameters = pytrades_lib.pytrades.system_parameters
logger.info('all_parameters')
logger.info('%s' %(' '.join(['%23.16E' %(x) for x in all_parameters])))

logger.info('')
logger.info('RUN PSO')
#pytrades_lib.pytrades.driver_run_pso(1,all_parameters,fitting_parameters)
fitting_parameters, fitness_pso = pytrades_lib.pytrades.pyrun_pso(nfit,1)
logger.info('')

lglhd, check_lg = pytrades_lib.pytrades.fortran_loglikelihood(fitting_parameters)
logger.info('loglikelihood = %23.16E with Check = %r' %(lglhd, check_lg))

logger.info('')

logger.info('%15s %23s %23s %23s %23s' %('name', 'min', 'max', 'parameter', 'fit'))
for ii in range(0,nfit):
  logger.info('%15s %23.16E %23.16E %23.16E %23.16E' %(parameter_names[ii], parameters_minmax[ii,0], parameters_minmax[ii,1],good_fit[ii],fitting_parameters[ii]))

logger.info('')

elapsed = time.time() - start
elapsed_d, elapsed_h, elapsed_m, elapsed_s = computation_time(elapsed)

logger.info('')
logger.info('pyTRADES FINISHED in %2d day %02d hour %02d min %.2f sec - bye bye' %(int(elapsed_d), int(elapsed_h), int(elapsed_m), elapsed_s))
pytrades_lib.pytrades.deallocate_variables()
logger.info('')

