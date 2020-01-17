#!/usr/bin/env python
# -*- coding: utf-8 -*-

 # no more "zero" integer division bugs!:P
import argparse
import os
import numpy as np # array
import h5py
import sys
import time
import glob

# import multiprocessing as mp
from multiprocessing import Pool
# from schwimmbad import JoblibPool as Pool

import emcee

from constants import Mjups, Msear
import ancillary as anc
import pytrades_lib

# =============================================================================

def get_args():
  parser = argparse.ArgumentParser(description='TRADES+EMCEE')

  # PATH FOLDER: full_path
  parser.add_argument('-p', '--path', 
                      action='store', dest='full_path', required=True,
                      help='The path (absolute or relative) with simulation files for TRADES.'
                      )

  # SUBFOLDER TO SAVE ALL SIMULATION DATA
  parser.add_argument('-s', '--sub-folder', '--sb',
                      action='store', dest='sub_folder', default='emcee_run',
                      help='Sub-folder name, without full path. Default = emcee_run'
                      )

  # NUMBER OF CPU TO USE WITH EMCE !!
  parser.add_argument('-c', '--cpu', '--nthreads',
                      action='store', dest='nthreads', default=1,
                      help='Number of threads to use. default nthreads = 1.'
                      )

  # NUMBER OF WALKERS TO USE WITH EMCEE
  parser.add_argument('-nw', '--nwalkers', '-np', '--npop',
                      action='store', dest='nwalkers', default=1,
                      help='Number of walkers (or number of chains) to use. default nwalkers = nfit*2'
                      )

  # NUMBER OF STEPS/RUNS TO DO FOR EACH WALKER OF EMCEE
  parser.add_argument('-nr', '--nruns', '-ns', '--nsteps',
                      action='store', dest='nruns', default=10000,
                      help='Number of runs/steps to use for each chain. default nruns = 10000.'
                      )

  # NUMBER OF STEPS/RUNS TO SAVE TEMPORARY EMCEE SIMULATION
  parser.add_argument('--isave', '--iter-save', '--iterations-save',
                      action='store', dest='nsave', default='False', 
                      help='Number of iterations to do for save temporary chain. No intermediate save.'
                      )

  # COMPUTE OR NOT CONSTANT TO ADD TO THE LGLIKELIHOOD: - (1/2dof) * SUM( ln 2pi * sigma_obs^2 )
  parser.add_argument('-l', '--ln-err', '--ln-err-const',
                      action='store', dest='ln_flag', default=True,
                       help='Computes or not constant to add to the lglikelihood: - (1/2dof) * SUM( ln 2pi * sigma_obs^2 ). Default = True. Set to False if not use this value.'
                       )

  parser.add_argument('-ds', '--d-sigma', '--delta-sigma',
                      action='store', dest='delta_sigma', default='1.e-4', type=str,
                      help='Value of the sigma to compute initial walkers from initial solution. Default=1.e-4'
                      )

  parser.add_argument('-e', '--emcee-previous',
                      action='store', dest='emcee_previous', default='None', type=str,
                      help='Provide an existing "emcee_summary.hdf5" file if you wanto to start from the last step of that simulation. Default is None, create new initial walkers.'
                      )

  parser.add_argument('--trades-fin', '--trades-final-previous', '--trades-final', 
                      action='store', dest='trades_previous', default='None', 
                      help='Define file from a previous TRADES simulation. File name and structure should be of type X_Y_finalNpar.dat. The parameters from this file will be the new original parameters'
                      )

  parser.add_argument('-seed', '--seed', 
                      action='store', dest='seed', default='None', 
                      help='Seed for random number generator. Default is None.'
                      )

  cli = parser.parse_args()

  cli.full_path = os.path.join(os.path.abspath(cli.full_path), '')
  cli.sub_folder = os.path.join(os.path.relpath(cli.sub_folder), '')

  cli.nthreads = int(cli.nthreads)

  cli.nwalkers = int(cli.nwalkers)
  cli.nruns = int(cli.nruns)

  cli.ln_flag = anc.set_bool_argument(cli.ln_flag)

  cli.emcee_previous = anc.set_adhoc_file(cli.emcee_previous)
  cli.trades_previous = anc.set_adhoc_file(cli.trades_previous)

  try:
    cli.seed = int(cli.seed)
    if(cli.seed <= 0): cli.seed = None
  except:
    cli.seed = None

  return cli

# =============================================================================

# INITIALISE FOLDER AND LOG FILE
init_folder = anc.init_folder

# ==============================================================================

def get_emcee_arguments(cli,nfit):

  # NUMBER OF WALKERS
  if (cli.nwalkers < nfit*2):
    nwalkers = nfit * 2
  else:
    nwalkers = cli.nwalkers
  if (nwalkers % 2) != 0:
    nwalkers += 1

  # NUMBER OF STEPS/RUNS FOR EACH WALKER
  if (cli.nruns < 1):
    nruns = 10000
  else:
    nruns = cli.nruns

  try:
    nsave = int(cli.nsave)
    if(nsave <= 0 or nsave >= nruns):
      nsave = False
  except:
    nsave = False
  
  npost=0

  return nwalkers, nruns, nsave, npost

# =============================================================================

def compute_proper_sigma(nfit, delta_sigma, parameter_names):

  delta_sigma_out = np.ones((nfit))*delta_sigma
  for ifit in range(0,nfit):
    if(delta_sigma > 1.e-6):
      if('Ms' in parameter_names[ifit]):
        delta_sigma_out[ifit] = delta_sigma * 1.e-3
      if(delta_sigma > 1.e-2):
        if('P' in parameter_names[ifit]):
          delta_sigma_out[ifit] = delta_sigma * 1.e-2
        elif('mA' in parameter_names[ifit] or 'lambda' in parameter_names[ifit]):
          delta_sigma_out[ifit] = 1.e-4

  return delta_sigma_out

# =============================================================================

def compute_initial_walkers(lnprob_sq, nfit, nwalkers, fitting_parameters, parameters_minmax, parameter_names, delta_sigma, of_run):
  anc.print_both(' Inititializing walkers with delta_sigma = %s' %(str(delta_sigma).strip()), of_run)
  p0 = []
  i_p0 = 0

  anc.print_both(' good p0:', of_run)

  # 2017-02-03 LUCA --0--
  try:
    d_sigma = np.float64(delta_sigma)
  except:
    d_sigma = np.float64(1.e-4)

  delta_sigma_out = compute_proper_sigma(nfit, d_sigma, parameter_names)
  print(' ', end=' ')
  # init all initial walkers
  while True:
    test_p0 = np.array([fitting_parameters[ifit] + np.random.normal(loc=0., scale=delta_sigma_out[ifit]) for ifit in range(0,nfit)], dtype=np.float64)
    test_lg = lnprob_sq(test_p0, parameter_names)
    if(not np.isinf(test_lg)):
      i_p0 +=1
      p0.append(test_p0)
      print(i_p0, end=' ')
      if(i_p0 == nwalkers): break
  p0[-1] = fitting_parameters # I want the original fitting paramameters in the initial walkers
  print()
  # if 'random' opt ==> create other Gaussian starting points (<->nwalkers)
  if('ran' in str(delta_sigma).strip().lower()):
    delta_parameters = np.abs(parameters_minmax[:,1] - parameters_minmax[:,0]) # DELTA BETWEEN MAX AND MIN OF BOUNDARIES
    nw_min = 30
    n_gpts = int((nwalkers-nw_min)/nw_min) # a new Gaussian starting point each nw_min walkers, keeping at least nw_min walkers Gaussian to the original fitting parameters
    print(' new gaussian starting points: ',n_gpts)
    if(n_gpts > 0):
      print(' doing random-gaussian points ... ')
      for i_gpt in range(0, n_gpts):
        # create new starting point, but check if lnL != -inf
        new_start = fitting_parameters.copy()
        sel_fit = int(np.random.random()*(nfit-1)) # change only parameter...
        print('gpt ',i_gpt+1)
        print('selected sel_fit = ',sel_fit,' ==> ',parameter_names[sel_fit])
        print('val = ', new_start[sel_fit],' with min = ',parameters_minmax[sel_fit,0],' and delta = ',delta_parameters[sel_fit])
        while True:
          new_start[sel_fit] = parameters_minmax[sel_fit,0] + delta_parameters[sel_fit]*np.random.random()
          test_lg = lnprob_sq(new_start, parameter_names)
          if(not np.isinf(test_lg)): break
        i_pos = nw_min * i_gpt
        print('i_pos = ', end=' ')
        while True:
          test_p0 = np.array([new_start[ifit] + np.random.normal(loc=0., scale=delta_sigma_out[ifit]) for ifit in range(0,nfit)], dtype=np.float64)
          test_lg = lnprob_sq(test_p0, parameter_names)
          if(not np.isinf(test_lg)):
            p0[i_pos] = test_p0
            print(i_pos, end=' ')
            i_pos +=1
            if(i_pos%nw_min == 0): break
      print()
    print()

  anc.print_both(' done initial walkers.', of_run)

  return p0

# =============================================================================
# =============================================================================
# MAIN SCRIPT - NOT IN FUNCTION DUE TO ISSUE WITH PARALLEL AND PICKLE FUNCTION OF LNPROB..
# =============================================================================
# =============================================================================

# def main():

# MAIN -- TRADES + EMCEE
# READ COMMAND LINE ARGUMENTS
cli = get_args()

# STARTING TIME
start = time.time()

# RENAME
working_path = cli.full_path
nthreads = cli.nthreads
np.random.RandomState(cli.seed)

# INITIALISE TRADES WITH SUBROUTINE WITHIN TRADES_LIB -> PARAMETER NAMES, MINMAX, INTEGRATION ARGS, READ DATA ...
pytrades_lib.pytrades.initialize_trades(working_path, cli.sub_folder, nthreads)

# RETRIEVE DATA AND VARIABLES FROM TRADES_LIB MODULE
n_bodies = pytrades_lib.pytrades.n_bodies # NUMBER OF TOTAL BODIES OF THE SYSTEM
n_planets = n_bodies - 1 # NUMBER OF PLANETS IN THE SYSTEM
ndata = pytrades_lib.pytrades.ndata # TOTAL NUMBER OF DATA AVAILABLE
nfit  = pytrades_lib.pytrades.nfit # NUMBER OF PARAMETERS TO FIT
nfree  = pytrades_lib.pytrades.nfree # NUMBER OF FREE PARAMETERS (ie nrvset)
dof   = pytrades_lib.pytrades.dof # NUMBER OF DEGREES OF FREEDOM = NDATA - NFIT
global inv_dof
inv_dof = pytrades_lib.pytrades.inv_dof

# READ THE NAMES OF THE PARAMETERS FROM THE TRADES_LIB AND CONVERT IT TO PYTHON STRINGS
str_len = pytrades_lib.pytrades.str_len
temp_names = pytrades_lib.pytrades.get_parameter_names(nfit,str_len)
trades_names = anc.convert_fortran_charray2python_strararray(temp_names)
parameter_names = anc.trades_names_to_emcee(trades_names)

if(cli.trades_previous is not None):
  temp_names, trades_parameters = anc.read_fitted_file(cli.trades_previous)
  if(nfit != np.shape(trades_parameters)[0]):
    anc.print_both(' NUMBER OF PARAMETERS (%d) IN TRADES-PREVIOUS FILE DOES NOT' \
                'MATCH THE CURRENT CONFIGURATION nfit=%d\nSTOP' \
                %(np.shape(trades_parameters)[0], nfit)
              )
    sys.exit()
  del temp_names
else:
  # INITIAL PARAMETER SET (NEEDED ONLY TO HAVE THE PROPER ARRAY/VECTOR)
  trades_parameters = pytrades_lib.pytrades.fitting_parameters

# save initial_fitting parameters into array
original_fit_parameters = trades_parameters.copy()
fitting_parameters = anc.e_to_sqrte_fitting(trades_parameters, trades_names)

trades_minmax = pytrades_lib.pytrades.parameters_minmax # PARAMETER BOUNDARIES
parameters_minmax = anc.e_to_sqrte_boundaries(trades_minmax, trades_names)

# RADIAL VELOCITIES SET
n_rv = pytrades_lib.pytrades.nrv
n_set_rv = pytrades_lib.pytrades.nrvset

# TRANSITS SET
n_t0 = pytrades_lib.pytrades.nt0
n_t0_sum = pytrades_lib.pytrades.ntts
n_set_t0 = 0
for i in range(0, n_bodies-1):
  if (n_t0[i] > 0): n_set_t0 += 1

# compute global constant for the loglhd
global ln_err_const
ln_err_const = pytrades_lib.pytrades.ln_err_const

# SET EMCEE PARAMETERS:
nwalkers, nruns, nsave, _ = get_emcee_arguments(cli,nfit)

# uses the ancillary.get_fitted(full_path)
# to obtain info for the conversion from fitted to physical parameters
# nfit, NB, bodies_file, id_fit, id_all, nfit_list, cols_list, case =  anc.get_fitted(working_path)
_, _, _, id_fit, _, _, cols_list, case =  anc.get_fitted(working_path)
# stellar mass in solar unit to earth, priors m in Earth masses
m_factor = pytrades_lib.pytrades.mr_star[0,0]*Msear 
# read priors file: priors.in
priors = anc.read_priors(working_path)
kep_elem = anc.all_parameters_to_kep_elem(pytrades_lib.pytrades.system_parameters, n_bodies)
# lnL_priors(p, priors, names_par, kep_elem, id_fit, case_list, cols_list, m_factor)
ln_prior = anc.lnL_priors(fitting_parameters, priors, parameter_names,
                          kep_elem,
                          id_fit, case, cols_list, 
                          m_factor
                          )

#
# LOGPROBABILITY FUNCTION NEEDED BY EMCEE
#
def lnprob(fitting_parameters):
  loglhd = 0.0
  check = 1
  loglhd, check = pytrades_lib.pytrades.fortran_loglikelihood(np.array(fitting_parameters,dtype=np.float64))
  loglhd = loglhd + ln_err_const # ln_err_const: global variable
  if ( check == 0 ):
    loglhd = -np.inf
  return loglhd

def lnprob_sq(fitting_parameters, names_par):

  fitting_trades = anc.sqrte_to_e_fitting(fitting_parameters, names_par)
  loglhd = lnprob(fitting_trades)
  if(np.isfinite(loglhd)):
    ln_prior = anc.lnL_priors(fitting_parameters, priors, parameter_names,
                              kep_elem,
                              id_fit, case, cols_list, m_factor
                              )
    loglhd += ln_prior


  return loglhd

# INITIALISE SCRIPT FOLDER/LOG FILE
working_folder, _, of_run = init_folder(working_path, cli.sub_folder)

anc.print_both('',of_run)
anc.print_both(' ======== ',of_run)
anc.print_both(' pyTRADES' ,of_run)
anc.print_both(' ======== ',of_run)
anc.print_both('',of_run)
anc.print_both(' WORKING PATH = %s' %(working_path),of_run)
anc.print_both(' NUMBER OF THREADS = %d' %(nthreads),of_run)
anc.print_both(' dof = ndata(%d) - nfit(%d) - nfree(%d) = %d' %(ndata, nfit, nfree, dof),of_run)
anc.print_both(' Total N_RV = %d for %d set(s)' %(n_rv, n_set_rv),of_run)
anc.print_both(' Total N_T0 = %d for %d out of %d planet(s)' %(n_t0_sum, n_set_t0, n_planets),of_run)
anc.print_both(' %s = %.7f' %('log constant error = ', ln_err_const),of_run)
anc.print_both(' %s = %.7f' %('IN FORTRAN log constant error = ', pytrades_lib.pytrades.ln_err_const),of_run)
anc.print_both(' seed = %s' %(str(cli.seed)), of_run)

if(cli.trades_previous is not None):
  anc.print_both('\n ******\n INITIAL FITTING PARAMETERS FROM PREVIOUS' \
            ' TRADES-EMCEE SIM IN FILE:\n %s\n ******\n' %(cli.trades_previous),
            of_run
            )

anc.print_both(' ORIGINAL PARAMETER VALUES -> 0000', of_run)
_, lgllhd_0000, _ = pytrades_lib.pytrades.write_summary_files(0, original_fit_parameters)
anc.print_both(' ', of_run)
anc.print_both(' TESTING LNPROB_SQ ...', of_run)

lgllhd_zero = lnprob(trades_parameters)
lgllhd_sq_zero = lnprob_sq(fitting_parameters, parameter_names)

anc.print_both(' ', of_run)
anc.print_both(' %15s %23s %23s %15s %23s' %('trades_names', 'original_trades', 'trades_par', 'emcee_names', 'emcee_par'), of_run)
for ifit in range(0, nfit):
  anc.print_both(' %15s %23.16e %23.16e %15s %23.16e' %(trades_names[ifit], original_fit_parameters[ifit], trades_parameters[ifit], parameter_names[ifit], fitting_parameters[ifit]), of_run)
anc.print_both(' ', of_run)
anc.print_both(' %15s %23.16e %23.16e %15s %23.16e' %('lnprob', lgllhd_0000, lgllhd_zero, 'lnprob_sq', lgllhd_sq_zero), of_run)
anc.print_both(' ', of_run)

# INITIALISES THE WALKERS
if(cli.emcee_previous is not None):
  anc.print_both(' Use a previous emcee simulation: %s' %(cli.emcee_previous), of_run)
  last_p0, old_nwalkers, last_done = anc.get_last_emcee_iteration(cli.emcee_previous, nwalkers)
  if(not last_done):
    anc.print_both('**STOP: USING A DIFFERENT NUMBER OF WALKERS (%d) W.R.T. PREVIOUS EMCEE SIMULATION (%d).' %(nwalkers, old_nwalkers), of_run)
    sys.exit()
  p0 = last_p0
else:
  p0 = compute_initial_walkers(lnprob_sq, nfit, nwalkers, fitting_parameters, parameters_minmax, parameter_names, cli.delta_sigma, of_run)

anc.print_both(' emcee chain: nwalkers = %d nruns = %d' %(nwalkers, nruns), of_run)
anc.print_both(' sampler ... ',of_run)

if(nthreads > 1):
  threads_pool = Pool(nthreads)
else:
  threads_pool = None

sampler = emcee.EnsembleSampler(nwalkers, nfit, lnprob_sq,
                                pool=threads_pool,
                                args=[parameter_names]
                                ) # needed to use sqrt(e) in emcee instead of e (in fortran)

anc.print_both(' A PRE-EMCEE OF {} STEPS'.format(int(0.05*nruns)), of_run)
p0 = sampler.run_mcmc(p0, int(0.05*nruns))
anc.print_both(' RESET OF THE SAMPLER', of_run)
sampler.reset()

anc.print_both(' ready to go', of_run)
anc.print_both(' with nsave = {}'.format(nsave), of_run)
sys.stdout.flush()

if (nsave != False):
  # save temporary sampling during emcee every nruns*10%
  if(os.path.exists(os.path.join(working_folder, 'emcee_summary.hdf5')) \
    and os.path.isfile(os.path.join(working_folder, 'emcee_summary.hdf5'))):
    os.remove(os.path.join(working_folder, 'emcee_summary.hdf5'))

  f_hdf5 = h5py.File(os.path.join(working_folder, 'emcee_summary.hdf5'), 'a')
  f_hdf5.create_dataset('parameter_names', data=anc.encode_list(parameter_names), dtype='S10')
  f_hdf5.create_dataset('boundaries', data=parameters_minmax, dtype=np.float64)
  f_hdf5.create_dataset('chains', (nwalkers, nruns, nfit), dtype=np.float64)
  f_hdf5['chains'].attrs['nwalkers'] = nwalkers
  f_hdf5['chains'].attrs['nruns']    = nruns
  f_hdf5['chains'].attrs['nfit']     = nfit
  f_hdf5['chains'].attrs['nfree']    = nfree
  f_hdf5.create_dataset('lnprobability', (nwalkers, nruns), dtype=np.float64)
  f_hdf5['lnprobability'].attrs['ln_err_const'] = ln_err_const
  f_hdf5.create_dataset('acceptance_fraction', data=np.zeros((nfit)), dtype=np.float64)
  f_hdf5.create_dataset('autocor_time', data=np.zeros((nfit)), dtype=np.float64)
  f_hdf5.close()

  pos = p0
  niter_save = int(nruns/nsave)
  anc.print_both(' Running emcee with temporary saving', of_run)
  sys.stdout.flush()
  for i in range(0, niter_save):
    anc.print_both('', of_run)
    anc.print_both(' iter: {0:6d} '.format(i+1), of_run)
    aaa = i*nsave
    bbb = aaa+nsave
    pos = sampler.run_mcmc(pos, nsave)
    anc.print_both('completed {0:d} steps of {1:d}'.format(bbb, nruns), of_run)

    f_hdf5 = h5py.File(os.path.join(working_folder, 'emcee_summary.hdf5'), 'a')
    temp_dset = f_hdf5['chains']
    temp_dset[:,aaa:bbb,:] = sampler.chain[:, aaa:bbb, :]
    temp_dset.attrs['completed_steps'] = bbb
    temp_lnprob = f_hdf5['lnprobability']
    try:
      temp_lnprob[:, aaa:bbb] = sampler.lnprobability[:, aaa:bbb]
    except:
      temp_lnprob[:, aaa:bbb] = sampler.lnprobability.T[:, aaa:bbb]
    mean_acceptance_fraction = np.mean(sampler.acceptance_fraction)
    acor_time = anc.compute_acor_time(sampler, steps_done=bbb)
    temp_acor = f_hdf5['autocor_time']
    temp_acor[...] = acor_time
    f_hdf5.close()
    sys.stdout.flush()

  anc.print_both('', of_run)
  anc.print_both('...done with saving temporary total shape = {}'.format(str(np.shape(sampler.chain))), of_run)
  anc.print_both('', of_run)
  sys.stdout.flush()

else:
  # GOOD COMPLETE SINGLE RUNNING OF EMCEE, WITHOUT REMOVING THE BURN-IN
  anc.print_both(' Running full emcee ...', of_run)
  sys.stdout.flush()
  sampler.run_mcmc(p0, nruns)
  anc.print_both('done', of_run)
  anc.print_both('', of_run)
  sys.stdout.flush()

  mean_acceptance_fraction = np.mean(sampler.acceptance_fraction)
  acor_time = anc.compute_acor_time(sampler)
  lnprobability = sampler.lnprobability

  # save chains with original shape as hdf5 file
  f_hdf5 = h5py.File(os.path.join(working_folder, 'emcee_summary.hdf5'), 'w')
  f_hdf5.create_dataset('chains', data=sampler.chain, dtype=np.float64)
  f_hdf5['chains'].attrs['nwalkers']        = nwalkers
  f_hdf5['chains'].attrs['nruns']           = nruns
  f_hdf5['chains'].attrs['nfit']            = nfit
  f_hdf5['chains'].attrs['nfree']           = nfree
  f_hdf5['chains'].attrs['completed_steps'] = nruns
  f_hdf5.create_dataset('parameter_names', data=anc.encode_list(parameter_names), dtype='S10')
  f_hdf5.create_dataset('boundaries', data=parameters_minmax, dtype=np.float64)
  f_hdf5.create_dataset('acceptance_fraction', data=sampler.acceptance_fraction, dtype=np.float64)
  f_hdf5['acceptance_fraction'].attrs['mean_acceptance_fraction'] = mean_acceptance_fraction
  f_hdf5.create_dataset('autocor_time', data=acor_time, dtype=np.float64)
  f_hdf5.create_dataset('lnprobability', data=lnprobability, dtype=np.float64)
  f_hdf5['lnprobability'].attrs['ln_err_const'] = ln_err_const
  f_hdf5.close()

anc.print_both(" Mean_acceptance_fraction should be between [0.25-0.5] = {0:.6f}".format(mean_acceptance_fraction),
                of_run)
anc.print_both('', of_run)

if(threads_pool is not None):
  # close the pool of threads
  threads_pool.close()
  # threads_pool.terminate()
  threads_pool.join()

anc.print_both('COMPLETED EMCEE', of_run)

elapsed = time.time() - start
elapsed_d, elapsed_h, elapsed_m, elapsed_s = anc.computation_time(elapsed)

anc.print_both('', of_run)
anc.print_both(' pyTRADES: EMCEE FINISHED in {0:2d} day {1:02d} hour {2:02d} min {3:.2f} sec - bye bye'.format(
               int(elapsed_d), int(elapsed_h), int(elapsed_m), elapsed_s),
               of_run
              )

anc.print_both('', of_run)
of_run.close()
pytrades_lib.pytrades.deallocate_variables()

  # return

# ==============================================================================
# ==============================================================================

# if __name__ == "__main__":
#   main()
