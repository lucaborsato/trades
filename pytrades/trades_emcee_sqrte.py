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
#from emcee.utils import MPIPool
from constants import Mjups
import ancillary as anc
import glob
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
  parser = argparse.ArgumentParser(description='TRADES+EMCEE')
  
  # PATH FOLDER: full_path
  parser.add_argument('-p', '--path', action='store', dest='full_path', required=True, help='The path (absolute or relative) with simulation files for TRADES.')

  # SUBFOLDER TO SAVE ALL SIMULATION DATA
  parser.add_argument('-s', '--sub-folder', '--sb', action='store', dest='sub_folder', default='emcee_run', help='Sub-folder name, without full path. Default = emcee_run')
  
  # NUMBER OF CPU TO USE WITH EMCE !!
  parser.add_argument('-c', '--cpu', '--nthreads', action='store', dest='nthreads', default=1, help='Number of threads to use. default nthreads = 1.')
  
  # NUMBER OF WALKERS TO USE WITH EMCEE
  parser.add_argument('-nw', '--nwalkers', '-np', '--npop', action='store', dest='nwalkers', default=1, help='Number of walkers (or number of chains) to use. default nwalkers = nfit*2')
  
  # NUMBER OF STEPS/RUNS TO DO FOR EACH WALKER OF EMCEE
  parser.add_argument('-nr', '--nruns', '-ns', '--nsteps', action='store', dest='nruns', default=10000, help='Number of runs/steps to use for each chain. default nruns = 10000.')
  
  # NUMBER OF STEPS/RUNS TO SAVE TEMPORARY EMCEE SIMULATION
  parser.add_argument('--isave', '--iter-save', '--iterations-save', action='store', dest='nsave', default='False', help='Number of iterations to do for save temporary chain. default each 0.1 of nruns.')
  
  # NUMBER OF BURNIN TO SKIP FOR PRELIMINARY ANALYSIS
  parser.add_argument('-nb', '--nburn', '--npost', action='store', dest='npost', default=1000, help='Number of burn in, or number of posterior to discard at the beginning of the each chain. default npost = 1000.')
  
  # COMPUTE OR NOT CONSTANT TO ADD TO THE LGLIKELIHOOD: - (1/2dof) * SUM( ln 2pi * sigma_obs^2 )
  parser.add_argument('-l', '--ln-err', '--ln-err-const', action='store', dest='ln_flag', default=True, help='Computes or not constant to add to the lglikelihood: - (1/2dof) * SUM( ln 2pi * sigma_obs^2 ). Default = True. Set to False if not use this value.')
  
  parser.add_argument('-ds', '--d-sigma', '--delta-sigma', action='store', dest='delta_sigma', default='1.e-4', type=str, help='Value of the sigma to compute initial walkers from initial solution. Default=1.e-4')
  
  parser.add_argument('-e', '--emcee-previous', action='store', dest='emcee_previous', default='None', type=str, help='Provide an existing "emcee_summary.hdf5" file if you wanto to start from the last step of that simulation. Default is None, create new initial walkers.')
  
  parser.add_argument('--trades-fin', '--trades-final-previous', '--trades-final', action='store', dest='trades_previous', default='None', help='Define file from a previous TRADES simulation. File name and structure should be of type X_Y_finalNpar.dat. The parameters from this file will be the new original parameters')
  
  cli = parser.parse_args()
  
  cli.full_path = os.path.join(os.path.abspath(cli.full_path), '')
  cli.sub_folder = os.path.join(os.path.relpath(cli.sub_folder), '')
  
  cli.nthreads = int(cli.nthreads)
  
  cli.nwalkers = int(cli.nwalkers)
  cli.nruns = int(cli.nruns)
  cli.npost = int(cli.npost)
  
  cli.ln_flag = anc.set_bool_argument(cli.ln_flag)
  
  cli.emcee_previous = anc.set_adhoc_file(cli.emcee_previous)
  cli.trades_previous = anc.set_adhoc_file(cli.trades_previous)
  #try:
    #cli.delta_sigma = np.float64(cli.delta_sigma)
  #except:
    #cli.delta_sigma = np.float64(1.e-4)
    
  return cli


# 
# LOGPROBABILITY FUNCTION NEEDED BY EMCEE
#
def lnprob(fitting_parameters):
  loglhd = 0.
  check = 1
  loglhd, check = pytrades_lib.pytrades.fortran_loglikelihood(np.array(fitting_parameters,dtype=np.float64))
  #print loglhd, ln_err_const
  loglhd = loglhd + ln_err_const # ln_err_const: global variable
  #print loglhd
  if ( check == 0 ):
    loglhd = -np.inf
  return loglhd

def lnprob_sq(fitting_parameters, names_par):
  
  fitting_trades = anc.sqrte_to_e_fitting(fitting_parameters, names_par)
  loglhd = lnprob(fitting_trades)
  
  return loglhd

#
# INITIALISE FOLDER AND LOG FILE
#
def init_folder(working_path, sub_folder):
  working_folder = os.path.join(working_path, sub_folder)
  if (not os.path.isdir(working_folder)):
      os.makedirs(working_folder)
      
  #arg_file = os.path.join(working_path, 'arg.in')
  #shutil.copy(arg_file, os.path.join(working_folder,''))
  #bodies_file = os.path.join(working_path, 'bodies.lst')
  #shutil.copy(bodies_file, os.path.join(working_folder,''))
  #obd = open(bodies_file, 'r')
  #for line in obd.readlines():
    #shutil.copy(os.path.join(working_path, line.strip().split()[0]), os.path.join(working_folder,''))
  #obd.close()
  #t0files = glob.glob(os.path.join(working_path,'NB*_observations.dat'))
  #for t0f in t0files:
    #shutil.copy(t0f, os.path.join(working_folder,''))
  #if(os.path.exists(os.path.join(working_path,'obsRV.dat'))):
    #shutil.copy(os.path.join(working_path,'obsRV.dat'), os.path.join(working_folder,''))
  
  # copy files
  anc.copy_simulation_files(working_path, working_folder)
  
  run_log = os.path.join(working_folder, "trades_run.log")
  of_run = open(run_log, 'w')
  anc.print_both("# pyTRADES LOG FILE", of_run)
  anc.print_both("# working_path = %s" %(working_path), of_run)
  anc.print_both("# working_folder = %s" %(working_folder), of_run)
  anc.print_both("# run_log = %s" %(run_log), of_run)
  
  return working_folder, run_log, of_run

def compute_ln_err_const(dof, e_RVo, e_T0o, ln_flag=False):
  
  if (ln_flag):
    eRV = e_RVo[e_RVo > 0.]
    eT0 = e_T0o[e_T0o > 0.]
    
    ln_e_RVo = np.sum(np.log(eRV*eRV))
    ln_e_T0o = np.sum(np.log(eT0*eT0))

    ln_err_const = - 0.5 * dof * np.log(2.*np.pi) - 0.5 * ( ln_e_RVo + ln_e_T0o)
  else:
    ln_err_const = 0.
  
  return ln_err_const


def get_emcee_arguments(cli,nfit):
  
  #
  # NUMBER OF WALKERS
  if (cli.nwalkers < nfit*2):
    nwalkers = nfit * 2
  else:
    nwalkers = cli.nwalkers
  if (nwalkers % 2) != 0:
    #print ' Provided odd nwalkers = %d ==> using nwalkers+1 = %d' %(nwalkers, nwalkers+1)
    nwalkers += 1

  # NUMBER OF STEPS/RUNS FOR EACH WALKER
  if (cli.nruns < 1):
    nruns = 10000
  else:
    nruns = cli.nruns

  # NUMBER OF SAVE STEPS
  if (cli.nsave != 'False'):
    if (int(cli.nsave) > 0 and int(cli.nsave) < nruns):
      nsave = int(cli.nsave)
    elif (int(cli.nsave) <= 0):
      nsave = False
    else:
      nsave = nruns/10
  else:
    nsave = False

  # NUMBER OF BURNIN/POSTERIOR TO DISCARD
  if (cli.npost < 0):
    npost = 1000
  else:
    npost = cli.npost
  #print nwalkers, nruns, npost
  
  return nwalkers, nruns, nsave, npost

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
      

def compute_initial_walkers(nfit, nwalkers, fitting_parameters, parameters_minmax, parameter_names, delta_sigma, of_run):
  # initial walkers as input fitting_parameters + N(loc=0.,sigma=1.,size=nwalkers)*delta_sigma
  #p0 = [parameters_minmax[:,0] + np.random.random(nfit)*delta_parameters for i in range(0, nwalkers)]
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
  print ' ',
  # init all initial walkers
  while True:
      test_p0 = np.array([fitting_parameters[ifit] + np.random.normal(loc=0., scale=delta_sigma_out[ifit]) for ifit in range(0,nfit)], dtype=np.float64)
      #test_lg = lnprob(test_p0)
      test_lg = lnprob_sq(test_p0, parameter_names)
      if(not np.isinf(test_lg)):
        i_p0 +=1
        p0.append(test_p0)
        print i_p0,
        if(i_p0 == nwalkers): break
  p0[-1] = fitting_parameters # I want the original fitting paramameters in the initial walkers
  print
  # if 'random' opt ==> create other Gaussian starting points (<->nwalkers)
  if('ran' in str(delta_sigma).strip().lower()):
    delta_parameters = np.abs(parameters_minmax[:,1] - parameters_minmax[:,0]) # DELTA BETWEEN MAX AND MIN OF BOUNDARIES
    nw_min = 30
    n_gpts = int((nwalkers-nw_min)/nw_min) # a new Gaussian starting point each nw_min walkers, keeping at least nw_min walkers Gaussian to the original fitting parameters
    print ' new gaussian starting points: ',n_gpts
    if(n_gpts > 0):
      print ' doing random-gaussian points ... '
      for i_gpt in range(0, n_gpts):
        # create new starting point, but check if lnL != -inf
        new_start = fitting_parameters.copy()
        sel_fit = int(np.random.random()*(nfit-1)) # change only parameter...
        print 'gpt ',i_gpt+1
        print 'selected sel_fit = ',sel_fit,' ==> ',parameter_names[sel_fit]
        print 'val = ', new_start[sel_fit],' with min = ',parameters_minmax[sel_fit,0],' and delta = ',delta_parameters[sel_fit]
        while True:
          new_start[sel_fit] = parameters_minmax[sel_fit,0] + delta_parameters[sel_fit]*np.random.random()
          test_lg = lnprob_sq(new_start, parameter_names)
          if(not np.isinf(test_lg)): break
        i_pos = nw_min * i_gpt
        print 'i_pos = ',
        while True:
          test_p0 = np.array([new_start[ifit] + np.random.normal(loc=0., scale=delta_sigma_out[ifit]) for ifit in range(0,nfit)], dtype=np.float64)
          test_lg = lnprob_sq(test_p0, parameter_names)
          if(not np.isinf(test_lg)):
            p0[i_pos] = test_p0
            print i_pos,
            i_pos +=1
            if(i_pos%nw_min == 0): break
      print
    print
  
  anc.print_both(' done initial walkers.', of_run)

  return p0


def main():
  # MAIN -- TRADES + EMCEE
  # READ COMMAND LINE ARGUMENTS
  cli = get_args()

  # STARTING TIME
  start = time.time()

  # RENAME 
  working_path = cli.full_path
  nthreads=cli.nthreads

  # INITIALISE TRADES WITH SUBROUTINE WITHIN TRADES_LIB -> PARAMETER NAMES, MINMAX, INTEGRATION ARGS, READ DATA ...
  pytrades_lib.pytrades.initialize_trades(working_path, cli.sub_folder, nthreads)

  # RETRIEVE DATA AND VARIABLES FROM TRADES_LIB MODULE
  
  #global n_bodies, n_planets, ndata, npar, nfit, dof, inv_dof
  n_bodies = pytrades_lib.pytrades.n_bodies # NUMBER OF TOTAL BODIES OF THE SYSTEM
  n_planets = n_bodies - 1 # NUMBER OF PLANETS IN THE SYSTEM
  ndata = pytrades_lib.pytrades.ndata # TOTAL NUMBER OF DATA AVAILABLE
  npar  = pytrades_lib.pytrades.npar # NUMBER OF TOTAL PARAMATERS ~n_planets X 6
  nfit  = pytrades_lib.pytrades.nfit # NUMBER OF PARAMETERS TO FIT
  nfree  = pytrades_lib.pytrades.nfree # NUMBER OF FREE PARAMETERS (ie nrvset)
  dof   = pytrades_lib.pytrades.dof # NUMBER OF DEGREES OF FREEDOM = NDATA - NFIT
  global inv_dof
  inv_dof = np.float64(1.0 / dof)

  # READ THE NAMES OF THE PARAMETERS FROM THE TRADES_LIB AND CONVERT IT TO PYTHON STRINGS
  #reshaped_names = pytrades_lib.pytrades.parameter_names.reshape((10,nfit), order='F').T
  #parameter_names = [''.join(reshaped_names[i,:]).strip() for i in range(0,nfit)]
  
  #parameter_names = anc.convert_fortran2python_strarray(pytrades_lib.pytrades.parameter_names, nfit, str_len=10)
  #trades_names = anc.convert_fortran2python_strarray(pytrades_lib.pytrades.parameter_names,
                                                     #nfit, str_len=10
                                                    #)
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
    #fitting_parameters = pytrades_lib.pytrades.fitting_parameters
    trades_parameters = pytrades_lib.pytrades.fitting_parameters
    # save initial_fitting parameters into array
  original_fit_parameters = trades_parameters.copy()
  fitting_parameters = anc.e_to_sqrte_fitting(trades_parameters, trades_names)
  
  trades_minmax = pytrades_lib.pytrades.parameters_minmax # PARAMETER BOUNDARIES
  #parameters_minmax = trades_minmax.copy()
  #parameters_minmax[:,0] = anc.e_to_sqrte_fitting(trades_minmax[:,0], trades_names)
  #parameters_minmax[:,1] = anc.e_to_sqrte_fitting(trades_minmax[:,1], trades_names)
  parameters_minmax = anc.e_to_sqrte_boundaries(trades_minmax, trades_names)

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
    # fortran variable RV in python will be rv!!!
    e_RVo = np.array(pytrades_lib.pytrades.ervobs[:], dtype=np.float64)
  except:
    e_RVo = np.array([0.], dtype=np.float64)
  try:
    e_T0o = np.array(pytrades_lib.pytrades.et0obs[:,:], dtype=np.float64).reshape((-1))
  except:
    e_T0o = np.array([0.], dtype=np.float64)
  ln_err_const = compute_ln_err_const(dof, e_RVo, e_T0o, cli.ln_flag)

  # SET EMCEE PARAMETERS:
  nwalkers, nruns, nsave, npost = get_emcee_arguments(cli,nfit)

  # INITIALISE SCRIPT FOLDER/LOG FILE
  working_folder, run_log, of_run = init_folder(working_path, cli.sub_folder)

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

  if(cli.trades_previous is not None):
    anc.print_both('\n ******\n INITIAL FITTING PARAMETERS FROM PREVIOUS' \
              ' TRADES-EMCEE SIM IN FILE:\n %s\n ******\n' %(cli.trades_previous),
              of_run
              )
    
  anc.print_both(' ORIGINAL PARAMETER VALUES -> 0000', of_run)
  fitness_0000, lgllhd_0000, check_0000 = pytrades_lib.pytrades.write_summary_files(0, original_fit_parameters)
  anc.print_both(' ', of_run)
  anc.print_both(' TESTING LNPROB_SQ ...', of_run)
  
  lgllhd_zero = lnprob(trades_parameters)
  lgllhd_sq_zero = lnprob_sq(fitting_parameters, parameter_names)

  anc.print_both(' ', of_run)
  anc.print_both(' %15s %23s %23s %15s %23s' %('trades_names', 'original_trades', 'trades_par', 'emcee_names', 'emcee_par'), of_run)
  for ifit in range(0, nfit):
    anc.print_both(' %15s %23.16e %23.16e %15s %23.16e' %(trades_names[ifit], original_fit_parameters[ifit], trades_parameters[ifit], parameter_names[ifit], fitting_parameters[ifit]), of_run)
  anc.print_both(' ', of_run)
  anc.print_both(' %15s %23.16e %23.16e %15s %23.16e' %('lnprob', lgllhd_0000,lgllhd_zero, 'lnprob_sq', lgllhd_sq_zero), of_run)
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
    p0 = compute_initial_walkers(nfit, nwalkers, fitting_parameters, parameters_minmax, parameter_names, cli.delta_sigma, of_run)

  anc.print_both(' emcee chain: nwalkers = %d nruns = %d' %(nwalkers, nruns), of_run)
  anc.print_both(' sampler ... ',of_run)
  
  # old version with threads
  #sampler = emcee.EnsembleSampler(nwalkers, nfit, lnprob, threads=nthreads)
  #sampler = emcee.EnsembleSampler(nwalkers, nfit, lnprob_sq, threads=nthreads, args=[parameter_names]) # needed to use sqrt(e) in emcee instead of e (in fortran)
  
  threads_pool = emcee.interruptible_pool.InterruptiblePool(nthreads)
  #sampler = emcee.EnsembleSampler(nwalkers, nfit, lnprob, pool=threads_pool)
  sampler = emcee.EnsembleSampler(nwalkers, nfit, lnprob_sq, pool=threads_pool, args=[parameter_names]) # needed to use sqrt(e) in emcee instead of e (in fortran)
  
  anc.print_both(' ready to go', of_run)
  anc.print_both(' with nsave = %r' %(nsave), of_run)
  sys.stdout.flush()

  #sys.exit()

  if (nsave != False):
    # save temporary sampling during emcee every nruns*10%
    #if(os.path.exists(os.path.join(working_folder, 'emcee_temp.hdf5')) and os.path.isfile(os.path.join(working_folder, 'emcee_temp.hdf5'))):
      #os.remove(os.path.join(working_folder, 'emcee_temp.hdf5'))
    if(os.path.exists(os.path.join(working_folder, 'emcee_summary.hdf5')) and os.path.isfile(os.path.join(working_folder, 'emcee_summary.hdf5'))):
      os.remove(os.path.join(working_folder, 'emcee_summary.hdf5'))
    f_hdf5 = h5py.File(os.path.join(working_folder, 'emcee_summary.hdf5'), 'a')
    f_hdf5.create_dataset('parameter_names', data=parameter_names, dtype='S10')
    f_hdf5.create_dataset('boundaries', data=parameters_minmax, dtype=np.float64)
    temp_dset = f_hdf5.create_dataset('chains', (nwalkers, nruns, nfit), dtype=np.float64)
    temp_lnprob = f_hdf5.create_dataset('lnprobability', (nwalkers, nruns), dtype=np.float64)
    temp_lnprob.attrs['ln_err_const'] = ln_err_const
    temp_acceptance = f_hdf5.create_dataset('acceptance_fraction', data=np.zeros((nfit)), dtype=np.float64)
    temp_acor = f_hdf5.create_dataset('autocor_time', data=np.zeros((nfit)), dtype=np.float64)
    f_hdf5.close()
    pos = p0
    nchains = int(nruns/nsave)
    state=None
    anc.print_both(' Running emcee with temporary saving', of_run)
    sys.stdout.flush()
    for i in range(0, nchains):
      anc.print_both('', of_run)
      anc.print_both(' iter: %6d ' %(i+1), of_run)
      aaa = i*nsave
      bbb = aaa+nsave
      pos, prob, state = sampler.run_mcmc(pos, N=nsave, rstate0=state)
      anc.print_both('completed %d steps of %d' %(bbb, nruns), of_run)
      f_hdf5 = h5py.File(os.path.join(working_folder, 'emcee_summary.hdf5'), 'a')
      temp_dset = f_hdf5['chains'] #[:,:,:]
      temp_dset[:,aaa:bbb,:] = sampler.chain[:, aaa:bbb, :]
      #f_hdf5['chains'].attrs['completed_steps'] = bbb
      temp_dset.attrs['completed_steps'] = bbb
      temp_lnprob = f_hdf5['lnprobability'] #[:,:]
      temp_lnprob[:, aaa:bbb] = sampler.lnprobability[:, aaa:bbb]
      shape_lnprob = sampler.lnprobability.shape
      
      acceptance_fraction = sampler.acceptance_fraction
      temp_acceptance = f_hdf5['acceptance_fraction']
      temp_acceptance = acceptance_fraction
      #f_hdf5.create_dataset('acceptance_fraction', data=acceptance_fraction, dtype=np.float64)
      mean_acceptance_fraction = np.mean(acceptance_fraction)
    
      #temp_chains_T = np.zeros((bbb, nwalkers, nfit))
      #for ifit in range(0,nfit):
        #temp_chains_T[:,:,ifit] = sampler.chain[:, :bbb, ifit].T
      #acor_time = anc.compute_autocor_time(temp_chains_T, walkers_transposed=True)
      acor_time = anc.compute_acor_time(sampler, steps_done=bbb)
      temp_acor = f_hdf5['autocor_time']
      temp_acor[...] = acor_time
      
      #f_hdf5.create_dataset('autocor_time', data=np.array(acor_temp, dtype=np.float64), dtype=np.float64)
      #f_hdf5.create_dataset('autocor_time', data=np.array(sampler.acor, dtype=np.float64), dtype=np.float64) # not working
      #print 'aaa = %6d bbb = %6d -> sampler.lnprobability.shape = (%6d , %6d)' %(aaa, bbb, shape_lnprob[0], shape_lnprob[1])
      f_hdf5.close()
      sys.stdout.flush()
    anc.print_both('', of_run)
    anc.print_both('...done with saving temporary total shape = %s' %(str(np.shape(sampler.chain))), of_run)
    anc.print_both('', of_run)
    sys.stdout.flush()

  # RUN EMCEE AND RESET AFTER REMOVE BURN-IN
  #pos, prob, state = sampler.run_mcmc(p0, npost)
  #sampler.reset()
  #sampler.run_mcmc(pos, nruns, rstate0=state)
  else:
    # GOOD COMPLETE SINGLE RUNNING OF EMCEE, WITHOUT REMOVING THE BURN-IN
    anc.print_both(' Running full emcee ...', of_run)
    sys.stdout.flush()
    sampler.run_mcmc(p0, nruns)
    anc.print_both('done', of_run)
    anc.print_both('', of_run)
    sys.stdout.flush()
    flatchains = sampler.chain[:, :, :].reshape((nwalkers*nruns, nfit)) # full chain values
    acceptance_fraction = sampler.acceptance_fraction
    mean_acceptance_fraction = np.mean(acceptance_fraction)
    #autocor_time = sampler.acor
    temp_chains_T = np.zeros((bbb, nwalkers, nfit))
    for ifit in range(0,nfit):
      temp_chains_T[:,:,ifit] = sampler.chain[:, :, ifit].T
    #acor_time = anc.compute_autocor_time(temp_chains_T, walkers_transposed=True)
    acor_time = anc.compute_acor_time(sampler)
    lnprobability = sampler.lnprobability
    # save chains with original shape as hdf5 file
    f_hdf5 = h5py.File(os.path.join(working_folder, 'emcee_summary.hdf5'), 'w')
    f_hdf5.create_dataset('chains', data=sampler.chain, dtype=np.float64)
    f_hdf5['chains'].attrs['completed_steps'] = nruns
    f_hdf5.create_dataset('parameter_names', data=parameter_names, dtype='S10')
    f_hdf5.create_dataset('boundaries', data=parameters_minmax, dtype=np.float64)
    f_hdf5.create_dataset('acceptance_fraction', data=acceptance_fraction, dtype=np.float64)
    f_hdf5.create_dataset('autocor_time', data=acor_time, dtype=np.float64)
    f_hdf5.create_dataset('lnprobability', data=lnprobability, dtype=np.float64)
    f_hdf5['lnprobability'].attrs['ln_err_const'] = ln_err_const
    f_hdf5.close()

  anc.print_both(" Mean_acceptance_fraction should be between [0.25-0.5] = %.6f" %(mean_acceptance_fraction), of_run)
  anc.print_both('', of_run)

  # close the pool of threads
  threads_pool.close()
  threads_pool.terminate()
  threads_pool.join()

  anc.print_both('COMPLETED EMCEE', of_run)

  elapsed = time.time() - start
  elapsed_d, elapsed_h, elapsed_m, elapsed_s = anc.computation_time(elapsed)

  anc.print_both('', of_run)
  anc.print_both(' pyTRADES: EMCEE FINISHED in %2d day %02d hour %02d min %.2f sec - bye bye' %(int(elapsed_d), int(elapsed_h), int(elapsed_m), elapsed_s), of_run)
  anc.print_both('', of_run)
  of_run.close()
  pytrades_lib.pytrades.deallocate_variables()

  return

if __name__ == "__main__":
  main()
